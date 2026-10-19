"""End-to-end processing of one dataset and aggregation across datasets.

Per dataset (``process_dataset``):
  1. QC (cell/gene filters, mito %, Scrublet per sample)
  2. CP10K normalisation (+ log2(CP10K+1))
  3. per-gene noise metrics for every patient x compartment group, for the
     pooled compartments (patient = 'ALL') and for patient x cell-type groups
  4. the same metrics after downsampling all cells to a common UMI depth
  5. per-cell DEPTH / DEPTH2 and pseudo-bulk (patient x compartment) DEPTH / DEPTH2
  6. patient-level paired tests (malignant vs non-malignant) with BH correction

Outputs are Parquet tables in ``out_dir/<dataset_id>/`` plus a JSON summary.
``aggregate`` concatenates the per-dataset tables into the compact tables the
web app reads (``data/app``).
"""

from __future__ import annotations

import gc
import json
import logging
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats as sps

from . import __version__
from .depth import depth2_matrix, depth_matrix, depth_scores_from_counts, pseudobulk_log_cpm
from .io import MALIGNANT, NON_MALIGNANT
from .metrics import METRIC_COLUMNS, gene_noise_table, log2_fold_change, mean_trend_residual, trend_excess
from .normalize import normalize_cp10k
from .qc import add_qc_metrics, common_depth_target, downsample_counts, filter_cells, filter_genes, run_scrublet
from .stats import bh_adjust, compartment_medians, one_sample_wilcoxon, paired_noise_tests, paired_score_test

log = logging.getLogger(__name__)

PATIENT_TABLE_COLUMNS = [
    "dataset",
    "cancer_type",
    "patient",
    "compartment",
    "gene",
    "n_cells",
    "mean",
    "detection_rate",
    "sd",
    "cv",
    "cv2",
    "mad",
    "fano",
    "cv2_resid",
    "mean_ds",
    "cv_ds",
    "mad_ds",
    "fano_ds",
    "cv2_resid_ds",
]


def _group_metrics(X: sp.csr_matrix, genes: np.ndarray, frac: float) -> pd.DataFrame:
    t = gene_noise_table(X, genes)
    t["cv2_resid"], _ = mean_trend_residual(t["mean"].to_numpy(), t["cv2"].to_numpy(), frac=frac)
    return t


def _metrics_for_groups(norm: sp.csr_matrix, groups: pd.core.groupby.DataFrameGroupBy, genes: np.ndarray, min_cells: int, frac: float, keys: list[str]) -> pd.DataFrame:
    frames = []
    for key, idx in groups.indices.items():
        if len(idx) < min_cells:
            continue
        t = _group_metrics(norm[idx], genes, frac)
        if not isinstance(key, tuple):
            key = (key,)
        for k, v in zip(keys, key, strict=True):
            t[k] = v
        frames.append(t)
    if not frames:
        return pd.DataFrame(columns=["gene", *METRIC_COLUMNS, "cv2_resid", *keys])
    return pd.concat(frames, ignore_index=True)


CONTRAST_PATIENT_COLUMNS = ["patient", "reference", "n_cells_malignant", "n_cells_reference", "n_genes", "log2fc_cv", "log2fc_cv_ds", "excess_cv2_ds"]
CONTRAST_SUMMARY_COLUMNS = [
    "reference",
    "n_patients",
    "median_log2fc_cv",
    "p_log2fc_cv",
    "median_log2fc_cv_ds",
    "p_log2fc_cv_ds",
    "median_excess_cv2_ds",
    "p_excess_cv2_ds",
    "rank_biserial_excess",
    "frac_patients_excess_positive",
]
CONTRAST_Q_COLUMNS = ["q_log2fc_cv", "q_log2fc_cv_ds", "q_excess_cv2_ds"]


def noise_contrasts(per_patient: pd.DataFrame, per_celltype: pd.DataFrame, min_det: float, frac: float, min_patients: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Patient-level summaries of how noisy malignant cells are relative to a reference group.

    References are all non-malignant cells of the patient and, separately, each
    non-malignant cell type with enough cells. For every patient and reference,
    over genes detected in >= ``min_det`` of cells in both groups:

    * ``log2fc_cv``     median over genes of log2(CV_malignant / CV_reference) (TCNA Eq. 1)
    * ``log2fc_cv_ds``  the same after downsampling all cells to a common UMI depth
    * ``excess_cv2_ds`` median over genes of the malignant log10(CV^2) residual from the
                        reference's mean-CV^2 trend, at matched depth (mean-corrected)

    Each per-patient value is then tested against 0 across patients (Wilcoxon
    signed-rank), so the patient is the unit of replication. scTCNA addition.
    """
    mal = per_patient[per_patient["compartment"] == MALIGNANT].set_index(["patient", "gene"])
    refs = [(NON_MALIGNANT, per_patient[per_patient["compartment"] == NON_MALIGNANT])]
    lineages = per_celltype[per_celltype["cell_type"] != "Malignant"]
    for ct, df in lineages.groupby("cell_type", observed=True):
        refs.append((str(ct), df))
    rows = []
    for ref_name, ref_df in refs:
        for patient, r in ref_df.groupby("patient", observed=True):
            if patient not in mal.index.get_level_values(0):
                continue
            m = mal.loc[patient]
            r = r.set_index("gene")
            shared = m.index.intersection(r.index)
            m, rr = m.loc[shared], r.loc[shared]
            keep = (m["detection_rate"] >= min_det) & (rr["detection_rate"] >= min_det)
            m_e, r_e = m[keep], rr[keep]
            fc = log2_fold_change(m_e["cv"], r_e["cv"])
            fc_ds = log2_fold_change(m_e["cv_ds"], r_e["cv_ds"])
            excess = trend_excess(r["mean_ds"], r["cv_ds"] ** 2, m_e["mean_ds"], m_e["cv_ds"] ** 2, frac=frac)
            rows.append(
                {
                    "patient": patient,
                    "reference": ref_name,
                    "n_cells_malignant": int(m["n_cells"].iloc[0]) if len(m) else 0,
                    "n_cells_reference": int(r["n_cells"].iloc[0]) if len(r) else 0,
                    "n_genes": int(keep.sum()),
                    "log2fc_cv": float(np.nanmedian(fc)) if np.isfinite(fc).any() else np.nan,
                    "log2fc_cv_ds": float(np.nanmedian(fc_ds)) if np.isfinite(fc_ds).any() else np.nan,
                    "excess_cv2_ds": float(np.nanmedian(excess)) if np.isfinite(excess).any() else np.nan,
                }
            )
    per = pd.DataFrame(rows, columns=CONTRAST_PATIENT_COLUMNS)
    out = []
    for ref_name, g in per.groupby("reference", observed=True, sort=False):
        t_fc = one_sample_wilcoxon(g["log2fc_cv"], min_patients)
        t_ds = one_sample_wilcoxon(g["log2fc_cv_ds"], min_patients)
        t_ex = one_sample_wilcoxon(g["excess_cv2_ds"], min_patients)
        out.append(
            {
                "reference": ref_name,
                "n_patients": int(g["patient"].nunique()),
                "median_log2fc_cv": t_fc["median"],
                "p_log2fc_cv": t_fc["p_value"],
                "median_log2fc_cv_ds": t_ds["median"],
                "p_log2fc_cv_ds": t_ds["p_value"],
                "median_excess_cv2_ds": t_ex["median"],
                "p_excess_cv2_ds": t_ex["p_value"],
                "rank_biserial_excess": t_ex["rank_biserial"],
                "frac_patients_excess_positive": t_ex["frac_positive"],
            }
        )
    summ = pd.DataFrame(out, columns=CONTRAST_SUMMARY_COLUMNS)
    # BH over every (reference, measure) test of the dataset
    p_cols = ["p_log2fc_cv", "p_log2fc_cv_ds", "p_excess_cv2_ds"]
    q = bh_adjust(summ[p_cols].to_numpy(dtype=float).ravel()).reshape(len(summ), len(p_cols)) if len(summ) else np.empty((0, 3))
    for j, c in enumerate(p_cols):
        summ["q" + c[1:]] = q[:, j]
    return per, summ


def process_dataset(adata: ad.AnnData, params: dict, meta: dict, out_dir: str | Path) -> dict:
    """Run the full per-dataset pipeline and write Parquet outputs. Returns the summary dict."""
    t0 = time.time()
    out_dir = Path(out_dir) / meta["id"]
    out_dir.mkdir(parents=True, exist_ok=True)
    qc, grp, norm_p, met_p, ds_p, depth_p = (params["qc"], params["groups"], params["normalisation"], params["metrics"], params["downsampling"], params["depth"])
    summary: dict = {"dataset": meta["id"], "cancer_type": meta["cancer_type"], "sctcna_version": __version__, "params": params}

    # new AnnData with copied annotations; the count matrix is shared (it is only read, never modified in place)
    X_in = sp.csr_matrix(adata.X)
    if not np.issubdtype(X_in.dtype, np.integer) and not np.all(np.mod(X_in.data, 1) == 0):
        raise ValueError(f"{meta['id']}: the expression matrix is not integer UMI counts; the pipeline expects raw counts")
    adata = ad.AnnData(X=sp.csr_matrix(X_in, dtype=np.int32), obs=adata.obs.copy(), var=adata.var.copy(), uns=dict(adata.uns))
    del X_in
    add_qc_metrics(adata, qc["mito_prefix"])
    adata, cell_stats = filter_cells(adata, qc["min_genes_per_cell"], qc["min_counts_per_cell"], qc["max_pct_mito"])
    adata, gene_stats = filter_genes(adata, qc["min_cells_per_gene"])
    doublet_stats = {"doublets_removed": 0, "scrublet_skipped_samples": [], "scrublet_failed_samples": []}
    if qc.get("run_scrublet", True):
        adata, doublet_stats = run_scrublet(adata, "sample", qc["scrublet_min_cells_per_sample"], params["sampling"]["seed"])
    else:
        adata.obs["doublet_score"] = np.nan
    sizes = adata.obs.groupby(["patient", "compartment"], observed=True)["cell_name"].transform("size")
    adata = adata[(sizes >= grp["min_cells_per_group"]).to_numpy()].copy()
    # per-cell totals over the genes that survived QC (used for the downsampling target and exported)
    X_qc = sp.csr_matrix(adata.X)
    adata.obs["n_counts"] = np.asarray(X_qc.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.diff(X_qc.indptr)
    summary["qc"] = {**cell_stats, **gene_stats, **doublet_stats, "cells_final": int(adata.n_obs), "genes_final": int(adata.n_vars)}
    comp_counts = adata.obs.groupby(["patient", "compartment"], observed=True).size().unstack(fill_value=0)
    summary["n_patients"] = int(adata.obs["patient"].nunique())
    summary["n_patients_paired"] = int(((comp_counts.get(MALIGNANT, 0) > 0) & (comp_counts.get(NON_MALIGNANT, 0) > 0)).sum())
    summary["n_cells"] = {c: int(n) for c, n in adata.obs["compartment"].value_counts().items()}
    log.info("%s: %d cells, %d genes, %d patients after QC (%.0fs)", meta["id"], adata.n_obs, adata.n_vars, summary["n_patients"], time.time() - t0)

    counts = sp.csr_matrix(adata.X)
    norm = normalize_cp10k(counts, norm_p["target_sum"])
    genes = adata.var_names.to_numpy()

    frac = met_p["lowess_frac"]
    by_pc = adata.obs.groupby(["patient", "compartment"], observed=True)
    per_patient = _metrics_for_groups(norm, by_pc, genes, grp["min_cells_per_group"], frac, ["patient", "compartment"])
    pooled = _metrics_for_groups(norm, adata.obs.groupby("compartment", observed=True), genes, grp["min_cells_per_group"], frac, ["compartment"])
    pooled["patient"] = "ALL"
    by_ct = adata.obs.groupby(["patient", "cell_type"], observed=True)
    per_celltype = _metrics_for_groups(norm, by_ct, genes, grp["min_cells_per_celltype_group"], frac, ["patient", "cell_type"])

    del norm
    gc.collect()

    if ds_p.get("enabled", True):
        target = common_depth_target(adata.obs["n_counts"].to_numpy(), ds_p["depth_quantile"], ds_p["min_depth"])
        counts_ds = downsample_counts(counts, target, params["sampling"]["seed"])
        norm_ds = normalize_cp10k(counts_ds, norm_p["target_sum"])
        del counts_ds
        pp_ds = _metrics_for_groups(norm_ds, by_pc, genes, grp["min_cells_per_group"], frac, ["patient", "compartment"])
        pp_ds = pp_ds[["gene", "patient", "compartment", "mean", "cv", "mad", "fano", "cv2_resid"]].rename(
            columns={"mean": "mean_ds", "cv": "cv_ds", "mad": "mad_ds", "fano": "fano_ds", "cv2_resid": "cv2_resid_ds"}
        )
        per_patient = per_patient.merge(pp_ds, on=["gene", "patient", "compartment"], how="left")
        ct_ds = _metrics_for_groups(norm_ds, by_ct, genes, grp["min_cells_per_celltype_group"], frac, ["patient", "cell_type"])
        ct_ds = ct_ds[["gene", "patient", "cell_type", "mean", "cv"]].rename(columns={"mean": "mean_ds", "cv": "cv_ds"})
        per_celltype = per_celltype.merge(ct_ds, on=["gene", "patient", "cell_type"], how="left")
        del norm_ds
        gc.collect()
        summary["downsampling"] = {"target_depth": int(target), "cells_above_target": int((adata.obs["n_counts"] > target).sum())}
        ok = per_patient[["cv", "cv_ds"]].dropna()
        summary["downsampling"]["spearman_cv_vs_cv_ds"] = float(sps.spearmanr(ok["cv"], ok["cv_ds"]).statistic) if len(ok) > 10 else np.nan
        okr = per_patient[["cv2_resid", "cv2_resid_ds"]].dropna()
        summary["downsampling"]["spearman_cv2resid_vs_ds"] = float(sps.spearmanr(okr["cv2_resid"], okr["cv2_resid_ds"]).statistic) if len(okr) > 10 else np.nan
    else:
        for c in ("mean_ds", "cv_ds", "mad_ds", "fano_ds", "cv2_resid_ds"):
            per_patient[c] = np.nan
        per_celltype["mean_ds"] = np.nan
        per_celltype["cv_ds"] = np.nan

    for df in (per_patient, pooled, per_celltype):
        df["dataset"] = meta["id"]
        df["cancer_type"] = meta["cancer_type"]
    per_patient = per_patient[PATIENT_TABLE_COLUMNS]
    pooled = pooled[[c for c in PATIENT_TABLE_COLUMNS if c in pooled.columns]]
    per_celltype = per_celltype[
        ["dataset", "cancer_type", "patient", "cell_type", "gene", "n_cells", "mean", "detection_rate", "sd", "cv", "cv2", "mad", "fano", "cv2_resid", "mean_ds", "cv_ds"]
    ]

    det = np.bincount(counts.indices, minlength=adata.n_vars) / adata.n_obs
    depth_genes = det >= depth_p["min_detection_rate_genes"]
    is_non_mal = (adata.obs["compartment"] == NON_MALIGNANT).to_numpy()
    is_mal = (adata.obs["compartment"] == MALIGNANT).to_numpy()
    groups_pc = (adata.obs["patient"].astype(str) + "|" + adata.obs["compartment"].astype(str)).to_numpy()
    # DEPTH: malignant cells against the mean of the non-malignant cells; DEPTH2: z-scores within patient x compartment
    cell_depth, cell_depth2 = depth_scores_from_counts(counts, depth_genes, groups_pc, is_mal, is_non_mal, norm_p["target_sum"])
    cell_scores = adata.obs[["cell_name", "sample", "patient", "cell_type", "compartment", "n_counts", "n_genes", "pct_mito", "doublet_score"]].copy()
    cell_scores["depth"] = cell_depth
    cell_scores["depth2"] = cell_depth2
    cell_scores["dataset"] = meta["id"]
    cell_scores["cancer_type"] = meta["cancer_type"]
    cell_scores = cell_scores.reset_index(drop=True)

    labels = adata.obs["patient"].astype(str) + "|" + adata.obs["compartment"].astype(str)
    pb, pb_n = pseudobulk_log_cpm(counts, labels, min_cells=grp["min_cells_per_group"], columns=depth_genes)
    pb_meta = pd.DataFrame({"label": pb.columns})
    pb_meta[["patient", "compartment"]] = pb_meta["label"].str.split("|", expand=True)
    pb_meta["n_cells"] = pb_n.reindex(pb.columns).to_numpy()
    pb_meta["depth2_pseudobulk"] = depth2_matrix(pb).reindex(pb.columns).to_numpy()
    normal_flag = pd.Series(pb_meta["compartment"].eq(NON_MALIGNANT).to_numpy(), index=pb.columns)
    d_tumour = depth_matrix(pb, normal_flag)
    d_all = depth_matrix(pb, pd.Series(False, index=pb.columns)) if not normal_flag.any() else None
    pb_meta["depth_pseudobulk"] = d_tumour.reindex(pb.columns).to_numpy()
    if d_all is not None:
        pb_meta["depth_pseudobulk"] = d_all.reindex(pb.columns).to_numpy()
    cell_med = (
        cell_scores.groupby(["patient", "compartment"], observed=True)[["depth", "depth2"]]
        .median()
        .rename(columns={"depth": "depth_cell_median", "depth2": "depth2_cell_median"})
        .reset_index()
    )
    sample_scores = pb_meta.drop(columns="label").merge(cell_med, on=["patient", "compartment"], how="left")
    sample_scores["dataset"] = meta["id"]
    sample_scores["cancer_type"] = meta["cancer_type"]
    summary["depth"] = {"n_genes_used": int(depth_genes.sum())}
    # DEPTH is only defined for malignant samples (the non-malignant cells are its reference), so only the two
    # DEPTH2 summaries can be compared between compartments.
    for score in ("depth2_pseudobulk", "depth2_cell_median"):
        summary["depth"][f"paired_test_{score}"] = paired_score_test(sample_scores, score, min_patients=grp["min_patients_per_test"])

    medians = compartment_medians(per_patient)
    # Genes enter the paired tests (and hence the BH family) only when detected in >= min_detection_rate of
    # cells in both compartments (median across patients); noise of sporadically detected genes is unreliable.
    det_cols = [c for c in (f"detection_rate_{MALIGNANT}", f"detection_rate_{NON_MALIGNANT}") if c in medians]
    medians["expressed"] = (medians[det_cols] >= met_p["min_detection_rate"]).all(axis=1) if len(det_cols) == 2 else False
    tested = per_patient[per_patient["gene"].isin(set(medians.loc[medians["expressed"], "gene"]))]
    tests_cv = paired_noise_tests(tested, "cv", grp["min_patients_per_test"])
    tests_mad = paired_noise_tests(tested, "mad", grp["min_patients_per_test"])[["dataset", "gene", "log2fc_mad_median", "p_value", "q_value"]].rename(
        columns={"p_value": "p_value_mad", "q_value": "q_value_mad"}
    )
    tests_fano = paired_noise_tests(tested, "fano", grp["min_patients_per_test"])[["dataset", "gene", "log2fc_fano_median", "p_value", "q_value"]].rename(
        columns={"p_value": "p_value_fano", "q_value": "q_value_fano"}
    )
    # mean-corrected CV^2: paired per-patient differences of the LOWESS residuals (already on a log10 scale)
    tests_resid = paired_noise_tests(tested, "cv2_resid", grp["min_patients_per_test"], difference=True)[
        ["dataset", "gene", "delta_cv2_resid_median", "p_value", "q_value"]
    ].rename(columns={"delta_cv2_resid_median": "delta_cv2_resid", "p_value": "p_value_cv2_resid", "q_value": "q_value_cv2_resid"})
    gene_summary = medians
    for t in (tests_cv, tests_mad, tests_fano, tests_resid):
        gene_summary = gene_summary.merge(t, on=["dataset", "gene"], how="left")
    gene_summary["cancer_type"] = meta["cancer_type"]
    # Eq. 1 on the cancer-type (median across patients) values, for display next to the per-patient version
    gene_summary["log2fc_cv_of_medians"] = log2_fold_change(gene_summary.get(f"cv_{MALIGNANT}"), gene_summary.get(f"cv_{NON_MALIGNANT}"))
    # pooled-vs-patient comparison (pseudoreplication diagnostic): CV from pooled cells vs median across patients
    pooled_cv = pooled.pivot_table(index="gene", columns="compartment", values="cv", aggfunc="first", observed=True)
    pooled_cv.columns = [f"cv_pooled_{c}" for c in pooled_cv.columns]
    gene_summary = gene_summary.merge(pooled_cv.reset_index(), on="gene", how="left")
    summary["n_genes_expressed"] = int(gene_summary["expressed"].sum())
    sig = gene_summary[gene_summary["expressed"]]
    summary["tests"] = {
        "n_genes_tested": int(sig["p_value"].notna().sum()),
        "n_sig_q05": int((sig["q_value"] < 0.05).sum()),
        "n_sig_q05_noisier_malignant": int(((sig["q_value"] < 0.05) & (sig["log2fc_cv_median"] > 0)).sum()),
        "n_sig_q05_noisier_non_malignant": int(((sig["q_value"] < 0.05) & (sig["log2fc_cv_median"] < 0)).sum()),
        "median_log2fc_cv_expressed_genes": float(sig["log2fc_cv_median"].median()) if len(sig) else np.nan,
        "median_delta_cv2_resid_expressed_genes": float(sig["delta_cv2_resid"].median()) if "delta_cv2_resid" in sig and len(sig) else np.nan,
        "n_sig_q05_cv2_resid": int((sig["q_value_cv2_resid"] < 0.05).sum()),
    }

    contrast_patient, contrast_summary = noise_contrasts(per_patient, per_celltype, met_p["min_detection_rate"], frac, grp["min_patients_per_test"])
    for df in (contrast_patient, contrast_summary):
        df.insert(0, "cancer_type", meta["cancer_type"])
        df.insert(0, "dataset", meta["id"])
    allref = contrast_summary[contrast_summary["reference"] == NON_MALIGNANT]
    summary["contrast_vs_non_malignant"] = allref.drop(columns=["dataset", "cancer_type", "reference"]).iloc[0].to_dict() if len(allref) else {}

    contrast_patient.to_parquet(out_dir / "noise_contrasts_patient.parquet", index=False, compression="zstd")
    contrast_summary.to_parquet(out_dir / "noise_contrasts.parquet", index=False, compression="zstd")
    per_patient.to_parquet(out_dir / "gene_metrics_patient.parquet", index=False, compression="zstd")
    pooled.to_parquet(out_dir / "gene_metrics_pooled.parquet", index=False, compression="zstd")
    per_celltype.to_parquet(out_dir / "gene_metrics_celltype.parquet", index=False, compression="zstd")
    gene_summary.to_parquet(out_dir / "gene_summary.parquet", index=False, compression="zstd")
    cell_scores.to_parquet(out_dir / "cell_scores.parquet", index=False, compression="zstd")
    sample_scores.to_parquet(out_dir / "sample_scores.parquet", index=False, compression="zstd")
    summary["runtime_seconds"] = round(time.time() - t0, 1)
    with open(out_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=_json_default)
    log.info("%s: done in %.0fs", meta["id"], summary["runtime_seconds"])
    return summary


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


APP_MIN_DETECTION = 0.01


def _compact(df: pd.DataFrame, categorical: list[str], float_columns: list[str] | None = None) -> pd.DataFrame:
    """Shrink an app table: repeated strings become categoricals and display-only floats float32.

    The web app holds every table in memory (512 MB on the free hosting tier),
    so per-cell and per-patient tables are stored compactly; values used for
    statistics (p/q-values in gene_summary) stay float64.
    """
    df = df.copy()
    for c in categorical:
        if c in df:
            df[c] = df[c].astype(str).astype("category")
    for c in float_columns or []:
        if c in df:
            df[c] = df[c].astype(np.float32)
    return df


def heterogeneity_tests(sample_scores: pd.DataFrame, min_patients: int = 6) -> pd.DataFrame:
    """Paired malignant vs non-malignant tests of the per-patient DEPTH2 summaries, per cancer type.

    DEPTH has no non-malignant values (those cells are its reference), so it is not tested.
    """
    rows = []
    for ct, df in sample_scores.groupby("cancer_type", observed=True):
        for level, column in (("cell", "depth2_cell_median"), ("sample", "depth2_pseudobulk")):
            rows.append({"cancer_type": ct, "score": "depth2", "level": level, **paired_score_test(df, column, min_patients=min_patients)})
    return pd.DataFrame(rows, columns=["cancer_type", "score", "level", "n_pairs", "median_diff", "p_value", "rank_biserial"])


def aggregate(processed_dir: str | Path, app_dir: str | Path, registry: dict, manifest: dict | None = None) -> dict:
    """Combine per-dataset outputs into the tables served by the web app."""
    processed_dir, app_dir = Path(processed_dir), Path(app_dir)
    app_dir.mkdir(parents=True, exist_ok=True)
    (app_dir / "patient").mkdir(exist_ok=True)
    summaries, gene_summaries, cell_scores, sample_scores, celltype_rows = [], [], [], [], []
    contrasts, contrasts_patient = [], []
    for entry in registry["datasets"]:
        d = processed_dir / entry["id"]
        if not (d / "summary.json").exists():
            continue
        s = json.load(open(d / "summary.json"))
        summaries.append((entry, s))
        gs = pd.read_parquet(d / "gene_summary.parquet")
        gene_summaries.append(gs)
        cell_scores.append(pd.read_parquet(d / "cell_scores.parquet"))
        sample_scores.append(pd.read_parquet(d / "sample_scores.parquet"))
        pp = pd.read_parquet(d / "gene_metrics_patient.parquet")
        keep_genes = set(gs.loc[gs["expressed"], "gene"])
        pp = pp[pp["gene"].isin(keep_genes)]
        cols = ["dataset", "cancer_type", "patient", "compartment", "gene", "n_cells", "mean", "detection_rate", "cv", "mad", "fano", "cv2_resid", "cv_ds", "cv2_resid_ds"]
        pp = pp[cols].copy()
        pp = _compact(pp, ["dataset", "cancer_type", "patient", "compartment", "gene"], [c for c in cols[6:]])
        pp.to_parquet(app_dir / "patient" / f"{entry['id']}.parquet", index=False, compression="zstd")
        contrasts.append(pd.read_parquet(d / "noise_contrasts.parquet"))
        contrasts_patient.append(pd.read_parquet(d / "noise_contrasts_patient.parquet"))
        ct = pd.read_parquet(d / "gene_metrics_celltype.parquet")
        ct = ct[ct["gene"].isin(keep_genes)]
        ct_med = (
            ct.groupby(["dataset", "cancer_type", "cell_type", "gene"], observed=True)
            .agg(
                n_patients=("patient", "nunique"),
                mean=("mean", "median"),
                detection_rate=("detection_rate", "median"),
                cv=("cv", "median"),
                mad=("mad", "median"),
                fano=("fano", "median"),
                cv2_resid=("cv2_resid", "median"),
            )
            .reset_index()
        )
        celltype_rows.append(ct_med)
    if not summaries:
        raise RuntimeError(f"no processed datasets found under {processed_dir}")
    rows = []
    for entry, s in summaries:
        rows.append(
            {
                "dataset": entry["id"],
                "cancer_type": entry["cancer_type"],
                "cancer_name": entry["cancer_name"],
                "tcna_site": entry.get("tcna_site"),
                "source": entry["source"],
                "study": entry["study"],
                "citation": entry["citation"],
                "technology": entry["technology"],
                "site": entry.get("site"),
                "data_url": entry["data_url"],
                "meta_url": entry.get("meta_url"),
                "notes": entry.get("notes"),
                "n_patients": s["n_patients"],
                "n_patients_paired": s["n_patients_paired"],
                "n_cells_malignant": s["n_cells"].get(MALIGNANT, 0),
                "n_cells_non_malignant": s["n_cells"].get(NON_MALIGNANT, 0),
                "n_genes": s["qc"]["genes_final"],
                "n_genes_expressed": s["n_genes_expressed"],
                "doublets_removed": s["qc"]["doublets_removed"],
                "downsampling_target_depth": s.get("downsampling", {}).get("target_depth"),
                "spearman_cv_vs_downsampled": s.get("downsampling", {}).get("spearman_cv_vs_cv_ds"),
                "n_sig_q05": s["tests"]["n_sig_q05"],
                "n_sig_noisier_malignant": s["tests"]["n_sig_q05_noisier_malignant"],
                "n_sig_noisier_non_malignant": s["tests"]["n_sig_q05_noisier_non_malignant"],
                "median_log2fc_cv": s["tests"]["median_log2fc_cv_expressed_genes"],
                "median_delta_cv2_resid": s["tests"].get("median_delta_cv2_resid_expressed_genes"),
                "depth2_pseudobulk_median_diff": s["depth"]["paired_test_depth2_pseudobulk"]["median_diff"],
                "depth2_pseudobulk_p": s["depth"]["paired_test_depth2_pseudobulk"]["p_value"],
                "depth2_cell_median_diff": s["depth"]["paired_test_depth2_cell_median"]["median_diff"],
                "depth2_cell_p": s["depth"]["paired_test_depth2_cell_median"]["p_value"],
                "log2fc_cv_ds_vs_non_malignant": s.get("contrast_vs_non_malignant", {}).get("median_log2fc_cv_ds"),
                "excess_cv2_ds_vs_non_malignant": s.get("contrast_vs_non_malignant", {}).get("median_excess_cv2_ds"),
                "excess_cv2_ds_p": s.get("contrast_vs_non_malignant", {}).get("p_excess_cv2_ds"),
                "download_date": (manifest or {}).get(entry["id"], {}).get("download_date"),
                "sha256": (manifest or {}).get(entry["id"], {}).get("sha256"),
                "sctcna_version": s["sctcna_version"],
                "runtime_seconds": s["runtime_seconds"],
            }
        )
    datasets = pd.DataFrame(rows)
    datasets.to_parquet(app_dir / "datasets.parquet", index=False)
    # app copy: genes detected in >= APP_MIN_DETECTION of cells (median over patients) in at least one compartment,
    # float32 display values and categorical labels, so the served tables fit a 512 MB instance; the full
    # per-dataset tables stay in data/processed/<dataset>/gene_summary.parquet
    gs_all = pd.concat(gene_summaries, ignore_index=True)
    detected = gs_all[["detection_rate_malignant", "detection_rate_non_malignant"]].max(axis=1) >= APP_MIN_DETECTION
    gs_app = gs_all[detected | gs_all["expressed"].fillna(False).astype(bool)].reset_index(drop=True)
    stat_columns = [c for c in gs_app.columns if c.startswith(("p_value", "q_value"))]
    float_columns = [c for c in gs_app.select_dtypes("float").columns if c not in stat_columns]
    _compact(gs_app, ["dataset", "cancer_type"], float_columns).to_parquet(app_dir / "gene_summary.parquet", index=False, compression="zstd")
    cs = _compact(pd.concat(cell_scores, ignore_index=True), ["dataset", "cancer_type", "patient", "cell_type", "compartment"], ["n_counts", "pct_mito", "depth", "depth2"])
    cs[["dataset", "cancer_type", "patient", "cell_type", "compartment", "n_counts", "n_genes", "pct_mito", "depth", "depth2"]].to_parquet(
        app_dir / "cell_scores.parquet", index=False, compression="zstd"
    )
    ss = pd.concat(sample_scores, ignore_index=True)
    ss.to_parquet(app_dir / "sample_scores.parquet", index=False, compression="zstd")
    min_patients = registry.get("params", {}).get("groups", {}).get("min_patients_per_test", 6)
    heterogeneity_tests(ss, min_patients).to_parquet(app_dir / "heterogeneity_tests.parquet", index=False)
    cts = _compact(pd.concat(celltype_rows, ignore_index=True), ["dataset", "cancer_type", "cell_type", "gene"], ["mean", "detection_rate", "cv", "mad", "fano", "cv2_resid"])
    cts.to_parquet(app_dir / "celltype_summary.parquet", index=False, compression="zstd")
    pd.concat(contrasts, ignore_index=True).to_parquet(app_dir / "noise_contrasts.parquet", index=False, compression="zstd")
    pd.concat(contrasts_patient, ignore_index=True).to_parquet(app_dir / "noise_contrasts_patient.parquet", index=False, compression="zstd")
    with open(app_dir / "build.json", "w") as fh:
        json.dump({"sctcna_version": __version__, "datasets": [e["id"] for e, _ in summaries], "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, fh, indent=2)
    return {"datasets": len(summaries), "app_dir": str(app_dir)}
