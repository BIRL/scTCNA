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
from .depth import depth2_matrix, depth2_per_cell, depth_matrix, depth_per_cell, gene_mean_sd, pseudobulk_log_cpm
from .io import MALIGNANT, NON_MALIGNANT
from .metrics import METRIC_COLUMNS, gene_noise_table, log2_fold_change, mean_trend_residual
from .normalize import log2p1, normalize_cp10k, select_columns
from .qc import add_qc_metrics, common_depth_target, downsample_counts, filter_cells, filter_genes, run_scrublet
from .stats import compartment_medians, paired_noise_tests, paired_score_test

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


def process_dataset(adata: ad.AnnData, params: dict, meta: dict, out_dir: str | Path) -> dict:
    """Run the full per-dataset pipeline and write Parquet outputs. Returns the summary dict."""
    t0 = time.time()
    out_dir = Path(out_dir) / meta["id"]
    out_dir.mkdir(parents=True, exist_ok=True)
    qc, grp, norm_p, met_p, ds_p, depth_p = (params["qc"], params["groups"], params["normalisation"], params["metrics"], params["downsampling"], params["depth"])
    summary: dict = {"dataset": meta["id"], "cancer_type": meta["cancer_type"], "sctcna_version": __version__, "params": params}

    # 1. QC -----------------------------------------------------------------
    adata = adata.copy()
    adata.X = sp.csr_matrix(adata.X, dtype=np.int32)
    add_qc_metrics(adata, qc["mito_prefix"])
    adata, cell_stats = filter_cells(adata, qc["min_genes_per_cell"], qc["min_counts_per_cell"], qc["max_pct_mito"])
    adata, gene_stats = filter_genes(adata, qc["min_cells_per_gene"])
    doublet_stats = {"doublets_removed": 0, "scrublet_skipped_samples": [], "scrublet_failed_samples": []}
    if qc.get("run_scrublet", True):
        adata, doublet_stats = run_scrublet(adata, "sample", qc["scrublet_min_cells_per_sample"], params["sampling"]["seed"])
    else:
        adata.obs["doublet_score"] = np.nan
    # drop patient x compartment groups that are too small after QC
    sizes = adata.obs.groupby(["patient", "compartment"], observed=True)["cell_name"].transform("size")
    adata = adata[(sizes >= grp["min_cells_per_group"]).to_numpy()].copy()
    summary["qc"] = {**cell_stats, **gene_stats, **doublet_stats, "cells_final": int(adata.n_obs), "genes_final": int(adata.n_vars)}
    comp_counts = adata.obs.groupby(["patient", "compartment"], observed=True).size().unstack(fill_value=0)
    summary["n_patients"] = int(adata.obs["patient"].nunique())
    summary["n_patients_paired"] = int(((comp_counts.get(MALIGNANT, 0) > 0) & (comp_counts.get(NON_MALIGNANT, 0) > 0)).sum())
    summary["n_cells"] = {c: int(n) for c, n in adata.obs["compartment"].value_counts().items()}
    log.info("%s: %d cells, %d genes, %d patients after QC (%.0fs)", meta["id"], adata.n_obs, adata.n_vars, summary["n_patients"], time.time() - t0)

    # 2. normalisation --------------------------------------------------------
    counts = sp.csr_matrix(adata.X)
    norm = normalize_cp10k(counts, norm_p["target_sum"])
    genes = adata.var_names.to_numpy()

    # 3. per-gene metrics -------------------------------------------------------
    frac = met_p["lowess_frac"]
    by_pc = adata.obs.groupby(["patient", "compartment"], observed=True)
    per_patient = _metrics_for_groups(norm, by_pc, genes, grp["min_cells_per_group"], frac, ["patient", "compartment"])
    pooled = _metrics_for_groups(norm, adata.obs.groupby("compartment", observed=True), genes, grp["min_cells_per_group"], frac, ["compartment"])
    pooled["patient"] = "ALL"
    by_ct = adata.obs.groupby(["patient", "cell_type"], observed=True)
    per_celltype = _metrics_for_groups(norm, by_ct, genes, grp["min_cells_per_celltype_group"], frac, ["patient", "cell_type"])

    # 4. downsampled robustness check -----------------------------------------------
    if ds_p.get("enabled", True):
        target = common_depth_target(adata.obs["n_counts"].to_numpy(), ds_p["depth_quantile"], ds_p["min_depth"])
        counts_ds = downsample_counts(counts, target, params["sampling"]["seed"])
        norm_ds = normalize_cp10k(counts_ds, norm_p["target_sum"])
        pp_ds = _metrics_for_groups(norm_ds, by_pc, genes, grp["min_cells_per_group"], frac, ["patient", "compartment"])
        pp_ds = pp_ds[["gene", "patient", "compartment", "mean", "cv", "mad", "fano", "cv2_resid"]].rename(
            columns={"mean": "mean_ds", "cv": "cv_ds", "mad": "mad_ds", "fano": "fano_ds", "cv2_resid": "cv2_resid_ds"}
        )
        per_patient = per_patient.merge(pp_ds, on=["gene", "patient", "compartment"], how="left")
        summary["downsampling"] = {"target_depth": int(target), "cells_above_target": int((adata.obs["n_counts"] > target).sum())}
        ok = per_patient[["cv", "cv_ds"]].dropna()
        summary["downsampling"]["spearman_cv_vs_cv_ds"] = float(sps.spearmanr(ok["cv"], ok["cv_ds"]).statistic) if len(ok) > 10 else np.nan
        okr = per_patient[["cv2_resid", "cv2_resid_ds"]].dropna()
        summary["downsampling"]["spearman_cv2resid_vs_ds"] = float(sps.spearmanr(okr["cv2_resid"], okr["cv2_resid_ds"]).statistic) if len(okr) > 10 else np.nan
    else:
        for c in ("mean_ds", "cv_ds", "mad_ds", "fano_ds", "cv2_resid_ds"):
            per_patient[c] = np.nan

    for df in (per_patient, pooled, per_celltype):
        df["dataset"] = meta["id"]
        df["cancer_type"] = meta["cancer_type"]
    per_patient = per_patient[PATIENT_TABLE_COLUMNS]
    pooled = pooled[[c for c in PATIENT_TABLE_COLUMNS if c in pooled.columns]]
    per_celltype = per_celltype[["dataset", "cancer_type", "patient", "cell_type", "gene", "n_cells", "mean", "detection_rate", "sd", "cv", "cv2", "mad", "fano", "cv2_resid"]]

    # 5. DEPTH / DEPTH2 ---------------------------------------------------------------
    logX = log2p1(norm)
    det = np.bincount(sp.csr_matrix(counts).indices, minlength=adata.n_vars) / adata.n_obs
    depth_genes = det >= depth_p["min_detection_rate_genes"]
    logX_d = select_columns(logX, depth_genes)
    is_non_mal = (adata.obs["compartment"] == NON_MALIGNANT).to_numpy()
    mean_all, sd_all = gene_mean_sd(logX_d)
    cell_depth2 = depth2_per_cell(logX_d, mean_all, sd_all)
    ref_mean = np.asarray(logX_d[is_non_mal].mean(axis=0)).ravel() if is_non_mal.any() else mean_all
    cell_depth = depth_per_cell(logX_d, ref_mean)
    cell_scores = adata.obs[["cell_name", "sample", "patient", "cell_type", "compartment", "n_counts", "n_genes", "pct_mito", "doublet_score"]].copy()
    cell_scores["depth"] = cell_depth
    cell_scores["depth2"] = cell_depth2
    cell_scores["dataset"] = meta["id"]
    cell_scores["cancer_type"] = meta["cancer_type"]
    cell_scores = cell_scores.reset_index(drop=True)

    labels = adata.obs["patient"].astype(str) + "|" + adata.obs["compartment"].astype(str)
    pb, pb_n = pseudobulk_log_cpm(select_columns(counts, depth_genes), labels, min_cells=grp["min_cells_per_group"])
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
    # per-patient medians of the per-cell scores
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
    # DEPTH on pseudo-bulks is only defined for malignant samples (non-malignant pseudo-bulks are the reference),
    # so a paired test is possible for the other three scores only.
    for score in ("depth2_pseudobulk", "depth2_cell_median", "depth_cell_median"):
        summary["depth"][f"paired_test_{score}"] = paired_score_test(sample_scores, score, min_patients=3)

    # 6. gene summary + paired tests -----------------------------------------------------
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
    gene_summary = (
        medians.merge(tests_cv, on=["dataset", "gene"], how="left").merge(tests_mad, on=["dataset", "gene"], how="left").merge(tests_fano, on=["dataset", "gene"], how="left")
    )
    gene_summary["cancer_type"] = meta["cancer_type"]
    # Eq. 1 on the cancer-type (median across patients) values, for display next to the per-patient version
    gene_summary["log2fc_cv_of_medians"] = log2_fold_change(gene_summary.get(f"cv_{MALIGNANT}"), gene_summary.get(f"cv_{NON_MALIGNANT}"))
    if f"cv2_resid_{MALIGNANT}" in gene_summary and f"cv2_resid_{NON_MALIGNANT}" in gene_summary:
        gene_summary["delta_cv2_resid"] = gene_summary[f"cv2_resid_{MALIGNANT}"] - gene_summary[f"cv2_resid_{NON_MALIGNANT}"]
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
    }

    # write ----------------------------------------------------------------------------
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


def aggregate(processed_dir: str | Path, app_dir: str | Path, registry: dict, manifest: dict | None = None) -> dict:
    """Combine per-dataset outputs into the tables served by the web app."""
    processed_dir, app_dir = Path(processed_dir), Path(app_dir)
    app_dir.mkdir(parents=True, exist_ok=True)
    (app_dir / "patient").mkdir(exist_ok=True)
    summaries, gene_summaries, cell_scores, sample_scores, celltype_rows = [], [], [], [], []
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
        for c in cols[5:]:
            if c != "n_cells":
                pp[c] = pp[c].astype(np.float32)
        pp.to_parquet(app_dir / "patient" / f"{entry['id']}.parquet", index=False, compression="zstd")
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
                "depth_cell_median_diff": s["depth"]["paired_test_depth_cell_median"]["median_diff"],
                "depth_cell_p": s["depth"]["paired_test_depth_cell_median"]["p_value"],
                "download_date": (manifest or {}).get(entry["id"], {}).get("download_date"),
                "sha256": (manifest or {}).get(entry["id"], {}).get("sha256"),
                "sctcna_version": s["sctcna_version"],
                "runtime_seconds": s["runtime_seconds"],
            }
        )
    datasets = pd.DataFrame(rows)
    datasets.to_parquet(app_dir / "datasets.parquet", index=False)
    pd.concat(gene_summaries, ignore_index=True).to_parquet(app_dir / "gene_summary.parquet", index=False, compression="zstd")
    cs = pd.concat(cell_scores, ignore_index=True)
    cs[["dataset", "cancer_type", "patient", "cell_type", "compartment", "n_counts", "n_genes", "pct_mito", "depth", "depth2"]].to_parquet(
        app_dir / "cell_scores.parquet", index=False, compression="zstd"
    )
    pd.concat(sample_scores, ignore_index=True).to_parquet(app_dir / "sample_scores.parquet", index=False, compression="zstd")
    pd.concat(celltype_rows, ignore_index=True).to_parquet(app_dir / "celltype_summary.parquet", index=False, compression="zstd")
    with open(app_dir / "build.json", "w") as fh:
        json.dump({"sctcna_version": __version__, "datasets": [e["id"] for e, _ in summaries], "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, fh, indent=2)
    return {"datasets": len(summaries), "app_dir": str(app_dir)}
