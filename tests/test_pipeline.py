"""Integration tests: the synthetic datasets through process_dataset and aggregate (conftest fixtures), plus the CLI."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml
from scipy import stats as sps

from sctcna import __version__
from sctcna.config import load_config
from sctcna.depth import depth2_matrix, depth2_per_cell, depth_matrix, gene_mean_sd, pseudobulk_log_cpm
from sctcna.metrics import gene_noise_table
from sctcna.normalize import log2p1, normalize_cp10k
from sctcna.pipeline import CONTRAST_PATIENT_COLUMNS, CONTRAST_SUMMARY_COLUMNS, PATIENT_TABLE_COLUMNS, process_dataset
from sctcna.qc import add_qc_metrics, filter_cells, filter_genes
from sctcna.synthetic import make_synthetic_adata

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATASETS = ["SYN_Test", "SYN_Two"]  # ids and patient counts of the conftest fixtures
N_PATIENTS = {"SYN_Test": 6, "SYN_Two": 5}
FIRST = "SYN_Test"
META = {"id": "SYN_Test", "cancer_type": "SYN"}

GENE_SUMMARY_COLUMNS = [
    "dataset", "gene", "mean_malignant", "mean_non_malignant", "detection_rate_malignant", "detection_rate_non_malignant",
    "sd_malignant", "sd_non_malignant", "cv_malignant", "cv_non_malignant", "cv2_malignant", "cv2_non_malignant",
    "mad_malignant", "mad_non_malignant", "fano_malignant", "fano_non_malignant", "cv2_resid_malignant", "cv2_resid_non_malignant",
    "n_patients_malignant", "n_patients_non_malignant", "n_pairs", "log2fc_cv_median", "log2fc_cv_mean", "wilcoxon_stat",
    "p_value", "rank_biserial", "q_value", "log2fc_mad_median", "p_value_mad", "q_value_mad", "log2fc_fano_median",
    "p_value_fano", "q_value_fano", "cancer_type", "log2fc_cv_of_medians", "delta_cv2_resid", "expressed",
    "cv_pooled_malignant", "cv_pooled_non_malignant",
]  # fmt: skip
SAMPLE_SCORE_COLUMNS = ["patient", "compartment", "n_cells", "depth2_pseudobulk", "depth_pseudobulk", "depth_cell_median", "depth2_cell_median", "dataset", "cancer_type"]
APP_COLUMNS = {
    "datasets.parquet": [
        "dataset", "cancer_type", "cancer_name", "tcna_site", "source", "study", "citation", "technology", "site", "data_url", "meta_url", "notes",
        "n_patients", "n_patients_paired", "n_cells_malignant", "n_cells_non_malignant", "n_genes", "n_genes_expressed", "doublets_removed",
        "downsampling_target_depth", "spearman_cv_vs_downsampled", "n_sig_q05", "n_sig_noisier_malignant", "n_sig_noisier_non_malignant",
        "median_log2fc_cv", "median_delta_cv2_resid", "depth2_pseudobulk_median_diff", "depth2_pseudobulk_p", "depth2_cell_median_diff",
        "depth2_cell_p", "depth_cell_median_diff", "depth_cell_p", "log2fc_cv_ds_vs_non_malignant", "excess_cv2_ds_vs_non_malignant",
        "excess_cv2_ds_p", "download_date", "sha256", "sctcna_version", "runtime_seconds",
    ],
    "gene_summary.parquet": GENE_SUMMARY_COLUMNS,
    "cell_scores.parquet": ["dataset", "cancer_type", "patient", "cell_type", "compartment", "n_counts", "n_genes", "pct_mito", "depth", "depth2"],
    "sample_scores.parquet": SAMPLE_SCORE_COLUMNS,
    "celltype_summary.parquet": ["dataset", "cancer_type", "cell_type", "gene", "n_patients", "mean", "detection_rate", "cv", "mad", "fano", "cv2_resid"],
    "noise_contrasts.parquet": ["dataset", "cancer_type", *CONTRAST_SUMMARY_COLUMNS],
    "noise_contrasts_patient.parquet": ["dataset", "cancer_type", *CONTRAST_PATIENT_COLUMNS],
}  # fmt: skip
PATIENT_APP_COLUMNS = ["dataset", "cancer_type", "patient", "compartment", "gene", "n_cells", "mean", "detection_rate", "cv", "mad", "fano", "cv2_resid", "cv_ds", "cv2_resid_ds"]
PROCESSED_COLUMNS = {
    "gene_metrics_patient.parquet": PATIENT_TABLE_COLUMNS,
    "gene_metrics_pooled.parquet": [c for c in PATIENT_TABLE_COLUMNS if not c.endswith("_ds")],
    "gene_metrics_celltype.parquet": ["dataset", "cancer_type", "patient", "cell_type", "gene", "n_cells", "mean", "detection_rate", "sd", "cv", "cv2", "mad", "fano", "cv2_resid", "mean_ds", "cv_ds"],
    "noise_contrasts.parquet": ["dataset", "cancer_type", *CONTRAST_SUMMARY_COLUMNS],
    "noise_contrasts_patient.parquet": ["dataset", "cancer_type", *CONTRAST_PATIENT_COLUMNS],
    "gene_summary.parquet": GENE_SUMMARY_COLUMNS,
    "cell_scores.parquet": ["cell_name", "sample", "patient", "cell_type", "compartment", "n_counts", "n_genes", "pct_mito", "doublet_score", "depth", "depth2", "dataset", "cancer_type"],
    "sample_scores.parquet": SAMPLE_SCORE_COLUMNS,
}  # fmt: skip


def _assert_columns(table: pd.DataFrame, columns: list[str], name: str) -> None:
    """The documented columns, all present and nothing else (the order is not part of the contract)."""
    assert not table.columns.duplicated().any(), name
    assert set(table.columns) == set(columns), f"{name}: {sorted(set(table.columns) ^ set(columns))}"


def _qc_like_pipeline(adata, params: dict):
    """Replay the QC of process_dataset (Scrublet off) so per-cell values can be recomputed from the fixture."""
    a = adata.copy()
    a.X = sp.csr_matrix(a.X, dtype=np.int32)
    qc = params["qc"]
    add_qc_metrics(a, qc["mito_prefix"])
    a, _ = filter_cells(a, qc["min_genes_per_cell"], qc["min_counts_per_cell"], qc["max_pct_mito"])
    a, _ = filter_genes(a, qc["min_cells_per_gene"])
    sizes = a.obs.groupby(["patient", "compartment"], observed=True)["cell_name"].transform("size")
    return a[(sizes >= params["groups"]["min_cells_per_group"]).to_numpy()].copy()


def _depth_gene_mask(adata, params: dict) -> np.ndarray:
    counts = sp.csr_matrix(adata.X)
    detection = np.bincount(counts.indices, minlength=adata.n_vars) / adata.n_obs
    return detection >= params["depth"]["min_detection_rate_genes"]


# ---------------------------------------------------------------------------
# process_dataset + aggregate on the conftest fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def gene_summary(app_data_dir, synthetic_adata) -> pd.DataFrame:
    """Gene summary of the first dataset with the planted gene class and mean count of every gene."""
    gs = pd.read_parquet(app_data_dir / "gene_summary.parquet")
    gs = gs[gs["dataset"] == FIRST].copy()
    gs["truth"] = gs["gene"].map(synthetic_adata.var["truth"])
    gs["mu"] = gs["gene"].map(synthetic_adata.var["mu"])
    assert gs["truth"].notna().all()
    return gs


def test_processed_outputs_exist_with_documented_columns(processed_dir):
    for ds in DATASETS:
        d = processed_dir / ds
        for name, cols in PROCESSED_COLUMNS.items():
            assert (d / name).exists(), name
            table = pd.read_parquet(d / name)
            _assert_columns(table, cols, name)
            assert len(table) > 0, name
            assert set(table["dataset"]) == {ds}, name
        pooled = pd.read_parquet(d / "gene_metrics_pooled.parquet")
        assert set(pooled["patient"]) == {"ALL"}
        cells = pd.read_parquet(d / "cell_scores.parquet")
        assert cells["doublet_score"].isna().all()  # Scrublet is off in the fixture parameters
        summary = json.loads((d / "summary.json").read_text())
        assert summary["dataset"] == ds and summary["sctcna_version"] == __version__
        for key in ("params", "qc", "n_patients", "n_patients_paired", "n_cells", "downsampling", "depth", "tests", "n_genes_expressed", "runtime_seconds"):
            assert key in summary, key
        for score in ("depth2_pseudobulk", "depth2_cell_median", "depth_cell_median"):
            assert set(summary["depth"][f"paired_test_{score}"]) == {"n_pairs", "median_diff", "p_value", "rank_biserial"}
        assert summary["qc"]["doublets_removed"] == 0
        assert summary["n_patients"] == summary["n_patients_paired"] == N_PATIENTS[ds]


def test_app_tables_have_documented_columns(app_data_dir):
    for name, cols in APP_COLUMNS.items():
        table = pd.read_parquet(app_data_dir / name)
        _assert_columns(table, cols, name)
        assert set(table["dataset"]) == set(DATASETS), name
    for ds in DATASETS:
        patient = pd.read_parquet(app_data_dir / "patient" / f"{ds}.parquet")
        _assert_columns(patient, PATIENT_APP_COLUMNS, f"patient/{ds}")
        assert set(patient["dataset"]) == {ds}
    build = json.loads((app_data_dir / "build.json").read_text())
    assert build["sctcna_version"] == __version__
    assert build["datasets"] == DATASETS
    assert build["built"]


def test_datasets_table(app_data_dir, processed_dir):
    datasets = pd.read_parquet(app_data_dir / "datasets.parquet").set_index("dataset")
    assert len(datasets) == 2 and list(datasets.index) == DATASETS
    assert datasets["cancer_type"].tolist() == ["SYN", "SY2"]
    assert datasets["n_patients"].tolist() == [6, 5]
    assert datasets["n_patients_paired"].tolist() == [6, 5]
    assert (datasets["doublets_removed"] == 0).all()
    assert (datasets["n_genes_expressed"] <= datasets["n_genes"]).all()
    assert (datasets["spearman_cv_vs_downsampled"] > 0.9).all()
    assert (datasets["sctcna_version"] == __version__).all()
    assert datasets["download_date"].isna().all() and datasets["sha256"].isna().all()  # no manifest in the fixture
    cells = pd.read_parquet(app_data_dir / "cell_scores.parquet")
    per_group = cells.groupby(["dataset", "compartment"]).size()
    for ds in DATASETS:
        assert datasets.loc[ds, "n_cells_malignant"] == per_group[(ds, "malignant")]
        assert datasets.loc[ds, "n_cells_non_malignant"] == per_group[(ds, "non_malignant")]
        summary = json.loads((processed_dir / ds / "summary.json").read_text())
        assert datasets.loc[ds, "downsampling_target_depth"] == summary["downsampling"]["target_depth"]
        assert datasets.loc[ds, "n_genes"] == summary["qc"]["genes_final"]
        assert datasets.loc[ds, "n_sig_q05"] == summary["tests"]["n_sig_q05"]
        np.testing.assert_equal(datasets.loc[ds, "depth2_pseudobulk_p"], summary["depth"]["paired_test_depth2_pseudobulk"]["p_value"])


def test_planted_noisy_genes_are_the_top_differential_noise_hits(gene_summary):
    noisy = gene_summary[gene_summary["truth"] == "noisy_malignant"]
    others = gene_summary[gene_summary["truth"] != "noisy_malignant"]
    assert len(noisy) == 20
    tested = noisy[noisy["expressed"]]  # paired statistics exist for expressed genes only
    assert (tested["log2fc_cv_median"] > 0).all()
    assert (tested["rank_biserial"] > 0).all()
    # a planted mean of >= 1 count per cell gives a clear signal: expected log2(CV_mal / CV_non-mal) = 0.5 * log2(1 + mu)
    strong = noisy[noisy["mu"] >= 1.0]
    assert len(strong) >= 10
    assert strong["expressed"].all()
    assert strong["log2fc_cv_median"].min() > others["log2fc_cv_median"].max()
    assert (strong["rank_biserial"] == 1.0).all()
    assert (strong["p_value"] == 2 / 64).all()  # six patients, every pair positive: the smallest two-sided Wilcoxon p
    assert strong["p_value"].max() <= gene_summary["p_value"].min()
    assert (strong["q_value"] <= others["q_value"].min()).all()
    top = gene_summary.nlargest(len(strong), "log2fc_cv_median")
    assert (top["truth"] == "noisy_malignant").all()
    assert (strong["log2fc_cv_of_medians"] > 0).all()
    assert (strong["delta_cv2_resid"] > 0).all()
    assert (strong["log2fc_mad_median"] > 0).all() and (strong["log2fc_fano_median"] > 0).all()


def test_constant_genes_have_zero_raw_noise_and_the_lowest_normalised_noise(gene_summary, synthetic_adata):
    const = gene_summary[gene_summary["truth"] == "constant"]
    others = gene_summary[gene_summary["truth"] != "constant"]
    assert len(const) == 5
    raw = gene_noise_table(synthetic_adata.X, synthetic_adata.var_names).set_index("gene")
    assert (raw.loc[const["gene"], ["sd", "cv", "cv2", "mad", "fano"]] == 0.0).all().all()
    # after per-cell CP10K scaling a constant count varies with the library size, so the CV is small but not zero
    for comp in ("malignant", "non_malignant"):
        assert (const[f"cv_{comp}"] < 0.1).all()
        assert (const[f"cv_{comp}"] < others[f"cv_{comp}"].min()).all()
        assert (const[f"fano_{comp}"] < others[f"fano_{comp}"].min()).all()
        assert (const[f"mad_{comp}"] < others[f"mad_{comp}"].min()).all()
        assert const[f"cv_{comp}"].nunique() == 1  # identical genes give identical values
        assert (const[f"detection_rate_{comp}"] == 1.0).all()


def test_poisson_and_overdispersed_genes(gene_summary, synthetic_adata):
    pois = gene_summary[gene_summary["truth"] == "poisson"]
    od = gene_summary[gene_summary["truth"] == "overdispersed"]
    od_strong = od[od["mu"] >= 2.0]
    assert len(pois) > 100 and len(od) == 40 and len(od_strong) >= 10
    # raw counts (same distribution in every group, so cells can be pooled): Poisson Fano ~ 1, NB with theta = 2 has Fano = 1 + mu/2
    raw = gene_noise_table(synthetic_adata.X, synthetic_adata.var_names).set_index("gene")
    assert abs(raw.loc[pois["gene"], "fano"].median() - 1.0) < 0.15
    assert raw.loc[od["gene"], "fano"].median() > 1.5
    assert (raw.loc[od_strong["gene"], "fano"] > 1.5).all()
    # CP10K scaling multiplies Fano by about 1e4 / library size, so the classes are compared with each other
    for comp in ("malignant", "non_malignant"):
        assert od[f"fano_{comp}"].median() > 1.5 * pois[f"fano_{comp}"].median()
        assert (od_strong[f"fano_{comp}"] > 1.5 * pois[f"fano_{comp}"].median()).all()
    # neither class is differentially noisy between the compartments
    assert abs(pois["log2fc_cv_median"].median()) < 0.1
    assert abs(od["log2fc_cv_median"].median()) < 0.1


def test_downsampled_metrics(app_data_dir, processed_dir):
    for ds in DATASETS:
        patient = pd.read_parquet(app_data_dir / "patient" / f"{ds}.parquet")
        assert np.isfinite(patient[["cv", "cv_ds", "cv2_resid", "cv2_resid_ds"]].to_numpy()).all()
        assert sps.spearmanr(patient["cv"], patient["cv_ds"]).statistic > 0.9
        summary = json.loads((processed_dir / ds / "summary.json").read_text())
        assert summary["downsampling"]["target_depth"] >= 200
        assert summary["downsampling"]["cells_above_target"] > 0
        assert summary["downsampling"]["spearman_cv_vs_cv_ds"] > 0.9
        full = pd.read_parquet(processed_dir / ds / "gene_metrics_patient.parquet")
        expressed = full[full["gene"].isin(patient["gene"])]
        assert np.isfinite(expressed[["mean_ds", "cv_ds", "mad_ds", "fano_ds", "cv2_resid_ds"]].to_numpy()).all()


def test_expressed_flag_logic(app_data_dir, fast_params):
    gs = pd.read_parquet(app_data_dir / "gene_summary.parquet")
    threshold = fast_params["metrics"]["min_detection_rate"]
    expected = (gs["detection_rate_malignant"] >= threshold) & (gs["detection_rate_non_malignant"] >= threshold)
    assert gs["expressed"].dtype == bool
    assert (gs["expressed"] == expected).all()
    assert 0 < gs["expressed"].sum() < len(gs)
    # paired tests (and the BH family) are restricted to expressed genes
    assert gs.loc[gs["expressed"], ["n_pairs", "p_value", "q_value"]].notna().all().all()
    assert gs.loc[~gs["expressed"], ["n_pairs", "p_value", "q_value", "p_value_mad", "p_value_fano"]].isna().all().all()
    datasets = pd.read_parquet(app_data_dir / "datasets.parquet").set_index("dataset")
    celltype = pd.read_parquet(app_data_dir / "celltype_summary.parquet")
    for ds in DATASETS:
        expressed = set(gs.loc[(gs["dataset"] == ds) & gs["expressed"], "gene"])
        assert datasets.loc[ds, "n_genes_expressed"] == len(expressed)
        patient = pd.read_parquet(app_data_dir / "patient" / f"{ds}.parquet")
        assert set(patient["gene"]) == expressed
        assert set(celltype.loc[celltype["dataset"] == ds, "gene"]) <= expressed


def test_sample_scores_one_row_per_patient_and_compartment(app_data_dir):
    ss = pd.read_parquet(app_data_dir / "sample_scores.parquet")
    cells = pd.read_parquet(app_data_dir / "cell_scores.parquet")
    for ds in DATASETS:
        sub = ss[ss["dataset"] == ds]
        assert len(sub) == 2 * N_PATIENTS[ds]
        assert not sub.duplicated(["patient", "compartment"]).any()
        assert set(sub["compartment"]) == {"malignant", "non_malignant"}
        n_cells = cells[cells["dataset"] == ds].groupby(["patient", "compartment"]).size()
        indexed = sub.set_index(["patient", "compartment"])
        assert (indexed["n_cells"] == n_cells.reindex(indexed.index)).all()
    non_malignant = ss["compartment"] == "non_malignant"
    assert ss.loc[non_malignant, "depth_pseudobulk"].isna().all()
    assert np.isfinite(ss.loc[~non_malignant, "depth_pseudobulk"]).all()
    for col in ("depth2_pseudobulk", "depth_cell_median", "depth2_cell_median"):
        assert np.isfinite(ss[col]).all() and (ss[col] > 0).all(), col
    medians = cells.groupby(["dataset", "patient", "compartment"])[["depth", "depth2"]].median()
    joined = ss.set_index(["dataset", "patient", "compartment"]).join(medians)
    np.testing.assert_allclose(joined["depth_cell_median"], joined["depth"], rtol=1e-6)
    np.testing.assert_allclose(joined["depth2_cell_median"], joined["depth2"], rtol=1e-6)


def test_cell_scores_match_depth2_recomputed_from_the_fixture(processed_dir, app_data_dir, synthetic_adata, fast_params):
    cells = pd.read_parquet(processed_dir / FIRST / "cell_scores.parquet")
    assert np.isfinite(cells[["depth", "depth2"]]).all().all()
    assert (cells[["depth", "depth2"]] > 0).all().all()
    adata = _qc_like_pipeline(synthetic_adata, fast_params)
    assert cells["cell_name"].tolist() == adata.obs["cell_name"].tolist()
    np.testing.assert_array_equal(cells["n_counts"].to_numpy(), adata.obs["n_counts"].to_numpy())
    used = _depth_gene_mask(adata, fast_params)
    summary = json.loads((processed_dir / FIRST / "summary.json").read_text())
    assert summary["depth"]["n_genes_used"] == int(used.sum())
    log_used = sp.csr_matrix(log2p1(normalize_cp10k(sp.csr_matrix(adata.X), fast_params["normalisation"]["target_sum"]))[:, used])
    mean, sd = gene_mean_sd(log_used)
    expected = depth2_per_cell(log_used, mean, sd)
    np.testing.assert_allclose(cells["depth2"].to_numpy(), expected, rtol=0, atol=1e-8)
    app_cells = pd.read_parquet(app_data_dir / "cell_scores.parquet")
    np.testing.assert_allclose(app_cells.loc[app_cells["dataset"] == FIRST, "depth2"].to_numpy(), expected, rtol=0, atol=1e-8)


def test_pseudobulk_scores_match_recomputation_from_counts(processed_dir, synthetic_adata, fast_params):
    adata = _qc_like_pipeline(synthetic_adata, fast_params)
    used = _depth_gene_mask(adata, fast_params)
    labels = adata.obs["patient"].astype(str) + "|" + adata.obs["compartment"].astype(str)
    table, n_cells = pseudobulk_log_cpm(sp.csr_matrix(adata.X)[:, used], labels, min_cells=fast_params["groups"]["min_cells_per_group"])
    expected_depth2 = depth2_matrix(table)
    expected_depth = depth_matrix(table, pd.Series(table.columns.str.endswith("|non_malignant"), index=table.columns))
    ss = pd.read_parquet(processed_dir / FIRST / "sample_scores.parquet")
    ss.index = ss["patient"] + "|" + ss["compartment"]
    assert (ss.loc[table.columns, "n_cells"].to_numpy() == n_cells.to_numpy()).all()
    np.testing.assert_allclose(ss.loc[table.columns, "depth2_pseudobulk"].to_numpy(), expected_depth2.to_numpy(), rtol=0, atol=1e-8)
    np.testing.assert_allclose(ss.loc[expected_depth.index, "depth_pseudobulk"].to_numpy(), expected_depth.to_numpy(), rtol=0, atol=1e-8)


def test_celltype_summary(app_data_dir):
    ct = pd.read_parquet(app_data_dir / "celltype_summary.parquet")
    assert set(ct["cell_type"]) == {"Malignant", "T_cell", "Fibroblast", "Epithelial"}
    assert not ct.duplicated(["dataset", "cell_type", "gene"]).any()
    assert (ct["n_patients"] >= 1).all()
    assert ct.loc[ct["dataset"] == FIRST, "n_patients"].max() == N_PATIENTS[FIRST]
    assert np.isfinite(ct[["mean", "detection_rate", "mad"]].to_numpy()).all()
    assert np.isfinite(ct.loc[ct["mean"] > 0, ["cv", "fano"]].to_numpy()).all()


# ---------------------------------------------------------------------------
# process_dataset on its own and the command line
# ---------------------------------------------------------------------------
def _small_params(**overrides) -> dict:
    """Default parameters shrunk for a few dozen cells per group; every sample stays below the Scrublet minimum, so no doublet detection runs."""
    params = copy.deepcopy(load_config()["params"])
    params["qc"]["scrublet_min_cells_per_sample"] = 10**9
    params["qc"]["min_cells_per_gene"] = 3
    params["groups"].update(min_cells_per_group=10, min_patients_per_test=3, min_cells_per_celltype_group=10)
    params["downsampling"]["min_depth"] = 100
    for section, values in overrides.items():
        params[section].update(values)
    return params


def _write_config(path: Path, params: dict) -> Path:
    cfg = load_config()
    path.write_text(yaml.safe_dump({"sources": cfg["sources"], "datasets": [], "params": params}))
    return path


def _run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "sctcna", *args], cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=300)


def test_cli_run_and_aggregate_synthetic(tmp_path):
    config = _write_config(tmp_path / "config.yaml", _small_params())
    out = tmp_path / "processed"
    res = _run_cli("run", "--synthetic", "--synthetic-patients", "3", "--synthetic-cells", "40", "--out", str(out), "--config", str(config))
    assert res.returncode == 0, res.stderr[-2000:]
    summary = json.loads((out / "SYN_Test" / "summary.json").read_text())
    assert summary["dataset"] == "SYN_Test" and summary["cancer_type"] == "SYN"
    assert summary["n_patients"] == summary["n_patients_paired"] == 3
    assert summary["qc"]["doublets_removed"] == 0
    assert summary["qc"]["scrublet_skipped_samples"] == ["P1", "P2", "P3"]
    assert summary["params"]["qc"]["scrublet_min_cells_per_sample"] == 10**9  # the config file was used
    for name in PROCESSED_COLUMNS:
        assert (out / "SYN_Test" / name).exists(), name
    app = tmp_path / "app"
    res = _run_cli("aggregate", "--config", str(config), "--raw", str(tmp_path / "raw"), "--processed", str(out), "--app", str(app), "--include-synthetic")
    assert res.returncode == 0, res.stderr[-2000:]
    assert json.loads(res.stdout.strip().splitlines()[-1]) == {"datasets": 1, "app_dir": str(app)}
    datasets = pd.read_parquet(app / "datasets.parquet")
    assert len(datasets) == 1 and datasets.loc[0, "dataset"] == "SYN_Test"
    _assert_columns(datasets, APP_COLUMNS["datasets.parquet"], "datasets.parquet")
    assert (app / "patient" / "SYN_Test.parquet").exists() and (app / "build.json").exists()


def test_process_dataset_with_scrublet_disabled(tmp_path):
    params = _small_params(qc={"run_scrublet": False})
    summary = process_dataset(make_synthetic_adata(n_patients=3, cells_per_group=30, seed=0), params, META, tmp_path)
    assert summary["qc"]["doublets_removed"] == 0 and summary["qc"]["scrublet_skipped_samples"] == []
    assert (tmp_path / "SYN_Test" / "summary.json").exists()
    cells = pd.read_parquet(tmp_path / "SYN_Test" / "cell_scores.parquet")
    assert cells["doublet_score"].isna().all()
    assert len(cells) == summary["qc"]["cells_final"]


def test_process_dataset_without_celltype_groups_writes_an_empty_celltype_table(tmp_path):
    params = _small_params(groups={"min_cells_per_celltype_group": 10**6})
    process_dataset(make_synthetic_adata(n_patients=3, cells_per_group=30, seed=0), params, META, tmp_path)
    table = pd.read_parquet(tmp_path / "SYN_Test" / "gene_metrics_celltype.parquet")
    assert table.empty
    _assert_columns(table, PROCESSED_COLUMNS["gene_metrics_celltype.parquet"], "gene_metrics_celltype.parquet")


def test_noise_contrast_tables(processed_dir, app_data_dir):
    per = pd.read_parquet(processed_dir / "SYN_Test" / "noise_contrasts_patient.parquet")
    summ = pd.read_parquet(processed_dir / "SYN_Test" / "noise_contrasts.parquet")
    assert {"non_malignant", "T_cell", "Fibroblast", "Epithelial"} <= set(summ["reference"])
    allref = per[per["reference"] == "non_malignant"]
    assert allref["patient"].nunique() == 6 and (allref["n_genes"] > 100).all()
    # simulated compartments share means and depth, so the median contrast over all genes is close to zero
    assert allref["excess_cv2_ds"].abs().max() < 0.1
    row = summ.set_index("reference").loc["non_malignant"]
    assert row["n_patients"] == 6 and 0 <= row["p_excess_cv2_ds"] <= 1
    app = pd.read_parquet(app_data_dir / "noise_contrasts.parquet")
    assert set(app["dataset"]) == {"SYN_Test", "SYN_Two"}
    assert "excess_cv2_ds_vs_non_malignant" in pd.read_parquet(app_data_dir / "datasets.parquet").columns
