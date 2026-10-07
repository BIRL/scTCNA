"""Quality control: cell/gene filters, mitochondrial fraction, doublets, depth downsampling.

All thresholds are scTCNA assumptions (the TCNA manuscript concerns bulk data);
defaults live in configs/datasets.yaml and are documented in METHODS.md. 3CA
and TISCH2 both apply their own QC upstream, so these filters are a light,
uniform re-filter.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp

from .normalize import select_columns

log = logging.getLogger(__name__)


def add_qc_metrics(adata: ad.AnnData, mito_prefix: str = "MT-") -> ad.AnnData:
    X = sp.csr_matrix(adata.X)
    adata.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.diff(X.indptr)
    mito = adata.var_names.str.upper().str.startswith(mito_prefix.upper())
    adata.var["mito"] = mito
    mito_counts = np.asarray(select_columns(X, mito).sum(axis=1)).ravel() if mito.any() else np.zeros(adata.n_obs)
    with np.errstate(divide="ignore", invalid="ignore"):
        adata.obs["pct_mito"] = np.where(adata.obs["n_counts"] > 0, 100.0 * mito_counts / adata.obs["n_counts"], 0.0)
    return adata


def filter_cells(adata: ad.AnnData, min_genes: int = 200, min_counts: int = 500, max_pct_mito: float = 20.0) -> tuple[ad.AnnData, dict]:
    if "pct_mito" not in adata.obs:
        add_qc_metrics(adata)
    keep = (adata.obs["n_genes"] >= min_genes) & (adata.obs["n_counts"] >= min_counts) & (adata.obs["pct_mito"] <= max_pct_mito)
    stats = {
        "cells_before": int(adata.n_obs),
        "removed_low_genes": int((adata.obs["n_genes"] < min_genes).sum()),
        "removed_low_counts": int((adata.obs["n_counts"] < min_counts).sum()),
        "removed_high_mito": int((adata.obs["pct_mito"] > max_pct_mito).sum()),
        "cells_after": int(keep.sum()),
    }
    return adata[keep.to_numpy()].copy(), stats


def filter_genes(adata: ad.AnnData, min_cells: int = 10) -> tuple[ad.AnnData, dict]:
    X = sp.csr_matrix(adata.X)
    n_cells = np.bincount(X.indices, minlength=adata.n_vars)
    keep = n_cells >= min_cells
    stats = {"genes_before": int(adata.n_vars), "genes_after": int(keep.sum())}
    return adata[:, keep].copy(), stats


def run_scrublet(adata: ad.AnnData, batch_key: str = "sample", min_cells: int = 100, seed: int = 0) -> tuple[ad.AnnData, dict]:
    """Scrublet doublet detection per sample; predicted doublets are removed.

    Samples with fewer than ``min_cells`` cells are skipped (kNN simulation is
    unreliable there) and their cells are kept.
    """
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    skipped, failed = [], []
    for batch, idx in adata.obs.groupby(batch_key, observed=True).indices.items():
        if len(idx) < min_cells:
            skipped.append(str(batch))
            continue
        sub = adata[idx].copy()
        try:
            sc.pp.scrublet(sub, random_state=seed, verbose=False)
            adata.obs.iloc[idx, adata.obs.columns.get_loc("doublet_score")] = sub.obs["doublet_score"].to_numpy()
            adata.obs.iloc[idx, adata.obs.columns.get_loc("predicted_doublet")] = sub.obs["predicted_doublet"].to_numpy().astype(bool)
        except Exception as exc:  # pragma: no cover - depends on data
            failed.append(str(batch))
            log.warning("scrublet failed for %s=%s (%s); cells kept", batch_key, batch, exc)
    n_doublets = int(adata.obs["predicted_doublet"].sum())
    keep = ~adata.obs["predicted_doublet"].to_numpy(dtype=bool)
    stats = {"doublets_removed": n_doublets, "scrublet_skipped_samples": skipped, "scrublet_failed_samples": failed}
    return adata[keep].copy(), stats


def common_depth_target(n_counts: np.ndarray, quantile: float = 0.1, min_depth: int = 1000) -> int:
    return int(max(min_depth, np.quantile(n_counts, quantile)))


def downsample_counts(counts: sp.spmatrix, target: int, seed: int = 0) -> sp.csr_matrix:
    """Downsample every cell with more than ``target`` counts to exactly ``target`` (cells below are unchanged)."""
    tmp = ad.AnnData(X=sp.csr_matrix(counts, dtype=np.float32))
    sc.pp.downsample_counts(tmp, counts_per_cell=target, random_state=seed, replace=False)
    out = sp.csr_matrix(tmp.X, dtype=np.int32)
    out.eliminate_zeros()
    return out
