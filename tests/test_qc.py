"""Tests for the QC step (sctcna.qc): per-cell metrics, filters, depth downsampling and the Scrublet wrapper."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sctcna.qc import add_qc_metrics, common_depth_target, downsample_counts, filter_cells, filter_genes, run_scrublet
from sctcna.synthetic import make_synthetic_adata


def _adata(X, genes: list[str]) -> ad.AnnData:
    X = np.asarray(X)
    return ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(index=[f"c{i}" for i in range(X.shape[0])]), var=pd.DataFrame(index=genes))


@pytest.fixture
def small() -> ad.AnnData:
    # genes: MT-CO1, mt-nd1 (the prefix match ignores case), ACTB, GAPDH
    return _adata([[5, 0, 2, 1], [0, 0, 0, 0], [1, 1, 1, 1], [10, 10, 0, 0]], ["MT-CO1", "mt-nd1", "ACTB", "GAPDH"])


def test_add_qc_metrics(small):
    add_qc_metrics(small, "MT-")
    assert small.obs["n_counts"].tolist() == [8, 0, 4, 20]
    assert small.obs["n_genes"].tolist() == [3, 0, 4, 2]
    np.testing.assert_allclose(small.obs["pct_mito"], [62.5, 0.0, 50.0, 100.0])
    assert small.var["mito"].tolist() == [True, True, False, False]
    # no mitochondrial genes at all -> 0 % everywhere
    other = _adata([[1, 2], [3, 4]], ["A", "B"])
    add_qc_metrics(other)
    assert (other.obs["pct_mito"] == 0).all()
    assert not other.var["mito"].any()


def test_filter_cells_thresholds_are_inclusive(small):
    kept, stats = filter_cells(small, min_genes=2, min_counts=4, max_pct_mito=50.0)
    assert kept.obs_names.tolist() == ["c2"]
    assert stats == {"cells_before": 4, "removed_low_genes": 1, "removed_low_counts": 1, "removed_high_mito": 2, "cells_after": 1}
    kept, stats = filter_cells(small, min_genes=0, min_counts=0, max_pct_mito=100.0)
    assert kept.n_obs == 4 and stats["cells_after"] == 4
    assert "pct_mito" in kept.obs  # computed on demand


def test_filter_genes(small):
    kept, stats = filter_genes(small, min_cells=2)
    assert kept.var_names.tolist() == ["MT-CO1", "mt-nd1", "ACTB", "GAPDH"]
    assert stats == {"genes_before": 4, "genes_after": 4}
    kept, stats = filter_genes(small, min_cells=3)
    assert kept.var_names.tolist() == ["MT-CO1"]
    assert stats == {"genes_before": 4, "genes_after": 1}
    assert np.array_equal(kept.X.toarray().ravel(), [5, 0, 1, 10])


def test_common_depth_target_uses_floor():
    assert common_depth_target(np.array([100, 200, 300, 400, 5000]), 0.1, 1000) == 1000
    n_counts = np.arange(1, 10001)
    target = common_depth_target(n_counts, 0.1, 10)
    assert target == 1000 and isinstance(target, int)
    assert common_depth_target(n_counts, 0.5, 10) == 5000


@pytest.fixture
def poisson_counts() -> tuple[np.ndarray, int]:
    rng = np.random.default_rng(0)
    dense = rng.poisson(2.0, size=(60, 40))
    target = int(np.median(dense.sum(axis=1)))
    totals = dense.sum(axis=1)
    assert (totals > target).any() and (totals <= target).any()
    return dense, target


def test_downsample_counts_known_properties(poisson_counts):
    dense, target = poisson_counts
    totals = dense.sum(axis=1)
    out = downsample_counts(sp.csr_matrix(dense), target, seed=0)
    assert out.dtype == np.int32 and out.format == "csr"
    assert (out.data != 0).all()
    ds = out.toarray()
    above = totals > target
    assert (ds[above].sum(axis=1) == target).all()
    assert np.array_equal(ds[~above], dense[~above])
    assert (ds <= dense).all()
    assert np.array_equal(downsample_counts(sp.csr_matrix(dense), target, seed=0).toarray(), ds)
    assert np.array_equal(downsample_counts(dense, target, seed=0).toarray(), ds)
    assert not np.array_equal(downsample_counts(sp.csr_matrix(dense), target, seed=1).toarray(), ds)
    # a target above every cell leaves the matrix unchanged
    assert np.array_equal(downsample_counts(sp.csr_matrix(dense), int(totals.max()), seed=0).toarray(), dense)


def test_downsample_counts_leaves_input_intact(poisson_counts):
    dense, target = poisson_counts
    counts = sp.csr_matrix(dense, dtype=np.int32)
    downsample_counts(counts, target, seed=0)
    assert np.array_equal(counts.toarray(), dense)


def test_run_scrublet_skips_small_samples_and_keeps_all_cells():
    adata = make_synthetic_adata(n_patients=2, cells_per_group=30, seed=0)
    adata.X = sp.csr_matrix(adata.X, dtype=np.int32)
    out, stats = run_scrublet(adata, "sample", min_cells=10_000, seed=0)
    assert out.n_obs == adata.n_obs
    assert out.obs_names.tolist() == adata.obs_names.tolist()
    assert stats == {"doublets_removed": 0, "scrublet_skipped_samples": ["P1", "P2"], "scrublet_failed_samples": []}
    assert out.obs["doublet_score"].isna().all()
    assert not out.obs["predicted_doublet"].any()


def test_run_scrublet_smoke_on_one_sample():
    # one sample of 200 simulated cells: a structural check of the wrapper, not of Scrublet's calls
    adata = make_synthetic_adata(n_patients=1, cells_per_group=100, seed=0)
    adata.X = sp.csr_matrix(adata.X, dtype=np.int32)
    out, stats = run_scrublet(adata, "sample", min_cells=100, seed=0)
    assert stats["scrublet_skipped_samples"] == [] and stats["scrublet_failed_samples"] == []
    assert out.n_obs == adata.n_obs - stats["doublets_removed"]
    assert out.obs["doublet_score"].between(0.0, 1.0).all()
    assert out.obs["predicted_doublet"].dtype == bool
    assert not out.obs["predicted_doublet"].any()
