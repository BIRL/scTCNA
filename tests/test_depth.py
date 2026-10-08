"""Known-answer tests for DEPTH and DEPTH2 (sctcna.depth): the sparse per-cell code against the dense R replicas."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sctcna.depth import depth2_matrix, depth2_per_cell, depth_matrix, depth_per_cell, gene_mean_sd, pseudobulk_log_cpm

# Worked example from the DEPTH2 package README (10 genes x 5 samples) and the scores it reports.
README_EXPR = pd.DataFrame(
    [
        [8.23, 8.68, 8.19, 9.06, 8.70],
        [9.82, 8.50, 8.82, 8.78, 9.73],
        [11.78, 12.42, 11.98, 11.42, 11.55],
        [10.52, 9.68, 10.03, 9.21, 9.99],
        [6.68, 7.97, 7.69, 8.52, 8.09],
        [9.89, 7.53, 8.25, 9.33, 7.57],
        [10.30, 10.05, 9.81, 9.09, 10.14],
        [4.85, 10.08, 6.55, 8.82, 5.50],
        [10.63, 10.53, 9.88, 9.29, 9.19],
        [6.16, 3.41, 5.24, 4.96, 4.20],
    ],
    columns=["S1", "S2", "S3", "S4", "S5"],
)
README_SCORES = [0.40, 0.47, 0.28, 0.44, 0.28]


@pytest.fixture
def log_expression() -> np.ndarray:
    """40 cells x 25 genes of log2(count + 1) with a constant non-zero gene, a never-detected gene and many zeros."""
    rng = np.random.default_rng(7)
    X = np.log2(rng.poisson(1.5, size=(40, 25)) + 1.0)
    X[:, 3] = 2.0
    X[:, 7] = 0.0
    return X


def _columns(n: int) -> list[str]:
    return [f"c{i}" for i in range(n)]


def test_depth2_reproduces_readme_example():
    scores = depth2_matrix(README_EXPR)
    assert scores.name == "DEPTH2"
    assert list(scores.index) == list(README_EXPR.columns)
    np.testing.assert_allclose(np.round(scores.to_numpy(), 2), README_SCORES, rtol=0, atol=1e-9)
    per_cell = depth2_per_cell(sp.csr_matrix(README_EXPR.to_numpy().T))
    np.testing.assert_allclose(per_cell, scores.to_numpy(), rtol=0, atol=1e-10)


def test_depth2_per_cell_matches_dense_replica(log_expression):
    X = log_expression
    dense = depth2_matrix(pd.DataFrame(X.T, columns=_columns(X.shape[0])))
    np.testing.assert_allclose(depth2_per_cell(sp.csr_matrix(X)), dense.to_numpy(), rtol=0, atol=1e-10)


def test_depth_per_cell_matches_dense_replica_for_tumour_samples(log_expression):
    X = log_expression
    cols = _columns(X.shape[0])
    normal = np.zeros(X.shape[0], dtype=bool)
    normal[:15] = True
    expr = pd.DataFrame(X.T, columns=cols)
    dense = depth_matrix(expr, pd.Series(normal, index=cols))
    assert dense.name == "DEPTH"
    assert list(dense.index) == [c for c, is_normal in zip(cols, normal, strict=True) if not is_normal]
    sparse = depth_per_cell(sp.csr_matrix(X), X[normal].mean(axis=0))
    np.testing.assert_allclose(sparse[~normal], dense.to_numpy(), rtol=0, atol=1e-10)
    # without normal samples the reference is the mean of the tumour samples themselves
    dense_all = depth_matrix(expr)
    assert list(dense_all.index) == cols
    np.testing.assert_allclose(depth_per_cell(sp.csr_matrix(X), X.mean(axis=0)), dense_all.to_numpy(), rtol=0, atol=1e-10)


def test_zero_variance_genes_count_as_z_zero(log_expression):
    X = log_expression
    mean = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    assert (sd == 0).sum() == 2
    # the R package: scale() -> abs() -> NaN to 0 -> sd across genes (ddof = 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.abs((X - mean) / sd)
    z[~np.isfinite(z)] = 0.0
    expected = z.std(axis=1, ddof=1)
    np.testing.assert_allclose(depth2_per_cell(sp.csr_matrix(X)), expected, rtol=0, atol=1e-10)
    np.testing.assert_allclose(depth2_matrix(pd.DataFrame(X.T)).to_numpy(), expected, rtol=0, atol=1e-10)
    # the constant genes still count towards the number of genes: dropping them changes the scores
    assert not np.allclose(depth2_per_cell(sp.csr_matrix(X[:, sd > 0])), expected)


def test_reference_moments_can_come_from_a_larger_set(log_expression):
    X = log_expression
    mean, sd = gene_mean_sd(sp.csr_matrix(X))
    np.testing.assert_allclose(mean, X.mean(axis=0), rtol=0, atol=1e-12)
    np.testing.assert_allclose(sd, X.std(axis=0, ddof=1), rtol=0, atol=1e-12)
    full = depth2_per_cell(sp.csr_matrix(X), mean, sd)
    np.testing.assert_allclose(full, depth2_per_cell(sp.csr_matrix(X)), rtol=0, atol=1e-12)
    np.testing.assert_allclose(depth2_per_cell(sp.csr_matrix(X[5:9]), mean, sd), full[5:9], rtol=0, atol=1e-10)
    np.testing.assert_allclose(depth2_per_cell(sp.csr_matrix(X[:1]), mean, sd), full[:1], rtol=0, atol=1e-10)


def test_single_gene_gives_nan():
    X = sp.csr_matrix(np.array([[1.0], [2.0], [4.0]]))
    assert np.isnan(depth2_per_cell(X)).all()
    assert np.isnan(depth_per_cell(X, np.array([2.0]))).all()


def test_single_cell_without_reference_gives_nan():
    scores = depth2_per_cell(sp.csr_matrix(np.array([[1.0, 2.0, 3.0]])))
    assert scores.shape == (1,)
    assert np.isnan(scores[0])


def test_depth_per_cell_checks_reference_length():
    with pytest.raises(ValueError):
        depth_per_cell(sp.csr_matrix(np.ones((2, 3))), np.ones(2))


def test_pseudobulk_log_cpm_sums_counts_per_group():
    counts = sp.csr_matrix(np.array([[1, 2, 0], [3, 0, 1], [0, 0, 5], [2, 2, 2]]))
    table, n_cells = pseudobulk_log_cpm(counts, pd.Series(["b", "a", "b", "a"]))
    assert list(table.columns) == ["a", "b"]
    assert list(n_cells.index) == ["a", "b"] and n_cells.tolist() == [2, 2]
    for group, sums in {"a": np.array([5, 2, 3]), "b": np.array([1, 2, 5])}.items():
        np.testing.assert_allclose(table[group].to_numpy(), np.log2(sums / sums.sum() * 1e6 + 1.0), rtol=0, atol=1e-12)
    # groups below min_cells are dropped from both outputs
    table, n_cells = pseudobulk_log_cpm(counts, pd.Series(["b", "a", "b", "c"]), min_cells=2)
    assert list(table.columns) == ["b"] and n_cells.tolist() == [2]
    # a group without counts stays at zero instead of dividing by zero
    table, _ = pseudobulk_log_cpm(sp.csr_matrix(np.array([[0, 0], [1, 3]])), pd.Series(["z", "y"]))
    assert (table["z"] == 0).all()
