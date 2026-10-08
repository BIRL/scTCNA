"""Known-answer tests for the per-gene noise metrics (sctcna.metrics)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sctcna.metrics import METRIC_COLUMNS, gene_moments, gene_noise_table, log2_fold_change, mean_absolute_deviation, mean_trend_residual

# Three cells x two genes, worked by hand (sample SD, ddof = 1):
#   gene A = (1, 2, 6): mean 3, var 7, sd sqrt(7), cv sqrt(7)/3, mad (2 + 1 + 3)/3 = 2, fano 7/3, detected in 3/3 cells
#   gene B = (2, 0, 4): mean 2, var 4, sd 2,       cv 1,         mad (0 + 2 + 2)/3,     fano 2,   detected in 2/3 cells
HAND = np.array([[1.0, 2.0], [2.0, 0.0], [6.0, 4.0]])
HAND_EXPECTED = {
    "mean": [3.0, 2.0],
    "var": [7.0, 4.0],
    "sd": [np.sqrt(7.0), 2.0],
    "cv": [np.sqrt(7.0) / 3.0, 1.0],
    "cv2": [7.0 / 9.0, 1.0],
    "mad": [2.0, 4.0 / 3.0],
    "fano": [7.0 / 3.0, 2.0],
    "detection_rate": [1.0, 2.0 / 3.0],
}


def _assert_table(table: pd.DataFrame, expected: dict[str, list[float]], atol: float = 1e-12) -> None:
    for col, values in expected.items():
        np.testing.assert_allclose(table[col].to_numpy(), values, rtol=0, atol=atol, err_msg=col)


def _hand_csr_with_explicit_zeros(dtype) -> sp.csr_matrix:
    """HAND as CSR with every entry stored, including the zero of gene B."""
    X = sp.csr_matrix((HAND.ravel().astype(dtype), np.tile([0, 1], 3), np.array([0, 2, 4, 6])), shape=HAND.shape)
    assert X.nnz == HAND.size
    return X


def test_constant_gene_has_exactly_zero_noise():
    table = gene_noise_table(np.full((25, 1), 3.0))
    assert table.loc[0, "mean"] == 3.0
    assert table.loc[0, "detection_rate"] == 1.0
    for col in ("sd", "var", "cv", "cv2", "mad", "fano"):
        assert table.loc[0, col] == 0.0, col


def test_all_zero_gene_is_undefined_not_zero():
    table = gene_noise_table(np.zeros((10, 1)))
    assert table.loc[0, "mean"] == 0.0
    assert table.loc[0, "detection_rate"] == 0.0
    assert table.loc[0, "sd"] == 0.0
    assert table.loc[0, "mad"] == 0.0
    assert np.isnan(table.loc[0, "cv"])
    assert np.isnan(table.loc[0, "cv2"])
    assert np.isnan(table.loc[0, "fano"])


def test_hand_computed_dense_matrix():
    table = gene_noise_table(HAND, gene_names=["A", "B"])
    assert list(table.columns) == ["gene", *METRIC_COLUMNS]
    assert table["gene"].tolist() == ["A", "B"]
    assert (table["n_cells"] == 3).all()
    _assert_table(table, HAND_EXPECTED)
    np.testing.assert_allclose(mean_absolute_deviation(HAND), HAND_EXPECTED["mad"], rtol=0, atol=1e-12)


def test_sparse_input_with_explicit_zeros_matches_dense_formulas():
    _assert_table(gene_noise_table(_hand_csr_with_explicit_zeros(np.float64)), HAND_EXPECTED)
    np.testing.assert_allclose(mean_absolute_deviation(_hand_csr_with_explicit_zeros(np.float64)), HAND_EXPECTED["mad"], rtol=0, atol=1e-12)
    # integer counts, the usual input
    _assert_table(gene_noise_table(sp.csr_matrix(HAND.astype(np.int32))), HAND_EXPECTED)


def test_gene_noise_table_leaves_float32_input_intact():
    X = _hand_csr_with_explicit_zeros(np.float32)
    gene_noise_table(X)
    assert np.array_equal(X.toarray(), HAND)


def test_detection_rate_counts_cells_with_non_zero_values():
    X = np.array([[0, 1, 5], [0, 0, 2], [0, 3, 0], [0, 0, 1]], dtype=float)
    np.testing.assert_allclose(gene_noise_table(X)["detection_rate"], [0.0, 0.5, 0.75])
    np.testing.assert_allclose(gene_moments(X)["detection_rate"], [0.0, 0.5, 0.75])


def test_single_cell_has_undefined_variance():
    table = gene_noise_table(np.array([[1.0, 4.0]]))
    np.testing.assert_allclose(table["mean"], [1.0, 4.0])
    assert table["sd"].isna().all()
    assert table["cv"].isna().all()
    assert table["fano"].isna().all()
    np.testing.assert_allclose(table["mad"], [0.0, 0.0])


@pytest.fixture(scope="module")
def poisson_nb_table() -> tuple[pd.DataFrame, np.ndarray]:
    """100 Poisson genes and 30 negative-binomial genes (theta = 1) with matched means, 5000 cells."""
    rng = np.random.default_rng(20240607)
    n_cells = 5000
    mu = rng.uniform(2.0, 20.0, size=100)
    poisson = rng.poisson(mu, size=(n_cells, mu.size))
    # NB with theta = 1: Var = mu + mu^2, so Fano = 1 + mu >= 3 for these means
    nb = rng.poisson(rng.gamma(shape=1.0, scale=mu[:30], size=(n_cells, 30)))
    is_nb = np.r_[np.zeros(100, dtype=bool), np.ones(30, dtype=bool)]
    return gene_noise_table(sp.csr_matrix(np.hstack([poisson, nb]))), is_nb


def test_poisson_genes_have_fano_near_one(poisson_nb_table):
    table, is_nb = poisson_nb_table
    assert abs(table.loc[~is_nb, "fano"].median() - 1.0) < 0.1


def test_negative_binomial_genes_rank_above_poisson_genes(poisson_nb_table):
    table, is_nb = poisson_nb_table
    fano = table["fano"].to_numpy()
    assert (fano[is_nb] > 2.0).all()
    assert fano[is_nb].min() > fano[~is_nb].max()
    resid, fitted = mean_trend_residual(table["mean"].to_numpy(), table["cv2"].to_numpy())
    assert np.isfinite(resid).all() and np.isfinite(fitted).all()
    assert (resid[is_nb] > 0).all()
    assert resid[is_nb].min() > resid[~is_nb].max()


def test_mean_trend_residual_is_nan_where_undefined_and_zero_on_an_exact_trend():
    mean = np.geomspace(0.1, 100.0, 80)
    cv2 = 1.0 / mean  # exact Poisson relation: log10(cv2) = -log10(mean)
    mean = np.r_[mean, 0.0, 5.0]
    cv2 = np.r_[cv2, 1.0, 0.0]
    resid, fitted = mean_trend_residual(mean, cv2)
    assert resid.shape == fitted.shape == mean.shape
    assert np.isnan(resid[-2:]).all() and np.isnan(fitted[-2:]).all()
    assert np.abs(resid[:-2]).max() < 1e-8
    np.testing.assert_allclose(fitted[:-2], np.log10(cv2[:-2]), rtol=0, atol=1e-8)


def test_mean_trend_residual_needs_enough_genes():
    mean = np.geomspace(1.0, 10.0, 10)
    resid, fitted = mean_trend_residual(mean, 1.0 / mean, min_points=50)
    assert np.isnan(resid).all() and np.isnan(fitted).all()
    resid, _ = mean_trend_residual(mean, 1.0 / mean, min_points=5)
    assert np.isfinite(resid).all()


def test_log2_fold_change_known_values():
    assert log2_fold_change(2.0, 1.0) == 1.0
    out = log2_fold_change([2.0, 1.0, 0.0, 4.0, np.nan, 3.0], [1.0, 0.0, 1.0, 8.0, 1.0, np.nan])
    np.testing.assert_array_equal(np.isnan(out), [False, True, True, False, True, True])
    assert out[0] == 1.0
    assert out[3] == -1.0


def test_log2_fold_change_vectorised_shapes():
    a = np.array([[2.0, 4.0], [8.0, 1.0]])
    out = log2_fold_change(a, np.ones((2, 2)))
    assert out.shape == (2, 2)
    np.testing.assert_allclose(out, [[1.0, 2.0], [3.0, 0.0]])
    assert log2_fold_change(np.array(2.0), np.array(1.0)).shape == ()
    assert log2_fold_change(pd.Series([2.0, 4.0]), pd.Series([1.0, 1.0])).tolist() == [1.0, 2.0]
    assert log2_fold_change([], []).shape == (0,)


def test_log2_fold_change_broadcasts_scalar_against_vector():
    np.testing.assert_allclose(log2_fold_change([2.0, 4.0], 1.0), [1.0, 2.0])
