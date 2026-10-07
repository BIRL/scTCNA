"""DEPTH and DEPTH2 tumour-heterogeneity scores.

Both follow the published R packages used by TCNA:

* DEPTH  (Li et al. 2020, WangX-Lab/DEPTH):  for each tumour sample the score is
  the SD across genes of (x_g - mean_g^normal)^2, with the tumour mean used when
  no normal samples exist.
* DEPTH2 (Song & Wang 2022, WangX-Lab/DEPTH2): each gene is z-scored across
  samples (R ``scale``: sample SD, NaN -> 0); the score of a sample is the SD
  across genes of |z|.

TCNA applies them to log2(x + 1) expression. Here "sample" is either a single
cell (per-cell scores, computed sparsely and exactly) or a pseudo-bulk
patient x compartment profile (dense, the bulk-like equivalent). Zero-variance
genes contribute z = 0 and still count towards the number of genes, as in R.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _row_sums(values: np.ndarray, indptr: np.ndarray) -> np.ndarray:
    cs = np.concatenate([[0.0], np.cumsum(values, dtype=np.float64)])
    return cs[indptr[1:]] - cs[indptr[:-1]]


def _sd_from_sums(s1: np.ndarray, s2: np.ndarray, g: int) -> np.ndarray:
    if g < 2:
        return np.full(s1.shape, np.nan)
    var = (s2 - s1**2 / g) / (g - 1)
    return np.sqrt(np.where(var < 0, 0.0, var))


def gene_mean_sd(X: sp.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    X = sp.csr_matrix(X, dtype=np.float64)
    n, g = X.shape
    s1 = np.asarray(X.sum(axis=0)).ravel()
    s2 = np.asarray(X.multiply(X).sum(axis=0)).ravel()
    mean = s1 / n
    var = (s2 - n * mean**2) / (n - 1) if n > 1 else np.full(g, np.nan)
    return mean, np.sqrt(np.where(var < 0, 0.0, var))


def depth2_per_cell(X: sp.spmatrix, mean: np.ndarray | None = None, sd: np.ndarray | None = None) -> np.ndarray:
    """DEPTH2 score of every row (cell) of a cells x genes log-expression matrix.

    Genes are z-scored across the rows of ``X`` unless ``mean``/``sd`` are given
    (e.g. computed on a larger reference set).
    """
    X = sp.csr_matrix(X, dtype=np.float64)
    X.eliminate_zeros()
    n, g = X.shape
    if mean is None or sd is None:
        mean, sd = gene_mean_sd(X)
    valid = np.isfinite(sd) & (sd > 0)
    sd_safe = np.where(valid, sd, 1.0)
    a = np.where(valid, np.abs(mean) / sd_safe, 0.0)  # |z| of an implicit zero entry
    idx = X.indices
    z = np.where(valid[idx], (X.data - mean[idx]) / sd_safe[idx], 0.0)
    t1 = np.abs(z) - a[idx]
    t2 = z**2 - a[idx] ** 2
    s1 = a.sum() + _row_sums(t1, X.indptr)
    s2 = (a**2).sum() + _row_sums(t2, X.indptr)
    return _sd_from_sums(s1, s2, g)


def depth_per_cell(X: sp.spmatrix, reference_mean: np.ndarray) -> np.ndarray:
    """DEPTH score of every row (cell): SD across genes of (x_g - reference_mean_g)^2."""
    X = sp.csr_matrix(X, dtype=np.float64)
    X.eliminate_zeros()
    n, g = X.shape
    m = np.asarray(reference_mean, dtype=np.float64)
    if m.shape[0] != g:
        raise ValueError("reference_mean length must equal the number of genes")
    idx = X.indices
    d2 = (X.data - m[idx]) ** 2
    t1 = d2 - m[idx] ** 2
    t2 = d2**2 - m[idx] ** 4
    s1 = (m**2).sum() + _row_sums(t1, X.indptr)
    s2 = (m**4).sum() + _row_sums(t2, X.indptr)
    return _sd_from_sums(s1, s2, g)


def depth2_matrix(expr: pd.DataFrame) -> pd.Series:
    """Reference (dense) DEPTH2 on a genes x samples table, replicating the R package."""
    M = expr.to_numpy(dtype=np.float64)
    mean = M.mean(axis=1, keepdims=True)
    sd = M.std(axis=1, ddof=1, keepdims=True) if M.shape[1] > 1 else np.full_like(mean, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.abs((M - mean) / sd)
    z[~np.isfinite(z)] = 0.0
    return pd.Series(z.std(axis=0, ddof=1), index=expr.columns, name="DEPTH2")


def depth_matrix(expr: pd.DataFrame, is_normal: pd.Series | None = None) -> pd.Series:
    """Reference (dense) DEPTH on a genes x samples table, replicating the R package.

    Scores are returned for the tumour (non-normal) samples; the reference mean
    is taken over normal samples when any exist, otherwise over tumour samples.
    """
    M = expr.to_numpy(dtype=np.float64)
    if is_normal is None:
        is_normal = pd.Series(False, index=expr.columns)
    normal = is_normal.reindex(expr.columns).fillna(False).to_numpy(dtype=bool)
    tumour = ~normal
    ref = M[:, normal].mean(axis=1, keepdims=True) if normal.any() else M[:, tumour].mean(axis=1, keepdims=True)
    score = (M[:, tumour] - ref) ** 2
    return pd.Series(score.std(axis=0, ddof=1), index=expr.columns[tumour], name="DEPTH")


def pseudobulk_log_cpm(counts: sp.spmatrix, groups: pd.Series, min_cells: int = 1) -> tuple[pd.DataFrame, pd.Series]:
    """Sum raw counts per group, scale to counts per million and log2(x + 1).

    Returns (genes x groups table, cells per group). scTCNA assumption: the
    bulk-like sample used for sample-level DEPTH/DEPTH2.
    """
    counts = sp.csr_matrix(counts, dtype=np.float64)
    groups = pd.Series(np.asarray(groups), index=range(counts.shape[0]))
    codes, uniques = pd.factorize(groups, sort=True)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes)))
    sums = np.asarray((ind @ counts).todense())  # groups x genes
    n_cells = pd.Series(np.bincount(codes, minlength=len(uniques)), index=uniques)
    keep = n_cells.to_numpy() >= min_cells
    sums = sums[keep]
    totals = sums.sum(axis=1, keepdims=True)
    cpm = np.where(totals > 0, sums / np.where(totals > 0, totals, 1.0) * 1e6, 0.0)
    table = pd.DataFrame(np.log2(cpm + 1.0).T, columns=np.asarray(uniques)[keep])
    return table, n_cells[keep]
