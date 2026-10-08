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

CHUNK_ROWS = 10_000  # cells scored at a time; bounds temporary memory on large datasets


def _row_sums(values: np.ndarray, indptr: np.ndarray) -> np.ndarray:
    """Per-row sums of CSR-ordered values (sequential per row, so independent of how rows are chunked)."""
    rows = np.repeat(np.arange(len(indptr) - 1), np.diff(indptr))
    return np.bincount(rows, weights=values, minlength=len(indptr) - 1)


def _sd_from_sums(s1: np.ndarray, s2: np.ndarray, g: int) -> np.ndarray:
    if g < 2:
        return np.full(s1.shape, np.nan)
    var = (s2 - s1**2 / g) / (g - 1)
    return np.sqrt(np.where(var < 0, 0.0, var))


def gene_mean_sd(X: sp.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    """Per-gene mean and sample SD (ddof = 1) across the rows of ``X``."""
    X = sp.csr_matrix(X)
    n, g = X.shape
    s1, s2 = np.zeros(g), np.zeros(g)
    for a in range(0, n, CHUNK_ROWS):
        lo, hi = X.indptr[a], X.indptr[min(n, a + CHUNK_ROWS)]
        d, idx = X.data[lo:hi].astype(np.float64), X.indices[lo:hi]
        s1 += np.bincount(idx, weights=d, minlength=g)
        s2 += np.bincount(idx, weights=d * d, minlength=g)
    mean = s1 / n
    var = (s2 - n * mean**2) / (n - 1) if n > 1 else np.full(g, np.nan)
    return mean, np.sqrt(np.where(var <= 1e-12 * mean**2, 0.0, var))


def depth_scores_from_counts(
    counts: sp.spmatrix,
    gene_mask: np.ndarray,
    groups: np.ndarray,
    scored_rows: np.ndarray,
    reference_rows: np.ndarray,
    target_sum: float = 1e4,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell DEPTH and DEPTH2 on log2(CP10K + 1), one group of cells at a time.

    * DEPTH is computed for ``scored_rows`` only (the malignant cells, the
      "tumour samples" of DEPTH) against the mean over ``reference_rows`` (the
      non-malignant cells, DEPTH's "normal" reference; the scored rows' own mean
      when there are none). Other rows get NaN.
    * DEPTH2 z-scores every gene across the cells of the cell's own group
      (patient x compartment), so a cell's score measures how far its profile
      departs from the other cells of the same tumour (or the same patient's
      non-malignant cells). Groups with fewer than two cells get NaN.

    scTCNA assumptions: the reference set for DEPTH and the z-scoring cohort for
    DEPTH2. Scoring both compartments against one pooled reference would turn
    differences in mean expression between compartments into apparent
    heterogeneity differences.
    """
    from .normalize import log2p1, normalize_cp10k, select_columns

    counts = sp.csr_matrix(counts)
    n = counts.shape[0]
    cols = np.flatnonzero(gene_mask)
    groups = np.asarray(groups)
    scored_rows = np.asarray(scored_rows, dtype=bool)
    reference_rows = np.asarray(reference_rows, dtype=bool)
    ref = reference_rows if reference_rows.any() else scored_rows

    def log_block(rows: np.ndarray) -> sp.csr_matrix:
        return select_columns(log2p1(normalize_cp10k(counts[rows], target_sum)), cols)

    # reference mean, accumulated block by block
    r1 = np.zeros(cols.size)
    ref_idx = np.flatnonzero(ref)
    for a in range(0, ref_idx.size, CHUNK_ROWS):
        L = log_block(ref_idx[a : a + CHUNK_ROWS])
        r1 += np.bincount(L.indices, weights=L.data.astype(np.float64), minlength=cols.size)
    ref_mean = r1 / max(ref_idx.size, 1)

    depth = np.full(n, np.nan)
    depth2 = np.full(n, np.nan)
    for _, idx in pd.Series(np.arange(n)).groupby(groups, sort=False).groups.items():
        idx = np.asarray(idx)
        L = log_block(idx)
        if idx.size >= 2:
            mean, sd = gene_mean_sd(L)
            depth2[idx] = depth2_per_cell(L, mean, sd)
        scored = scored_rows[idx]
        if scored.any():
            depth[idx[scored]] = depth_per_cell(L[scored], ref_mean)
    return depth, depth2


def depth2_per_cell(X: sp.spmatrix, mean: np.ndarray | None = None, sd: np.ndarray | None = None) -> np.ndarray:
    """DEPTH2 score of every row (cell) of a cells x genes log-expression matrix.

    Genes are z-scored across the rows of ``X`` unless ``mean``/``sd`` are given
    (e.g. computed on a larger reference set).
    """
    X = sp.csr_matrix(X, dtype=np.float64, copy=True)
    X.eliminate_zeros()
    n, g = X.shape
    if mean is None or sd is None:
        if n < 2:  # z-scores need at least two samples; the score is undefined
            return np.full(n, np.nan)
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
    X = sp.csr_matrix(X, dtype=np.float64, copy=True)
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
    if M.shape[1] < 2:
        return pd.Series(np.nan, index=expr.columns, name="DEPTH2")
    mean = M.mean(axis=1, keepdims=True)
    sd = M.std(axis=1, ddof=1, keepdims=True)
    sd = np.where(sd <= 1e-12 * np.abs(mean), 0.0, sd)  # constant genes: rounding must not create |z| ~ 1
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


def pseudobulk_log_cpm(counts: sp.spmatrix, groups: pd.Series, min_cells: int = 1, columns=None) -> tuple[pd.DataFrame, pd.Series]:
    """Sum raw counts per group, scale to counts per million and log2(x + 1).

    ``columns`` (boolean mask or indices) restricts the genes after summing;
    CPM totals are taken over the kept genes only.

    Returns (genes x groups table, cells per group). scTCNA assumption: the
    bulk-like sample used for sample-level DEPTH/DEPTH2.
    """
    counts = sp.csr_matrix(counts)
    groups = pd.Series(np.asarray(groups), index=range(counts.shape[0]))
    codes, uniques = pd.factorize(groups, sort=True)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes)))
    sums = np.asarray((ind @ counts).todense(), dtype=np.float64)  # groups x genes; float64 indicator, no copy of counts
    if columns is not None:
        sums = sums[:, np.asarray(columns)]
    n_cells = pd.Series(np.bincount(codes, minlength=len(uniques)), index=uniques)
    keep = n_cells.to_numpy() >= min_cells
    sums = sums[keep]
    totals = sums.sum(axis=1, keepdims=True)
    cpm = np.where(totals > 0, sums / np.where(totals > 0, totals, 1.0) * 1e6, 0.0)
    table = pd.DataFrame(np.log2(cpm + 1.0).T, columns=np.asarray(uniques)[keep])
    return table, n_cells[keep]
