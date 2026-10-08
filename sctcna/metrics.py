"""Per-gene noise metrics.

Definitions follow the TCNA manuscript (Shah et al. 2026) and the TCNA code:

* mean           : arithmetic mean across cells
* SD             : sample standard deviation (ddof = 1)
* CV             : SD / mean  (TCNA's web code reports CV x 100; the ratio is stored here)
* CV^2           : CV squared
* MAD            : mean absolute deviation around the mean, mean(|x - mean|)
* Differential noise (TCNA Eq. 1): logFC = log2(CV_tumour / CV_normal)

Single-cell additions (scTCNA assumptions, see METHODS.md):

* detection_rate : fraction of cells with a non-zero value
* Fano           : variance / mean on the CP10K scale (for Poisson counts at library
                   size N this is about 10^4 / N, so it depends on sequencing depth)
* cv2_resid      : residual of log10(CV^2) from a LOWESS trend on log10(mean),
                   fitted across genes within the same group of cells, i.e.
                   noise corrected for its dependence on mean expression.

All functions accept a cells x genes matrix (sparse CSR or dense) and return
NaN where a metric is undefined (mean = 0), rather than 0.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from statsmodels.nonparametric.smoothers_lowess import lowess

METRIC_COLUMNS = ["n_cells", "mean", "detection_rate", "sd", "var", "cv", "cv2", "mad", "fano"]


CHUNK_ROWS = 20_000  # rows processed at a time; bounds temporary memory on large groups of cells


def _as_csr(X) -> sp.csr_matrix:
    """CSR view of ``X`` (dense or sparse) with duplicates summed; the caller's matrix is never modified."""
    X = sp.csr_matrix(X)
    if not X.has_canonical_format:
        X = X.copy()
        X.sum_duplicates()
    return X


def _chunks(X: sp.csr_matrix):
    """Yield (data as float64, column indices) of consecutive row blocks without copying the matrix."""
    for a in range(0, X.shape[0], CHUNK_ROWS):
        b = min(X.shape[0], a + CHUNK_ROWS)
        lo, hi = X.indptr[a], X.indptr[b]
        yield X.data[lo:hi].astype(np.float64), X.indices[lo:hi]


def gene_moments(X) -> dict[str, np.ndarray]:
    """Per-gene n, mean, sample variance (ddof=1), SD and detection rate.

    Sums are accumulated per row block, so memory stays bounded for large
    groups. Explicitly stored zeros do not count as detected.
    """
    X = _as_csr(X)
    n, g = X.shape
    s1, s2, nnz = np.zeros(g), np.zeros(g), np.zeros(g)
    for d, idx in _chunks(X):
        s1 += np.bincount(idx, weights=d, minlength=g)
        s2 += np.bincount(idx, weights=d * d, minlength=g)
        nnz += np.bincount(idx, weights=(d != 0).astype(np.float64), minlength=g)
    mean = s1 / n if n else np.full(g, np.nan)
    if n > 1:
        var = (s2 - n * mean**2) / (n - 1)
        # one-pass variance: rounding can leave a tiny positive or negative value for constant genes
        var = np.where(var <= 1e-12 * mean**2, 0.0, var)
    else:
        var = np.full(g, np.nan)
    return {"n_cells": np.full(g, n, dtype=np.int64), "mean": mean, "var": var, "sd": np.sqrt(var), "detection_rate": nnz / n if n else np.full(g, np.nan), "_nnz": nnz}


def mean_absolute_deviation(X, mean: np.ndarray | None = None) -> np.ndarray:
    """MAD = mean(|x - mean|) per gene (TCNA definition: mean absolute deviation about the mean)."""
    X = _as_csr(X)
    n, g = X.shape
    if n == 0:
        return np.full(g, np.nan)
    if mean is None:
        mean = gene_moments(X)["mean"]
    sum_dev, stored = np.zeros(g), np.zeros(g)
    for d, idx in _chunks(X):
        sum_dev += np.bincount(idx, weights=np.abs(d - mean[idx]), minlength=g)
        stored += np.bincount(idx, minlength=g)
    # every entry not stored is a zero, contributing |0 - mean|
    return (sum_dev + (n - stored) * np.abs(mean)) / n


def gene_noise_table(X, gene_names=None) -> pd.DataFrame:
    """All per-gene metrics for one group of cells."""
    X = _as_csr(X)
    mom = gene_moments(X)
    mean, var, sd = mom["mean"], mom["var"], mom["sd"]
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(mean > 0, sd / mean, np.nan)
        fano = np.where(mean > 0, var / mean, np.nan)
    mad = mean_absolute_deviation(X, mean)
    df = pd.DataFrame(
        {
            "n_cells": mom["n_cells"],
            "mean": mean,
            "detection_rate": mom["detection_rate"],
            "sd": sd,
            "var": var,
            "cv": cv,
            "cv2": cv**2,
            "mad": mad,
            "fano": fano,
        }
    )
    if gene_names is not None:
        df.insert(0, "gene", np.asarray(gene_names))
    return df


def mean_trend_residual(mean: np.ndarray, cv2: np.ndarray, frac: float = 0.3, min_points: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Residual of log10(CV^2) from a LOWESS fit on log10(mean) across genes.

    Returns (residual, fitted) on the log10 scale; NaN where undefined or where
    the fit could not be made. scTCNA assumption: corrects the strong
    mean-variance dependence of single-cell noise (Brennecke et al. 2013 style).
    """
    mean = np.asarray(mean, dtype=float)
    cv2 = np.asarray(cv2, dtype=float)
    ok = np.isfinite(mean) & np.isfinite(cv2) & (mean > 0) & (cv2 > 0)
    resid = np.full(mean.shape, np.nan)
    fitted = np.full(mean.shape, np.nan)
    if ok.sum() < min_points:
        return resid, fitted
    x = np.log10(mean[ok])
    y = np.log10(cv2[ok])
    delta = 0.01 * (x.max() - x.min())
    fit = lowess(y, x, frac=frac, it=2, delta=delta, return_sorted=False)
    fitted[ok] = fit
    resid[ok] = y - fit
    return resid, fitted


def log2_fold_change(numerator, denominator) -> np.ndarray:
    """TCNA Eq. 1: log2(noise_tumour / noise_normal); NaN where either side is not positive."""
    a, b = np.broadcast_arrays(np.asarray(numerator, dtype=float), np.asarray(denominator, dtype=float))
    out = np.full(a.shape, np.nan)
    ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    out[ok] = np.log2(a[ok] / b[ok])
    return out


def trend_excess(ref_mean, ref_cv2, query_mean, query_cv2, frac: float = 0.3, min_points: int = 50) -> np.ndarray:
    """Excess noise of query genes over a reference group's mean-CV^2 trend.

    A LOWESS curve of log10(CV^2) on log10(mean) is fitted across the genes of
    the reference group; each query gene's residual is log10(CV^2_query) minus
    the curve evaluated at the query gene's mean. Query means outside the
    reference's fitted range give NaN (no extrapolation). Positive values mean
    the query cells are noisier than reference cells at the same expression
    level. scTCNA addition, used to compare malignant cells with non-malignant
    cells (or one lineage) without the mean dependence of CV.
    """
    rm, rc = np.asarray(ref_mean, dtype=float), np.asarray(ref_cv2, dtype=float)
    qm, qc = np.asarray(query_mean, dtype=float), np.asarray(query_cv2, dtype=float)
    out = np.full(qm.shape, np.nan)
    ok = np.isfinite(rm) & np.isfinite(rc) & (rm > 0) & (rc > 0)
    if ok.sum() < min_points:
        return out
    x, y = np.log10(rm[ok]), np.log10(rc[ok])
    fit = lowess(y, x, frac=frac, it=2, delta=0.01 * (x.max() - x.min()), return_sorted=True)
    xs, ys = fit[:, 0], fit[:, 1]
    q_ok = np.isfinite(qm) & np.isfinite(qc) & (qm > 0) & (qc > 0)
    lx = np.log10(np.where(q_ok, qm, 1.0))
    inside = q_ok & (lx >= xs[0]) & (lx <= xs[-1])
    out[inside] = np.log10(qc[inside]) - np.interp(lx[inside], xs, ys)
    return out
