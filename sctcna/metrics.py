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
* Fano           : variance / mean (Poisson noise gives Fano = 1)
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


def _as_csr(X) -> sp.csr_matrix:
    X = sp.csr_matrix(X, dtype=np.float64)
    X.sum_duplicates()
    X.eliminate_zeros()
    return X


def gene_moments(X) -> dict[str, np.ndarray]:
    """Per-gene n, mean, sample variance (ddof=1), SD and detection rate."""
    X = _as_csr(X)
    n, g = X.shape
    s1 = np.asarray(X.sum(axis=0)).ravel()
    s2 = np.asarray(X.multiply(X).sum(axis=0)).ravel()
    nnz = np.bincount(X.indices, minlength=g).astype(np.float64)
    mean = s1 / n if n else np.full(g, np.nan)
    if n > 1:
        var = (s2 - n * mean**2) / (n - 1)
        var = np.where(var < 0, 0.0, var)  # numerical noise on constant genes
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
        mean = np.asarray(X.sum(axis=0)).ravel() / n
    dev = np.abs(X.data - mean[X.indices])
    sum_dev = np.bincount(X.indices, weights=dev, minlength=g)
    nnz = np.bincount(X.indices, minlength=g)
    return (sum_dev + (n - nnz) * np.abs(mean)) / n


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
    a = np.asarray(numerator, dtype=float)
    b = np.asarray(denominator, dtype=float)
    out = np.full(np.broadcast(a, b).shape, np.nan)
    ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    out[ok] = np.log2(a[ok] / b[ok])
    return out
