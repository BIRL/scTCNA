"""Normalisation.

TCNA computes gene noise on normalised bulk expression (TPM/FPKM, linear scale)
and DEPTH/DEPTH2 on log2(x + 1) values. The single-cell equivalents used here
are counts per 10,000 (CP10K, library-size normalisation; scTCNA assumption)
and log2(CP10K + 1).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def normalize_cp10k(counts: sp.spmatrix, target_sum: float = 1e4, chunk_rows: int = 20_000) -> sp.csr_matrix:
    """Scale each cell (row) to ``target_sum`` total counts. Zero-count cells stay zero.

    Rows are scaled in place on one float32 copy, block by block, so the peak
    memory is about one extra matrix rather than two.
    """
    X = sp.csr_matrix(counts, dtype=np.float32, copy=True)
    totals = np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
    scale = np.where(totals > 0, target_sum / np.where(totals > 0, totals, 1.0), 0.0).astype(np.float32)
    for a in range(0, X.shape[0], chunk_rows):
        b = min(X.shape[0], a + chunk_rows)
        lo, hi = X.indptr[a], X.indptr[b]
        X.data[lo:hi] *= np.repeat(scale[a:b], np.diff(X.indptr[a : b + 1]))
    return X


def log2p1(X: sp.spmatrix) -> sp.csr_matrix:
    """log2(x + 1) applied to the stored entries (zeros stay zero)."""
    out = sp.csr_matrix(X, dtype=np.float32, copy=True)
    out.data = np.log2(out.data + 1.0).astype(np.float32)
    return out


def select_columns(X: sp.spmatrix, mask_or_index) -> sp.csr_matrix:
    """Column subset of a sparse matrix returned as canonical CSR.

    Goes through CSC so that the subset is a plain row selection on the
    transposed structure, which avoids the CSR fancy-column-indexing kernel
    (observed to crash on large matrices with some SciPy/NumPy builds).
    """
    idx = np.flatnonzero(mask_or_index) if np.asarray(mask_or_index).dtype == bool else np.asarray(mask_or_index, dtype=np.int64)
    Xc = sp.csc_matrix(X)
    Xc.sum_duplicates()
    out = sp.csr_matrix(Xc[:, idx])
    out.sum_duplicates()
    out.sort_indices()
    return out
