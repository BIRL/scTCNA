"""Normalisation.

TCNA computes gene noise on normalised bulk expression (TPM/FPKM, linear scale)
and DEPTH/DEPTH2 on log2(x + 1) values. The single-cell equivalents used here
are counts per 10,000 (CP10K, library-size normalisation; scTCNA assumption)
and log2(CP10K + 1).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def normalize_cp10k(counts: sp.spmatrix, target_sum: float = 1e4) -> sp.csr_matrix:
    """Scale each cell (row) to ``target_sum`` total counts. Zero-count cells stay zero."""
    X = sp.csr_matrix(counts, dtype=np.float32)
    totals = np.asarray(X.sum(axis=1)).ravel()
    scale = np.where(totals > 0, target_sum / np.where(totals > 0, totals, 1.0), 0.0).astype(np.float32)
    out = sp.diags(scale) @ X
    return sp.csr_matrix(out, dtype=np.float32)


def log2p1(X: sp.spmatrix) -> sp.csr_matrix:
    """log2(x + 1) applied to the stored entries (zeros stay zero)."""
    out = sp.csr_matrix(X, dtype=np.float32, copy=True)
    out.data = np.log2(out.data + 1.0).astype(np.float32)
    return out
