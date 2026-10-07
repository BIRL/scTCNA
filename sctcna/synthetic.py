"""Synthetic single-cell data with known noise structure, used by tests and demos."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _nb_counts(rng, mu, theta, size):
    """Negative binomial with mean mu and dispersion theta (Var = mu + mu^2/theta)."""
    if np.isinf(theta):
        return rng.poisson(mu, size=size)
    lam = rng.gamma(shape=theta, scale=mu / theta, size=size)
    return rng.poisson(lam)


def make_synthetic_adata(
    n_patients: int = 4,
    cells_per_group: int = 200,
    n_genes: int = 300,
    n_noisy_in_malignant: int = 20,
    seed: int = 0,
    dataset_id: str = "SYN_Test",
    cancer_type: str = "SYN",
) -> ad.AnnData:
    """Two compartments x n_patients with Poisson, over-dispersed and malignant-noisy genes.

    Gene classes (stored in ``adata.var['truth']``):
    * ``poisson``       : Fano = 1 in both compartments
    * ``overdispersed`` : NB with theta = 2 in both compartments (Fano >> 1)
    * ``noisy_malignant``: Poisson in non-malignant cells, NB theta = 1 in malignant cells (same mean)
    * ``constant``      : the same count in every cell (zero noise)
    * ``mito``          : Poisson genes named MT-* for the mitochondrial fraction
    """
    rng = np.random.default_rng(seed)
    classes = ["constant"] * 5 + ["mito"] * 5 + ["noisy_malignant"] * n_noisy_in_malignant + ["overdispersed"] * 40 + ["poisson"] * (n_genes - 50 - n_noisy_in_malignant)
    n_genes = len(classes)
    names = []
    for i, c in enumerate(classes):
        prefix = {"constant": "CONST", "mito": "MT-G", "noisy_malignant": "NOISY", "overdispersed": "OD", "poisson": "POIS"}[c]
        names.append(f"{prefix}{i}")
    mu = rng.lognormal(mean=0.5, sigma=1.0, size=n_genes)  # mean counts per cell
    obs_rows, blocks = [], []
    for p in range(n_patients):
        for comp in ("malignant", "non_malignant"):
            block = np.zeros((cells_per_group, n_genes), dtype=np.int32)
            for g, c in enumerate(classes):
                if c == "constant":
                    block[:, g] = 3
                elif c == "overdispersed":
                    block[:, g] = _nb_counts(rng, mu[g], 2.0, cells_per_group)
                elif c == "noisy_malignant":
                    block[:, g] = _nb_counts(rng, mu[g], 1.0 if comp == "malignant" else np.inf, cells_per_group)
                else:
                    block[:, g] = rng.poisson(mu[g], size=cells_per_group)
            blocks.append(block)
            sample = f"P{p + 1}_{'T' if comp == 'malignant' else 'N'}"
            for i in range(cells_per_group):
                obs_rows.append(
                    {
                        "cell_name": f"{sample}_c{i}",
                        "sample": f"P{p + 1}",
                        "patient": f"P{p + 1}",
                        "cell_type": "Malignant" if comp == "malignant" else rng.choice(["T_cell", "Fibroblast", "Epithelial"]),
                        "compartment": comp,
                        "dataset": dataset_id,
                        "cancer_type": cancer_type,
                    }
                )
    X = sp.csr_matrix(np.vstack(blocks))
    obs = pd.DataFrame(obs_rows)
    obs.index = obs["cell_name"].values
    var = pd.DataFrame({"truth": classes, "mu": mu}, index=names)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["source"] = {"atlas": "synthetic", "values": "UMI"}
    return adata
