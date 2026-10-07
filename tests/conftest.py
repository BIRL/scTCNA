"""Shared fixtures: a small synthetic dataset run through the pipeline and aggregated into app tables."""

from __future__ import annotations

import copy
from pathlib import Path

import pytest

from sctcna.config import load_config
from sctcna.pipeline import aggregate, process_dataset
from sctcna.synthetic import make_synthetic_adata

SYN_META = {
    "id": "SYN_Test",
    "cancer_type": "SYN",
    "cancer_name": "Synthetic test data",
    "tcna_site": "Synthetic",
    "source": "synthetic",
    "study": "synthetic",
    "citation": "n/a",
    "technology": "simulated",
    "site": "n/a",
    "data_url": "",
    "enabled": True,
}
SYN2_META = {**SYN_META, "id": "SYN_Two", "cancer_type": "SY2", "cancer_name": "Second synthetic cancer"}


@pytest.fixture(scope="session")
def fast_params() -> dict:
    cfg = load_config()
    p = copy.deepcopy(cfg["params"])
    p["qc"]["run_scrublet"] = False  # Scrublet is slow and not meaningful on simulated data
    p["qc"]["min_cells_per_gene"] = 3
    p["groups"]["min_cells_per_group"] = 10
    p["groups"]["min_patients_per_test"] = 4
    p["groups"]["min_cells_per_celltype_group"] = 10
    p["downsampling"]["min_depth"] = 200
    return p


@pytest.fixture(scope="session")
def synthetic_adata():
    return make_synthetic_adata(n_patients=6, cells_per_group=120, seed=1)


@pytest.fixture(scope="session")
def processed_dir(tmp_path_factory, synthetic_adata, fast_params) -> Path:
    out = tmp_path_factory.mktemp("processed")
    process_dataset(synthetic_adata, fast_params, SYN_META, out)
    second = make_synthetic_adata(n_patients=5, cells_per_group=100, seed=2, dataset_id="SYN_Two", cancer_type="SY2")
    process_dataset(second, fast_params, SYN2_META, out)
    return out


@pytest.fixture(scope="session")
def app_data_dir(tmp_path_factory, processed_dir) -> Path:
    app_dir = tmp_path_factory.mktemp("app")
    registry = {"datasets": [SYN_META, SYN2_META]}
    aggregate(processed_dir, app_dir, registry, manifest={})
    return app_dir
