"""Loader tests on tiny 3CA and TISCH2 bundles written to tmp_path (sctcna.io)."""

from __future__ import annotations

import shutil
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse as sp

from sctcna.io import (
    MALIGNANT,
    NON_MALIGNANT,
    OBS_COLUMNS,
    UNASSIGNED,
    find_3ca_dir,
    load_3ca,
    load_dataset,
    load_tisch2,
    open_3ca_matrix,
    read_3ca_cells,
    read_mtx_cells_subset,
    read_mtx_header,
    read_tisch2_h5,
    recover_counts_from_lognorm,
    select_cells,
)

GENES = ["TP53", "MT-CO1", "ACTB", "GAPDH"]
# 3CA leaves the cell_type empty for cells it did not assign
CELLS = pd.DataFrame(
    {
        "cell_name": [f"cell{i}" for i in range(7)],
        "sample": ["S1", "S1", "S2", "S2", "S3", "S3", "S3"],
        "cell_type": ["Malignant", "T_cell", "", "Malignant", "Fibroblast", "Malignant", "Malignant"],
    }
)
SAMPLES = pd.DataFrame({"sample": ["S1", "S2", "S3"], "patient": ["PA", "PA", "PB"], "technology": ["10x", "10x", "10x"]})
CELL_PATIENTS = ["PA", "PA", "PA", "PA", "PB", "PB", "PB"]
CELL_COMPARTMENTS = [MALIGNANT, NON_MALIGNANT, UNASSIGNED, MALIGNANT, NON_MALIGNANT, MALIGNANT, MALIGNANT]


# ---------------------------------------------------------------------------
# 3CA
# ---------------------------------------------------------------------------
@pytest.fixture
def counts_genes_by_cells() -> np.ndarray:
    rng = np.random.default_rng(11)
    M = rng.integers(0, 6, size=(len(GENES), len(CELLS)))
    M[0, 0] = 0
    return M


@pytest.fixture
def bundle_3ca(tmp_path, counts_genes_by_cells) -> Path:
    """A study directory with the bundle one level down, as in the 3CA tarballs."""
    bundle = tmp_path / "Data_Study" / "Study"
    bundle.mkdir(parents=True)
    (bundle / "Genes.txt").write_text("".join(f'"{g}"\n' for g in GENES))
    CELLS.to_csv(bundle / "Cells.csv", index=False)
    SAMPLES.to_csv(bundle / "Samples.csv", index=False)
    scipy.io.mmwrite(bundle / "Exp_data_UMIcounts.mtx", sp.coo_matrix(counts_genes_by_cells), field="integer")
    return tmp_path / "Data_Study"


def test_find_3ca_dir_handles_nested_layout(bundle_3ca, tmp_path):
    assert find_3ca_dir(bundle_3ca) == bundle_3ca / "Study"
    assert find_3ca_dir(bundle_3ca / "Study") == bundle_3ca / "Study"
    with pytest.raises(FileNotFoundError):
        find_3ca_dir(tmp_path / "missing")


def test_read_3ca_cells_assigns_patient_and_compartment(bundle_3ca):
    cells = read_3ca_cells(bundle_3ca / "Study")
    assert cells["patient"].tolist() == CELL_PATIENTS
    assert cells["compartment"].tolist() == CELL_COMPARTMENTS
    # without Samples.csv the sample id doubles as patient id
    (bundle_3ca / "Study" / "Samples.csv").unlink()
    assert read_3ca_cells(bundle_3ca / "Study")["patient"].tolist() == CELLS["sample"].tolist()


def test_load_3ca_returns_counts_for_kept_cells(bundle_3ca, counts_genes_by_cells):
    adata = load_3ca(bundle_3ca, "DS", "CT")
    kept = (CELLS["cell_type"] != "").to_numpy()
    assert adata.shape == (6, 4)
    assert adata.obs_names.tolist() == CELLS.loc[kept, "cell_name"].tolist()
    assert adata.var_names.tolist() == GENES  # quotes stripped
    assert adata.X.dtype == np.int32 and adata.X.format == "csr"
    assert np.array_equal(adata.X.toarray(), counts_genes_by_cells[:, kept].T)
    assert list(adata.obs.columns) == OBS_COLUMNS
    assert adata.obs["compartment"].tolist() == [c for c in CELL_COMPARTMENTS if c != UNASSIGNED]
    assert adata.obs["patient"].tolist() == [p for p, c in zip(CELL_PATIENTS, CELL_COMPARTMENTS, strict=True) if c != UNASSIGNED]
    assert adata.obs["dataset"].unique().tolist() == ["DS"]
    assert adata.obs["cancer_type"].unique().tolist() == ["CT"]
    assert adata.uns["source"]["atlas"] == "3CA" and adata.uns["source"]["values"] == "UMI"
    full = load_3ca(bundle_3ca, "DS", "CT", drop_unassigned=False)
    assert full.n_obs == 7
    assert full.obs["compartment"].tolist() == CELL_COMPARTMENTS
    assert np.array_equal(full.X.toarray(), counts_genes_by_cells.T)


def test_load_3ca_checks_dimensions(bundle_3ca):
    (bundle_3ca / "Study" / "Genes.txt").write_text('"ONLY_ONE"\n')
    with pytest.raises(ValueError):
        load_3ca(bundle_3ca, "DS", "CT")


def test_select_cells_caps_groups_deterministically(bundle_3ca):
    cells = read_3ca_cells(bundle_3ca / "Study")
    assert select_cells(cells, None).all()
    assert select_cells(cells, 0).all()
    keep = select_cells(cells, 1, seed=0)
    assert keep.dtype == bool and keep.shape == (7,)
    assert np.array_equal(keep, select_cells(cells, 1, seed=0))
    sizes = cells[keep].groupby(["patient", "compartment"]).size()
    assert len(sizes) == 5 and (sizes == 1).all()
    assert keep[[1, 2, 4]].all()  # single-cell groups are never touched
    capped = load_3ca(bundle_3ca, "DS", "CT", max_cells_per_group=1, seed=3)
    assert capped.n_obs == 4
    assert capped.obs.groupby(["patient", "compartment"], observed=True).size().max() == 1


def test_read_mtx_cells_subset(bundle_3ca, counts_genes_by_cells):
    mtx = bundle_3ca / "Study" / "Exp_data_UMIcounts.mtx"
    n_genes, n_cells, nnz, is_integer, n_header = read_mtx_header(mtx)
    assert (n_genes, n_cells, nnz, is_integer) == (4, 7, int((counts_genes_by_cells != 0).sum()), True)
    assert n_header >= 2
    empty = read_mtx_cells_subset(mtx, np.zeros(7, dtype=bool))
    assert empty.shape == (0, 4) and empty.nnz == 0 and empty.format == "csr"
    with pytest.raises(ValueError):
        read_mtx_cells_subset(mtx, np.ones(6, dtype=bool))
    mask = np.array([True, False, True, False, False, False, True])
    X = read_mtx_cells_subset(mtx, mask, chunksize=3)  # forces several chunks
    assert X.dtype == np.int32
    assert np.array_equal(X.toarray(), counts_genes_by_cells[:, mask].T)
    # an open text stream works the same way and is left open
    with open(mtx) as fh:
        assert np.array_equal(read_mtx_cells_subset(fh, mask).toarray(), X.toarray())
        assert not fh.closed


def test_load_3ca_streams_the_matrix_from_the_archive(bundle_3ca, counts_genes_by_cells):
    bundle = bundle_3ca / "Study"
    mtx = bundle / "Exp_data_UMIcounts.mtx"
    archive = bundle_3ca / "Data_Study.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        tf.add(mtx, arcname="Study/Exp_data_UMIcounts.mtx")
    mtx.unlink()
    with pytest.raises(FileNotFoundError):
        load_3ca(bundle_3ca, "DS", "CT", archive=bundle_3ca / "missing.tar.gz")
    with open_3ca_matrix(bundle, archive) as (stream, name):
        assert name == "Exp_data_UMIcounts.mtx"
        assert _read_header(stream) == (4, 7)
    kept = (CELLS["cell_type"] != "").to_numpy()
    adata = load_3ca(bundle_3ca, "DS", "CT")  # the Data_*.tar.gz next to the bundle is found on its own
    assert np.array_equal(adata.X.toarray(), counts_genes_by_cells[:, kept].T)
    assert adata.uns["source"]["matrix"] == "Exp_data_UMIcounts.mtx" and adata.uns["source"]["values"] == "UMI"
    explicit = load_3ca(bundle_3ca, "DS", "CT", archive=archive)
    assert np.array_equal(explicit.X.toarray(), adata.X.toarray())


def _read_header(stream) -> tuple[int, int]:
    """(n_genes, n_cells) from the dimension line of an open MatrixMarket stream."""
    for line in stream:
        if not line.startswith("%"):
            n_genes, n_cells, _ = (int(x) for x in line.split())
            return n_genes, n_cells
    raise AssertionError("no dimension line")


# ---------------------------------------------------------------------------
# TISCH2
# ---------------------------------------------------------------------------
def _random_counts(rng: np.random.Generator, n_cells: int, n_genes: int) -> np.ndarray:
    """Random counts where every cell has at least one count of exactly 1, which the inversion is anchored on."""
    counts = rng.poisson(0.8, size=(n_cells, n_genes))
    counts[np.arange(n_cells), rng.integers(0, n_genes, n_cells)] = 1
    return counts


def _log_normalise(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    totals = counts.sum(axis=1, keepdims=True)
    return np.log1p(counts / np.where(totals > 0, totals, 1) * target_sum)


def _write_tisch2_h5(path: Path, log_norm_cells_by_genes: np.ndarray, barcodes: list[str], genes: list[str]) -> Path:
    """10x-style HDF5 layout: genes x cells CSC under /matrix with a /matrix/features group."""
    Xc = sp.csc_matrix(log_norm_cells_by_genes.T.astype(np.float32))
    with h5py.File(path, "w") as fh:
        g = fh.create_group("matrix")
        g.create_dataset("data", data=Xc.data)
        g.create_dataset("indices", data=Xc.indices)
        g.create_dataset("indptr", data=Xc.indptr)
        g.create_dataset("shape", data=np.array(Xc.shape, dtype=np.int64))
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        f = g.create_group("features")
        f.create_dataset("name", data=np.array(genes, dtype="S"))
        f.create_dataset("id", data=np.array([f"ENSG{i:05d}" for i in range(len(genes))], dtype="S"))
        f.create_dataset("feature_type", data=np.array(["Gene Expression"] * len(genes), dtype="S"))
        f.create_dataset("genome", data=np.array(["GRCh38"] * len(genes), dtype="S"))
        f.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))
    return path


@pytest.fixture
def tisch2_bundle(tmp_path) -> tuple[Path, Path, np.ndarray, list[str], list[str]]:
    rng = np.random.default_rng(5)
    n_cells, n_genes = 12, 40
    counts = _random_counts(rng, n_cells, n_genes)
    barcodes = [f"BC{i}@S{i % 2 + 1}" for i in range(n_cells)]
    genes = [f"GENE{i}" for i in range(n_genes)]
    h5 = _write_tisch2_h5(tmp_path / "expression.h5", _log_normalise(counts), barcodes, genes)
    malignant = [i % 3 == 0 for i in range(n_cells)]
    meta = pd.DataFrame(
        {
            "Cell": barcodes,
            "Celltype (malignancy)": ["Malignant cells" if m else "Immune cells" for m in malignant],
            "Celltype (major-lineage)": ["Malignant" if m else "CD8T" for m in malignant],
            "Celltype (minor-lineage)": ["Malignant" if m else "CD8Tex" for m in malignant],
            "Patient": [f"P{i % 2 + 1}" for i in range(n_cells)],
            "Sample": [f"S{i % 2 + 1}" for i in range(n_cells)],
        }
    )
    meta.to_csv(tmp_path / "CellMetainfo_table.tsv", sep="\t", index=False)
    return h5, tmp_path / "CellMetainfo_table.tsv", counts, barcodes, genes


def test_read_tisch2_h5_returns_cells_by_genes(tisch2_bundle):
    h5, _, counts, barcodes, genes = tisch2_bundle
    X, bc, gn = read_tisch2_h5(h5)
    assert X.shape == (12, 40) and X.format == "csr"
    assert bc == barcodes and gn == genes
    np.testing.assert_allclose(X.toarray(), _log_normalise(counts), rtol=1e-6, atol=1e-6)


def test_load_tisch2_recovers_integer_counts_and_compartments(tisch2_bundle):
    h5, meta, counts, barcodes, genes = tisch2_bundle
    adata = load_tisch2(h5, meta, "T2", "CT")
    assert adata.obs_names.tolist() == barcodes
    assert adata.var_names.tolist() == genes
    assert adata.X.dtype == np.int32
    assert np.array_equal(adata.X.toarray(), counts)
    assert adata.uns["source"]["fraction_cells_integer"] == 1.0
    assert adata.uns["source"]["atlas"] == "TISCH2"
    assert list(adata.obs.columns) == OBS_COLUMNS
    assert adata.obs["compartment"].tolist() == [MALIGNANT if i % 3 == 0 else NON_MALIGNANT for i in range(12)]
    assert adata.obs["cell_type"].tolist() == ["Malignant" if i % 3 == 0 else "CD8T" for i in range(12)]
    assert adata.obs["patient"].tolist() == [f"P{i % 2 + 1}" for i in range(12)]
    assert adata.obs["sample"].tolist() == [f"S{i % 2 + 1}" for i in range(12)]
    assert adata.obs["dataset"].unique().tolist() == ["T2"] and adata.obs["cancer_type"].unique().tolist() == ["CT"]
    capped = load_tisch2(h5, meta, "T2", "CT", max_cells_per_group=2, seed=1)
    assert capped.obs.groupby(["patient", "compartment"], observed=True).size().max() == 2
    assert capped.n_obs == 8


def test_recover_counts_from_lognorm_roundtrip():
    rng = np.random.default_rng(1)
    for n_cells, n_genes in ((5, 30), (50, 120), (200, 25)):
        counts = _random_counts(rng, n_cells, n_genes)
        recovered, frac = recover_counts_from_lognorm(sp.csr_matrix(_log_normalise(counts)))
        assert frac == 1.0
        assert recovered.dtype == np.int32 and recovered.format == "csr"
        assert np.array_equal(recovered.toarray(), counts)
    counts = _random_counts(rng, 20, 50)
    recovered, frac = recover_counts_from_lognorm(sp.csr_matrix(_log_normalise(counts, 1e6)), target_sum=1e6)
    assert frac == 1.0 and np.array_equal(recovered.toarray(), counts)


def test_recover_counts_flags_cells_without_a_unit_count():
    rng = np.random.default_rng(2)
    counts = _random_counts(rng, 10, 30)
    counts[0] = 0
    counts[0, :3] = [2, 3, 2]  # smallest non-zero count is 2: the inferred total is half the true one
    counts[1] = 0  # an empty cell stays empty and is not counted as a failure
    recovered, frac = recover_counts_from_lognorm(sp.csr_matrix(_log_normalise(counts)))
    assert frac == 0.9
    assert np.array_equal(recovered.toarray()[2:], counts[2:])
    assert not recovered[1].nnz
    assert not np.array_equal(recovered.toarray()[0], counts[0])


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
def test_load_dataset_prefers_cached_h5ad(tmp_path, bundle_3ca):
    raw = tmp_path / "raw"
    ds_dir = raw / "3ca" / "DS"
    ds_dir.mkdir(parents=True)
    entry = {"id": "DS", "source": "3CA", "cancer_type": "CT", "data_url": "https://example.org/Data_Study.tar.gz"}
    with pytest.raises(FileNotFoundError):
        load_dataset(entry, raw)
    for p in (bundle_3ca / "Study").iterdir():
        shutil.copy(p, ds_dir / p.name)
    from_bundle = load_dataset(entry, raw, max_cells_per_group=None, seed=0)
    assert from_bundle.n_obs == 6 and from_bundle.uns["source"]["atlas"] == "3CA"
    cached = from_bundle[:2].copy()
    cached.uns["source"] = {"atlas": "cache"}
    cached.write_h5ad(ds_dir / "DS.h5ad")
    from_cache = load_dataset(entry, raw)
    assert from_cache.n_obs == 2 and from_cache.uns["source"]["atlas"] == "cache"
    assert from_cache.obs_names.tolist() == from_bundle.obs_names[:2].tolist()
    with pytest.raises(ValueError):
        load_dataset({**entry, "source": "OTHER"}, raw)


def test_load_dataset_dispatches_tisch2(tmp_path, tisch2_bundle):
    h5, meta, counts, _, _ = tisch2_bundle
    raw = tmp_path / "raw"
    ds_dir = raw / "tisch2" / "T2"
    ds_dir.mkdir(parents=True)
    shutil.copy(h5, ds_dir / h5.name)
    shutil.copy(meta, ds_dir / meta.name)
    entry = {
        "id": "T2",
        "source": "TISCH2",
        "cancer_type": "CT",
        "data_url": f"https://example.org/static/data/T2/{h5.name}?download=1",
        "meta_url": f"https://example.org/static/data/T2/{meta.name}",
    }
    adata = load_dataset(entry, raw)
    assert adata.uns["source"]["atlas"] == "TISCH2"
    assert np.array_equal(adata.X.toarray(), counts)
