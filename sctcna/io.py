"""Loaders that turn 3CA and TISCH2 study bundles into AnnData objects.

The returned AnnData is cells x genes with raw UMI counts in ``X`` (CSR, int32)
and the following ``obs`` columns: ``cell_name``, ``sample``, ``patient``,
``cell_type``, ``compartment`` (``malignant`` / ``non_malignant``), ``dataset``
and ``cancer_type``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger(__name__)

MALIGNANT = "malignant"
NON_MALIGNANT = "non_malignant"
UNASSIGNED = "unassigned"

OBS_COLUMNS = ["cell_name", "sample", "patient", "cell_type", "compartment", "dataset", "cancer_type"]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def select_cells(obs: pd.DataFrame, max_cells_per_group: int | None, seed: int = 0) -> np.ndarray:
    """Boolean mask capping cells per (patient, compartment) group.

    scTCNA assumption: groups larger than the cap are randomly subsampled so
    that very large studies stay tractable. The cap is applied before QC.
    """
    keep = np.ones(len(obs), dtype=bool)
    if not max_cells_per_group:
        return keep
    rng = np.random.default_rng(seed)
    for _, idx in obs.groupby(["patient", "compartment"], observed=True).indices.items():
        if len(idx) > max_cells_per_group:
            drop = rng.choice(idx, size=len(idx) - max_cells_per_group, replace=False)
            keep[drop] = False
    return keep


def _finalise(X: sp.csr_matrix, obs: pd.DataFrame, genes: list[str], dataset_id: str, cancer_type: str) -> ad.AnnData:
    obs = obs.copy()
    obs["dataset"] = dataset_id
    obs["cancer_type"] = cancer_type
    obs.index = obs["cell_name"].astype(str).values
    obs.index.name = None
    var = pd.DataFrame(index=pd.Index(genes, name=None).astype(str))
    adata = ad.AnnData(X=X, obs=obs[OBS_COLUMNS + [c for c in obs.columns if c not in OBS_COLUMNS]], var=var)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


# ---------------------------------------------------------------------------
# 3CA
# ---------------------------------------------------------------------------
def find_3ca_dir(study_dir: str | Path) -> Path:
    """Locate the folder holding the 3CA expression files (handles the nested tar layout)."""
    study_dir = Path(study_dir)
    candidates = [study_dir, *[p for p in study_dir.iterdir() if p.is_dir()]] if study_dir.is_dir() else []
    for c in candidates:
        if (c / "Cells.csv").exists() and (c / "Genes.txt").exists():
            return c
    raise FileNotFoundError(f"no 3CA bundle (Cells.csv + Genes.txt) found under {study_dir}")


def _3ca_matrix_path(bundle: Path) -> Path:
    for name in ("Exp_data_UMIcounts.mtx", "Exp_data_TPM.mtx"):
        if (bundle / name).exists():
            return bundle / name
    raise FileNotFoundError(f"no Exp_data_*.mtx in {bundle}")


def read_3ca_cells(bundle: str | Path) -> pd.DataFrame:
    """Read Cells.csv (+ Samples.csv) and derive ``patient`` and ``compartment``.

    3CA labels malignant cells as ``Malignant``; every other non-empty
    ``cell_type`` is a curated non-malignant lineage. Cells with an empty
    ``cell_type`` were left unassigned by 3CA and are marked ``unassigned``.
    """
    bundle = Path(bundle)
    cells = pd.read_csv(bundle / "Cells.csv", dtype=str, keep_default_na=False)
    cells.columns = [c.strip() for c in cells.columns]
    if "patient" not in cells.columns or (cells["patient"] == "").all():
        samples_path = bundle / "Samples.csv"
        if samples_path.exists():
            samples = pd.read_csv(samples_path, dtype=str, keep_default_na=False)
            s2p = dict(zip(samples["sample"], samples.get("patient", samples["sample"]), strict=True))
            cells["patient"] = cells["sample"].map(s2p).fillna(cells["sample"])
        else:
            cells["patient"] = cells["sample"]
    cells["patient"] = cells["patient"].replace("", pd.NA).fillna(cells["sample"]).astype(str)
    ct = cells["cell_type"].astype(str).str.strip()
    cells["cell_type"] = ct
    cells["compartment"] = np.where(ct == "Malignant", MALIGNANT, np.where(ct == "", UNASSIGNED, NON_MALIGNANT))
    return cells


def read_mtx_header(path: str | Path) -> tuple[int, int, int, bool, int]:
    """Return (n_rows, n_cols, nnz, is_integer, n_header_lines) of a coordinate MatrixMarket file."""
    n_header = 0
    is_integer = False
    with open(path) as fh:
        for line in fh:
            n_header += 1
            if line.startswith("%"):
                if n_header == 1:
                    is_integer = "integer" in line
                continue
            nrow, ncol, nnz = (int(x) for x in line.split())
            return nrow, ncol, nnz, is_integer, n_header
    raise ValueError(f"{path}: no dimension line found")


def read_mtx_cells_subset(path: str | Path, keep_cells: np.ndarray, chunksize: int = 20_000_000) -> sp.csr_matrix:
    """Stream a genes x cells coordinate MatrixMarket file, keeping only ``keep_cells``.

    Returns a cells x genes CSR matrix (rows ordered as the kept cells appear in
    the file's column order). Streaming keeps memory proportional to the kept
    entries rather than the whole file.
    """
    n_genes, n_cells, _nnz, is_integer, n_header = read_mtx_header(path)
    keep_cells = np.asarray(keep_cells, dtype=bool)
    if keep_cells.shape[0] != n_cells:
        raise ValueError(f"keep_cells has {keep_cells.shape[0]} entries but matrix has {n_cells} cells")
    new_index = np.cumsum(keep_cells) - 1
    n_keep = int(keep_cells.sum())
    val_dtype = np.int32 if is_integer else np.float32
    rows, cols, vals = [], [], []
    reader = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["g", "c", "v"],
        skiprows=n_header,
        dtype={"g": np.int32, "c": np.int32, "v": val_dtype},
        chunksize=chunksize,
        engine="c",
    )
    for chunk in reader:
        c = chunk["c"].to_numpy() - 1
        m = keep_cells[c]
        if not m.any():
            continue
        rows.append(new_index[c[m]])
        cols.append(chunk["g"].to_numpy()[m] - 1)
        vals.append(chunk["v"].to_numpy()[m])
    if rows:
        r = np.concatenate(rows)
        cidx = np.concatenate(cols)
        v = np.concatenate(vals)
    else:
        r = np.zeros(0, dtype=np.int32)
        cidx = np.zeros(0, dtype=np.int32)
        v = np.zeros(0, dtype=val_dtype)
    X = sp.coo_matrix((v, (r, cidx)), shape=(n_keep, n_genes)).tocsr()
    X.sum_duplicates()
    X.eliminate_zeros()
    return X


def load_3ca(
    study_dir: str | Path,
    dataset_id: str,
    cancer_type: str,
    max_cells_per_group: int | None = None,
    seed: int = 0,
    drop_unassigned: bool = True,
) -> ad.AnnData:
    """Load a 3CA study bundle into AnnData (raw UMI counts)."""
    bundle = find_3ca_dir(study_dir)
    cells = read_3ca_cells(bundle)
    genes = [g.strip().strip('"') for g in open(bundle / "Genes.txt") if g.strip()]
    mtx = _3ca_matrix_path(bundle)
    n_genes, n_cells, *_ = read_mtx_header(mtx)
    if n_cells != len(cells) or n_genes != len(genes):
        raise ValueError(f"{mtx}: matrix is {n_genes}x{n_cells} but Genes.txt/Cells.csv have {len(genes)}/{len(cells)} entries")
    keep = np.ones(len(cells), dtype=bool)
    if drop_unassigned:
        keep &= (cells["compartment"] != UNASSIGNED).to_numpy()
    keep &= select_cells(cells, max_cells_per_group, seed)
    log.info("%s: keeping %d of %d cells (%s)", dataset_id, keep.sum(), len(cells), mtx.name)
    X = read_mtx_cells_subset(mtx, keep)
    if mtx.name == "Exp_data_TPM.mtx":
        log.warning("%s: 3CA bundle provides TPM (non-UMI technology); counts-based steps are approximate", dataset_id)
    obs = cells.loc[keep].reset_index(drop=True)
    adata = _finalise(X, obs, genes, dataset_id, cancer_type)
    adata.uns["source"] = {"atlas": "3CA", "bundle": str(bundle), "matrix": mtx.name, "values": "UMI" if "UMI" in mtx.name else "TPM"}
    return adata


# ---------------------------------------------------------------------------
# TISCH2
# ---------------------------------------------------------------------------
def recover_counts_from_lognorm(X: sp.csr_matrix, target_sum: float = 1e4, tol: float = 0.02) -> tuple[sp.csr_matrix, float]:
    """Invert Seurat-style log1p(counts / total * target_sum) normalisation.

    For a cell whose smallest non-zero count is 1, the smallest non-zero
    normalised value equals target_sum / total, which gives the total and hence
    every count. The fraction of cells whose recovered values are all within
    ``tol`` of an integer is returned so callers can verify the inversion.

    scTCNA assumption: used only for TISCH2 matrices, which are distributed
    without raw counts.
    """
    X = sp.csr_matrix(X, dtype=np.float64)
    X.eliminate_zeros()
    cp = np.expm1(X.data)
    n = X.shape[0]
    mins = np.full(n, np.inf)
    rows = np.repeat(np.arange(n), np.diff(X.indptr))
    np.minimum.at(mins, rows, cp)
    totals = np.where(np.isfinite(mins) & (mins > 0), target_sum / mins, 0.0)
    counts = cp * totals[rows] / target_sum
    rounded = np.rint(counts)
    ok = np.abs(counts - rounded) <= tol
    cell_ok = np.ones(n, dtype=bool)
    np.logical_and.at(cell_ok, rows, ok)
    frac_ok = float(cell_ok.mean()) if n else 1.0
    out = sp.csr_matrix((rounded.astype(np.int32), X.indices.copy(), X.indptr.copy()), shape=X.shape)
    out.eliminate_zeros()
    return out, frac_ok


def read_tisch2_h5(path: str | Path) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """Read a TISCH2 10x-style HDF5 file. Returns (cells x genes CSR, barcodes, genes)."""
    with h5py.File(path, "r") as fh:
        g = fh["matrix"]
        shape = tuple(int(x) for x in g["shape"][:])  # (genes, cells)
        X = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        barcodes = [b.decode() if isinstance(b, bytes) else str(b) for b in g["barcodes"][:]]
        genes = [b.decode() if isinstance(b, bytes) else str(b) for b in g["features"]["name"][:]]
    return X.T.tocsr(), barcodes, genes


def load_tisch2(
    h5_path: str | Path,
    meta_path: str | Path,
    dataset_id: str,
    cancer_type: str,
    max_cells_per_group: int | None = None,
    seed: int = 0,
    target_sum: float = 1e4,
) -> ad.AnnData:
    """Load a TISCH2 dataset (expression .h5 + CellMetainfo table) into AnnData."""
    X, barcodes, genes = read_tisch2_h5(h5_path)
    meta = pd.read_csv(meta_path, sep="\t", dtype=str, keep_default_na=False)
    meta = meta.set_index("Cell").reindex(barcodes)
    mal = meta["Celltype (malignancy)"].astype(str)
    obs = pd.DataFrame(
        {
            "cell_name": barcodes,
            "sample": meta.get("Sample", pd.Series(index=meta.index, dtype=str)).astype(str).values,
            "patient": meta.get("Patient", pd.Series(index=meta.index, dtype=str)).astype(str).values,
            "cell_type": meta["Celltype (major-lineage)"].astype(str).values,
        }
    )
    obs["patient"] = obs["patient"].replace({"": pd.NA, "nan": pd.NA}).fillna(obs["sample"])
    obs["compartment"] = np.where(mal.values == "Malignant cells", MALIGNANT, NON_MALIGNANT)
    keep = select_cells(obs, max_cells_per_group, seed)
    X = X[keep]
    obs = obs.loc[keep].reset_index(drop=True)
    counts, frac_ok = recover_counts_from_lognorm(X, target_sum=target_sum)
    if frac_ok < 0.95:
        log.warning("%s: only %.1f%% of cells inverted to integer counts; treat count-based steps with caution", dataset_id, 100 * frac_ok)
    adata = _finalise(counts, obs, genes, dataset_id, cancer_type)
    adata.uns["source"] = {"atlas": "TISCH2", "h5": str(h5_path), "values": "UMI recovered from log1p(CP10K)", "fraction_cells_integer": frac_ok}
    return adata


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
def load_dataset(entry: dict, raw_dir: str | Path, max_cells_per_group: int | None = None, seed: int = 0) -> ad.AnnData:
    """Load a registry entry from ``raw_dir/<source>/<id>/`` (or its cached .h5ad)."""
    raw_dir = Path(raw_dir)
    ds_dir = raw_dir / entry["source"].lower() / entry["id"]
    cached = ds_dir / f"{entry['id']}.h5ad"
    if cached.exists():
        return ad.read_h5ad(cached)
    if entry["source"] == "3CA":
        return load_3ca(ds_dir, entry["id"], entry["cancer_type"], max_cells_per_group, seed)
    if entry["source"] == "TISCH2":
        h5 = ds_dir / Path(entry["data_url"].split("?")[0]).name
        meta = ds_dir / Path(entry["meta_url"].split("?")[0]).name
        return load_tisch2(h5, meta, entry["id"], entry["cancer_type"], max_cells_per_group, seed)
    raise ValueError(f"unknown source {entry['source']}")
