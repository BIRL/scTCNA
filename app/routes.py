"""JSON API routes and CSV downloads.

Every JSON body goes through ``clean`` so NaN/inf and numpy scalars never reach
the client; ``json.dumps(allow_nan=False)`` then guarantees a hard failure
instead of an invalid document should anything slip through.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterator
from typing import Annotated, Any

import numpy as np
import pandas as pd
from fastapi import APIRouter, Depends, HTTPException, Query, Request
from fastapi.responses import Response, StreamingResponse

from .store import Store

router = APIRouter(prefix="/api")

CSV_TABLES = ("gene_summary", "top_genes", "sample_scores", "cell_scores", "datasets", "gene_patients", "heatmap", "noise_contrasts")


def clean(obj: Any) -> Any:
    """Recursively convert a response payload to JSON-safe Python values (NaN/inf -> None)."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return clean(obj.tolist())
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if math.isfinite(obj) else None
    if obj is None or obj is pd.NA or obj is pd.NaT:
        return None
    if isinstance(obj, np.str_):
        return str(obj)
    return obj


def json_response(payload: Any, status_code: int = 200) -> Response:
    body = json.dumps(clean(payload), allow_nan=False, separators=(",", ":"))
    return Response(content=body, status_code=status_code, media_type="application/json")


def get_store(request: Request) -> Store:
    store = request.app.state.store
    if store is None:
        raise HTTPException(status_code=503, detail=request.app.state.data_error)
    return store


StoreDep = Annotated[Store, Depends(get_store)]


# ---------------------------------------------------------------------------- overview
@router.api_route("/health", methods=["GET", "HEAD"])
def health(request: Request) -> Response:
    """200 with build info when the tables are loaded; 503 otherwise, so a host's health check fails a deploy without data."""
    store = request.app.state.store
    if store is None:
        return json_response({"status": "no-data", "n_datasets": 0, "n_genes": 0, "build": None, "detail": request.app.state.data_error}, status_code=503)
    return json_response({"status": "ok", **store.health()})


@router.get("/cancer-types")
def cancer_types(store: StoreDep) -> Response:
    return json_response(store.cancer_type_rows())


@router.get("/datasets")
def datasets(store: StoreDep) -> Response:
    return json_response(store.dataset_rows())


# ---------------------------------------------------------------------------- genes
@router.get("/genes/search")
def search_genes(store: StoreDep, q: str = "", limit: Annotated[int, Query(ge=1, le=100)] = 20) -> Response:
    return json_response({"query": q, "results": store.search_genes(q, limit)})


@router.get("/genes/{gene}")
def gene(store: StoreDep, gene: str) -> Response:
    return json_response(store.gene(gene))


@router.get("/genes/{gene}/patients")
def gene_patients(store: StoreDep, gene: str, cancer_type: str) -> Response:
    return json_response(store.gene_patients(gene, cancer_type))


@router.get("/genes/{gene}/celltypes")
def gene_celltypes(store: StoreDep, gene: str, cancer_type: str) -> Response:
    return json_response(store.gene_celltypes(gene, cancer_type))


# ---------------------------------------------------------------------------- per cancer type
@router.get("/cancer-types/{cancer_type}/top-genes")
def top_genes(store: StoreDep, cancer_type: str, metric: str, direction: str = "noisy", n: Annotated[int, Query(ge=1, le=500)] = 50, expressed_only: bool = True) -> Response:
    return json_response(store.top_genes(cancer_type, metric, direction, n, expressed_only))


@router.get("/cancer-types/{cancer_type}/comparison")
def comparison(store: StoreDep, cancer_type: str, metric: str = "cv", expressed_only: bool = True) -> Response:
    return json_response(store.comparison(cancer_type, metric, expressed_only))


@router.get("/cancer-types/{cancer_type}/tumor-heterogeneity")
def tumor_heterogeneity(store: StoreDep, cancer_type: str, score: str = "depth2", level: str = "cell") -> Response:
    return json_response(store.tumor_heterogeneity(cancer_type, score, level))


@router.get("/heatmap")
def heatmap(store: StoreDep, genes: str = "", metric: str = "log2fc_cv_median", expressed_only: bool = False) -> Response:
    return json_response(store.heatmap(genes, metric, expressed_only))


@router.get("/tumor-heterogeneity/summary")
def heterogeneity_summary(store: StoreDep) -> Response:
    return json_response(store.heterogeneity_summary())


@router.get("/noise-contrasts")
def noise_contrasts(store: StoreDep, cancer_type: str | None = None) -> Response:
    return json_response(store.contrast_summary(cancer_type))


@router.get("/cancer-types/{cancer_type}/noise-contrasts")
def noise_contrasts_patients(store: StoreDep, cancer_type: str) -> Response:
    return json_response(store.contrast_patients(cancer_type))


# ---------------------------------------------------------------------------- downloads
def _csv_chunks(df: pd.DataFrame, rows_per_chunk: int = 50_000) -> Iterator[str]:
    yield df.iloc[:0].to_csv(index=False)
    for start in range(0, len(df), rows_per_chunk):
        yield df.iloc[start : start + rows_per_chunk].to_csv(index=False, header=False)


def csv_response(df: pd.DataFrame, filename: str) -> StreamingResponse:
    return StreamingResponse(_csv_chunks(df), media_type="text/csv", headers={"Content-Disposition": f'attachment; filename="{filename}"'})


def _require(value: str | None, name: str, table: str) -> str:
    if value is None or not value.strip():
        raise HTTPException(status_code=422, detail=f"query parameter '{name}' is required for the {table} download")
    return value


@router.get("/download/{table}.csv")
def download(
    store: StoreDep,
    table: str,
    cancer_type: str | None = None,
    gene: str | None = None,
    genes: str | None = None,
    metric: str | None = None,
    direction: str = "noisy",
    n: Annotated[int, Query(ge=1, le=500)] = 50,
    expressed_only: bool | None = None,
) -> StreamingResponse:
    if table not in CSV_TABLES:
        raise HTTPException(status_code=404, detail=f"unknown table: {table}; expected one of {', '.join(CSV_TABLES)}")
    suffix = f"_{cancer_type}" if cancer_type else ""
    if table == "datasets":
        return csv_response(store.datasets, "datasets.csv")
    if table == "gene_summary":
        return csv_response(store.gene_summary_frame(cancer_type), f"gene_summary{suffix}.csv")
    if table == "sample_scores":
        return csv_response(store.sample_scores_frame(cancer_type), f"sample_scores{suffix}.csv")
    if table == "noise_contrasts":
        return csv_response(store.contrast_frame(cancer_type), f"noise_contrasts{suffix}.csv")
    if table == "cell_scores":
        _require(cancer_type, "cancer_type", table)
        return csv_response(store.cell_scores_frame(cancer_type), f"cell_scores{suffix}.csv")
    if table == "top_genes":
        _require(cancer_type, "cancer_type", table)
        _require(metric, "metric", table)
        expressed = True if expressed_only is None else expressed_only
        return csv_response(store.top_genes_frame(cancer_type, metric, direction, n, expressed), f"top_genes{suffix}_{metric}_{direction}.csv")
    if table == "gene_patients":
        _require(gene, "gene", table)
        _require(cancer_type, "cancer_type", table)
        df = store.gene_patients_frame(gene, cancer_type)
        return csv_response(df, f"gene_patients_{df['gene'].iloc[0]}{suffix}.csv")
    _require(genes, "genes", table)
    _require(metric, "metric", table)
    expressed = False if expressed_only is None else expressed_only
    return csv_response(store.heatmap_frame(genes, metric, expressed), f"heatmap_{metric}.csv")
