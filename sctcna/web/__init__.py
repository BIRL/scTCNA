"""Web layer of scTCNA: the in-memory store of app tables and the JSON API.

``create_api(data_dir)`` returns a self-contained FastAPI application with
routes such as ``/health`` and ``/genes/{gene}``. scTCNA's own server mounts
it at ``/api`` (``app.main``); another FastAPI application, such as the TCNA
portal, can mount the same object at any prefix, for example ``/api/sc``.
Only pandas, numpy, pyarrow and FastAPI are needed at runtime.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from fastapi import FastAPI, Request

from .. import __version__
from .routes import json_response, router
from .store import Store, StoreError

__all__ = ["Store", "StoreError", "create_api", "load_store", "resolve_data_dir"]

log = logging.getLogger(__name__)

DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "app"


def resolve_data_dir(data_dir: str | Path | None = None) -> Path:
    """``data_dir`` argument, else the ``SCTCNA_DATA_DIR`` environment variable, else ``<repo>/data/app``."""
    return Path(data_dir or os.environ.get("SCTCNA_DATA_DIR") or DEFAULT_DATA_DIR)


def load_store(data_dir: Path) -> tuple[Store | None, str | None]:
    """Load the tables; on failure return (None, reason) so the API can answer 503 instead of crashing."""
    try:
        store = Store(data_dir)
        log.info("loaded scTCNA tables from %s (%d datasets)", data_dir.name, store.health()["n_datasets"])
        return store, None
    except (FileNotFoundError, ValueError) as exc:
        reason = f"app data not available: {exc}"
        log.warning("%s; the API will answer 503 until the tables exist and the server is restarted", reason)
        return None, reason


def create_api(data_dir: str | Path | None = None, store: Store | None = None, data_error: str | None = None) -> FastAPI:
    """FastAPI application with the scTCNA JSON API and CSV downloads (no static files)."""
    if store is None and data_error is None:
        store, data_error = load_store(resolve_data_dir(data_dir))
    api = FastAPI(title="scTCNA API", version=__version__, description="Single-cell transcriptional noise atlas for cancer: JSON API over the aggregated tables.")
    api.state.store = store
    api.state.data_error = data_error

    @api.exception_handler(StoreError)
    async def store_error(request: Request, exc: StoreError):
        return json_response({"detail": str(exc)}, status_code=exc.status_code)

    api.include_router(router)
    return api
