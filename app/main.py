"""FastAPI application factory.

``create_app`` loads the app tables from ``data_dir`` (argument, else the
``SCTCNA_DATA_DIR`` environment variable, else ``<repo>/data/app``) and mounts
the static frontend. A missing data directory does not break start-up: the
API then answers 503 on every route except ``/api/health``, which reports
``no-data``, so the server can come up before the tables are built.

Run with ``uvicorn app.main:app``.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles

from sctcna import __version__

from .routes import json_response, router
from .store import Store, StoreError

log = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
STATIC_DIR = Path(__file__).resolve().parent / "static"


def resolve_data_dir(data_dir: str | Path | None = None) -> Path:
    return Path(data_dir or os.environ.get("SCTCNA_DATA_DIR") or REPO_ROOT / "data" / "app")


def create_app(data_dir: str | Path | None = None) -> FastAPI:
    data_dir = resolve_data_dir(data_dir)
    app = FastAPI(title="scTCNA", version=__version__, description="Single-cell transcriptional noise atlas for cancer: JSON API over the aggregated tables.")
    app.state.data_dir = data_dir
    app.state.store = None
    app.state.data_error = None
    try:
        app.state.store = Store(data_dir)
        log.info("loaded app tables from %s (%d datasets)", data_dir, app.state.store.health()["n_datasets"])
    except FileNotFoundError as exc:
        app.state.data_error = f"app data not available: {exc}"
        log.warning("%s; the API will answer 503 until the tables exist and the server is restarted", app.state.data_error)

    @app.exception_handler(StoreError)
    async def store_error(request: Request, exc: StoreError):
        return json_response({"detail": str(exc)}, status_code=exc.status_code)

    app.include_router(router)
    if STATIC_DIR.is_dir():
        app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="static")
    return app


app = create_app()
