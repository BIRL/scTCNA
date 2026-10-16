"""scTCNA's standalone web server: the JSON API at ``/api`` plus the static frontend.

The API itself lives in ``sctcna.web`` so that other applications (such as the
TCNA portal, which serves it at ``/api/sc``) can mount it. ``create_app`` loads
the tables from ``data_dir`` (argument, else ``SCTCNA_DATA_DIR``, else
``<repo>/data/app``). A missing data directory does not break start-up: the API
then answers 503 on every route, including ``/api/health``.

Run with ``uvicorn app.main:app``.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from sctcna import __version__
from sctcna.web import create_api, load_store, resolve_data_dir

STATIC_DIR = Path(__file__).resolve().parent / "static"

__all__ = ["create_app", "resolve_data_dir"]


def create_app(data_dir: str | Path | None = None) -> FastAPI:
    data_dir = resolve_data_dir(data_dir)
    store, data_error = load_store(data_dir)
    app = FastAPI(title="scTCNA", version=__version__)
    app.state.data_dir = data_dir
    app.state.store = store
    app.state.data_error = data_error
    app.mount("/api", create_api(store=store, data_error=data_error), name="api")
    if STATIC_DIR.is_dir():
        app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="static")
    return app


app = create_app()
