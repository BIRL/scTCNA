"""Web application for scTCNA: a FastAPI JSON API over the aggregated tables plus the static frontend.

``app.store.Store`` loads the Parquet tables written by ``sctcna aggregate``;
``app.routes`` exposes them as ``/api/...`` endpoints and CSV downloads;
``app.main.create_app`` wires both together (``uvicorn app.main:app``).
"""
