"""Registry and parameter loading."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

DEFAULT_CONFIG = Path(__file__).resolve().parent.parent / "configs" / "datasets.yaml"


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load the YAML registry (datasets, sources, pipeline parameters)."""
    path = Path(path) if path else DEFAULT_CONFIG
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    for key in ("datasets", "params", "sources"):
        if key not in cfg:
            raise ValueError(f"config {path} is missing the '{key}' section")
    return cfg


def enabled_datasets(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    return [d for d in cfg["datasets"] if d.get("enabled", False)]


def get_dataset(cfg: dict[str, Any], dataset_id: str) -> dict[str, Any]:
    for d in cfg["datasets"]:
        if d["id"] == dataset_id:
            return d
    raise KeyError(f"unknown dataset id: {dataset_id}")


def params(cfg: dict[str, Any]) -> dict[str, Any]:
    return cfg["params"]
