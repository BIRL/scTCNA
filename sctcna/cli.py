"""Command line interface: ``python -m sctcna <command>``."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from . import __version__
from .config import enabled_datasets, get_dataset, load_config, params
from .io import load_dataset
from .pipeline import aggregate, process_dataset
from .synthetic import make_synthetic_adata


def _run(args) -> int:
    cfg = load_config(args.config)
    p = params(cfg)
    if args.synthetic:
        adata = make_synthetic_adata(n_patients=args.synthetic_patients, cells_per_group=args.synthetic_cells, seed=p["sampling"]["seed"])
        meta = {
            "id": "SYN_Test",
            "cancer_type": "SYN",
            "cancer_name": "Synthetic test data",
            "source": "synthetic",
            "study": "synthetic",
            "citation": "n/a",
            "technology": "simulated",
            "data_url": "",
            "enabled": True,
        }
        process_dataset(adata, p, meta, args.out)
        return 0
    entries = [get_dataset(cfg, d) for d in args.datasets] if args.datasets else enabled_datasets(cfg)
    failures = 0
    for entry in entries:
        try:
            # passed straight through so the loaded matrix is not kept alive by this frame
            process_dataset(load_dataset(entry, args.raw, p["sampling"]["max_cells_per_patient_compartment"], p["sampling"]["seed"]), p, entry, args.out)
        except FileNotFoundError as exc:
            logging.error("%s: data not found (%s). Run scripts/download_data.py first.", entry["id"], exc)
            failures += 1
        except Exception:
            logging.exception("%s: processing failed", entry["id"])
            failures += 1
            if args.fail_fast:
                raise
    return 1 if failures else 0


def _aggregate(args) -> int:
    cfg = load_config(args.config)
    manifest = {}
    mpath = Path(args.raw) / "manifest.json"
    if mpath.exists():
        manifest = json.load(open(mpath))
    if args.include_synthetic:
        cfg = dict(cfg)
        cfg["datasets"] = list(cfg["datasets"]) + [
            {
                "id": "SYN_Test",
                "cancer_type": "SYN",
                "cancer_name": "Synthetic test data",
                "source": "synthetic",
                "study": "synthetic",
                "citation": "n/a",
                "technology": "simulated",
                "data_url": "",
            }
        ]
    res = aggregate(args.processed, args.app, cfg, manifest)
    print(json.dumps(res))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="sctcna", description="scTCNA pipeline")
    ap.add_argument("--version", action="version", version=__version__)
    ap.add_argument("-v", "--verbose", action="store_true")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="process datasets into Parquet tables")
    r.add_argument("--config", default=None)
    r.add_argument("--raw", default="data/raw", help="directory with downloaded datasets")
    r.add_argument("--out", default="data/processed")
    r.add_argument("--datasets", nargs="*", help="dataset ids (default: all enabled)")
    r.add_argument("--synthetic", action="store_true", help="process a synthetic dataset instead")
    r.add_argument("--synthetic-patients", type=int, default=6)
    r.add_argument("--synthetic-cells", type=int, default=300)
    r.add_argument("--fail-fast", action="store_true")
    r.set_defaults(func=_run)

    a = sub.add_parser("aggregate", help="build the app tables from processed datasets")
    a.add_argument("--config", default=None)
    a.add_argument("--raw", default="data/raw")
    a.add_argument("--processed", default="data/processed")
    a.add_argument("--app", default="data/app")
    a.add_argument("--include-synthetic", action="store_true")
    a.set_defaults(func=_aggregate)

    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
