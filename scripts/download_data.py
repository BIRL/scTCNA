#!/usr/bin/env python
"""Download scTCNA source datasets listed in configs/datasets.yaml.

3CA bundles (Dropbox) and TISCH2 files (tisch.compbio.cn) are fetched into
``data/raw/<source>/<dataset_id>/``. Each download is recorded in
``data/raw/manifest.json`` with URL, size, SHA-256 and UTC download date.
With ``--convert`` the bundle is turned into ``<dataset_id>.h5ad`` (raw UMI
counts, curated cell annotations). The multi-gigabyte expression matrix is
streamed straight out of the .tar.gz (never extracted) and the archive is
deleted afterwards unless ``--keep-archives``; the pipeline reads the .h5ad.

Usage:
    python scripts/download_data.py --list
    python scripts/download_data.py --datasets LUAD_Kim2020 COAD_Lee2020 --convert
    python scripts/download_data.py --all --convert
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
import tarfile
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sctcna.config import enabled_datasets, get_dataset, load_config, params  # noqa: E402

log = logging.getLogger("download_data")
UA = {"User-Agent": "scTCNA-downloader/0.1 (+https://github.com/)"}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, timeout: int = 60) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, headers=UA, timeout=timeout, allow_redirects=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length") or 0)
        done = 0
        t0 = time.time()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                done += len(chunk)
                if total and done % (50 << 20) < (1 << 20):
                    log.info("  %s: %.0f/%.0f MB (%.1f MB/s)", dest.name, done / 1e6, total / 1e6, done / 1e6 / max(time.time() - t0, 1e-6))
    tmp.rename(dest)
    log.info("  %s: %.1f MB in %.0fs", dest.name, dest.stat().st_size / 1e6, time.time() - t0)
    return dest


def load_manifest(raw_dir: Path) -> dict:
    p = raw_dir / "manifest.json"
    return json.load(open(p)) if p.exists() else {}


def save_manifest(raw_dir: Path, manifest: dict) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    with open(raw_dir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)


def fetch_3ca(entry: dict, ds_dir: Path, keep_archives: bool) -> dict:
    rec = {"files": {}}
    for key, url in (("data", entry["data_url"]), ("meta", entry.get("meta_url"))):
        if not url:
            continue
        name = entry["archive"] if key == "data" else "Meta-" + entry["archive"].replace("Data_", "", 1)
        arc = ds_dir / name
        if not arc.exists():
            log.info("downloading %s %s", entry["id"], key)
            download(url, arc)
        rec["files"][key] = {"url": url, "file": name, "bytes": arc.stat().st_size, "sha256": sha256_of(arc)}
        with tarfile.open(arc) as tf:
            # the expression matrix (gigabytes of text) stays inside the archive and is streamed at conversion time
            members = [m for m in tf.getmembers() if not Path(m.name).name.startswith("Exp_data_")]
            tf.extractall(ds_dir, members=members, filter="data")
        if not keep_archives and key == "meta":
            arc.unlink()
    # 3CA ships the curated cell table in both bundles; prefer the Meta-data copy when present
    meta_dirs = [p for p in ds_dir.iterdir() if p.is_dir() and p.name.startswith("Meta-data")]
    data_dirs = [p for p in ds_dir.iterdir() if p.is_dir() and p.name.startswith("Data_")]
    if meta_dirs and data_dirs:
        for fname in ("Cells.csv", "Samples.csv"):
            src = meta_dirs[0] / fname
            if src.exists():
                shutil.copy2(src, data_dirs[0] / fname)
    return rec


def fetch_tisch2(entry: dict, ds_dir: Path) -> dict:
    rec = {"files": {}}
    for key, url in (("data", entry["data_url"]), ("meta", entry["meta_url"])):
        name = Path(url.split("?")[0]).name
        dest = ds_dir / name
        if not dest.exists():
            log.info("downloading %s %s", entry["id"], key)
            download(url, dest)
        rec["files"][key] = {"url": url, "file": name, "bytes": dest.stat().st_size, "sha256": sha256_of(dest)}
    return rec


def convert(entry: dict, ds_dir: Path, p: dict, delete_raw: bool) -> Path:
    from sctcna.io import load_3ca, load_tisch2

    out = ds_dir / f"{entry['id']}.h5ad"
    if out.exists():
        log.info("%s: %s exists, skipping conversion", entry["id"], out.name)
        return out
    cap = p["sampling"]["max_cells_per_patient_compartment"]
    seed = p["sampling"]["seed"]
    if entry["source"] == "3CA":
        adata = load_3ca(ds_dir, entry["id"], entry["cancer_type"], cap, seed, archive=ds_dir / entry["archive"])
    else:
        h5 = ds_dir / Path(entry["data_url"].split("?")[0]).name
        meta = ds_dir / Path(entry["meta_url"].split("?")[0]).name
        adata = load_tisch2(h5, meta, entry["id"], entry["cancer_type"], cap, seed)
    adata.write_h5ad(out, compression="gzip")
    log.info("%s: wrote %s (%d cells x %d genes)", entry["id"], out.name, adata.n_obs, adata.n_vars)
    if delete_raw:
        for big in [*ds_dir.rglob("Exp_data_*.mtx"), *ds_dir.glob("Data_*.tar.gz")]:
            big.unlink()
            log.info("%s: deleted %s to free disk space", entry["id"], big.name)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--raw", default=str(ROOT / "data" / "raw"))
    ap.add_argument("--datasets", nargs="*", help="dataset ids to fetch")
    ap.add_argument("--all", action="store_true", help="fetch every enabled dataset")
    ap.add_argument("--list", action="store_true", help="print the registry and exit")
    ap.add_argument("--convert", action="store_true", help="convert to .h5ad after download")
    ap.add_argument("--keep-archives", action="store_true", help="keep .tar.gz archives after conversion")
    ap.add_argument("--keep-raw", action="store_true", help="alias of --keep-archives")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    cfg = load_config(args.config)
    if args.list:
        for d in cfg["datasets"]:
            print(f"{'*' if d.get('enabled') else ' '} {d['id']:28s} {d['source']:7s} {d['cancer_type']:6s} {d['study']:28s} {d['data_url']}")
        print("* = enabled by default")
        return 0
    entries = [get_dataset(cfg, d) for d in args.datasets] if args.datasets else (enabled_datasets(cfg) if args.all else [])
    if not entries:
        ap.error("give --datasets ID [ID ...] or --all (or --list)")
    raw_dir = Path(args.raw)
    manifest = load_manifest(raw_dir)
    failures = 0
    for entry in entries:
        ds_dir = raw_dir / entry["source"].lower() / entry["id"]
        ds_dir.mkdir(parents=True, exist_ok=True)
        try:
            if (ds_dir / f"{entry['id']}.h5ad").exists() and entry["id"] in manifest:
                log.info("%s: already converted, skipping", entry["id"])
                continue
            rec = fetch_3ca(entry, ds_dir, args.keep_archives) if entry["source"] == "3CA" else fetch_tisch2(entry, ds_dir)
            rec.update(dataset=entry["id"], source=entry["source"], download_date=time.strftime("%Y-%m-%d", time.gmtime()))
            rec["sha256"] = rec["files"]["data"]["sha256"]
            manifest[entry["id"]] = rec
            save_manifest(raw_dir, manifest)
            if args.convert:
                convert(entry, ds_dir, params(cfg), delete_raw=not (args.keep_raw or args.keep_archives))
        except Exception:
            log.exception("%s: failed", entry["id"])
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
