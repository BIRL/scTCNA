#!/usr/bin/env python
"""Measure the web app's peak memory under a worst-case request mix.

Starts ``uvicorn app.main:app`` as a child process on the given data directory,
requests every page and every heavy endpoint for every cancer type (all
comparisons, all top-gene metrics, per-cell violin data, per-patient tables,
200-gene heatmaps for every metric), stops the server and reports the child's
peak resident memory. Used to check that the app fits the 512 MB of Render's
free tier.

Usage:
    python scripts/check_memory.py --data data/app --out docs/memory_check.txt
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import socket
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RANK_METRICS = [
    "log2fc_cv_median", "log2fc_mad_median", "log2fc_fano_median", "delta_cv2_resid", "cv_malignant", "cv_non_malignant",
    "mad_malignant", "mad_non_malignant", "fano_malignant", "fano_non_malignant", "cv2_resid_malignant", "cv2_resid_non_malignant",
]  # fmt: skip
PAGES = ["/", "/summary.html", "/top-genes.html", "/tumor.html", "/heatmap.html", "/about.html"]


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _get(base: str, path: str, params: dict | None = None) -> bytes:
    url = base + path + ("?" + urllib.parse.urlencode(params) if params else "")
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=str(ROOT / "data" / "app"))
    ap.add_argument("--out", default=None, help="also write the report to this file")
    args = ap.parse_args(argv)

    port = _free_port()
    env = {**os.environ, "SCTCNA_DATA_DIR": str(Path(args.data).resolve())}
    proc = subprocess.Popen([sys.executable, "-m", "uvicorn", "app.main:app", "--host", "127.0.0.1", "--port", str(port)], cwd=ROOT, env=env)
    base = f"http://127.0.0.1:{port}"
    for _ in range(300):
        try:
            health = json.loads(_get(base, "/api/health"))
            break
        except OSError:
            time.sleep(0.2)
    else:
        proc.kill()
        raise SystemExit("server did not start")

    t0 = time.time()
    n_requests = 0
    try:
        for page in PAGES:
            _get(base, page)
            n_requests += 1
        cancers = [c["cancer_type"] for c in json.loads(_get(base, "/api/cancer-types"))]
        genes_by_cancer = {}
        for ct in cancers:
            for metric in ("cv", "mad", "fano", "cv2_resid", "sd"):
                _get(base, f"/api/cancer-types/{ct}/comparison", {"metric": metric})
                n_requests += 1
            for metric in RANK_METRICS:
                top = json.loads(_get(base, f"/api/cancer-types/{ct}/top-genes", {"metric": metric, "n": 500}))
                n_requests += 1
            genes_by_cancer[ct] = [r["gene"] for r in top["rows"]]
            for score in ("depth2", "depth"):
                for level in ("cell", "sample"):
                    _get(base, f"/api/cancer-types/{ct}/tumor-heterogeneity", {"score": score, "level": level})
                    n_requests += 1
            _get(base, f"/api/cancer-types/{ct}/noise-contrasts")
            gene = genes_by_cancer[ct][0]
            _get(base, f"/api/genes/{gene}")
            _get(base, f"/api/genes/{gene}/patients", {"cancer_type": ct})
            _get(base, f"/api/genes/{gene}/celltypes", {"cancer_type": ct})
            n_requests += 4
        pool = sorted({g for gs in genes_by_cancer.values() for g in gs})[:200]
        for metric in RANK_METRICS:
            for expressed in ("true", "false"):
                _get(base, "/api/heatmap", {"genes": ",".join(pool), "metric": metric, "expressed_only": expressed})
                n_requests += 1
        _get(base, "/api/tumor-heterogeneity/summary")
        _get(base, "/api/noise-contrasts")
        _get(base, "/api/download/gene_summary.csv")
        n_requests += 3
    finally:
        proc.terminate()
        proc.wait(timeout=30)

    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak_mb = peak / 1e6 if platform.system() == "Darwin" else peak / 1e3  # bytes on macOS, KiB on Linux
    report = (
        f"scTCNA app memory check ({time.strftime('%Y-%m-%d %H:%M', time.localtime())})\n"
        f"data: {args.data} ({health.get('n_datasets')} datasets, {health.get('n_genes')} genes)\n"
        f"requests: {n_requests} (every page, every heavy endpoint for every cancer type) in {time.time() - t0:.0f} s\n"
        f"peak resident memory of the uvicorn process: {peak_mb:.0f} MB (limit on Render free tier: 512 MB)\n"
        f"platform: {platform.system()} {platform.machine()}, Python {platform.python_version()}\n"
    )
    print(report, end="")
    if args.out:
        Path(args.out).write_text(report)
    return 0 if peak_mb < 512 else 1


if __name__ == "__main__":
    sys.exit(main())
