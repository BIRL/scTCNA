#!/usr/bin/env python
"""Print the per-cancer results tables (Markdown) used in README.md from the app tables.

Usage:
    python scripts/summarise_results.py --app data/app
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def _f(x, digits: int = 2) -> str:
    return "n/a" if pd.isna(x) else f"{x:.{digits}f}"


def _p(x) -> str:
    return "n/a" if pd.isna(x) else f"{x:.1e}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--app", default=str(ROOT / "data" / "app"))
    args = ap.parse_args(argv)
    app = Path(args.app)
    ds = pd.read_parquet(app / "datasets.parquet").sort_values("cancer_type")
    nc = pd.read_parquet(app / "noise_contrasts.parquet")
    allref = nc[nc["reference"] == "non_malignant"].set_index("cancer_type")
    tests = pd.read_parquet(app / "heterogeneity_tests.parquet").set_index(["cancer_type", "level"])

    print(
        "| Cancer | Study | Patients (paired) | Cells mal. / non-mal. | Genes tested | q < 0.05 (noisier mal. / non-mal.) | Median log2 CV ratio, raw | Depth-matched | Excess noise, log10 CV² (q) |"
    )
    print("|---|---|---|---|---|---|---|---|---|")
    for _, r in ds.iterrows():
        c = allref.loc[r["cancer_type"]] if r["cancer_type"] in allref.index else None
        print(
            f"| {r['cancer_type']} | {r['study']} | {r['n_patients']} ({r['n_patients_paired']}) | {r['n_cells_malignant']:,} / {r['n_cells_non_malignant']:,} | "
            f"{r['n_genes_expressed']:,} | {r['n_sig_q05']:,} ({r['n_sig_noisier_malignant']:,} / {r['n_sig_noisier_non_malignant']:,}) | "
            f"{_f(c['median_log2fc_cv'] if c is not None else None)} | {_f(c['median_log2fc_cv_ds'] if c is not None else None)} | "
            f"{_f(c['median_excess_cv2_ds'] if c is not None else None, 3)} ({_p(c['q_excess_cv2_ds'] if c is not None else None)}) |"
        )
    print()
    print("| Cancer | DEPTH2 per cell, median diff. mal. - non-mal. (p) | DEPTH2 pseudo-bulk, median diff. (p) | Spearman CV vs downsampled CV |")
    print("|---|---|---|---|")
    for _, r in ds.iterrows():
        cell = tests.loc[(r["cancer_type"], "cell")] if (r["cancer_type"], "cell") in tests.index else None
        pb = tests.loc[(r["cancer_type"], "sample")] if (r["cancer_type"], "sample") in tests.index else None
        print(
            f"| {r['cancer_type']} | {_f(cell['median_diff'] if cell is not None else None, 3)} ({_p(cell['p_value'] if cell is not None else None)}) | "
            f"{_f(pb['median_diff'] if pb is not None else None, 3)} ({_p(pb['p_value'] if pb is not None else None)}) | {_f(r['spearman_cv_vs_downsampled'])} |"
        )
    print()
    lin = nc[(nc["reference"] != "non_malignant") & (nc["n_patients"] >= 6)]
    print("Single-lineage references with >= 6 patients (median excess noise of malignant cells, log10 CV²; q):")
    for ct, g in lin.groupby("cancer_type"):
        parts = [f"{r.reference} {r.median_excess_cv2_ds:+.3f} ({_p(r.q_excess_cv2_ds)})" for r in g.sort_values("reference").itertuples()]
        print(f"- {ct}: " + "; ".join(parts))
    return 0


if __name__ == "__main__":
    sys.exit(main())
