"""Statistics with the patient as the unit of replication.

Every per-gene metric is first computed within one patient x compartment group
of cells. Cancer-type values are medians across patients, and malignant vs
non-malignant comparisons are paired across patients (Wilcoxon signed-rank on
the per-patient log2 CV ratio, which is TCNA Eq. 1 evaluated per patient).
Multiple testing is controlled per cancer type with Benjamini-Hochberg.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps
from statsmodels.stats.multitest import multipletests

SUMMARY_METRICS = ["mean", "detection_rate", "sd", "cv", "cv2", "mad", "fano", "cv2_resid"]


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN p-values stay NaN and do not count towards m."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q


def rank_biserial(d: np.ndarray) -> float:
    """Matched-pairs rank-biserial correlation of paired differences (Kerby 2014).

    +1 when every non-zero difference is positive, -1 when every one is negative.
    """
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d) & (d != 0)]
    if d.size == 0:
        return np.nan
    r = sps.rankdata(np.abs(d))
    pos = r[d > 0].sum()
    neg = r[d < 0].sum()
    return float((pos - neg) / (pos + neg))


def compartment_medians(metrics_long: pd.DataFrame, metrics: list[str] | None = None, min_patients: int = 1) -> pd.DataFrame:
    """Median of each metric across patients, per (dataset, gene, compartment); wide layout.

    Columns are ``<metric>_<compartment>`` plus ``n_patients_<compartment>``.
    """
    metrics = metrics or [m for m in SUMMARY_METRICS if m in metrics_long.columns]
    g = metrics_long.groupby(["dataset", "gene", "compartment"], observed=True)
    med = g[metrics].median()
    med["n_patients"] = g["patient"].nunique()
    wide = med.unstack("compartment")
    wide.columns = [f"{m}_{c}" for m, c in wide.columns]
    wide = wide.reset_index()
    return wide


def paired_noise_tests(
    metrics_long: pd.DataFrame,
    metric: str = "cv",
    min_patients: int = 6,
    treat: str = "malignant",
    control: str = "non_malignant",
    difference: bool = False,
) -> pd.DataFrame:
    """Per (dataset, gene) paired malignant vs non-malignant test across patients.

    The per-patient contrast is log2(metric_treat / metric_control) (TCNA Eq. 1
    per patient when metric = 'cv'), or the plain difference when
    ``difference`` is set (for metrics already on a log scale, such as the
    LOWESS residual). Returns n_pairs, the median and mean per-patient contrast
    (``log2fc_<metric>_*`` or ``delta_<metric>_*``), the Wilcoxon signed-rank
    statistic and p-value, the matched-pairs rank-biserial effect size, and BH
    q-values per dataset.
    """
    prefix = f"delta_{metric}" if difference else f"log2fc_{metric}"
    sub = metrics_long[metrics_long["compartment"].isin([treat, control])]
    wide = sub.pivot_table(index=["dataset", "gene", "patient"], columns="compartment", values=metric, aggfunc="first", observed=True)
    if treat not in wide.columns or control not in wide.columns:
        return pd.DataFrame(columns=["dataset", "gene", "n_pairs", f"{prefix}_median", f"{prefix}_mean", "wilcoxon_stat", "p_value", "rank_biserial", "q_value"])
    a = wide[treat].to_numpy(dtype=float)
    b = wide[control].to_numpy(dtype=float)
    d = np.full(a.shape, np.nan)
    if difference:
        ok = np.isfinite(a) & np.isfinite(b)
        d[ok] = a[ok] - b[ok]
    else:
        ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
        d[ok] = np.log2(a[ok] / b[ok])
    pairs = pd.DataFrame({"d": d}, index=wide.index).dropna().reset_index()
    rows = []
    for (ds, gene), grp in pairs.groupby(["dataset", "gene"], observed=True, sort=False):
        dd = grp["d"].to_numpy()
        n = dd.size
        rec = {"dataset": ds, "gene": gene, "n_pairs": n, f"{prefix}_median": float(np.median(dd)), f"{prefix}_mean": float(dd.mean())}
        if n >= min_patients:
            nz = dd[dd != 0]
            if nz.size == 0:
                rec.update(wilcoxon_stat=0.0, p_value=1.0)
            else:
                res = sps.wilcoxon(nz, zero_method="wilcox", alternative="two-sided")
                rec.update(wilcoxon_stat=float(res.statistic), p_value=float(res.pvalue))
            rec["rank_biserial"] = rank_biserial(dd)
        else:
            rec.update(wilcoxon_stat=np.nan, p_value=np.nan, rank_biserial=np.nan)
        rows.append(rec)
    out = pd.DataFrame(rows)
    if out.empty:
        out["q_value"] = []
        return out
    out["q_value"] = np.nan
    for _ds, idx in out.groupby("dataset", observed=True).indices.items():
        out.loc[out.index[idx], "q_value"] = bh_adjust(out.loc[out.index[idx], "p_value"].to_numpy())
    return out


def paired_score_test(sample_scores: pd.DataFrame, score: str, treat: str = "malignant", control: str = "non_malignant", min_patients: int = 6) -> dict:
    """Paired Wilcoxon of a per-patient score (e.g. DEPTH2) between compartments."""
    wide = sample_scores.pivot_table(index="patient", columns="compartment", values=score, aggfunc="first", observed=True)
    if treat not in wide.columns or control not in wide.columns:
        return {"n_pairs": 0, "median_diff": np.nan, "p_value": np.nan, "rank_biserial": np.nan}
    pairs = wide[[treat, control]].dropna()
    d = (pairs[treat] - pairs[control]).to_numpy(dtype=float)
    rec = {"n_pairs": int(d.size), "median_diff": float(np.median(d)) if d.size else np.nan, "p_value": np.nan, "rank_biserial": np.nan}
    if d.size >= min_patients:
        nz = d[d != 0]
        rec["p_value"] = 1.0 if nz.size == 0 else float(sps.wilcoxon(nz, zero_method="wilcox").pvalue)
        rec["rank_biserial"] = rank_biserial(d)
    return rec


def one_sample_wilcoxon(values, min_n: int = 6) -> dict:
    """Wilcoxon signed-rank test of per-patient values against 0 (two-sided) with the rank-biserial effect size."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    rec = {
        "n": int(v.size),
        "median": float(np.median(v)) if v.size else np.nan,
        "frac_positive": float((v > 0).mean()) if v.size else np.nan,
        "p_value": np.nan,
        "rank_biserial": np.nan,
    }
    if v.size >= min_n:
        nz = v[v != 0]
        rec["p_value"] = 1.0 if nz.size == 0 else float(sps.wilcoxon(nz, zero_method="wilcox").pvalue)
        rec["rank_biserial"] = rank_biserial(v)
    return rec
