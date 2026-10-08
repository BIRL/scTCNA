"""In-memory access to the aggregated scTCNA app tables.

A ``Store`` reads every Parquet table of an app data directory once (the
tables are written by ``sctcna.pipeline.aggregate``), keeps a few lookup
structures around and answers API questions with plain Python values.
Request-time work is limited to filtering and sorting; the heavier views
(compartment comparison, heatmap, tumour heterogeneity) are memoised per
parameter combination because the tables never change while serving.

Validation problems are raised as ``NotFound`` / ``BadRequest`` so the web
layer can map them to HTTP status codes without the store importing FastAPI.
Returned structures may still contain NaN floats; the web layer turns those
into JSON null.
"""

from __future__ import annotations

import json
import re
from bisect import bisect_left
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

MALIGNANT = "malignant"
NON_MALIGNANT = "non_malignant"
COMPARTMENTS = (MALIGNANT, NON_MALIGNANT)

# Metrics that rank genes (top-genes, heatmap) and the q-value column that goes with each; None means the
# metric is a plain per-compartment summary without a paired test.
RANK_METRICS: dict[str, str | None] = {
    "log2fc_cv_median": "q_value",
    "log2fc_mad_median": "q_value_mad",
    "log2fc_fano_median": "q_value_fano",
    "delta_cv2_resid": "q_value",
    "cv_malignant": None,
    "cv_non_malignant": None,
    "mad_malignant": None,
    "mad_non_malignant": None,
    "fano_malignant": None,
    "fano_non_malignant": None,
    "cv2_resid_malignant": None,
    "cv2_resid_non_malignant": None,
}
# Metrics of the malignant vs non-malignant comparison: (per-gene fold-change column or None when it is
# computed from the two compartment medians, q-value column).
COMPARISON_METRICS: dict[str, tuple[str | None, str]] = {
    "cv": ("log2fc_cv_median", "q_value"),
    "mad": ("log2fc_mad_median", "q_value_mad"),
    "fano": ("log2fc_fano_median", "q_value_fano"),
    "cv2_resid": ("delta_cv2_resid", "q_value"),
    "sd": (None, "q_value"),
}
DIRECTIONS = ("noisy", "stable")
SCORES = ("depth2", "depth")
LEVELS = ("cell", "sample")
MAX_HEATMAP_GENES = 200
MAX_PATIENT_TABLES = 3
MAX_CELL_VALUES = 5000
SIGNIFICANCE = 0.05

CANCER_TYPE_COLUMNS = [
    "cancer_type",
    "cancer_name",
    "tcna_site",
    "dataset",
    "study",
    "citation",
    "technology",
    "site",
    "n_patients",
    "n_patients_paired",
    "n_cells_malignant",
    "n_cells_non_malignant",
    "n_genes",
    "n_genes_expressed",
]
TOP_GENE_COLUMNS = [
    "gene",
    "mean_malignant",
    "mean_non_malignant",
    "detection_rate_malignant",
    "detection_rate_non_malignant",
    "cv_malignant",
    "cv_non_malignant",
    "mad_malignant",
    "mad_non_malignant",
    "fano_malignant",
    "fano_non_malignant",
    "cv2_resid_malignant",
    "cv2_resid_non_malignant",
    "log2fc_cv_median",
    "delta_cv2_resid",
    "p_value",
    "q_value",
    "rank_biserial",
    "n_pairs",
]
PATIENT_ROW_COLUMNS = ["patient", "compartment", "n_cells", "mean", "detection_rate", "cv", "mad", "fano", "cv2_resid", "cv_ds", "cv2_resid_ds"]
CELLTYPE_ROW_COLUMNS = ["cell_type", "n_patients", "mean", "detection_rate", "cv", "mad", "fano", "cv2_resid"]
SAMPLE_ROW_COLUMNS = ["patient", "compartment", "n_cells", "depth2_pseudobulk", "depth_pseudobulk", "depth_cell_median", "depth2_cell_median"]
CELL_ROW_COLUMNS = ["dataset", "cancer_type", "patient", "cell_type", "compartment", "n_counts", "n_genes", "pct_mito", "depth", "depth2"]


class StoreError(Exception):
    """A request the store cannot answer; ``status_code`` is the HTTP status to use."""

    status_code = 400


class NotFound(StoreError):
    status_code = 404


class BadRequest(StoreError):
    status_code = 400


def _records(df: pd.DataFrame) -> list[dict[str, Any]]:
    df = df.copy()
    for c in df.columns:
        if isinstance(df[c].dtype, pd.CategoricalDtype):
            df[c] = df[c].astype(object)
    return df.to_dict("records")


def _log2_ratio(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """log2(a / b), NaN where either side is not a positive finite number (TCNA Eq. 1)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    out = np.full(a.shape, np.nan)
    ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    out[ok] = np.log2(a[ok] / b[ok])
    return out


def _floats(values: np.ndarray) -> list[float]:
    return np.asarray(values, dtype=float).tolist()


def _nanmedian(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else float("nan")


class Store:
    """All app tables of one data directory, loaded at construction."""

    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)
        if not (self.data_dir / "datasets.parquet").exists():
            raise FileNotFoundError(f"no app tables found in {self.data_dir} (run `python -m sctcna aggregate` first)")
        self.datasets = pd.read_parquet(self.data_dir / "datasets.parquet")
        self.gene_summary = pd.read_parquet(self.data_dir / "gene_summary.parquet")
        self.cell_scores = pd.read_parquet(self.data_dir / "cell_scores.parquet")
        self.sample_scores = pd.read_parquet(self.data_dir / "sample_scores.parquet")
        self.celltype_summary = pd.read_parquet(self.data_dir / "celltype_summary.parquet")
        build_path = self.data_dir / "build.json"
        self.build = json.load(open(build_path)) if build_path.exists() else {}

        # one dataset per cancer type in this release; keep the first row should the table ever hold two
        self.datasets = self.datasets.sort_values(["cancer_type", "dataset"]).reset_index(drop=True)
        self._cancers: dict[str, dict[str, Any]] = {r["cancer_type"]: r for r in _records(self.datasets.drop_duplicates("cancer_type"))}
        self.cancer_types = sorted(self._cancers)

        # gene symbols: case-insensitive lookup of the canonical spelling, plus a sorted key list for prefix search
        self._canonical: dict[str, str] = {}
        for gene in self.gene_summary["gene"].drop_duplicates():
            self._canonical.setdefault(gene.upper(), gene)
        self._search_keys = sorted(self._canonical)
        self._n_cancer_types = {g: int(n) for g, n in self.gene_summary.groupby("gene")["cancer_type"].nunique().items()}

        # row lookups into the summary tables: positional indices rather than per-cancer copies keep memory low
        self._summary_by_gene = self.gene_summary.groupby("gene").indices
        self._summary_idx_by_cancer = {str(ct): idx for ct, idx in self.gene_summary.groupby("cancer_type", observed=True).indices.items()}
        self._cells_idx_by_cancer = {str(ct): idx for ct, idx in self.cell_scores.groupby("cancer_type", observed=True).indices.items()}
        self._samples_by_cancer = {str(ct): df.reset_index(drop=True) for ct, df in self.sample_scores.groupby("cancer_type", observed=True)}
        ct_key = self.celltype_summary["cancer_type"].astype(str) + "\t" + self.celltype_summary["gene"].astype(str)
        order = np.argsort(ct_key.to_numpy(dtype=object), kind="stable")
        self.celltype_summary = self.celltype_summary.iloc[order].reset_index(drop=True)
        self._celltype_keys = ct_key.to_numpy(dtype=object)[order]

        tests_path = self.data_dir / "heterogeneity_tests.parquet"
        tests = pd.read_parquet(tests_path) if tests_path.exists() else pd.DataFrame(columns=["cancer_type", "score", "level"])
        self._heterogeneity_tests = {(str(r["cancer_type"]), r["score"], r["level"]): r for r in tests.to_dict("records")}
        contrasts_path = self.data_dir / "noise_contrasts.parquet"
        self.noise_contrasts = pd.read_parquet(contrasts_path) if contrasts_path.exists() else pd.DataFrame(columns=["dataset", "cancer_type", "reference"])
        contrasts_patient_path = self.data_dir / "noise_contrasts_patient.parquet"
        self.noise_contrasts_patient = (
            pd.read_parquet(contrasts_patient_path) if contrasts_patient_path.exists() else pd.DataFrame(columns=["dataset", "cancer_type", "patient", "reference"])
        )

        self._patient_tables: dict[str, tuple[pd.DataFrame, dict[str, np.ndarray]]] = {}
        self._cache: dict[tuple, Any] = {}

    # ------------------------------------------------------------------ helpers
    def _memo(self, key: tuple, compute: Callable[[], Any], max_entries: int = 512) -> Any:
        """Memoise ``compute()`` under ``key``; the oldest entry goes when the cache is full."""
        if key not in self._cache:
            if len(self._cache) >= max_entries:
                del self._cache[next(iter(self._cache))]
            self._cache[key] = compute()
        return self._cache[key]

    def cancer(self, cancer_type: str) -> dict[str, Any]:
        """Dataset row of a cancer type."""
        try:
            return self._cancers[cancer_type]
        except KeyError:
            raise NotFound(f"unknown cancer type: {cancer_type}") from None

    def gene_symbol(self, gene: str) -> str:
        """Canonical spelling of a gene symbol (case-insensitive)."""
        try:
            return self._canonical[gene.strip().upper()]
        except KeyError:
            raise NotFound(f"gene not found in any dataset: {gene}") from None

    def _patient_table(self, dataset: str) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
        """Per-patient table of one dataset, read on first use; at most MAX_PATIENT_TABLES stay in memory."""
        if dataset not in self._patient_tables:
            path = self.data_dir / "patient" / f"{dataset}.parquet"
            df = pd.read_parquet(path) if path.exists() else pd.DataFrame(columns=["gene", *PATIENT_ROW_COLUMNS])
            if len(self._patient_tables) >= MAX_PATIENT_TABLES:
                del self._patient_tables[next(iter(self._patient_tables))]
            self._patient_tables[dataset] = (df, {str(g): idx for g, idx in df.groupby("gene", observed=True).indices.items()})
        return self._patient_tables[dataset]

    def _cancer_summary(self, cancer_type: str, expressed_only: bool) -> pd.DataFrame:
        self.cancer(cancer_type)
        idx = self._summary_idx_by_cancer.get(cancer_type)
        if idx is None:
            return self.gene_summary.iloc[:0]
        df = self.gene_summary.iloc[idx]
        if expressed_only:
            df = df[df["expressed"].fillna(False).to_numpy(dtype=bool)]
        return df

    @staticmethod
    def _check_rank_metric(metric: str) -> None:
        if metric not in RANK_METRICS:
            raise BadRequest(f"unknown metric: {metric}; expected one of {', '.join(RANK_METRICS)}")

    # ------------------------------------------------------------------ overview
    def health(self) -> dict[str, Any]:
        return {"n_datasets": int(self.datasets["dataset"].nunique()), "n_genes": len(self._canonical), "build": self.build}

    def cancer_type_rows(self) -> list[dict[str, Any]]:
        return [{c: self._cancers[ct].get(c) for c in CANCER_TYPE_COLUMNS} for ct in self.cancer_types]

    def dataset_rows(self) -> list[dict[str, Any]]:
        return _records(self.datasets)

    # ------------------------------------------------------------------ genes
    def search_genes(self, q: str, limit: int = 20) -> list[dict[str, Any]]:
        """Case-insensitive prefix matches; the exact symbol sorts first, then alphabetically."""
        key = q.strip().upper()
        if not key:
            raise BadRequest("q must not be empty")
        # the exact key is the smallest string with that prefix, so the sorted order already puts it first
        i = bisect_left(self._search_keys, key)
        results = []
        while i < len(self._search_keys) and len(results) < limit and self._search_keys[i].startswith(key):
            gene = self._canonical[self._search_keys[i]]
            results.append({"gene": gene, "n_cancer_types": self._n_cancer_types.get(gene, 0)})
            i += 1
        return results

    def gene(self, gene: str) -> dict[str, Any]:
        symbol = self.gene_symbol(gene)
        rows = self.gene_summary.iloc[self._summary_by_gene[symbol]].sort_values("cancer_type")
        cancers = []
        for rec in _records(rows):
            rec["cancer_name"] = self._cancers.get(rec["cancer_type"], {}).get("cancer_name")
            cancers.append(rec)
        return {"gene": symbol, "cancers": cancers}

    def gene_patients(self, gene: str, cancer_type: str) -> dict[str, Any]:
        symbol = self.gene_symbol(gene)
        dataset = self.cancer(cancer_type)["dataset"]
        table, by_gene = self._patient_table(dataset)
        idx = by_gene.get(symbol)
        if idx is None or len(idx) == 0:
            raise NotFound(f"no per-patient values for {symbol} in {cancer_type} (per-patient tables only hold genes expressed in both compartments)")
        rows = table.iloc[idx].sort_values(["patient", "compartment"])[PATIENT_ROW_COLUMNS]
        return {"gene": symbol, "cancer_type": cancer_type, "dataset": dataset, "rows": _records(rows)}

    def gene_celltypes(self, gene: str, cancer_type: str) -> dict[str, Any]:
        symbol = self.gene_symbol(gene)
        self.cancer(cancer_type)
        key = f"{cancer_type}\t{symbol}"
        lo = int(np.searchsorted(self._celltype_keys, key, side="left"))
        hi = int(np.searchsorted(self._celltype_keys, key, side="right"))
        if hi <= lo:
            raise NotFound(f"no cell-type values for {symbol} in {cancer_type}")
        rows = self.celltype_summary.iloc[lo:hi].sort_values("cell_type")[CELLTYPE_ROW_COLUMNS]
        return {"gene": symbol, "cancer_type": cancer_type, "rows": _records(rows)}

    # ------------------------------------------------------------------ per cancer type
    def top_genes(self, cancer_type: str, metric: str, direction: str = "noisy", n: int = 50, expressed_only: bool = True) -> dict[str, Any]:
        self._check_rank_metric(metric)
        if direction not in DIRECTIONS:
            raise BadRequest(f"unknown direction: {direction}; expected noisy or stable")
        df = self._cancer_summary(cancer_type, expressed_only)
        df = df[np.isfinite(df[metric].to_numpy(dtype=float))]
        df = df.sort_values([metric, "gene"], ascending=[direction == "stable", True]).head(n)
        rows = df[TOP_GENE_COLUMNS].copy()
        rows.insert(1, "value", df[metric].to_numpy(dtype=float))
        return {"cancer_type": cancer_type, "metric": metric, "direction": direction, "n": n, "rows": _records(rows)}

    def comparison(self, cancer_type: str, metric: str = "cv", expressed_only: bool = True) -> dict[str, Any]:
        """Per-gene malignant vs non-malignant values of one noise metric (genes with both values only)."""
        if metric not in COMPARISON_METRICS:
            raise BadRequest(f"unknown metric: {metric}; expected one of {', '.join(COMPARISON_METRICS)}")
        self.cancer(cancer_type)
        return self._memo(("comparison", cancer_type, metric, expressed_only), lambda: self._comparison(cancer_type, metric, expressed_only))

    def _comparison(self, cancer_type: str, metric: str, expressed_only: bool) -> dict[str, Any]:
        fc_column, q_column = COMPARISON_METRICS[metric]
        df = self._cancer_summary(cancer_type, expressed_only)
        a = df[f"{metric}_{MALIGNANT}"].to_numpy(dtype=float)
        b = df[f"{metric}_{NON_MALIGNANT}"].to_numpy(dtype=float)
        keep = np.isfinite(a) & np.isfinite(b)
        df, a, b = df[keep], a[keep], b[keep]
        log2fc = _log2_ratio(a, b) if fc_column is None else df[fc_column].to_numpy(dtype=float)
        q = df[q_column].to_numpy(dtype=float)
        with np.errstate(invalid="ignore"):
            sig = q < SIGNIFICANCE
        return {
            "cancer_type": cancer_type,
            "metric": metric,
            "n_genes": int(len(df)),
            "genes": df["gene"].tolist(),
            "malignant": _floats(a),
            "non_malignant": _floats(b),
            "log2fc": _floats(log2fc),
            "q_value": _floats(q),
            "summary": {
                "median_log2fc": _nanmedian(log2fc),
                "n_sig_q05": int(sig.sum()),
                "n_sig_noisier_malignant": int((sig & (log2fc > 0)).sum()),
                "n_sig_noisier_non_malignant": int((sig & (log2fc < 0)).sum()),
                "n_patients_paired": int(self.cancer(cancer_type)["n_patients_paired"]),
            },
        }

    # ------------------------------------------------------------------ heatmap
    @staticmethod
    def parse_gene_list(text: str) -> list[str]:
        """Split a user-supplied gene list on commas/whitespace, dropping case-insensitive duplicates."""
        seen: set[str] = set()
        genes = []
        for token in re.split(r"[,\s]+", text or ""):
            if token and token.upper() not in seen:
                seen.add(token.upper())
                genes.append(token)
        return genes

    def _pivot(self, metric: str, expressed_only: bool) -> tuple[pd.DataFrame, pd.DataFrame | None]:
        """Gene x cancer-type matrices of a metric and of its q-value (None when the metric has no test)."""

        def build():
            df = self.gene_summary
            if expressed_only:
                df = df[df["expressed"].fillna(False).to_numpy(dtype=bool)]
            df = df.drop_duplicates(["gene", "cancer_type"])
            values = df.pivot(index="gene", columns="cancer_type", values=metric).reindex(columns=self.cancer_types)
            q_column = RANK_METRICS[metric]
            q_values = df.pivot(index="gene", columns="cancer_type", values=q_column).reindex(columns=self.cancer_types) if q_column else None
            return values, q_values

        return self._memo(("pivot", metric, expressed_only), build)

    def heatmap(self, genes: str, metric: str, expressed_only: bool = False) -> dict[str, Any]:
        self._check_rank_metric(metric)
        requested = self.parse_gene_list(genes)
        if not requested:
            raise BadRequest("genes must list at least one gene symbol")
        if len(requested) > MAX_HEATMAP_GENES:
            raise BadRequest(f"at most {MAX_HEATMAP_GENES} genes per heatmap (got {len(requested)})")
        found, missing = [], []
        for token in requested:
            symbol = self._canonical.get(token.upper())
            (found if symbol else missing).append(symbol or token)
        if not found:
            raise BadRequest("none of the requested genes is present in any dataset")

        def build():
            values, q_values = self._pivot(metric, expressed_only)
            return {
                "metric": metric,
                "genes": found,
                "missing": missing,
                "cancer_types": self.cancer_types,
                "values": values.reindex(found).to_numpy(dtype=float).tolist(),
                "q_values": q_values.reindex(found).to_numpy(dtype=float).tolist() if q_values is not None else None,
            }

        return self._memo(("heatmap", tuple(found), tuple(missing), metric, expressed_only), build)

    # ------------------------------------------------------------------ tumour heterogeneity
    def tumor_heterogeneity(self, cancer_type: str, score: str = "depth2", level: str = "cell") -> dict[str, Any]:
        if score not in SCORES:
            raise BadRequest(f"unknown score: {score}; expected depth2 or depth")
        if level not in LEVELS:
            raise BadRequest(f"unknown level: {level}; expected cell or sample")
        self.cancer(cancer_type)
        return self._memo(("heterogeneity", cancer_type, score, level), lambda: self._tumor_heterogeneity(cancer_type, score, level))

    def _tumor_heterogeneity(self, cancer_type: str, score: str, level: str) -> dict[str, Any]:
        samples = self._samples_by_cancer.get(cancer_type, self.sample_scores.iloc[:0])
        out: dict[str, Any] = {"cancer_type": cancer_type, "score": score, "level": level}
        if level == "sample":
            rows = samples.sort_values(["patient", "compartment"])[SAMPLE_ROW_COLUMNS]
            out["rows"] = _records(rows)
            out["paired_test"] = self._paired_test(cancer_type, score, level)
            return out
        cells = self.cell_scores.iloc[self._cells_idx_by_cancer.get(cancer_type, [])]
        compartments = {}
        for comp in COMPARTMENTS:
            values = cells.loc[cells["compartment"] == comp, score].to_numpy(dtype=float)
            shown = values
            if values.size > MAX_CELL_VALUES:
                # fixed seed: the same cells are shown on every request
                shown = values[np.sort(np.random.default_rng(0).choice(values.size, MAX_CELL_VALUES, replace=False))]
            compartments[comp] = {"values": _floats(shown), "n": int(values.size), "median": _nanmedian(values)}
        out["compartments"] = compartments
        patients = samples.sort_values(["patient", "compartment"])[["patient", "compartment", "n_cells", f"{score}_cell_median"]]
        out["patients"] = _records(patients.rename(columns={f"{score}_cell_median": "median"}))
        out["paired_test"] = self._paired_test(cancer_type, score, level)
        return out

    def _paired_test(self, cancer_type: str, score: str, level: str) -> dict[str, Any]:
        """Pre-computed paired Wilcoxon test (written by ``sctcna aggregate``) of malignant vs non-malignant per-patient scores."""
        rec = self._heterogeneity_tests.get((cancer_type, score, level))
        if rec is None:
            return {"n_pairs": 0, "median_diff": float("nan"), "p_value": float("nan"), "rank_biserial": float("nan")}
        return {k: rec[k] for k in ("n_pairs", "median_diff", "p_value", "rank_biserial")}

    # ------------------------------------------------------------------ noise contrasts
    def contrast_summary(self, cancer_type: str | None = None) -> dict[str, Any]:
        """Malignant vs reference-group noise contrasts (depth-matched, mean-corrected), one row per cancer type x reference."""
        df = self.noise_contrasts
        if cancer_type is not None:
            self.cancer(cancer_type)
            df = df[df["cancer_type"] == cancer_type]
        df = df.sort_values(["cancer_type", "reference"])
        return {"rows": _records(df)}

    def contrast_patients(self, cancer_type: str) -> dict[str, Any]:
        self.cancer(cancer_type)
        df = self.noise_contrasts_patient
        df = df[df["cancer_type"] == cancer_type].sort_values(["reference", "patient"])
        return {"cancer_type": cancer_type, "summary": self.contrast_summary(cancer_type)["rows"], "rows": _records(df)}

    def heterogeneity_summary(self) -> dict[str, Any]:
        """One row per cancer type: the pre-computed paired tests plus the median across patients of the per-patient DEPTH2 cell medians."""

        def build():
            rows = []
            for ct in self.cancer_types:
                info = self._cancers[ct]
                samples = self._samples_by_cancer.get(ct, self.sample_scores.iloc[:0])
                per_comp = {comp: _nanmedian(samples.loc[samples["compartment"] == comp, "depth2_cell_median"].to_numpy(dtype=float)) for comp in COMPARTMENTS}
                rows.append(
                    {
                        "cancer_type": ct,
                        "cancer_name": info.get("cancer_name"),
                        "n_patients_paired": info.get("n_patients_paired"),
                        "depth2_cell_median_malignant": per_comp[MALIGNANT],
                        "depth2_cell_median_non_malignant": per_comp[NON_MALIGNANT],
                        **{
                            c: info.get(c)
                            for c in ("depth2_cell_median_diff", "depth2_cell_p", "depth2_pseudobulk_median_diff", "depth2_pseudobulk_p", "depth_cell_median_diff", "depth_cell_p")
                        },
                    }
                )
            return {"rows": rows}

        return self._memo(("heterogeneity_summary",), build)

    # ------------------------------------------------------------------ frames for CSV downloads
    def gene_summary_frame(self, cancer_type: str | None = None) -> pd.DataFrame:
        return self._cancer_summary(cancer_type, expressed_only=False) if cancer_type else self.gene_summary

    def sample_scores_frame(self, cancer_type: str | None = None) -> pd.DataFrame:
        if cancer_type is None:
            return self.sample_scores
        self.cancer(cancer_type)
        return self._samples_by_cancer.get(cancer_type, self.sample_scores.iloc[:0])

    def cell_scores_frame(self, cancer_type: str) -> pd.DataFrame:
        self.cancer(cancer_type)
        return self.cell_scores.iloc[self._cells_idx_by_cancer.get(cancer_type, [])][CELL_ROW_COLUMNS]

    def top_genes_frame(self, cancer_type: str, metric: str, direction: str, n: int, expressed_only: bool) -> pd.DataFrame:
        return pd.DataFrame(self.top_genes(cancer_type, metric, direction, n, expressed_only)["rows"], columns=["gene", "value", *TOP_GENE_COLUMNS[1:]])

    def gene_patients_frame(self, gene: str, cancer_type: str) -> pd.DataFrame:
        res = self.gene_patients(gene, cancer_type)
        df = pd.DataFrame(res["rows"], columns=PATIENT_ROW_COLUMNS)
        df.insert(0, "gene", res["gene"])
        df.insert(1, "cancer_type", cancer_type)
        return df

    def heatmap_frame(self, genes: str, metric: str, expressed_only: bool) -> pd.DataFrame:
        """Wide gene x cancer-type table; ``<cancer_type>_q`` columns are added when the metric has a q-value."""
        res = self.heatmap(genes, metric, expressed_only)
        df = pd.DataFrame(res["values"], columns=res["cancer_types"])
        if res["q_values"] is not None:
            for j, ct in enumerate(res["cancer_types"]):
                df[f"{ct}_q"] = [row[j] for row in res["q_values"]]
        df.insert(0, "gene", res["genes"])
        return df
