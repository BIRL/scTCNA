"""API tests against the synthetic app tables built by the ``app_data_dir`` fixture (cancer types SYN and SY2)."""

from __future__ import annotations

import importlib
import json
import math

import pytest
from fastapi.testclient import TestClient

from app.main import create_app

RANK_METRICS = [
    "log2fc_cv_median",
    "log2fc_mad_median",
    "log2fc_fano_median",
    "delta_cv2_resid",
    "cv_malignant",
    "cv_non_malignant",
    "mad_malignant",
    "mad_non_malignant",
    "fano_malignant",
    "fano_non_malignant",
    "cv2_resid_malignant",
    "cv2_resid_non_malignant",
]
TOP_ROW_KEYS = {
    "gene",
    "value",
    "value_q",
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
}
PATIENT_ROW_KEYS = {"patient", "compartment", "n_cells", "mean", "detection_rate", "cv", "mad", "fano", "cv2_resid", "cv_ds", "cv2_resid_ds"}


def _reject_constant(name: str):
    raise AssertionError(f"non-finite JSON token leaked into the response: {name}")


def body(response, status=200):
    """Parsed JSON body; fails if the document contains NaN/Infinity tokens."""
    assert response.status_code == status, response.text
    assert response.headers["content-type"].startswith("application/json")
    return json.loads(response.content, parse_constant=_reject_constant)


def _walk(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _walk(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _walk(v)
    else:
        yield obj


def assert_finite(obj):
    for v in _walk(obj):
        assert v is None or isinstance(v, (str, bool, int)) or (isinstance(v, float) and math.isfinite(v)), v


@pytest.fixture(scope="module")
def client(app_data_dir):
    return TestClient(create_app(app_data_dir))


@pytest.fixture(scope="module")
def a_gene(client):
    """A gene that has per-patient rows in SYN (i.e. expressed in both compartments)."""
    rows = body(client.get("/api/cancer-types/SYN/top-genes", params={"metric": "log2fc_cv_median", "n": 1}))["rows"]
    return rows[0]["gene"]


# ------------------------------------------------------------------ health / overview
def test_health(client):
    data = body(client.get("/api/health"))
    assert data["status"] == "ok"
    assert data["n_datasets"] == 2
    assert data["n_genes"] > 0
    assert set(data["build"]["datasets"]) == {"SYN_Test", "SYN_Two"}
    assert "sctcna_version" in data["build"]


def test_cancer_types(client):
    rows = body(client.get("/api/cancer-types"))
    assert [r["cancer_type"] for r in rows] == ["SY2", "SYN"]
    expected = {
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
    }
    assert set(rows[0]) == expected
    syn = next(r for r in rows if r["cancer_type"] == "SYN")
    assert syn["dataset"] == "SYN_Test" and syn["n_patients"] == 6 and syn["n_patients_paired"] == 6


def test_datasets(client):
    rows = body(client.get("/api/datasets"))
    assert len(rows) == 2
    for col in ("dataset", "cancer_type", "depth2_cell_p", "sha256", "runtime_seconds", "sctcna_version"):
        assert col in rows[0]
    assert_finite(rows)


# ------------------------------------------------------------------ genes
def test_gene_search(client):
    data = body(client.get("/api/genes/search", params={"q": "pois", "limit": 5}))
    assert data["query"] == "pois"
    assert len(data["results"]) == 5
    assert all(r["gene"].startswith("POIS") for r in data["results"])
    assert all(r["n_cancer_types"] == 2 for r in data["results"])
    genes = [r["gene"] for r in data["results"]]
    assert genes == sorted(genes)


def test_gene_search_exact_match_first(client):
    # over-dispersed genes are named OD30..OD69, so "od3" matches many of them and "od30" is itself a symbol
    data = body(client.get("/api/genes/search", params={"q": "od3"}))
    assert len(data["results"]) == 10 and all(r["gene"].startswith("OD3") for r in data["results"])
    data = body(client.get("/api/genes/search", params={"q": "od30"}))
    assert data["results"][0]["gene"] == "OD30"


def test_gene_search_bad_input(client):
    assert client.get("/api/genes/search", params={"q": "   "}).status_code == 400
    assert client.get("/api/genes/search", params={"q": ""}).status_code == 400
    assert client.get("/api/genes/search").status_code == 400
    assert client.get("/api/genes/search", params={"q": "P", "limit": 0}).status_code == 422
    assert client.get("/api/genes/search", params={"q": "P", "limit": 101}).status_code == 422
    assert body(client.get("/api/genes/search", params={"q": "ZZZNOPE"}))["results"] == []


def test_gene_detail(client, a_gene):
    data = body(client.get(f"/api/genes/{a_gene.lower()}"))
    assert data["gene"] == a_gene
    assert [c["cancer_type"] for c in data["cancers"]] == ["SY2", "SYN"]
    row = data["cancers"][1]
    for col in ("cv_malignant", "cv_non_malignant", "log2fc_cv_median", "p_value", "q_value", "expressed", "cancer_name", "dataset", "cv_pooled_malignant"):
        assert col in row
    assert row["dataset"] == "SYN_Test" and row["cancer_name"] == "Synthetic test data"
    assert_finite(data)


def test_gene_detail_unknown(client):
    r = client.get("/api/genes/NOSUCHGENE")
    assert r.status_code == 404
    assert "detail" in r.json()


def test_gene_patients(client, a_gene):
    data = body(client.get(f"/api/genes/{a_gene}/patients", params={"cancer_type": "SYN"}))
    assert data["gene"] == a_gene and data["cancer_type"] == "SYN" and data["dataset"] == "SYN_Test"
    assert len(data["rows"]) == 12  # 6 patients x 2 compartments
    assert set(data["rows"][0]) == PATIENT_ROW_KEYS
    assert {r["compartment"] for r in data["rows"]} == {"malignant", "non_malignant"}
    assert_finite(data)


def test_gene_patients_errors(client, a_gene):
    assert client.get(f"/api/genes/{a_gene}/patients").status_code == 422
    assert client.get(f"/api/genes/{a_gene}/patients", params={"cancer_type": "NOPE"}).status_code == 404
    assert client.get("/api/genes/NOSUCHGENE/patients", params={"cancer_type": "SYN"}).status_code == 404
    # per-patient tables only hold expressed genes, so a known but non-expressed gene is a 404 as well
    all_genes = {r["gene"] for r in body(client.get("/api/cancer-types/SYN/top-genes", params={"metric": "cv_malignant", "n": 500, "expressed_only": "false"}))["rows"]}
    expressed = {r["gene"] for r in body(client.get("/api/cancer-types/SYN/top-genes", params={"metric": "cv_malignant", "n": 500}))["rows"]}
    not_expressed = sorted(all_genes - expressed)
    if not_expressed:
        assert client.get(f"/api/genes/{not_expressed[0]}/patients", params={"cancer_type": "SYN"}).status_code == 404


def test_gene_celltypes(client, a_gene):
    data = body(client.get(f"/api/genes/{a_gene}/celltypes", params={"cancer_type": "SYN"}))
    assert data["gene"] == a_gene and data["cancer_type"] == "SYN"
    cell_types = {r["cell_type"] for r in data["rows"]}
    assert "Malignant" in cell_types and len(cell_types) >= 2
    assert set(data["rows"][0]) == {"cell_type", "n_patients", "mean", "detection_rate", "cv", "mad", "fano", "cv2_resid"}
    assert client.get(f"/api/genes/{a_gene}/celltypes").status_code == 422
    assert client.get(f"/api/genes/{a_gene}/celltypes", params={"cancer_type": "NOPE"}).status_code == 404
    assert client.get("/api/genes/NOSUCHGENE/celltypes", params={"cancer_type": "SYN"}).status_code == 404


# ------------------------------------------------------------------ top genes
@pytest.mark.parametrize("metric", RANK_METRICS)
def test_top_genes_every_metric(client, metric):
    data = body(client.get("/api/cancer-types/SYN/top-genes", params={"metric": metric, "n": 10}))
    assert data["metric"] == metric and data["direction"] == "noisy" and data["n"] == 10
    assert 0 < len(data["rows"]) <= 10
    assert set(data["rows"][0]) == TOP_ROW_KEYS
    values = [r["value"] for r in data["rows"]]
    assert all(v is not None for v in values)
    assert values == sorted(values, reverse=True)


def test_top_genes_stable_and_filters(client):
    stable = body(client.get("/api/cancer-types/SYN/top-genes", params={"metric": "log2fc_cv_median", "direction": "stable", "n": 5}))
    values = [r["value"] for r in stable["rows"]]
    assert values == sorted(values)
    noisy = body(client.get("/api/cancer-types/SYN/top-genes", params={"metric": "log2fc_cv_median", "n": 500}))
    everything = body(client.get("/api/cancer-types/SYN/top-genes", params={"metric": "log2fc_cv_median", "n": 500, "expressed_only": "false"}))
    assert len(everything["rows"]) >= len(noisy["rows"])
    # simulated NOISY genes are over-dispersed only in malignant cells, so they should lead the ranking
    assert any(r["gene"].startswith("NOISY") for r in noisy["rows"][:25])


def test_top_genes_errors(client):
    assert client.get("/api/cancer-types/NOPE/top-genes", params={"metric": "log2fc_cv_median"}).status_code == 404
    assert client.get("/api/cancer-types/SYN/top-genes", params={"metric": "bogus"}).status_code == 400
    assert client.get("/api/cancer-types/SYN/top-genes", params={"metric": "log2fc_cv_median", "direction": "sideways"}).status_code == 400
    assert client.get("/api/cancer-types/SYN/top-genes", params={"metric": "log2fc_cv_median", "n": 0}).status_code == 422
    assert client.get("/api/cancer-types/SYN/top-genes", params={"metric": "log2fc_cv_median", "n": 501}).status_code == 422
    assert client.get("/api/cancer-types/SYN/top-genes").status_code == 422


# ------------------------------------------------------------------ comparison
@pytest.mark.parametrize("metric", ["cv", "mad", "fano", "cv2_resid", "sd"])
def test_comparison(client, metric):
    data = body(client.get("/api/cancer-types/SYN/comparison", params={"metric": metric}))
    assert data["cancer_type"] == "SYN" and data["metric"] == metric
    n = data["n_genes"]
    assert n > 0
    for key in ("genes", "malignant", "non_malignant", "log2fc", "q_value"):
        assert len(data[key]) == n
    assert all(v is not None for v in data["malignant"]) and all(v is not None for v in data["non_malignant"])
    summary = data["summary"]
    assert set(summary) == {"median_log2fc", "n_sig_q05", "n_sig_noisier_malignant", "n_sig_noisier_non_malignant", "n_patients_paired"}
    assert summary["n_patients_paired"] == 6
    assert summary["n_sig_q05"] >= summary["n_sig_noisier_malignant"] + summary["n_sig_noisier_non_malignant"]
    if metric == "sd":
        i = 0
        assert data["log2fc"][i] == pytest.approx(math.log2(data["malignant"][i] / data["non_malignant"][i]))
    assert_finite(data)


def test_comparison_errors(client):
    assert client.get("/api/cancer-types/NOPE/comparison", params={"metric": "cv"}).status_code == 404
    assert client.get("/api/cancer-types/SYN/comparison", params={"metric": "variance"}).status_code == 400
    # repeated calls hit the cache and stay identical
    a = client.get("/api/cancer-types/SYN/comparison", params={"metric": "cv"}).content
    b = client.get("/api/cancer-types/SYN/comparison", params={"metric": "cv"}).content
    assert a == b


# ------------------------------------------------------------------ heatmap
def test_heatmap(client):
    data = body(client.get("/api/heatmap", params={"genes": "od30, pois100\nNOISY10 CONST0 od30 NOPE1", "metric": "log2fc_cv_median"}))
    assert data["metric"] == "log2fc_cv_median"
    assert data["genes"] == ["OD30", "POIS100", "NOISY10", "CONST0"]
    assert data["missing"] == ["NOPE1"]
    assert data["cancer_types"] == ["SY2", "SYN"]
    assert len(data["values"]) == 4 and all(len(row) == 2 for row in data["values"])
    assert len(data["q_values"]) == 4 and all(len(row) == 2 for row in data["q_values"])
    assert_finite(data)


def test_heatmap_no_qvalues_for_compartment_metric(client):
    data = body(client.get("/api/heatmap", params={"genes": "OD30,OD31", "metric": "cv_malignant"}))
    assert data["q_values"] is None
    assert all(v is not None for row in data["values"] for v in row)


def test_heatmap_expressed_only(client):
    full = body(client.get("/api/heatmap", params={"genes": "OD30", "metric": "cv_malignant"}))
    filtered = body(client.get("/api/heatmap", params={"genes": "OD30", "metric": "cv_malignant", "expressed_only": "true"}))
    assert len(filtered["values"]) == 1
    assert full["cancer_types"] == filtered["cancer_types"]


def test_heatmap_errors(client):
    assert client.get("/api/heatmap", params={"genes": "OD30", "metric": "bogus"}).status_code == 400
    assert client.get("/api/heatmap", params={"genes": "", "metric": "log2fc_cv_median"}).status_code == 400
    assert client.get("/api/heatmap", params={"genes": "NOPE1,NOPE2", "metric": "log2fc_cv_median"}).status_code == 400
    too_many = ",".join(f"G{i}" for i in range(201))
    assert client.get("/api/heatmap", params={"genes": too_many, "metric": "log2fc_cv_median"}).status_code == 400
    at_limit = ",".join(f"G{i}" for i in range(199)) + ",OD30"
    assert client.get("/api/heatmap", params={"genes": at_limit, "metric": "log2fc_cv_median"}).status_code == 200


# ------------------------------------------------------------------ tumour heterogeneity
@pytest.mark.parametrize("score", ["depth2", "depth"])
def test_tumor_heterogeneity_cell_level(client, score):
    data = body(client.get("/api/cancer-types/SYN/tumor-heterogeneity", params={"score": score, "level": "cell"}))
    assert data["cancer_type"] == "SYN" and data["score"] == score and data["level"] == "cell"
    comps = data["compartments"]
    assert set(comps) == {"malignant", "non_malignant"}
    scored = ["malignant", "non_malignant"] if score == "depth2" else ["malignant"]  # DEPTH: non-malignant cells are the reference
    for name, comp in comps.items():
        if name in scored:
            assert 0 < len(comp["values"]) <= 5000
            assert comp["n"] >= len(comp["values"])
            assert comp["median"] is not None
        else:
            assert comp["values"] == [] and comp["n"] == 0 and comp["median"] is None
    assert len(data["patients"]) == 12
    assert set(data["patients"][0]) == {"patient", "compartment", "n_cells", "median"}
    test = data["paired_test"]
    assert set(test) == {"n_pairs", "median_diff", "p_value", "rank_biserial"}
    assert test["n_pairs"] == (6 if score == "depth2" else 0)
    assert_finite(data)
    # deterministic subsample: identical bodies on repeat
    assert client.get("/api/cancer-types/SYN/tumor-heterogeneity", params={"score": score, "level": "cell"}).json() == data


def test_tumor_heterogeneity_sample_level(client):
    data = body(client.get("/api/cancer-types/SYN/tumor-heterogeneity", params={"score": "depth2", "level": "sample"}))
    assert data["level"] == "sample"
    assert len(data["rows"]) == 12
    assert set(data["rows"][0]) == {"patient", "compartment", "n_cells", "depth2_pseudobulk", "depth_pseudobulk", "depth_cell_median", "depth2_cell_median"}
    non_mal = [r for r in data["rows"] if r["compartment"] == "non_malignant"]
    assert all(r["depth_pseudobulk"] is None for r in non_mal)  # NaN in the table must come out as null
    assert data["paired_test"]["n_pairs"] == 6
    depth = body(client.get("/api/cancer-types/SYN/tumor-heterogeneity", params={"score": "depth", "level": "sample"}))
    assert depth["paired_test"]["n_pairs"] == 0 and depth["paired_test"]["p_value"] is None
    assert_finite(depth)


def test_tumor_heterogeneity_errors(client):
    assert client.get("/api/cancer-types/NOPE/tumor-heterogeneity").status_code == 404
    assert client.get("/api/cancer-types/SYN/tumor-heterogeneity", params={"score": "depth3"}).status_code == 400
    assert client.get("/api/cancer-types/SYN/tumor-heterogeneity", params={"level": "patient"}).status_code == 400


def test_tumor_heterogeneity_summary(client):
    rows = body(client.get("/api/tumor-heterogeneity/summary"))["rows"]
    assert [r["cancer_type"] for r in rows] == ["SY2", "SYN"]
    expected = {
        "cancer_type",
        "cancer_name",
        "n_patients_paired",
        "depth2_cell_median_malignant",
        "depth2_cell_median_non_malignant",
        "depth2_cell_median_diff",
        "depth2_cell_p",
        "depth2_pseudobulk_median_diff",
        "depth2_pseudobulk_p",
    }
    assert set(rows[0]) == expected
    assert rows[1]["depth2_cell_median_malignant"] is not None
    assert_finite(rows)


# ------------------------------------------------------------------ downloads
def _csv(response, status=200):
    assert response.status_code == status, response.text
    assert response.headers["content-type"].startswith("text/csv")
    assert response.headers["content-disposition"].startswith("attachment; filename=")
    return response.text.splitlines()


def test_download_datasets_and_gene_summary(client, a_gene):
    lines = _csv(client.get("/api/download/datasets.csv"))
    assert lines[0].startswith("dataset,cancer_type,") and len(lines) == 3
    lines = _csv(client.get("/api/download/gene_summary.csv"))
    assert "log2fc_cv_median" in lines[0] and len(lines) > 300
    syn_only = _csv(client.get("/api/download/gene_summary.csv", params={"cancer_type": "SYN"}))
    assert len(syn_only) < len(lines)
    assert client.get("/api/download/gene_summary.csv", params={"cancer_type": "NOPE"}).status_code == 404


def test_download_scores(client):
    lines = _csv(client.get("/api/download/sample_scores.csv"))
    assert len(lines) == 1 + 12 + 10
    assert len(_csv(client.get("/api/download/sample_scores.csv", params={"cancer_type": "SY2"}))) == 11
    cells = _csv(client.get("/api/download/cell_scores.csv", params={"cancer_type": "SY2"}))
    assert cells[0] == "dataset,cancer_type,patient,cell_type,compartment,n_counts,n_genes,pct_mito,depth,depth2"
    assert len(cells) > 100
    assert client.get("/api/download/cell_scores.csv").status_code == 422
    assert client.get("/api/download/cell_scores.csv", params={"cancer_type": "NOPE"}).status_code == 404


def test_download_top_genes(client):
    lines = _csv(client.get("/api/download/top_genes.csv", params={"cancer_type": "SYN", "metric": "log2fc_cv_median", "n": 7}))
    assert lines[0].startswith("gene,value,value_q,mean_malignant,") and len(lines) == 8
    assert client.get("/api/download/top_genes.csv", params={"metric": "log2fc_cv_median"}).status_code == 422
    assert client.get("/api/download/top_genes.csv", params={"cancer_type": "SYN"}).status_code == 422
    assert client.get("/api/download/top_genes.csv", params={"cancer_type": "SYN", "metric": "bogus"}).status_code == 400
    assert client.get("/api/download/top_genes.csv", params={"cancer_type": "SYN", "metric": "log2fc_cv_median", "n": 0}).status_code == 422


def test_download_gene_patients_and_heatmap(client, a_gene):
    lines = _csv(client.get("/api/download/gene_patients.csv", params={"gene": a_gene.lower(), "cancer_type": "SYN"}))
    assert lines[0] == "gene,cancer_type,patient,compartment,n_cells,mean,detection_rate,cv,mad,fano,cv2_resid,cv_ds,cv2_resid_ds"
    assert len(lines) == 13 and lines[1].startswith(f"{a_gene},SYN,")
    assert client.get("/api/download/gene_patients.csv", params={"gene": a_gene}).status_code == 422
    assert client.get("/api/download/gene_patients.csv", params={"gene": "NOSUCHGENE", "cancer_type": "SYN"}).status_code == 404
    lines = _csv(client.get("/api/download/heatmap.csv", params={"genes": "OD30,OD31,NOPE", "metric": "log2fc_cv_median"}))
    assert lines[0] == "gene,SY2,SYN,SY2_q,SYN_q" and len(lines) == 3
    assert client.get("/api/download/heatmap.csv", params={"metric": "log2fc_cv_median"}).status_code == 422
    assert client.get("/api/download/heatmap.csv", params={"genes": "OD30", "metric": "bogus"}).status_code == 400


def test_download_unknown_table(client):
    r = client.get("/api/download/secrets.csv")
    assert r.status_code == 404
    assert "detail" in r.json()


# ------------------------------------------------------------------ missing data directory
def test_missing_data_dir_gives_503(tmp_path):
    client = TestClient(create_app(tmp_path / "does-not-exist"))
    health = body(client.get("/api/health"))
    assert health["status"] == "no-data" and health["n_datasets"] == 0 and health["build"] is None
    for path in ("/api/cancer-types", "/api/datasets", "/api/genes/OD30", "/api/heatmap?genes=OD30&metric=cv_malignant", "/api/download/datasets.csv"):
        r = client.get(path)
        assert r.status_code == 503, path
        assert "app data not available" in r.json()["detail"]


def test_module_level_app_survives_missing_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("SCTCNA_DATA_DIR", str(tmp_path / "nothing-here"))
    import app.main as main_module

    reloaded = importlib.reload(main_module)
    assert reloaded.app.state.store is None
    assert TestClient(reloaded.app).get("/api/health").json()["status"] == "no-data"


# ------------------------------------------------------------------ noise contrasts
def test_noise_contrasts_all_and_filtered(client):
    rows = body(client.get("/api/noise-contrasts"))["rows"]
    assert {r["cancer_type"] for r in rows} == {"SYN", "SY2"}
    refs = {r["reference"] for r in rows if r["cancer_type"] == "SYN"}
    assert "non_malignant" in refs and refs - {"non_malignant"}  # pooled reference plus single lineages
    keys = {"reference", "n_patients", "median_log2fc_cv", "p_log2fc_cv", "median_log2fc_cv_ds", "median_excess_cv2_ds", "p_excess_cv2_ds", "frac_patients_excess_positive"}
    assert keys <= set(rows[0])
    syn = body(client.get("/api/noise-contrasts", params={"cancer_type": "SYN"}))["rows"]
    assert syn and all(r["cancer_type"] == "SYN" for r in syn)
    assert client.get("/api/noise-contrasts", params={"cancer_type": "NOPE"}).status_code == 404


def test_noise_contrasts_per_patient(client):
    data = body(client.get("/api/cancer-types/SYN/noise-contrasts"))
    assert data["cancer_type"] == "SYN" and data["summary"]
    allref = [r for r in data["rows"] if r["reference"] == "non_malignant"]
    assert len({r["patient"] for r in allref}) == 6
    assert client.get("/api/cancer-types/NOPE/noise-contrasts").status_code == 404


def test_download_noise_contrasts(client):
    lines = _csv(client.get("/api/download/noise_contrasts.csv", params={"cancer_type": "SY2"}))
    header = lines[0].split(",")
    col = header.index("cancer_type")
    assert len(lines) > 1 and all(line.split(",")[col] == "SY2" for line in lines[1:])


def test_heterogeneity_tests_are_precomputed(client):
    store = client.app.state.store
    assert ("SYN", "depth2", "cell") in store._heterogeneity_tests
    data = body(client.get("/api/cancer-types/SYN/tumor-heterogeneity", params={"score": "depth2", "level": "cell"}))
    assert data["paired_test"]["n_pairs"] == 6
