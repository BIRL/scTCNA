"""Browser smoke tests: every page loads against the synthetic app tables and renders its plots.

Run with ``pytest -m e2e tests/test_e2e.py`` (needs ``python -m playwright install chromium``).
A uvicorn server is started on a free port with SCTCNA_DATA_DIR pointing at the synthetic
tables built by the ``app_data_dir`` fixture (cancer types SYN and SY2; genes NOISY*, POIS*, OD*).
"""

from __future__ import annotations

import os
import socket
import subprocess
import sys
import time
from pathlib import Path

import httpx
import pytest

pytestmark = pytest.mark.e2e
PLOT = ".js-plotly-plot"
REPO = Path(__file__).resolve().parent.parent


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="module")
def server(app_data_dir):
    port = _free_port()
    env = {**os.environ, "SCTCNA_DATA_DIR": str(app_data_dir)}
    proc = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "app.main:app", "--host", "127.0.0.1", "--port", str(port)],
        cwd=REPO,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    base = f"http://127.0.0.1:{port}"
    deadline = time.time() + 60
    while time.time() < deadline:
        try:
            if httpx.get(base + "/api/health", timeout=2).status_code == 200:
                break
        except httpx.HTTPError:
            pass
        time.sleep(0.3)
    else:
        proc.kill()
        pytest.fail("uvicorn did not start within 60 s")
    yield base
    proc.terminate()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        proc.kill()


@pytest.fixture
def errors(page):
    """Console errors and uncaught exceptions raised while a test drives the page."""
    found: list[str] = []
    page.on("console", lambda m: found.append(m.text) if m.type == "error" and "favicon" not in m.text else None)
    page.on("pageerror", lambda e: found.append(f"pageerror: {e}"))
    return found


def test_home_page_loads_with_sample_counts(server, page, errors):
    page.goto(server + "/")
    assert "scTCNA" in page.title()
    page.wait_for_selector(f"#plot-cells{PLOT}", timeout=15_000)
    page.wait_for_selector(f"#plot-patients{PLOT}", timeout=15_000)
    assert page.locator("#cancer-form").is_hidden()  # cancer-specific inputs stay hidden in pan-cancer mode
    assert errors == []


def test_gene_search_renders_plots_and_table(server, page, errors):
    page.goto(server + "/")
    page.fill("#gene-input", "NOISY")
    page.wait_for_selector(".autocomplete-list li", timeout=10_000)
    gene = page.locator(".autocomplete-list li").first.inner_text().split()[0]
    page.locator(".autocomplete-list li").first.click()
    page.wait_for_selector(f"#plot-cv{PLOT}", timeout=15_000)
    assert page.locator(f"#pan-results {PLOT}").count() >= 4
    assert gene in page.locator("#pan-title").inner_text()
    assert page.locator("#pan-table tbody tr").count() == 2  # one row per cancer type
    assert errors == []


def test_cancer_specific_mode(server, page, errors):
    page.goto(server + "/")
    page.locator("label[for=mode-cancer]").click()
    assert page.locator("#pan-form").is_hidden()
    page.wait_for_selector("#cancer-select option", state="attached", timeout=10_000)
    page.fill("#genes-text", "NOISY10, POIS100\nOD60")
    page.locator("#cancer-form button[type=submit]").click()
    page.wait_for_selector(f"#plot-cs-cv{PLOT}", timeout=15_000)
    assert page.locator("#cs-table tbody tr").count() == 3
    assert errors == []


def test_summary_page(server, page, errors):
    page.goto(server + "/summary.html")
    page.wait_for_selector(f"#plot-all{PLOT}", timeout=15_000)
    page.wait_for_selector(f"#plot-lineages{PLOT}", timeout=15_000)
    page.wait_for_selector(f"#plot-patients{PLOT}", timeout=15_000)
    assert page.locator("#table-all tbody tr").count() == 2
    assert errors == []


def test_top_genes_page(server, page, errors):
    page.goto(server + "/top-genes.html")
    page.wait_for_selector("#table tbody tr", timeout=15_000)
    assert page.locator("#table tbody tr").count() > 10
    page.wait_for_selector(f"#plot-top{PLOT}", timeout=15_000)
    assert errors == []


def test_tumor_page_violins(server, page, errors):
    page.goto(server + "/tumor.html")
    page.wait_for_selector(f"#plot-main{PLOT} .violinlayer path", state="attached", timeout=15_000)
    page.wait_for_function("document.getElementById('stats').textContent.includes('patient pairs')", timeout=15_000)
    assert "Paired Wilcoxon" in page.locator("#stats").inner_text()
    assert errors == []


def test_heatmap_page(server, page, errors):
    page.goto(server + "/heatmap.html")
    page.fill("#genes", "NOISY10 NOISY11 POIS100 NOTAGENE")
    page.locator("#form button[type=submit]").click()
    page.wait_for_selector(f"#plot-heatmap{PLOT}", timeout=15_000)
    assert "NOTAGENE" in page.locator("#missing").inner_text()
    assert errors == []


def test_about_page_lists_datasets(server, page, errors):
    page.goto(server + "/about.html")
    page.wait_for_selector("#datasets-table tbody tr", timeout=15_000)
    assert page.locator("#datasets-table tbody tr").count() == 2
    assert "Cancer Noise Atlas" in page.content()
    assert errors == []
