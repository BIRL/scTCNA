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
    assert page.locator("#pan-table thead th").count() == 3  # metrics in rows: the corner cell plus one column per cancer type
    assert errors == []


def test_cancer_specific_mode(server, page, errors):
    page.goto(server + "/")
    page.locator("label[for=mode-cancer]").click()
    assert page.locator("#pan-form").is_hidden()
    page.wait_for_selector("#cancer-select option", state="attached", timeout=10_000)
    page.fill("#genes-text", "NOISY10, POIS100\nOD60")
    page.locator("#cancer-form button[type=submit]").click()
    page.wait_for_selector(f"#plot-cs-cv{PLOT}", timeout=15_000)
    assert page.locator("#cs-table thead th").count() == 4  # corner cell plus one column per gene
    assert errors == []


def test_summary_page(server, page, errors):
    page.goto(server + "/summary.html")
    page.wait_for_selector(f"#plot-all{PLOT}", timeout=15_000)
    page.wait_for_selector(f"#plot-lineages{PLOT}", timeout=15_000)
    page.wait_for_selector(f"#plot-patients{PLOT}", timeout=15_000)
    assert page.locator("#table-all thead th").count() == 3  # corner cell plus one column per cancer type
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


BAR_PATHS = "ps => ps.map(p => p.getAttribute('d')).join('|')"
SVG_TEXT = "ts => ts.map(t => t.textContent)"  # inner_text() is None for SVG nodes


def _pan_gene(page, gene):
    page.fill("#gene-input", gene)
    page.locator("#pan-form button[type=submit]").click()
    page.wait_for_function(f"document.getElementById('pan-title').textContent.startsWith('{gene} ')", timeout=15_000)
    page.wait_for_selector(f"#plot-cv{PLOT} .main-svg", timeout=15_000)


def _cancer_mode(page, server, cancer="SYN"):
    page.goto(server + "/?mode=cancer")
    page.wait_for_selector("#cancer-select option", state="attached", timeout=15_000)
    page.select_option("#cancer-select", cancer)


def _analyse(page, genes):
    page.fill("#genes-text", genes)
    page.locator("#cancer-form button[type=submit]").click()


def test_gene_plots_update_on_every_new_selection(server, page, errors):
    """A1: the second and third analyses redraw the plots on screen, not in a detached node."""
    page.goto(server + "/")
    seen = []
    for gene in ("NOISY10", "OD60", "POIS100"):
        _pan_gene(page, gene)
        page.wait_for_function("document.querySelector('#plot-cv').data[0].y.length === 2 && document.querySelector('#pan-table thead').textContent.length > 0", timeout=15_000)
        paths = page.locator("#plot-cv .bars path").evaluate_all(BAR_PATHS)
        assert paths and paths not in seen
        seen.append(paths)
    _cancer_mode(page, server)
    for genes, first in (("NOISY10, POIS100", "NOISY10"), ("OD60 OD61 OD62", "OD60")):
        _analyse(page, genes)
        page.wait_for_function(f"Array.from(document.querySelectorAll('#plot-cs-cv .xtick text')).some(t => t.textContent === '{first}')", timeout=15_000)
        ticks = page.locator("#plot-cs-cv .xtick text").evaluate_all(SVG_TEXT)
        assert ticks == genes.replace(",", "").split()
    assert errors == []


def test_summary_patient_plot_follows_selection(server, page, errors):
    """A1 on the summary page: changing the cancer type and then the measure redraws the per-patient plot."""
    page.goto(server + "/summary.html")
    page.wait_for_selector(f"#plot-patients{PLOT} .gtitle", timeout=15_000)
    first = page.locator("#plot-patients .gtitle").text_content()
    other = "SY2" if "SY2" not in first else "SYN"
    page.select_option("#cancer", other)
    page.wait_for_function(f"document.querySelector('#plot-patients .gtitle').textContent.includes('{other}')", timeout=15_000)
    ytitle = page.locator("#plot-patients .ytitle").text_content()
    assert page.locator("#measure").input_value() == "excess_cv2_ds"  # the default measure
    page.select_option("#measure", "log2fc_cv")
    page.wait_for_function(f"document.querySelector('#plot-patients .ytitle').textContent !== {ytitle!r}", timeout=15_000)
    assert "raw" in page.locator("#plot-patients .gtitle").text_content()
    assert errors == []


def test_pan_cancer_type_checkboxes(server, page, errors):
    """B3: unticking a cancer type drops its column; Clear explains the empty view; the choice is kept in the URL."""
    page.goto(server + "/")
    page.wait_for_selector("#pan-checks input", state="attached", timeout=15_000)
    assert page.locator("#pan-checks input:checked").count() == 2
    _pan_gene(page, "NOISY10")
    assert page.locator("#pan-table thead th").count() == 3
    page.locator("#pan-checks input[value=SY2]").uncheck()
    page.wait_for_function("document.querySelectorAll('#pan-table thead th').length === 2", timeout=15_000)
    assert "cancer_types=SYN" in page.url
    assert "1 of 2" in page.locator("#pan-count").inner_text()
    page.reload()
    page.wait_for_selector("#pan-table thead th", timeout=15_000)
    assert page.locator("#pan-checks input:checked").count() == 1
    page.locator("#pan-none").click()
    page.wait_for_function("document.getElementById('pan-note').textContent.includes('Select at least one cancer type')", timeout=15_000)
    assert page.locator("#plot-cv").inner_text() == "No cancer type selected."
    page.locator("#pan-all").click()
    page.wait_for_function("document.querySelectorAll('#pan-table thead th').length === 3", timeout=15_000)
    assert "cancer_types" not in page.url
    assert errors == []


def test_gene_notes_explain_missing_and_untested_genes(server, page, errors):
    """A2: unknown symbols and genes detected in < 10% of cells are explained inline."""
    page.goto(server + "/")
    _pan_gene(page, "NOISY24")  # detected in < 10% of cells in SYN only
    note = page.locator("#pan-note").inner_text()
    assert "NOISY24 is detected in < 10% of cells in SYN" in note
    _cancer_mode(page, server, "SYN")
    _analyse(page, "NOISY10 NOISY24 NOTAGENE")
    page.wait_for_selector(f"#plot-cs-diff{PLOT}", timeout=15_000)
    missing = page.locator("#cs-missing").inner_text()
    assert "Not found in any dataset: NOTAGENE." in missing
    assert "Detected in < 10% of cells" in missing and "NOISY24" in missing
    assert "Paired test" in page.locator("#cs-table").inner_text()
    _analyse(page, "NOTAGENE")
    page.wait_for_function("document.getElementById('genes-note').textContent.includes('NOTAGENE')", timeout=15_000)
    assert page.locator("#cancer-results").is_hidden()
    assert [e for e in errors if "status of 404" not in e] == []  # the browser logs the 404 of the unknown symbol


def test_gene_list_upload(server, page, errors):
    """B4: an uploaded .txt list is merged, upper-cased and de-duplicated; unknown symbols are listed."""
    upload = {"name": "genes.txt", "mimeType": "text/plain", "buffer": b"od60\nnoisy10, NOTAGENE od60\n"}
    _cancer_mode(page, server, "SYN")
    page.fill("#genes-text", "OD60")
    page.set_input_files("#genes-file", files=[upload])
    page.wait_for_function("document.getElementById('genes-note').textContent.includes('NOTAGENE')", timeout=15_000)
    assert page.locator("#genes-text").input_value().split() == ["OD60", "NOISY10", "NOTAGENE"]
    page.goto(server + "/heatmap.html")
    page.fill("#genes", "")
    page.set_input_files("#genes-file", files=[upload])
    page.wait_for_function("document.getElementById('genes-note').textContent.includes('NOTAGENE')", timeout=15_000)
    assert page.locator("#genes").input_value().split() == ["OD60", "NOISY10", "NOTAGENE"]
    assert errors == []


def test_many_genes_use_collapsible_metric_cards(server, page, errors):
    """B6: more than 8 genes show one collapsible card per metric, the first one open."""
    _cancer_mode(page, server, "SYN")
    _analyse(page, " ".join(f"OD{i}" for i in range(40, 50)))
    page.wait_for_selector("#cs-grid.collapsible", timeout=15_000)
    cards = page.locator("#cs-grid details.metric-card")
    assert cards.count() == 6
    assert cards.nth(0).get_attribute("open") is not None and cards.nth(1).get_attribute("open") is None
    cards.nth(1).locator("summary").click()
    page.wait_for_function("document.querySelectorAll('#cs-grid details.metric-card')[1].open", timeout=15_000)
    width = cards.nth(1).locator(".main-svg").first.bounding_box()["width"]
    assert width > 300
    assert page.locator("#cs-table thead th").count() == 11
    assert errors == []


def test_metric_table_row_label_sorts_columns(server, page, errors):
    """B5: clicking a metric row label sorts the columns by that row."""
    _cancer_mode(page, server, "SYN")
    _analyse(page, "OD60 OD61 OD62 OD63")
    page.wait_for_selector("#cs-table thead th", timeout=15_000)
    head = lambda: page.locator("#cs-table thead th").all_inner_texts()[1:]  # noqa: E731
    before = head()
    label = page.locator("#cs-table tbody th", has_text="log2 CV ratio (malignant / non-malignant), median").first
    label.click()
    first = head()
    label.click()
    second = head()
    assert sorted(first) == sorted(before) and second == first[::-1]
    assert errors == []


def test_info_button_opens_page_notes(server, page, errors):
    """C9: the "i" button next to the title opens the page notes; Escape closes them."""
    for path, text in (("/", "Pan-cancer"), ("/summary.html", "comparison group"), ("/about.html", "API")):
        page.goto(server + path)
        page.locator(".page-head .info-btn").click()
        panel = page.locator("#page-info-panel")
        assert panel.is_visible() and text in panel.inner_text()
        page.keyboard.press("Escape")
        assert panel.is_hidden()
    assert errors == []
