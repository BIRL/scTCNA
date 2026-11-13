# scTCNA — single-cell extension of The Cancer Noise Atlas (TCNA)

scTCNA applies the transcriptional-noise metrics of The Cancer Noise Atlas (TCNA) to
single-cell RNA-seq. Instead of tumour and normal bulk samples it compares malignant and
non-malignant cells of the same patients, in 10 cancer types, using one 10x Genomics study
per cancer type from the Curated Cancer Cell Atlas (3CA): 389,697 cells from 157 patients
after quality control (151,563 malignant, 238,134 non-malignant; 138 patients have both).
Per gene it computes TCNA's CV, CV², SD and MAD (plus Fano and detection rate) within each
patient x compartment group and TCNA's differential noise (Eq. 1, log2 of the CV ratio) per
patient; per cell and per pseudo-bulk profile it computes the tumour-heterogeneity scores
DEPTH and DEPTH2. Because single-cell noise depends on sequencing depth and mean expression,
it adds downsampling to a common UMI depth per dataset (the 10th percentile of per-cell UMI
totals, at least 1,000; cells below it are left unchanged), a mean-corrected CV² (LOWESS
residual), and paired tests with the patient as the unit (Wilcoxon signed-rank, BH correction
per cancer type). The results are precomputed tables served by a web app modelled on TCNA's
dashboards.

- TCNA portal: <https://tcna.lums.edu.pk>
- TCNA paper: Shah R, Faisal A, bin Faiz A, Bashir I, Jawaid A, Bashir K, Shoaib M,
  Chaudhary SU. The Cancer Noise Atlas (TCNA): A Webserver of Multiscale Transcriptional
  Noise in Cancer. *J Mol Biol* 2026. Code: <https://github.com/BIRL/TCNA>
- scTCNA demo: <https://sctcna.vercel.app>
- 3CA: <https://www.weizmann.ac.il/sites/3CA/>
- Results: [RESULTS.md](RESULTS.md); methods and assumptions: [METHODS.md](METHODS.md)

## Repository layout

```
sctcna/                 pipeline package (scanpy-based)
  io.py                 3CA / TISCH2 loaders -> AnnData with raw UMI counts
  qc.py                 cell and gene filters, mitochondrial fraction, Scrublet, depth downsampling
  normalize.py          CP10K and log2(x + 1)
  metrics.py            per-gene noise metrics, LOWESS residual, Eq. 1
  depth.py              DEPTH and DEPTH2, per cell and per pseudo-bulk
  stats.py              patient-level medians, paired Wilcoxon, BH, rank-biserial
  pipeline.py           process_dataset() and aggregate()
  synthetic.py          simulated data with known noise structure (tests, demo)
  cli.py                python -m sctcna run | aggregate
  web/                  JSON API over the app tables (create_api), mountable in other FastAPI apps
app/                    standalone web server (app.main:app) and the static frontend
configs/datasets.yaml   dataset registry and every pipeline parameter
scripts/                download_data.py, summarise_results.py, check_memory.py
tests/                  pytest suite (unit, API, pipeline integration) and Playwright browser tests
data/app/               aggregated tables served by the app (committed)
deploy/windows/         Windows Server 2019 set-up (IIS, NSSM)
docs/                   DEPLOY.md (Docker and Render), memory_check.txt
Dockerfile, render.yaml container image and Render Blueprint
requirements.txt        pinned full environment (pipeline, app, tests)
requirements-app.txt    pinned web-app dependencies only
METHODS.md, RESULTS.md  methods and assumptions; results per cancer type
```

`data/raw/` (downloads) and `data/processed/` (per-dataset outputs) are created by the
pipeline and are not committed.

## Quick start: the web app

The app serves the committed tables in `data/app`; no download or database is needed.
Python 3.12:

```bash
git clone https://github.com/BIRL/scTCNA.git && cd scTCNA
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements-app.txt
SCTCNA_DATA_DIR=data/app uvicorn app.main:app --host 127.0.0.1 --port 8000
curl http://127.0.0.1:8000/api/health    # {"status":"ok","n_datasets":10,...}
```

Then open <http://127.0.0.1:8000>. Without `SCTCNA_DATA_DIR` the app reads `data/app`.

- `pip install -e .` also works; it installs the full pipeline stack as well.
- A wheel (`pip wheel --no-deps .`) contains the `sctcna` package and its API
  (`sctcna.web`) but not `app/`; the standalone server runs from a checkout.
- With `make`: `make venv install` (full environment) and `make app`.

Pages: Gene Noise Analysis (pan-cancer and cancer-specific), Cross-cancer summary, Top genes,
Tumor Heterogeneity (DEPTH / DEPTH2), Heatmap, About / Methods.

## Reproduce the pipeline

Install the full environment (`pip install -r requirements.txt`, or `make venv install`), then:

```bash
python scripts/download_data.py --list              # show the registry
python scripts/download_data.py --all --convert     # download every enabled dataset to data/raw
python -m sctcna run                                # data/raw -> data/processed/<dataset>/
python -m sctcna aggregate                          # data/processed -> data/app
python scripts/summarise_results.py                 # the tables in RESULTS.md
```

The same steps are `make data`, `make pipeline` and `make aggregate`. Single datasets:
`--datasets LUAD_Kim2020 COAD_Lee2020` for both the download and `run`. Without any download,
`python -m sctcna run --synthetic --out data/processed_syn` processes simulated data.

`--convert` streams the expression matrix out of each 3CA archive into
`data/raw/<source>/<dataset>/<dataset>.h5ad` (raw UMI counts, capped at 3,000 cells per
patient x compartment) and deletes the archive unless `--keep-archives` is given. The 3CA
archives are 0.2-1.5 GB each; conversion needs about twice the archive size of free disk.
All thresholds are in the `params` section of `configs/datasets.yaml`.

## Deployment

| Target | Files | Instructions |
|---|---|---|
| Windows Server 2019 with IIS and NSSM | `deploy/windows/` | [deploy/windows/README_WINDOWS.md](deploy/windows/README_WINDOWS.md) |
| Any Docker host | `Dockerfile` | `make docker-build docker-run`, then <http://localhost:7860> |
| Render (free tier, 512 MB) | `render.yaml` | [docs/DEPLOY.md](docs/DEPLOY.md) |

The app's peak memory was 389 MB with all ten cancer types loaded
([docs/memory_check.txt](docs/memory_check.txt)).

The API can also run inside another FastAPI application. The TCNA portal mounts it at
`/api/sc` (branch `single-cell-mode` of <https://github.com/sabasiddique1/TCNA>):

```python
from sctcna.web import create_api
portal.mount("/api/sc", create_api("/path/to/data/app"))
```

## API

Routes of `sctcna/web/routes.py`, served under `/api` by the standalone app:

| Route | Returns |
|---|---|
| `GET /api/health` | build information, the number of loaded datasets and the package version; 503 without tables |
| `GET /api/cancer-types` | cancer types with dataset metadata |
| `GET /api/datasets` | the datasets table |
| `GET /api/genes/search?q=&cancer_type=` | gene name suggestions; with `cancer_type`, only genes measured in that cancer type |
| `GET /api/genes/lookup?genes=&cancer_type=` | sorts a list of up to 500 symbols into `found`, `unknown` and `not_in_cancer_type` |
| `GET /api/genes/{gene}` | per-cancer-type summary for one gene |
| `GET /api/genes/{gene}/patients?cancer_type=` | per-patient values for one gene |
| `GET /api/genes/{gene}/celltypes?cancer_type=` | per-cell-type values for one gene |
| `GET /api/cancer-types/{cancer_type}/top-genes?metric=` | genes ranked by a noise metric, differential (e.g. `log2fc_cv_median`) or per compartment (e.g. `cv_malignant`) |
| `GET /api/cancer-types/{cancer_type}/comparison` | malignant vs non-malignant values of one metric (`metric=`, default `cv`) per expressed gene; `expressed_only=false` for all genes |
| `GET /api/cancer-types/{cancer_type}/tumor-heterogeneity` | DEPTH or DEPTH2, per cell or per pseudo-bulk |
| `GET /api/heatmap?genes=` | genes x cancer types matrix (up to 200 genes) |
| `GET /api/tumor-heterogeneity/summary` | DEPTH2 summary per cancer type (per-cell and pseudo-bulk paired tests) |
| `GET /api/noise-contrasts` | malignant vs reference-group contrasts for every cancer type |
| `GET /api/cancer-types/{cancer_type}/noise-contrasts` | the same contrasts per patient |
| `GET /api/download/{table}.csv` | CSV export of a served table |

## Tests

```bash
make test        # 127 tests: unit, API and pipeline integration, on synthetic data
make test-e2e    # 16 Playwright browser tests (needs: python -m playwright install chromium)
make lint        # ruff check and ruff format --check
python scripts/check_memory.py --data data/app   # peak memory of the app
```

| File | Tests |
|---|---|
| `tests/test_api.py` | 57 |
| `tests/test_metrics.py` | 16 |
| `tests/test_pipeline.py` | 16 |
| `tests/test_io.py` | 13 |
| `tests/test_depth.py` | 9 |
| `tests/test_e2e.py` | 16 |
| `tests/test_qc.py` | 8 |
| `tests/test_stats.py` | 8 |

## Methods and assumptions

[METHODS.md](METHODS.md) describes data handling, QC, normalisation, every metric formula,
the statistics, DEPTH / DEPTH2 and the output tables, and lists which choices come from TCNA
and which are scTCNA assumptions. In short: every metric is computed per gene within one
patient x compartment group; the patient is the unit of replication; cancer-type values are
medians across patients; differential noise is TCNA Eq. 1 per patient, tested with a paired
Wilcoxon signed-rank test and BH correction per cancer type.

| | TCNA | scTCNA |
|---|---|---|
| Data | Bulk RNA-seq from GDC, Primary Tumor vs Solid Tissue Normal, 26 sites | 10x scRNA-seq from 3CA, 10 cancer types, one study each |
| Contrast | Tumour vs normal samples | Malignant vs non-malignant cells of the same patients |
| Expression scale | TPM / FPKM / FPKM-UQ | CP10K; log2(CP10K + 1) per cell and log2(CPM + 1) per pseudo-bulk for DEPTH / DEPTH2 |
| Gene noise metrics | CV, CV², SD, MAD | the same, plus Fano, detection rate and mean-corrected CV² |
| Statistics | differential noise as log2 fold change | paired Wilcoxon across patients (at least 6), rank-biserial effect size, BH per cancer type |
| Robustness checks | none | downsampling to a common UMI depth; depth-matched, mean-corrected excess noise |
| Pathway analysis, custom upload | yes | no |

Known limitations (detailed in METHODS.md): one study per cancer type; the non-malignant
compartment mixes several lineages; malignant calls come from 3CA and are not re-derived;
the melanoma study is metastatic; 10x UMI data only; cross-cancer comparisons are confounded
by study effects; 6-27 paired patients per cancer type limit the power of the paired tests.

## Data

Primary source: 3CA, the Curated Cancer Cell Atlas (Gavish et al., *Nature* 2023). Each
study bundle provides raw UMI counts, curated cell types (malignant vs annotated
non-malignant lineages) and patient identifiers. Secondary source: TISCH2 (Han et al.,
*Nucleic Acids Res* 2023); its two entries in the registry are disabled. Counts after QC,
from `data/app/datasets.parquet`:

| Dataset | Cancer type (TCNA site) | Study | Site | Patients (paired) | Cells | Malignant / non-malignant |
|---|---|---|---|---|---|---|
| `BRCA_Wu2021` | BRCA (Breast) | Wu et al. 2021 | primary | 26 (20) | 74,723 | 22,334 / 52,389 |
| `COAD_Lee2020` | COAD (Colorectal) | Lee et al. 2020 | primary | 15 (15) | 20,597 | 11,056 / 9,541 |
| `GBM_Wang2019` | GBM (Brain) | Wang et al. 2019 | primary | 9 (7) | 15,630 | 9,408 / 6,222 |
| `HNSC_Kurten2021` | HNSC (Head and Neck) | Kürten et al. 2021 | primary | 18 (14) | 54,127 | 11,390 / 42,737 |
| `KIRC_Zhang2021` | KIRC (Kidney) | Zhang et al. 2021 | primary | 9 (8) | 25,220 | 9,934 / 15,286 |
| `LIHC_Ma2019` | LIHC (Liver) | Ma et al. 2019 | primary | 8 (6) | 2,596 | 1,609 / 987 |
| `LUAD_Kim2020` | LUAD (Lung) | Kim et al. 2020 | primary | 14 (14) | 29,196 | 9,308 / 19,888 |
| `OV_Zhang2022` | OV (Ovary) | Zhang et al. 2022 | mixed | 11 (11) | 41,235 | 8,772 / 32,463 |
| `PDAC_Peng2019` | PDAC (Pancreas) | Peng et al. 2019 | primary | 16 (16) | 26,419 | 7,930 / 18,489 |
| `SKCM_Biermann2022` | SKCM (Skin) | Biermann et al. 2022 | metastatic | 31 (27) | 99,954 | 59,822 / 40,132 |
| Total | | | | 157 (138) | 389,697 | 151,563 / 238,134 |

Full citations and download URLs are in `configs/datasets.yaml`; `scripts/download_data.py`
records checksums and download dates in `data/raw/manifest.json`. The registry does not
record licences: the datasets keep the licences and terms of their original publishers (3CA,
TISCH2 and the individual studies), and `data/app` contains summary tables derived from them.

## Citation

scTCNA has no publication yet (manuscript in preparation). Please cite TCNA and the data
sources:

- Shah R, Faisal A, bin Faiz A, Bashir I, Jawaid A, Bashir K, Shoaib M, Chaudhary SU.
  The Cancer Noise Atlas (TCNA): A Webserver of Multiscale Transcriptional Noise in Cancer.
  *J Mol Biol* 2026.
- Gavish A, Tyler M, Greenwald AC, et al. Hallmarks of transcriptional intratumour
  heterogeneity across a thousand tumours. *Nature* 618, 598-606 (2023).
- The individual studies in the Data table (citations in `configs/datasets.yaml`).
- DEPTH: Li M, Zhang Z, Li L, Wang X. *Commun Biol* 3, 505 (2020). DEPTH2: Song D, Wang X.
  *J Transl Med* 20, 150 (2022).

## License

MIT, see [LICENSE](LICENSE).

## Author

Saba Siddique, Biomedical Informatics & Engineering Research Laboratory (BIRL), Lahore
University of Management Sciences (LUMS). Contact: Dr. Safee Ullah Chaudhary, BIRL, LUMS.
