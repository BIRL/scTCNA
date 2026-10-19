# scTCNA: a single-cell Cancer Noise Atlas

scTCNA is a single-cell RNA-seq counterpart of **The Cancer Noise Atlas (TCNA)**
(Shah R, Faisal A, bin Faiz A, Bashir I, Jawaid A, Bashir K, Shoaib M, Chaudhary SU;
*J Mol Biol* 2026; portal <https://tcna.lums.edu.pk/>, code <https://github.com/BIRL/TCNA>).

TCNA quantifies transcriptional noise from bulk RNA-seq (GDC, tumour vs solid tissue
normal) at the gene, pathway and tumour scales. scTCNA keeps TCNA's metric definitions
and dashboards but replaces bulk samples by single cells and the tumour-vs-normal
contrast by **malignant vs non-malignant cells of the same patients**, using curated
10x Genomics studies from the Curated Cancer Cell Atlas (3CA).

## What it is

- A **pipeline** (`sctcna/`, scanpy-based) that turns a 3CA (or TISCH2) study bundle
  into per-gene noise metrics (mean, SD, CV, CV^2, MAD, Fano, detection rate,
  mean-corrected CV^2), differential noise between compartments with patient-level
  statistics, and tumour-heterogeneity scores (DEPTH, DEPTH2) per cell and per
  pseudo-bulk profile.
- A set of **precomputed Parquet tables** (`data/app/`) that summarise ten cancer types.
- A **web app** (`app/`, FastAPI + static frontend) with the same dashboards as TCNA's
  Gene Noise Analysis and Tumor Heterogeneity pages, plus a Cross-cancer summary
  (is malignant-cell noise higher, patient by patient?), a Top genes view, a
  cross-cancer heatmap and an About/Methods page.

Metric vocabulary used throughout the app: mean (CP10K), SD, CV = SD/mean, CV^2,
MAD = mean absolute deviation, Fano = variance/mean, mean-corrected CV^2 (LOWESS
residual), differential noise log2(CV_mal/CV_non-mal), DEPTH and DEPTH2. Compartments
are **Malignant** and **Non-malignant**, the single-cell stand-ins for TCNA's Tumor and
Normal.

## How it differs from TCNA

| | TCNA | scTCNA |
|---|---|---|
| Data | Bulk RNA-seq from GDC: Primary Tumor vs Solid Tissue Normal, 26 anatomical sites | 10x scRNA-seq from 3CA: 10 cancer types, one study per cancer type |
| Unit of observation | Sample | Cell, grouped by patient and compartment |
| Contrast | Tumour vs normal samples | Malignant vs non-malignant cells within the same patients |
| Expression scale for gene noise | TPM / FPKM / FPKM-UQ (user-selectable) | CP10K (counts per 10,000) |
| Expression scale for DEPTH / DEPTH2 | log2(x + 1) | log2(CP10K + 1) per cell; log2(CPM + 1) per pseudo-bulk |
| Gene noise metrics | CV, CV^2, SD, MAD | CV, CV^2, SD, MAD, plus Fano, detection rate and mean-corrected CV^2 |
| Differential noise | log2(CV_tumour / CV_normal), Eq. 1 | Eq. 1 evaluated per patient and summarised as the median across patients; also for MAD, Fano and the difference of mean-corrected CV^2 |
| Statistics | Differential noise shown as logFC; no per-gene test described | Paired Wilcoxon signed-rank across patients (>= 6), rank-biserial effect size, BH correction per cancer type; patient-level malignant vs reference contrasts |
| Robustness checks | none | Downsampling to a common UMI depth; depth-matched, mean-corrected excess noise; pooled-cells vs per-patient CV diagnostic |
| Tumour heterogeneity | DEPTH and DEPTH2 per sample | DEPTH (malignant cells vs the non-malignant mean) and DEPTH2 (z-scores within each patient x compartment) per cell, plus both per patient x compartment pseudo-bulk |
| Lineage resolution | none | Per cell type (patient x cell type) metrics, so malignant cells can be compared with the matched normal lineage |
| Pathway analysis | STRINGdb enrichment dashboard | Not replicated |
| Custom data upload | Supported | Not supported |
| Computation | On demand (FastAPI + MariaDB) | Precomputed Parquet tables served by FastAPI |
| Cross-cancer view | Differential noise across sites | Heatmap of within-dataset values (ranks, z-scores); no cross-dataset integration |

## Data

Primary source: **3CA**, the Curated Cancer Cell Atlas (Gavish et al., *Nature* 2023;
<https://www.weizmann.ac.il/sites/3CA/>). Each study bundle (hosted on Dropbox) ships
`Exp_data_UMIcounts.mtx`, `Genes.txt`, `Cells.csv` and `Samples.csv`. 3CA labels
malignant cells as `Malignant` and gives every other cell a curated non-malignant
lineage; cells with an empty `cell_type` are unassigned and are excluded.

Ten 10x Genomics studies are enabled in `configs/datasets.yaml` (patient and cell
counts are those reported by 3CA, before scTCNA's per-patient cap and QC). Sun et al.
2021 (liver) was planned but 3CA distributes it only as normalised TPM values, not UMI
counts; the pipeline refuses non-integer input, so Ma et al. 2019 is used instead.

| Dataset id | Cancer type (TCNA site) | Study | Patients | Cells |
|---|---|---|---|---|
| `BRCA_Wu2021` | BRCA (Breast) | Wu et al. 2021, *Nat Genet* | 26 | 100,064 |
| `LUAD_Kim2020` | LUAD (Lung) | Kim et al. 2020, *Nat Commun* | 14 | 32,493 |
| `COAD_Lee2020` | COAD (Colorectal) | Lee et al. 2020, *Nat Genet* | 15 | 21,657 |
| `HNSC_Kurten2021` | HNSC (Head and Neck) | Kürten et al. 2021, *Nat Commun* | 18 | 62,700 |
| `GBM_Wang2019` | GBM (Brain) | Wang et al. 2019, *Cancer Discov* | 9 | 16,317 |
| `LIHC_Ma2019` | LIHC (Liver: HCC and iCCA) | Ma et al. 2019, *Cancer Cell* | 8 | 3,018 |
| `PDAC_Peng2019` | PDAC (Pancreas) | Peng et al. 2019, *Cell Res* | 16 | 27,953 |
| `SKCM_Biermann2022` | SKCM (Skin) | Biermann et al. 2022, *Cell* | 31 | 136,973 (all metastatic) |
| `OV_Zhang2022` | OV (Ovary) | Zhang et al. 2022, *Sci Adv* | 11 | 51,786 |
| `KIRC_Zhang2021` | KIRC (Kidney) | Zhang et al. 2021, *PNAS* | 9 | 29,474 |

Exact download URLs, SHA-256 checksums and download dates are written by
`scripts/download_data.py` to `data/raw/manifest.json` and surfaced in the app's
datasets table.

Secondary source: **TISCH2** (Han et al., *NAR* 2023; <https://tisch.compbio.cn/>).
TISCH2 distributes log1p(CP10K) matrices rather than raw counts; the loader inverts the
per-cell scaling to recover integer UMI counts and verifies integrality. Two TISCH2
entries (BRCA and PDAC re-processings of the same studies) are in the registry but
disabled.

Note on availability: at development time the original 3CA website (weizmann.ac.il) and
the old TISCH2 domain were not reachable from the development machine; the Dropbox
bundles and tisch.compbio.cn were.

## Results summary

All ten cancer types were processed by the same pipeline code (`sctcna/` as of commit
`a01a0d1`); the numbers below are produced by `python scripts/summarise_results.py`
from the committed `data/app` tables. In total 389,697 cells (151,563 malignant,
238,134 non-malignant) from 157 patients passed QC; 138 patients have both compartments.

**Per-gene differential noise and patient-level contrasts against all non-malignant
cells.** "q < 0.05" counts expressed genes whose per-patient log2(CV_mal/CV_non-mal) is
significant (paired Wilcoxon, BH per cancer type). The last three columns are medians over
patients of the median over genes (Cross-cancer summary page); negative values mean
malignant cells are *less* noisy.

| Cancer | Study | Patients (paired) | Cells mal. / non-mal. | Genes tested | q < 0.05 (noisier mal. / non-mal.) | Median log2 CV ratio, raw | Depth-matched | Excess noise, log10 CV² (q) |
|---|---|---|---|---|---|---|---|---|
| BRCA | Wu et al. 2021 | 26 (20) | 22,334 / 52,389 | 2,919 | 2,456 (79 / 2,377) | -0.58 | -0.03 | -0.032 (2.9e-05) |
| COAD | Lee et al. 2020 | 15 (15) | 11,056 / 9,541 | 4,412 | 4,099 (96 / 4,003) | -0.81 | -0.05 | -0.051 (4.6e-04) |
| GBM | Wang et al. 2019 | 9 (7) | 9,408 / 6,222 | 5,917 | 0 (0 / 0) | -0.10 | 0.02 | -0.011 (3.8e-01) |
| HNSC | Kürten et al. 2021 | 18 (14) | 11,390 / 42,737 | 3,471 | 2,761 (77 / 2,684) | -0.58 | -0.03 | -0.037 (8.5e-04) |
| KIRC | Zhang et al. 2021 | 9 (8) | 9,934 / 15,286 | 5,392 | 954 (70 / 884) | -0.31 | -0.01 | -0.021 (3.1e-02) |
| LIHC | Ma et al. 2019 | 8 (6) | 1,609 / 987 | 4,016 | 0 (0 / 0) | -0.21 | 0.12 | -0.020 (4.7e-02) |
| LUAD | Kim et al. 2020 | 14 (14) | 9,308 / 19,888 | 4,628 | 3,715 (92 / 3,623) | -0.52 | -0.08 | -0.039 (1.3e-03) |
| OV | Zhang et al. 2022 | 11 (11) | 8,772 / 32,463 | 3,729 | 3,444 (37 / 3,407) | -0.99 | 0.02 | -0.037 (4.7e-03) |
| PDAC | Peng et al. 2019 | 16 (16) | 7,930 / 18,489 | 6,888 | 5,093 (186 / 4,907) | -0.30 | -0.04 | -0.028 (3.2e-04) |
| SKCM | Biermann et al. 2022 | 31 (27) | 59,822 / 40,132 | 5,033 | 4,290 (177 / 4,113) | -0.42 | 0.02 | -0.056 (1.3e-07) |

**Tumour heterogeneity (DEPTH2).** Per cell, DEPTH2 measures how atypical a cell is within
its own patient x compartment group; per pseudo-bulk, how much a patient's compartment
profile differs from the other patients'. Median paired difference malignant minus
non-malignant, with the Wilcoxon p-value. The last column is the Spearman correlation
between raw and depth-downsampled per-patient CVs (all genes and groups).

| Cancer | DEPTH2 per cell, median diff. mal. - non-mal. (p) | DEPTH2 pseudo-bulk, median diff. (p) | Spearman CV vs downsampled CV |
|---|---|---|---|
| BRCA | -0.075 (3.8e-06) | 0.123 (5.7e-06) | 0.92 |
| COAD | -0.107 (6.1e-05) | 0.062 (1.2e-01) | 0.91 |
| GBM | -0.034 (3.8e-01) | 0.027 (3.0e-01) | 0.95 |
| HNSC | -0.097 (1.2e-04) | 0.195 (2.4e-04) | 0.96 |
| KIRC | -0.020 (4.6e-01) | 0.132 (7.8e-02) | 0.96 |
| LIHC | 0.033 (2.2e-01) | -0.058 (1.6e-01) | 0.91 |
| LUAD | -0.084 (1.2e-03) | 0.177 (2.4e-04) | 0.95 |
| OV | -0.148 (9.8e-04) | 0.139 (9.8e-04) | 0.92 |
| PDAC | -0.039 (2.1e-04) | 0.141 (3.1e-05) | 0.96 |
| SKCM | -0.022 (1.2e-01) | 0.021 (1.3e-01) | 0.96 |

**Malignant cells against single non-malignant lineages** with >= 6 patients (median
depth-matched excess noise of malignant cells, log10 CV²; q):
- BRCA: B_cell -0.019 (8.8e-02); Endothelial -0.004 (1.5e-01); Fibroblast -0.005 (3.1e-02); Macrophage -0.006 (5.2e-02); Monocyte -0.000 (1.0e+00); NK_cell -0.003 (6.2e-01); Pericyte -0.006 (1.0e-01); Plasmablast -0.037 (8.8e-02); T_cell -0.018 (1.4e-03)
- COAD: B_cell -0.046 (9.8e-03); Fibroblast -0.006 (9.8e-03); Macrophage -0.009 (9.8e-03); T_cell -0.039 (9.2e-04)
- HNSC: B_cell -0.013 (4.8e-02); Endothelial +0.002 (6.3e-01); Fibroblast -0.002 (3.3e-01); Lymphovascular -0.005 (6.3e-01); Macrophage -0.017 (4.8e-02); T_cell -0.033 (1.9e-03)
- KIRC: Endothelial -0.002 (5.6e-01); Macrophage -0.016 (3.1e-02); Myocyte -0.004 (4.0e-01)
- LUAD: B_cell -0.026 (2.5e-02); Epithelial -0.005 (3.6e-01); Fibroblast -0.004 (4.2e-01); Macrophage -0.010 (5.1e-03); Mast +0.003 (4.6e-01); T_cell -0.027 (4.1e-03)
- OV: B_cell +0.001 (1.0e+00); Dendritic -0.004 (5.9e-02); Fibroblast -0.011 (1.7e-02); Macrophage -0.016 (6.7e-03); Mesothelial -0.002 (7.3e-02); NK_cell -0.016 (6.6e-02); T_cell -0.021 (4.7e-03)
- PDAC: Endothelial -0.002 (9.0e-01); Epithelial -0.002 (2.9e-01); Fibroblast -0.002 (6.7e-01); Macrophage -0.011 (2.7e-02); Stellate -0.007 (5.8e-01); T_cell -0.015 (2.7e-02)
- SKCM: B_cell -0.018 (2.8e-01); Fibroblast -0.022 (4.9e-02); Myeloid -0.029 (1.4e-06); Plasma -0.014 (6.3e-01); T_NK_cell -0.045 (2.3e-05)

### Key findings

1. With TCNA's raw measure (Eq. 1 per patient), malignant cells look less noisy than the
   pooled non-malignant cells in all ten cancer types (significant in nine); thousands of
   genes per cancer reach q < 0.05, almost all in that direction.
2. Most of that gap appears to track sequencing depth. Malignant cells have more UMIs per
   cell than non-malignant cells in 134 of 138 patients, and after downsampling every cell
   to the same depth the median CV ratio lies within ±0.1 in nine of ten cancer types
   (LIHC +0.12 with 6 patients; still significant in LUAD and COAD). Whether the higher
   depth is technical or reflects more mRNA per malignant cell cannot be told apart here.
3. Controlling for both depth and mean expression, malignant cells are slightly *less*
   variable than non-malignant cells at the same expression level in nine of ten cancer
   types (q < 0.05; median excess -0.01 to -0.06 log10 CV², GBM not significant). Against
   single lineages the deficit is largest for T, B and myeloid cells and near zero for
   fibroblasts, endothelial and normal epithelial cells.
4. Per-cell DEPTH2 is lower in malignant cells than in the patient's pooled non-malignant
   cells in nine of ten cancer types (p < 0.05 in six); the non-malignant group mixes several
   lineages, which raises its within-group heterogeneity. Pseudo-bulk DEPTH2 is higher for
   malignant profiles in nine of ten (p < 0.05 in five), i.e. malignant profiles differ more
   between patients. This is consistent with between-patient differences contributing much
   of the tumour heterogeneity seen in bulk data, but scTCNA was not compared with TCNA's
   bulk scores directly.
5. Genes most often noisier in malignant cells across cancers include MALAT1, HLA-B, VIM,
   S100A4, ARHGDIB and DUSP2 (q < 0.05 in five to seven cancer types).

## Repository layout

```
README.md, METHODS.md          this file and the methods description
pyproject.toml, requirements.txt
configs/datasets.yaml          dataset registry and every pipeline parameter
sctcna/                        pipeline package
  config.py                    registry loading
  io.py                        3CA / TISCH2 loaders -> AnnData (raw UMI counts)
  qc.py                        cell/gene filters, mito fraction, Scrublet, depth downsampling
  normalize.py                 CP10K and log2(x + 1)
  metrics.py                   per-gene noise metrics, LOWESS residual, Eq. 1
  depth.py                     DEPTH and DEPTH2 (per cell and per pseudo-bulk)
  stats.py                     patient-level medians, paired Wilcoxon, BH, rank-biserial
  pipeline.py                  process_dataset() and aggregate()
  synthetic.py                 simulated data with known noise structure
  cli.py                       python -m sctcna run | aggregate
scripts/
  download_data.py             fetch 3CA / TISCH2 bundles, write data/raw/manifest.json
  check_memory.py              peak memory of the app under a worst-case request mix
app/
  main.py                      FastAPI app (create_app; reads SCTCNA_DATA_DIR)
  static/                      frontend pages and local assets (no CDN at runtime)
tests/                         pytest suite (unit, API, integration) and Playwright e2e
docs/                          TCNA manuscript (docx + extracted text), DEPLOY.md, memory_check.txt
Dockerfile, render.yaml        container image and Render Blueprint used for deployment
data/
  raw/                         downloaded bundles + manifest.json (git-ignored)
  processed/<dataset>/         per-dataset pipeline outputs (git-ignored)
  app/                         aggregated tables served by the app (committed)
```

## Install

Python 3.12. Dependencies are pinned in `requirements.txt`.

```bash
uv venv --python 3.12 .venv && source .venv/bin/activate && uv pip install -r requirements.txt && python -m playwright install chromium
```

Without `uv`:

```bash
python3.12 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt && python -m playwright install chromium
```

The Playwright Chromium install is only needed for the end-to-end browser tests.

## Download data

```bash
python scripts/download_data.py --list                                   # show the registry
python scripts/download_data.py --all --convert                          # every enabled dataset
python scripts/download_data.py --datasets LUAD_Kim2020 COAD_Lee2020 --convert
```

Bundles land in `data/raw/<source>/<dataset_id>/`. `--convert` writes
`<dataset_id>.h5ad` (raw UMI counts with the curated annotations, already capped at
3,000 cells per patient x compartment). The multi-gigabyte expression matrix is streamed
straight out of the `.tar.gz` and never extracted; the archive is deleted after a
successful conversion unless `--keep-archives` is given. The pipeline reads the `.h5ad`
when it exists. The 3CA archives are 0.2-1.5 GB each; conversion needs roughly twice
the archive size of free disk at a time.

## Run the pipeline

```bash
python -m sctcna run                                  # all enabled datasets -> data/processed/<id>/
python -m sctcna run --datasets LUAD_Kim2020          # one or more dataset ids
python -m sctcna run --synthetic --out data/processed_syn   # simulated data, no download needed
python -m sctcna aggregate                            # data/processed -> data/app
```

`run` writes per-dataset Parquet tables plus a `summary.json` under
`data/processed/<id>/`; `aggregate` combines them into the compact tables the app reads
(`datasets`, `gene_summary`, `cell_scores`, `sample_scores`, `celltype_summary`,
`patient/<id>` and `build.json`). Useful flags: `run --config PATH --raw DIR --out DIR
--fail-fast`, `aggregate --processed DIR --app DIR --include-synthetic` (adds the
synthetic dataset), and `python -m sctcna -v run` for verbose logging, for example

```bash
python -m sctcna aggregate --processed data/processed_syn --app data/app_syn --include-synthetic
```

All thresholds live in the `params` section of `configs/datasets.yaml`.

## Run the app

```bash
SCTCNA_DATA_DIR=data/app python -m uvicorn app.main:app --reload
```

Then open <http://127.0.0.1:8000>. `SCTCNA_DATA_DIR` points at an aggregated data
directory (for example `data/app_syn` for a quick look at synthetic data).

Pages: Gene Noise Analysis (pan-cancer and cancer-specific modes), Cross-cancer summary,
Top genes, Tumor Heterogeneity (DEPTH / DEPTH2), Heatmap, About / Methods.

API endpoints:

| Endpoint | Returns |
|---|---|
| `GET /api/health` | build info and loaded datasets |
| `GET /api/cancer-types` | cancer types with dataset metadata |
| `GET /api/datasets` | the datasets table |
| `GET /api/genes/search?q=` | gene name autocomplete |
| `GET /api/genes/{gene}` | per-cancer-type summary for one gene |
| `GET /api/genes/{gene}/patients` | per-patient values for one gene |
| `GET /api/genes/{gene}/celltypes` | per-cell-type values for one gene |
| `GET /api/cancer-types/{ct}/top-genes` | ranked differential-noise genes |
| `GET /api/cancer-types/{ct}/comparison` | malignant vs non-malignant comparison for a cancer type |
| `GET /api/heatmap` | cross-cancer matrix for the heatmap |
| `GET /api/cancer-types/{ct}/tumor-heterogeneity` | DEPTH / DEPTH2 per cell and per pseudo-bulk |
| `GET /api/tumor-heterogeneity/summary` | per-cancer-type DEPTH / DEPTH2 summary |
| `GET /api/noise-contrasts` | malignant vs reference-group contrasts for every cancer type |
| `GET /api/cancer-types/{ct}/noise-contrasts` | the same with per-patient values |
| `GET /api/download/{table}.csv` | CSV export of a served table |

### Embedding the API in another FastAPI application

The API is also importable on its own, without the pipeline stack. `sctcna.web.create_api`
returns a FastAPI application that another server can mount at any prefix; the TCNA portal
serves it at `/api/sc`:

```python
from sctcna.web import create_api
portal.mount("/api/sc", create_api("/path/to/data/app"))
```

It needs only pandas (2.2 or 3.x), numpy, pyarrow and FastAPI, so it can be installed into
another environment with `pip install --no-deps <path to this repository>`.

## Run the tests

```bash
pytest -m "not e2e"                 # unit, API and pipeline integration tests (synthetic data)
pytest -m e2e tests/test_e2e.py     # Playwright browser smoke tests (needs Playwright Chromium)
ruff check . && ruff format --check .
python scripts/check_memory.py --data data/app   # peak memory of the app, must stay under 512 MB
```

Plain `pytest` runs both sets.

## Deploy

The app is deployed as one free Docker web service on **Render** (512 MB RAM, no card
needed); `render.yaml` is the Blueprint. Push the repository to GitHub, then on
<https://dashboard.render.com> choose **New > Blueprint**, pick the repository and click
**Apply**. Step-by-step instructions, the memory measurement and troubleshooting are in
[docs/DEPLOY.md](docs/DEPLOY.md).

Locally:

```bash
docker build -t sctcna . && docker run --rm -p 7860:7860 sctcna
```

## Methods

The full description of data handling, QC, normalisation, metric formulas, statistics,
DEPTH / DEPTH2, output tables and the list of what is taken from the TCNA manuscript
versus what is an scTCNA assumption is in [METHODS.md](METHODS.md). In short: every
metric is computed per gene within one patient x compartment group of cells; the patient
is the unit of replication; cancer-type values are medians across patients; differential
noise is TCNA Eq. 1 per patient tested with a paired Wilcoxon signed-rank test and BH
correction per cancer type.

Known limitations (detailed in METHODS.md): one study per cancer type; the non-malignant
compartment is a mixture of lineages; malignant calls come from the atlas and are not
re-derived; the melanoma study is metastatic; no pathway module and no custom upload;
10x UMI data only; cross-cancer comparisons are confounded by study effects; patient
numbers (6-27 paired patients per cancer type) limit the power of the paired
tests.

## Citation and acknowledgements

If you use scTCNA, please cite TCNA and the data sources:

- Shah R, Faisal A, bin Faiz A, Bashir I, Jawaid A, Bashir K, Shoaib M, Chaudhary SU.
  The Cancer Noise Atlas (TCNA): A Webserver of Multiscale Transcriptional Noise in
  Cancer. *J Mol Biol* 2026. <https://tcna.lums.edu.pk/>
- Gavish A, Tyler M, Greenwald AC, et al. Hallmarks of transcriptional intratumour
  heterogeneity across a thousand tumours. *Nature* 618, 598-606 (2023). (3CA)
- Han Y, Wang Y, Dong X, et al. TISCH2: expanded datasets and new tools for single-cell
  transcriptome analyses of the tumor microenvironment. *Nucleic Acids Res* 51,
  D1425-D1431 (2023).
- The individual studies listed in the Data table (full citations in
  `configs/datasets.yaml`).
- DEPTH: Li et al. 2020 (<https://github.com/WangX-Lab/DEPTH>); DEPTH2: Song & Wang
  2022 (<https://github.com/WangX-Lab/DEPTH2>).

The pipeline is built on scanpy, anndata, Scrublet, statsmodels, SciPy, pandas and
pyarrow; the app on FastAPI and uvicorn.

## Licence

MIT, see [LICENSE](LICENSE). The underlying datasets keep the licences and terms of their
original publishers (3CA, TISCH2 and the individual studies).
