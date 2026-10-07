---
title: scTCNA
colorFrom: indigo
colorTo: blue
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: Single-cell transcriptional noise atlas for cancer
---

# scTCNA

A single-cell RNA-seq counterpart of The Cancer Noise Atlas (TCNA). For each
cancer type the atlas reports per-gene expression noise in malignant and
non-malignant cells, computed per patient and summarised as medians across
patients:

- mean (CP10K), SD, CV = SD/mean, CV^2, MAD = mean absolute deviation,
  Fano = variance/mean, and mean-corrected CV^2 (LOWESS residual);
- differential noise log2(CV_mal/CV_non-mal) with paired Wilcoxon tests across
  patients and Benjamini-Hochberg correction per cancer type;
- DEPTH and DEPTH2 tumour-heterogeneity scores per cell and per patient
  pseudo-bulk;
- the same metrics after downsampling all cells to a common UMI depth.

Compartments are "Malignant" and "Non-malignant", the single-cell stand-in for
TCNA's Tumor vs Normal comparison.

This Space serves precomputed tables; nothing is recomputed at request time.
The pipeline that builds them (QC, doublet removal, CP10K normalisation,
per-patient metrics, downsampling, DEPTH/DEPTH2 and the paired statistics),
the dataset registry and the full documentation live in the project
repository. Datasets come from the Curated Cancer Cell Atlas (3CA; Gavish et
al., Nature 2023) and are listed, with citations, on the app's datasets page.

Built from a Docker image (`python:3.12-slim`, FastAPI + uvicorn, Parquet via
pyarrow). Health check: `/api/health`.
