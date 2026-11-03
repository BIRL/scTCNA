# scTCNA results

All ten cancer types were processed by the same pipeline code (`sctcna/` as of commit
`ec1e423`); the numbers below are produced by `python scripts/summarise_results.py`
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
