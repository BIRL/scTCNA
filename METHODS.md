# scTCNA methods

scTCNA re-implements the noise metrics of The Cancer Noise Atlas (TCNA; Shah et al.,
*J Mol Biol* 2026) on single-cell RNA-seq. Where the TCNA manuscript or code defines a
quantity it is used unchanged; every single-cell specific choice is an scTCNA assumption
and is listed as such in the table near the end. Parameter values are the defaults in
`configs/datasets.yaml`.

Notation: a *group* is one patient x compartment set of cells; `x_ig` is the normalised
expression of gene `g` in cell `i`; `n` is the number of cells in the group.

## Data

**3CA** (Curated Cancer Cell Atlas; Gavish et al. 2023). A study bundle holds
`Exp_data_UMIcounts.mtx` (genes x cells, raw UMI counts), `Genes.txt`, `Cells.csv` and
`Samples.csv`; the matching Meta-data bundle carries the curated cell-type calls and is
preferred when present. Patients are taken from the `patient` column of `Cells.csv`,
otherwise from the sample-to-patient map in `Samples.csv`, otherwise the sample id is
used. Compartment: `cell_type == "Malignant"` is **malignant**; every other non-empty
`cell_type` (the curated non-malignant lineages) is pooled into **non-malignant**; cells
with an empty `cell_type` were left unassigned by 3CA and are excluded. Malignant calls
are taken from the atlas (3CA derives them from inferred copy-number aberrations) and are
not re-derived. All ten enabled datasets are 10x Genomics UMI data; the loader warns if a
bundle only provides TPM, and the pipeline stops if the matrix is not integer counts.

**TISCH2** (Han et al. 2023) is supported as a second source. Its expression `.h5` files
are log1p(CP10K); for each cell the smallest non-zero value corresponds to a count of 1,
which gives the cell total and hence every count. The fraction of cells whose recovered
values are all within 0.02 of an integer is recorded and a warning is raised below 95%.
`Celltype (malignancy) == "Malignant cells"` defines the malignant compartment and
`Celltype (major-lineage)` the cell type. The two TISCH2 registry entries are disabled.

**Cell cap.** At load, every patient x compartment group larger than 3,000 cells is
randomly subsampled to 3,000 (without replacement, seed 0). The cap is applied before QC.

Downloads are recorded in `data/raw/manifest.json` (URL, size, SHA-256, UTC date).

## Quality control

Applied uniformly after the atlas's own processing, in this order:

1. Cells: >= 200 detected genes, >= 500 UMIs, mitochondrial fraction <= 20% (genes with
   the `MT-` prefix).
2. Genes: detected (count > 0) in >= 10 cells.
3. Doublets: Scrublet per sample (scanpy's `pp.scrublet`, seed 0); predicted doublets
   are removed. Samples with fewer than 100 cells are skipped and kept, since the kNN
   simulation is unreliable there.
4. Patient x compartment groups with fewer than 20 cells are dropped.

Per-cell `n_counts`, `n_genes`, `pct_mito` and `doublet_score` are kept and exported.

## Normalisation

- Gene noise metrics: **CP10K**, each cell scaled to 10,000 total counts (linear scale).
  This is the single-cell stand-in for TCNA's TPM / FPKM.
- DEPTH and DEPTH2: **log2(CP10K + 1)** per cell, mirroring TCNA's log2(x + 1) before the
  heterogeneity metrics. Pseudo-bulk profiles use log2(CPM + 1) of summed counts (below).

No global low-expression filter is applied to the metrics themselves (as in TCNA, the
mean is reported next to every noise value); a detection-rate threshold is used only to
decide which genes enter tests and top lists.

## Noise metrics

Computed per gene within every group of cells, with `ddof = 1` for the variance.
Undefined values (mean = 0) are NaN rather than 0.

```
mean_g   = (1/n) * sum_i x_ig
SD_g     = sqrt( (1/(n-1)) * sum_i (x_ig - mean_g)^2 )
CV_g     = SD_g / mean_g
CV2_g    = CV_g^2
MAD_g    = (1/n) * sum_i |x_ig - mean_g|          (mean absolute deviation, TCNA definition)
Fano_g   = SD_g^2 / mean_g                        (scTCNA addition; on the CP10K scale, see below)
det_g    = (1/n) * sum_i [x_ig > 0]               (detection rate, scTCNA addition)
```

TCNA's web code reports CV x 100; scTCNA stores the plain ratio. Because all metrics are
computed on CP10K values, they depend on sequencing depth: for Poisson counts at library
size N the Fano factor is about 10^4 / N (not 1), and CV^2 has a Poisson floor of about
10^4 / (N x mean). Compartments sequenced at different depths therefore differ in raw CV,
MAD and Fano for technical reasons; the downsampled metrics and the depth-matched
contrasts below remove that difference.

**Differential noise** is TCNA Eq. 1, `logFC_g = log2(CV_g^tumour / CV_g^normal)`,
evaluated per patient `p`:

```
D_gp = log2( CV_gp^malignant / CV_gp^non-malignant )
```

and summarised per cancer type as the median of `D_gp` across patients
(`log2fc_cv_median`; the mean is stored as `log2fc_cv_mean`). The same ratio is computed
for MAD (`log2fc_mad_median`) and Fano (`log2fc_fano_median`), each with its own paired
test and q-values. For the mean-corrected CV^2 the per-patient difference of residuals is
used (`delta_cv2_resid` = median over patients of `cv2_resid_malignant - cv2_resid_non_malignant`,
tested in the same paired way: `p_value_cv2_resid`, `q_value_cv2_resid`), since residuals
are already on a log scale and can be negative. The residuals are relative to each group's
own trend, so their median over genes is close to zero by construction: they rank genes
within a group but do not say which compartment is noisier overall (see the contrasts
below for that). For display next to the per-patient version,
Eq. 1 is also evaluated on the cancer-type (median across patients) CVs
(`log2fc_cv_of_medians`).

Metrics are computed for three kinds of groups: patient x compartment (the primary
unit), compartment with all patients' cells pooled (`patient = "ALL"`, a diagnostic, see
Statistics) and patient x cell type (lineage level, groups with >= 50 cells).

## Mean-dependence correction

In single-cell data CV^2 falls steeply with the mean (sampling noise dominates lowly
expressed genes). Within each group a LOWESS curve `f` (statsmodels, `frac = 0.3`, two
robustifying iterations) is fitted to `log10(CV2_g)` against `log10(mean_g)` across all
genes with `mean > 0` and `CV2 > 0`, and the residual

```
cv2_resid_g = log10(CV2_g) - f(log10(mean_g))
```

is the **mean-corrected CV^2 (LOWESS residual)**: positive values are noisier than genes
of the same mean, negative values quieter. The fit needs at least 50 usable genes;
otherwise the residual is NaN for the group. This follows the mean-variance trend
approach of Brennecke et al. (2013) and is an scTCNA addition.

## Depth robustness

CV and MAD also depend on sequencing depth through dropout. As a check, every cell with
more than a common target depth is downsampled (without replacement, seed 0) to exactly
that depth; the target is the 10th percentile of per-cell total UMIs after QC, with a
floor of 1,000. Cells below the target are left unchanged. Counts are re-normalised to
CP10K and the per-patient metrics recomputed (`mean_ds`, `cv_ds`, `mad_ds`, `fano_ds`,
`cv2_resid_ds`). The Spearman correlation between CV and downsampled CV across all gene
x patient x compartment rows is reported per dataset (`spearman_cv_vs_downsampled`; the
same for the LOWESS residual is kept in `summary.json`).

## Replication unit and statistics

The **patient** is the unit of replication; cells within a patient are not independent
observations of a cancer type.

- Cancer-type values of every metric are **medians across patients** of the per-patient
  values, per compartment (`<metric>_malignant`, `<metric>_non_malignant`), together with
  the number of contributing patients.
- **Paired test.** For every gene with at least 6 patients having both compartments, the
  per-patient `D_gp` are tested with a two-sided Wilcoxon signed-rank test (zero
  differences discarded; if all differences are zero, p = 1). Effect sizes are the median
  per-patient log2 ratio and the matched-pairs rank-biserial correlation
  `r = (R+ - R-) / (R+ + R-)`, where `R+` and `R-` are the rank sums of positive and
  negative differences (Kerby 2014; +1 when every patient is noisier in malignant cells).
- **Multiple testing.** Only expressed genes (detected in >= 10% of cells in both
  compartments, median across patients) enter the paired tests; Benjamini-Hochberg
  q-values are computed per cancer type over that family of tested genes.
- **Expressed genes.** A gene is `expressed` when its detection rate (median across
  patients) is >= 10% in both compartments. Only expressed genes are counted in the
  per-dataset significance tallies, listed in top-gene views and kept in the per-patient
  app tables.
- **Pseudoreplication diagnostic.** CVs computed on all patients' cells pooled per
  compartment (`cv_pooled_malignant`, `cv_pooled_non_malignant`) are stored next to the
  median-of-patients CVs. The pooled value treats cells as independent and absorbs
  between-patient differences; the gap between the two shows how much apparent noise is
  inter-patient rather than intra-tumour.
- **Heterogeneity scores.** Per-patient pseudo-bulk DEPTH2 and per-patient medians of
  the per-cell DEPTH2 are compared between compartments with the same paired Wilcoxon
  test (>= 6 patients); the median difference and p-value are reported per dataset. DEPTH
  is not compared, because only malignant profiles have a DEPTH score.

## Malignant versus reference noise contrasts

The per-gene tests ask which genes differ. To ask whether malignant cells are noisier
overall, every patient `p` is compared with a reference group `R` (all non-malignant
cells of that patient, or one non-malignant cell type with >= 50 cells). Over the genes
`G_pR` detected in >= 10% of cells of both groups:

```
log2fc_cv_pR     = median_{g in G_pR} log2( CV_gp^mal / CV_gp^R )                  (TCNA Eq. 1)
log2fc_cv_ds_pR  = the same with CVs after downsampling to the common UMI depth
excess_cv2_ds_pR = median_{g in G_pR} [ log10 CV2_gp^mal,ds - f_pR( log10 mean_gp^mal,ds ) ]
```

where `f_pR` is a LOWESS curve (frac 0.3) of log10 CV^2 on log10 mean fitted across all
genes of the reference group after downsampling; malignant genes whose mean falls outside
the reference's fitted range are skipped. The excess therefore compares malignant genes
with reference genes of the same expression level at the same sequencing depth. Each of
the three per-patient values is tested against zero across patients with a two-sided
Wilcoxon signed-rank test (>= 6 patients; with 5 the smallest two-sided p is 0.0625), with the matched-pairs rank-biserial
correlation and the fraction of patients above zero as effect sizes, and
Benjamini-Hochberg q-values over all reference x measure tests of a dataset
(`noise_contrasts.parquet`; per-patient values in `noise_contrasts_patient.parquet`). This
contrast is an scTCNA addition. Because the pooled non-malignant compartment mixes
lineages, its raw CV includes between-lineage differences; the single-lineage references
remove most of that composition effect.

## Batch effects

Every comparison is made within one study and within patients, so study- and
patient-level batch effects cancel in the paired contrast. Cross-cancer views use
within-dataset quantities only (ranks, z-scores of per-dataset values); there is no
cross-dataset integration or batch correction, and differences between cancer types are
therefore confounded with differences between studies.

## DEPTH and DEPTH2

Both follow the R packages used by TCNA (WangX-Lab/DEPTH, WangX-Lab/DEPTH2) and are
applied to log2-transformed expression, as in TCNA. "Sample" is a cell or a pseudo-bulk
profile. Genes enter the scores when detected in >= 1% of the dataset's cells.

```
DEPTH_s  = SD_g[ (x_gs - mu_g^ref)^2 ]        mu_g^ref = mean of gene g over reference (normal) samples,
                                               or over the tumour samples when there are no normals
DEPTH2_s = SD_g[ |z_gs| ]                     z_gs = (x_gs - mu_g) / sigma_g across samples,
                                               sigma_g = 0 (constant gene) -> z = 0, gene still counted
```

- **Per cell:** `x` is log2(CP10K + 1). DEPTH2 z-scores every gene across the cells of
  the cell's own patient x compartment group, so the score says how atypical a cell is
  within its own tumour (or within the patient's non-malignant cells), the single-cell
  analogue of DEPTH2's cohort of tumours. DEPTH is computed for malignant cells only,
  against the mean over all non-malignant cells of the dataset (DEPTH's normal
  reference); non-malignant cells have no DEPTH. Scoring both compartments against one
  pooled reference or one pooled z-score would turn differences in mean expression
  between compartments into apparent heterogeneity differences. Scores are in
  `cell_scores` (`depth`, `depth2`) and summarised per patient x compartment as medians
  (`depth_cell_median`, `depth2_cell_median`).
- **Per pseudo-bulk:** counts of each patient x compartment group (>= 20 cells) are
  summed over the DEPTH genes, scaled to counts per million of those genes and
  log2(x + 1)-transformed, which is the closest
  analogue of a bulk sample. DEPTH2 is computed across all pseudo-bulks of the dataset
  (`depth2_pseudobulk`). DEPTH is computed for malignant pseudo-bulks against the mean of
  the non-malignant pseudo-bulks (`depth_pseudobulk`; null for non-malignant profiles, by
  DEPTH's definition). If a dataset had no non-malignant pseudo-bulk, all profiles would be
  scored against their own mean.

## Outputs

Per dataset, `python -m sctcna run` writes to `data/processed/<dataset>/`:
`gene_metrics_patient.parquet`, `gene_metrics_pooled.parquet`,
`gene_metrics_celltype.parquet`, `gene_summary.parquet`, `cell_scores.parquet`,
`sample_scores.parquet` and `summary.json` (QC counts, downsampling target and
concordance, test tallies, DEPTH tests, runtime, parameters).

`python -m sctcna aggregate` combines them into the tables the app serves
(`data/app/`). `noise_contrasts.parquet` and `noise_contrasts_patient.parquet` are also written per dataset. A cancer type maps to exactly one dataset in this release.

| Table | One row per | Columns |
|---|---|---|
| `datasets.parquet` | dataset | `dataset, cancer_type, cancer_name, tcna_site, source, study, citation, technology, site, data_url, meta_url, notes, n_patients, n_patients_paired, n_cells_malignant, n_cells_non_malignant, n_genes, n_genes_expressed, doublets_removed, downsampling_target_depth, spearman_cv_vs_downsampled, n_sig_q05, n_sig_noisier_malignant, n_sig_noisier_non_malignant, median_log2fc_cv, median_delta_cv2_resid, depth2_pseudobulk_median_diff, depth2_pseudobulk_p, depth2_cell_median_diff, depth2_cell_p, log2fc_cv_ds_vs_non_malignant, excess_cv2_ds_vs_non_malignant, excess_cv2_ds_p, download_date, sha256, sctcna_version, runtime_seconds` |
| `gene_summary.parquet` | dataset x gene | `dataset, gene, mean_malignant, mean_non_malignant, detection_rate_malignant, detection_rate_non_malignant, sd_malignant, sd_non_malignant, cv_malignant, cv_non_malignant, cv2_malignant, cv2_non_malignant, mad_malignant, mad_non_malignant, fano_malignant, fano_non_malignant, cv2_resid_malignant, cv2_resid_non_malignant, n_patients_malignant, n_patients_non_malignant, n_pairs, log2fc_cv_median, log2fc_cv_mean, wilcoxon_stat, p_value, rank_biserial, q_value, log2fc_mad_median, p_value_mad, q_value_mad, log2fc_fano_median, p_value_fano, q_value_fano, delta_cv2_resid, p_value_cv2_resid, q_value_cv2_resid, cancer_type, log2fc_cv_of_medians, expressed, cv_pooled_malignant, cv_pooled_non_malignant`. Metric columns are medians across patients of the per-patient values. |
| `cell_scores.parquet` | cell | `dataset, cancer_type, patient, cell_type, compartment, n_counts, n_genes, pct_mito, depth, depth2` |
| `sample_scores.parquet` | patient x compartment | `patient, compartment, n_cells, depth2_pseudobulk, depth_pseudobulk, depth_cell_median, depth2_cell_median, dataset, cancer_type` (`depth_pseudobulk` and `depth_cell_median` are null for non-malignant profiles) |
| `celltype_summary.parquet` | dataset x cell type x gene | `dataset, cancer_type, cell_type, gene, n_patients, mean, detection_rate, cv, mad, fano, cv2_resid` (medians across patients; `Malignant` is one of the cell types; expressed genes only) |
| `patient/<dataset>.parquet` | patient x compartment x gene | `dataset, cancer_type, patient, compartment, gene, n_cells, mean, detection_rate, cv, mad, fano, cv2_resid, cv_ds, cv2_resid_ds` (expressed genes only; `*_ds` after downsampling) |
| `noise_contrasts.parquet` | dataset x reference group | `dataset, cancer_type, reference, n_patients, median_log2fc_cv, p_log2fc_cv, median_log2fc_cv_ds, p_log2fc_cv_ds, median_excess_cv2_ds, p_excess_cv2_ds, rank_biserial_excess, frac_patients_excess_positive, q_log2fc_cv, q_log2fc_cv_ds, q_excess_cv2_ds` (q: BH over every reference x measure test of the dataset) |
| `noise_contrasts_patient.parquet` | patient x reference group | `dataset, cancer_type, patient, reference, n_cells_malignant, n_cells_reference, n_genes, log2fc_cv, log2fc_cv_ds, excess_cv2_ds` |
| `heterogeneity_tests.parquet` | cancer type x level | `cancer_type, score, level, n_pairs, median_diff, p_value, rank_biserial` (paired DEPTH2 tests per cell and per pseudo-bulk, precomputed so the app needs no SciPy) |
| `build.json` | build | `sctcna_version, datasets, built` |

In the app tables, repeated strings are stored as categoricals and display-only values
(per-cell scores, per-patient and per-cell-type metrics) as float32 to keep the served
tables within the 512 MB of the free hosting tier; `gene_summary` keeps float64.

## Assumptions: manuscript vs scTCNA

| Item | Origin |
|---|---|
| CV = SD / mean, CV^2, SD, MAD = mean absolute deviation about the mean | TCNA manuscript / code |
| Differential noise = log2(CV_tumour / CV_normal) (Eq. 1) | TCNA manuscript |
| log2(x + 1) before DEPTH / DEPTH2 | TCNA manuscript (after Song & Wang 2022) |
| DEPTH and DEPTH2 definitions | DEPTH / DEPTH2 R packages used by TCNA |
| No global low-expression filter; mean shown next to every noise value | TCNA manuscript |
| CV stored as a ratio (TCNA web code reports CV x 100) | scTCNA storage convention |
| Tumour vs normal replaced by malignant vs non-malignant cells of the same patient | scTCNA assumption |
| Non-malignant = all curated non-malignant lineages pooled; unassigned cells excluded | scTCNA assumption |
| Malignant calls taken from the atlas, not re-derived | scTCNA assumption |
| Cap of 3,000 cells per patient x compartment (random, seed 0, before QC) | scTCNA assumption |
| QC thresholds: >= 200 genes, >= 500 UMIs, <= 20% mito, genes in >= 10 cells, Scrublet per sample, groups >= 20 cells | scTCNA assumption |
| CP10K for gene noise (TCNA: TPM / FPKM / FPKM-UQ); log2(CP10K + 1) and log2(CPM + 1) for DEPTH metrics | scTCNA assumption |
| Fano factor and detection rate | scTCNA addition |
| Mean-corrected CV^2 (LOWESS residual, frac 0.3) and its compartment difference | scTCNA addition |
| Downsampling to a common UMI depth (10th percentile, floor 1,000) | scTCNA addition |
| Eq. 1 evaluated per patient and summarised as the median across patients; same ratio for MAD and Fano | scTCNA assumption |
| Patient as unit of replication; paired Wilcoxon (>= 6 patients); rank-biserial; BH per cancer type | scTCNA assumption |
| Per-cell DEPTH for malignant cells only (reference: all non-malignant cells); per-cell DEPTH2 z-scored within patient x compartment | scTCNA assumption |
| Expressed = detection rate >= 10% in both compartments (median over patients) for tests / top lists | scTCNA assumption |
| Pooled-cell CV stored as a pseudoreplication diagnostic | scTCNA addition |
| Malignant vs reference contrasts (median over genes per patient; depth-matched excess over the reference mean-CV^2 trend; one-sample Wilcoxon across patients) | scTCNA addition |
| Lineage-level (patient x cell type, >= 50 cells) metrics | scTCNA addition |
| Cells as DEPTH / DEPTH2 samples; z-scores across all cells of a dataset; DEPTH reference = non-malignant cells | scTCNA assumption |
| Pseudo-bulk = summed counts per patient x compartment, log2(CPM + 1); DEPTH2 across all pseudo-bulks; DEPTH vs non-malignant pseudo-bulk mean | scTCNA assumption |
| Genes for DEPTH / DEPTH2 detected in >= 1% of cells | scTCNA assumption |

## Limitations

- One study per cancer type, so cancer-type and study effects cannot be separated;
  cross-cancer comparisons are confounded by study design, tissue handling and
  sequencing depth.
- The non-malignant compartment is a mixture of lineages (immune, stromal, normal
  epithelium); composition differences between patients contribute to its noise. The
  lineage-level table mitigates this by allowing malignant cells to be compared with one
  matched normal lineage.
- Malignant calls come from the atlases (CNA-based) and are not re-derived here.
- The melanoma study (Biermann et al. 2022) consists of metastases only. The liver dataset
  (Ma et al. 2019) mixes hepatocellular carcinoma and intrahepatic cholangiocarcinoma, as
  TCNA's Liver site does, and is small (6 paired patients, 2,596 cells after QC). Sun et
  al. 2021 was replaced because 3CA provides it only as normalised TPM values.
- The per-patient cap and the uniform QC differ from the original studies' processing,
  so cell numbers do not match the publications.
- TCNA's Pathway dashboard (STRING enrichment) is not replicated, and there is no custom
  data upload.
- Only 10x UMI data are used; the TISCH2 loader is provided but its entries are disabled.
- Statistical power of the paired tests is limited by patient numbers (6-27
  paired patients per cancer type).
- Dropout and zero inflation inflate CV and MAD for lowly expressed genes; the
  mean-corrected and downsampled versions are provided for this reason and should be read
  alongside the raw metrics.

## References

- Shah R, Faisal A, bin Faiz A, Bashir I, Jawaid A, Bashir K, Shoaib M, Chaudhary SU.
  The Cancer Noise Atlas (TCNA): A Webserver of Multiscale Transcriptional Noise in
  Cancer. *J Mol Biol* 2026. <https://tcna.lums.edu.pk/>; code
  <https://github.com/BIRL/TCNA>.
- Gavish A, Tyler M, Greenwald AC, et al. Hallmarks of transcriptional intratumour
  heterogeneity across a thousand tumours. *Nature* 618, 598-606 (2023).
  <https://www.weizmann.ac.il/sites/3CA/>
- Han Y, Wang Y, Dong X, et al. TISCH2: expanded datasets and new tools for single-cell
  transcriptome analyses of the tumor microenvironment. *Nucleic Acids Res* 51,
  D1425-D1431 (2023). <https://tisch.compbio.cn/>
- Li M, Zhang Z, Li L, Wang X. An algorithm to quantify intratumor heterogeneity based
  on alterations of gene expression profiles. *Commun Biol* 3, 505 (2020).
  <https://github.com/WangX-Lab/DEPTH>
- Song D, Wang X. DEPTH2: an mRNA-based algorithm to evaluate intratumor heterogeneity
  without reference to normal controls. *J Transl Med* 20, 150 (2022).
  <https://github.com/WangX-Lab/DEPTH2>
- Brennecke P, Anders S, Kim JK, et al. Accounting for technical noise in single-cell
  RNA-seq experiments. *Nat Methods* 10, 1093-1095 (2013).
- Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets
  in single-cell transcriptomic data. *Cell Syst* 8, 281-291 (2019).
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data
  analysis. *Genome Biol* 19, 15 (2018).
- Cleveland WS. Robust locally weighted regression and smoothing scatterplots.
  *J Am Stat Assoc* 74, 829-836 (1979).
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and
  powerful approach to multiple testing. *J R Stat Soc B* 57, 289-300 (1995).
- Kerby DS. The simple difference formula: an approach to teaching nonparametric
  correlation. *Compr Psychol* 3, 11.IT.3.1 (2014).
- Study citations for the ten datasets are in `configs/datasets.yaml`.
