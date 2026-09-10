# Single-Cell RNA-seq Analysis in AML

A reproducible single-cell RNA-seq workflow in **R / Seurat v5** using a public AML dataset.

## Dataset

This portfolio uses **GSE145410**, a public GEO dataset containing 8 10X Genomics scRNA-seq samples from an AML bone-marrow experiment:

- DMSO A / B
- INCB059872 A / B
- AZA A / B
- INCB059872 + AZA A / B

The dataset is publicly available from NCBI GEO.

The workflow downloads the public matrices automatically; no clinical or private CRCM data are included.

## Analysis workflow

```text
GSE145410
   |
   v
Download 10X matrices
   |
   v
Create Seurat objects
   |
   v
QC metrics
   |
   v
Adaptive QC with scuttle::isOutlier()
   |
   v
Normalization + 2,000 HVGs
   |
   v
PCA
   |
   v
Seurat v5 CCA integration
   |
   v
Clustering + UMAP
   |
   +--------------------+
   |                    |
   v                    v
SingleR              scPred
   |                    |
   +---------+----------+
             |
             v
          scAnnoX
             |
             v
   Annotation comparison
             |
             v
 Differential expression
             |
             v
       Final figures
```

## Annotation methods

The annotation comparison is deliberately limited to the three approaches used in the research workflow represented by this portfolio:

- **SingleR**
- **scPred**
- **scAnnoX**

ACTINN is not part of this repository.

### Reference datasets

SingleR uses the hematopoietic **Novershtern** reference available through `celldex`.

scPred uses its public PBMC reference (`scPred::pbmc_1`) to train supervised classifiers.

scAnnoX uses the same public scPred PBMC reference when the local Seurat-5-compatible scAnnoX package is available.

Because these references do not contain every AML leukemic state, cells that do not match the reference populations may remain unassigned. This is preferable to forcing an inappropriate cell identity.

## Important note about scAnnoX

The original scAnnoX package was developed before the current Seurat v5 layer system. The original package is therefore **not silently substituted** here.

The repository expects the Seurat-v5-compatible version of scAnnoX used in the original research workflow to be placed in:

```text
external/scAnnoX/
```

The script `08_annotation_scAnnoX.R` checks for this directory and stops with a clear message if it is absent.

This keeps the portfolio reproducible and makes the Seurat-v5 compatibility modification explicit.

## Project structure

```text
scRNAseq-AML-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── data/
│   └── README.md
├── external/
│   └── scAnnoX/
├── scripts/
│   ├── 00_install_packages.R
│   ├── 01_download_GSE145410.R
│   ├── 02_import_and_QC.R
│   ├── 03_normalization_HVG.R
│   ├── 04_PCA_clustering_UMAP.R
│   ├── 05_Seurat5_integration.R
│   ├── 06_annotation_SingleR.R
│   ├── 07_annotation_scPred.R
│   ├── 08_annotation_scAnnoX.R
│   ├── 09_compare_annotations.R
│   ├── 10_differential_expression.R
│   └── 11_final_figures.R
├── R/
│   ├── functions.R
│   └── plotting_functions.R
├── figures/
│   ├── qc/
│   ├── integration/
│   ├── annotation/
│   └── differential_expression/
├── results/
│   ├── qc/
│   ├── annotation/
│   └── differential_expression/
├── docs/
│   └── workflow.md
└── environment/
```

## Reproducibility

Run the scripts in this order:

```text
00_install_packages.R
01_download_GSE145410.R
02_import_and_QC.R
03_normalization_HVG.R
04_PCA_clustering_UMAP.R
05_Seurat5_integration.R
06_annotation_SingleR.R
07_annotation_scPred.R
08_annotation_scAnnoX.R
09_compare_annotations.R
10_differential_expression.R
11_final_figures.R
```

## Skills demonstrated

- R
- Seurat v5
- Bioconductor
- scuttle
- SingleR
- scPred
- scAnnoX
- PCA
- CCA integration
- clustering
- UMAP
- differential expression
- reproducible workflow organisation
- Git / GitHub

## Background

This portfolio repository was designed from the computational workflow developed during a research experience in single-cell RNA-seq analysis of AML treatment response and resistance.

No patient-level or confidential research data are included.
