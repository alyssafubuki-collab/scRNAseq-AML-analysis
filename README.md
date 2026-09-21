# Single-Cell RNA-seq Analysis in AML

A reproducible single-cell RNA-seq workflow in **R / Seurat v5** using a public AML dataset.

## Dataset

This portfolio uses **GSE145410**, a public GEO dataset containing 8 10X Genomics scRNA-seq samples from an AML bone-marrow experiment:

- DMSO A / B
- INCB059872 A / B
- AZA A / B
- INCB059872 + AZA A / B

The dataset is publicly available from NCBI GEO.

The workflow automatically downloads the public matrices; it includes no clinical or private data.

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
   Annotation comparison
             |
             v
 Differential expression
             |
             v
       Final figures
```

## Annotation methods

The annotation comparison is deliberately limited to the two approaches used in the research workflow represented by this portfolio:

- **SingleR**
- **scPred**

### Reference datasets

SingleR uses the hematopoietic **Novershtern** reference available through `celldex`.

scPred uses its public PBMC reference (`scPred::pbmc_1`) to train supervised classifiers.

Because these references do not contain every AML leukemic state, cells that do not match the reference populations may remain unassigned. This is preferable to forcing an inappropriate cell identity.

## Project structure

```text
scRNAseq-AML-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── data/
│   └── README.md
├── scripts/
│   ├── 00_install_packages.R
│   ├── 01_download_GSE145410.R
│   ├── 02_import_and_QC.R
│   ├── 03_normalization_HVG.R
│   ├── 04_PCA_clustering_UMAP.R
│   ├── 05_Seurat5_integration.R
│   ├── 06_annotation_SingleR.R
│   ├── 07_annotation_scPred.R
│   ├── 08b_consensus_annotations.R
│   ├── 09b_differential_expression.R
│   └── 10b_final_figures.R
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
08b_consensus_annotations.R
09b_differential_expression.R
10b_final_figures.R
```

## Skills demonstrated

- R
- Seurat v5
- Bioconductor
- scuttle
- SingleR
- scPred
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
