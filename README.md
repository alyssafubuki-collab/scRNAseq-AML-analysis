# Single-Cell RNA-seq Analysis in AML

A reproducible single-cell RNA-seq analysis workflow developed in R
using Seurat and Bioconductor.

## Overview

This repository presents a complete workflow for the analysis of
single-cell RNA sequencing data.

The workflow covers the main stages of a typical scRNA-seq analysis:

- Quality control
- Normalization
- Highly variable gene selection
- Principal component analysis
- Dataset integration
- Clustering
- UMAP visualization
- Cell-type annotation
- Differential gene expression
- Biological interpretation

The project is designed as a bioinformatics portfolio demonstrating
practical experience with R, Seurat, Bioconductor and single-cell
transcriptomics.

> This repository is a public portfolio version of analytical
> workflows developed during my research experience.
>
> No patient-level or confidential research data are included.

---

## Biological context

The original research experience involved the analysis of
single-cell RNA-seq data in the context of treatment response
and resistance in acute myeloid leukemia (AML).

The public version of this repository focuses on the computational
workflow rather than confidential clinical data.

---

## Workflow

```text
Single-cell RNA-seq data
          |
          v
   Quality control
          |
          v
    Normalization
          |
          v
 Highly variable genes
          |
          v
        PCA
          |
          v
 Dataset integration
          |
          v
      Clustering
          |
          v
        UMAP
          |
          v
   Cell annotation
          |
          v
 Differential expression
          |
          v
 Biological interpretation
```
---

Main analyses
---
## 1. Quality control

The workflow evaluates:

- Number of detected genes
- Total RNA counts
- Mitochondrial gene percentage
- Distribution of QC metrics
- Filtering of low-quality cells

Example filtering strategy:
```text
nFeature_RNA > 300
nFeature_RNA < 2500
percent.mt < 10
```
These thresholds are examples and should be adapted to the
characteristics of each dataset.

## 2. Normalization

Gene expression counts are normalized using the Seurat workflow.

The analysis includes:

- Log normalization
- Highly variable gene identification
- Scaling
- Principal component analysis

## 3. Dimensionality reduction

Principal component analysis (PCA) is used to reduce the
dimensionality of the expression matrix.

The first principal components are then used for:

- Nearest-neighbor graph construction
- Clustering
- UMAP visualization
  
## 4. Dataset integration

Multiple samples can be integrated to reduce technical
variation while preserving biological differences.

The workflow supports Seurat v5 integration using
CCAIntegration.

## 5. Cell-type annotation

Several annotation approaches can be compared:

- SingleR
- scPred
- scAnnoX

The purpose is to evaluate the consistency of predicted
cell identities across different approaches.

## 6. Differential expression

Differentially expressed genes can be identified between:

- Cell populations
- Experimental conditions
- Treatment time points
- Biological response groups

The workflow uses Seurat differential expression functions.

---
## Technologies
---
## Programming
- R
- Bash
- Linux
  
## Bioinformatics
- Seurat
- SeuratObject
- Bioconductor
- SingleR
- scPred
- scAnnoX
  
## Data analysis
- PCA
- Clustering
- UMAP
- Differential expression
- Dimensionality reduction

## Visualization
- ggplot2
- Seurat visualization tools

## Reproducibility
- Git
- GitHub
- Session information
- Script-based workflow

---
Repository structure
---
```text
scRNAseq-AML-analysis/
│
├── scripts/
│   ├── 00_setup.R
│   ├── 01_quality_control.R
│   ├── 02_normalization_hvg.R
│   ├── 03_pca_clustering_umap.R
│   ├── 04_integration.R
│   ├── 05_cell_annotation.R
│   ├── 06_differential_expression.R
│   └── 07_visualization.R
│
├── R/
│   ├── functions.R
│   └── plotting_functions.R
│
├── data/
├── figures/
├── results/
├── docs/
└── environment/
```
---
Reproducibility
---
The analysis is organized as a sequence of numbered scripts.

The recommended execution order is:
```text
00_setup.R
01_quality_control.R
02_normalization_hvg.R
03_pca_clustering_umap.R
04_integration.R
05_cell_annotation.R
06_differential_expression.R
07_visualization.R
```

Package versions and session information can be stored in:
```text
environment/packages.R
environment/sessionInfo.txt
```
---
Data availability
---
No clinical or patient-level data are included in this repository.

To reproduce the workflow, a public scRNA-seq dataset can be placed
in the data/ directory.

The repository is designed so that the analysis scripts can be
adapted to different scRNA-seq datasets.

---
My contribution
---
I developed and implemented the bioinformatics workflows presented
in this repository.

This includes:

- Data preprocessing
- Quality control
- Normalization
- Feature selection
- Dimensionality reduction
- Dataset integration
- Clustering
- Cell-type annotation
- Differential expression analysis
- Data visualization

The workflow was developed from practical experience analyzing
single-cell RNA-seq data during a research internship.

---
Background
---

This work was developed during my research experience at the
Centre de Recherche en Cancérologie de Marseille (CRCM),
Aix-Marseille Université.

The original research project involved single-cell analysis of
treatment response and resistance in acute myeloid leukemia.

Only non-confidential computational material is presented here.

---
Author
---
Alyssa Chellal

Master Biologie-Santé parcours Biomarkers & Artificial Intelligence

Aix-Marseille Université

Bioinformatics | Single-cell RNA-seq | Multi-omics | Data analysis

---
License
---
This project is released under the MIT License.

