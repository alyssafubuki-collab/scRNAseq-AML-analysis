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
