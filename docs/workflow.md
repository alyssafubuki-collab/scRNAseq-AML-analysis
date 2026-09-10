# Analysis workflow

## Overview

The analysis follows a standard single-cell RNA-seq workflow.

```text
Raw counts
    |
    v
Quality control
    |
    v
Cell filtering
    |
    v
Normalization
    |
    v
Highly variable genes
    |
    v
Scaling
    |
    v
PCA
    |
    v
Integration
    |
    v
Nearest-neighbor graph
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

Quality control
```

Cells are filtered according to:

Number of detected genes
Total RNA counts
Mitochondrial RNA percentage

Example thresholds:
```text
nFeature_RNA > 300
nFeature_RNA < 2500
percent.mt < 10
```
These thresholds are dataset-dependent.

Normalization

Expression counts are normalized using Seurat's
LogNormalize method.

Highly variable genes are identified using the variance-stabilizing
transformation (vst).

Dimensionality reduction

PCA is performed using the selected highly variable genes.

The first 30 principal components are used for downstream
neighborhood analysis and UMAP.

Integration

For multi-sample analysis, Seurat v5 integration is performed using
CCA integration.

The integrated representation is then used for:

Neighbor detection
Clustering
UMAP
Annotation

Cell identities can be assigned using:

Marker gene expression
SingleR
scPred
ACTINN

Using several methods allows comparison of annotation consistency.

Differential expression

Differentially expressed genes are identified between cell
populations or experimental groups.

The analysis can be adapted to:

Cell type
Treatment
Time point
Biological response
Patient groups
