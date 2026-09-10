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
