# Methods

## Software

The workflow was developed in R using:

- Seurat
- SeuratObject
- Bioconductor
- ggplot2
- dplyr
- patchwork

## Quality control

Cells with fewer than 300 detected genes or more than 2500 detected
genes can be excluded.

Cells with more than 10% mitochondrial transcripts can also be
excluded.

These thresholds are examples and should be adapted according to
dataset-specific QC distributions.

## Normalization

Gene expression counts are normalized using the Seurat
LogNormalize approach with a scale factor of 10,000.

## Highly variable genes

The top 2,000 highly variable genes are selected using the
variance-stabilizing transformation method.

## PCA

Principal component analysis is performed on the selected variable
genes.

The first 30 principal components are used for downstream analysis.

## Clustering

A shared nearest-neighbor graph is constructed from the selected
principal components.

Clusters are identified using Seurat's graph-based clustering
algorithm.

## UMAP

UMAP is performed using the first 30 principal components.

## Integration

For multi-sample datasets, Seurat v5 CCA integration can be used to
reduce technical differences between datasets.

## Annotation

Cell identities can be assigned based on:

1. Canonical marker genes
2. Reference-based annotation
3. Comparison between annotation methods

## Differential expression

Differential expression is performed using Seurat's
`FindMarkers()` and `FindAllMarkers()` functions.

Genes can be prioritized according to:

- Adjusted p-value
- Log fold change
- Percentage of cells expressing the gene
