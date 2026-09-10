# Workflow

## 1. Public data

GSE145410 is downloaded from GEO as the supplementary 10X matrix archive.

## 2. Import

Each GEO sample is converted into a Seurat object and receives explicit metadata:

- sample
- treatment
- replicate

## 3. Quality control

QC is adaptive rather than based on fixed thresholds.

For each sample, `scuttle::isOutlier()` is applied to:

- `nCount_RNA`: lower outliers, log-transformed
- `nFeature_RNA`: lower outliers, log-transformed
- `percent.mt`: upper outliers

The sample identifier is passed as the `batch` argument so that QC thresholds are calculated independently for each sample.

This follows the MAD-based outlier strategy used in the original workflow.

## 4. Normalization

Seurat LogNormalize is used with a scale factor of 10,000.

The top 2,000 variable features are selected with the `vst` method.

## 5. PCA and clustering

PCA is performed using the variable features.

The first 30 PCs are used for nearest-neighbor graph construction, clustering and UMAP.

## 6. Integration

Seurat v5 CCA integration is performed across the eight samples.

The integrated representation is used for clustering and UMAP.

## 7. Annotation

Three annotation approaches are compared:

1. SingleR
2. scPred
3. scAnnoX

No ACTINN workflow is included.

## 8. Differential expression

Differential expression is performed after joining the RNA layers.

The analysis includes:

- cluster markers
- treatment comparisons

The treatment comparison is exploratory because GSE145410 contains ex-vivo treatment replicates from a single AML patient rather than an independent patient cohort.
