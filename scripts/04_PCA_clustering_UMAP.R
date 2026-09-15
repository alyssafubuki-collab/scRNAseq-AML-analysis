# ============================================================
# 04_PCA_clustering_UMAP.R
# Unintegrated PCA, neighbors, clustering and UMAP
#
# This step is intentionally kept as a pre-integration diagnostic.
# Only the 2,000 HVGs are scaled.
# ============================================================

source("R/functions.R")

library(Seurat)
library(ggplot2)

project_dir <- get_project_dir()

object <- readRDS(
  file.path(
    project_dir,
    "results",
    "objects",
    "03_normalized.rds"
  )
)

DefaultAssay(object) <- "RNA"

if (length(VariableFeatures(object)) == 0) {
  stop("No variable features found.")
}

# ============================================================
# SCALE ONLY HVGs
# ============================================================

message("Scaling 2,000 HVGs for pre-integration PCA...")

object <- ScaleData(
  object,
  features = VariableFeatures(object),
  verbose = FALSE
)

# ============================================================
# PCA
# ============================================================

object <- RunPCA(
  object,
  features = VariableFeatures(object),
  npcs = 20,
  verbose = FALSE
)

pdf(
  file.path(
    project_dir,
    "figures",
    "integration",
    "PCA_elbow_plot.pdf"
  ),
  width = 7,
  height = 5
)

print(
  ElbowPlot(
    object,
    ndims = 20
  )
)

dev.off()

# ============================================================
# NEIGHBORS + CLUSTERS
# ============================================================

object <- FindNeighbors(
  object,
  dims = 1:20,
  verbose = FALSE
)

object <- FindClusters(
  object,
  resolution = 0.4,
  verbose = FALSE
)

# ============================================================
# UMAP
# ============================================================

object <- RunUMAP(
  object,
  dims = 1:20,
  reduction.name = "umap.unintegrated",
  reduction.key = "unintegratedUMAP_",
  verbose = FALSE
)

# ============================================================
# PLOTS
# ============================================================

p1 <- DimPlot(
  object,
  reduction = "umap.unintegrated",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "integration",
    "UMAP_unintegrated_clusters.png"
  ),
  p1,
  width = 8,
  height = 6,
  dpi = 300
)

p2 <- DimPlot(
  object,
  reduction = "umap.unintegrated",
  group.by = "sample"
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "integration",
    "UMAP_unintegrated_samples.png"
  ),
  p2,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# SAVE
# ============================================================

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "04_unintegrated.rds"
  )
)

message("PCA, clustering and unintegrated UMAP completed.")
