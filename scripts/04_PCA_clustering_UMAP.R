# ============================================================
# 04_PCA_clustering_UMAP.R
# Unintegrated PCA and UMAP
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
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

# ============================================================
# PCA
# ============================================================

message("Running unintegrated PCA...")

object <- ScaleData(
  object,
  features = VariableFeatures(object),
  verbose = FALSE
)

object <- RunPCA(
  object,
  features = VariableFeatures(object),
  npcs = 30,
  verbose = FALSE
)

gc()

# ============================================================
# ELBOW
# ============================================================

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
    ndims = 30
  )
)

dev.off()

# ============================================================
# NEIGHBORS
# ============================================================

object <- FindNeighbors(
  object,
  dims = 1:30,
  verbose = FALSE
)

# ============================================================
# CLUSTERS
# ============================================================

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
  dims = 1:30,
  reduction = "pca",
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

rm(p1, p2)
gc()

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

message("Unintegrated PCA/UMAP completed.")
