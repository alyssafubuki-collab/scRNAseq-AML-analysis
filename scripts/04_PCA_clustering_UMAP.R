# ============================================================
# 04_PCA_clustering_UMAP.R
# PCA, neighbors, clustering and UMAP before integration
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

object <- RunPCA(
  object,
  features = VariableFeatures(object),
  npcs = 30,
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
    ndims = 30
  )
)

dev.off()

object <- FindNeighbors(
  object,
  dims = 1:30,
  verbose = FALSE
)

object <- FindClusters(
  object,
  resolution = 0.4,
  verbose = FALSE
)

object <- RunUMAP(
  object,
  dims = 1:30,
  verbose = FALSE
)

p1 <- DimPlot(
  object,
  reduction = "umap",
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
  reduction = "umap",
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

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "04_unintegrated.rds"
  )
)

message("PCA, clustering and UMAP completed.")
