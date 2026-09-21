# ============================================================
# 04_PCA_clustering_UMAP.R
#
# Unintegrated PCA / clustering / UMAP
#
# This step is a reference/control analysis before integration.
#
# Pipeline:
#   03_normalized.rds
#       ↓
#   ScaleData on 2,000 HVGs
#       ↓
#   PCA
#       ↓
#   Neighbors
#       ↓
#   Clustering
#       ↓
#   UMAP
#
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
library(ggplot2)

project_dir <- get_project_dir()

input <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
)

if (!file.exists(input)) {
  stop(
    "Input object not found: ",
    input
  )
}

# ============================================================
# LOAD OBJECT
# ============================================================

message("============================================================")
message("Loading normalized object")
message("============================================================")

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop(
    "Input is not a Seurat object."
  )
}

DefaultAssay(object) <- "RNA"

message(
  "Cells: ",
  ncol(object)
)

message(
  "Features: ",
  nrow(object)
)

# ============================================================
# CHECK HVGs
# ============================================================

hvg <- VariableFeatures(object)

message("============================================================")
message("Variable features")
message("============================================================")

message(
  "Number of HVGs: ",
  length(hvg)
)

if (length(hvg) < 100) {
  stop(
    "Too few variable features for PCA."
  )
}

# ============================================================
# SCALE HVGs
# ============================================================

message("============================================================")
message("Scaling HVGs for unintegrated PCA")
message("============================================================")

object <- ScaleData(
  object,
  features = hvg,
  verbose = FALSE
)

gc()

# ============================================================
# VALIDATE SCALE.DATA
# ============================================================

scale_layers <- Layers(
  object[["RNA"]],
  search = "^scale\\.data$"
)

message("============================================================")
message("Scale layer")
message("============================================================")

print(scale_layers)

if (length(scale_layers) == 0) {
  stop(
    "scale.data layer was not created by ScaleData()."
  )
}

# ============================================================
# PCA
# ============================================================

message("============================================================")
message("Running unintegrated PCA")
message("============================================================")

npcs <- 30

object <- RunPCA(
  object,
  features = hvg,
  npcs = npcs,
  verbose = FALSE
)

if (!"pca" %in% Reductions(object)) {
  stop(
    "PCA reduction was not created."
  )
}

message(
  "PCA dimensions: ",
  ncol(
    Embeddings(
      object,
      reduction = "pca"
    )
  )
)

# ============================================================
# ELBOW PLOT
# ============================================================

message("============================================================")
message("Generating elbow plot")
message("============================================================")

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
    ndims = npcs
  )
)

dev.off()

# ============================================================
# DIMENSIONS
# ============================================================

dims_use <- 1:npcs

# ============================================================
# NEIGHBORS
# ============================================================

message("============================================================")
message("Finding unintegrated neighbors")
message("============================================================")

object <- FindNeighbors(
  object,
  reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)

# ============================================================
# CLUSTERING
# ============================================================

message("============================================================")
message("Finding unintegrated clusters")
message("============================================================")

object <- FindClusters(
  object,
  resolution = 0.4,
  verbose = FALSE
)

# ============================================================
# UMAP
# ============================================================

message("============================================================")
message("Running unintegrated UMAP")
message("============================================================")

object <- RunUMAP(
  object,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap.unintegrated",
  reduction.key = "unintegratedUMAP_",
  verbose = FALSE
)

if (!"umap.unintegrated" %in% Reductions(object)) {
  stop(
    "Unintegrated UMAP was not created."
  )
}

# ============================================================
# UMAP BY CLUSTER
# ============================================================

message("============================================================")
message("Generating cluster UMAP")
message("============================================================")

p_cluster <- DimPlot(
  object,
  reduction = "umap.unintegrated",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Unintegrated PCA / UMAP - clusters"
  )

ggsave(
  file.path(
    project_dir,
    "results",
    "figures",
    "integration",
    "UMAP_unintegrated_clusters.png"
  ),
  p_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# UMAP BY SAMPLE
# ============================================================

message("============================================================")
message("Generating sample UMAP")
message("============================================================")

p_sample <- DimPlot(
  object,
  reduction = "umap.unintegrated",
  group.by = "sample"
) +
  ggtitle(
    "Unintegrated PCA / UMAP - samples"
  )

ggsave(
  file.path(
    project_dir,
    "results",
    "figures",
    "integration",
    "UMAP_unintegrated_samples.png"
  ),
  p_sample,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# UMAP BY TREATMENT
# ============================================================

if (
  "treatment" %in%
  colnames(object@meta.data)
) {

  message("============================================================")
  message("Generating treatment UMAP")
  message("============================================================")

  p_treatment <- DimPlot(
    object,
    reduction = "umap.unintegrated",
    group.by = "treatment"
  ) +
    ggtitle(
      "Unintegrated PCA / UMAP - treatment"
    )

  ggsave(
    file.path(
      project_dir,
      "results",
      "figures",
      "integration",
      "UMAP_unintegrated_treatment.png"
    ),
    p_treatment,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# CELL COUNTS
# ============================================================

cluster_counts <- as.data.frame(
  table(
    object$seurat_clusters
  )
)

colnames(cluster_counts) <- c(
  "cluster",
  "cells"
)

write.csv(
  cluster_counts,
  file.path(
    project_dir,
    "results",
    "qc",
    "unintegrated_cluster_counts.csv"
  ),
  row.names = FALSE
)

# ============================================================
# SAVE
# ============================================================

output <- file.path(
  project_dir,
  "results",
  "objects",
  "04_unintegrated.rds"
)

saveRDS(
  object,
  output
)

message("============================================================")
message("Unintegrated PCA / clustering / UMAP completed")
message("============================================================")

message(
  "Cells: ",
  ncol(object)
)

message(
  "Clusters: ",
  length(
    unique(
      object$seurat_clusters
    )
  )
)

message(
  "Output: ",
  output
)

message("============================================================")
