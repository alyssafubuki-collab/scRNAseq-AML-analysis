# ============================================================
# 04_PCA_clustering_UMAP.R
#
# Unintegrated PCA / clustering / UMAP
# Used as a pre-integration reference.
# ============================================================

source("R/functions.R")

library(Seurat)
library(ggplot2)

project_dir <- get_project_dir()

input <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
)

if (!file.exists(input)) {
  stop("Input object not found: ", input)
}

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop("Input is not a Seurat object.")
}

DefaultAssay(object) <- "RNA"

hvg <- VariableFeatures(object)

if (length(hvg) < 20) {
  stop("Too few variable features for PCA.")
}

dims_use <- 1:min(
  30,
  length(hvg) - 1
)

message("============================================================")
message("Running unintegrated PCA")
message("============================================================")

object <- RunPCA(
  object,
  features = hvg,
  npcs = 30,
  verbose = FALSE
)

# ============================================================
# ELBOW PLOT
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
  dims = dims_use,
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
  dims = dims_use,
  reduction = "pca",
  reduction.name = "umap.unintegrated",
  reduction.key = "unintegratedUMAP_",
  verbose = FALSE
)

# ============================================================
# FIGURES
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
message("Output: ", output)
message("============================================================")
