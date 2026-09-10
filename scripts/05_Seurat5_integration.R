# ============================================================
# 05_Seurat5_integration.R
# Seurat v5 CCA integration
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

# ------------------------------------------------------------
# Split RNA assay into sample-specific layers.
# This is the Seurat v5 layer-based workflow.
# ------------------------------------------------------------

object[["RNA"]] <- split(
  object[["RNA"]],
  f = object$sample
)

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

object <- ScaleData(
  object,
  verbose = FALSE
)

object <- RunPCA(
  object,
  features = VariableFeatures(object),
  npcs = 30,
  verbose = FALSE
)

# ------------------------------------------------------------
# CCA integration
# ------------------------------------------------------------

object <- IntegrateLayers(
  object = object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
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
  reduction = "integrated.cca",
  dims = 1:30,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  verbose = FALSE
)

p1 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "sample"
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "integration",
    "UMAP_integrated_samples.png"
  ),
  p1,
  width = 8,
  height = 6,
  dpi = 300
)

p2 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "treatment"
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "integration",
    "UMAP_integrated_treatment.png"
  ),
  p2,
  width = 8,
  height = 6,
  dpi = 300
)

p3 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "integration",
    "UMAP_integrated_clusters.png"
  ),
  p3,
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
    "05_integrated.rds"
  )
)

message("Seurat v5 CCA integration completed.")
