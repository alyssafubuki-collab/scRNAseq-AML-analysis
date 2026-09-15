# ============================================================
# 05_Seurat5_integration.R
# Seurat v5 CCA integration
#
# Layer policy:
#   1. 02 merges samples and joins layers.
#   2. 03 keeps RNA unsplit.
#   3. 05 is the ONLY script that splits RNA by sample.
#
# Resource-conscious settings:
#   - 2,000 HVGs
#   - 20 PCs
#   - 20 integrated dimensions
#   - k.weight = 50
#   - ScaleData only on HVGs
#
# No scAnnoX.
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
library(ggplot2)

project_dir <- get_project_dir()

input_file <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
)

if (!file.exists(input_file)) {
  stop("Input object not found: ", input_file)
}

object <- readRDS(input_file)

if (!inherits(object, "Seurat")) {
  stop("03_normalized.rds does not contain a Seurat object.")
}

DefaultAssay(object) <- "RNA"

if (!"sample" %in% colnames(object[[]])) {
  stop("Metadata column 'sample' is missing.")
}

if (any(is.na(object$sample))) {
  stop("The 'sample' metadata contains NA values.")
}

# ============================================================
# VALIDATE PRE-INTEGRATION LAYERS
# ============================================================

rna_layers <- Layers(object[["RNA"]])

message("============================================================")
message("RNA layers before integration")
message("============================================================")
print(rna_layers)

if (any(grepl("SeuratProject", rna_layers))) {
  stop(
    "Invalid SeuratProject layers detected before integration."
  )
}

if (length(rna_layers) > 2L) {
  stop(
    "RNA is already split before script 05. ",
    "This indicates a broken upstream layer state."
  )
}

# ============================================================
# SPLIT ONCE BY SAMPLE
# ============================================================

message("============================================================")
message("Splitting RNA assay by sample")
message("============================================================")

object[["RNA"]] <- split(
  object[["RNA"]],
  f = object$sample
)

split_layers <- Layers(object[["RNA"]])

print(split_layers)

if (length(split_layers) < 3L) {
  stop(
    "RNA assay was not split into sample-specific layers."
  )
}

# ============================================================
# HVG SELECTION
# ============================================================

message("============================================================")
message("Selecting 2,000 HVGs")
message("============================================================")

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# ============================================================
# SCALE ONLY HVGs
# ============================================================

message("============================================================")
message("Scaling 2,000 HVGs")
message("============================================================")

object <- ScaleData(
  object,
  features = VariableFeatures(object),
  verbose = FALSE
)

gc()

# ============================================================
# PCA
# ============================================================

message("============================================================")
message("Running PCA: 20 dimensions")
message("============================================================")

object <- RunPCA(
  object,
  features = VariableFeatures(object),
  npcs = 20,
  verbose = FALSE
)

# ============================================================
# CCA INTEGRATION
# ============================================================

message("============================================================")
message("Running Seurat v5 CCA integration")
message("============================================================")
message("Parameters: dims=1:20, dims.to.integrate=20, k.weight=50")

object <- IntegrateLayers(
  object = object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:20,
  dims.to.integrate = 20,
  k.weight = 50,
  verbose = FALSE
)

gc()

# ============================================================
# INTEGRATED NEIGHBORS
# ============================================================

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  verbose = FALSE
)

object <- FindClusters(
  object,
  resolution = 0.4,
  verbose = FALSE
)

# ============================================================
# INTEGRATED UMAP
# ============================================================

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  verbose = FALSE
)

# ============================================================
# FIGURES
# ============================================================

dir.create(
  file.path(
    project_dir,
    "figures",
    "integration"
  ),
  recursive = TRUE,
  showWarnings = FALSE
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

# ============================================================
# SAVE
# ============================================================

output_file <- file.path(
  project_dir,
  "results",
  "objects",
  "05_integrated.rds"
)

saveRDS(
  object,
  output_file
)

if (!file.exists(output_file)) {
  stop("Failed to save integrated Seurat object.")
}

message("============================================================")
message("Integration completed successfully")
message("Output: ", output_file)
message("Cells: ", ncol(object))
message(
  "Clusters: ",
  length(unique(object$seurat_clusters))
)
message("Available reductions:")
print(Reductions(object))
message("============================================================")
