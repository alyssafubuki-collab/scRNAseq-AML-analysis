# ============================================================
# 05_Seurat5_integration.R
#
# Seurat v5 CCA integration
#
# Pipeline:
#   clean RNA layers
#       ↓
#   split by sample
#       ↓
#   NormalizeData per layer
#       ↓
#   FindVariableFeatures
#       ↓
#   ScaleData on HVGs
#       ↓
#   PCA
#       ↓
#   CCAIntegration
#
# No scAnnoX.
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
# LOAD
# ============================================================

message("============================================================")
message("Loading normalized Seurat object")
message("============================================================")

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop(
    "03_normalized.rds does not contain a Seurat object."
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
# METADATA VALIDATION
# ============================================================

if (!"sample" %in% colnames(object[[]])) {
  stop(
    "Metadata column 'sample' is missing."
  )
}

if (any(is.na(object$sample))) {
  stop(
    "sample metadata contains NA values."
  )
}

sample_counts <- table(object$sample)

message("============================================================")
message("Cells per sample")
message("============================================================")

print(sample_counts)

if (length(sample_counts) != 8) {
  warning(
    "Expected 8 samples, found ",
    length(sample_counts)
  )
}

if (any(sample_counts < 50)) {
  stop(
    "At least one sample contains fewer than 50 cells. ",
    "CCA integration cannot be safely performed."
  )
}

# ============================================================
# INITIAL LAYERS
# ============================================================

message("============================================================")
message("RNA layers before integration")
message("============================================================")

initial_layers <- Layers(
  object[["RNA"]]
)

print(initial_layers)

if (any(grepl("SeuratProject", initial_layers))) {
  stop(
    "Invalid SeuratProject layers detected."
  )
}

# We expect the object to be joined before integration.
if (!all(
  c("counts", "data") %in% initial_layers
)) {

  message(
    "RNA layers are not in the expected joined state."
  )

  message(
    "Attempting JoinLayers() once before splitting."
  )

  object[["RNA"]] <- JoinLayers(
    object[["RNA"]]
  )
}

joined_layers <- Layers(
  object[["RNA"]]
)

message("RNA layers after validation:")
print(joined_layers)

if (!"counts" %in% joined_layers) {
  stop(
    "A counts layer is required before integration."
  )
}

# ============================================================
# SPLIT BY SAMPLE
# ============================================================

message("============================================================")
message("Splitting RNA assay by sample")
message("============================================================")

object[["RNA"]] <- split(
  object[["RNA"]],
  f = object$sample
)

split_layers <- Layers(
  object[["RNA"]]
)

message("RNA layers after split:")
print(split_layers)

data_layers <- Layers(
  object[["RNA"]],
  search = "^data\\."
)

counts_layers <- Layers(
  object[["RNA"]],
  search = "^counts\\."
)

message("Normalized data layers:")
print(data_layers)

message("Counts layers:")
print(counts_layers)

if (length(data_layers) != 8) {
  stop(
    "Expected 8 normalized data layers after split, found ",
    length(data_layers)
  )
}

if (length(counts_layers) != 8) {
  stop(
    "Expected 8 counts layers after split, found ",
    length(counts_layers)
  )
}

# ============================================================
# NORMALIZE EACH SAMPLE LAYER
#
# This follows the Seurat v5 integration workflow:
# split → NormalizeData → FindVariableFeatures → ScaleData → PCA
# ============================================================

message("============================================================")
message("Re-normalizing split RNA layers")
message("============================================================")

object <- NormalizeData(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

gc()

# ============================================================
# VARIABLE FEATURES
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

hvg <- VariableFeatures(object)

message(
  "Number of HVGs: ",
  length(hvg)
)

if (length(hvg) < 2000) {
  warning(
    "Fewer than 2,000 HVGs were returned: ",
    length(hvg)
  )
}

# ============================================================
# SCALE HVGs ONLY
# ============================================================

message("============================================================")
message("Scaling 2,000 HVGs")
message("============================================================")

object <- ScaleData(
  object,
  features = hvg,
  verbose = FALSE
)

gc()

# ============================================================
# VALIDATE SCALE LAYERS
# ============================================================

scale_layers <- Layers(
  object[["RNA"]],
  search = "^scale\\.data\\."
)

message("Scale layers:")
print(scale_layers)

if (length(scale_layers) != 8) {
  stop(
    "Expected 8 scale.data layers, found ",
    length(scale_layers)
  )
}

# ============================================================
# PCA
# ============================================================

message("============================================================")
message("Running PCA: 20 dimensions")
message("============================================================")

object <- RunPCA(
  object,
  features = hvg,
  npcs = 20,
  verbose = FALSE
)

if (!"pca" %in% Reductions(object)) {
  stop(
    "PCA reduction was not created."
  )
}

pca_cells <- nrow(
  Embeddings(
    object,
    reduction = "pca"
  )
)

if (pca_cells != ncol(object)) {
  stop(
    "PCA cell count does not match Seurat object."
  )
}

# ============================================================
# CCA INTEGRATION
# ============================================================

message("============================================================")
message("Running Seurat v5 CCA integration")
message("============================================================")

message(
  "Parameters:",
  " dims=1:20",
  ", k.filter=50",
  ", k.weight=50"
)

gc()

object <- tryCatch(

  IntegrateLayers(
    object = object,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    assay = "RNA",
    features = hvg,
    layers = data_layers,
    scale.layer = scale_layers,
    dims = 1:20,
    k.filter = 50,
    k.weight = 50,
    verbose = TRUE
  ),

  error = function(e) {

    message("============================================================")
    message("CCA INTEGRATION FAILED")
    message("============================================================")

    message(
      "Error message:"
    )

    message(
      conditionMessage(e)
    )

    message("------------------------------------------------------------")
    message("RNA layers:")
    print(
      Layers(object[["RNA"]])
    )

    message("------------------------------------------------------------")
    message("Data layers:")
    print(
      Layers(
        object[["RNA"]],
        search = "^data\\."
      )
    )

    message("------------------------------------------------------------")
    message("Scale layers:")
    print(
      Layers(
        object[["RNA"]],
        search = "^scale\\.data\\."
      )
    )

    message("------------------------------------------------------------")
    message("PCA:")
    print(
      object[["pca"]]
    )

    stop(
      "Seurat CCAIntegration failed: ",
      conditionMessage(e)
    )
  }
)

# ============================================================
# VALIDATE CCA
# ============================================================

if (!"integrated.cca" %in% Reductions(object)) {
  stop(
    "CCA integration did not create integrated.cca."
  )
}

cca_embeddings <- Embeddings(
  object,
  reduction = "integrated.cca"
)

message("============================================================")
message("CCA integration successful")
message("============================================================")

message(
  "Integrated cells: ",
  nrow(cca_embeddings)
)

message(
  "Integrated dimensions: ",
  ncol(cca_embeddings)
)

if (nrow(cca_embeddings) != ncol(object)) {
  stop(
    "integrated.cca contains a different number of cells ",
    "than the Seurat object."
  )
}

# ============================================================
# NEIGHBORS
# ============================================================

message("============================================================")
message("Finding integrated neighbors")
message("============================================================")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  verbose = FALSE
)

# ============================================================
# CLUSTERS
# ============================================================

message("============================================================")
message("Finding clusters")
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
message("Running integrated UMAP")
message("============================================================")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  verbose = FALSE
)

# ============================================================
# FIGURE DIRECTORY
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

# ============================================================
# UMAP BY SAMPLE
# ============================================================

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

# ============================================================
# UMAP BY TREATMENT
# ============================================================

if ("treatment" %in% colnames(object[[]])) {

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
}

# ============================================================
# UMAP BY CLUSTER
# ============================================================

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

output <- file.path(
  project_dir,
  "results",
  "objects",
  "05_integrated.rds"
)

saveRDS(
  object,
  output
)

message("============================================================")
message("Integration completed successfully")
message("============================================================")

message(
  "Output: ",
  output
)

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
  "Available reductions:"
)

print(
  Reductions(object)
)
