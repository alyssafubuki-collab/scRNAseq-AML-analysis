# ============================================================
# 05_Seurat5_integration.R
# Seurat v5 CCA integration
#
# Robust layer handling:
#   - joins existing split layers before re-splitting
#   - prevents repeated split() errors
#   - Seurat v5 workflow
#
# No scAnnoX
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
library(ggplot2)

project_dir <- get_project_dir()


# ============================================================
# Input object
# ============================================================

input_file <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
)

if (!file.exists(input_file)) {

  stop(
    "Input object not found: ",
    input_file
  )
}

message("============================================================")
message("Loading normalized Seurat object")
message("============================================================")

object <- readRDS(input_file)

if (!inherits(object, "Seurat")) {

  stop(
    "03_normalized.rds does not contain a Seurat object."
  )
}

message(
  "Cells: ",
  ncol(object)
)

message(
  "Features: ",
  nrow(object)
)


# ============================================================
# Check RNA assay
# ============================================================

if (!"RNA" %in% Assays(object)) {

  stop(
    "RNA assay not found in Seurat object."
  )
}

DefaultAssay(object) <- "RNA"

rna <- object[["RNA"]]


# ============================================================
# Inspect current layers
# ============================================================

message("============================================================")
message("RNA layers before integration")
message("============================================================")

rna_layers <- Layers(rna)

print(rna_layers)

message(
  "Number of RNA layers: ",
  length(rna_layers)
)


# ============================================================
# Check sample metadata
# ============================================================

if (!"sample" %in% colnames(object[[]])) {

  stop(
    "Metadata column 'sample' is missing."
  )
}

sample_values <- object$sample

if (any(is.na(sample_values))) {

  stop(
    "The 'sample' metadata contains NA values."
  )
}

message("============================================================")
message("Sample information")
message("============================================================")

print(
  table(object$sample)
)


# ============================================================
# Join existing split layers
#
# Seurat v5 does not allow split() on an assay whose layers
# are already split.
#
# We therefore normalize the layer state first.
# ============================================================

if (length(rna_layers) > 2) {

  message("============================================================")
  message("Existing split layers detected")
  message("============================================================")

  message(
    "Joining RNA layers before splitting by sample..."
  )

  object[["RNA"]] <- JoinLayers(
    object[["RNA"]]
  )

  rna_layers_after_join <- Layers(
    object[["RNA"]]
  )

  message(
    "RNA layers after JoinLayers():"
  )

  print(
    rna_layers_after_join
  )

} else {

  message(
    "RNA assay does not appear to contain multiple split layers."
  )
}


# ============================================================
# Split RNA assay by sample
# ============================================================

message("============================================================")
message("Splitting RNA assay by sample")
message("============================================================")

object[["RNA"]] <- split(
  object[["RNA"]],
  f = object$sample
)

message(
  "RNA layers after split:"
)

print(
  Layers(object[["RNA"]])
)


# ============================================================
# Variable features
# ============================================================

message("============================================================")
message("Finding variable features")
message("============================================================")

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)


# ============================================================
# Scaling
# ============================================================

message("============================================================")
message("Scaling data")
message("============================================================")

object <- ScaleData(
  object,
  verbose = FALSE
)


# ============================================================
# PCA
# ============================================================

message("============================================================")
message("Running PCA")
message("============================================================")

object <- RunPCA(
  object,
  features = VariableFeatures(object),
  npcs = 30,
  verbose = FALSE
)


# ============================================================
# CCA integration
# ============================================================

message("============================================================")
message("Running Seurat v5 CCA integration")
message("============================================================")

object <- IntegrateLayers(
  object = object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)


# ============================================================
# Neighbors
# ============================================================

message("============================================================")
message("Finding neighbors")
message("============================================================")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)


# ============================================================
# Clustering
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
message("Running UMAP")
message("============================================================")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:30,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  verbose = FALSE
)


# ============================================================
# Create output directories
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
# UMAP by sample
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
# UMAP by treatment
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

} else {

  warning(
    "Metadata column 'treatment' not found. ",
    "Treatment UMAP will not be generated."
  )
}


# ============================================================
# UMAP by cluster
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
# Save integrated object
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


# ============================================================
# Final validation
# ============================================================

if (!file.exists(output_file)) {

  stop(
    "Failed to save integrated Seurat object."
  )
}

message("============================================================")
message("Integration completed successfully")
message("============================================================")

message(
  "Output: ",
  output_file
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
