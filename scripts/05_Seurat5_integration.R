# ============================================================
# 05_Seurat5_integration.R
#
# Seurat v5 Harmony integration
#
# Workflow based directly on the previous working workflow:
#
# QC object
#   -> split RNA layers by sample
#   -> NormalizeData
#   -> FindVariableFeatures
#   -> ScaleData
#   -> RunPCA
#   -> IntegrateLayers / HarmonyIntegration
#   -> integrated UMAP
#   -> integrated clustering
#
# No SelectIntegrationFeatures()
# No scAnnoX
#
# ------------------------------------------------------------
# NOTE (see chat discussion): CCAIntegration reproducibly failed
# with "Error in data.use1[anchors1, ] : subscript out of bounds"
# inside Seurat's anchor-weighting step, both with the default
# all-pairs merge tree and with an explicit DMSO reference. This
# happens when a sample pair shares fewer mutual-nearest-neighbor
# anchors than k.weight. Given AZA (demethylating agent) and
# INCB059872 (LSD1 inhibitor) are both expected to cause large
# transcriptional shifts relative to DMSO, at least one pair
# likely lacks enough shared anchors for CCA regardless of
# k.weight or reference choice.
#
# Fix: use HarmonyIntegration instead. Harmony integrates on the
# PCA embedding directly with no anchor-finding step, so it is
# not subject to this failure mode and is standard practice for
# designs with strong expected treatment effects across samples.
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
library(ggplot2)

set.seed(123)

# ------------------------------------------------------------
# MEMORY / PARALLELIZATION
# ------------------------------------------------------------

if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

options(future.globals.maxSize = 4 * 1024^3)

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

project_dir <- get_project_dir()

input_file <- file.path(
  project_dir,
  "results",
  "objects",
  "02_QC_filtered.rds"
)

output_file <- file.path(
  project_dir,
  "results",
  "objects",
  "05_integrated.rds"
)

fig_dir <- file.path(
  project_dir,
  "results",
  "figures"
)

dir.create(
  file.path(project_dir, "results", "objects"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file
  )
}

# ------------------------------------------------------------
# LOAD QC OBJECT
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LOADING QC OBJECT\n")
cat("============================================================\n\n")

object <- readRDS(input_file)

DefaultAssay(object) <- "RNA"

cat(
  "Cells: ",
  ncol(object),
  "\n",
  sep = ""
)

cat(
  "Features: ",
  nrow(object),
  "\n",
  sep = ""
)

if (!"sample" %in% colnames(object[[]])) {
  stop(
    "Metadata column 'sample' is missing from the Seurat object."
  )
}

cat("\nSamples:\n")
print(table(object$sample))

# ------------------------------------------------------------
# CHECK INITIAL RNA LAYERS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INITIAL RNA LAYERS\n")
cat("============================================================\n\n")

initial_layers <- Layers(object[["RNA"]])

print(initial_layers)

# The object coming from script 02 should be unsplit.
# If it is already split, stop instead of silently modifying it.

if (any(grepl("^counts\\.", initial_layers))) {
  stop(
    paste0(
      "The input QC object is already split into sample-specific ",
      "layers. Expected an unsplit RNA assay before script 05."
    )
  )
}

# ------------------------------------------------------------
# SPLIT
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SPLITTING RNA LAYERS BY SAMPLE\n")
cat("============================================================\n\n")

object[["RNA"]] <- split(
  object[["RNA"]],
  f = object$sample
)

split_layers <- Layers(object[["RNA"]])

cat("RNA layers after split:\n")
print(split_layers)

# ------------------------------------------------------------
# CHECK SPLIT STRUCTURE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CHECKING SPLIT STRUCTURE\n")
cat("============================================================\n\n")

count_layers <- Layers(
  object[["RNA"]],
  search = "^counts\\."
)

cat("RNA layers after split:\n")
print(Layers(object[["RNA"]]))

cat("\nNumber of counts layers: ")
cat(length(count_layers))
cat("\n")

if (length(count_layers) != 8) {
  stop(
    "Expected 8 sample-specific counts layers, found ",
    length(count_layers),
    "."
  )
}

cat("\nSplit structure check: OK\n")

# ------------------------------------------------------------
# NORMALIZATION
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("NORMALIZATION\n")
cat("============================================================\n\n")

object <- NormalizeData(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("Normalization completed.\n")

# ------------------------------------------------------------
# CHECK NORMALIZED DATA LAYERS
# ------------------------------------------------------------

data_layers <- Layers(
  object[["RNA"]],
  search = "^data\\."
)

cat("\nNumber of data layers after normalization: ")
cat(length(data_layers))
cat("\n")

print(data_layers)

if (length(data_layers) != 8) {
  stop(
    "Expected 8 sample-specific data layers after normalization, found ",
    length(data_layers),
    "."
  )
}

cat("\nNormalized data layer check: OK\n")

# ------------------------------------------------------------
# VARIABLE FEATURES
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("HIGHLY VARIABLE FEATURES\n")
cat("============================================================\n\n")

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

hvg <- VariableFeatures(object)

cat(
  "Number of variable features: ",
  length(hvg),
  "\n",
  sep = ""
)

if (length(hvg) != 2000) {
  warning(
    "Expected 2000 variable features but found ",
    length(hvg),
    "."
  )
}

# ------------------------------------------------------------
# SCALE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SCALING\n")
cat("============================================================\n\n")

object <- ScaleData(
  object,
  features = hvg,
  verbose = FALSE
)

cat("Scaling completed.\n")

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PCA\n")
cat("============================================================\n\n")

object <- RunPCA(
  object,
  features = hvg,
  npcs = 30,
  verbose = FALSE
)

cat("PCA completed.\n")

# ------------------------------------------------------------
# CCA PRE-FLIGHT DIAGNOSTICS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INTEGRATION PRE-FLIGHT DIAGNOSTICS\n")
cat("============================================================\n\n")

cat("RNA layers:\n")
print(Layers(object[["RNA"]]))

cat("\nPCA dimensions:\n")
print(dim(Embeddings(object[["pca"]])))

cat("\nNumber of PCA cells:\n")
cat(nrow(Embeddings(object[["pca"]])))
cat("\n")

cat("\nNumber of object cells:\n")
cat(ncol(object))
cat("\n")

# Check that PCA contains exactly the cells of the Seurat object.

pca_cells <- rownames(
  Embeddings(object[["pca"]])
)

object_cells <- colnames(object)

if (!setequal(pca_cells, object_cells)) {
  stop(
    "PCA cell names do not match Seurat object cell names."
  )
}

cat("\nPCA cell identity check: OK\n")

# Check that each sample has cells.

sample_counts <- table(object$sample)

cat("\nCells per sample:\n")
print(sample_counts)

if (any(sample_counts == 0)) {
  stop(
    "At least one sample contains zero cells."
  )
}

# Check sample-specific layer names.

sample_order <- sub(
  "^counts\\.",
  "",
  count_layers
)

cat("\nSample order in RNA layers:\n")
print(sample_order)

if (!setequal(sample_order, unique(object$sample))) {
  stop(
    "RNA layer sample names do not match object$sample."
  )
}

cat("\nSample/layer identity check: OK\n")

# Check for duplicate cell barcodes across samples. 10X barcodes can
# repeat across GEM wells; if any upstream step merged objects without
# enforcing unique cell names, integration indexing can break silently.

n_dup_cells <- sum(duplicated(colnames(object)))

cat("\nDuplicate cell barcode check: ")
cat(n_dup_cells)
cat(" duplicates found\n")

if (n_dup_cells > 0) {
  stop(
    "Found ",
    n_dup_cells,
    " duplicate cell barcodes across samples. ",
    "Cell names must be unique before integration."
  )
}

# ------------------------------------------------------------
# HARMONY INTEGRATION
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("HARMONY INTEGRATION\n")
cat("============================================================\n\n")

# CCAIntegration (with and without an explicit reference) reproducibly
# hit "Error in data.use1[anchors1, ] : subscript out of bounds" inside
# Seurat's anchor-weighting step (FindIntegrationMatrix). This happens
# when a pair of samples being merged shares fewer mutual-nearest-
# neighbor anchors than k.weight. Given the biology here -- AZA
# (a demethylating agent) and INCB059872 (an LSD1 inhibitor) are both
# expected to cause large transcriptional shifts relative to DMSO --
# at least one sample pair likely doesn't have enough shared anchors
# for CCA, regardless of k.weight or reference choice.
#
# Harmony integrates directly on the PCA embedding with no anchor-
# finding step, so it isn't subject to this failure mode and is a
# standard, robust choice for designs with strong expected treatment
# effects across many samples.

if (!requireNamespace("harmony", quietly = TRUE)) {
  stop(
    "Package 'harmony' is required for HarmonyIntegration but is not ",
    "installed. Install it with install.packages('harmony')."
  )
}

cat(
  "Running HarmonyIntegration across all 8 samples.\n"
)

cat(
  "Using dimensions 1:30.\n\n"
)

object <- IntegrateLayers(
  object = object,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  dims = 1:30,
  verbose = FALSE
)

cat("\n")
cat("Harmony integration completed successfully.\n")

# ------------------------------------------------------------
# CHECK INTEGRATED REDUCTION
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CHECKING INTEGRATED REDUCTION\n")
cat("============================================================\n\n")

if (!"integrated.harmony" %in% Reductions(object)) {
  stop(
    "integrated.harmony reduction was not created."
  )
}

integrated_embeddings <- Embeddings(
  object[["integrated.harmony"]]
)

cat(
  "Integrated reduction dimensions: ",
  nrow(integrated_embeddings),
  " cells x ",
  ncol(integrated_embeddings),
  " dimensions\n",
  sep = ""
)

if (nrow(integrated_embeddings) != ncol(object)) {
  stop(
    "Number of cells in integrated.harmony does not match Seurat object."
  )
}

cat("Integrated reduction check: OK\n")

# ------------------------------------------------------------
# INTEGRATED NEIGHBORS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INTEGRATED NEIGHBORS\n")
cat("============================================================\n\n")

object <- FindNeighbors(
  object,
  reduction = "integrated.harmony",
  dims = 1:30,
  verbose = FALSE
)

# ------------------------------------------------------------
# INTEGRATED CLUSTERS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INTEGRATED CLUSTERS\n")
cat("============================================================\n\n")

object <- FindClusters(
  object,
  resolution = 0.4,
  cluster.name = "cluster.harmony",
  verbose = FALSE
)

cat(
  "Number of Harmony clusters: ",
  length(unique(object$cluster.harmony)),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# INTEGRATED UMAP
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INTEGRATED UMAP\n")
cat("============================================================\n\n")

object <- RunUMAP(
  object,
  reduction = "integrated.harmony",
  dims = 1:30,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  return.model = FALSE,
  verbose = FALSE
)

p_integrated_sample <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "sample"
) +
  ggtitle("Integrated UMAP - samples")

ggsave(
  filename = file.path(
    fig_dir,
    "05_UMAP_integrated_sample.png"
  ),
  plot = p_integrated_sample,
  width = 10,
  height = 7,
  dpi = 300
)

p_integrated_cluster <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "cluster.harmony",
  label = TRUE
) +
  ggtitle("Integrated UMAP - Harmony clusters")

ggsave(
  filename = file.path(
    fig_dir,
    "05_UMAP_integrated_clusters.png"
  ),
  plot = p_integrated_cluster,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# UNINTEGRATED UMAP
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("UNINTEGRATED UMAP\n")
cat("============================================================\n\n")

object <- RunUMAP(
  object,
  dims = 1:30,
  reduction.name = "umap.unintegrated",
  reduction.key = "unintegratedUMAP_",
  return.model = FALSE,
  verbose = FALSE
)

p_unintegrated <- DimPlot(
  object,
  reduction = "umap.unintegrated",
  group.by = "sample"
) +
  ggtitle("Unintegrated UMAP - samples")

ggsave(
  filename = file.path(
    fig_dir,
    "05_UMAP_unintegrated_sample.png"
  ),
  plot = p_unintegrated,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SAVING INTEGRATED OBJECT\n")
cat("============================================================\n\n")

saveRDS(
  object,
  output_file
)

cat(
  "Saved integrated object to:\n",
  output_file,
  "\n"
)

cat("\n")
cat("============================================================\n")
cat("SCRIPT 05 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
