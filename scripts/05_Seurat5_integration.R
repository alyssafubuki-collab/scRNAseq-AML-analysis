# ============================================================
# 05_Seurat5_integration.R
#
# Seurat v5 CCA integration
#
# Based on the previous working integration workflow:
# split -> NormalizeData -> FindVariableFeatures ->
# ScaleData -> RunPCA -> IntegrateLayers(CCAIntegration)
#
# No SelectIntegrationFeatures()
# No scAnnoX
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
library(ggplot2)

set.seed(123)

# Avoid parallel workers consuming too much RAM on GitHub Actions
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

options(future.globals.maxSize = 4 * 1024^3)

# ============================================================
# PATHS
# ============================================================

project_dir <- get_project_dir()

input_file <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
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

# ============================================================
# LOAD OBJECT
# ============================================================

cat("\n============================================================\n")
cat("LOADING NORMALIZED OBJECT\n")
cat("============================================================\n\n")

object <- readRDS(input_file)

DefaultAssay(object) <- "RNA"

cat("Cells:", ncol(object), "\n")
cat("Features:", nrow(object), "\n")

if (!"sample" %in% colnames(object[[]])) {
  stop(
    "Metadata column 'sample' is missing from the Seurat object."
  )
}

cat("\nSamples:\n")
print(table(object$sample))

# ============================================================
# SPLIT RNA LAYERS
# ============================================================

cat("\n============================================================\n")
cat("SPLITTING RNA LAYERS BY SAMPLE\n")
cat("============================================================\n\n")

rna_layers_before <- Layers(object[["RNA"]])

cat("RNA layers before split:\n")
print(rna_layers_before)

# Only split if the assay is currently not sample-split.
#
# A normal unsplit Seurat v5 assay has layers such as:
# counts / data
#
# A split assay has:
# counts.SampleA / counts.SampleB / ...
# data.SampleA / data.SampleB / ...

already_split <- any(
  grepl("^counts\\.", rna_layers_before)
)

if (!already_split) {

  object[["RNA"]] <- split(
    object[["RNA"]],
    f = object$sample
  )

  cat("\nRNA layers after split:\n")
  print(Layers(object[["RNA"]]))

} else {

  cat(
    "\nRNA assay is already split by sample. ",
    "Skipping split().\n",
    sep = ""
  )
}

# ============================================================
# NORMALIZATION
# ============================================================

cat("\n============================================================\n")
cat("NORMALIZATION\n")
cat("============================================================\n\n")

object <- NormalizeData(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ============================================================
# VARIABLE FEATURES
# ============================================================

cat("\n============================================================\n")
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
  "Number of variable features:",
  length(hvg),
  "\n"
)

# ============================================================
# SCALE DATA
# ============================================================

cat("\n============================================================\n")
cat("SCALING\n")
cat("============================================================\n\n")

object <- ScaleData(
  object,
  features = hvg,
  verbose = FALSE
)

# ============================================================
# PCA
# ============================================================

cat("\n============================================================\n")
cat("PCA\n")
cat("============================================================\n\n")

object <- RunPCA(
  object,
  features = hvg,
  npcs = 30,
  verbose = FALSE
)

cat("PCA completed.\n")

# ============================================================
# DEFINE REFERENCE SAMPLES
# ============================================================

cat("\n============================================================\n")
cat("CCA REFERENCE SAMPLES\n")
cat("============================================================\n\n")

# The public GSE145410 dataset does not contain the
# pre/post-treatment timepoint structure used in the
# original CRCM workflow.
#
# Here DMSO_A and DMSO_B are used as the reference samples,
# corresponding to the control condition.

reference_samples <- c(
  "DMSO_A",
  "DMSO_B"
)

# Get the order of the sample-specific counts layers.
count_layers <- Layers(
  object[["RNA"]],
  search = "^counts\\."
)

if (length(count_layers) == 0) {
  stop(
    "No sample-specific counts layers were found after splitting."
  )
}

sample_order <- sub(
  "^counts\\.",
  "",
  count_layers
)

cat("Sample order used by Seurat:\n")
print(sample_order)

reference_indices <- which(
  sample_order %in% reference_samples
)

cat("\nReference samples:\n")
print(reference_samples)

cat("\nReference indices:\n")
print(reference_indices)

if (length(reference_indices) != length(reference_samples)) {
  stop(
    "Could not identify all reference samples. Expected: ",
    paste(reference_samples, collapse = ", "),
    ". Found in layers: ",
    paste(sample_order, collapse = ", ")
  )
}

# ============================================================
# CCA INTEGRATION
# ============================================================

cat("\n============================================================\n")
cat("CCA INTEGRATION\n")
cat("============================================================\n\n")

object <- IntegrateLayers(
  object = object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  reference = reference_indices,
  dims = 1:30,
  dims.to.integrate = 30,
  k.weight = 50,
  verbose = FALSE
)

cat("\nCCA integration completed successfully.\n")

# ============================================================
# INTEGRATED NEIGHBORS / CLUSTERS\n# ============================================================

cat("\n============================================================\n")
cat("INTEGRATED CLUSTERING\n")
cat("============================================================\n\n")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)

object <- FindClusters(
  object,
  resolution = 0.4,
  cluster.name = "cluster.cca",
  verbose = FALSE
)

cat("Integrated clustering completed.\n")

# ============================================================
# INTEGRATED UMAP
# ============================================================

cat("\n============================================================\n")
cat("INTEGRATED UMAP\n")
cat("============================================================\n\n")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:30,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  return.model = FALSE,
  verbose = FALSE
)

p_integrated <- DimPlot(
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
  plot = p_integrated,
  width = 10,
  height = 7,
  dpi = 300
)

p_cluster <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "cluster.cca",
  label = TRUE
) +
  ggtitle("Integrated UMAP - CCA clusters")

ggsave(
  filename = file.path(
    fig_dir,
    "05_UMAP_integrated_clusters.png"
  ),
  plot = p_cluster,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# UNINTEGRATED UMAP
# ============================================================

cat("\n============================================================\n")
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

# ============================================================
# SAVE
# ============================================================

cat("\n============================================================\n")
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

cat("\n============================================================\n")
cat("SCRIPT 05 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
