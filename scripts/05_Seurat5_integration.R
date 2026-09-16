# ============================================================
# 05_Seurat5_integration.R
#
# Seurat v5 CCA integration
#
# Workflow:
# split -> NormalizeData -> FindVariableFeatures ->
# ScaleData -> RunPCA -> IntegrateLayers(CCAIntegration)
#
# No SelectIntegrationFeatures()
# No scAnnoX
# No reference integration in this version
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

# ------------------------------------------------------------
# LOAD
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LOADING NORMALIZED OBJECT\n")
cat("============================================================\n\n")

object <- readRDS(input_file)

DefaultAssay(object) <- "RNA"

cat("Cells: ", ncol(object), "\n", sep = "")
cat("Features: ", nrow(object), "\n", sep = "")

if (!"sample" %in% colnames(object[[]])) {
  stop(
    "Metadata column 'sample' is missing."
  )
}

cat("\nSamples:\n")
print(table(object$sample))

# ------------------------------------------------------------
# SPLIT
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SPLITTING RNA LAYERS BY SAMPLE\n")
cat("============================================================\n\n")

rna_layers <- Layers(object[["RNA"]])

cat("RNA layers before split:\n")
print(rna_layers)

already_split <- any(
  grepl("^counts\\.", rna_layers)
)

if (!already_split) {

  object[["RNA"]] <- split(
    object[["RNA"]],
    f = object$sample
  )

} else {

  cat("RNA assay already split. Skipping split().\n")
}

cat("\nRNA layers after split:\n")
print(Layers(object[["RNA"]]))

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

if (length(hvg) < 1000) {
  stop(
    "Too few variable features were selected: ",
    length(hvg)
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
# CCA INTEGRATION
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CCA INTEGRATION\n")
cat("============================================================\n\n")

cat("Running CCA integration across all samples.\n")
cat("Reference integration is disabled for this test.\n\n")

object <- IntegrateLayers(
  object = object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:30,
  dims.to.integrate = 30,
  k.weight = 50,
  verbose = FALSE
)

cat("\nCCA integration completed successfully.\n")

# ------------------------------------------------------------
# INTEGRATED NEIGHBORS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("INTEGRATED NEIGHBORS / CLUSTERS\n")
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

# ------------------------------------------------------------
# INTEGRATED UMAP
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
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
  group.by = "cluster.cca",
  label = TRUE
) +
  ggtitle("Integrated UMAP - CCA clusters")

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
