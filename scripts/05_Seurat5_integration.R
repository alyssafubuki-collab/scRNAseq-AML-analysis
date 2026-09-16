# ============================================================
# 05_Seurat5_integration.R
#
# Seurat v5 CCA anchor-based integration
#
# Based directly on the original integration workflow:
#
#   1. Load normalized object
#   2. Split RNA assay by sample
#   3. Normalize sample-specific layers
#   4. Select variable features
#   5. Scale HVGs
#   6. PCA
#   7. CCA integration using reference samples
#   8. Integrated UMAP
#   9. Integrated neighbors / clusters
#  10. Unintegrated UMAP for comparison
#  11. Join RNA layers
#  12. Save integrated object
#
# Dataset:
# GSE145410
#
# Reference:
# DMSO_A + DMSO_B
#
# IMPORTANT:
# The reference is based on the biological control condition
# rather than on arbitrary layer ordering.
# ============================================================


# ============================================================
# 0. INITIALIZATION
# ============================================================

source("R/functions.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

set.seed(123)

project_dir <- get_project_dir()

options(
  future.globals.maxSize = 4 * 1024^3
)


# ============================================================
# 1. PATHS
# ============================================================

input <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
)

output <- file.path(
  project_dir,
  "results",
  "objects",
  "05_integrated.rds"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "integration"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. LOAD OBJECT
# ============================================================

message("============================================================")
message("05 - SEURAT V5 CCA INTEGRATION")
message("============================================================")

if (!file.exists(input)) {
  stop(
    "Input file not found: ",
    input
  )
}

message("Loading normalized Seurat object")

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop(
    "Input file does not contain a Seurat object."
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
# 3. CHECK SAMPLE METADATA
# ============================================================

if (!"sample" %in% colnames(object@meta.data)) {
  stop(
    "Metadata column 'sample' is missing."
  )
}

object$sample <- as.character(object$sample)

samples <- unique(object$sample)

message("============================================================")
message("SAMPLES")
message("============================================================")

print(
  table(object$sample)
)

message("Number of samples: ", length(samples))

if (length(samples) != 8) {
  warning(
    "Expected 8 GSE145410 samples, found ",
    length(samples),
    "."
  )
}


# ============================================================
# 4. CLEAN EXISTING RNA LAYERS
# ============================================================
#
# The object produced by script 03 should contain:
#
#   counts
#   data
#
# If sample-specific layers already exist, they are first joined.
#
# This prevents the repeated:
#
#   "The following layers are already split"
#
# error encountered previously.
# ============================================================

message("============================================================")
message("INITIAL RNA LAYERS")
message("============================================================")

initial_layers <- Layers(
  object[["RNA"]]
)

print(initial_layers)


sample_layer_pattern <- paste0(
  "\\.",
  paste(
    gsub(
      "([.|()\\+])",
      "\\\\\\1",
      samples
    ),
    collapse = "|"
  ),
  "$"
)

has_sample_layers <- any(
  grepl(
    "^counts\\.",
    initial_layers
  )
) ||
  any(
    grepl(
      "^data\\.",
      initial_layers
    )
  )

if (has_sample_layers) {

  message(
    "Existing sample-specific RNA layers detected."
  )

  message(
    "Joining RNA layers before re-splitting."
  )

  object[["RNA"]] <- JoinLayers(
    object[["RNA"]]
  )

  gc()

} else {

  message(
    "RNA assay contains unsplit counts/data layers."
  )
}

message("RNA layers after JoinLayers check:")

print(
  Layers(object[["RNA"]]
  )
)


# ============================================================
# 5. SPLIT RNA ASSAY BY SAMPLE
# ============================================================

message("============================================================")
message("SPLITTING RNA ASSAY BY SAMPLE")
message("============================================================")

object[["RNA"]] <- split(
  object[["RNA"]],
  f = object$sample
)

gc()

split_layers <- Layers(
  object[["RNA"]]
)

message("RNA layers after split:")

print(split_layers)


# ============================================================
# 6. VERIFY 8 COUNTS + 8 DATA LAYERS
# ============================================================

count_layers <- split_layers[
  grepl(
    "^counts\\.",
    split_layers
  )
]

data_layers <- split_layers[
  grepl(
    "^data\\.",
    split_layers
  )
]

message(
  "Counts layers: ",
  length(count_layers)
)

message(
  "Data layers: ",
  length(data_layers)
)

if (length(count_layers) != length(samples)) {
  stop(
    "Unexpected number of counts layers. Expected ",
    length(samples),
    ", found ",
    length(count_layers),
    "."
  )
}

if (length(data_layers) != length(samples)) {
  stop(
    "Unexpected number of data layers. Expected ",
    length(samples),
    ", found ",
    length(data_layers),
    "."
  )
}


# ============================================================
# 7. NORMALIZATION
# ============================================================
#
# Normalization is repeated here intentionally.
#
# This follows the original workflow where normalization is
# performed after splitting the assay into sample-specific layers.
# ============================================================

message("============================================================")
message("NORMALIZATION")
message("============================================================")

object <- NormalizeData(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

gc()


# ============================================================
# 8. SELECT INTEGRATION FEATURES
# ============================================================
#
# We deliberately avoid manually calling SelectIntegrationFeatures()
# here.
#
# With Seurat v5 layer-based integration, FindVariableFeatures()
# operates on the split layers and provides the feature set used
# by the PCA/integration workflow.
#
# This avoids the S4 coercion error encountered previously from
# calling SelectIntegrationFeatures() directly on this object.
# ============================================================

message("============================================================")
message("HIGHLY VARIABLE FEATURES")
message("============================================================")

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

hvg <- VariableFeatures(
  object
)

message(
  "Number of HVGs: ",
  length(hvg)
)

if (length(hvg) < 1000) {
  stop(
    "Too few variable features detected: ",
    length(hvg)
  )
}


# ============================================================
# 9. SCALE HVGs
# ============================================================

message("============================================================")
message("SCALING HVGs")
message("============================================================")

object <- ScaleData(
  object,
  features = hvg,
  verbose = TRUE
)

gc()


# ============================================================
# 10. PCA
# ============================================================

message("============================================================")
message("PCA")
message("============================================================")

object <- RunPCA(
  object,
  features = hvg,
  npcs = 30,
  verbose = FALSE
)

gc()

if (!"pca" %in% Reductions(object)) {
  stop(
    "PCA reduction was not successfully created."
  )
}

message(
  "PCA completed."
)


# ============================================================
# 11. DEFINE REFERENCE SAMPLES
# ============================================================
#
# GSE145410 samples:
#
# DMSO_A
# DMSO_B
# INCB059872_A
# INCB059872_B
# AZA_A
# AZA_B
# INCB059872_AZA_A
# INCB059872_AZA_B
#
# DMSO_A and DMSO_B are the untreated/control samples.
#
# They are therefore used as the integration reference.
# ============================================================

reference_samples <- c(
  "DMSO_A",
  "DMSO_B"
)

missing_reference <- setdiff(
  reference_samples,
  samples
)

if (length(missing_reference) > 0) {
  stop(
    "Reference sample(s) missing: ",
    paste(
      missing_reference,
      collapse = ", "
    )
  )
}

message("============================================================")
message("INTEGRATION REFERENCE")
message("============================================================")

message(
  "Reference samples: ",
  paste(
    reference_samples,
    collapse = ", "
  )
)

message("All samples:")

print(
  samples
)


# ============================================================
# 12. DETERMINE REFERENCE INDICES
# ============================================================
#
# IntegrateLayers() expects indices corresponding to the dataset
# order represented by the split layers.
#
# We obtain the dataset/sample order directly from the RNA layer
# names rather than relying on alphabetical factor levels.
# ============================================================

counts_layers <- Layers(
  object[["RNA"]],
  search = "^counts\\."
)

if (length(counts_layers) != length(samples)) {
  stop(
    "Could not determine the 8 sample-specific counts layers."
  )
}

layer_samples <- sub(
  "^counts\\.",
  "",
  counts_layers
)

message("Dataset order used by Seurat:")

print(
  data.frame(
    index = seq_along(layer_samples),
    sample = layer_samples
  )
)

reference_indices <- which(
  layer_samples %in% reference_samples
)

message(
  "Reference indices: ",
  paste(
    reference_indices,
    collapse = ", "
  )
)

if (length(reference_indices) != 2) {
  stop(
    "Expected exactly 2 reference datasets (DMSO_A and DMSO_B), found ",
    length(reference_indices),
    "."
  )
}


# ============================================================
# 13. CCA INTEGRATION
# ============================================================

message("============================================================")
message("CCA INTEGRATION")
message("============================================================")

message(
  "Using DMSO_A + DMSO_B as reference."
)

message(
  "Dimensions: 1:30"
)

message(
  "Method: CCAIntegration"
)

integration_error <- NULL

object <- tryCatch(

  {

    IntegrateLayers(
      object = object,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      reference = reference_indices,
      dims = 1:30,
      k.weight = 50,
      verbose = TRUE
    )

  },

  error = function(e) {

    integration_error <<- conditionMessage(e)

    NULL
  }
)

if (is.null(object)) {

  message("============================================================")
  message("CCA INTEGRATION FAILED")
  message("============================================================")

  message(
    integration_error
  )

  stop(
    "CCA integration failed. Original error: ",
    integration_error
  )
}

gc()

if (!"integrated.cca" %in% Reductions(object)) {
  stop(
    "CCA integration completed without creating ",
    "'integrated.cca'."
  )
}

message("CCA integration completed successfully.")


# ============================================================
# 14. INTEGRATED NEIGHBORS
# ============================================================

message("============================================================")
message("INTEGRATED NEIGHBORS")
message("============================================================")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)

gc()


# ============================================================
# 15. INTEGRATED CLUSTERING
# ============================================================

message("============================================================")
message("INTEGRATED CLUSTERING")
message("============================================================")

object <- FindClusters(
  object,
  resolution = 0.4,
  cluster.name = "cluster.cca",
  verbose = FALSE
)

gc()


# ============================================================
# 16. INTEGRATED UMAP
# ============================================================

message("============================================================")
message("INTEGRATED UMAP")
message("============================================================")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:30,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  verbose = FALSE
)

gc()


# ============================================================
# 17. UNINTEGRATED UMAP
# ============================================================
#
# This reproduces the comparison in the original workflow:
#
#   integrated UMAP
#   versus
#   non-integrated PCA/UMAP
#
# The unintegrated reduction remains useful for assessing whether
# integration substantially changes the structure.
# ============================================================

message("============================================================")
message("UNINTEGRATED UMAP")
message("============================================================")

object <- FindNeighbors(
  object,
  reduction = "pca",
  dims = 1:30,
  graph.name = "pca_snn_unintegrated",
  verbose = FALSE
)

object <- RunUMAP(
  object,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.unintegrated",
  reduction.key = "unintegratedUMAP_",
  verbose = FALSE
)

gc()


# ============================================================
# 18. INTEGRATED UMAP BY SAMPLE
# ============================================================

p1 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "sample"
) +
  ggtitle(
    "Integrated UMAP - Samples"
  )

ggsave(
  file.path(
    figure_dir,
    "UMAP_integrated_samples.png"
  ),
  p1,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 19. INTEGRATED UMAP BY TREATMENT
# ============================================================

if ("treatment" %in% colnames(object@meta.data)) {

  p2 <- DimPlot(
    object,
    reduction = "umap.integrated",
    group.by = "treatment"
  ) +
    ggtitle(
      "Integrated UMAP - Treatment"
    )

  ggsave(
    file.path(
      figure_dir,
      "UMAP_integrated_treatment.png"
    ),
    p2,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 20. INTEGRATED UMAP BY CLUSTER
# ============================================================

p3 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "cluster.cca",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Integrated UMAP - CCA clusters"
  )

ggsave(
  file.path(
    figure_dir,
    "UMAP_integrated_clusters.png"
  ),
  p3,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 21. UNINTEGRATED UMAP BY SAMPLE
# ============================================================

p4 <- DimPlot(
  object,
  reduction = "umap.unintegrated",
  group.by = "sample"
) +
  ggtitle(
    "Unintegrated UMAP - Samples"
  )

ggsave(
  file.path(
    figure_dir,
    "UMAP_unintegrated_samples.png"
  ),
  p4,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 22. JOIN RNA LAYERS
# ============================================================
#
# Integration has now been performed.
#
# We can safely rejoin the RNA layers for downstream annotation
# and differential expression.
# ============================================================

message("============================================================")
message("JOINING RNA LAYERS")
message("============================================================")

object <- JoinLayers(
  object[["RNA"]]
)

gc()

message("RNA layers after JoinLayers:")

print(
  Layers(
    object[["RNA"]]
  )
)


# ============================================================
# 23. FINAL VALIDATION
# ============================================================

message("============================================================")
message("FINAL VALIDATION")
message("============================================================")

if (!"integrated.cca" %in% Reductions(object)) {
  stop(
    "Final object does not contain integrated.cca."
  )
}

if (!"umap.integrated" %in% Reductions(object)) {
  stop(
    "Final object does not contain umap.integrated."
  )
}

if (!"cluster.cca" %in% colnames(object@meta.data)) {
  stop(
    "Final object does not contain cluster.cca."
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

message(
  "Reductions:"
)

print(
  Reductions(object)
)


# ============================================================
# 24. SAVE OBJECT
# ============================================================

message("============================================================")
message("SAVING OBJECT")
message("============================================================")

saveRDS(
  object,
  output
)

if (!file.exists(output)) {
  stop(
    "Output object was not created."
  )
}

message(
  "Output: ",
  output
)

message("============================================================")
message("05 - SEURAT V5 CCA INTEGRATION COMPLETED")
message("============================================================")
