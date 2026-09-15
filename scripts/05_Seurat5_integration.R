# ============================================================
# 05_Seurat5_integration.R
# GSE145410 AML scRNA-seq
#
# Seurat v5 CCA integration - V2
#
# Input:
#   results/objects/03_normalized.rds
#
# Output:
#   results/objects/05_integrated.rds
#   results/tables/05_*.csv
#   figures/05_integration/*
#
# IMPORTANT
# ----------
# - Seurat v5 layer-based integration
# - CCAIntegration
# - Explicit balanced sample tree
# - k.weight = 20
# - k.filter = 50
# - 20 PCs
# - 2,000 HVGs
# - NO scAnnoX
# - NO manual Harmony integration
#
# ============================================================


# ============================================================
# 0. OPTIONS
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================================\n")
cat("05 - Seurat v5 CCA integration V2\n")
cat("============================================================\n")
cat("\n")


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})


# ============================================================
# 2. PARAMETERS
# ============================================================

INPUT_FILE <- "results/objects/03_normalized.rds"

OUTPUT_FILE <- "results/objects/05_integrated.rds"

N_HVG <- 2000

N_PCS <- 20

CCA_K_FILTER <- 50

CCA_K_WEIGHT <- 20

CLUSTER_RESOLUTION <- 0.4


# ============================================================
# 3. OUTPUT DIRECTORIES
# ============================================================

dir.create(
  "results/objects",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/05_integration",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. CHECK INPUT
# ============================================================

if (!file.exists(INPUT_FILE)) {

  stop(
    paste0(
      "\nERROR: Input file does not exist:\n",
      INPUT_FILE,
      "\n\n",
      "Run 03_normalization_HVG.R first.\n"
    )
  )
}


# ============================================================
# 5. LOAD OBJECT
# ============================================================

cat("Loading input object:\n")
cat("  ", INPUT_FILE, "\n\n", sep = "")

object <- readRDS(INPUT_FILE)


if (!inherits(object, "Seurat")) {

  stop(
    "ERROR: Input object is not a Seurat object."
  )
}


if (!"RNA" %in% Assays(object)) {

  stop(
    "ERROR: RNA assay is missing."
  )
}


if (!"sample" %in% colnames(object[[]])) {

  stop(
    "ERROR: metadata column 'sample' is missing."
  )
}


# ============================================================
# 6. BASIC OBJECT INFORMATION
# ============================================================

cat("============================================================\n")
cat("INPUT OBJECT\n")
cat("============================================================\n\n")

cat("Cells :", ncol(object), "\n")
cat("Genes :", nrow(object), "\n\n")

cat("Assays:\n")
print(Assays(object))

cat("\nMetadata columns:\n")
print(colnames(object[[]]))


# ============================================================
# 7. SAMPLE INFORMATION
# ============================================================

sample_ids <- sort(
  unique(
    as.character(object$sample)
  )
)

cat("\n")
cat("============================================================\n")
cat("SAMPLES\n")
cat("============================================================\n\n")

cat("Number of samples:", length(sample_ids), "\n\n")

print(sample_ids)


# ============================================================
# 8. EXPECTED SAMPLE ORDER
# ============================================================
#
# This order is used to construct an explicit balanced
# integration tree.
#
# 1  DMSO_A
# 2  DMSO_B
# 3  INCB059872_A
# 4  INCB059872_B
# 5  AZA_A
# 6  AZA_B
# 7  INCB059872_AZA_A
# 8  INCB059872_AZA_B
#
# ============================================================

expected_samples <- c(
  "DMSO_A",
  "DMSO_B",
  "INCB059872_A",
  "INCB059872_B",
  "AZA_A",
  "AZA_B",
  "INCB059872_AZA_A",
  "INCB059872_AZA_B"
)


# ============================================================
# 9. VALIDATE SAMPLE SET
# ============================================================

missing_samples <- setdiff(
  expected_samples,
  sample_ids
)

extra_samples <- setdiff(
  sample_ids,
  expected_samples
)


if (length(missing_samples) > 0) {

  stop(
    paste0(
      "\nERROR: Expected sample(s) missing:\n",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}


if (length(extra_samples) > 0) {

  stop(
    paste0(
      "\nERROR: Unexpected sample(s) detected:\n",
      paste(
        extra_samples,
        collapse = ", "
      )
    )
  )
}


# Force exact expected order
sample_ids <- expected_samples


# ============================================================
# 10. CELLS PER SAMPLE
# ============================================================

cells_per_sample <- table(
  factor(
    object$sample,
    levels = sample_ids
  )
)

cat("\n")
cat("Cells per sample:\n\n")

print(cells_per_sample)

write.csv(
  as.data.frame(cells_per_sample),
  "results/tables/05_cells_per_sample.csv",
  row.names = FALSE
)


# ============================================================
# 11. INITIAL RNA LAYERS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("INITIAL RNA LAYERS\n")
cat("============================================================\n\n")

initial_layers <- Layers(
  object[["RNA"]]
)

print(initial_layers)


# ============================================================
# 12. JOIN PREVIOUSLY SPLIT LAYERS IF NECESSARY
# ============================================================
#
# We start from 03_normalized.rds.
#
# If the object contains sample-specific layers already,
# JoinLayers() is used before creating the clean split required
# for this integration.
#
# ============================================================

sample_layer_pattern <- paste0(
  "\\.(",
  paste(
    sample_ids,
    collapse = "|"
  ),
  ")$"
)

has_sample_specific_layers <- any(
  grepl(
    sample_layer_pattern,
    initial_layers
  )
)


if (has_sample_specific_layers) {

  cat("\n")
  cat("Sample-specific layers detected.\n")
  cat("Joining RNA layers before splitting again...\n\n")

  object[["RNA"]] <- JoinLayers(
    object[["RNA"]]
  )

} else {

  cat("\n")
  cat("No sample-specific layers detected.\n")
  cat("RNA assay is already suitable for splitting.\n")
}


# ============================================================
# 13. RNA LAYERS AFTER JOIN
# ============================================================

cat("\n")
cat("RNA layers after JoinLayers check:\n\n")

layers_after_join <- Layers(
  object[["RNA"]]
)

print(layers_after_join)


# ============================================================
# 14. SPLIT RNA BY SAMPLE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SPLITTING RNA ASSAY BY SAMPLE\n")
cat("============================================================\n\n")

object[["RNA"]] <- split(
  object[["RNA"]],
  f = factor(
    object$sample,
    levels = sample_ids
  )
)


# ============================================================
# 15. VALIDATE SPLIT
# ============================================================

split_layers <- Layers(
  object[["RNA"]]
)

cat("RNA layers after split:\n\n")
print(split_layers)


counts_layers <- grep(
  "^counts\\.",
  split_layers,
  value = TRUE
)

data_layers <- grep(
  "^data\\.",
  split_layers,
  value = TRUE
)


cat("\n")
cat("Counts layers:", length(counts_layers), "\n")
print(counts_layers)

cat("\n")
cat("Data layers:", length(data_layers), "\n")
print(data_layers)


if (length(counts_layers) != length(sample_ids)) {

  stop(
    paste0(
      "\nERROR: Expected ",
      length(sample_ids),
      " counts layers but found ",
      length(counts_layers),
      ".\n\n",
      "RNA layers are:\n",
      paste(
        split_layers,
        collapse = "\n"
      )
    )
  )
}


if (length(data_layers) != length(sample_ids)) {

  stop(
    paste0(
      "\nERROR: Expected ",
      length(sample_ids),
      " data layers but found ",
      length(data_layers),
      ".\n\n",
      "RNA layers are:\n",
      paste(
        split_layers,
        collapse = "\n"
      )
    )
  }


# ============================================================
# 16. NORMALIZATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("NORMALIZATION\n")
cat("============================================================\n\n")

object <- NormalizeData(
  object,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)


# ============================================================
# 17. VARIABLE FEATURES
# ============================================================

cat("\n")
cat("============================================================\n")
cat("VARIABLE FEATURES\n")
cat("============================================================\n\n")

object <- FindVariableFeatures(
  object,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = N_HVG,
  verbose = FALSE
)

hvg <- VariableFeatures(
  object
)


cat(
  "Number of HVGs:",
  length(hvg),
  "\n"
)


if (length(hvg) == 0) {

  stop(
    "ERROR: No variable features were identified."
  )
}


if (length(hvg) < 1000) {

  warning(
    "Fewer than 1,000 HVGs were identified: ",
    length(hvg)
  )
}


write.csv(
  data.frame(
    rank = seq_along(hvg),
    feature = hvg
  ),
  "results/tables/05_variable_features.csv",
  row.names = FALSE
)


# ============================================================
# 18. SCALE DATA
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCALING VARIABLE FEATURES\n")
cat("============================================================\n\n")

object <- ScaleData(
  object,
  assay = "RNA",
  features = hvg,
  verbose = FALSE
)


# ============================================================
# 19. CHECK SCALE DATA
# ============================================================

cat("\n")
cat("RNA layers after ScaleData:\n\n")

layers_after_scale <- Layers(
  object[["RNA"]]
)

print(layers_after_scale)


scale_layers <- grep(
  "^scale\\.data",
  layers_after_scale,
  value = TRUE
)


cat(
  "\nScale layer(s) detected:",
  length(scale_layers),
  "\n"
)

print(scale_layers)


if (length(scale_layers) == 0) {

  stop(
    paste0(
      "\nERROR: No scale.data layer was created.\n\n",
      "Current RNA layers:\n",
      paste(
        layers_after_scale,
        collapse = "\n"
      )
    )
  )
}


# ============================================================
# 20. PCA
# ============================================================

cat("\n")
cat("============================================================\n")
cat("PCA\n")
cat("============================================================\n\n")

object <- RunPCA(
  object,
  assay = "RNA",
  features = hvg,
  npcs = N_PCS,
  verbose = FALSE
)


# ============================================================
# 21. CHECK PCA
# ============================================================

if (!"pca" %in% Reductions(object)) {

  stop(
    "ERROR: PCA reduction was not generated."
  )
}


pca_embeddings <- Embeddings(
  object,
  reduction = "pca"
)

pca_dims <- ncol(
  pca_embeddings
)

cat(
  "PCA dimensions:",
  pca_dims,
  "\n"
)


if (pca_dims < N_PCS) {

  warning(
    "Only ",
    pca_dims,
    " PCA dimensions were generated."
  )
}


# ============================================================
# 22. ELBOW PLOT
# ============================================================

pdf(
  "figures/05_integration/05_elbowplot.pdf",
  width = 7,
  height = 5
)

print(
  ElbowPlot(
    object,
    reduction = "pca",
    ndims = min(
      N_PCS,
      pca_dims
    )
  )
)

dev.off()


# ============================================================
# 23. EXPLICIT BALANCED SAMPLE TREE
# ============================================================
#
# Dataset numbering:
#
#   1 = DMSO_A
#   2 = DMSO_B
#   3 = INCB059872_A
#   4 = INCB059872_B
#   5 = AZA_A
#   6 = AZA_B
#   7 = INCB059872_AZA_A
#   8 = INCB059872_AZA_B
#
# Integration steps:
#
#   row 1: dataset 1 + dataset 2
#   row 2: dataset 3 + dataset 4
#   row 3: dataset 5 + dataset 6
#   row 4: dataset 7 + dataset 8
#
#   row 5: result 1 + result 2
#   row 6: result 3 + result 4
#
#   row 7: result 5 + result 6
#
# This avoids the automatically generated unbalanced tree
# that previously produced:
#
#   Merging dataset 5 into 2
#
# immediately before the failing integration-vector step.
#
# ============================================================

sample_tree <- matrix(
  c(
    -1, -2,
    -3, -4,
    -5, -6,
    -7, -8,
     1,  2,
     3,  4,
     5,  6
  ),
  ncol = 2,
  byrow = TRUE
)


cat("\n")
cat("============================================================\n")
cat("EXPLICIT SAMPLE TREE\n")
cat("============================================================\n\n")

print(sample_tree)

cat("\nDataset mapping:\n")

for (i in seq_along(sample_ids)) {

  cat(
    "  ",
    i,
    " = ",
    sample_ids[i],
    "\n",
    sep = ""
  )
}


write.csv(
  sample_tree,
  "results/tables/05_sample_tree.csv",
  row.names = FALSE
)


# ============================================================
# 24. CCA INTEGRATION PARAMETERS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CCA INTEGRATION PARAMETERS\n")
cat("============================================================\n\n")

cat("Method            : CCAIntegration\n")
cat("Assay             : RNA\n")
cat("HVGs              :", length(hvg), "\n")
cat("Dimensions        : 1:", N_PCS, "\n", sep = "")
cat("k.filter          :", CCA_K_FILTER, "\n")
cat("k.weight          :", CCA_K_WEIGHT, "\n")
cat("dims.to.integrate :", N_PCS, "\n")
cat("preserve.order    : TRUE\n")
cat("sample.tree       : explicit balanced tree\n")
cat("\n")


# ============================================================
# 25. RUN CCA INTEGRATION
# ============================================================

cat("============================================================\n")
cat("RUNNING CCA INTEGRATION\n")
cat("============================================================\n\n")

integration_start <- Sys.time()

integration_error <- NULL

integration_traceback <- NULL


object_integrated <- tryCatch(

  {

    IntegrateLayers(
      object = object,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      assay = "RNA",
      features = hvg,
      normalization.method = "LogNormalize",
      dims = 1:N_PCS,
      k.filter = CCA_K_FILTER,
      k.weight = CCA_K_WEIGHT,
      dims.to.integrate = N_PCS,
      sample.tree = sample_tree,
      preserve.order = TRUE,
      verbose = TRUE
    )

  },

  error = function(e) {

    integration_error <<- conditionMessage(e)

    integration_traceback <<- capture.output(
      traceback()
    )

    NULL
  }
)


integration_end <- Sys.time()

integration_time <- difftime(
  integration_end,
  integration_start,
  units = "mins"
)


cat("\n")
cat("CCA runtime:",
    round(
      as.numeric(integration_time),
      2
    ),
    "minutes\n"
)


# ============================================================
# 26. HANDLE CCA FAILURE
# ============================================================

if (is.null(object_integrated)) {

  cat("\n")
  cat("============================================================\n")
  cat("CCA INTEGRATION FAILED\n")
  cat("============================================================\n\n")

  cat("Error message:\n")
  cat(
    integration_error,
    "\n\n"
  )


  # ----------------------------------------------------------
  # Save diagnostics
  # ----------------------------------------------------------

  diagnostic_file <- paste0(
    "results/tables/",
    "05_CCA_failure_diagnostics.txt"
  )

  diagnostic_lines <- c(
    "============================================================",
    "CCA INTEGRATION FAILURE",
    "============================================================",
    "",
    paste0(
      "Input file: ",
      INPUT_FILE
    ),
    paste0(
      "Cells: ",
      ncol(object)
    ),
    paste0(
      "Genes: ",
      nrow(object)
    ),
    paste0(
      "Samples: ",
      length(sample_ids)
    ),
    "",
    "Samples:",
    paste(
      paste0(
        seq_along(sample_ids),
        " = ",
        sample_ids
      ),
      collapse = "\n"
    ),
    "",
    "RNA layers:",
    paste(
      Layers(object[["RNA"]]),
      collapse = "\n"
    ),
    "",
    "PCA dimensions:",
    as.character(pca_dims),
    "",
    "HVGs:",
    as.character(length(hvg)),
    "",
    "k.filter:",
    as.character(CCA_K_FILTER),
    "",
    "k.weight:",
    as.character(CCA_K_WEIGHT),
    "",
    "Sample tree:",
    paste(
      apply(
        sample_tree,
        1,
        paste,
        collapse = " "
      ),
      collapse = "\n"
    ),
    "",
    "ERROR:",
    integration_error,
    "",
    "TRACEBACK:",
    integration_traceback
  )

  writeLines(
    diagnostic_lines,
    diagnostic_file
  )


  # ----------------------------------------------------------
  # Save preprocessing object for debugging
  # ----------------------------------------------------------

  debug_file <- paste0(
    "results/objects/",
    "05_CCA_preintegration_debug.rds"
  )

  saveRDS(
    object,
    debug_file,
    compress = TRUE
  )


  cat(
    "Diagnostic file saved to:\n  ",
    diagnostic_file,
    "\n",
    sep = ""
  )

  cat(
    "Debug object saved to:\n  ",
    debug_file,
    "\n",
    sep = ""
  )


  stop(
    paste0(
      "\nCCA integration failed.\n\n",
      "Error:\n",
      integration_error,
      "\n\n",
      "The pre-integration object and diagnostic report were saved."
    )
  )
}


# ============================================================
# 27. RENAME OBJECT
# ============================================================

object <- object_integrated

rm(object_integrated)

gc()


# ============================================================
# 28. CHECK INTEGRATED REDUCTION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CHECKING INTEGRATED REDUCTION\n")
cat("============================================================\n\n")

available_reductions <- Reductions(
  object
)

cat("Available reductions:\n")
print(available_reductions)


if (!"integrated.cca" %in% available_reductions) {

  stop(
    paste0(
      "\nERROR: integrated.cca reduction was not created.\n\n",
      "Available reductions:\n",
      paste(
        available_reductions,
        collapse = ", "
      )
    )
  )
}


integrated_embeddings <- Embeddings(
  object,
  reduction = "integrated.cca"
)


integrated_ncells <- nrow(
  integrated_embeddings
)

integrated_ndims <- ncol(
  integrated_embeddings
)


cat(
  "\nIntegrated cells:",
  integrated_ncells,
  "\n"
)

cat(
  "Integrated dimensions:",
  integrated_ndims,
  "\n"
)


if (integrated_ncells != ncol(object)) {

  stop(
    paste0(
      "\nERROR: integrated.cca does not contain the same number ",
      "of cells as the Seurat object.\n\n",
      "Integrated cells: ",
      integrated_ncells,
      "\n",
      "Object cells: ",
      ncol(object)
    )
  )
}


if (integrated_ndims < N_PCS) {

  warning(
    "integrated.cca contains only ",
    integrated_ndims,
    " dimensions."
  )
}


# ============================================================
# 29. FIND NEIGHBORS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINDING NEIGHBORS\n")
cat("============================================================\n\n")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:min(
    N_PCS,
    integrated_ndims
  ),
  verbose = FALSE
)


# ============================================================
# 30. FIND CLUSTERS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINDING CLUSTERS\n")
cat("============================================================\n\n")

object <- FindClusters(
  object,
  resolution = CLUSTER_RESOLUTION,
  verbose = FALSE
)


# ============================================================
# 31. CLUSTER COUNTS
# ============================================================

cluster_counts <- as.data.frame(
  table(
    cluster = Idents(object)
  )
)

colnames(cluster_counts) <- c(
  "cluster",
  "n_cells"
)


write.csv(
  cluster_counts,
  "results/tables/05_cluster_counts.csv",
  row.names = FALSE
)


cat("\nCluster counts:\n\n")
print(cluster_counts)


# ============================================================
# 32. UMAP
# ============================================================

cat("\n")
cat("============================================================\n")
cat("RUNNING UMAP\n")
cat("============================================================\n\n")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:min(
    N_PCS,
    integrated_ndims
  ),
  reduction.name = "umap.cca",
  reduction.key = "UMAPCCA_",
  verbose = FALSE
)


# ============================================================
# 33. UMAP - CLUSTERS
# ============================================================

p_cluster <- DimPlot(
  object,
  reduction = "umap.cca",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "GSE145410 - CCA integrated clusters"
  ) +
  theme_classic()


ggsave(
  "figures/05_integration/05_UMAP_CCA_clusters.pdf",
  p_cluster,
  width = 8,
  height = 6
)


ggsave(
  "figures/05_integration/05_UMAP_CCA_clusters.png",
  p_cluster,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 34. UMAP - SAMPLE
# ============================================================

p_sample <- DimPlot(
  object,
  reduction = "umap.cca",
  group.by = "sample"
) +
  ggtitle(
    "GSE145410 - CCA integrated by sample"
  ) +
  theme_classic()


ggsave(
  "figures/05_integration/05_UMAP_CCA_sample.pdf",
  p_sample,
  width = 9,
  height = 6
)


ggsave(
  "figures/05_integration/05_UMAP_CCA_sample.png",
  p_sample,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 35. UMAP - TREATMENT
# ============================================================

if ("treatment" %in% colnames(object[[]])) {

  p_treatment <- DimPlot(
    object,
    reduction = "umap.cca",
    group.by = "treatment"
  ) +
    ggtitle(
      "GSE145410 - CCA integrated by treatment"
    ) +
    theme_classic()


  ggsave(
    "figures/05_integration/05_UMAP_CCA_treatment.pdf",
    p_treatment,
    width = 9,
    height = 6
  )


  ggsave(
    "figures/05_integration/05_UMAP_CCA_treatment.png",
    p_treatment,
    width = 9,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 36. UMAP - REPLICATE
# ============================================================

if ("replicate" %in% colnames(object[[]])) {

  p_replicate <- DimPlot(
    object,
    reduction = "umap.cca",
    group.by = "replicate"
  ) +
    ggtitle(
      "GSE145410 - CCA integrated by replicate"
    ) +
    theme_classic()


  ggsave(
    "figures/05_integration/05_UMAP_CCA_replicate.pdf",
    p_replicate,
    width = 8,
    height = 6
  )


  ggsave(
    "figures/05_integration/05_UMAP_CCA_replicate.png",
    p_replicate,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 37. INTEGRATION SUMMARY
# ============================================================

integration_summary <- data.frame(

  parameter = c(
    "input_cells",
    "input_genes",
    "number_of_samples",
    "number_of_HVGs",
    "PCA_dimensions",
    "integrated_CCA_dimensions",
    "CCA_k_filter",
    "CCA_k_weight",
    "cluster_resolution",
    "integration_runtime_minutes"
  ),

  value = c(
    ncol(object),
    nrow(object),
    length(sample_ids),
    length(hvg),
    pca_dims,
    integrated_ndims,
    CCA_K_FILTER,
    CCA_K_WEIGHT,
    CLUSTER_RESOLUTION,
    round(
      as.numeric(integration_time),
      2
    )
  )
)


write.csv(
  integration_summary,
  "results/tables/05_integration_summary.csv",
  row.names = FALSE
)


# ============================================================
# 38. SAVE INTEGRATED OBJECT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SAVING INTEGRATED OBJECT\n")
cat("============================================================\n\n")

saveRDS(
  object,
  OUTPUT_FILE,
  compress = TRUE
)


cat(
  "Saved:\n  ",
  OUTPUT_FILE,
  "\n",
  sep = ""
)


# ============================================================
# 39. FINAL DIAGNOSTICS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL OBJECT DIAGNOSTICS\n")
cat("============================================================\n\n")

cat("Cells:", ncol(object), "\n")
cat("Genes:", nrow(object), "\n")

cat("\nSamples:\n")
print(
  table(object$sample)
)

cat("\nRNA layers:\n")
print(
  Layers(object[["RNA"]])
)

cat("\nReductions:\n")
print(
  Reductions(object)
)

cat("\nClusters:\n")
print(
  table(Idents(object))
)


# ============================================================
# 40. SESSION INFORMATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SESSION INFORMATION\n")
cat("============================================================\n\n")

print(
  sessionInfo()
)


# ============================================================
# 41. COMPLETION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("05_Seurat5_integration.R V2 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Input:\n")
cat("  ", INPUT_FILE, "\n\n", sep = "")

cat("Output:\n")
cat("  ", OUTPUT_FILE, "\n\n", sep = "")

cat("Integration:\n")
cat("  Method            : CCAIntegration\n")
cat("  HVGs              :", length(hvg), "\n")
cat("  PCA dimensions    :", N_PCS, "\n")
cat("  k.filter          :", CCA_K_FILTER, "\n")
cat("  k.weight          :", CCA_K_WEIGHT, "\n")
cat("  dims.to.integrate :", N_PCS, "\n")
cat("  preserve.order    : TRUE\n")
cat("  sample.tree       : explicit balanced tree\n")
cat("  cluster resolution:", CLUSTER_RESOLUTION, "\n")

cat("\n")
cat("No scAnnoX.\n")
cat("No manual Harmony integration.\n")
cat("\n")
