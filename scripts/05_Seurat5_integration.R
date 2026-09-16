# ============================================================
# 05_Seurat5_integration.R
# GSE145410 AML scRNA-seq
#
# Seurat v5 CCA integration - V2.1
#
# Input:
#   results/objects/03_normalized.rds
#
# Output:
#   results/objects/05_integrated.rds
#   results/tables/05_*.csv
#   figures/05_integration/*
#
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================================\n")
cat("05 - Seurat v5 CCA integration V2.1\n")
cat("============================================================\n\n")


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
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
# 3. DIRECTORIES
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
# 4. INPUT
# ============================================================

if (!file.exists(INPUT_FILE)) {
  stop(
    paste0(
      "Input file not found: ",
      INPUT_FILE
    )
  )
}

cat("Loading object:\n")
cat(INPUT_FILE, "\n\n")

object <- readRDS(INPUT_FILE)


# ============================================================
# 5. VALIDATION
# ============================================================

if (!inherits(object, "Seurat")) {
  stop("Input is not a Seurat object.")
}

if (!"RNA" %in% Assays(object)) {
  stop("RNA assay not found.")
}

if (!"sample" %in% colnames(object[[]])) {
  stop("Metadata column 'sample' not found.")
}


cat("Cells:", ncol(object), "\n")
cat("Genes:", nrow(object), "\n\n")


# ============================================================
# 6. SAMPLE DEFINITIONS
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

sample_ids <- sort(
  unique(
    as.character(object$sample)
  )
)

cat("Samples detected:\n")
print(sample_ids)
cat("\n")


# ============================================================
# 7. CHECK SAMPLE SET
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
      "Missing expected samples: ",
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
      "Unexpected samples detected: ",
      paste(
        extra_samples,
        collapse = ", "
      )
    )
  )
}

sample_ids <- expected_samples


# ============================================================
# 8. CELLS PER SAMPLE
# ============================================================

cells_per_sample <- table(
  factor(
    object$sample,
    levels = sample_ids
  )
)

cat("Cells per sample:\n")
print(cells_per_sample)
cat("\n")

write.csv(
  as.data.frame(cells_per_sample),
  "results/tables/05_cells_per_sample.csv",
  row.names = FALSE
)


# ============================================================
# 9. INITIAL RNA LAYERS
# ============================================================

cat("============================================================\n")
cat("INITIAL RNA LAYERS\n")
cat("============================================================\n\n")

initial_layers <- Layers(
  object[["RNA"]]
)

print(initial_layers)
cat("\n")


# ============================================================
# 10. JOIN EXISTING SAMPLE LAYERS
# ============================================================

sample_pattern <- paste(
  sample_ids,
  collapse = "|"
)

sample_layer_pattern <- paste0(
  "\\.(",
  sample_pattern,
  ")$"
)

has_sample_layers <- any(
  grepl(
    sample_layer_pattern,
    initial_layers
  )
)

if (has_sample_layers) {

  cat("Sample-specific layers detected.\n")
  cat("Joining RNA layers before splitting again.\n\n")

  object[["RNA"]] <- JoinLayers(
    object[["RNA"]]
  )

} else {

  cat("No sample-specific layers detected.\n")
  cat("RNA assay is already suitable for splitting.\n\n")
}


# ============================================================
# 11. SPLIT RNA BY SAMPLE
# ============================================================

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
# 12. CHECK SPLIT LAYERS
# ============================================================

split_layers <- Layers(
  object[["RNA"]]
)

cat("RNA layers after split:\n\n")
print(split_layers)
cat("\n")


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

cat("Counts layers:", length(counts_layers), "\n")
print(counts_layers)

cat("\nData layers:", length(data_layers), "\n")
print(data_layers)

cat("\n")


if (length(counts_layers) != length(sample_ids)) {
  stop(
    paste0(
      "Expected ",
      length(sample_ids),
      " counts layers, found ",
      length(counts_layers)
    )
  )
}

if (length(data_layers) != length(sample_ids)) {
  stop(
    paste0(
      "Expected ",
      length(sample_ids),
      " data layers, found ",
      length(data_layers)
    )
  )


cat("\nSplit validation successful.\n\n")


# ============================================================
# 13. NORMALIZATION
# ============================================================

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
# 14. VARIABLE FEATURES
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
  "\n\n"
)

if (length(hvg) == 0) {
  stop("No variable features detected.")
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
# 15. SCALE DATA
# ============================================================

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
# 16. CHECK SCALE DATA
# ============================================================

scale_layers <- grep(
  "^scale\\.data",
  Layers(object[["RNA"]]),
  value = TRUE
)

cat("\nScale layer(s):\n")
print(scale_layers)
cat("\n")

if (length(scale_layers) == 0) {
  stop(
    "ScaleData did not produce a scale.data layer."
  )
}


# ============================================================
# 17. PCA
# ============================================================

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

if (!"pca" %in% Reductions(object)) {
  stop("PCA reduction was not created.")
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
  "\n\n"
)


# ============================================================
# 18. ELBOW PLOT
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
# 19. EXPLICIT BALANCED SAMPLE TREE
# ============================================================
#
# Dataset numbering:
#
# 1 = DMSO_A
# 2 = DMSO_B
# 3 = INCB059872_A
# 4 = INCB059872_B
# 5 = AZA_A
# 6 = AZA_B
# 7 = INCB059872_AZA_A
# 8 = INCB059872_AZA_B
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

cat("============================================================\n")
cat("SAMPLE TREE\n")
cat("============================================================\n\n")

print(sample_tree)

cat("\nDataset mapping:\n")

for (i in seq_along(sample_ids)) {
  cat(
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
# 20. CCA PARAMETERS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CCA PARAMETERS\n")
cat("============================================================\n\n")

cat("HVGs              :", N_HVG, "\n")
cat("PCA dimensions    :", N_PCS, "\n")
cat("k.filter          :", CCA_K_FILTER, "\n")
cat("k.weight          :", CCA_K_WEIGHT, "\n")
cat("dims.to.integrate :", N_PCS, "\n")
cat("preserve.order    : TRUE\n")
cat("\n")


# ============================================================
# 21. CCA INTEGRATION
# ============================================================

cat("============================================================\n")
cat("RUNNING CCA INTEGRATION\n")
cat("============================================================\n\n")

integration_start <- Sys.time()

integration_error <- NULL

object_integrated <- tryCatch(

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
  ),

  error = function(e) {

    integration_error <<- conditionMessage(e)

    return(NULL)
  }
)

integration_end <- Sys.time()

integration_minutes <- as.numeric(
  difftime(
    integration_end,
    integration_start,
    units = "mins"
  )
)

cat("\n")
cat(
  "CCA runtime:",
  round(
    integration_minutes,
    2
  ),
  "minutes\n"
)


# ============================================================
# 22. HANDLE INTEGRATION ERROR
# ============================================================

if (is.null(object_integrated)) {

  cat("\n")
  cat("============================================================\n")
  cat("CCA INTEGRATION FAILED\n")
  cat("============================================================\n\n")

  cat("Error:\n")
  cat(
    integration_error,
    "\n\n"
  )


  diagnostic_file <- paste0(
    "results/tables/",
    "05_CCA_failure_diagnostics.txt"
  )

  diagnostic_text <- c(
    "CCA INTEGRATION FAILURE",
    "",
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
    paste0(
      "HVGs: ",
      length(hvg)
    ),
    paste0(
      "PCA dimensions: ",
      pca_dims
    ),
    paste0(
      "k.filter: ",
      CCA_K_FILTER
    ),
    paste0(
      "k.weight: ",
      CCA_K_WEIGHT
    ),
    "",
    "Samples:",
    paste(
      sample_ids,
      collapse = "\n"
    ),
    "",
    "RNA layers:",
    paste(
      Layers(object[["RNA"]]),
      collapse = "\n"
    ),
    "",
    "Error:",
    integration_error
  )

  writeLines(
    diagnostic_text,
    diagnostic_file
  )


  saveRDS(
    object,
    "results/objects/05_CCA_preintegration_debug.rds",
    compress = TRUE
  )


  stop(
    paste0(
      "\nCCA integration failed: ",
      integration_error
    )
  )
}


# ============================================================
# 23. FINAL INTEGRATED OBJECT
# ============================================================

object <- object_integrated

rm(object_integrated)

gc()


# ============================================================
# 24. CHECK INTEGRATED REDUCTION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CHECKING INTEGRATED REDUCTION\n")
cat("============================================================\n\n")

print(
  Reductions(object)
)

if (!"integrated.cca" %in% Reductions(object)) {
  stop(
    "integrated.cca reduction was not created."
  )
}

integrated_embeddings <- Embeddings(
  object,
  reduction = "integrated.cca"
)

integrated_cells <- nrow(
  integrated_embeddings
)

integrated_dims <- ncol(
  integrated_embeddings
)

cat("\nIntegrated cells:", integrated_cells, "\n")
cat("Integrated dimensions:", integrated_dims, "\n")

if (integrated_cells != ncol(object)) {
  stop(
    "Number of integrated cells does not match Seurat object."
  )
}


# ============================================================
# 25. NEIGHBORS
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
    integrated_dims
  ),
  verbose = FALSE
)


# ============================================================
# 26. CLUSTERS
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

print(cluster_counts)


# ============================================================
# 27. UMAP
# ============================================================

cat("\n")
cat("============================================================\n")
cat("UMAP\n")
cat("============================================================\n\n")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:min(
    N_PCS,
    integrated_dims
  ),
  reduction.name = "umap.cca",
  reduction.key = "UMAPCCA_",
  verbose = FALSE
)


# ============================================================
# 28. UMAP - CLUSTERS
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
# 29. UMAP - SAMPLE
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
# 30. UMAP - TREATMENT
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
# 31. UMAP - REPLICATE
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
# 32. SUMMARY
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
    integrated_dims,
    CCA_K_FILTER,
    CCA_K_WEIGHT,
    CLUSTER_RESOLUTION,
    round(
      integration_minutes,
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
# 33. SAVE
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
  "Saved:",
  OUTPUT_FILE,
  "\n"
)


# ============================================================
# 34. FINAL DIAGNOSTICS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL DIAGNOSTICS\n")
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
# 35. SESSION INFO
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SESSION INFO\n")
cat("============================================================\n\n")

print(
  sessionInfo()
)


# ============================================================
# 36. SUCCESS
# ============================================================

cat("\n")
cat("============================================================\n")
cat("05_Seurat5_integration.R V2.1 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Output:\n")
cat(
  OUTPUT_FILE,
  "\n\n"
)

cat("CCA integration:\n")
cat("  HVGs       :", N_HVG, "\n")
cat("  PCs        :", N_PCS, "\n")
cat("  k.filter   :", CCA_K_FILTER, "\n")
cat("  k.weight   :", CCA_K_WEIGHT, "\n")
cat("  resolution :", CLUSTER_RESOLUTION, "\n")

cat("\nNo scAnnoX.\n")
cat("No manual Harmony integration.\n\n")
