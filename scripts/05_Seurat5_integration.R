# ============================================================
# 05_Seurat5_integration.R
# GSE145410 AML scRNA-seq
#
# Seurat v5 CCA integration
#
# Input:
#   results/objects/03_normalized.rds
#
# Output:
#   results/objects/05_integrated.rds
#   results/tables/05_cluster_counts.csv
#   figures/05_integration/
#
# Important:
#   - Seurat v5 layer-based integration
#   - Integration by sample
#   - CCAIntegration
#   - NO scAnnoX
#   - NO manual Harmony integration
#   - scPred/Harmony is handled later in the annotation workflow
# ============================================================


# ============================================================
# 0. Setup
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================================\n")
cat("05 - Seurat v5 CCA integration\n")
cat("============================================================\n\n")


# ============================================================
# 1. Load packages
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
# 2. Create output directories
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
# 3. Input
# ============================================================

input_file <- "results/objects/03_normalized.rds"

if (!file.exists(input_file)) {
  stop(
    paste0(
      "\nERROR: Input file not found:\n",
      input_file,
      "\n\nRun 03_normalization_HVG.R first.\n"
    )
  )
}

cat("Loading:", input_file, "\n")

object <- readRDS(input_file)


# ============================================================
# 4. Basic validation
# ============================================================

if (!inherits(object, "Seurat")) {
  stop("ERROR: input object is not a Seurat object.")
}

if (!"RNA" %in% Assays(object)) {
  stop("ERROR: RNA assay not found.")
}

if (!"sample" %in% colnames(object[[]])) {
  stop(
    "ERROR: metadata column 'sample' is missing. ",
    "The integration must be performed by sample."
  )
}

cat("\nObject successfully loaded.\n")
cat("Cells:", ncol(object), "\n")
cat("Genes:", nrow(object), "\n")

cat("\nMetadata columns:\n")
print(colnames(object[[]]))


# ============================================================
# 5. Define samples
# ============================================================

sample_ids <- sort(unique(as.character(object$sample)))

cat("\n============================================================\n")
cat("Samples detected\n")
cat("============================================================\n")

print(sample_ids)

cat("\nNumber of samples:", length(sample_ids), "\n")


# Expected GSE145410 samples
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

missing_samples <- setdiff(expected_samples, sample_ids)
extra_samples <- setdiff(sample_ids, expected_samples)

if (length(missing_samples) > 0) {
  warning(
    "Expected samples missing: ",
    paste(missing_samples, collapse = ", ")
  )
}

if (length(extra_samples) > 0) {
  warning(
    "Additional sample labels detected: ",
    paste(extra_samples, collapse = ", ")
  )
}


# ============================================================
# 6. Cells per sample
# ============================================================

cells_per_sample <- table(object$sample)

cat("\n============================================================\n")
cat("Cells per sample\n")
cat("============================================================\n")

print(cells_per_sample)

write.csv(
  as.data.frame(cells_per_sample),
  "results/tables/05_cells_per_sample_before_integration.csv",
  row.names = FALSE
)


# ============================================================
# 7. Check current RNA layers
# ============================================================

cat("\n============================================================\n")
cat("RNA layers before preparation\n")
cat("============================================================\n")

rna_layers_before <- Layers(object[["RNA"]])

print(rna_layers_before)


# ============================================================
# 8. Rejoin RNA layers if necessary
#
# We start from 03_normalized.rds.
#
# Integration requires the assay to be split by sample.
# If the object already contains sample-specific layers,
# first join them to avoid nested/repeated split layers.
# ============================================================

has_sample_layers <- any(
  grepl(
    paste0(
      "\\.",
      paste(sample_ids, collapse = "|"),
      "$"
    ),
    rna_layers_before
  )
)

if (has_sample_layers) {

  cat("\nSample-specific RNA layers already detected.\n")
  cat("Joining RNA layers before the new split...\n\n")

  object[["RNA"]] <- JoinLayers(
    object[["RNA"]]
  )

} else {

  cat("\nRNA assay is already in joined form.\n")
}


# ============================================================
# 9. Check layers after JoinLayers
# ============================================================

cat("\n============================================================\n")
cat("RNA layers after JoinLayers check\n")
cat("============================================================\n")

print(
  Layers(object[["RNA"]])
)


# ============================================================
# 10. Split RNA assay by sample
#
# Seurat v5 integration works with sample-specific layers:
#
# counts.SAMPLE
# data.SAMPLE
#
# We keep all cells in a single Seurat object.
# ============================================================

cat("\n============================================================\n")
cat("Splitting RNA assay by sample\n")
cat("============================================================\n\n")

object[["RNA"]] <- split(
  object[["RNA"]],
  f = object$sample
)


# ============================================================
# 11. Inspect split layers
# ============================================================

cat("\nRNA layers after split:\n\n")

rna_layers_split <- Layers(
  object[["RNA"]]
)

print(rna_layers_split)


# ============================================================
# 12. Validate sample-specific counts/data layers
# ============================================================

counts_layers <- grep(
  "^counts\\.",
  rna_layers_split,
  value = TRUE
)

data_layers <- grep(
  "^data\\.",
  rna_layers_split,
  value = TRUE
)

cat("\nCounts layers detected:", length(counts_layers), "\n")
print(counts_layers)

cat("\nData layers detected:", length(data_layers), "\n")
print(data_layers)


if (length(counts_layers) != length(sample_ids)) {

  stop(
    paste0(
      "\nERROR: Expected ",
      length(sample_ids),
      " sample-specific counts layers but found ",
      length(counts_layers),
      ".\n\n",
      "Current RNA layers:\n",
      paste(rna_layers_split, collapse = "\n")
    )
  )
}

if (length(data_layers) != length(sample_ids)) {

  stop(
    paste0(
      "\nERROR: Expected ",
      length(sample_ids),
      " sample-specific data layers but found ",
      length(data_layers),
      ".\n\n",
      "Current RNA layers:\n",
      paste(rna_layers_split, collapse = "\n")
    )
  )
}


# ============================================================
# 13. Normalization
#
# The input was already normalized in script 03.
#
# We nevertheless call NormalizeData here because Seurat v5
# performs normalization independently on the split layers.
# This follows the Seurat v5 integration workflow.
# ============================================================

cat("\n============================================================\n")
cat("Normalizing split RNA layers\n")
cat("============================================================\n\n")

object <- NormalizeData(
  object,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)


# ============================================================
# 14. Find variable features
#
# 2,000 HVGs, as used in the rest of the project.
#
# With split layers, Seurat determines the variable features
# for the individual layers and builds a consensus feature set.
# ============================================================

cat("\n============================================================\n")
cat("Finding variable features\n")
cat("============================================================\n\n")

object <- FindVariableFeatures(
  object,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

hvg <- VariableFeatures(object)

cat("Number of variable features:", length(hvg), "\n")

if (length(hvg) == 0) {
  stop("ERROR: No variable features were identified.")
}

if (length(hvg) < 1000) {
  warning(
    "Less than 1,000 variable features were identified: ",
    length(hvg)
  )
}

cat("\nFirst variable features:\n")
print(head(hvg, 20))


# ============================================================
# 15. Save HVG table
# ============================================================

hvg_table <- data.frame(
  feature = hvg,
  rank = seq_along(hvg)
)

write.csv(
  hvg_table,
  "results/tables/05_variable_features.csv",
  row.names = FALSE
)


# ============================================================
# 16. Scale data
#
# IMPORTANT:
# Do NOT require eight scale.data.sample layers here.
#
# Seurat v5's IntegrateLayers() accepts:
#
#   scale.layer = "scale.data"
#
# and internally handles the split assay layers.
#
# The previous script incorrectly expected:
#
#   scale.data.DMSO_A
#   scale.data.DMSO_B
#   ...
#
# which caused the CI failure:
#
#   Expected 8 scale.data layers, found 0
#
# We therefore only verify that Seurat produced a scale.data
# layer usable by the integration workflow.
# ============================================================

cat("\n============================================================\n")
cat("Scaling variable features\n")
cat("============================================================\n\n")

object <- ScaleData(
  object,
  assay = "RNA",
  features = hvg,
  verbose = FALSE
)


# ============================================================
# 17. Inspect scale.data layers
# ============================================================

cat("\n============================================================\n")
cat("Checking scale.data after ScaleData\n")
cat("============================================================\n\n")

rna_layers_scaled <- Layers(
  object[["RNA"]]
)

print(rna_layers_scaled)

scale_layers <- grep(
  "^scale\\.data",
  rna_layers_scaled,
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
      "\nERROR: ScaleData() did not create a usable scale.data layer.\n\n",
      "Current RNA layers:\n",
      paste(rna_layers_scaled, collapse = "\n"),
      "\n\n",
      "The integration cannot continue without scaled data."
    )
  )
}


# ============================================================
# 18. PCA
#
# PCA is computed on the same 2,000 HVGs.
#
# 20 PCs are used for the CCA integration to reduce memory
# consumption on GitHub Actions while retaining sufficient
# dimensional information for this dataset.
# ============================================================

cat("\n============================================================\n")
cat("Running PCA\n")
cat("============================================================\n\n")

object <- RunPCA(
  object,
  assay = "RNA",
  features = hvg,
  npcs = 20,
  verbose = FALSE
)


# ============================================================
# 19. Validate PCA
# ============================================================

if (!"pca" %in% Reductions(object)) {
  stop("ERROR: PCA reduction was not created.")
}

pca_dims <- ncol(
  Embeddings(
    object,
    reduction = "pca"
  )
)

cat(
  "PCA dimensions available:",
  pca_dims,
  "\n"
)

if (pca_dims < 20) {
  warning(
    "Fewer than 20 PCA dimensions were generated: ",
    pca_dims
  )
}


# ============================================================
# 20. Save PCA plot
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
    ndims = min(20, pca_dims)
  )
)

dev.off()


# ============================================================
# 21. CCA integration
#
# Official Seurat v5 workflow:
#
# IntegrateLayers(
#   object = object,
#   method = CCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "integrated.cca"
# )
#
# We explicitly provide:
#   - assay = RNA
#   - features = HVGs
#   - layers = normalized sample-specific layers
#   - scale.layer = "scale.data"
#   - dims = 1:20
#   - k.filter = 50
#   - k.weight = 50
#   - dims.to.integrate = 20
#
# No Harmony integration is performed here.
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Running Seurat v5 CCA integration\n")
cat("============================================================\n\n")

cat("Integration parameters:\n")
cat("  Method              : CCAIntegration\n")
cat("  Assay               : RNA\n")
cat("  Features            :", length(hvg), "\n")
cat("  Dimensions          : 1:20\n")
cat("  k.filter            : 50\n")
cat("  k.weight            : 50\n")
cat("  dims.to.integrate   : 20\n")
cat("  scale.layer         : scale.data\n")
cat("\n")


# ------------------------------------------------------------
# Integration wrapped in tryCatch to provide useful diagnostics
# ------------------------------------------------------------

integration_error <- NULL

object <- tryCatch(

  {

    IntegrateLayers(
      object = object,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      assay = "RNA",
      features = hvg,
      layers = data_layers,
      scale.layer = "scale.data",
      normalization.method = "LogNormalize",
      dims = 1:20,
      k.filter = 50,
      k.weight = 50,
      dims.to.integrate = 20,
      verbose = TRUE
    )

  },

  error = function(e) {

    integration_error <<- conditionMessage(e)

    cat("\n")
    cat("============================================================\n")
    cat("CCA INTEGRATION ERROR\n")
    cat("============================================================\n")
    cat("\n")
    cat(conditionMessage(e))
    cat("\n\n")

    NULL
  }
)


# ============================================================
# 22. Stop cleanly if integration failed
# ============================================================

if (is.null(object)) {

  cat("\n============================================================\n")
  cat("Integration diagnostics\n")
  cat("============================================================\n\n")

  cat("RNA layers:\n")
  print(Layers(object = readRDS(input_file)[["RNA"]]))

  cat("\nSamples:\n")
  print(table(object$sample))

  cat("\nOriginal integration error:\n")
  cat(integration_error, "\n")

  stop(
    paste0(
      "\nCCA integration failed.\n\n",
      "Original error:\n",
      integration_error,
      "\n\n",
      "See the GitHub Actions log for the detailed Seurat traceback."
    )
  )
}


# ============================================================
# 23. Validate integrated reduction
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Validating integrated reduction\n")
cat("============================================================\n\n")

available_reductions <- Reductions(object)

cat("Available reductions:\n")
print(available_reductions)


if (!"integrated.cca" %in% available_reductions) {

  stop(
    paste0(
      "\nERROR: integrated.cca reduction was not created.\n\n",
      "Available reductions:\n",
      paste(available_reductions, collapse = ", ")
    )
  )
}


integrated_embeddings <- Embeddings(
  object,
  reduction = "integrated.cca"
)

cat(
  "\nIntegrated cells:",
  nrow(integrated_embeddings),
  "\n"
)

cat(
  "Integrated dimensions:",
  ncol(integrated_embeddings),
  "\n"
)


if (nrow(integrated_embeddings) != ncol(object)) {

  stop(
    paste0(
      "\nERROR: Number of cells in integrated reduction does not ",
      "match the Seurat object.\n",
      "Integrated cells: ",
      nrow(integrated_embeddings),
      "\n",
      "Object cells: ",
      ncol(object)
    )
  )
}


# ============================================================
# 24. Neighbors using integrated CCA space
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Computing neighbors from integrated CCA space\n")
cat("============================================================\n\n")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  verbose = FALSE
)


# ============================================================
# 25. Clustering
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Clustering\n")
cat("============================================================\n\n")

object <- FindClusters(
  object,
  resolution = 0.4,
  verbose = FALSE
)


# ============================================================
# 26. Cluster statistics
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

cat("\nCluster counts:\n")
print(cluster_counts)


# ============================================================
# 27. UMAP from integrated CCA
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Running UMAP on integrated CCA space\n")
cat("============================================================\n\n")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  reduction.name = "umap.cca",
  reduction.key = "UMAPCCA_",
  verbose = FALSE
)


# ============================================================
# 28. UMAP - clusters
# ============================================================

p_cluster <- DimPlot(
  object,
  reduction = "umap.cca",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("GSE145410 - CCA integrated clusters") +
  theme_classic()

ggsave(
  filename = "figures/05_integration/05_UMAP_CCA_clusters.pdf",
  plot = p_cluster,
  width = 8,
  height = 6
)

ggsave(
  filename = "figures/05_integration/05_UMAP_CCA_clusters.png",
  plot = p_cluster,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 29. UMAP - sample
# ============================================================

p_sample <- DimPlot(
  object,
  reduction = "umap.cca",
  group.by = "sample"
) +
  ggtitle("GSE145410 - CCA integrated by sample") +
  theme_classic()

ggsave(
  filename = "figures/05_integration/05_UMAP_CCA_sample.pdf",
  plot = p_sample,
  width = 9,
  height = 6
)

ggsave(
  filename = "figures/05_integration/05_UMAP_CCA_sample.png",
  plot = p_sample,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 30. UMAP - treatment
# ============================================================

if ("treatment" %in% colnames(object[[]])) {

  p_treatment <- DimPlot(
    object,
    reduction = "umap.cca",
    group.by = "treatment"
  ) +
    ggtitle("GSE145410 - CCA integrated by treatment") +
    theme_classic()

  ggsave(
    filename = "figures/05_integration/05_UMAP_CCA_treatment.pdf",
    plot = p_treatment,
    width = 9,
    height = 6
  )

  ggsave(
    filename = "figures/05_integration/05_UMAP_CCA_treatment.png",
    plot = p_treatment,
    width = 9,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 31. UMAP - replicate
# ============================================================

if ("replicate" %in% colnames(object[[]])) {

  p_replicate <- DimPlot(
    object,
    reduction = "umap.cca",
    group.by = "replicate"
  ) +
    ggtitle("GSE145410 - CCA integrated by replicate") +
    theme_classic()

  ggsave(
    filename = "figures/05_integration/05_UMAP_CCA_replicate.pdf",
    plot = p_replicate,
    width = 8,
    height = 6
  )

  ggsave(
    filename = "figures/05_integration/05_UMAP_CCA_replicate.png",
    plot = p_replicate,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 32. Integration summary
# ============================================================

integration_summary <- data.frame(
  parameter = c(
    "input_cells",
    "input_genes",
    "number_of_samples",
    "number_of_HVGs",
    "PCA_dimensions",
    "CCA_dimensions",
    "CCA_k_filter",
    "CCA_k_weight",
    "cluster_resolution"
  ),
  value = c(
    ncol(object),
    nrow(object),
    length(sample_ids),
    length(hvg),
    pca_dims,
    ncol(integrated_embeddings),
    50,
    50,
    0.4
  )
)

write.csv(
  integration_summary,
  "results/tables/05_integration_summary.csv",
  row.names = FALSE
)


# ============================================================
# 33. Save integrated object
#
# IMPORTANT:
# We keep the RNA layers split here.
#
# JoinLayers() will be performed later before differential
# expression, after annotation/integration analysis.
# ============================================================

output_file <- "results/objects/05_integrated.rds"

cat("\n")
cat("============================================================\n")
cat("Saving integrated object\n")
cat("============================================================\n\n")

saveRDS(
  object,
  output_file,
  compress = TRUE
)

cat("Saved:", output_file, "\n")


# ============================================================
# 34. Final object diagnostics
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL OBJECT\n")
cat("============================================================\n\n")

cat("Cells:", ncol(object), "\n")
cat("Genes:", nrow(object), "\n")

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
# 35. Session information
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Session information\n")
cat("============================================================\n\n")

print(
  sessionInfo()
)


# ============================================================
# 36. Completion message
# ============================================================

cat("\n")
cat("============================================================\n")
cat("05_Seurat5_integration.R COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
cat("\n")

cat("Output:\n")
cat("  ", output_file, "\n", sep = "")

cat("\nFigures:\n")
cat("  figures/05_integration/\n")

cat("\nTables:\n")
cat("  results/tables/05_*.csv\n")

cat("\nIntegration:\n")
cat("  Method       : CCAIntegration\n")
cat("  PCs          : 20\n")
cat("  HVGs         : 2000\n")
cat("  k.filter     : 50\n")
cat("  k.weight     : 50\n")
cat("  Resolution   : 0.4\n")

cat("\nNo scAnnoX used.\n")
cat("No manual Harmony integration used.\n")
cat("\n")
