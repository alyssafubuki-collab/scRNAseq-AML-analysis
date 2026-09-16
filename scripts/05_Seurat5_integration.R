# ============================================================
# 05_Seurat5_integration.R
# Seurat v5 integration using CCA
#
# Input:
#   results/03_normalized.rds
#
# Output:
#   results/05_integrated.rds
#   results/05_integration_umap.png
#   results/05_integration_elbow.png
#
# Annotation is NOT performed here.
# scAnnoX is NOT used.
# ============================================================

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  future.globals.maxSize = 8 * 1024^3
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

cat("\n")
cat("============================================================\n")
cat("05 - SEURAT V5 CCA INTEGRATION\n")
cat("============================================================\n\n")


# ============================================================
# 1. Paths
# ============================================================

project_dir <- getwd()

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

umap_file <- file.path(
  project_dir,
  "figures",
  "integration",
  "05_integration_umap.png"
)

elbow_file <- file.path(
  project_dir,
  "figures",
  "integration",
  "05_integration_elbow.png"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ",
      input_file
    )
  )
}


# ============================================================
# 2. Load object
# ============================================================

cat("Loading:", input_file, "\n")

object <- readRDS(input_file)

cat("Object loaded successfully.\n\n")

cat("Cells:", ncol(object), "\n")
cat("Genes:", nrow(object), "\n\n")


# ============================================================
# 3. Basic validation
# ============================================================

if (!"RNA" %in% Assays(object)) {
  stop("RNA assay not found in Seurat object.")
}

DefaultAssay(object) <- "RNA"

cat("Default assay:", DefaultAssay(object), "\n\n")


# ============================================================
# 4. Inspect initial RNA layers
# ============================================================

cat("============================================================\n")
cat("INITIAL RNA LAYERS\n")
cat("============================================================\n")

initial_layers <- Layers(object[["RNA"]])

print(initial_layers)

cat("\n")


# ============================================================
# 5. Detect whether object is already split
# ============================================================

counts_layers <- grep(
  "^counts\\.",
  initial_layers,
  value = TRUE
)

data_layers <- grep(
  "^data\\.",
  initial_layers,
  value = TRUE
)

sample_ids <- character(0)

if (length(counts_layers) > 0) {

  sample_ids <- sub(
    "^counts\\.",
    "",
    counts_layers
  )

}

if (length(counts_layers) == 0) {

  cat("No sample-specific counts layers detected.\n")
  cat("RNA layers will be split by sample.\n\n")

} else {

  cat("Sample-specific counts layers already detected.\n")

  cat("Counts layers:", length(counts_layers), "\n")
  print(counts_layers)

  cat("\n")

  cat("Data layers:", length(data_layers), "\n")
  print(data_layers)

  cat("\n")

}


# ============================================================
# 6. Split RNA assay if necessary
# ============================================================

if (length(counts_layers) == 0) {

  if (!"sample" %in% colnames(object@meta.data)) {
    stop(
      "Metadata column 'sample' is required to split the RNA assay."
    )
  }

  sample_ids <- unique(
    as.character(object$sample)
  )

  sample_ids <- sample_ids[
    !is.na(sample_ids) &
    sample_ids != ""
  ]

  if (length(sample_ids) != 8) {

    stop(
      paste0(
        "Expected 8 samples, found ",
        length(sample_ids),
        "."
      )
    )

  }

  cat("Samples detected:\n")
  print(sample_ids)

  cat("\n")

  cat("Splitting RNA assay by sample...\n\n")

  object[["RNA"]] <- split(
    object[["RNA"]],
    f = object$sample
  )

  cat("RNA layers after split:\n")

  split_layers <- Layers(object[["RNA"]])

  print(split_layers)

  cat("\n")

}


# ============================================================
# 7. Validate split layers
# ============================================================

rna_layers <- Layers(object[["RNA"]])

counts_layers <- grep(
  "^counts\\.",
  rna_layers,
  value = TRUE
)

data_layers <- grep(
  "^data\\.",
  rna_layers,
  value = TRUE
)

cat("Counts layers:", length(counts_layers), "\n")
print(counts_layers)

cat("\n")

cat("Data layers:", length(data_layers), "\n")
print(data_layers)

cat("\n")


if (length(counts_layers) != length(sample_ids)) {

  stop(
    paste0(
      "Expected ",
      length(sample_ids),
      " counts layers, found ",
      length(counts_layers),
      "."
    )
  )

}


if (length(data_layers) != length(sample_ids)) {

  stop(
    paste0(
      "Expected ",
      length(sample_ids),
      " data layers, found ",
      length(data_layers),
      "."
    )
  )

}


cat("Split validation successful.\n\n")


# ============================================================
# 8. Normalize sample-specific layers
# ============================================================

cat("============================================================\n")
cat("NORMALIZATION\n")
cat("============================================================\n\n")

cat("Running NormalizeData on sample-specific RNA layers...\n")

object <- NormalizeData(
  object,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

cat("\nNormalization completed.\n\n")


# ============================================================
# 9. Select integration features
# ============================================================

cat("============================================================\n")
cat("INTEGRATION FEATURES\n")
cat("============================================================\n\n")

integration_features <- SelectIntegrationFeatures(
  object = object,
  nfeatures = 2000
)

cat(
  "Number of integration features:",
  length(integration_features),
  "\n\n"
)


# ============================================================
# 10. Scale integration features
# ============================================================

cat("============================================================\n")
cat("SCALING\n")
cat("============================================================\n\n")

cat("Scaling integration features...\n")

object <- ScaleData(
  object,
  assay = "RNA",
  features = integration_features,
  verbose = TRUE
)

cat("\nScaling completed.\n\n")


# ============================================================
# 11. PCA
# ============================================================

cat("============================================================\n")
cat("PCA\n")
cat("============================================================\n\n")

cat("Running PCA...\n")

object <- RunPCA(
  object,
  assay = "RNA",
  features = integration_features,
  npcs = 30,
  verbose = TRUE
)

cat("\nPCA completed.\n\n")


# ============================================================
# 12. Elbow plot
# ============================================================

cat("Saving PCA elbow plot...\n")

png(
  filename = elbow_file,
  width = 1200,
  height = 900,
  res = 150
)

print(
  ElbowPlot(
    object,
    ndims = 30
  )
)

dev.off()

cat(
  "Saved:",
  elbow_file,
  "\n\n"
)


# ============================================================
# 13. CCA integration
# ============================================================

cat("============================================================\n")
cat("SEURAT V5 CCA INTEGRATION\n")
cat("============================================================\n\n")

cat("Integration parameters:\n")
cat("  method            = CCAIntegration\n")
cat("  dimensions        = 1:20\n")
cat("  dimensions output = 20\n")
cat("  k.weight          = 20\n")
cat("  nfeatures         = 2000\n\n")


# ------------------------------------------------------------
# Balanced integration tree
#
# 8 datasets:
#
#   1 + 2
#   3 + 4
#   5 + 6
#   7 + 8
#
#   then pair the four merged groups.
# ------------------------------------------------------------

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

cat("Integration tree:\n")
print(sample_tree)
cat("\n")


# ============================================================
# 14. Run IntegrateLayers
# ============================================================

integration_error <- NULL

integrated_object <- tryCatch(

  {

    IntegrateLayers(
      object = object,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      dims = 1:20,
      k.weight = 20,
      sample.tree = sample_tree,
      preserve.order = TRUE,
      verbose = TRUE
    )

  },

  error = function(e) {

    integration_error <<- e

    NULL

  }

)


# ============================================================
# 15. Integration error diagnostics
# ============================================================

if (is.null(integrated_object)) {

  cat("\n")
  cat("============================================================\n")
  cat("CCA INTEGRATION ERROR\n")
  cat("============================================================\n\n")

  cat(
    conditionMessage(integration_error),
    "\n\n"
  )

  cat("Diagnostics before integration:\n\n")

  cat("Number of cells:", ncol(object), "\n")
  cat("Number of genes:", nrow(object), "\n\n")

  cat("RNA layers:\n")
  print(Layers(object[["RNA"]]))

  cat("\n")

  if ("sample" %in% colnames(object@meta.data)) {

    cat("Cells per sample:\n")
    print(table(object$sample))

    cat("\n")

  }

  cat("PCA dimensions:\n")

  if ("pca" %in% Reductions(object)) {

    print(
      dim(
        Embeddings(
          object,
          "pca"
        )
      )
    )

  } else {

    cat("PCA reduction not found.\n")

  }

  cat("\n")

  stop(
    "CCA integration failed."
  )

}


# ============================================================
# 16. Use integrated object
# ============================================================

object <- integrated_object

cat("\n")
cat("CCA integration completed successfully.\n\n")


# ============================================================
# 17. Check integrated reduction
# ============================================================

if (!"integrated.cca" %in% Reductions(object)) {

  stop(
    "integrated.cca reduction was not created."
  )

}

cat("Integrated reduction found.\n")

cat(
  "Integrated dimensions:",
  ncol(
    Embeddings(
      object,
      "integrated.cca"
    )
  ),
  "\n\n"
)


# ============================================================
# 18. Neighbors
# ============================================================

cat("============================================================\n")
cat("NEIGHBOR GRAPH\n")
cat("============================================================\n\n")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  graph.name = c(
    "integrated.cca.nn",
    "integrated.cca.snn"
  ),
  verbose = TRUE
)

cat("\nNeighbor graph completed.\n\n")


# ============================================================
# 19. Clustering
# ============================================================

cat("============================================================\n")
cat("CLUSTERING\n")
cat("============================================================\n\n")

object <- FindClusters(
  object,
  graph.name = "integrated.cca.snn",
  resolution = 0.4,
  verbose = TRUE
)

cat("\nClustering completed.\n\n")


# ============================================================
# 20. UMAP
# ============================================================

cat("============================================================\n")
cat("UMAP\n")
cat("============================================================\n\n")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:20,
  reduction.name = "umap.integrated",
  reduction.key = "UMAPINT_",
  verbose = TRUE
)

cat("\nUMAP completed.\n\n")


# ============================================================
# 21. Save UMAP
# ============================================================

cat("Saving integrated UMAP...\n")

p_umap <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "sample",
  label = FALSE
) +
  ggtitle(
    "Seurat v5 CCA integration - samples"
  )

ggsave(
  filename = umap_file,
  plot = p_umap,
  width = 10,
  height = 7,
  dpi = 150
)

cat(
  "Saved:",
  umap_file,
  "\n\n"
)


# ============================================================
# 22. Save integrated object
# ============================================================

cat("============================================================\n")
cat("SAVING\n")
cat("============================================================\n\n")

saveRDS(
  object,
  output_file
)

cat(
  "Saved:",
  output_file,
  "\n\n"
)


# ============================================================
# 23. Final summary
# ============================================================

cat("============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n\n")

cat("Cells:", ncol(object), "\n")
cat("Genes:", nrow(object), "\n")

cat(
  "Assays:",
  paste(
    Assays(object),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Reductions:",
  paste(
    Reductions(object),
    collapse = ", "
  ),
  "\n\n"
)

cat("Clusters:\n")
print(
  table(
    Idents(object)
  )
)

cat("\n")

if ("sample" %in% colnames(object@meta.data)) {

  cat("Cells per sample:\n")
  print(
    table(
      object$sample
    )
  )

  cat("\n")

}

cat("============================================================\n")
cat("05 - COMPLETE\n")
cat("============================================================\n")
