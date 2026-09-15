# ============================================================
# 05_Seurat5_integration.R
# Memory-conscious Seurat v5 CCA integration
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

output <- file.path(
  project_dir,
  "results",
  "objects",
  "05_integrated.rds"
)

# ============================================================
# SETTINGS
# ============================================================

N_HVG <- 2000

N_PCS <- 20

CLUSTER_RESOLUTION <- 0.4

# ============================================================
# MEMORY SETTINGS
# ============================================================

options(
  future.globals.maxSize = 8 * 1024^3
)

set.seed(1234)

# ============================================================
# LOAD
# ============================================================

message("============================================================")
message("Loading normalized Seurat object")
message("============================================================")

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop("Input is not a Seurat object.")
}

message("Cells: ", ncol(object))
message("Features: ", nrow(object))

# ============================================================
# SAMPLE CHECK
# ============================================================

if (!"sample" %in% colnames(object@meta.data)) {
  stop(
    "Metadata column 'sample' is missing."
  )
}

object$sample <- factor(
  object$sample
)

message("============================================================")
message("Samples")
message("============================================================")

print(
  table(object$sample)
)

# ============================================================
# LAYER CHECK
# ============================================================

rna_layers <- Layers(
  object[["RNA"]]
)

message("============================================================")
message("RNA layers before split")
message("============================================================")

print(rna_layers)

if (any(grepl("SeuratProject", rna_layers))) {

  stop(
    paste0(
      "Invalid RNA layers detected.\n",
      paste(rna_layers, collapse = "\n"),
      "\n\nRegenerate 02_QC_filtered.rds using the corrected ",
      "02_import_and_QC.R."
    )
  )
}

# ============================================================
# SPLIT ONCE
# ============================================================

message("============================================================")
message("Splitting RNA by sample")
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
# HVGs
# ============================================================

message("============================================================")
message("Finding variable features")
message("============================================================")

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = N_HVG,
  verbose = FALSE
)

# ============================================================
# SCALE ONLY HVGs
# ============================================================

message("============================================================")
message("Scaling only variable features")
message("============================================================")

object <- ScaleData(
  object,
  features = VariableFeatures(object),
  verbose = FALSE
)

gc()

# ============================================================
# PCA
# ============================================================

message("============================================================")
message("Running PCA")
message("============================================================")

object <- RunPCA(
  object,
  features = VariableFeatures(object),
  npcs = N_PCS,
  verbose = FALSE
)

gc()

# ============================================================
# CCA INTEGRATION
# ============================================================

message("============================================================")
message("Starting Seurat v5 CCA integration")
message("============================================================")

message(
  "Dimensions: ",
  N_PCS
)

message(
  "HVGs: ",
  length(VariableFeatures(object))
)

message(
  "Cells: ",
  ncol(object)
)

gc()

# ------------------------------------------------------------
# IMPORTANT:
#
# CCAIntegration is performed in low-dimensional space.
#
# k.weight is reduced from the default 100 to 50 to reduce
# memory consumption during anchor weighting.
#
# dims.to.integrate is limited to 20.
# ------------------------------------------------------------

object <- IntegrateLayers(
  object = object,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:N_PCS,
  dims.to.integrate = N_PCS,
  k.weight = 50,
  verbose = TRUE
)

gc()

message("CCA integration completed.")

# ============================================================
# NEIGHBORS
# ============================================================

message("============================================================")
message("Finding neighbors on integrated space")
message("============================================================")

object <- FindNeighbors(
  object,
  reduction = "integrated.cca",
  dims = 1:N_PCS,
  verbose = FALSE
)

gc()

# ============================================================
# CLUSTERS
# ============================================================

message("============================================================")
message("Clustering")
message("============================================================")

object <- FindClusters(
  object,
  resolution = CLUSTER_RESOLUTION,
  verbose = FALSE
)

gc()

# ============================================================
# UMAP
# ============================================================

message("============================================================")
message("Running integrated UMAP")
message("============================================================")

object <- RunUMAP(
  object,
  reduction = "integrated.cca",
  dims = 1:N_PCS,
  reduction.name = "umap.integrated",
  reduction.key = "integratedUMAP_",
  verbose = FALSE
)

gc()

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

rm(p1)
gc()

# ============================================================
# UMAP BY TREATMENT
# ============================================================

if ("treatment" %in% colnames(object@meta.data)) {

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

  rm(p2)
  gc()
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

rm(p3)
gc()

# ============================================================
# SAVE
# ============================================================

message("============================================================")
message("Saving integrated object")
message("============================================================")

saveRDS(
  object,
  output,
  compress = TRUE
)

# ============================================================
# VALIDATION
# ============================================================

if (!file.exists(output)) {
  stop(
    "Integrated object was not created."
  )
}

message("============================================================")
message("CCA integration completed successfully")
message("Output: ", output)
message("============================================================")
