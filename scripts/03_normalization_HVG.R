# ============================================================
# 03_normalization_HVG.R
# Normalization and variable feature selection
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
  "02_QC_filtered.rds"
)

# ============================================================
# LOAD
# ============================================================

message("============================================================")
message("Loading QC-filtered object")
message("============================================================")

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop("Input is not a Seurat object.")
}

message("Cells: ", ncol(object))
message("Features: ", nrow(object))

# ============================================================
# RNA LAYERS
# ============================================================

rna_layers <- Layers(
  object[["RNA"]]
)

message("RNA layers:")
print(rna_layers)

if (any(grepl("SeuratProject", rna_layers))) {
  stop(
    "Invalid SeuratProject layers detected before normalization."
  )
}

# ============================================================
# NORMALIZATION
# ============================================================

message("============================================================")
message("LogNormalize")
message("============================================================")

object <- NormalizeData(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

gc()

# ============================================================
# HVG
# ============================================================

message("============================================================")
message("Selecting 2,000 highly variable genes")
message("============================================================")

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

message(
  "HVGs: ",
  length(VariableFeatures(object))
)

# ============================================================
# PLOT
# ============================================================

top10 <- head(
  VariableFeatures(object),
  10
)

p <- VariableFeaturePlot(
  object
)

p <- LabelPoints(
  plot = p,
  points = top10,
  repel = TRUE
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "qc",
    "highly_variable_genes.png"
  ),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

write.csv(
  data.frame(
    gene = VariableFeatures(object)
  ),
  file.path(
    project_dir,
    "results",
    "qc",
    "highly_variable_genes.csv"
  ),
  row.names = FALSE
)

# ============================================================
# SAVE
# ============================================================

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "03_normalized.rds"
  )
)

message("============================================================")
message("Normalization and HVG selection completed")
message("============================================================")
