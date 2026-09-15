# ============================================================
# 03_normalization_HVG.R
#
# LogNormalize + highly variable genes
#
# No ScaleData here.
# Scaling is performed after layer splitting in script 05.
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

if (!file.exists(input)) {
  stop("Input object not found: ", input)
}

message("============================================================")
message("Loading QC-filtered object")
message("============================================================")

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop("Input is not a Seurat object.")
}

DefaultAssay(object) <- "RNA"

message("Cells: ", ncol(object))
message("Features: ", nrow(object))

# ============================================================
# VALIDATE RNA LAYERS
# ============================================================

layers <- Layers(object[["RNA"]])

message("RNA layers before normalization:")
print(layers)

if (any(grepl("SeuratProject", layers))) {
  stop(
    "Invalid SeuratProject layers detected."
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
# HIGHLY VARIABLE GENES
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

hvg <- VariableFeatures(object)

message(
  "Number of HVGs: ",
  length(hvg)
)

if (length(hvg) != 2000) {
  warning(
    "Expected 2,000 HVGs but obtained ",
    length(hvg)
  )
}

# ============================================================
# HVG PLOT
# ============================================================

top10 <- head(
  hvg,
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
    gene = hvg
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

output <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
)

saveRDS(
  object,
  output
)

message("============================================================")
message("Normalization and HVG selection completed")
message("Output: ", output)
message("============================================================")
