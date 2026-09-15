# ============================================================
# 03_normalization_HVG.R
# Log-normalization and highly variable gene selection
#
# Deliberately NO ScaleData() here.
# Scaling is performed only on the 2,000 HVGs in scripts 04/05.
# This avoids creating a large all-gene scale.data matrix.
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

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop("Input is not a Seurat object.")
}

DefaultAssay(object) <- "RNA"

message("============================================================")
message("Loading QC-filtered object")
message("============================================================")
message("Cells: ", ncol(object))
message("Features: ", nrow(object))

rna_layers <- Layers(object[["RNA"]])

message("RNA layers:")
print(rna_layers)

if (any(grepl("SeuratProject", rna_layers))) {
  stop(
    "Invalid SeuratProject layers detected before normalization."
  )
}

if (length(rna_layers) > 2L) {
  stop(
    "RNA assay is already split before integration. ",
    "The layer split must occur only in script 05."
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
  "HVGs selected: ",
  length(VariableFeatures(object))
)

# ============================================================
# HVG PLOT
# ============================================================

top10 <- head(
  VariableFeatures(object),
  10
)

p <- VariableFeaturePlot(object)

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

output_file <- file.path(
  project_dir,
  "results",
  "objects",
  "03_normalized.rds"
)

saveRDS(
  object,
  output_file
)

message("============================================================")
message("Normalization and HVG selection completed")
message("Output: ", output_file)
message("============================================================")
