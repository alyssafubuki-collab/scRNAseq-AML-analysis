# ============================================================
# 03_normalization_HVG.R
# Normalization and variable feature selection
# ============================================================

source("R/functions.R")

library(Seurat)
library(ggplot2)

project_dir <- get_project_dir()

input <- file.path(
  project_dir,
  "results",
  "objects",
  "02_QC_filtered.rds"
)

object <- readRDS(input)

object <- NormalizeData(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# Scale all genes for the downstream PCA workflow.
object <- ScaleData(
  object,
  features = rownames(object),
  verbose = FALSE
)

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

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "03_normalized.rds"
  )
)

message("Normalization and HVG selection completed.")
