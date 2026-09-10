# ============================================================
# 07_annotation_scPred.R
# scPred supervised annotation
# ============================================================

source("R/functions.R")

library(Seurat)
library(scPred)
library(ggplot2)

project_dir <- get_project_dir()

object <- readRDS(
  file.path(
    project_dir,
    "results",
    "objects",
    "05_integrated.rds"
  )
)

# ------------------------------------------------------------
# Public scPred reference
# ------------------------------------------------------------

reference <- scPred::pbmc_1

# The reference already contains cell_type labels.
# Re-run the standard preprocessing required by scPred.
reference <- NormalizeData(
  reference,
  verbose = FALSE
)

reference <- FindVariableFeatures(
  reference,
  verbose = FALSE
)

reference <- ScaleData(
  reference,
  verbose = FALSE
)

reference <- RunPCA(
  reference,
  npcs = 30,
  verbose = FALSE
)

# Build the supervised feature space.
reference <- getFeatureSpace(
  reference,
  "cell_type"
)

# Train classifiers.
# tuneLength=1 keeps this reproducible and lighter for a portfolio repo.
reference <- trainModel(
  reference,
  tuneLength = 1,
  number = 3,
  seed = 66,
  allowParallel = FALSE
)

# Query must use the same normalization strategy.
object <- NormalizeData(
  object,
  verbose = FALSE
)

# scPred uses the trained reference models and Harmony alignment.
object <- scPredict(
  object,
  reference,
  threshold = 0.55
)

object$scPred_label <- object$scpred_prediction

write.csv(
  data.frame(
    cell = colnames(object),
    label = object$scPred_label
  ),
  file.path(
    project_dir,
    "results",
    "annotation",
    "scPred_predictions.csv"
  ),
  row.names = FALSE
)

p <- DimPlot(
  object,
  reduction = "scpred",
  group.by = "scPred_label",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "annotation",
    "scPred_UMAP.png"
  ),
  p,
  width = 10,
  height = 7,
  dpi = 300
)

saveRDS(
  reference,
  file.path(
    project_dir,
    "results",
    "annotation",
    "scPred_reference_trained.rds"
  )
)

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "07_scPred.rds"
  )
)

message("scPred annotation completed.")
