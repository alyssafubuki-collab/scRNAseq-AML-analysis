# ============================================================
# 07_annotation_scPred.R
#
# scPred supervised annotation
#
# Input:
#   06_SingleR.rds
#
# scPredict performs its own Harmony-based alignment.
#
# ------------------------------------------------------------
# NOTE (see chat discussion): scPred's project_query() internally
# calls GetAssayData(new, "data"), passing "data" as the SECOND
# POSITIONAL argument. scPred was written against an older
# Seurat/SeuratObject API where that position was `slot`. In
# current SeuratObject (v5), the second positional argument is
# `assay`, not `slot`/`layer`, so "data" gets validated as an
# assay name instead -- which fails because the object only has
# an assay named "RNA":
#   Error: `assay` must be one of "RNA", not "data".
#
# Workaround (carried over from a previous analysis): duplicate
# the RNA assay under the literal name "data" and make it the
# default assay before calling scPredict(). This makes "data" a
# valid assay name, so GetAssayData(new, "data") now resolves
# correctly instead of erroring -- no need to patch scPred's
# source. The duplicate assay is removed again right after
# scPredict() finishes so it isn't carried into the saved object.
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
library(scPred)
library(harmony)
library(ggplot2)

project_dir <- get_project_dir()

input <- file.path(
  project_dir,
  "results",
  "objects",
  "06_SingleR.rds"
)

if (!file.exists(input)) {
  stop(
    "Input object not found: ",
    input
  )
}

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop(
    "Input is not a Seurat object."
  )
}

# ============================================================
# VALIDATE SINGLE R
# ============================================================

if (!"SingleR_label" %in% colnames(object@meta.data)) {
  stop(
    "SingleR_label is missing.",
    "\nThe scPred script must receive 06_SingleR.rds."
  )
}

# ============================================================
# REFERENCE
# ============================================================

message("============================================================")
message("Loading scPred reference")
message("============================================================")

reference <- scPred::pbmc_1

# ============================================================
# PREPROCESS REFERENCE
# ============================================================

message("============================================================")
message("Preparing scPred reference")
message("============================================================")

reference <- NormalizeData(
  reference,
  verbose = FALSE
)

reference <- FindVariableFeatures(
  reference,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

reference <- ScaleData(
  reference,
  features = VariableFeatures(reference),
  verbose = FALSE
)

reference <- RunPCA(
  reference,
  features = VariableFeatures(reference),
  npcs = 30,
  verbose = FALSE
)

# ============================================================
# BUILD FEATURE SPACE
# ============================================================

message("============================================================")
message("Building scPred feature space")
message("============================================================")

reference <- getFeatureSpace(
  reference,
  "cell_type"
)

# ============================================================
# TRAIN CLASSIFIERS
# ============================================================

message("============================================================")
message("Training scPred classifiers")
message("============================================================")

reference <- trainModel(
  reference,
  tuneLength = 1,
  number = 3,
  seed = 66,
  allowParallel = FALSE
)

# ============================================================
# PREPARE QUERY
# ============================================================

message("============================================================")
message("Preparing query object")
message("============================================================")

DefaultAssay(object) <- "RNA"

object <- NormalizeData(
  object,
  verbose = FALSE
)

# ============================================================
# SCPRED COMPATIBILITY WORKAROUND
# ============================================================
# See note at top of file: duplicate the RNA assay under the name
# "data" and set it as default so scPred's positional
# GetAssayData(new, "data") call resolves as a valid assay name.

object[["data"]] <- object[["RNA"]]

DefaultAssay(object) <- "data"

# ============================================================
# SCPRED
# ============================================================

message("============================================================")
message("Running scPred")
message("============================================================")

object <- scPredict(
  object,
  reference,
  threshold = 0.55,
  max.iter.harmony = 20,
  recompute_alignment = TRUE,
  seed = 66
)

# ============================================================
# REMOVE WORKAROUND ASSAY
# ============================================================
# Drop the duplicate "data" assay now that scPredict has run, so it
# isn't carried into the saved object and doesn't double the object's
# size on disk. Restore RNA as the default assay for downstream scripts.

object[["data"]] <- NULL

DefaultAssay(object) <- "RNA"

# ============================================================
# VALIDATE SCPRED OUTPUT
# ============================================================

if (!"scpred_prediction" %in% colnames(object@meta.data)) {
  stop(
    "scPred did not generate scpred_prediction."
  )
}

if (!"scpred" %in% Reductions(object)) {
  warning(
    "scpred reduction was not found."
  )
}

# ============================================================
# STORE LABEL
# ============================================================

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

# ============================================================
# SCPRED UMAP
# ============================================================

if ("scpred" %in% Reductions(object)) {

  p <- DimPlot(
    object,
    reduction = "scpred",
    group.by = "scPred_label",
    label = TRUE,
    repel = TRUE,
    na.value = "grey80"
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
}

# ============================================================
# SAVE TRAINED REFERENCE
# ============================================================

saveRDS(
  reference,
  file.path(
    project_dir,
    "results",
    "annotation",
    "scPred_reference_trained.rds"
  )
)

# ============================================================
# SAVE QUERY
# ============================================================

output <- file.path(
  project_dir,
  "results",
  "objects",
  "07_scPred.rds"
)

saveRDS(
  object,
  output
)

message("============================================================")
message("scPred annotation completed")
message("Output: ", output)
message("============================================================")
