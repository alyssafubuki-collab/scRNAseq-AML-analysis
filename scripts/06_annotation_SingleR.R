# ============================================================
# 06_annotation_SingleR.R
#
# SingleR annotation using the
# Novershtern hematopoietic reference
# ============================================================

source("R/functions.R")

library(Seurat)
library(SeuratObject)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)

project_dir <- get_project_dir()

input <- file.path(
  project_dir,
  "results",
  "objects",
  "05_integrated.rds"
)

if (!file.exists(input)) {
  stop("Input object not found: ", input)
}

object <- readRDS(input)

if (!inherits(object, "Seurat")) {
  stop("Input is not a Seurat object.")
}

DefaultAssay(object) <- "RNA"

# ============================================================
# JOIN RNA LAYERS
#
# scPred will operate on the complete query object.
# Joining after integration is supported by Seurat v5.
# ============================================================

message("============================================================")
message("Joining RNA layers")
message("============================================================")

object[["RNA"]] <- JoinLayers(
  object[["RNA"]]
)

message("RNA layers:")
print(
  Layers(object[["RNA"]])
)

# ============================================================
# CONVERT TO SINGLECELLEXPERIMENT
# ============================================================

message("============================================================")
message("Preparing SingleR input")
message("============================================================")

sce <- as.SingleCellExperiment(
  object,
  assay = "RNA"
)

# ============================================================
# REFERENCE
# ============================================================

message("============================================================")
message("Loading Novershtern hematopoietic reference")
message("============================================================")

ref <- celldex::NovershternHematopoieticData(
  cell.ont = "all"
)

# ============================================================
# SINGLER
# ============================================================

message("============================================================")
message("Running SingleR")
message("============================================================")

pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

# ============================================================
# STORE RESULTS
# ============================================================

object$SingleR_label <- pred$pruned.labels
object$SingleR_label_raw <- pred$labels
object$SingleR_delta <- pred$delta.next

prediction_table <- data.frame(
  cell = rownames(pred),
  label = pred$pruned.labels,
  raw_label = pred$labels,
  delta = pred$delta.next
)

write.csv(
  prediction_table,
  file.path(
    project_dir,
    "results",
    "annotation",
    "SingleR_predictions.csv"
  ),
  row.names = FALSE
)

# ============================================================
# UMAP
# ============================================================

if ("umap.integrated" %in% Reductions(object)) {

  p <- DimPlot(
    object,
    reduction = "umap.integrated",
    group.by = "SingleR_label",
    label = TRUE,
    repel = TRUE,
    na.value = "grey80"
  )

  ggsave(
    file.path(
      project_dir,
      "results",
      "figures",
      "annotation",
      "SingleR_UMAP.png"
    ),
    p,
    width = 10,
    height = 7,
    dpi = 300
  )
}

# ============================================================
# SAVE
# ============================================================

output <- file.path(
  project_dir,
  "results",
  "objects",
  "06_SingleR.rds"
)

saveRDS(
  object,
  output
)

message("============================================================")
message("SingleR annotation completed")
message("Output: ", output)
message("============================================================")
