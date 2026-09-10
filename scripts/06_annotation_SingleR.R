# ============================================================
# 06_annotation_SingleR.R
# SingleR annotation using the Novershtern hematopoietic reference
# ============================================================

source("R/functions.R")

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
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

# Join RNA layers so the expression matrix can be exported cleanly.
object <- JoinLayers(
  object,
  assay = "RNA"
)

sce <- as.SingleCellExperiment(
  object,
  assay = "RNA"
)

# The Novershtern reference is specifically hematopoietic and contains
# HSCs, CMPs, MEPs, GMPs, granulocytes, monocytes, NK and T/B populations.
ref <- celldex::NovershternHematopoieticData(
  cell.ont = "all"
)

pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

object$SingleR_label <- pred$pruned.labels
object$SingleR_label_raw <- pred$labels
object$SingleR_delta <- pred$delta.next

write.csv(
  data.frame(
    cell = rownames(pred),
    label = pred$pruned.labels,
    raw_label = pred$labels,
    delta = pred$delta.next
  ),
  file.path(
    project_dir,
    "results",
    "annotation",
    "SingleR_predictions.csv"
  ),
  row.names = FALSE
)

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
    "figures",
    "annotation",
    "SingleR_UMAP.png"
  ),
  p,
  width = 10,
  height = 7,
  dpi = 300
)

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "06_SingleR.rds"
  )
)

message("SingleR annotation completed.")
