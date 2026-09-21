# ============================================================
# 10b_final_figures.R
#
# Final portfolio figures
#
# ============================================================

source("R/functions.R")

library(Seurat)
library(ggplot2)
library(patchwork)

project_dir <- get_project_dir()

input <- file.path(
  project_dir,
  "results",
  "objects",
  "09_DE_ready.rds"
)

if (!file.exists(input)) {
  stop(
    "Input object not found: ",
    input
  )
}

object <- readRDS(input)

# ============================================================
# INTEGRATED CLUSTERS
# ============================================================

if ("umap.integrated" %in% Reductions(object)) {

  p_cluster <- DimPlot(
    object,
    reduction = "umap.integrated",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(
      "Integrated scRNA-seq clusters"
    )

  ggsave(
    file.path(
      project_dir,
      "figures",
      "integration",
      "FINAL_integrated_clusters.png"
    ),
    p_cluster,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# TREATMENT
# ============================================================

if (
  "treatment" %in%
  colnames(object@meta.data) &&
  "umap.integrated" %in% Reductions(object)
) {

  p_treatment <- DimPlot(
    object,
    reduction = "umap.integrated",
    group.by = "treatment"
  ) +
    ggtitle(
      "Treatment"
    )

  ggsave(
    file.path(
      project_dir,
      "figures",
      "integration",
      "FINAL_integrated_treatment.png"
    ),
    p_treatment,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# ANNOTATION COMPARISON
# ============================================================

if (
  all(
    c(
      "SingleR_label",
      "scPred_label",
      "annotation_consensus"
    ) %in%
      colnames(object@meta.data)
  ) &&
  "umap.integrated" %in% Reductions(object)
) {

  p1 <- DimPlot(
    object,
    reduction = "umap.integrated",
    group.by = "SingleR_label",
    label = TRUE,
    repel = TRUE,
    na.value = "grey80"
  ) +
    ggtitle("SingleR")

  p2 <- DimPlot(
    object,
    reduction = "umap.integrated",
    group.by = "scPred_label",
    label = TRUE,
    repel = TRUE,
    na.value = "grey80"
  ) +
    ggtitle("scPred")

  p3 <- DimPlot(
    object,
    reduction = "umap.integrated",
    group.by = "annotation_consensus",
    label = TRUE,
    repel = TRUE,
    na.value = "grey80"
  ) +
    ggtitle(
      "SingleR + scPred consensus"
    )

  annotation_panel <- (
    p1 | p2
  ) / p3

  ggsave(
    file.path(
      project_dir,
      "figures",
      "annotation",
      "FINAL_annotation_comparison.png"
    ),
    annotation_panel,
    width = 14,
    height = 11,
    dpi = 300
  )
}

# ============================================================
# QC
# ============================================================

if (
  all(
    c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "sample"
    ) %in%
      colnames(object@meta.data)
  )
) {

  p_qc <- VlnPlot(
    object,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    group.by = "sample",
    ncol = 3,
    pt.size = 0.03
  )

  ggsave(
    file.path(
      project_dir,
      "figures",
      "qc",
      "FINAL_QC.png"
    ),
    p_qc,
    width = 15,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# CONSENSUS CELL COUNTS
# ============================================================

if (
  "annotation_consensus" %in%
  colnames(object@meta.data)
) {

  counts <- as.data.frame(
    table(
      object$annotation_consensus
    )
  )

  colnames(counts) <- c(
    "cell_type",
    "cells"
  )

  p_counts <- ggplot(
    counts,
    aes(
      x = reorder(cell_type, cells),
      y = cells
    )
  ) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Consensus cell-type composition",
      x = "Cell type",
      y = "Number of cells"
    )

  ggsave(
    file.path(
      project_dir,
      "figures",
      "annotation",
      "FINAL_consensus_cell_counts.png"
    ),
    p_counts,
    width = 9,
    height = 6,
    dpi = 300
  )
}

# ============================================================
# SESSION INFO
# ============================================================

save_session_info(
  project_dir
)

message("============================================================")
message("Final figures generated successfully")
message("============================================================")
