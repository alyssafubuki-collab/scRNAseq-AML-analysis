# ============================================================
# 11_final_figures.R
# Final portfolio figures
# ============================================================

source("R/functions.R")

library(Seurat)
library(ggplot2)
library(patchwork)

project_dir <- get_project_dir()

object <- readRDS(
  file.path(
    project_dir,
    "results",
    "objects",
    "10_DE_ready.rds"
  )
)

# ------------------------------------------------------------
# Integrated UMAP
# ------------------------------------------------------------

p_cluster <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Integrated scRNA-seq clusters")

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

# ------------------------------------------------------------
# Treatment
# ------------------------------------------------------------

p_treatment <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "treatment"
) +
  ggtitle("Treatment")

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

# ------------------------------------------------------------
# Annotation comparison
# ------------------------------------------------------------

p1 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "SingleR_label"
) +
  ggtitle("SingleR")

p2 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "scPred_label"
) +
  ggtitle("scPred")

p3 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "scAnnoX_label"
) +
  ggtitle("scAnnoX")

p4 <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "annotation_consensus"
) +
  ggtitle("Consensus")

annotation_panel <- (
  p1 | p2
) / (
  p3 | p4
)

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

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

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

save_session_info(project_dir)

message("Final figures generated.")
