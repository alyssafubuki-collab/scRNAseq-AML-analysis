# ============================================================
# 07_visualization.R
# Final visualization
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

input_file <- file.path(
  results_dir,
  "05_annotated_seurat_object.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Annotated object not found."
  )
}

seurat_obj <- readRDS(input_file)

# ------------------------------------------------------------
# UMAP
# ------------------------------------------------------------

p_umap <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "Single-cell RNA-seq clustering"
  )

ggsave(
  file.path(
    figures_dir,
    "annotation",
    "final_umap_clusters.png"
  ),
  p_umap,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

p_qc <- VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3,
  pt.size = 0.05
)

ggsave(
  file.path(
    figures_dir,
    "qc",
    "final_qc.png"
  ),
  p_qc,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Common marker visualization
# ------------------------------------------------------------

markers <- intersect(
  c(
    "CD3D",
    "NKG7",
    "LYZ",
    "CD79A",
    "CD34",
    "MPO"
  ),
  rownames(seurat_obj)
)

if (length(markers) > 0) {

  p_dot <- DotPlot(
    seurat_obj,
    features = markers
  ) +
    RotatedAxis()

  ggsave(
    file.path(
      figures_dir,
      "annotation",
      "final_marker_dotplot.png"
    ),
    p_dot,
    width = 10,
    height = 6,
    dpi = 300
  )

}

message(
  "Final visualization completed."
)
