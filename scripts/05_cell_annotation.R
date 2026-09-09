# ============================================================
# 05_cell_annotation.R
# Cell-type annotation
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

input_file <- file.path(
  results_dir,
  "03_clustered_seurat_object.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Clustered object not found. ",
    "Run 03_pca_clustering_umap.R first."
  )
}

seurat_obj <- readRDS(input_file)

# ------------------------------------------------------------
# Marker genes
# ------------------------------------------------------------

marker_genes <- list(

  T_cell = c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "IL7R"
  ),

  NK_cell = c(
    "NKG7",
    "GNLY",
    "KLRD1"
  ),

  B_cell = c(
    "CD79A",
    "MS4A1",
    "CD37"
  ),

  Myeloid = c(
    "LYZ",
    "S100A8",
    "S100A9",
    "CTSS"
  ),

  HSC_like = c(
    "HOPX",
    "MALAT1",
    "GATA2",
    "AVP"
  ),

  Progenitor = c(
    "MPO",
    "CD34",
    "KIT"
  )

)

# ------------------------------------------------------------
# DotPlot
# ------------------------------------------------------------

available_markers <- unique(
  unlist(marker_genes)
)

available_markers <- available_markers[
  available_markers %in%
    rownames(seurat_obj)
]

p_markers <- DotPlot(
  seurat_obj,
  features = available_markers
) +
  RotatedAxis()

ggsave(
  file.path(
    figures_dir,
    "annotation",
    "marker_dotplot.png"
  ),
  p_markers,
  width = 12,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# Feature plots
# ------------------------------------------------------------

selected_markers <- intersect(
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

if (length(selected_markers) > 0) {

  p_features <- FeaturePlot(
    seurat_obj,
    features = selected_markers,
    ncol = 3
  )

  ggsave(
    file.path(
      figures_dir,
      "annotation",
      "marker_featureplots.png"
    ),
    p_features,
    width = 12,
    height = 8,
    dpi = 300
  )

}

# ------------------------------------------------------------
# Manual annotation
# ------------------------------------------------------------

# Example:
#
# cluster_annotations <- c(
#   "HSC",
#   "T cell",
#   "NK cell",
#   "Myeloid",
#   "Progenitor"
# )
#
# names(cluster_annotations) <- levels(seurat_obj)
#
# seurat_obj <- RenameIdents(
#   seurat_obj,
#   cluster_annotations
# )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file.path(
    results_dir,
    "05_annotated_seurat_object.rds"
  )
)

message(
  "Annotation marker analysis completed."
)
