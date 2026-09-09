# ============================================================
# 03_pca_clustering_umap.R
# PCA, clustering and UMAP
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

input_file <- file.path(
  results_dir,
  "02_normalized_seurat_object.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Normalized object not found. ",
    "Run 02_normalization_hvg.R first."
  )
}

seurat_obj <- readRDS(input_file)

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = 30,
  verbose = FALSE
)

# ------------------------------------------------------------
# PCA diagnostics
# ------------------------------------------------------------

pdf(
  file.path(
    figures_dir,
    "integration",
    "pca_elbow_plot.pdf"
  ),
  width = 7,
  height = 5
)

print(
  ElbowPlot(
    seurat_obj,
    ndims = 30
  )
)

dev.off()

# ------------------------------------------------------------
# Neighbors
# ------------------------------------------------------------

seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:30,
  verbose = FALSE
)

# ------------------------------------------------------------
# Clustering
# ------------------------------------------------------------

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.4,
  verbose = FALSE
)

# ------------------------------------------------------------
# UMAP
# ------------------------------------------------------------

seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:30,
  verbose = FALSE
)

# ------------------------------------------------------------
# UMAP plots
# ------------------------------------------------------------

p_cluster <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(
    figures_dir,
    "integration",
    "umap_clusters.png"
  ),
  p_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

p_feature <- FeaturePlot(
  seurat_obj,
  features = c(
    "MALAT1",
    "LYZ"
  ),
  ncol = 2
)

ggsave(
  file.path(
    figures_dir,
    "integration",
    "umap_marker_genes.png"
  ),
  p_feature,
  width = 10,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file.path(
    results_dir,
    "03_clustered_seurat_object.rds"
  )
)

message(
  "PCA, clustering and UMAP completed."
)
