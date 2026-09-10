# ============================================================
# plotting_functions.R
# ============================================================

plot_qc <- function(object) {

  Seurat::VlnPlot(
    object,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    group.by = "sample",
    ncol = 3,
    pt.size = 0.05
  )
}


plot_umap_clusters <- function(object, reduction = "umap") {

  Seurat::DimPlot(
    object,
    reduction = reduction,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  )
}


plot_umap_group <- function(
    object,
    group,
    reduction = "umap"
) {

  Seurat::DimPlot(
    object,
    reduction = reduction,
    group.by = group
  )
}
