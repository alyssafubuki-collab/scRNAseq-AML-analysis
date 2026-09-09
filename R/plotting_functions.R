# ============================================================
# plotting_functions.R
# Reusable visualization functions
# ============================================================

plot_qc_metrics <- function(
    seurat_object
) {

  VlnPlot(
    seurat_object,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3,
    pt.size = 0.1
  )
}


plot_clusters <- function(
    seurat_object,
    reduction = "umap"
) {

  DimPlot(
    seurat_object,
    reduction = reduction,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  )
}


plot_sample_distribution <- function(
    seurat_object,
    reduction = "umap"
) {

  if (!"sample" %in% colnames(seurat_object@meta.data)) {

    stop(
      "The Seurat object does not contain a 'sample' metadata column."
    )

  }

  DimPlot(
    seurat_object,
    reduction = reduction,
    group.by = "sample"
  )
}


plot_marker_genes <- function(
    seurat_object,
    genes
) {

  genes <- genes[
    genes %in% rownames(seurat_object)
  ]

  if (length(genes) == 0) {

    stop(
      "None of the requested genes are present in the dataset."
    )

  }

  DotPlot(
    seurat_object,
    features = genes
  ) +
    RotatedAxis()
}
