# ============================================================
# functions.R
# Reusable helper functions
# ============================================================

calculate_qc_summary <- function(seurat_object) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  data.frame(
    cells = ncol(seurat_object),
    genes = nrow(seurat_object),
    median_features =
      median(seurat_object$nFeature_RNA),
    median_counts =
      median(seurat_object$nCount_RNA),
    median_percent_mt =
      median(seurat_object$percent.mt)
  )
}


save_seurat_object <- function(
    seurat_object,
    output_file
) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  saveRDS(
    seurat_object,
    output_file
  )

  message(
    "Seurat object saved to: ",
    output_file
  )
}


load_seurat_object <- function(
    input_file
) {

  if (!file.exists(input_file)) {
    stop(
      "File not found: ",
      input_file
    )
  }

  object <- readRDS(input_file)

  if (!inherits(object, "Seurat")) {
    stop(
      "The file does not contain a Seurat object."
    )
  }

  return(object)
}
