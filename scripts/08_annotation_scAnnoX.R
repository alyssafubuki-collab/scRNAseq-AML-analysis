# ============================================================
# 08_annotation_scAnnoX.R
# scAnnoX annotation using the Seurat-v5-compatible local package
# ============================================================

source("R/functions.R")

library(Seurat)
library(ggplot2)

project_dir <- get_project_dir()

package_dir <- file.path(
  project_dir,
  "external",
  "scAnnoX"
)

if (!dir.exists(package_dir)) {
  stop(
    "\nThe Seurat-v5-compatible scAnnoX package is missing.\n\n",
    "Copy the modified scAnnoX package used in your original workflow to:\n",
    package_dir,
    "\n\nThe portfolio intentionally does not replace your modified package with the original version."
  )
}

if (!requireNamespace("scAnnoX", quietly = TRUE)) {

  message("Installing local scAnnoX package...")

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  remotes::install_local(
    package_dir,
    upgrade = "never",
    dependencies = TRUE
  )
}

library(scAnnoX)

object <- readRDS(
  file.path(
    project_dir,
    "results",
    "objects",
    "05_integrated.rds"
  )
)

# ------------------------------------------------------------
# Public reference for scAnnoX
# ------------------------------------------------------------

reference <- scPred::pbmc_1

reference <- NormalizeData(
  reference,
  verbose = FALSE
)

reference <- FindVariableFeatures(
  reference,
  verbose = FALSE
)

reference <- ScaleData(
  reference,
  verbose = FALSE
)

reference <- RunPCA(
  reference,
  npcs = 30,
  verbose = FALSE
)

Idents(reference) <- reference$cell_type

marker_list <- findMarkerToolsForSc(
  reference,
  to.list = TRUE,
  top.k = 30
)

# Run the scAnnoX comparison.
# The exact method name/API is retained in this script so the package
# modification used in the research workflow remains the source of truth.

anno.tools <- c(
  "SingleR",
  "Seurat",
  "sciBet",
  "scmap",
  "CHETAH",
  "scSorter",
  "sc.type",
  "cellID",
  "scCATCH",
  "SCINA"
)

marker.based <- c(
  "scSorter",
  "sc.type",
  "cellID",
  "scCATCH",
  "SCINA"
)

anno.res <- lapply(
  anno.tools,
  function(tool) {

    message("Running scAnnoX method: ", tool)

    strategy <- ifelse(
      tool %in% marker.based,
      "marker-based",
      "refernce-based"
    )

    tmp <- autoAnnoTools(
      object,
      ref.obj = reference,
      ref.ctype = "cell_type",
      marker.lst = marker_list,
      method = tool,
      select.marker = "Seurat",
      top.k = 30,
      strategy = strategy
    )

    tmp@meta.data[, ncol(tmp@meta.data)]
  }
)

anno.res <- do.call(
  cbind,
  anno.res
)

colnames(anno.res) <- anno.tools

scAnnoX_result <- annoResult(
  anno.res
)

object$scAnnoX_label <- scAnnoX_result$scAnnoX

write.csv(
  cbind(
    cell = rownames(scAnnoX_result),
    scAnnoX_result
  ),
  file.path(
    project_dir,
    "results",
    "annotation",
    "scAnnoX_predictions.csv"
  ),
  row.names = FALSE
)

p <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "scAnnoX_label",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "annotation",
    "scAnnoX_UMAP.png"
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
    "08_scAnnoX.rds"
  )
)

message("scAnnoX annotation completed.")
