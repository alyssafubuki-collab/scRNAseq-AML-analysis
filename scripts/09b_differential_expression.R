# ============================================================
# 09_differential_expression.R
#
# Differential expression after SingleR + scPred consensus.
# ============================================================

source("R/functions.R")

library(Seurat)
library(dplyr)
library(ggplot2)

project_dir <- get_project_dir()

object <- readRDS(
  file.path(
    project_dir,
    "results",
    "objects",
    "08_annotation_consensus.rds"
  )
)

# ------------------------------------------------------------
# JOIN SEURAT v5 RNA LAYERS
# ------------------------------------------------------------

object <- JoinLayers(
  object,
  assay = "RNA"
)

# ------------------------------------------------------------
# CLUSTER MARKERS
# ------------------------------------------------------------

Idents(object) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(
  object,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  cluster_markers,
  file.path(
    project_dir,
    "results",
    "differential_expression",
    "cluster_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# TOP CLUSTER MARKERS
# ------------------------------------------------------------

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

write.csv(
  top_markers,
  file.path(
    project_dir,
    "results",
    "differential_expression",
    "top_cluster_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# TREATMENT DE
# ------------------------------------------------------------

Idents(object) <- "treatment"

treatments <- levels(
  factor(object$treatment)
)

if (length(treatments) >= 2) {

  for (i in seq_len(length(treatments) - 1)) {

    for (j in (i + 1):length(treatments)) {

      a <- treatments[i]
      b <- treatments[j]

      message(
        "DE: ",
        a,
        " vs ",
        b
      )

      de <- FindMarkers(
        object,
        ident.1 = a,
        ident.2 = b,
        min.pct = 0.10,
        logfc.threshold = 0.25
      )

      de$gene <- rownames(de)

      filename <- paste0(
        "DE_",
        make.names(a),
        "_vs_",
        make.names(b),
        ".csv"
      )

      write.csv(
        de,
        file.path(
          project_dir,
          "results",
          "differential_expression",
          filename
        ),
        row.names = FALSE
      )
    }
  }
}

# ------------------------------------------------------------
# CONSENSUS ANNOTATION MARKERS
# ------------------------------------------------------------

if ("annotation_consensus" %in%
    colnames(object@meta.data)) {

  valid_cells <- (
    !is.na(object$annotation_consensus) &
    object$annotation_consensus != "unassigned"
  )

  annotation_object <- subset(
    object,
    cells = colnames(object)[valid_cells]
  )

  if (length(unique(
    annotation_object$annotation_consensus
  )) >= 2) {

    Idents(annotation_object) <-
      "annotation_consensus"

    annotation_markers <- FindAllMarkers(
      annotation_object,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )

    write.csv(
      annotation_markers,
      file.path(
        project_dir,
        "results",
        "differential_expression",
        "consensus_annotation_markers.csv"
      ),
      row.names = FALSE
    )
  }
}

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "09_DE_ready.rds"
  )
)

message(
  "Differential expression completed."
)
