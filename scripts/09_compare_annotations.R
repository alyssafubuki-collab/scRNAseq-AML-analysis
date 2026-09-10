# ============================================================
# 09_compare_annotations.R
# Compare SingleR, scPred and scAnnoX
# ============================================================

source("R/functions.R")

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

project_dir <- get_project_dir()

object <- readRDS(
  file.path(
    project_dir,
    "results",
    "objects",
    "08_scAnnoX.rds"
  )
)

meta <- object@meta.data

required <- c(
  "SingleR_label",
  "scPred_label",
  "scAnnoX_label"
)

missing <- setdiff(
  required,
  colnames(meta)
)

if (length(missing) > 0) {
  stop(
    "Missing annotation columns: ",
    paste(missing, collapse = ", "),
    "\nRun scripts 06, 07 and 08 first."
  )
}

annotation_table <- meta %>%
  select(
    SingleR = SingleR_label,
    scPred = scPred_label,
    scAnnoX = scAnnoX_label
  )

write.csv(
  annotation_table,
  file.path(
    project_dir,
    "results",
    "annotation",
    "annotation_comparison_per_cell.csv"
  ),
  row.names = TRUE
)

# Pairwise agreement
agreement <- data.frame(
  comparison = c(
    "SingleR_vs_scPred",
    "SingleR_vs_scAnnoX",
    "scPred_vs_scAnnoX"
  ),
  agreement = c(
    mean(
      annotation_table$SingleR ==
        annotation_table$scPred,
      na.rm = TRUE
    ),
    mean(
      annotation_table$SingleR ==
        annotation_table$scAnnoX,
      na.rm = TRUE
    ),
    mean(
      annotation_table$scPred ==
        annotation_table$scAnnoX,
      na.rm = TRUE
    )
  )
)

write.csv(
  agreement,
  file.path(
    project_dir,
    "results",
    "annotation",
    "pairwise_annotation_agreement.csv"
  ),
  row.names = FALSE
)

# Consensus label:
# retain a label only when at least two methods agree.
annotation_table$consensus <- apply(
  annotation_table,
  1,
  function(x) {

    x <- x[
      !is.na(x) &
      x != "" &
      x != "unassigned"
    ]

    if (length(x) < 2) {
      return("unassigned")
    }

    tab <- sort(
      table(x),
      decreasing = TRUE
    )

    if (tab[1] >= 2) {
      names(tab)[1]
    } else {
      "unassigned"
    }
  }
)

object$annotation_consensus <- annotation_table$consensus

consensus_counts <- as.data.frame(
  table(
    object$annotation_consensus
  )
)

colnames(consensus_counts) <- c(
  "cell_type",
  "cells"
)

write.csv(
  consensus_counts,
  file.path(
    project_dir,
    "results",
    "annotation",
    "consensus_annotation_counts.csv"
  ),
  row.names = FALSE
)

p <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "annotation_consensus",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "annotation",
    "annotation_consensus_UMAP.png"
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
    "09_annotation_comparison.rds"
  )
)

message("Annotation comparison completed.")
