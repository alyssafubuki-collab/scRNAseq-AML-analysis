# ============================================================
# 08_consensus_annotation.R
#
# Consensus annotation using:
#   1. SingleR
#   2. scPred
#
# scAnnoX intentionally excluded.
# ============================================================

source("R/functions.R")

library(Seurat)
library(dplyr)
library(ggplot2)

project_dir <- get_project_dir()

input <- file.path(
  project_dir,
  "results",
  "objects",
  "07_scPred.rds"
)

object <- readRDS(input)

# ------------------------------------------------------------
# CHECK REQUIRED ANNOTATIONS
# ------------------------------------------------------------

required <- c(
  "SingleR_label",
  "scPred_label"
)

missing <- setdiff(
  required,
  colnames(object@meta.data)
)

if (length(missing) > 0) {
  stop(
    "Missing annotation columns: ",
    paste(missing, collapse = ", "),
    "\nRun scripts 06 and 07 first."
  )
}

# ------------------------------------------------------------
# CREATE ANNOTATION TABLE
# ------------------------------------------------------------

annotation_table <- object@meta.data %>%
  select(
    SingleR = SingleR_label,
    scPred = scPred_label
  )

# Convert NA to explicit unassigned state
annotation_table$SingleR[
  is.na(annotation_table$SingleR)
] <- "unassigned"

annotation_table$scPred[
  is.na(annotation_table$scPred)
] <- "unassigned"

# ------------------------------------------------------------
# PAIRWISE AGREEMENT
# ------------------------------------------------------------

valid <- (
  annotation_table$SingleR != "unassigned" &
  annotation_table$scPred != "unassigned"
)

agreement <- mean(
  annotation_table$SingleR[valid] ==
    annotation_table$scPred[valid]
)

agreement_table <- data.frame(
  comparison = "SingleR_vs_scPred",
  agreement = agreement,
  cells_compared = sum(valid)
)

write.csv(
  agreement_table,
  file.path(
    project_dir,
    "results",
    "annotation",
    "pairwise_annotation_agreement.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# CONSENSUS
#
# With only two annotation methods:
#
# - same prediction -> consensus
# - different predictions -> unassigned
# - one method unassigned -> use the other only if desired
#
# Here we use the conservative approach:
# disagreement = unassigned.
# ------------------------------------------------------------

annotation_table$consensus <- "unassigned"

same_prediction <- (
  annotation_table$SingleR ==
    annotation_table$scPred &
  annotation_table$SingleR != "unassigned"
)

annotation_table$consensus[
  same_prediction
] <- annotation_table$SingleR[
  same_prediction
]

# ------------------------------------------------------------
# OPTIONAL SINGLE-METHOD LABEL
#
# Keep the individual annotations available for inspection.
# ------------------------------------------------------------

object$annotation_consensus <-
  annotation_table$consensus

# ------------------------------------------------------------
# SAVE PER-CELL ANNOTATIONS
# ------------------------------------------------------------

write.csv(
  data.frame(
    cell = rownames(annotation_table),
    annotation_table
  ),
  file.path(
    project_dir,
    "results",
    "annotation",
    "annotation_comparison_per_cell.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# CONSENSUS COUNTS
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# SINGLE METHOD COUNTS
# ------------------------------------------------------------

singleR_counts <- as.data.frame(
  table(object$SingleR_label)
)

colnames(singleR_counts) <- c(
  "cell_type",
  "cells"
)

write.csv(
  singleR_counts,
  file.path(
    project_dir,
    "results",
    "annotation",
    "SingleR_annotation_counts.csv"
  ),
  row.names = FALSE
)

scPred_counts <- as.data.frame(
  table(object$scPred_label)
)

colnames(scPred_counts) <- c(
  "cell_type",
  "cells"
)

write.csv(
  scPred_counts,
  file.path(
    project_dir,
    "results",
    "annotation",
    "scPred_annotation_counts.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# CONSENSUS UMAP
# ------------------------------------------------------------

p <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "annotation_consensus",
  label = TRUE,
  repel = TRUE,
  na.value = "grey80"
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

# ------------------------------------------------------------
# SAVE OBJECT
# ------------------------------------------------------------

saveRDS(
  object,
  file.path(
    project_dir,
    "results",
    "objects",
    "08_annotation_consensus.rds"
  )
)

message(
  "Consensus annotation completed."
)

message(
  "SingleR/scPred agreement: ",
  round(agreement * 100, 2),
  "%"
)
