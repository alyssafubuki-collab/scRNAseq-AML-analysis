# ============================================================
# 08_consensus_annotation.R
# Consensus annotation using SingleR + scPred
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

if (!file.exists(input)) {
  stop("Input object not found: ", input)
}

object <- readRDS(input)

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

annotation_table <- object@meta.data %>%
  select(
    SingleR = SingleR_label,
    scPred = scPred_label
  )

annotation_table$SingleR[
  is.na(annotation_table$SingleR)
] <- "unassigned"

annotation_table$scPred[
  is.na(annotation_table$scPred)
] <- "unassigned"

valid <- (
  annotation_table$SingleR != "unassigned" &
  annotation_table$scPred != "unassigned"
)

if (sum(valid) == 0) {
  agreement <- NA_real_
} else {
  agreement <- mean(
    annotation_table$SingleR[valid] ==
      annotation_table$scPred[valid]
  )
}

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

# Conservative two-method consensus:
# identical non-unassigned labels are retained;
# disagreements remain unassigned.
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

object$annotation_consensus <-
  annotation_table$consensus

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

consensus_counts <- as.data.frame(
  table(object$annotation_consensus)
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
    "08_annotation_consensus.rds"
  )
)

message(
  "Consensus annotation completed."
)

message(
  "SingleR/scPred agreement: ",
  ifelse(
    is.na(agreement),
    "NA",
    paste0(round(agreement * 100, 2), "%")
  )
)
