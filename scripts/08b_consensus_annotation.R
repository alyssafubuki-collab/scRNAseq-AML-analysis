# ============================================================
# 08b_consensus_annotation.R
#
# Consensus annotation:
#   SingleR + scPred
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
  stop(
    "Input object not found: ",
    input
  )
}

object <- readRDS(input)

# ============================================================
# VALIDATION
# ============================================================

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
    paste(missing, collapse = ", ")
  )
}

# ============================================================
# ANNOTATION TABLE
# ============================================================

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

# ============================================================
# AGREEMENT
# ============================================================

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

# ============================================================
# CONSERVATIVE CONSENSUS
#
# Agreement → label
# Disagreement → unassigned
# ============================================================

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

# ============================================================
# PER-CELL TABLE
# ============================================================

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

# ============================================================
# COUNTS
# ============================================================

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

# SingleR counts
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

# scPred counts
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

# ============================================================
# UMAP
# ============================================================

if ("umap.integrated" %in% Reductions(object)) {

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
}

# ============================================================
# SAVE
# ============================================================

output <- file.path(
  project_dir,
  "results",
  "objects",
  "08_annotation_consensus.rds"
)

saveRDS(
  object,
  output
)

message("============================================================")
message("Consensus annotation completed")
message("============================================================")

message(
  "SingleR/scPred agreement: ",
  ifelse(
    is.na(agreement),
    "NA",
    paste0(
      round(agreement * 100, 2),
      "%"
    )
  )
)

message(
  "Output: ",
  output
)
