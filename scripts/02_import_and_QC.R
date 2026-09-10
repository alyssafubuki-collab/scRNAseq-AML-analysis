# ============================================================
# 02_import_and_QC.R
# Import all 8 GSE145410 samples and perform adaptive QC
# ============================================================

source("R/functions.R")
source("R/plotting_functions.R")

library(Seurat)
library(SeuratObject)
library(scuttle)
library(ggplot2)
library(dplyr)

project_dir <- get_project_dir()
make_project_dirs(project_dir)

metadata <- get_gse145410_metadata()

extract_dir <- file.path(
  project_dir,
  "data",
  "GSE145410",
  "extracted"
)

objects <- list()

for (i in seq_len(nrow(metadata))) {

  gsm <- metadata$gsm[i]
  sample_name <- metadata$sample[i]

  sample_dir <- file.path(
    extract_dir,
    gsm
  )

  if (!dir.exists(sample_dir)) {
    stop(
      "Sample directory not found: ",
      sample_dir,
      "\nRun scripts/01_download_GSE145410.R first."
    )
  }

  message("Reading ", gsm, " - ", sample_name)

  counts <- read_gse145410_sample(sample_dir)

  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE145410",
    min.cells = 3,
    min.features = 0
  )

  obj$gsm <- gsm
  obj$sample <- sample_name
  obj$treatment <- metadata$treatment[i]
  obj$replicate <- metadata$replicate[i]

  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  objects[[gsm]] <- obj
}

# Merge all samples before QC.
# QC thresholds are nevertheless calculated per sample via batch=sample.
combined <- Reduce(
  function(x, y) merge(x, y),
  objects
)

combined$sample <- factor(combined$sample)

# ------------------------------------------------------------
# QC before filtering
# ------------------------------------------------------------

qc_before <- combined@meta.data %>%
  count(sample, name = "cells_before")

write.csv(
  qc_before,
  file.path(
    project_dir,
    "results",
    "qc",
    "cells_before_QC.csv"
  ),
  row.names = FALSE
)

p_before <- VlnPlot(
  combined,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  ncol = 3,
  pt.size = 0.05
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "qc",
    "QC_before_filtering.png"
  ),
  p_before,
  width = 15,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Adaptive QC with isOutlier()
# ------------------------------------------------------------

# Thresholds are calculated independently per sample.
qc_counts <- isOutlier(
  combined$nCount_RNA,
  nmads = 3,
  type = "lower",
  log = TRUE,
  batch = combined$sample
)

qc_features <- isOutlier(
  combined$nFeature_RNA,
  nmads = 3,
  type = "lower",
  log = TRUE,
  batch = combined$sample
)

qc_mito <- isOutlier(
  combined$percent.mt,
  nmads = 3,
  type = "higher",
  batch = combined$sample
)

combined$discard_counts <- qc_counts
combined$discard_features <- qc_features
combined$discard_mito <- qc_mito

combined$discard <- (
  qc_counts |
  qc_features |
  qc_mito
)

# ------------------------------------------------------------
# QC threshold table
# ------------------------------------------------------------

get_thresholds <- function(x, metric_name) {

  thresholds <- attr(x, "threshold")

  data.frame(
    metric = metric_name,
    lower = thresholds[1],
    upper = thresholds[2]
  )
}

thresholds <- bind_rows(
  get_thresholds(qc_counts, "nCount_RNA"),
  get_thresholds(qc_features, "nFeature_RNA"),
  get_thresholds(qc_mito, "percent.mt")
)

write.csv(
  thresholds,
  file.path(
    project_dir,
    "results",
    "qc",
    "QC_thresholds.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Filtering
# ------------------------------------------------------------

combined <- subset(
  combined,
  subset = discard == FALSE
)

combined$discard <- NULL
combined$discard_counts <- NULL
combined$discard_features <- NULL
combined$discard_mito <- NULL

# ------------------------------------------------------------
# QC after filtering
# ------------------------------------------------------------

qc_after <- combined@meta.data %>%
  count(sample, name = "cells_after")

write.csv(
  qc_after,
  file.path(
    project_dir,
    "results",
    "qc",
    "cells_after_QC.csv"
  ),
  row.names = FALSE
)

p_after <- VlnPlot(
  combined,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "sample",
  ncol = 3,
  pt.size = 0.05
)

ggsave(
  file.path(
    project_dir,
    "figures",
    "qc",
    "QC_after_filtering.png"
  ),
  p_after,
  width = 15,
  height = 6,
  dpi = 300
)

saveRDS(
  combined,
  file.path(
    project_dir,
    "results",
    "objects",
    "02_QC_filtered.rds"
  )
)

save_session_info(project_dir)

message(
  "QC completed. Cells retained: ",
  ncol(combined)
)
