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

# ============================================================
# IMPORT EACH SAMPLE SEPARATELY
# ============================================================

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

  message("============================================================")
  message("Reading ", gsm, " - ", sample_name)
  message("============================================================")

  counts <- read_gse145410_sample(sample_dir)

  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 0
  )

  # ----------------------------------------------------------
  # Make cell IDs unique before merging
  # ----------------------------------------------------------

  obj <- RenameCells(
    obj,
    add.cell.id = sample_name
  )

  # ----------------------------------------------------------
  # Metadata
  # ----------------------------------------------------------

  obj$gsm <- gsm
  obj$sample <- sample_name
  obj$treatment <- metadata$treatment[i]
  obj$replicate <- metadata$replicate[i]

  # ----------------------------------------------------------
  # Mitochondrial percentage
  # ----------------------------------------------------------

  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  objects[[sample_name]] <- obj
}

# ============================================================
# MERGE ALL SAMPLES
# ============================================================

message("============================================================")
message("Merging all samples")
message("============================================================")

combined <- merge(
  x = objects[[1]],
  y = objects[-1],
  merge.data = FALSE,
  project = "GSE145410"
)

# ============================================================
# IMPORTANT SEURAT v5 LAYER CLEANUP
# ============================================================
#
# merge() can create several counts.* layers.
# We immediately collapse them into a single RNA counts layer.
#
# This is intentional:
#
# 02 -> ONE clean RNA layer
# 03 -> normalization, no split
# 05 -> ONE split operation before integration
#
# This prevents repeated:
# .SeuratProject.SeuratProject.SeuratProject...
# layer names.
# ============================================================

message("============================================================")
message("RNA layers immediately after merge")
message("============================================================")

rna_layers <- Layers(combined[["RNA"]])

print(rna_layers)

if (length(rna_layers) > 1) {

  message("Multiple RNA layers detected after merge.")
  message("Joining RNA layers into a single clean assay...")

  combined[["RNA"]] <- JoinLayers(
    combined[["RNA"]]
  )
}

message("============================================================")
message("RNA layers after cleanup")
message("============================================================")

print(
  Layers(combined[["RNA"]])
)

# ============================================================
# CHECK SAMPLE METADATA
# ============================================================

combined$sample <- factor(
  combined$sample
)

message("============================================================")
message("Sample information")
message("============================================================")

print(
  table(combined$sample)
)

# ============================================================
# QC BEFORE FILTERING
# ============================================================

qc_before <- combined@meta.data %>%
  count(
    sample,
    name = "cells_before"
  )

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

# ============================================================
# ADAPTIVE QC WITH isOutlier()
# ============================================================
#
# Thresholds are calculated independently for each sample.
# This preserves the actual adaptive QC strategy.
# ============================================================

message("============================================================")
message("Calculating adaptive QC thresholds")
message("============================================================")

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

# ============================================================
# QC THRESHOLD TABLE
# ============================================================

get_thresholds <- function(
  x,
  metric_name
) {

  thresholds <- attr(
    x,
    "threshold"
  )

  data.frame(
    metric = metric_name,
    lower = thresholds[1],
    upper = thresholds[2]
  )
}

thresholds <- bind_rows(
  get_thresholds(
    qc_counts,
    "nCount_RNA"
  ),
  get_thresholds(
    qc_features,
    "nFeature_RNA"
  ),
  get_thresholds(
    qc_mito,
    "percent.mt"
  )
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

# ============================================================
# FILTER CELLS
# ============================================================

message("============================================================")
message("Filtering low-quality cells")
message("============================================================")

cells_before <- ncol(combined)

combined <- subset(
  combined,
  subset = discard == FALSE
)

cells_after <- ncol(combined)

message(
  "Cells before QC: ",
  cells_before
)

message(
  "Cells after QC: ",
  cells_after
)

# ============================================================
# REMOVE TEMPORARY QC FLAGS
# ============================================================

combined$discard <- NULL
combined$discard_counts <- NULL
combined$discard_features <- NULL
combined$discard_mito <- NULL

# ============================================================
# QC AFTER FILTERING
# ============================================================

qc_after <- combined@meta.data %>%
  count(
    sample,
    name = "cells_after"
  )

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

# ============================================================
# FINAL RNA LAYER VALIDATION
# ============================================================

message("============================================================")
message("Final RNA layer validation")
message("============================================================")

final_layers <- Layers(
  combined[["RNA"]]
)

print(final_layers)

if (length(final_layers) > 3) {
  warning(
    "More than three RNA layers remain after QC: ",
    paste(final_layers, collapse = ", ")
  )
}

# ============================================================
# SAVE
# ============================================================

saveRDS(
  combined,
  file.path(
    project_dir,
    "results",
    "objects",
    "02_QC_filtered.rds"
  )
)

save_session_info(
  project_dir
)

message("============================================================")
message("QC completed successfully")
message("Cells retained: ", ncol(combined))
message("RNA layers: ", paste(final_layers, collapse = ", "))
message("============================================================")
