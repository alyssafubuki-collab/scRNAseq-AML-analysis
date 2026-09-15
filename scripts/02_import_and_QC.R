# ============================================================
# 02_import_and_QC.R
# Import all 8 GSE145410 samples and perform adaptive QC
#
# Seurat v5 layer policy:
#   - create one Seurat object per sample
#   - merge samples
#   - immediately JoinLayers() once
#   - keep a clean, unsplit RNA assay until integration
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
# IMPORT EACH SAMPLE
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

  # Make cell names globally unique before merging.
  obj <- RenameCells(
    obj,
    add.cell.id = sample_name
  )

  obj$gsm <- gsm
  obj$sample <- sample_name
  obj$treatment <- metadata$treatment[i]
  obj$replicate <- metadata$replicate[i]

  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  objects[[sample_name]] <- obj
}

# ============================================================
# MERGE
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

rm(objects)
gc()

# ============================================================
# NORMALIZE SEURAT v5 LAYER STATE
# ============================================================

message("============================================================")
message("RNA layers immediately after merge")
message("============================================================")

print(Layers(combined[["RNA"]]))

message("Joining RNA layers once after merge...")

combined[["RNA"]] <- JoinLayers(
  combined[["RNA"]]
)

gc()

final_layers <- Layers(combined[["RNA"]])

message("============================================================")
message("RNA layers after JoinLayers()")
message("============================================================")

print(final_layers)

if (length(final_layers) > 2L) {
  stop(
    "RNA assay still contains split layers after JoinLayers():\n",
    paste(final_layers, collapse = "\n")
  )
}

if (any(grepl("SeuratProject", final_layers))) {
  stop(
    "Invalid SeuratProject layer names remain after JoinLayers()."
  )
}

# ============================================================
# SAMPLE METADATA
# ============================================================

combined$sample <- factor(
  combined$sample,
  levels = metadata$sample
)

message("============================================================")
message("Sample information")
message("============================================================")

print(table(combined$sample))

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
# ADAPTIVE QC WITH scuttle::isOutlier()
#
# Thresholds are calculated independently by sample.
# No arbitrary fixed nFeature/nCount/mitochondrial cutoffs.
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
# QC THRESHOLDS
# ============================================================

get_thresholds <- function(x, metric_name) {

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

# ============================================================
# FILTER
# ============================================================

cells_before <- ncol(combined)

combined <- subset(
  combined,
  subset = discard == FALSE
)

cells_after <- ncol(combined)

message("Cells before QC: ", cells_before)
message("Cells after QC: ", cells_after)

# ============================================================
# REMOVE QC FLAGS
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
# FINAL LAYER VALIDATION
# ============================================================

final_layers <- Layers(
  combined[["RNA"]]
)

message("============================================================")
message("Final RNA layers")
message("============================================================")

print(final_layers)

if (length(final_layers) > 2L) {
  stop(
    "Unexpected split RNA layers remain in 02_QC_filtered.rds."
  )
}

if (any(grepl("SeuratProject", final_layers))) {
  stop(
    "Invalid SeuratProject layer names remain."
  )
}

# ============================================================
# SAVE
# ============================================================

output_file <- file.path(
  project_dir,
  "results",
  "objects",
  "02_QC_filtered.rds"
)

saveRDS(
  combined,
  output_file
)

save_session_info(project_dir)

message("============================================================")
message("QC completed successfully")
message("Cells retained: ", ncol(combined))
message("Output: ", output_file)
message("============================================================")
