# ============================================================
# 02_import_and_QC.R
#
# Import GSE145410
# Adaptive QC using scuttle::isOutlier()
# Clean Seurat v5 RNA layers
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

  # Unique cell names before merging
  obj <- RenameCells(
    obj,
    add.cell.id = sample_name
  )

  # Metadata
  obj$gsm <- gsm
  obj$sample <- sample_name
  obj$treatment <- metadata$treatment[i]
  obj$replicate <- metadata$replicate[i]

  # Mitochondrial content
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
# CLEAN SEURAT v5 LAYERS
# ============================================================

message("============================================================")
message("RNA layers after merge")
message("============================================================")

print(
  Layers(combined[["RNA"]])
)

combined[["RNA"]] <- JoinLayers(
  combined[["RNA"]]
)

gc()

message("============================================================")
message("RNA layers after JoinLayers")
message("============================================================")

print(
  Layers(combined[["RNA"]])
)

layer_names <- Layers(
  combined[["RNA"]]
)

if (any(grepl("SeuratProject", layer_names))) {
  stop(
    "Invalid SeuratProject layers detected after JoinLayers:\n",
    paste(layer_names, collapse = "\n")
  )
}

# At this point the object must have one counts layer.
if (!"counts" %in% layer_names) {
  stop(
    "Expected a single 'counts' layer after JoinLayers()."
  )
}

# ============================================================
# SAMPLE METADATA
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
    "results",
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
# ADAPTIVE QC
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
# SAVE QC THRESHOLDS
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

message(
  "Cells before QC: ",
  cells_before
)

message(
  "Cells after QC: ",
  cells_after
)

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
    "results",
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
# FINAL VALIDATION
# ============================================================

final_layers <- Layers(
  combined[["RNA"]]
)

message("============================================================")
message("Final RNA layers")
message("============================================================")

print(final_layers)

if (length(final_layers) != 2 ||
    !all(c("counts", "data") %in% final_layers)) {

  # At this stage data should normally not exist yet.
  # We explicitly require counts only.
  if (!identical(final_layers, "counts")) {
    stop(
      "Unexpected RNA layers after QC: ",
      paste(final_layers, collapse = ", ")
    )
  }
}

# ============================================================
# SAVE
# ============================================================

output <- file.path(
  project_dir,
  "results",
  "objects",
  "02_QC_filtered.rds"
)

saveRDS(
  combined,
  output
)

save_session_info(
  project_dir
)

message("============================================================")
message("QC completed successfully")
message("Cells retained: ", ncol(combined))
message("Output: ", output)
message("============================================================")
