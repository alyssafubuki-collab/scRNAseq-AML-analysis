# ============================================================
# 01_quality_control.R
# Single-cell RNA-seq quality control
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------

# Expected input:
# data/raw_matrix/
#
# For 10X data, the directory should contain:
# barcodes.tsv.gz
# features.tsv.gz
# matrix.mtx.gz

input_dir <- file.path(
  data_dir,
  "raw_matrix"
)

if (!dir.exists(input_dir)) {

  stop(
    "Input directory not found: ",
    input_dir,
    "\nPlease place your public 10X dataset in data/raw_matrix/"
  )

}

# ------------------------------------------------------------
# Read 10X data
# ------------------------------------------------------------

counts <- Read10X(
  data.dir = input_dir
)

# If multiple modalities are returned,
# keep the RNA expression matrix.
if (is.list(counts)) {

  if ("Gene Expression" %in% names(counts)) {

    counts <- counts[["Gene Expression"]]

  } else {

    counts <- counts[[1]]

  }

}

# ------------------------------------------------------------
# Create Seurat object
# ------------------------------------------------------------

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "AML_scRNAseq",
  min.cells = 3,
  min.features = 200
)

# ------------------------------------------------------------
# Mitochondrial genes
# ------------------------------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

qc_summary <- data.frame(
  metric = c(
    "Cells",
    "Genes"
  ),
  value = c(
    ncol(seurat_obj),
    nrow(seurat_obj)
  )
)

write.csv(
  qc_summary,
  file.path(results_dir, "qc_summary_before_filtering.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# QC plots
# ------------------------------------------------------------

p1 <- VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3,
  pt.size = 0.1
)

ggsave(
  filename = file.path(
    figures_dir,
    "qc",
    "qc_before_filtering.png"
  ),
  plot = p1,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Filtering
# ------------------------------------------------------------

library(scuttle)

qc.lib <- isOutlier(
    seurat_obj$nCount_RNA,
    nmads = 3,
    type = "lower",
    log = TRUE
)

qc.features <- isOutlier(
    seurat_obj$nFeature_RNA,
    nmads = 3,
    type = "lower",
    log = TRUE
)

qc.mito <- isOutlier(
    seurat_obj$percent.mt,
    nmads = 3,
    type = "higher"
)

discard <- qc.lib | qc.features | qc.mito

seurat_obj$discard <- discard

seurat_obj <- subset(
    seurat_obj,
    subset = discard == FALSE
)

# ------------------------------------------------------------
# QC after filtering
# ------------------------------------------------------------

qc_summary_after <- data.frame(
  metric = c(
    "Cells after filtering",
    "Genes"
  ),
  value = c(
    ncol(seurat_obj),
    nrow(seurat_obj)
  )
)

write.csv(
  qc_summary_after,
  file.path(results_dir, "qc_summary_after_filtering.csv"),
  row.names = FALSE
)

p2 <- VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3,
  pt.size = 0.1
)

ggsave(
  filename = file.path(
    figures_dir,
    "qc",
    "qc_after_filtering.png"
  ),
  plot = p2,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# Save object
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file.path(
    results_dir,
    "01_qc_seurat_object.rds"
  )
)

message(
  "QC completed: ",
  ncol(seurat_obj),
  " cells retained."
)
