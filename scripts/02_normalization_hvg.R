# ============================================================
# 02_normalization_hvg.R
# Normalization and highly variable genes
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

input_file <- file.path(
  results_dir,
  "01_qc_seurat_object.rds"
)

if (!file.exists(input_file)) {
  stop("QC object not found. Run 01_quality_control.R first.")
}

seurat_obj <- readRDS(input_file)

# ------------------------------------------------------------
# Normalize
# ------------------------------------------------------------

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

# ------------------------------------------------------------
# Find highly variable genes
# ------------------------------------------------------------

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

# ------------------------------------------------------------
# Plot HVGs
# ------------------------------------------------------------

top10 <- head(
  VariableFeatures(seurat_obj),
  10
)

variable_plot <- VariableFeaturePlot(
  seurat_obj
)

label_plot <- LabelPoints(
  plot = variable_plot,
  points = top10,
  repel = TRUE
)

ggsave(
  filename = file.path(
    figures_dir,
    "qc",
    "highly_variable_genes.png"
  ),
  plot = label_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Scaling
# ------------------------------------------------------------

all_genes <- rownames(seurat_obj)

seurat_obj <- ScaleData(
  seurat_obj,
  features = all_genes,
  verbose = FALSE
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file.path(
    results_dir,
    "02_normalized_seurat_object.rds"
  )
)

write.csv(
  data.frame(
    gene = VariableFeatures(seurat_obj)
  ),
  file.path(
    results_dir,
    "highly_variable_genes.csv"
  ),
  row.names = FALSE
)

message("Normalization and HVG selection completed.")
