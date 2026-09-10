# ============================================================
# 00_install_packages.R
# Install packages required by the portfolio workflow
# ============================================================

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork",
  "Matrix",
  "GEOquery",
  "remotes"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_packages <- c(
  "scuttle",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

BiocManager::install(
  bioc_packages,
  ask = FALSE,
  update = FALSE
)

# scPred requires Harmony and is installed from GitHub.
if (!requireNamespace("harmony", quietly = TRUE)) {
  remotes::install_github("immunogenomics/harmony")
}

if (!requireNamespace("scPred", quietly = TRUE)) {
  remotes::install_github("powellgenomicslab/scPred")
}

message("Core packages installed.")
message("IMPORTANT: copy your Seurat-v5-compatible scAnnoX package into:")
message("external/scAnnoX/")
