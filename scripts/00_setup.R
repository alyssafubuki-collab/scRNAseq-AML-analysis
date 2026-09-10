# ============================================================
# 00_setup.R
# Project setup and package loading
# ============================================================

# Clear workspace
rm(list = ls())

# ------------------------------------------------------------
# Project directories
# ------------------------------------------------------------

project_dir <- normalizePath(
  file.path(getwd()),
  winslash = "/",
  mustWork = FALSE
)

data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(
  file.path(figures_dir, "qc"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(figures_dir, "integration"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(figures_dir, "annotation"),
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  file.path(figures_dir, "differential_expression"),
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork",
  "Matrix"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {

  message(
    "The following packages are missing:\n",
    paste(missing_packages, collapse = ", ")
  )

}

# Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(Matrix)
})

message("Seurat version: ", packageVersion("Seurat"))
message("SeuratObject version: ", packageVersion("SeuratObject"))

# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

set.seed(1234)

# ------------------------------------------------------------
# Session information
# ------------------------------------------------------------

session_file <- file.path(
  project_dir,
  "environment",
  "sessionInfo.txt"
)

dir.create(
  dirname(session_file),
  showWarnings = FALSE,
  recursive = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  con = session_file
)

message("Project successfully initialized.")
