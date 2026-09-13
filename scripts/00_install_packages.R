# ============================================================
# 00_install_packages.R
#
# Installation of all packages required for the AML
# scRNA-seq analysis workflow.
#
# R 4.4.x
# Seurat v5
# SingleR
# scPred
#
# scAnnoX is NOT used.
# ============================================================

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

# ============================================================
# 1. CRAN
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

message("============================================================")
message("Installing CRAN packages")
message("============================================================")

for (pkg in cran_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    message("Installing: ", pkg)

    install.packages(
      pkg,
      dependencies = TRUE
    )

  } else {

    message(
      "Already installed: ",
      pkg,
      " ",
      packageVersion(pkg)
    )
  }
}

# ============================================================
# 2. BIOCMANAGER
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages("BiocManager")
}

# ============================================================
# 3. BIOCONDUCTOR
# ============================================================

bioc_packages <- c(
  "scuttle",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

message("============================================================")
message("Installing Bioconductor packages")
message("============================================================")

BiocManager::install(
  bioc_packages,
  ask = FALSE,
  update = FALSE
)

# ============================================================
# 4. HARMONY
# ============================================================

message("============================================================")
message("Checking Harmony")
message("============================================================")

if (!requireNamespace("harmony", quietly = TRUE)) {

  message("Harmony is not installed.")

  stop(
    "Harmony must be installed before running the analysis."
  )

} else {

  message(
    "Harmony already installed: ",
    packageVersion("harmony")
  )
}

# ============================================================
# 5. SCPRED
# ============================================================

message("============================================================")
message("Installing scPred")
message("============================================================")

if (!requireNamespace("scPred", quietly = TRUE)) {

  message("Downloading scPred from public GitHub repository")

  system2(
    "curl",
    args = c(
      "-L",
      "--fail",
      "--retry", "3",
      "-H", "Authorization:",
      "-H", "X-GitHub-Api-Version:",
      "https://github.com/powellgenomicslab/scPred/archive/refs/heads/master.tar.gz",
      "-o", "scPred.tar.gz"
    )
  )

  unlink("scPred-src", recursive = TRUE)

  dir.create(
    "scPred-src",
    recursive = TRUE
  )

  system2(
    "tar",
    args = c(
      "-xzf",
      "scPred.tar.gz",
      "--strip-components=1",
      "-C",
      "scPred-src"
    )
  )

  system2(
    "R",
    args = c(
      "CMD",
      "INSTALL",
      "--no-multiarch",
      "--with-keep.source",
      "scPred-src"
    )
  )

} else {

  message(
    "scPred already installed: ",
    packageVersion("scPred")
  )
}

# ============================================================
# 6. VERIFY
# ============================================================

message("============================================================")
message("Checking package installation")
message("============================================================")

required_packages <- c(
  cran_packages,
  "BiocManager",
  bioc_packages,
  "harmony",
  "scPred"
)

failed <- character(0)

for (pkg in required_packages) {

  if (requireNamespace(pkg, quietly = TRUE)) {

    message(
      "OK: ",
      pkg,
      " ",
      packageVersion(pkg)
    )

  } else {

    message(
      "FAILED: ",
      pkg
    )

    failed <- c(
      failed,
      pkg
    )
  }
}

# ============================================================
# 7. STOP IF A PACKAGE IS MISSING
# ============================================================

if (length(failed) > 0) {

  stop(
    paste(
      "Package installation failed:",
      paste(failed, collapse = ", ")
    )
  )
}

# ============================================================
# 8. FINAL CHECK
# ============================================================

message("============================================================")
message("ALL PACKAGES INSTALLED SUCCESSFULLY")
message("============================================================")

message(
  "R version: ",
  R.version.string
)

message(
  "Seurat: ",
  packageVersion("Seurat")
)

message(
  "SeuratObject: ",
  packageVersion("SeuratObject")
)

message(
  "SingleR: ",
  packageVersion("SingleR")
)

message(
  "scPred: ",
  packageVersion("scPred")
)

message(
  "Harmony: ",
  packageVersion("harmony")
)

message("scAnnoX: NOT USED")

