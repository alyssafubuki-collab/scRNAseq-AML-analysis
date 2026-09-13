```r
# ============================================================
# 00_install_packages.R
#
# Install packages required for the complete AML scRNA-seq
# analysis workflow.
#
# Compatible with:
#   R 4.4.x
#   Seurat v5
#   SingleR
#   scPred
#
# scAnnoX is NOT used.
# ============================================================

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

# ============================================================
# 1. CRAN PACKAGES
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

    message("Installing CRAN package: ", pkg)

    install.packages(
      pkg,
      dependencies = TRUE
    )
  } else {

    message("Already installed: ", pkg)
  }
}

# ============================================================
# 2. BIOCONDUCTOR
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  message("Installing BiocManager...")

  install.packages(
    "BiocManager"
  )
}

# Use the Bioconductor release compatible with R 4.4
bioc_version <- "3.20"

message("============================================================")
message("Configuring Bioconductor ", bioc_version)
message("============================================================")

BiocManager::install(
  version = bioc_version,
  ask = FALSE
)

# ============================================================
# 3. BIOCONDUCTOR PACKAGES
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
#
# Harmony is required by scPred.
#
# IMPORTANT:
# We explicitly disable GitHub credentials so that a stale
# GitHub PAT cannot cause an HTTP 401 authentication error.
#
# The GitHub repository is public.
# ============================================================

message("============================================================")
message("Installing Harmony")
message("============================================================")

if (!requireNamespace("harmony", quietly = TRUE)) {

  Sys.unsetenv("GITHUB_PAT")
  Sys.unsetenv("GITHUB_TOKEN")

  remotes::install_github(
    "immunogenomics/harmony",
    upgrade = "never",
    dependencies = TRUE,
    auth_token = NULL
  )

} else {

  message("Harmony already installed.")
}

# ============================================================
# 5. SCPRED
# ============================================================
#
# scPred is installed from its public GitHub repository.
#
# No personal GitHub PAT is required.
# ============================================================

message("============================================================")
message("Installing scPred")
message("============================================================")

if (!requireNamespace("scPred", quietly = TRUE)) {

  Sys.unsetenv("GITHUB_PAT")
  Sys.unsetenv("GITHUB_TOKEN")

  remotes::install_github(
    "powellgenomicslab/scPred",
    upgrade = "never",
    dependencies = TRUE,
    auth_token = NULL
  )

} else {

  message("scPred already installed.")
}

# ============================================================
# 6. VERIFY INSTALLATIONS
# ============================================================

message("============================================================")
message("Verifying installed packages")
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
      as.character(
        packageVersion(pkg)
      )
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
# 7. STOP IF SOMETHING IS MISSING
# ============================================================

if (length(failed) > 0) {

  stop(
    paste(
      "The following packages failed to install:",
      paste(failed, collapse = ", ")
    )
  )
}

# ============================================================
# 8. FINAL MESSAGE
# ============================================================

message("============================================================")
message("ALL REQUIRED PACKAGES INSTALLED SUCCESSFULLY")
message("============================================================")

message("R version: ", R.version.string)

message(
  "Seurat: ",
  as.character(packageVersion("Seurat"))
)

message(
  "SeuratObject: ",
  as.character(packageVersion("SeuratObject"))
)

message(
  "SingleR: ",
  as.character(packageVersion("SingleR"))
)

message(
  "scPred: ",
  as.character(packageVersion("scPred"))
)

message(
  "Harmony: ",
  as.character(packageVersion("harmony"))
)

message("scAnnoX is NOT used in this workflow.")
```
