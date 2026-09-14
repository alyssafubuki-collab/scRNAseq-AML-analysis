# ============================================================
# AML scRNA-seq analysis
# Package installation for GitHub Actions
#
# R version: 4.4.3
#
# No scAnnoX
# No GitHub PAT required
# Harmony 1.2.4
# scPred 1.9.2
# ============================================================

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

options(
  download.file.method = "libcurl"
)

options(
  warn = 1
)

message("============================================================")
message("AML scRNA-seq package installation")
message("============================================================")

message("R version: ", R.version.string)
message("Library paths:")
print(.libPaths())

# ============================================================
# Helper functions
# ============================================================

install_cran <- function(packages) {

  for (pkg in packages) {

    if (requireNamespace(pkg, quietly = TRUE)) {

      message(
        "Already installed: ",
        pkg,
        " ",
        as.character(packageVersion(pkg))
      )

    } else {

      message("Installing CRAN package: ", pkg)

      install.packages(
        pkg,
        dependencies = c(
          "Depends",
          "Imports",
          "LinkingTo"
        )
      )

      if (!requireNamespace(pkg, quietly = TRUE)) {

        stop(
          "Installation failed for CRAN package: ",
          pkg
        )
      }

      message(
        "Successfully installed: ",
        pkg,
        " ",
        as.character(packageVersion(pkg))
      )
    }
  }
}


install_bioc <- function(packages) {

  for (pkg in packages) {

    if (requireNamespace(pkg, quietly = TRUE)) {

      message(
        "Already installed: ",
        pkg,
        " ",
        as.character(packageVersion(pkg))
      )

    } else {

      message("Installing Bioconductor package: ", pkg)

      BiocManager::install(
        pkg,
        ask = FALSE,
        update = FALSE
      )

      if (!requireNamespace(pkg, quietly = TRUE)) {

        stop(
          "Installation failed for Bioconductor package: ",
          pkg
        )
      }

      message(
        "Successfully installed: ",
        pkg,
        " ",
        as.character(packageVersion(pkg))
      )
    }
  }
}


# ============================================================
# 1. Basic CRAN dependencies
# ============================================================

message("============================================================")
message("1/7 - Installing basic CRAN packages")
message("============================================================")

basic_cran <- c(
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork",
  "Matrix",
  "data.table",
  "future",
  "remotes",
  "RhpcBLASctl"
)

install_cran(basic_cran)


# ============================================================
# 2. Seurat
#
# Install separately so that a Seurat failure is isolated.
# ============================================================

message("============================================================")
message("2/7 - Installing Seurat")
message("============================================================")

if (requireNamespace("Seurat", quietly = TRUE)) {

  message(
    "Seurat already installed: ",
    as.character(packageVersion("Seurat"))
  )

} else {

  message("Seurat is not installed.")
  message("Installing Seurat separately...")

  install.packages(
    "Seurat",
    dependencies = c(
      "Depends",
      "Imports",
      "LinkingTo"
    )
  )
}

if (!requireNamespace("Seurat", quietly = TRUE)) {

  stop(
    paste(
      "Seurat installation failed.",
      "See the installation messages immediately above."
    )
  )
}

message(
  "OK: Seurat ",
  as.character(packageVersion("Seurat"))
)


# ============================================================
# 3. Bioconductor packages
# ============================================================

message("============================================================")
message("3/7 - Installing Bioconductor packages")
message("============================================================")

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages(
    "BiocManager",
    dependencies = c(
      "Depends",
      "Imports",
      "LinkingTo"
    )
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  stop(
    "BiocManager installation failed."
  )
}

message(
  "BiocManager ",
  as.character(packageVersion("BiocManager"))
)

bioc_packages <- c(
  "GEOquery",
  "scuttle",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

install_bioc(bioc_packages)


# ============================================================
# 4. scPred dependencies
# ============================================================

message("============================================================")
message("4/7 - Installing scPred dependencies")
message("============================================================")

scpred_dependencies <- c(
  "ggbeeswarm",
  "MLmetrics",
  "caret",
  "kernlab",
  "pROC"
)

install_cran(scpred_dependencies)


# ============================================================
# 5. Harmony 1.2.4
#
# Harmony 1.2.4 is used because the workflow requires
# HarmonyMatrix.
# ============================================================

message("============================================================")
message("5/7 - Installing Harmony 1.2.4")
message("============================================================")

harmony_required <- "1.2.4"

harmony_is_correct <- FALSE

if (requireNamespace("harmony", quietly = TRUE)) {

  current_harmony <- as.character(
    packageVersion("harmony")
  )

  message(
    "Current Harmony version: ",
    current_harmony
  )

  harmony_is_correct <- (
    current_harmony == harmony_required &&
    "HarmonyMatrix" %in%
      getNamespaceExports("harmony")
  )
}

if (!harmony_is_correct) {

  message(
    "Installing Harmony ",
    harmony_required,
    " from CRAN archive."
  )

  harmony_url <- paste0(
    "https://cran.r-project.org/src/contrib/Archive/",
    "harmony/harmony_1.2.4.tar.gz"
  )

  harmony_file <- tempfile(
    pattern = "harmony_",
    fileext = ".tar.gz"
  )

  download_status <- download.file(
    url = harmony_url,
    destfile = harmony_file,
    mode = "wb"
  )

  if (!identical(download_status, 0L)) {

    stop(
      "Failed to download Harmony 1.2.4."
    )
  }

  install_status <- system2(
    "R",
    args = c(
      "CMD",
      "INSTALL",
      "--no-multiarch",
      "--with-keep.source",
      shQuote(harmony_file)
    )
  )

  unlink(harmony_file)

  if (!identical(install_status, 0L)) {

    stop(
      "Harmony 1.2.4 installation failed."
    )
  }
}

# Verify Harmony

if (!requireNamespace("harmony", quietly = TRUE)) {

  stop(
    "Harmony is not installed after installation attempt."
  )
}

harmony_version <- as.character(
  packageVersion("harmony")
)

message(
  "Installed Harmony version: ",
  harmony_version
)

if (harmony_version != harmony_required) {

  stop(
    paste0(
      "Incorrect Harmony version. Expected ",
      harmony_required,
      " but found ",
      harmony_version,
      "."
    )
  )
}

if (
  !"HarmonyMatrix" %in%
  getNamespaceExports("harmony")
) {

  stop(
    "HarmonyMatrix is not available in Harmony 1.2.4."
  )
}

message(
  "OK: harmony::HarmonyMatrix available"
)


# ============================================================
# 6. scPred 1.9.2
#
# Direct GitHub archive installation.
# This avoids remotes::install_github() and GitHub PAT.
# ============================================================

message("============================================================")
message("6/7 - Installing scPred")
message("============================================================")

if (requireNamespace("scPred", quietly = TRUE)) {

  message(
    "scPred already installed: ",
    as.character(packageVersion("scPred"))
  )

} else {

  message(
    "Downloading scPred from public GitHub repository."
  )

  scpred_url <- paste0(
    "https://github.com/powellgenomicslab/scPred/",
    "archive/refs/heads/master.tar.gz"
  )

  scpred_tar <- tempfile(
    pattern = "scPred_",
    fileext = ".tar.gz"
  )

  scpred_source <- tempfile(
    pattern = "scPred_source_"
  )

  dir.create(
    scpred_source,
    recursive = TRUE
  )

  # ----------------------------------------------------------
  # Download
  # ----------------------------------------------------------

  download_status <- download.file(
    url = scpred_url,
    destfile = scpred_tar,
    mode = "wb"
  )

  if (!identical(download_status, 0L)) {

    stop(
      "Failed to download scPred from GitHub."
    )
  }

  # ----------------------------------------------------------
  # Extract
  # ----------------------------------------------------------

  extract_status <- system2(
    "tar",
    args = c(
      "-xzf",
      shQuote(scpred_tar),
      "--strip-components=1",
      "-C",
      shQuote(scpred_source)
    )
  )

  if (!identical(extract_status, 0L)) {

    stop(
      "Failed to extract scPred source archive."
    )
  }

  # ----------------------------------------------------------
  # Install
  # ----------------------------------------------------------

  message(
    "Installing scPred from local source."
  )

  install_status <- system2(
    "R",
    args = c(
      "CMD",
      "INSTALL",
      "--no-multiarch",
      "--with-keep.source",
      shQuote(scpred_source)
    )
  )

  if (!identical(install_status, 0L)) {

    stop(
      "scPred installation failed."
    )
  }

  unlink(scpred_tar)
  unlink(
    scpred_source,
    recursive = TRUE
  )
}

if (!requireNamespace("scPred", quietly = TRUE)) {

  stop(
    "scPred is not installed."
  )
}

message(
  "OK: scPred ",
  as.character(packageVersion("scPred"))
)


# ============================================================
# 7. Final verification
# ============================================================

message("============================================================")
message("7/7 - Final package verification")
message("============================================================")

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork",
  "Matrix",
  "data.table",
  "future",
  "GEOquery",
  "remotes",
  "BiocManager",
  "scuttle",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "RhpcBLASctl",
  "ggbeeswarm",
  "MLmetrics",
  "caret",
  "kernlab",
  "pROC",
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
      as.character(packageVersion(pkg))
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
# Compatibility checks
# ============================================================

message("============================================================")
message("Compatibility checks")
message("============================================================")

# HarmonyMatrix

if (
  requireNamespace("harmony", quietly = TRUE) &&
  "HarmonyMatrix" %in%
    getNamespaceExports("harmony")
) {

  message(
    "OK: harmony::HarmonyMatrix available"
  )

} else {

  failed <- unique(
    c(
      failed,
      "harmony::HarmonyMatrix"
    )
  )

  message(
    "FAILED: harmony::HarmonyMatrix"
  )
}


# scPred

if (requireNamespace("scPred", quietly = TRUE)) {

  message(
    "OK: scPred installed: ",
    as.character(packageVersion("scPred"))
  )

} else {

  failed <- unique(
    c(
      failed,
      "scPred"
    )
  )
}


# ============================================================
# Final result
# ============================================================

message("============================================================")
message("Installation summary")
message("============================================================")

if (length(failed) > 0) {

  message("FAILED PACKAGES:")

  for (pkg in unique(failed)) {

    message(
      " - ",
      pkg
    )
  }

  stop(
    "Package installation failed: ",
    paste(
      unique(failed),
      collapse = ", "
    )
  )
}

message(
  "ALL REQUIRED PACKAGES INSTALLED SUCCESSFULLY."
)

message("============================================================")
message("Installation complete")
message("============================================================")
