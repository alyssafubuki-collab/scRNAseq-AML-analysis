options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

options(
  download.file.method = "libcurl"
)

message("============================================================")
message("Starting package installation")
message("============================================================")

# ============================================================
# Helper
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

    }
  }
}

# ============================================================
# Core CRAN packages
# ============================================================

message("============================================================")
message("Installing core CRAN packages")
message("============================================================")

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork",
  "Matrix",
  "GEOquery",
  "data.table",
  "future",
  "remotes"
)

install_cran(cran_packages)

# ============================================================
# Bioconductor
# ============================================================

message("============================================================")
message("Installing Bioconductor")
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

# ============================================================
# scPred dependencies
# ============================================================

message("============================================================")
message("Installing scPred dependencies")
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
# Verify scPred dependencies
# ============================================================

message("============================================================")
message("Checking scPred dependencies")
message("============================================================")

required_scpred_dependencies <- c(
  "ggbeeswarm",
  "MLmetrics",
  "caret",
  "kernlab",
  "pROC"
)

for (pkg in required_scpred_dependencies) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    stop(
      "Required scPred dependency is missing: ",
      pkg
    )

  } else {

    message(
      "OK: ",
      pkg,
      " ",
      as.character(packageVersion(pkg))
    )
  }
}

# ============================================================
# Harmony
#
# scPred requires HarmonyMatrix in the workflow.
# Use the known-compatible Harmony 1.2.4 release.
# ============================================================

message("============================================================")
message("Installing/checking Harmony")
message("============================================================")

harmony_version_required <- "1.2.4"

harmony_ok <- FALSE

if (requireNamespace("harmony", quietly = TRUE)) {

  installed_harmony <- as.character(
    packageVersion("harmony")
  )

  message(
    "Harmony currently installed: ",
    installed_harmony
  )

  harmony_ok <- (
    installed_harmony == harmony_version_required &&
    "HarmonyMatrix" %in%
      getNamespaceExports("harmony")
  )
}

if (!harmony_ok) {

  message(
    "Installing Harmony ",
    harmony_version_required
  )

  harmony_url <- paste0(
    "https://cran.r-project.org/src/contrib/Archive/harmony/",
    "harmony_1.2.4.tar.gz"
  )

  harmony_tar <- "harmony_1.2.4.tar.gz"

  unlink(harmony_tar)

  download_ok <- download.file(
    url = harmony_url,
    destfile = harmony_tar,
    mode = "wb",
    quiet = FALSE
  )

  if (!identical(download_ok, 0L)) {

    stop(
      "Failed to download Harmony ",
      harmony_version_required
    )
  }

  system2(
    "R",
    args = c(
      "CMD",
      "INSTALL",
      "--no-multiarch",
      "--with-keep.source",
      harmony_tar
    )
  )

  unlink(harmony_tar)

} else {

  message(
    "Harmony ",
    harmony_version_required,
    " with HarmonyMatrix already installed."
  )
}

# ============================================================
# Verify Harmony
# ============================================================

message("============================================================")
message("Checking Harmony installation")
message("============================================================")

if (!requireNamespace("harmony", quietly = TRUE)) {

  stop("Harmony installation failed.")

}

harmony_version <- as.character(
  packageVersion("harmony")
)

message(
  "Harmony version: ",
  harmony_version
)

if (
  harmony_version != harmony_version_required
) {

  stop(
    "Wrong Harmony version. Expected ",
    harmony_version_required,
    " but found ",
    harmony_version
  )
}

if (
  !"HarmonyMatrix" %in%
  getNamespaceExports("harmony")
) {

  stop(
    "HarmonyMatrix is not exported by the installed Harmony package."
  )
}

message("OK: HarmonyMatrix is available.")

# ============================================================
# scPred
#
# Install directly from the public GitHub archive.
# This deliberately avoids remotes::install_github()
# and therefore avoids GitHub PAT authentication.
# ============================================================

message("============================================================")
message("Installing scPred")
message("============================================================")

if (requireNamespace("scPred", quietly = TRUE)) {

  message(
    "scPred already installed: ",
    as.character(packageVersion("scPred"))
  )

} else {

  message(
    "Downloading scPred from public GitHub repository"
  )

  scpred_url <- paste0(
    "https://github.com/powellgenomicslab/scPred/",
    "archive/refs/heads/master.tar.gz"
  )

  scpred_tar <- "scPred.tar.gz"

  scpred_src <- "scPred-src"

  unlink(scpred_tar)

  unlink(
    scpred_src,
    recursive = TRUE
  )

  # ----------------------------------------------------------
  # Download
  # ----------------------------------------------------------

  download_ok <- download.file(
    url = scpred_url,
    destfile = scpred_tar,
    mode = "wb",
    quiet = FALSE
  )

  if (!identical(download_ok, 0L)) {

    stop(
      "Failed to download scPred from GitHub."
    )
  }

  # ----------------------------------------------------------
  # Extract source
  # ----------------------------------------------------------

  dir.create(
    scpred_src,
    recursive = TRUE,
    showWarnings = FALSE
  )

  tar_status <- system2(
    "tar",
    args = c(
      "-xzf",
      scpred_tar,
      "--strip-components=1",
      "-C",
      scpred_src
    )
  )

  if (!identical(tar_status, 0L)) {

    stop(
      "Failed to extract scPred source archive."
    )
  }

  # ----------------------------------------------------------
  # Install from local source
  # ----------------------------------------------------------

  install_status <- system2(
    "R",
    args = c(
      "CMD",
      "INSTALL",
      "--no-multiarch",
      "--with-keep.source",
      scpred_src
    )
  )

  if (!identical(install_status, 0L)) {

    stop(
      "R CMD INSTALL failed for scPred."
    )
  }

  # ----------------------------------------------------------
  # Cleanup
  # ----------------------------------------------------------

  unlink(
    scpred_tar
  )

  unlink(
    scpred_src,
    recursive = TRUE
  )
}

# ============================================================
# Final package verification
# ============================================================

message("============================================================")
message("Checking package installation")
message("============================================================")

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork",
  "Matrix",
  "GEOquery",
  "remotes",
  "BiocManager",
  "scuttle",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "harmony",
  "ggbeeswarm",
  "MLmetrics",
  "caret",
  "kernlab",
  "pROC",
  "scPred"
)

failed_packages <- character(0)

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

    failed_packages <- c(
      failed_packages,
      pkg
    )
  }
}

# ============================================================
# Specific checks
# ============================================================

message("============================================================")
message("Running specific compatibility checks")
message("============================================================")

# Harmony
if (
  requireNamespace("harmony", quietly = TRUE) &&
  "HarmonyMatrix" %in%
    getNamespaceExports("harmony")
) {

  message(
    "OK: harmony::HarmonyMatrix available"
  )

} else {

  failed_packages <- unique(
    c(
      failed_packages,
      "harmony::HarmonyMatrix"
    )
  )

  message(
    "FAILED: harmony::HarmonyMatrix unavailable"
  )
}

# scPred
if (requireNamespace("scPred", quietly = TRUE)) {

  message(
    "OK: scPred installed"
  )

  # Check bundled reference dataset
  if (
    "pbmc_1" %in%
    getNamespaceExports("scPred")
  ) {

    message(
      "OK: scPred::pbmc_1 available"
    )

  } else {

    message(
      "WARNING: scPred::pbmc_1 not exported by this version"
    )
  }

} else {

  message(
    "FAILED: scPred not installed"
  )
}

# ============================================================
# Final result
# ============================================================

message("============================================================")
message("Installation summary")
message("============================================================")

if (length(failed_packages) > 0) {

  message(
    "FAILED PACKAGES:"
  )

  for (pkg in failed_packages) {
    message(
      " - ",
      pkg
    )
  }

  stop(
    "Package installation failed: ",
    paste(
      failed_packages,
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
