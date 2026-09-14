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
# Helper: CRAN installation
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
  "remotes",
  "RhpcBLASctl"
)

install_cran(cran_packages)

# ============================================================
# Bioconductor
# ============================================================

message("============================================================")
message("Installing Bioconductor packages")
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
# Explicit Seurat + GEOquery verification/reinstallation
# ============================================================

message("============================================================")
message("Ensuring Seurat and GEOquery are installed")
message("============================================================")

# Seurat
if (!requireNamespace("Seurat", quietly = TRUE)) {

  message("Seurat is missing. Installing Seurat explicitly.")

  install.packages(
    "Seurat",
    dependencies = c(
      "Depends",
      "Imports",
      "LinkingTo"
    )
  )
}

# GEOquery
if (!requireNamespace("GEOquery", quietly = TRUE)) {

  message("GEOquery is missing. Installing GEOquery explicitly.")

  BiocManager::install(
    "GEOquery",
    ask = FALSE,
    update = FALSE
  )
}

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
# Harmony 1.2.4
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

  download_status <- download.file(
    url = harmony_url,
    destfile = harmony_tar,
    mode = "wb",
    quiet = FALSE
  )

  if (!identical(download_status, 0L)) {

    stop(
      "Failed to download Harmony ",
      harmony_version_required
    )
  }

  install_status <- system2(
    "R",
    args = c(
      "CMD",
      "INSTALL",
      "--no-multiarch",
      "--with-keep.source",
      harmony_tar
    )
  )

  if (!identical(install_status, 0L)) {

    stop(
      "Harmony installation failed."
    )
  }

  unlink(harmony_tar)
}

# ============================================================
# Verify Harmony
# ============================================================

message("============================================================")
message("Checking Harmony")
message("============================================================")

if (!requireNamespace("harmony", quietly = TRUE)) {

  stop(
    "Harmony installation failed."
  )
}

harmony_version <- as.character(
  packageVersion("harmony")
)

message(
  "Harmony version: ",
  harmony_version
)

if (harmony_version != harmony_version_required) {

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
    "HarmonyMatrix is not available."
  )
}

message(
  "OK: harmony::HarmonyMatrix available"
)

# ============================================================
# scPred
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

  # Download
  download_status <- download.file(
    url = scpred_url,
    destfile = scpred_tar,
    mode = "wb",
    quiet = FALSE
  )

  if (!identical(download_status, 0L)) {

    stop(
      "Failed to download scPred."
    )
  }

  # Extract
  dir.create(
    scpred_src,
    recursive = TRUE,
    showWarnings = FALSE
  )

  extract_status <- system2(
    "tar",
    args = c(
      "-xzf",
      scpred_tar,
      "--strip-components=1",
      "-C",
      scpred_src
    )
  )

  if (!identical(extract_status, 0L)) {

    stop(
      "Failed to extract scPred source."
    )
  }

  # Install
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
      "scPred installation failed."
    )
  }

  unlink(scpred_tar)

  unlink(
    scpred_src,
    recursive = TRUE
  )
}

# ============================================================
# Final verification
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
# Specific compatibility checks
# ============================================================

message("============================================================")
message("Running compatibility checks")
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

  failed_packages <- unique(
    c(
      failed_packages,
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

  failed_packages <- unique(
    c(
      failed_packages,
      "scPred"
    )
  )

  message(
    "FAILED: scPred"
  )
}

# ============================================================
# Final result
# ============================================================

message("============================================================")
message("Installation summary")
message("============================================================")

if (length(failed_packages) > 0) {

  message("FAILED PACKAGES:")

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
