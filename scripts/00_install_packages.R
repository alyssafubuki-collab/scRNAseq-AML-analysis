# ============================================================
# AML scRNA-seq analysis
# GitHub Actions package installation
#
# R 4.4.3
# No scAnnoX
# No GitHub PAT
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

# ============================================================
# Helper
# ============================================================

install_cran <- function(packages) {

  for (pkg in packages) {

    if (requireNamespace(pkg, quietly = TRUE)) {

      message(
        "OK already installed: ",
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
          "Installation failed for package: ",
          pkg
        )
      }

      message(
        "OK installed: ",
        pkg,
        " ",
        as.character(packageVersion(pkg))
      )
    }
  }
}

# ============================================================
# 1. Seurat prerequisites
#
# Install problematic Seurat dependencies explicitly first.
# ============================================================

message("============================================================")
message("1/7 - Installing Seurat prerequisites")
message("============================================================")

seurat_prerequisites <- c(
  "rmarkdown",
  "shiny",
  "htmlwidgets",
  "plotly",
  "miniUI"
)

install_cran(seurat_prerequisites)

# ============================================================
# 2. Core CRAN packages
# ============================================================

message("============================================================")
message("2/7 - Installing core CRAN packages")
message("============================================================")

core_cran <- c(
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

install_cran(core_cran)

# ============================================================
# 3. Seurat
# ============================================================

message("============================================================")
message("3/7 - Installing Seurat")
message("============================================================")

if (requireNamespace("Seurat", quietly = TRUE)) {

  message(
    "Seurat already installed: ",
    as.character(packageVersion("Seurat"))
  )

} else {

  message("Seurat is not installed.")
  message("Installing Seurat explicitly...")

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
    "Seurat installation failed."
  )
}

message(
  "OK: Seurat ",
  as.character(packageVersion("Seurat"))
)

# ============================================================
# 4. Bioconductor
# ============================================================

message("============================================================")
message("4/7 - Installing Bioconductor packages")
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

bioc_packages <- c(
  "GEOquery",
  "scuttle",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

for (pkg in bioc_packages) {

  if (requireNamespace(pkg, quietly = TRUE)) {

    message(
      "OK already installed: ",
      pkg,
      " ",
      as.character(packageVersion(pkg))
    )

  } else {

    message(
      "Installing Bioconductor package: ",
      pkg
    )

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )

    if (!requireNamespace(pkg, quietly = TRUE)) {

      stop(
        "Bioconductor installation failed: ",
        pkg
      )
    }

    message(
      "OK installed: ",
      pkg,
      " ",
      as.character(packageVersion(pkg))
    )
  }
}

# ============================================================
# 5. scPred dependencies
# ============================================================

message("============================================================")
message("5/7 - Installing scPred dependencies")
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
# 6. Harmony 1.2.4
# ============================================================

message("============================================================")
message("6/7 - Installing Harmony 1.2.4")
message("============================================================")

harmony_required <- "1.2.4"

harmony_ok <- FALSE

if (requireNamespace("harmony", quietly = TRUE)) {

  current_harmony <- as.character(
    packageVersion("harmony")
  )

  message(
    "Current Harmony version: ",
    current_harmony
  )

  harmony_ok <- (
    current_harmony == harmony_required &&
    "HarmonyMatrix" %in%
      getNamespaceExports("harmony")
  )
}

if (!harmony_ok) {

  message(
    "Installing Harmony ",
    harmony_required
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
    harmony_url,
    harmony_file,
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

if (!requireNamespace("harmony", quietly = TRUE)) {

  stop(
    "Harmony is not installed."
  )
}

harmony_version <- as.character(
  packageVersion("harmony")
)

message(
  "Harmony version: ",
  harmony_version
)

if (harmony_version != harmony_required) {

  stop(
    "Wrong Harmony version."
  )
}

if (
  !"HarmonyMatrix" %in%
  getNamespaceExports("harmony")
) {

  stop(
    "HarmonyMatrix is unavailable."
  )
}

message(
  "OK: harmony::HarmonyMatrix available"
)

# ============================================================
# 7. scPred
# ============================================================

message("============================================================")
message("7/7 - Installing scPred")
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

  download_status <- download.file(
    scpred_url,
    scpred_tar,
    mode = "wb"
  )

  if (!identical(download_status, 0L)) {

    stop(
      "Failed to download scPred."
    )
  }

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
      "Failed to extract scPred."
    )
  }

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

  unlink(scpred_tar)

  unlink(
    scpred_source,
    recursive = TRUE
  )

  if (!identical(install_status, 0L)) {

    stop(
      "scPred installation failed."
    )
  }
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
# Final verification
# ============================================================

message("============================================================")
message("Final package verification")
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
  "rmarkdown",
  "shiny",
  "htmlwidgets",
  "plotly",
  "miniUI",
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
# Harmony compatibility
# ============================================================

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
