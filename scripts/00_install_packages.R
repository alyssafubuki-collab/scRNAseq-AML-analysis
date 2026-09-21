# ============================================================
# AML scRNA-seq analysis
# GitHub Actions package installation
#
# R 4.4.3
# No GitHub PAT
#
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

message(
  "R version: ",
  R.version.string
)

message("Library paths:")
print(.libPaths())


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

      message(
        "Installing CRAN package: ",
        pkg
      )

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
# 1. Dependencies required by Seurat
#
# Explicit order:
#
# sass
#   ↓
# bslib
#   ↓
# rmarkdown
#   ↓
# htmlwidgets
#   ↓
# plotly
#
# shiny
#   ↓
# miniUI
# ============================================================

message("============================================================")
message("1/8 - Installing Seurat web/reporting dependencies")
message("============================================================")

# ============================================================
# Seurat prerequisites
# ============================================================

cat("\n============================================================\n")
cat("Installing Seurat prerequisites\n")
cat("============================================================\n")

# fs is required by sass
install_cran(c("fs"))

# sass is required by bslib
install_cran(c("sass"))

# bslib is required by rmarkdown
install_cran(c("bslib"))

# Remaining rmarkdown dependencies
install_cran(c(
  "fontawesome",
  "jquerylib",
  "knitr"
))

# rmarkdown
install_cran(c("rmarkdown"))

# htmlwidgets / plotly
install_cran(c(
  "htmlwidgets",
  "plotly"
))

# shiny / miniUI
install_cran(c(
  "shiny",
  "miniUI"
))
# ------------------------------------------------------------
# rmarkdown dependencies
# ------------------------------------------------------------

install_cran(
  c(
    "fontawesome",
    "jquerylib",
    "knitr",
    "rmarkdown"
  )
)

# ------------------------------------------------------------
# htmlwidgets / plotly
# ------------------------------------------------------------

install_cran(
  c(
    "htmlwidgets",
    "plotly"
  )
)

# ------------------------------------------------------------
# shiny / miniUI
# ------------------------------------------------------------

install_cran(
  c(
    "shiny",
    "miniUI"
  )
)


# ============================================================
# 2. Core CRAN packages
# ============================================================

message("============================================================")
message("2/8 - Installing core CRAN packages")
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
message("3/8 - Installing Seurat")
message("============================================================")

if (requireNamespace("Seurat", quietly = TRUE)) {

  message(
    "Seurat already installed: ",
    as.character(packageVersion("Seurat"))
  )

} else {

  message(
    "Seurat is not installed."
  )

  message(
    "Installing Seurat explicitly..."
  )

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
# 4. Bioconductor packages
# ============================================================

message("============================================================")
message("4/8 - Installing Bioconductor packages")
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
message("5/8 - Installing scPred dependencies")
message("============================================================")

scpred_dependencies <- c(
  "ggbeeswarm",
  "MLmetrics",
  "caret",
  "kernlab",
  "pROC"
)

install_cran(
  scpred_dependencies
)


# ============================================================
# 6. Harmony 1.2.4
# ============================================================

message("============================================================")
message("6/8 - Installing Harmony 1.2.4")
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
    paste0(
      "Wrong Harmony version. Expected ",
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
    "HarmonyMatrix is unavailable."
  )
}

message(
  "OK: harmony::HarmonyMatrix available"
)


# ============================================================
# 7. scPred 1.9.2
# ============================================================

message("============================================================")
message("7/8 - Installing scPred")
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
    url = scpred_url,
    destfile = scpred_tar,
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

  message(
    "Installing scPred from local source..."
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
# 8. Final verification
# ============================================================

message("============================================================")
message("8/8 - Final package verification")
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
  "sass",
  "bslib",
  "rmarkdown",
  "fontawesome",
  "jquerylib",
  "knitr",
  "htmlwidgets",
  "plotly",
  "shiny",
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
# scPred compatibility
# ============================================================

if (requireNamespace("scPred", quietly = TRUE)) {

  message(
    "OK: scPred ",
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

  message(
    "FAILED PACKAGES:"
  )

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
