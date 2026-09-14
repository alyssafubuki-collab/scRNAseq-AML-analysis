# ============================================================
# 01_download_GSE145410.R
# Download and extract public GEO data
#
# Dataset:
#   GSE145410
#
# 8 public AML scRNA-seq samples
#
# No private data
# No scAnnoX
# ============================================================

source("R/functions.R")

project_dir <- get_project_dir()

data_dir <- file.path(
  project_dir,
  "data"
)

geo_dir <- file.path(
  data_dir,
  "GSE145410"
)

download_dir <- file.path(
  geo_dir,
  "download"
)

extract_dir <- file.path(
  geo_dir,
  "extracted"
)

dir.create(
  download_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  extract_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Check GEOquery
# ============================================================

if (!requireNamespace("GEOquery", quietly = TRUE)) {

  stop(
    paste(
      "GEOquery is not installed.",
      "Run scripts/00_install_packages.R first."
    )
  )
}

library(GEOquery)


# ============================================================
# GEO archive
# ============================================================

geo_tar <- file.path(
  download_dir,
  "GSE145410_RAW.tar"
)


# ============================================================
# Download GEO supplementary files
# ============================================================

if (!file.exists(geo_tar)) {

  message("============================================================")
  message("Downloading GSE145410 supplementary data")
  message("============================================================")

  getGEOSuppFiles(
    "GSE145410",
    makeDirectory = FALSE,
    baseDir = download_dir
  )

  downloaded <- list.files(
    download_dir,
    pattern = "^GSE145410_RAW\\.tar$",
    full.names = TRUE
  )

  if (length(downloaded) == 0) {

    stop(
      "GSE145410_RAW.tar was not found after GEO download."
    )
  }

  source_file <- normalizePath(
    downloaded[1],
    mustWork = FALSE
  )

  target_file <- normalizePath(
    geo_tar,
    mustWork = FALSE
  )

  # ----------------------------------------------------------
  # Avoid copying a file onto itself
  # ----------------------------------------------------------

  if (source_file != target_file) {

    message(
      "Moving downloaded GEO archive to target location..."
    )

    success <- file.copy(
      source_file,
      target_file,
      overwrite = TRUE
    )

    if (!success) {

      stop(
        "Failed to copy GSE145410_RAW.tar."
      )
    }

  } else {

    message(
      "GEO archive already located at target path."
    )
  }

} else {

  message(
    "GSE145410_RAW.tar already exists."
  )
}


# ============================================================
# Validate archive
# ============================================================

if (!file.exists(geo_tar)) {

  stop(
    "GSE145410_RAW.tar does not exist."
  )
}

archive_size <- file.info(
  geo_tar
)$size

message(
  "GSE145410_RAW.tar size: ",
  round(archive_size / 1024^2, 1),
  " MB"
)


# ============================================================
# GEO samples
# ============================================================

gsm_ids <- c(
  "GSM4317809",
  "GSM4317810",
  "GSM4317811",
  "GSM4317812",
  "GSM4317813",
  "GSM4317814",
  "GSM4317815",
  "GSM4317816"
)


# ============================================================
# Extract archive
# ============================================================

message("============================================================")
message("Extracting GEO archive")
message("============================================================")

untar(
  geo_tar,
  exdir = extract_dir
)


# ============================================================
# Find extracted files
# ============================================================

all_files <- list.files(
  extract_dir,
  recursive = TRUE,
  full.names = TRUE
)

if (length(all_files) == 0) {

  stop(
    "No files were extracted from GSE145410_RAW.tar."
  )
}

message(
  "Number of extracted files: ",
  length(all_files)
)


# ============================================================
# Reorganize files by GSM
# ============================================================

message("============================================================")
message("Organizing files by GSM")
message("============================================================")

for (gsm in gsm_ids) {

  message(
    "Processing ",
    gsm
  )

  gsm_files <- all_files[
    grepl(
      gsm,
      basename(all_files),
      fixed = TRUE
    )
  ]

  if (length(gsm_files) == 0) {

    warning(
      "No files found for ",
      gsm
    )

    next
  }

  target <- file.path(
    extract_dir,
    gsm
  )

  dir.create(
    target,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (f in gsm_files) {

    source_path <- normalizePath(
      f,
      mustWork = FALSE
    )

    target_path <- normalizePath(
      file.path(
        target,
        basename(f)
      ),
      mustWork = FALSE
    )

    # Do not copy a file onto itself
    if (source_path == target_path) {
      next
    }

    file.copy(
      f,
      target_path,
      overwrite = TRUE
    )
  }
}


# ============================================================
# Validate GSM directories
# ============================================================

message("============================================================")
message("Validating GSM files")
message("============================================================")

for (gsm in gsm_ids) {

  gsm_dir <- file.path(
    extract_dir,
    gsm
  )

  if (!dir.exists(gsm_dir)) {

    warning(
      "Missing GSM directory: ",
      gsm
    )

    next
  }

  gsm_files <- list.files(
    gsm_dir,
    full.names = TRUE
  )

  message(
    gsm,
    ": ",
    length(gsm_files),
    " file(s)"
  )
}


# ============================================================
# Dataset metadata
# ============================================================

dataset_info <- c(
  "GSE145410",
  "",
  "Public GEO AML single-cell RNA-seq dataset",
  "",
  "GEO samples:",
  "GSM4317809 DMSO A",
  "GSM4317810 DMSO B",
  "GSM4317811 INCB059872 A",
  "GSM4317812 INCB059872 B",
  "GSM4317813 AZA A",
  "GSM4317814 AZA B",
  "GSM4317815 INCB059872+AZA A",
  "GSM4317816 INCB059872+AZA B"
)

writeLines(
  dataset_info,
  file.path(
    geo_dir,
    "dataset_info.txt"
  )
)


# ============================================================
# Final validation
# ============================================================

message("============================================================")
message("Final GEO download validation")
message("============================================================")

if (!file.exists(geo_tar)) {

  stop(
    "Final validation failed: GEO archive is missing."
  )
}

existing_gsm <- gsm_ids[
  dir.exists(
    file.path(
      extract_dir,
      gsm_ids
    )
  )
]

message(
  "GSM directories detected: ",
  length(existing_gsm),
  "/",
  length(gsm_ids)
)

if (length(existing_gsm) != length(gsm_ids)) {

  warning(
    "Not all expected GSM directories were detected."
  )
}

message("============================================================")
message("GSE145410 download/extraction completed successfully.")
message("============================================================")
