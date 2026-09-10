# ============================================================
# 01_download_GSE145410.R
# Download and extract public GEO data
# ============================================================

source("R/functions.R")

project_dir <- get_project_dir()
data_dir <- file.path(project_dir, "data")
geo_dir <- file.path(data_dir, "GSE145410")
download_dir <- file.path(geo_dir, "download")
extract_dir <- file.path(geo_dir, "extracted")

dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("GEOquery is not installed. Run scripts/00_install_packages.R first.")
}

library(GEOquery)

geo_tar <- file.path(
  download_dir,
  "GSE145410_RAW.tar"
)

if (!file.exists(geo_tar)) {
  message("Downloading GSE145410 supplementary data...")
  getGEOSuppFiles(
    "GSE145410",
    makeDirectory = FALSE,
    baseDir = download_dir
  )

  downloaded <- list.files(
    download_dir,
    pattern = "GSE145410_RAW\\.tar$",
    full.names = TRUE
  )

  if (length(downloaded) == 0) {
    stop("GSE145410_RAW.tar was not found after GEO download.")
  }

  file.copy(downloaded[1], geo_tar, overwrite = TRUE)
}

message("Extracting GEO archive...")

untar(
  geo_tar,
  exdir = extract_dir
)

# GEO's archive contains the 10X triplets. We reorganize them by GSM.
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

all_files <- list.files(
  extract_dir,
  recursive = TRUE,
  full.names = TRUE
)

for (gsm in gsm_ids) {

  gsm_files <- all_files[
    grepl(gsm, basename(all_files), fixed = TRUE)
  ]

  if (length(gsm_files) == 0) {
    warning("No files found for ", gsm)
    next
  }

  target <- file.path(extract_dir, gsm)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)

  for (f in gsm_files) {
    if (normalizePath(dirname(f), mustWork = FALSE) !=
        normalizePath(target, mustWork = FALSE)) {
      file.copy(
        f,
        file.path(target, basename(f)),
        overwrite = TRUE
      )
    }
  }
}

writeLines(
  c(
    "GSE145410",
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
  ),
  file.path(geo_dir, "dataset_info.txt")
)

message("GSE145410 download/extraction completed.")
