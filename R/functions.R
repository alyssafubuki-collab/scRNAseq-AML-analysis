# ============================================================
# functions.R
# Reusable functions
# ============================================================

get_project_dir <- function() {

  # Works when scripts are launched from the repository root.
  # If launched from an RStudio project, getwd() is normally the root.

  project_dir <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )

  if (!dir.exists(file.path(project_dir, "scripts"))) {
    stop(
      "Project root not detected.\n",
      "Run the scripts from the repository root or open the repository as an RStudio project."
    )
  }

  project_dir
}


make_project_dirs <- function(project_dir) {

  dirs <- c(
    "figures/qc",
    "figures/integration",
    "figures/annotation",
    "figures/differential_expression",
    "results/qc",
    "results/annotation",
    "results/differential_expression",
    "results/objects",
    "environment"
  )

  invisible(
    lapply(
      dirs,
      function(x) {
        dir.create(
          file.path(project_dir, x),
          recursive = TRUE,
          showWarnings = FALSE
        )
      }
    )
  )
}


get_gse145410_metadata <- function() {

  data.frame(
    gsm = c(
      "GSM4317809",
      "GSM4317810",
      "GSM4317811",
      "GSM4317812",
      "GSM4317813",
      "GSM4317814",
      "GSM4317815",
      "GSM4317816"
    ),
    sample = c(
      "DMSO_A",
      "DMSO_B",
      "INCB059872_A",
      "INCB059872_B",
      "AZA_A",
      "AZA_B",
      "INCB059872_AZA_A",
      "INCB059872_AZA_B"
    ),
    treatment = c(
      "DMSO",
      "DMSO",
      "INCB059872",
      "INCB059872",
      "AZA",
      "AZA",
      "INCB059872_AZA",
      "INCB059872_AZA"
    ),
    replicate = c(
      "A", "B", "A", "B", "A", "B", "A", "B"
    ),
    stringsAsFactors = FALSE
  )
}


find_10x_file <- function(sample_dir, pattern) {

  x <- list.files(
    sample_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = pattern
  )

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[1]
}


read_gse145410_sample <- function(sample_dir) {

  matrix_file <- find_10x_file(
    sample_dir,
    "matrix\\.mtx(\\.gz)?$"
  )

  barcode_file <- find_10x_file(
    sample_dir,
    "barcodes\\.tsv(\\.gz)?$"
  )

  feature_file <- find_10x_file(
    sample_dir,
    "(features|genes)\\.tsv(\\.gz)?$"
  )

  if (any(is.na(c(matrix_file, barcode_file, feature_file)))) {
    stop(
      "Incomplete 10X files in: ",
      sample_dir
    )
  }

  # Read10X expects a directory containing the standard files.
  # We create a temporary standardized directory.
  tmp <- tempfile("10x_")
  dir.create(tmp)

  file.copy(matrix_file, file.path(tmp, "matrix.mtx.gz"))
  file.copy(barcode_file, file.path(tmp, "barcodes.tsv.gz"))
  file.copy(feature_file, file.path(tmp, "features.tsv.gz"))

  counts <- Seurat::Read10X(
    data.dir = tmp,
    gene.column = 2
  )

  unlink(tmp, recursive = TRUE)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  counts
}


save_session_info <- function(project_dir) {

  writeLines(
    capture.output(sessionInfo()),
    file.path(
      project_dir,
      "environment",
      "sessionInfo.txt"
    )
  )
}
