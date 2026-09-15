# ============================================================
# functions.R
# Reusable functions for the AML scRNA-seq portfolio
# ============================================================

get_project_dir <- function() {

  project_dir <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )

  if (!dir.exists(file.path(project_dir, "scripts"))) {
    stop(
      "Project root not detected.\n",
      "Run the scripts from the repository root or open ",
      "the repository as an RStudio project."
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


copy_10x_file <- function(
  source_file,
  destination_file
) {

  if (is.na(source_file) || !file.exists(source_file)) {
    stop(
      "10X source file not found: ",
      source_file
    )
  }

  source_is_gz <- grepl(
    "\\.gz$",
    source_file,
    ignore.case = TRUE
  )

  destination_is_gz <- grepl(
    "\\.gz$",
    destination_file,
    ignore.case = TRUE
  )

  if (source_is_gz && destination_is_gz) {

    ok <- file.copy(
      source_file,
      destination_file,
      overwrite = TRUE
    )

  } else if (!source_is_gz && !destination_is_gz) {

    ok <- file.copy(
      source_file,
      destination_file,
      overwrite = TRUE
    )

  } else if (!source_is_gz && destination_is_gz) {

    con_in <- file(
      source_file,
      open = "rb"
    )

    con_out <- gzfile(
      destination_file,
      open = "wb"
    )

    ok <- FALSE

    tryCatch(
      {
        repeat {
          chunk <- readBin(
            con_in,
            what = "raw",
            n = 1024^2
          )

          if (length(chunk) == 0) {
            break
          }

          writeBin(
            chunk,
            con_out
          )
        }

        ok <- TRUE
      },
      finally = {
        close(con_in)
        close(con_out)
      }
    )

  } else {

    stop(
      "Unsupported 10X compression conversion: ",
      source_file,
      " -> ",
      destination_file
    )
  }

  if (!isTRUE(ok)) {
    stop(
      "Failed to copy 10X file: ",
      source_file
    )
  }

  invisible(destination_file)
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

  if (any(is.na(
    c(
      matrix_file,
      barcode_file,
      feature_file
    )
  ))) {
    stop(
      "Incomplete 10X files in: ",
      sample_dir
    )
  }

  tmp <- tempfile("10x_")
  dir.create(tmp)

  on.exit(
    unlink(
      tmp,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  copy_10x_file(
    matrix_file,
    file.path(
      tmp,
      "matrix.mtx.gz"
    )
  )

  copy_10x_file(
    barcode_file,
    file.path(
      tmp,
      "barcodes.tsv.gz"
    )
  )

  copy_10x_file(
    feature_file,
    file.path(
      tmp,
      "features.tsv.gz"
    )
  )

  counts <- Seurat::Read10X(
    data.dir = tmp,
    gene.column = 2
  )

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
