# ============================================================
# Record installed package versions
# ============================================================

packages <- c(
  "Seurat",
  "SeuratObject",
  "scuttle",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "scPred",
  "harmony",
  "GEOquery",
  "ggplot2",
  "dplyr",
  "patchwork"
)

installed <- installed.packages()

out <- data.frame(
  package = packages,
  version = sapply(
    packages,
    function(x) {
      if (x %in% rownames(installed)) {
        as.character(packageVersion(x))
      } else {
        NA_character_
      }
    }
  ),
  stringsAsFactors = FALSE
)

dir.create(
  "environment",
  showWarnings = FALSE
)

write.csv(
  out,
  "environment/package_versions.csv",
  row.names = FALSE
)

print(out)
