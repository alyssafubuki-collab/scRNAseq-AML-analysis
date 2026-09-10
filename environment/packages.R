# ============================================================
# packages.R
# Package information
# ============================================================

packages <- c(
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "dplyr",
  "patchwork",
  "Matrix"
)

installed <- installed.packages()

package_information <- data.frame(
  package = packages,
  version = sapply(
    packages,
    function(x) {
      if (x %in% rownames(installed)) {
        as.character(
          packageVersion(x)
        )
      } else {
        NA
      }
    }
  )
)

print(package_information)

write.csv(
  package_information,
  "environment/package_versions.csv",
  row.names = FALSE
)
