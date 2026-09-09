# ============================================================
# 04_integration.R
# Seurat v5 dataset integration
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------

sample_dirs <- list.dirs(
  file.path(data_dir, "samples"),
  recursive = FALSE,
  full.names = TRUE
)

if (length(sample_dirs) < 2) {

  stop(
    "At least two sample directories are required.\n",
    "Expected structure:\n",
    "data/samples/sample1/\n",
    "data/samples/sample2/"
  )

}

# ------------------------------------------------------------
# Read samples
# ------------------------------------------------------------

sample_list <- list()

for (sample_dir in sample_dirs) {

  sample_name <- basename(sample_dir)

  message(
    "Reading sample: ",
    sample_name
  )

  counts <- Read10X(
    data.dir = sample_dir
  )

  if (is.list(counts)) {

    if ("Gene Expression" %in% names(counts)) {

      counts <- counts[["Gene Expression"]]

    } else {

      counts <- counts[[1]]

    }

  }

  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )

  obj$sample <- sample_name

  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )

  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 300 &
      nFeature_RNA < 2500 &
      percent.mt < 10
  )

  obj <- NormalizeData(
    obj,
    verbose = FALSE
  )

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )

  sample_list[[sample_name]] <- obj
}

# ------------------------------------------------------------
# Select integration features
# ------------------------------------------------------------

integration_features <- SelectIntegrationFeatures(
  object.list = sample_list,
  nfeatures = 3000
)

# ------------------------------------------------------------
# Prepare for integration
# ------------------------------------------------------------

sample_list <- lapply(
  sample_list,
  function(x) {

    x <- ScaleData(
      x,
      features = integration_features,
      verbose = FALSE
    )

    x <- RunPCA(
      x,
      features = integration_features,
      npcs = 30,
      verbose = FALSE
    )

    return(x)
  }
)

# ------------------------------------------------------------
# Seurat v5 integration
# ------------------------------------------------------------

seurat_obj <- IntegrateLayers(
  object = merge(
    x = sample_list[[1]],
    y = sample_list[-1]
  ),
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# ------------------------------------------------------------
# Integrated neighbors
# ------------------------------------------------------------

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "integrated.cca",
  dims = 1:30,
  verbose = FALSE
)

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.4,
  verbose = FALSE
)

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "integrated.cca",
  dims = 1:30,
  reduction.name = "umap.integrated",
  verbose = FALSE
)

# ------------------------------------------------------------
# Visualization
# ------------------------------------------------------------

p_sample <- DimPlot(
  seurat_obj,
  reduction = "umap.integrated",
  group.by = "sample"
)

ggsave(
  file.path(
    figures_dir,
    "integration",
    "umap_by_sample.png"
  ),
  p_sample,
  width = 8,
  height = 6,
  dpi = 300
)

p_cluster <- DimPlot(
  seurat_obj,
  reduction = "umap.integrated",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(
    figures_dir,
    "integration",
    "umap_integrated_clusters.png"
  ),
  p_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file.path(
    results_dir,
    "04_integrated_seurat_object.rds"
  )
)

message("Seurat v5 integration completed.")
