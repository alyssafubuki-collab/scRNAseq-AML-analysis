# ============================================================
# 06_differential_expression.R
# Differential gene expression
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

input_file <- file.path(
  results_dir,
  "05_annotated_seurat_object.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Annotated object not found. ",
    "Run 05_cell_annotation.R first."
  )
}

seurat_obj <- readRDS(input_file)

# ------------------------------------------------------------
# Identify markers for each cluster
# ------------------------------------------------------------

Idents(seurat_obj) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# ------------------------------------------------------------
# Save all markers
# ------------------------------------------------------------

write.csv(
  cluster_markers,
  file.path(
    results_dir,
    "cluster_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Top markers
# ------------------------------------------------------------

top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

write.csv(
  top_markers,
  file.path(
    results_dir,
    "top_cluster_markers.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Example pairwise comparison
# ------------------------------------------------------------

clusters <- levels(Idents(seurat_obj))

if (length(clusters) >= 2) {

  cluster_1 <- clusters[1]
  cluster_2 <- clusters[2]

  comparison <- FindMarkers(
    seurat_obj,
    ident.1 = cluster_1,
    ident.2 = cluster_2,
    min.pct = 0.1,
    logfc.threshold = 0.25
  )

  comparison$gene <- rownames(comparison)

  write.csv(
    comparison,
    file.path(
      results_dir,
      paste0(
        "DE_cluster_",
        cluster_1,
        "_vs_",
        cluster_2,
        ".csv"
      )
    ),
    row.names = FALSE
  )

}

# ------------------------------------------------------------
# Volcano-like plot
# ------------------------------------------------------------

if (exists("comparison")) {

  comparison <- comparison %>%
    mutate(
      significance = case_when(
        p_val_adj < 0.05 &
          avg_log2FC > 0.5 ~ "Up",
        p_val_adj < 0.05 &
          avg_log2FC < -0.5 ~ "Down",
        TRUE ~ "Not significant"
      )
    )

  p <- ggplot(
    comparison,
    aes(
      x = avg_log2FC,
      y = -log10(p_val_adj + 1e-300)
    )
  ) +
    geom_point(
      alpha = 0.6
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    theme_classic() +
    labs(
      title = paste(
        "Differential expression:",
        cluster_1,
        "vs",
        cluster_2
      ),
      x = "Average log2 fold change",
      y = "-log10 adjusted p-value"
    )

  ggsave(
    file.path(
      figures_dir,
      "differential_expression",
      "DE_comparison.png"
    ),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

message(
  "Differential expression analysis completed."
)
