# ============================================================
# 08_annotation_scAnnoX.R
#
# scAnnoX annotation
#
# scAnnoX is installed by GitHub Actions before this script.
# The package source is stored in:
#
# files/external/scAnnoX/
#
# Annotation framework:
#
# SingleR
# Seurat
# sciBet
# scmap
# CHETAH
# scSorter
# sc.type
# cellID
# scCATCH
# SCINA
#
# NOTE:
# scAnnoX is used here as an annotation-comparison framework.
# The independent annotation workflows used elsewhere in this
# repository remain SingleR and scPred.
# ============================================================


# ============================================================
# 1. LOAD PROJECT FUNCTIONS
# ============================================================

source("R/functions.R")


# ============================================================
# 2. LOAD LIBRARIES
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(scAnnoX)
})


# ============================================================
# 3. PROJECT DIRECTORY
# ============================================================

project_dir <- get_project_dir()


# ============================================================
# 4. CHECK scAnnoX INSTALLATION
# ============================================================

if (!requireNamespace(
  "scAnnoX",
  quietly = TRUE
)) {

  stop(
    paste(
      "\nscAnnoX is not installed.\n",
      "The GitHub Actions workflow must install scAnnoX before",
      "running 08_annotation_scAnnoX.R.\n"
    )
  )

}


# ============================================================
# 5. CHECK INPUT OBJECT
# ============================================================

input_file <- file.path(
  project_dir,
  "results",
  "objects",
  "05_integrated.rds"
)


if (!file.exists(input_file)) {

  stop(
    paste(
      "\nIntegrated Seurat object not found:\n",
      input_file,
      "\n\nRun 05_Seurat5_integration.R first.\n"
    )
  )

}


object <- readRDS(input_file)


# ============================================================
# 6. CHECK REQUIRED REDUCTION
# ============================================================

if (!"umap.integrated" %in% names(object@reductions)) {

  stop(
    paste(
      "\nThe Seurat object does not contain",
      "'umap.integrated'.\n",
      "Run 05_Seurat5_integration.R first.\n"
    )
  )

}


# ============================================================
# 7. CHECK OUTPUT DIRECTORIES
# ============================================================

dir.create(
  file.path(
    project_dir,
    "results",
    "annotation"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  file.path(
    project_dir,
    "results",
    "objects"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  file.path(
    project_dir,
    "figures",
    "annotation"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 8. PUBLIC REFERENCE
#
# scPred provides pbmc_1 as a public reference object.
# ============================================================

if (!requireNamespace(
  "scPred",
  quietly = TRUE
)) {

  stop(
    "scPred is required because pbmc_1 is used as the public reference."
  )

}


reference <- scPred::pbmc_1


# ============================================================
# 9. PREPARE REFERENCE
# ============================================================

reference <- NormalizeData(
  reference,
  verbose = FALSE
)


reference <- FindVariableFeatures(
  reference,
  verbose = FALSE
)


reference <- ScaleData(
  reference,
  verbose = FALSE
)


reference <- RunPCA(
  reference,
  npcs = 30,
  verbose = FALSE
)


# ============================================================
# 10. CHECK REFERENCE CELL TYPES
# ============================================================

if (!"cell_type" %in% colnames(reference@meta.data)) {

  stop(
    "The scPred::pbmc_1 reference does not contain 'cell_type'."
  )

}


Idents(reference) <- reference$cell_type


# ============================================================
# 11. BUILD MARKER LIST
# ============================================================

message(
  "Building marker list for scAnnoX..."
)


marker_list <- findMarkerToolsForSc(
  reference,
  to.list = TRUE,
  top.k = 30
)


if (is.null(marker_list)) {

  stop(
    "findMarkerToolsForSc() returned NULL."
  )

}


# ============================================================
# 12. METHODS USED BY scAnnoX
# ============================================================

anno_tools <- c(
  "SingleR",
  "Seurat",
  "sciBet",
  "scmap",
  "CHETAH",
  "scSorter",
  "sc.type",
  "cellID",
  "scCATCH",
  "SCINA"
)


marker_based <- c(
  "scSorter",
  "sc.type",
  "cellID",
  "scCATCH",
  "SCINA"
)


# ============================================================
# 13. CHECK REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "SingleR",
  "Seurat",
  "scibetR",
  "scmap",
  "CHETAH",
  "scSorter",
  "CelliD",
  "scCATCH",
  "SCINA"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0) {

  stop(
    paste(
      "The following scAnnoX runtime packages are missing:",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )

}


# ============================================================
# 14. RUN scAnnoX METHODS
# ============================================================

message("")
message("==============================================")
message("Running scAnnoX annotation comparison")
message("==============================================")
message("")


annotation_results <- list()


for (tool in anno_tools) {

  message("")
  message("----------------------------------------------")
  message(
    "Running scAnnoX method: ",
    tool
  )
  message("----------------------------------------------")


  strategy <- if (
    tool %in% marker_based
  ) {
    "marker-based"
  } else {
    "refernce-based"
  }


  result <- tryCatch(

    {

      autoAnnoTools(
        object,
        ref.obj = reference,
        ref.ctype = "cell_type",
        marker.lst = marker_list,
        method = tool,
        select.marker = "Seurat",
        top.k = 30,
        strategy = strategy
      )

    },

    error = function(e) {

      stop(
        paste(
          "\nscAnnoX failed for method:",
          tool,
          "\n\nError:",
          conditionMessage(e),
          "\n"
        )
      )

    }

  )


  if (is.null(result)) {

    stop(
      paste(
        "scAnnoX returned NULL for method:",
        tool
      )
    )

  }


  if (!inherits(
    result,
    "Seurat"
  )) {

    stop(
      paste(
        "Unexpected result type returned by scAnnoX method:",
        tool
      )
    )

  }


  if (ncol(result@meta.data) < 1) {

    stop(
      paste(
        "No annotation column returned by scAnnoX method:",
        tool
      )
    )

  }


  annotation_results[[tool]] <-
    result@meta.data[
      ,
      ncol(result@meta.data)
    ]


  message(
    "Completed: ",
    tool
  )

}


# ============================================================
# 15. COMBINE RESULTS
# ============================================================

anno_res <- do.call(
  cbind,
  annotation_results
)


anno_res <- as.data.frame(
  anno_res,
  stringsAsFactors = FALSE
)


colnames(anno_res) <- anno_tools


# ============================================================
# 16. ENSURE CELL NAMES ARE PRESERVED
# ============================================================

if (is.null(rownames(anno_res))) {

  stop(
    "scAnnoX annotation results do not contain cell names."
  )

}


# ============================================================
# 17. scAnnoX CONSENSUS RESULT
# ============================================================

message("")
message("Calculating scAnnoX consensus annotation...")


scAnnoX_result <- annoResult(
  anno_res
)


if (is.null(scAnnoX_result)) {

  stop(
    "annoResult() returned NULL."
  )

}


if (!"scAnnoX" %in% colnames(scAnnoX_result)) {

  stop(
    "annoResult() did not return a 'scAnnoX' column."
  )

}


# ============================================================
# 18. ADD FINAL ANNOTATION TO SEURAT OBJECT
# ============================================================

object$scAnnoX_label <-
  scAnnoX_result$scAnnoX


# ============================================================
# 19. SAVE RAW METHOD PREDICTIONS
# ============================================================

prediction_table <- cbind(
  cell = rownames(anno_res),
  anno_res
)


write.csv(
  prediction_table,
  file = file.path(
    project_dir,
    "results",
    "annotation",
    "scAnnoX_method_predictions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. SAVE CONSENSUS PREDICTIONS
# ============================================================

consensus_table <- cbind(
  cell = rownames(scAnnoX_result),
  scAnnoX_result
)


write.csv(
  consensus_table,
  file = file.path(
    project_dir,
    "results",
    "annotation",
    "scAnnoX_predictions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. SAVE ANNOTATED SEURAT OBJECT
# ============================================================

saveRDS(
  object,
  file = file.path(
    project_dir,
    "results",
    "objects",
    "08_scAnnoX.rds"
  )
)


# ============================================================
# 22. GENERATE UMAP
# ============================================================

p <- DimPlot(
  object,
  reduction = "umap.integrated",
  group.by = "scAnnoX_label",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "scAnnoX consensus annotation"
  )


ggsave(
  file.path(
    project_dir,
    "figures",
    "annotation",
    "scAnnoX_UMAP.png"
  ),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)


# ============================================================
# 23. SAVE METHOD COMPARISON HEATMAP DATA
# ============================================================

write.csv(
  prediction_table,
  file = file.path(
    project_dir,
    "results",
    "annotation",
    "scAnnoX_all_methods.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. SUMMARY
# ============================================================

message("")
message("==============================================")
message("scAnnoX annotation completed successfully.")
message("==============================================")
message("")

message(
  "Methods completed: ",
  paste(
    anno_tools,
    collapse = ", "
  )
)

message("")

message(
  "Consensus annotation:",
  paste(
    sort(
      unique(
        object$scAnnoX_label
      )
    ),
    collapse = ", "
  )
)

message("")

message(
  "Output object:",
  file.path(
    project_dir,
    "results",
    "objects",
    "08_scAnnoX.rds"
  )
)

message(
  "Prediction table:",
  file.path(
    project_dir,
    "results",
    "annotation",
    "scAnnoX_predictions.csv"
  )
)

message(
  "UMAP:",
  file.path(
    project_dir,
    "figures",
    "annotation",
    "scAnnoX_UMAP.png"
  )
)

message("")
