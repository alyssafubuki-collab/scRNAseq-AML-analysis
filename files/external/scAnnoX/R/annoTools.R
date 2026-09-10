#------------------------------------------------------------------------------------------
# Reference-based auto-annotation tools.

#' singleRAnno

#' Automated single-cell annotation using SingleR (Aran et al., 2019).
#' @param obj.seu Query Seurat object, which needs to be annotated.
#' @param ref.obj Reference Seurat object.
#' @param ref.ctype Please specify information about the cell type in the reference object, otherwise: defaults to NULL.
#' i.e. Idents are extracted set to the cell type.
#' @param ... More arguments can be assessed using the SingleR function in the SingleR package.
#' @return Annotated Seurat object, predicted results are embedded into meta.data slot.
#' @export singleRAnno
#'
#' @examples
#' obj.seu <- system.file("data/test", "test.obj.rds", package = "Biotools") %>% readRDS(.)
#' ref.obj <- system.file("data/test", "ref.obj.rds", package = "Biotools") %>% readRDS(.)
#' obj.seu <- singleRAnno(obj.seu, ref.obj)
singleRAnno <- function(obj.seu, ref.obj, ref.ctype = NULL, ...) {
    test.expr <- GetAssayData(obj.seu, layer = "data") %>% as.matrix()
    ref.expr <- GetAssayData(ref.obj, layer = "data") %>% as.matrix()

    if (!is.null(ref.ctype)) {
        label.ref <- ref.obj@meta.data[, ref.ctype]
    } else {
        label.ref <- Idents(ref.obj)
    }
    res.pred <- SingleR::SingleR(test = test.expr, ref = ref.expr, labels = label.ref, ...)
    obj.seu$SingleR <- res.pred$labels
    return(obj.seu)
}