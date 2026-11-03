#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/normalization.R — Normalization & Scaling Functions
# =============================================================================
# Extracted from: R_scripts/02_preprocessing/Preprocessing_functions.R
# =============================================================================

#' Z-score scaling of marker columns in a flowFrame
#'
#' @param flowframe A flowFrame object
#' @param marker_indices Integer vector of column indices to scale
#' @return Scaled flowFrame
z_score_scaling <- function(flowframe, marker_indices) {
  tf_list <- flowCore::transformList(
    flowCore::featureNames(flowframe[, marker_indices]), scale
  )
  flowCore::transform(flowframe, tf_list)
}

#' Quantile normalization across marker columns
#'
#' @param flowframe A flowFrame object
#' @param marker_indices Integer vector of column indices to normalize
#' @return Quantile-normalized flowFrame
quantile_norm <- function(flowframe, marker_indices) {
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("Package 'preprocessCore' required for quantile normalization. ",
         "Install with: BiocManager::install('preprocessCore')")
  }
  expr_data <- flowCore::exprs(flowframe)
  expr_data[, marker_indices] <- preprocessCore::normalize.quantiles(
    expr_data[, marker_indices]
  )
  flowCore::exprs(flowframe) <- expr_data
  flowframe
}

#' Normalize expression data using Harmony for batch correction
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param batch_vector Character/factor vector of batch labels
#' @param ... Additional arguments to harmony::HarmonyMatrix
#' @return Batch-corrected expression matrix
# Harmony's matrix entry point moved: HarmonyMatrix() (harmony < 2) was removed
# in harmony 2.x, where RunHarmony() takes a cells x features matrix directly.
# Uses whichever this installation has. do_pca = TRUE (HarmonyMatrix only)
# becomes an explicit PCA before RunHarmony.
.harmony_on_matrix <- function(mat, meta_data, vars_use, do_pca = FALSE, ...) {
  ns <- getNamespaceExports("harmony")
  if ("HarmonyMatrix" %in% ns) {
    return(harmony::HarmonyMatrix(mat, meta_data = meta_data,
                                  vars_use = vars_use, do_pca = do_pca, ...))
  }
  if (isTRUE(do_pca)) {
    mat <- stats::prcomp(mat, center = TRUE, scale. = TRUE)$x
  }
  harmony::RunHarmony(mat, meta_data = meta_data, vars_use = vars_use,
                      verbose = FALSE, ...)
}

harmony_normalize <- function(expr_matrix, batch_vector, ...) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Package 'harmony' required for batch correction. ",
         "Install with: install.packages('harmony')")
  }
  meta <- data.frame(batch = batch_vector)
  corrected <- .harmony_on_matrix(expr_matrix, meta_data = meta,
                                  vars_use = "batch", do_pca = FALSE, ...)
  corrected
}
