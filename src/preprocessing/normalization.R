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
harmony_normalize <- function(expr_matrix, batch_vector, ...) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Package 'harmony' required for batch correction. ",
         "Install with: install.packages('harmony')")
  }
  meta <- data.frame(batch = batch_vector)
  corrected <- harmony::HarmonyMatrix(
    expr_matrix,
    meta_data = meta,
    vars_use = "batch",
    do_pca = FALSE,
    ...
  )
  corrected
}
