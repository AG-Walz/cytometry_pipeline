#!/usr/bin/env Rscript
# =============================================================================
# harmonization/batch_correction.R — Batch Correction with Harmony
# =============================================================================
# Refactored from: R_scripts/03_normalization/Batch_correction.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Run Harmony batch correction on concatenated expression data
#'
#' @param flow_frames Named list of flowFrame objects (must be homogenized)
#' @param sample_names Character vector of sample names
#' @param batch_labels Character vector mapping each sample to a batch
#' @param markers_of_interest Character vector of markers to correct
#' @param do_pca Run PCA before Harmony (default FALSE for expression data)
#' @param ... Additional arguments to harmony::HarmonyMatrix
#' @return List with corrected_matrix, sample_labels, batch_labels
run_batch_correction <- function(flow_frames, sample_names, batch_labels,
                                  markers_of_interest = NULL,
                                  do_pca = FALSE, ...) {
  require(harmony)

  # Concatenate expression matrices
  expr_list <- list()
  sample_vec <- c()
  batch_vec <- c()

  for (i in seq_along(flow_frames)) {
    ff <- flow_frames[[i]]
    expr <- flowCore::exprs(ff)
    if (!is.null(markers_of_interest)) {
      marker_data <- explore_data(ff)
      cols <- which(marker_data$markers %in% markers_of_interest)
      expr <- expr[, cols, drop = FALSE]
      colnames(expr) <- marker_data$markers[cols]
    }
    expr_list[[i]] <- expr
    sample_vec <- c(sample_vec, rep(sample_names[i], nrow(expr)))
    batch_vec <- c(batch_vec, rep(batch_labels[i], nrow(expr)))
  }

  cat("Concatenated", length(expr_list), "samples (",
      sum(sapply(expr_list, nrow)), "total cells)\n")

  combined <- do.call(rbind, expr_list)
  meta <- data.frame(batch = batch_vec, sample = sample_vec)

  cat("Running Harmony batch correction...\n")
  corrected <- harmony::HarmonyMatrix(
    combined,
    meta_data = meta,
    vars_use = "batch",
    do_pca = do_pca,
    ...
  )
  cat("Batch correction complete.\n")

  list(
    corrected_matrix = corrected,
    sample_labels = sample_vec,
    batch_labels = batch_vec,
    markers = colnames(combined)
  )
}
