#!/usr/bin/env Rscript
# =============================================================================
# harmonization/batch_correction.R — Batch Correction with Harmony
# =============================================================================
# Refactored from: R_scripts/03_normalization/Batch_correction.R
# =============================================================================

# ---- cytometry src bootstrap (Nextflow-safe); see src/R_shared/paths.R ----
if (!exists(".cyto_paths_loaded", inherits = TRUE)) {
  local({
    root <- Sys.getenv("CYTOMETRY_SRC", "")
    if (!nzchar(root) || !dir.exists(file.path(root, "R_shared"))) {
      of <- NULL
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        f <- frames[[i]]$ofile
        if (!is.null(f) && nzchar(f)) { of <- f; break }
      }
      # Rscript src/foo/bar.R carries no source frame; recover from --file=.
      if (is.null(of)) {
        a <- commandArgs(FALSE); m <- grep("^--file=", a)
        if (length(m)) of <- sub("^--file=", "", a[m[[1]]])
      }
      d <- if (!is.null(of)) dirname(normalizePath(of, mustWork = FALSE)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

# Guarded, NOT unconditional: nextflow/bin/*_cli.R and scripts/*.R set shared_dir
# themselves (the --src_dir convention), and clobbering it here would send every
# source() below to whichever tree this file happens to live in rather than the
# one the caller chose.
if (!exists("shared_dir")) shared_dir <- cyto_shared_dir()
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

  # Validate all frames have identical marker columns
  ref_cols <- sort(colnames(flowCore::exprs(flow_frames[[1]])))
  for (i in seq_along(flow_frames)[-1]) {
    curr_cols <- sort(colnames(flowCore::exprs(flow_frames[[i]])))
    if (!identical(ref_cols, curr_cols)) {
      stop("Frame '", names(flow_frames)[i], "' has different markers than frame '",
           names(flow_frames)[1], "'. Use homogenize_flowframes() first.")
    }
  }

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
