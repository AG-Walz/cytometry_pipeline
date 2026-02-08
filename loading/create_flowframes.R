#!/usr/bin/env Rscript
# =============================================================================
# loading/create_flowframes.R — Create FlowFrames from Expression Matrices
# =============================================================================
# Refactored from: R_scripts/01_loading/Create_ffr.R
# Utilities moved to R_shared/flowframe_utils.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Create flowFrames from expression matrices
#'
#' @param expr_list Named list of expression matrices (or data.frames)
#' @param markers_of_interest Character vector of marker names
#' @param save_fcs Logical: save as FCS files
#' @param output_dir Output directory for FCS files
#' @return List of flowFrame objects
create_flowframes_from_matrices <- function(expr_list,
                                             markers_of_interest = NULL,
                                             save_fcs = FALSE,
                                             output_dir = NULL) {
  flow_frames <- list()

  for (nm in names(expr_list)) {
    cat("Creating flowFrame for:", nm, "...")
    mtx <- as.matrix(expr_list[[nm]])

    # Filter to markers of interest if specified
    if (!is.null(markers_of_interest)) {
      matched <- intersect(colnames(mtx), markers_of_interest)
      if (length(matched) > 0) mtx <- mtx[, matched, drop = FALSE]
    }

    # Generate cell IDs
    cell_ids <- paste0(nm, "_cell_", seq_len(nrow(mtx)))
    rownames(mtx) <- cell_ids

    # Convert to flowFrame
    ff <- from_mtx_to_ffr(mtx)
    flow_frames[[nm]] <- ff
    cat(" ", nrow(mtx), "cells x", ncol(mtx), "markers\n")

    # Save as FCS if requested
    if (save_fcs && !is.null(output_dir)) {
      make_outdir(output_dir)
      fcs_path <- file.path(output_dir, paste0(nm, ".fcs"))
      flowCore::write.FCS(ff, fcs_path)
      cat("  Saved:", fcs_path, "\n")
    }
  }

  flow_frames
}

#' Create a merged flowSet from individual flowFrames
#'
#' @param flow_frames Named list of flowFrame objects
#' @return flowSet object
create_merged_flowset <- function(flow_frames) {
  cat("Creating flowSet from", length(flow_frames), "flowFrames...\n")
  from_ffl_to_fSet(flow_frames, sample_names = names(flow_frames))
}
