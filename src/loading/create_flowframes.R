#!/usr/bin/env Rscript
# =============================================================================
# loading/create_flowframes.R — Create FlowFrames from Expression Matrices
# =============================================================================
# Refactored from: R_scripts/01_loading/Create_ffr.R
# Utilities moved to R_shared/flowframe_utils.R
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

#' Create flowFrames from expression matrices
#'
#' Bridge function: converts matrix/data.frame data into flowFrame objects,
#' enabling CSV expression data to flow into flowFrame-based scripts
#' (e.g., plot_umap(), plot_density(), preprocessing).
#'
#' Input from: load_csv_expression() (loading/load_csv_expression.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R)
#' Output to:  plot_umap() (visualization/umap.R),
#'             plot_density() (visualization/density_plots.R),
#'             homogenize_flowframes() (harmonization/panel_homogenization.R),
#'             downsample_flowframes() (preprocessing/downsampling.R)
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

  # Validate all frames have same number of columns
  ncols <- sapply(flow_frames, ncol)
  if (length(unique(ncols)) > 1) {
    stop("Cannot merge flowFrames with different numbers of columns. ",
         "Column counts: ", paste(unique(ncols), collapse = ", "),
         ". Use homogenize_flowframes() first.")
  }

  from_ffl_to_fSet(flow_frames, sample_names = names(flow_frames))
}
