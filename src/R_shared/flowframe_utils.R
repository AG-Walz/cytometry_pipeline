# =============================================================================
# R_shared/flowframe_utils.R — FlowFrame & FlowSet Utility Functions
# =============================================================================
# Deduplicated from: Functions.R, Create_ffr.R, data_loading_functions.R
# =============================================================================

#' Explore flow cytometry data from a flowFrame
#' @param flowframe A flowFrame object
#' @return List with dyes, markers, and mcols (marker column indices)
explore_data <- function(flowframe) {
  marker_cols <- c()
  dye_names <- character(ncol(flowCore::exprs(flowframe)))
  marker_names <- character(ncol(flowCore::exprs(flowframe)))

  for (i in seq_len(ncol(flowCore::exprs(flowframe)))) {
    dye_names[i] <- toString(flowCore::featureNames(flowframe[, i]))
    marker_names[i] <- toString(flowCore::markernames(flowframe[, i]))
  }
  for (j in seq_along(marker_names)) {
    if (marker_names[j] != "") marker_cols <- c(marker_cols, j)
  }
  list(dyes = dye_names, markers = marker_names, mcols = marker_cols)
}

#' Compare panels across a list of flowFrames
#' @param flowframes_list List of flowFrame objects
#' @return List of unique marker sets (each element = markers vector)
compare_panels <- function(flowframes_list) {
  marker_dict <- list()
  marker_dict[[1]] <- explore_data(flowframes_list[[1]])$markers
  cat("\nComparing panels.\n")

  for (i in 2:length(flowframes_list)) {
    current_markers <- explore_data(flowframes_list[[i]])$markers
    marker_found <- FALSE
    for (j in seq_along(marker_dict)) {
      if (identical(marker_dict[[j]], current_markers)) {
        marker_dict[[j]] <- c(marker_dict[[j]], i)
        marker_found <- TRUE
        break
      }
    }
    if (!marker_found) {
      marker_dict[[length(marker_dict) + 1]] <- current_markers
    }
  }
  marker_dict
}

#' Extract expression matrices from a list of flowFrames
#'
#' Output to: build_combined_expr_matrix() (this file)
#'
#' @param flow_frames_list List of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param markers_of_interest Character vector of marker names to extract
#' @return Named list of expression matrices (one per sample)
extract_exprmtx_list <- function(flow_frames_list, sample_names, markers_of_interest) {
  exprmtx_list <- list()
  for (i in seq_along(flow_frames_list)) {
    ff <- flow_frames_list[[i]]
    expr_data <- flowCore::exprs(ff)
    marker_data <- explore_data(ff)

    # Find columns matching markers of interest
    matched_cols <- which(marker_data$markers %in% markers_of_interest)
    if (length(matched_cols) > 0) {
      mtx <- expr_data[, matched_cols, drop = FALSE]
      colnames(mtx) <- marker_data$markers[matched_cols]
      exprmtx_list[[sample_names[i]]] <- mtx
    } else {
      warning("No matching markers found for sample: ", sample_names[i])
    }
  }
  exprmtx_list
}

#' Build a combined expression matrix with sample labels from flowFrames
#'
#' Calls extract_exprmtx_list() to get per-sample matrices, row-binds them,
#' and creates a parallel sample-label vector tracking which sample each row
#' belongs to.
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             downsample_flowframes() (preprocessing/downsampling.R),
#'             homogenize_flowframes() (harmonization/panel_homogenization.R)
#' Output to:  run_leiden_with_umap() (clustering/leiden.R),
#'             clusters_to_gating_csv() (clustering/leiden.R),
#'             plot_umap() (visualization/umap.R)
#'
#' @param flow_frames_list List of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param markers_of_interest Character vector of marker names to extract
#' @return List with components:
#'   \item{combined_matrix}{Numeric matrix (cells x markers)}
#'   \item{sample_labels}{Character vector of sample labels, one per row}
build_combined_expr_matrix <- function(flow_frames_list, sample_names, markers_of_interest) {
  expr_list <- extract_exprmtx_list(flow_frames_list, sample_names, markers_of_interest)
  combined_matrix <- do.call(rbind, expr_list)
  sample_labels <- rep(names(expr_list), sapply(expr_list, nrow))
  cat("Combined matrix:", nrow(combined_matrix), "cells x",
      ncol(combined_matrix), "markers from", length(expr_list), "samples\n")
  list(combined_matrix = combined_matrix, sample_labels = sample_labels)
}

#' Generate cell IDs for a flowFrame
#' @param flowframe A flowFrame object
#' @param sample_name Sample name prefix
#' @return Character vector of cell IDs
generate_cellID <- function(flowframe, sample_name) {
  n_cells <- nrow(flowCore::exprs(flowframe))
  paste0(sample_name, "_cell_", seq_len(n_cells))
}

#' Convert expression matrix to flowFrame
#' @param mtx Numeric matrix (rows = cells, columns = markers)
#' @param marker_names Character vector of marker names
#' @return flowFrame object
from_mtx_to_ffr <- function(mtx, marker_names = NULL) {
  if (!is.null(marker_names)) colnames(mtx) <- marker_names

  # Create minimal annotation
  params <- Biobase::AnnotatedDataFrame(
    data.frame(
      name = colnames(mtx),
      desc = colnames(mtx),
      range = apply(mtx, 2, function(x) diff(range(x, na.rm = TRUE))),
      minRange = apply(mtx, 2, min, na.rm = TRUE),
      maxRange = apply(mtx, 2, max, na.rm = TRUE),
      row.names = paste0("$P", seq_len(ncol(mtx))),
      stringsAsFactors = FALSE
    )
  )
  flowCore::flowFrame(exprs = mtx, parameters = params)
}

#' Convert a list of flowFrames to a flowSet
#' @param flow_frame_list List of flowFrame objects
#' @param sample_names Optional character vector of sample names
#' @return flowSet object
from_ffl_to_fSet <- function(flow_frame_list, sample_names = NULL) {
  if (!is.null(sample_names)) {
    names(flow_frame_list) <- sample_names
  }
  flowCore::flowSet(flow_frame_list)
}
