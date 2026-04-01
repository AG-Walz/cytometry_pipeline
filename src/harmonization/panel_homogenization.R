#!/usr/bin/env Rscript
# =============================================================================
# harmonization/panel_homogenization.R — Panel Homogenization
# =============================================================================
# Combined from: R_scripts/03_normalization/Homogenize.R + Unify.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Translate elements of a vector using a dictionary
#' @param vector Character vector
#' @param dictionary Named list (names = original, values = replacement)
#' @return Translated vector
translate_vector <- function(vector, dictionary) {
  sapply(vector, function(x) {
    if (x %in% names(dictionary)) dictionary[[x]] else "Unknown"
  })
}

#' Replace marker names in a flowFrame
#' @param flowframe A flowFrame object
#' @param replace_map Named list (names = old markers, values = new markers)
#' @return flowFrame with updated marker names
replace_markernames <- function(flowframe, replace_map) {
  marker_list <- flowCore::markernames(flowframe)
  marker_pos <- explore_data(flowframe)$mcols
  for (i in seq_along(marker_list)) {
    if (marker_list[i] %in% names(replace_map)) {
      new_marker <- replace_map[[marker_list[i]]]
      cat("Marker", marker_list[i], "->", new_marker, "\n")
      names(new_marker) <- flowCore::featureNames(flowframe)[marker_pos[i]]
      flowCore::markernames(flowframe)[i] <- new_marker
    }
  }
  flowframe
}

#' Reorder expression matrix to match target marker order
#' @param flowframe A flowFrame object
#' @param marker_name_list Character vector of target marker order
#' @return Reordered expression matrix
reorder_exprmtx <- function(flowframe, marker_name_list) {
  data <- explore_data(flowframe)
  result <- matrix(NA, nrow = nrow(flowCore::exprs(flowframe)),
                   ncol = length(marker_name_list))
  colnames(result) <- marker_name_list
  missing_markers <- c()
  for (j in seq_along(marker_name_list)) {
    col_idx <- which(data$markers == marker_name_list[j])
    if (length(col_idx) > 0) {
      result[, j] <- flowCore::exprs(flowframe)[, col_idx[1]]
    } else {
      missing_markers <- c(missing_markers, marker_name_list[j])
    }
  }
  if (length(missing_markers) > 0) {
    warning("Markers filled with NA (not found in flowFrame): ",
            paste(missing_markers, collapse = ", "))
  }
  result
}

#' Homogenize a list of flowFrames to a common panel
#'
#' @param flow_frames Named list of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param markers_of_interest Target marker list
#' @param markername_replacing Named list for marker renaming
#' @param save_csv Save reordered CSVs
#' @param save_fcs Save reordered FCS files
#' @param output_dir Output directory
#' @return Named list of homogenized flowFrame objects
homogenize_flowframes <- function(flow_frames, sample_names,
                                   markers_of_interest,
                                   markername_replacing = list(),
                                   save_csv = FALSE,
                                   save_fcs = FALSE,
                                   output_dir = NULL) {
  cat("Homogenizing", length(flow_frames), "samples to",
      length(markers_of_interest), "markers\n")

  if (!is.null(output_dir)) make_outdir(output_dir)

  hg_frames <- list()
  cell_ids <- generate_cellID_list(flow_frames)

  for (i in seq_along(flow_frames)) {
    sname <- sample_names[i]
    cat("  Homogenizing", sname, "...")

    ff <- flow_frames[[i]]

    # Rename markers if needed
    if (length(markername_replacing) > 0) {
      # Check for collisions: different old names mapping to the same new name
      new_names <- unlist(markername_replacing)
      if (any(duplicated(new_names))) {
        dupes <- unique(new_names[duplicated(new_names)])
        stop("markername_replacing collision: multiple markers map to the same name: ",
             paste(dupes, collapse = ", "))
      }
      ff <- replace_markernames(ff, markername_replacing)
    }

    # Reorder to target panel
    mtx <- reorder_exprmtx(ff, markers_of_interest)

    # Convert back to flowFrame
    hg_ff <- from_mtx_to_ffr(mtx, markers_of_interest)
    hg_frames[[sname]] <- hg_ff

    cat(" ", nrow(mtx), "cells x", ncol(mtx), "markers\n")

    # Save if requested
    if (save_csv && !is.null(output_dir)) {
      write.csv(mtx, file.path(output_dir, paste0("hg_", sname, ".csv")))
    }
    if (save_fcs && !is.null(output_dir)) {
      flowCore::write.FCS(hg_ff, file.path(output_dir, paste0("hg_", sname, ".fcs")))
    }
  }

  hg_frames
}

# Helper: generate cell IDs for a list of flowFrames
generate_cellID_list <- function(flow_frames) {
  total <- 0
  ids <- list()
  for (i in seq_along(flow_frames)) {
    n <- nrow(flow_frames[[i]])
    ids[[i]] <- (total + 1):(total + n)
    total <- total + n
  }
  ids
}
