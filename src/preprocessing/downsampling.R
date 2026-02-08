#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/downsampling.R — Downsampling Functions
# =============================================================================
# Merged from: R_scripts/02_preprocessing/Downsample.R + Downsampling_functions.R
# =============================================================================

#' Downsample a list of flowFrames
#'
#' @param flow_frames List of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param target_cells Target number of cells
#' @param method "Consistent" (same number per sample) or "Maximal" (proportional)
#' @param seed Random seed
#' @param save_fcs Save downsampled FCS files
#' @param output_dir Output directory for FCS files
#' @return List of downsampled flowFrame objects
downsample_flowframes <- function(flow_frames, sample_names,
                                   target_cells = 10000,
                                   method = "Maximal",
                                   seed = 42,
                                   save_fcs = FALSE,
                                   output_dir = NULL) {
  set.seed(seed)
  n_samples <- length(flow_frames)
  cell_counts <- sapply(flow_frames, nrow)
  total_cells <- sum(cell_counts)

  cat("Downsampling", n_samples, "samples (total:", total_cells, "cells)\n")
  cat("Method:", method, "| Target:", target_cells, "cells\n")

  ds_frames <- list()

  if (method == "Consistent") {
    # Each sample contributes the same number of cells
    per_sample <- min(target_cells %/% n_samples, min(cell_counts))
    cat("Consistent: taking", per_sample, "cells per sample\n")

    for (i in seq_along(flow_frames)) {
      idx <- sort(sample(nrow(flow_frames[[i]]), per_sample))
      ds_frames[[i]] <- flow_frames[[i]][idx, ]
      cat("  ", sample_names[i], ":", nrow(flow_frames[[i]]), "->", per_sample, "\n")
    }
  } else if (method == "Maximal") {
    # Proportional downsampling
    fractions <- cell_counts / total_cells
    per_sample <- pmin(round(fractions * target_cells), cell_counts)

    for (i in seq_along(flow_frames)) {
      n_take <- per_sample[i]
      if (n_take >= nrow(flow_frames[[i]])) {
        ds_frames[[i]] <- flow_frames[[i]]
      } else {
        idx <- sort(sample(nrow(flow_frames[[i]]), n_take))
        ds_frames[[i]] <- flow_frames[[i]][idx, ]
      }
      cat("  ", sample_names[i], ":", nrow(flow_frames[[i]]), "->",
          nrow(ds_frames[[i]]), "\n")
    }
  } else {
    stop("Unknown method: '", method, "'. Use 'Consistent' or 'Maximal'.")
  }

  names(ds_frames) <- sample_names

  if (save_fcs && !is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    for (i in seq_along(ds_frames)) {
      fcs_path <- file.path(output_dir, paste0("ds_", sample_names[i], ".fcs"))
      flowCore::write.FCS(ds_frames[[i]], fcs_path)
    }
    cat("Saved downsampled FCS files to:", output_dir, "\n")
  }

  cat("Total after downsampling:", sum(sapply(ds_frames, nrow)), "cells\n")
  ds_frames
}
