#!/usr/bin/env Rscript
# =============================================================================
# visualization/umap.R — UMAP Dimensionality Reduction Plots
# =============================================================================
# Merged from: R_scripts/06_visualization/UMAP.R + UMAP_filt.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Create UMAP plots with marker expression overlays
#'
#' @param flow_frames Named list of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param markers_of_interest Character vector of markers to plot
#' @param output_dir Output directory for plots
#' @param n_neighbors UMAP n_neighbors parameter
#' @param min_dist UMAP min_dist parameter
#' @param marker_selection_mode "all", "popup", or "terminal"
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @param seed Random seed
#' @return List with umap_coords, plots
plot_umap <- function(flow_frames, sample_names,
                      markers_of_interest = NULL,
                      output_dir,
                      n_neighbors = 15,
                      min_dist = 0.1,
                      marker_selection_mode = "all",
                      output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
                      plot_width = 8, plot_height = 6,
                      seed = 42) {
  require(ggplot2)
  make_outdir(output_dir)

  # Concatenate expression data
  expr_list <- list()
  sample_labels <- c()
  for (i in seq_along(flow_frames)) {
    ff <- flow_frames[[i]]
    data <- explore_data(ff)
    expr <- flowCore::exprs(ff)[, data$mcols, drop = FALSE]
    colnames(expr) <- data$markers[data$mcols]
    expr_list[[i]] <- expr
    sample_labels <- c(sample_labels, rep(sample_names[i], nrow(expr)))
  }
  combined <- do.call(rbind, expr_list)

  # Filter markers
  exclude_patterns <- c("^$", "Time", "FSC", "SSC", "index")
  all_markers <- colnames(combined)
  filtered_idx <- which(!grepl(paste(exclude_patterns, collapse = "|"),
                               all_markers, ignore.case = TRUE))
  filtered_markers <- all_markers[filtered_idx]

  if (is.null(markers_of_interest)) {
    markers_of_interest <- select_markers_by_mode(
      marker_selection_mode, filtered_markers, filtered_idx
    )
    markers_of_interest <- all_markers[markers_of_interest]
  }

  expr_subset <- combined[, markers_of_interest, drop = FALSE]

  # Compute UMAP
  cat("Computing UMAP on", nrow(expr_subset), "cells x",
      ncol(expr_subset), "markers...\n")
  umap_coords <- run_umap(expr_subset, n_neighbors = n_neighbors,
                           min_dist = min_dist, seed = seed)

  umap_df <- data.frame(UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2],
                         sample = sample_labels)

  # Plot by sample
  p_sample <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sample)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_minimal() + labs(title = "UMAP by Sample")
  save_plot_formats(p_sample, file.path(output_dir, "umap_sample"),
                    width = plot_width, height = plot_height,
                    formats = output_formats)

  # Plot per marker
  plots <- list(sample = p_sample)
  for (marker in markers_of_interest) {
    umap_df$expr <- expr_subset[, marker]
    p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = expr)) +
      geom_point(size = 0.3) +
      viridis::scale_color_viridis(name = marker) +
      theme_minimal() + labs(title = paste("UMAP -", marker))
    save_plot_formats(p, file.path(output_dir, paste0("umap_", create_safe_filename(marker))),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
    plots[[marker]] <- p
  }

  cat("UMAP plots saved to:", output_dir, "\n")

  list(umap_coords = umap_coords, sample_labels = sample_labels, plots = plots)
}
