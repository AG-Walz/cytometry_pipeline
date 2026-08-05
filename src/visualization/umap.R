#!/usr/bin/env Rscript
# =============================================================================
# visualization/umap.R — UMAP Dimensionality Reduction Plots
# =============================================================================
# Merged from: R_scripts/06_visualization/UMAP.R + UMAP_filt.R
# =============================================================================

# Respect a shared_dir the caller already set — the CLI convention used by
# nextflow/bin/*_cli.R and scripts/*.R — and only guess from the working
# directory when there is none. The unconditional assignment this replaces
# clobbered the caller's value, so sourcing this file from anywhere but the
# subproject directory sent every source() below to a path that does not exist.
# That is not hypothetical: run_dimred() source()s this file, and the FlowJo
# UMAP run died on "/tmp/src/R_shared/io_utils.R: No such file or directory".
# Same guard as loading/load_flowjo_wsp.R.
if (!exists("shared_dir")) {
  shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
}
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Create UMAP plots with marker expression overlays (flowFrame input)
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             downsample_flowframes() (preprocessing/downsampling.R),
#'             homogenize_flowframes() (harmonization/panel_homogenization.R),
#'             create_flowframes_from_matrices() (loading/create_flowframes.R)
#' Output to:  standalone visualization output
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
                      output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                      plot_width = 10, plot_height = 8,
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

#' Create UMAP plots from a pre-built expression matrix
#'
#' Matrix-based alternative to plot_umap() that accepts a numeric matrix
#' directly, without requiring flowFrame objects. This enables CSV expression
#' data and any matrix-producing script to generate UMAP visualizations.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             clusters_to_gating_csv() (clustering/leiden.R) [read back],
#'             run_leiden_with_umap() (clustering/leiden.R) [for cluster overlay]
#' Output to:  standalone visualization output
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param sample_labels Character vector of sample labels (one per row)
#' @param markers_of_interest Character vector of markers to plot (NULL = all columns)
#' @param cluster_labels Optional integer/character vector of cluster assignments
#' @param condition_labels Optional character vector of condition/group labels (one per row).
#'                          When non-NULL, generates a faceted UMAP plot split by condition.
#' @param output_dir Output directory for plots
#' @param n_neighbors UMAP n_neighbors parameter
#' @param min_dist UMAP min_dist parameter
#' @param umap_coords Optional pre-computed UMAP coordinates (2-column matrix).
#'                     If provided, skips UMAP computation.
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @param seed Random seed
#' @return List with umap_coords, sample_labels, plots
plot_umap_from_matrix <- function(expr_matrix, sample_labels = NULL,
                                   markers_of_interest = NULL,
                                   cluster_labels = NULL,
                                   condition_labels = NULL,
                                   output_dir,
                                   n_neighbors = 15,
                                   min_dist = 0.1,
                                   metric = "euclidean",
                                   n_threads = 1L, progress = FALSE,
                                   umap_coords = NULL,
                                   output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                                   plot_width = 10, plot_height = 8,
                                   seed = 42,
                                   draw_plots = TRUE) {
  require(ggplot2)
  # draw_plots = FALSE is the "compute the embedding, I will draw it myself"
  # path: a caller whose figures are its own (the FlowJo population UMAP colours
  # by gating subpopulation, by workspace and by every marker) still wants THIS
  # function's coordinates, so that there is one implementation of the embedding
  # rather than a private copy of run_umap() per caller. No output directory is
  # needed in that case, and creating one would leave an empty folder behind.
  if (draw_plots) make_outdir(output_dir)

  # Filter markers
  exclude_patterns <- c("^$", "Time", "FSC", "SSC", "index")
  all_markers <- colnames(expr_matrix)
  filtered_idx <- which(!grepl(paste(exclude_patterns, collapse = "|"),
                               all_markers, ignore.case = TRUE))
  filtered_markers <- all_markers[filtered_idx]

  if (is.null(markers_of_interest)) {
    markers_of_interest <- filtered_markers
  } else {
    markers_of_interest <- intersect(markers_of_interest, all_markers)
  }

  expr_subset <- expr_matrix[, markers_of_interest, drop = FALSE]

  # Compute or reuse UMAP coordinates
  if (is.null(umap_coords)) {
    cat("Computing UMAP on", nrow(expr_subset), "cells x",
        ncol(expr_subset), "markers...\n")
    umap_coords <- run_umap(expr_subset, n_neighbors = n_neighbors,
                             min_dist = min_dist, metric = metric, seed = seed,
                             n_threads = n_threads, progress = progress)
  } else {
    cat("Using pre-computed UMAP coordinates for", nrow(umap_coords), "cells\n")
  }

  if (is.null(sample_labels)) {
    sample_labels <- rep("all", nrow(expr_matrix))
  }

  if (!draw_plots) {
    return(list(umap_coords = umap_coords, sample_labels = sample_labels,
                plots = list()))
  }

  umap_df <- data.frame(UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2],
                         sample = sample_labels)

  # Plot by sample
  p_sample <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sample)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_minimal() + labs(title = "UMAP by Sample")
  save_plot_formats(p_sample, file.path(output_dir, "umap_sample"),
                    width = plot_width, height = plot_height,
                    formats = output_formats)

  plots <- list(sample = p_sample)

  # Plot by cluster if provided
  if (!is.null(cluster_labels)) {
    umap_df$cluster <- as.factor(cluster_labels)
    p_cluster <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
      geom_point(size = 0.3, alpha = 0.5) +
      theme_minimal() + labs(title = "UMAP by Cluster")
    save_plot_formats(p_cluster, file.path(output_dir, "umap_cluster"),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
    plots[["cluster"]] <- p_cluster
  }

  # Plot faceted by condition if provided
  if (!is.null(condition_labels)) {
    umap_df$condition <- condition_labels
    p_condition <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = sample)) +
      geom_point(size = 0.3, alpha = 0.5) +
      facet_wrap(~ condition) +
      theme_minimal() + labs(title = "UMAP by Condition")
    save_plot_formats(p_condition, file.path(output_dir, "umap_condition"),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
    plots[["condition"]] <- p_condition
  }

  # Plot per marker
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
