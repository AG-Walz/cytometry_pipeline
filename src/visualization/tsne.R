#!/usr/bin/env Rscript
# =============================================================================
# visualization/tsne.R — t-SNE Dimensionality Reduction Plots
# =============================================================================
# t-SNE visualization with dual input modes (flowFrame and matrix)
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
      d <- if (!is.null(of)) dirname(normalizePath(of)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Create t-SNE plots with marker expression overlays (flowFrame input)
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             downsample_flowframes() (preprocessing/downsampling.R),
#'             homogenize_flowframes() (harmonization/panel_homogenization.R),
#'             create_flowframes_from_matrices() (loading/create_flowframes.R)
#' Output to:  standalone visualization output
#'
#' @param flow_frames Named list of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param markers_of_interest Character vector of markers to plot (NULL = select by mode)
#' @param output_dir Output directory for plots
#' @param perplexity t-SNE perplexity parameter
#' @param max_iter Maximum number of t-SNE iterations
#' @param marker_selection_mode "all", "popup", or "terminal"
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @param seed Random seed
#' @return List with tsne_coords, sample_labels, plots
plot_tsne <- function(flow_frames, sample_names,
                      markers_of_interest = NULL,
                      output_dir,
                      perplexity = 30,
                      max_iter = 1000,
                      marker_selection_mode = "all",
                      output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                      plot_width = 10, plot_height = 8,
                      seed = 42) {
  require(ggplot2)
  require(Rtsne)
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

  # Compute t-SNE
  cat("Computing t-SNE on", nrow(expr_subset), "cells x",
      ncol(expr_subset), "markers (perplexity =", perplexity,
      ", max_iter =", max_iter, ")...\n")
  set.seed(seed)
  tsne_result <- Rtsne::Rtsne(expr_subset, dims = 2, perplexity = perplexity,
                                max_iter = max_iter, check_duplicates = FALSE)
  tsne_coords <- tsne_result$Y

  tsne_df <- data.frame(tSNE1 = tsne_coords[, 1], tSNE2 = tsne_coords[, 2],
                         sample = sample_labels)

  # Plot by sample
  p_sample <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = sample)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_minimal() + labs(title = "t-SNE by Sample")
  save_plot_formats(p_sample, file.path(output_dir, "tsne_sample"),
                    width = plot_width, height = plot_height,
                    formats = output_formats)

  # Plot per marker
  plots <- list(sample = p_sample)
  for (marker in markers_of_interest) {
    tsne_df$expr <- expr_subset[, marker]
    p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = expr)) +
      geom_point(size = 0.3) +
      viridis::scale_color_viridis(name = marker) +
      theme_minimal() + labs(title = paste("t-SNE -", marker))
    save_plot_formats(p, file.path(output_dir, paste0("tsne_", create_safe_filename(marker))),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
    plots[[marker]] <- p
  }

  cat("t-SNE plots saved to:", output_dir, "\n")

  list(tsne_coords = tsne_coords, sample_labels = sample_labels, plots = plots)
}

#' Create t-SNE plots from a pre-built expression matrix
#'
#' Matrix-based alternative to plot_tsne() that accepts a numeric matrix
#' directly, without requiring flowFrame objects. This enables CSV expression
#' data and any matrix-producing script to generate t-SNE visualizations.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             clusters_to_gating_csv() (clustering/leiden.R) [read back],
#'             run_leiden_with_umap() (clustering/leiden.R) [for cluster overlay]
#' Output to:  standalone visualization output
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param sample_labels Character vector of sample labels (one per row)
#' @param cluster_labels Optional integer/character vector of cluster assignments
#' @param output_dir Output directory for plots
#' @param perplexity t-SNE perplexity parameter
#' @param max_iter Maximum number of t-SNE iterations
#' @param tsne_coords Optional pre-computed t-SNE coordinates (2-column matrix).
#'                     If provided, skips t-SNE computation.
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @param seed Random seed
#' @return List with tsne_coords, sample_labels, plots
plot_tsne_from_matrix <- function(expr_matrix, sample_labels = NULL,
                                   cluster_labels = NULL,
                                   output_dir,
                                   perplexity = 30,
                                   max_iter = 1000,
                                   tsne_coords = NULL,
                                   output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                                   plot_width = 10, plot_height = 8,
                                   seed = 42) {
  require(ggplot2)
  require(Rtsne)
  make_outdir(output_dir)

  # Filter markers
  exclude_patterns <- c("^$", "Time", "FSC", "SSC", "index")
  all_markers <- colnames(expr_matrix)
  filtered_idx <- which(!grepl(paste(exclude_patterns, collapse = "|"),
                               all_markers, ignore.case = TRUE))
  markers_of_interest <- all_markers[filtered_idx]

  expr_subset <- expr_matrix[, markers_of_interest, drop = FALSE]

  # Compute or reuse t-SNE coordinates
  if (is.null(tsne_coords)) {
    cat("Computing t-SNE on", nrow(expr_subset), "cells x",
        ncol(expr_subset), "markers (perplexity =", perplexity,
        ", max_iter =", max_iter, ")...\n")
    set.seed(seed)
    tsne_result <- Rtsne::Rtsne(expr_subset, dims = 2, perplexity = perplexity,
                                  max_iter = max_iter, check_duplicates = FALSE)
    tsne_coords <- tsne_result$Y
  } else {
    cat("Using pre-computed t-SNE coordinates for", nrow(tsne_coords), "cells\n")
  }

  if (is.null(sample_labels)) {
    sample_labels <- rep("all", nrow(expr_matrix))
  }

  tsne_df <- data.frame(tSNE1 = tsne_coords[, 1], tSNE2 = tsne_coords[, 2],
                         sample = sample_labels)

  # Plot by sample
  p_sample <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = sample)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_minimal() + labs(title = "t-SNE by Sample")
  save_plot_formats(p_sample, file.path(output_dir, "tsne_sample"),
                    width = plot_width, height = plot_height,
                    formats = output_formats)

  plots <- list(sample = p_sample)

  # Plot by cluster if provided
  if (!is.null(cluster_labels)) {
    tsne_df$cluster <- as.factor(cluster_labels)
    p_cluster <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = cluster)) +
      geom_point(size = 0.3, alpha = 0.5) +
      theme_minimal() + labs(title = "t-SNE by Cluster")
    save_plot_formats(p_cluster, file.path(output_dir, "tsne_cluster"),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
    plots[["cluster"]] <- p_cluster
  }

  # Plot per marker
  for (marker in markers_of_interest) {
    tsne_df$expr <- expr_subset[, marker]
    p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = expr)) +
      geom_point(size = 0.3) +
      viridis::scale_color_viridis(name = marker) +
      theme_minimal() + labs(title = paste("t-SNE -", marker))
    save_plot_formats(p, file.path(output_dir, paste0("tsne_", create_safe_filename(marker))),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
    plots[[marker]] <- p
  }

  cat("t-SNE plots saved to:", output_dir, "\n")
  list(tsne_coords = tsne_coords, sample_labels = sample_labels, plots = plots)
}
