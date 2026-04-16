#!/usr/bin/env Rscript
# =============================================================================
# visualization/diffusion_map.R — Diffusion Map Visualization
# =============================================================================
# Diffusion map dimensionality reduction using the destiny package
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Create diffusion map plots from expression data
#'
#' Uses the destiny package to compute diffusion maps and visualize
#' cells in diffusion component space. Generates plots of DC1 vs DC2
#' and DC1 vs DC3, colored by sample, cluster, and individual markers.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             load_csv_expression() (loading/load_csv_expression.R),
#'             run_leiden_with_umap() (clustering/leiden.R) [for cluster overlay]
#' Output to:  standalone visualization output
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param sample_labels Character vector of sample labels (one per row)
#' @param cluster_labels Optional integer/character vector of cluster assignments
#' @param output_dir Output directory for plots
#' @param n_components Number of diffusion components to compute
#' @param sigma Diffusion kernel width: "local" for adaptive, or numeric value
#' @param seed Random seed
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @return List with dm_coords, eigenvalues, plots
plot_diffusion_map <- function(expr_matrix, sample_labels = NULL,
                                cluster_labels = NULL,
                                output_dir,
                                n_components = 3,
                                sigma = "local",
                                seed = 42,
                                output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                                plot_width = 10, plot_height = 8) {
  require(ggplot2)
  require(destiny)
  make_outdir(output_dir)

  # Filter markers (exclude scatter/time/index)
  exclude_patterns <- c("^$", "Time", "FSC", "SSC", "index")
  all_markers <- colnames(expr_matrix)
  filtered_idx <- which(!grepl(paste(exclude_patterns, collapse = "|"),
                               all_markers, ignore.case = TRUE))
  markers_of_interest <- all_markers[filtered_idx]

  expr_subset <- expr_matrix[, markers_of_interest, drop = FALSE]

  if (is.null(sample_labels)) {
    sample_labels <- rep("all", nrow(expr_matrix))
  }

  # Compute diffusion map
  cat("Computing diffusion map on", nrow(expr_subset), "cells x",
      ncol(expr_subset), "markers (n_components =", n_components, ")...\n")
  set.seed(seed)
  dm <- destiny::DiffusionMap(expr_subset, n_dcs = n_components, sigma = sigma)
  dm_coords <- as.matrix(destiny::eigenvectors(dm)[, seq_len(n_components)])
  colnames(dm_coords) <- paste0("DC", seq_len(n_components))
  eigenvalues <- destiny::eigenvalues(dm)[seq_len(n_components)]

  cat("Diffusion map computed. Top eigenvalues:",
      paste(round(eigenvalues, 4), collapse = ", "), "\n")

  plots <- list()

  # Define component pairs to plot
  component_pairs <- list(c(1, 2))
  if (n_components >= 3) {
    component_pairs <- c(component_pairs, list(c(1, 3)))
  }

  # Plot each component pair colored by sample
  for (pair in component_pairs) {
    dc_x <- paste0("DC", pair[1])
    dc_y <- paste0("DC", pair[2])
    pair_label <- paste0("dc", pair[1], "_vs_dc", pair[2])

    dm_df <- data.frame(
      x = dm_coords[, pair[1]],
      y = dm_coords[, pair[2]],
      sample = sample_labels
    )

    p_sample <- ggplot(dm_df, aes(x = x, y = y, color = sample)) +
      geom_point(size = 0.3, alpha = 0.5) +
      theme_minimal() +
      labs(title = paste("Diffusion Map by Sample -", dc_x, "vs", dc_y),
           x = dc_x, y = dc_y)
    save_plot_formats(p_sample,
                      file.path(output_dir, paste0("diffmap_sample_", pair_label)),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
    plots[[paste0("sample_", pair_label)]] <- p_sample

    # Plot by cluster if provided
    if (!is.null(cluster_labels)) {
      dm_df$cluster <- as.factor(cluster_labels)
      p_cluster <- ggplot(dm_df, aes(x = x, y = y, color = cluster)) +
        geom_point(size = 0.3, alpha = 0.5) +
        theme_minimal() +
        labs(title = paste("Diffusion Map by Cluster -", dc_x, "vs", dc_y),
             x = dc_x, y = dc_y)
      save_plot_formats(p_cluster,
                        file.path(output_dir, paste0("diffmap_cluster_", pair_label)),
                        width = plot_width, height = plot_height,
                        formats = output_formats)
      plots[[paste0("cluster_", pair_label)]] <- p_cluster
    }

    # Plot per marker
    for (marker in markers_of_interest) {
      dm_df$expr <- expr_subset[, marker]
      p_marker <- ggplot(dm_df, aes(x = x, y = y, color = expr)) +
        geom_point(size = 0.3) +
        viridis::scale_color_viridis(name = marker) +
        theme_minimal() +
        labs(title = paste("Diffusion Map -", marker, "-", dc_x, "vs", dc_y),
             x = dc_x, y = dc_y)
      save_plot_formats(p_marker,
                        file.path(output_dir,
                                  paste0("diffmap_", create_safe_filename(marker),
                                         "_", pair_label)),
                        width = plot_width, height = plot_height,
                        formats = output_formats)
      plots[[paste0(marker, "_", pair_label)]] <- p_marker
    }
  }

  cat("Diffusion map plots saved to:", output_dir, "\n")
  list(dm_coords = dm_coords, eigenvalues = eigenvalues, plots = plots)
}
