#!/usr/bin/env Rscript
# =============================================================================
# visualization/biaxial_plots.R — Biaxial Scatter / Pseudocolor Plots
# =============================================================================
# Classic flow cytometry biaxial plots with density, marker, or cluster coloring
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Create biaxial scatter/pseudocolor plots
#'
#' Classic flow cytometry biaxial plot showing two markers on x and y axes
#' with configurable color mapping. Supports density pseudocolor, marker
#' expression overlay, or cluster assignment coloring. Optionally draws
#' polygon gate boundaries.
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R) [via flowframe],
#'             build_combined_expr_matrix() (R_shared/flowframe_utils.R) [via expr_matrix],
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             run_leiden_with_umap() (clustering/leiden.R) [for cluster_labels]
#' Output to:  standalone visualization output
#'
#' @param flowframe Optional flowFrame object (extracts expression data)
#' @param expr_matrix Optional numeric matrix (rows = cells, cols = markers).
#'                     One of flowframe or expr_matrix must be provided.
#' @param x_marker Name of the marker for the x axis
#' @param y_marker Name of the marker for the y axis
#' @param color_by Coloring mode: "density" (pseudocolor via kde2d),
#'                 a marker name (viridis scale), or "cluster" (categorical)
#' @param sample_label Optional sample name for plot title
#' @param cluster_labels Optional integer/character vector of cluster assignments
#'                        (required when color_by = "cluster")
#' @param gate_coords Optional list with x and y numeric vectors defining a
#'                     polygon gate boundary to overlay
#' @param output_dir Output directory for plots
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @return Invisible ggplot object
plot_biaxial <- function(flowframe = NULL, expr_matrix = NULL,
                          x_marker, y_marker,
                          color_by = "density",
                          sample_label = NULL,
                          cluster_labels = NULL,
                          gate_coords = NULL,
                          output_dir,
                          output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                          plot_width = 10, plot_height = 8) {
  require(ggplot2)
  make_outdir(output_dir)

  # Extract expression data from flowframe or use matrix
  if (!is.null(flowframe)) {
    expr_data <- as.data.frame(flowCore::exprs(flowframe))
    # Try to map column names using marker names
    if (requireNamespace("flowCore", quietly = TRUE)) {
      marker_names <- flowCore::markernames(flowframe)
      named_idx <- which(marker_names != "" & !is.na(marker_names))
      if (length(named_idx) > 0) {
        colnames(expr_data)[named_idx] <- marker_names[named_idx]
      }
    }
  } else if (!is.null(expr_matrix)) {
    expr_data <- as.data.frame(expr_matrix)
  } else {
    stop("One of 'flowframe' or 'expr_matrix' must be provided")
  }

  # Validate markers exist
  if (!x_marker %in% colnames(expr_data)) {
    stop("x_marker '", x_marker, "' not found in data columns")
  }
  if (!y_marker %in% colnames(expr_data)) {
    stop("y_marker '", y_marker, "' not found in data columns")
  }

  plot_df <- data.frame(
    x = expr_data[[x_marker]],
    y = expr_data[[y_marker]]
  )

  # Build title
  title_parts <- paste(x_marker, "vs", y_marker)
  if (!is.null(sample_label)) {
    title_parts <- paste(title_parts, "-", sample_label)
  }

  # Build plot based on color_by mode
  if (color_by == "density") {
    require(MASS)
    cat("Computing 2D kernel density for pseudocolor...\n")
    kde <- MASS::kde2d(plot_df$x, plot_df$y, n = 100)
    # Interpolate density values at each data point
    ix <- findInterval(plot_df$x, kde$x)
    iy <- findInterval(plot_df$y, kde$y)
    ix <- pmax(1, pmin(ix, nrow(kde$z)))
    iy <- pmax(1, pmin(iy, ncol(kde$z)))
    plot_df$density <- mapply(function(i, j) kde$z[i, j], ix, iy)

    p <- ggplot(plot_df, aes(x = x, y = y, color = density)) +
      geom_point(size = 0.3) +
      viridis::scale_color_viridis(name = "Density") +
      theme_minimal() +
      labs(title = paste("Biaxial (density) -", title_parts),
           x = x_marker, y = y_marker)

  } else if (color_by == "cluster") {
    if (is.null(cluster_labels)) {
      stop("cluster_labels must be provided when color_by = 'cluster'")
    }
    plot_df$cluster <- as.factor(cluster_labels)
    p <- ggplot(plot_df, aes(x = x, y = y, color = cluster)) +
      geom_point(size = 0.3, alpha = 0.5) +
      theme_minimal() +
      labs(title = paste("Biaxial (cluster) -", title_parts),
           x = x_marker, y = y_marker)

  } else {
    # color_by is a marker name
    if (!color_by %in% colnames(expr_data)) {
      stop("color_by marker '", color_by, "' not found in data columns")
    }
    plot_df$marker_expr <- expr_data[[color_by]]
    p <- ggplot(plot_df, aes(x = x, y = y, color = marker_expr)) +
      geom_point(size = 0.3) +
      viridis::scale_color_viridis(name = color_by) +
      theme_minimal() +
      labs(title = paste("Biaxial (", color_by, ") -", title_parts),
           x = x_marker, y = y_marker)
  }

  # Overlay gate boundary if provided
  if (!is.null(gate_coords)) {
    gate_df <- data.frame(
      x = gate_coords$x,
      y = gate_coords$y
    )
    # Close the polygon
    gate_df <- rbind(gate_df, gate_df[1, ])
    p <- p + geom_path(data = gate_df, aes(x = x, y = y),
                        color = "red", linewidth = 0.8,
                        inherit.aes = FALSE)
    cat("Gate boundary overlaid on plot\n")
  }

  # Save plot
  safe_x <- create_safe_filename(x_marker)
  safe_y <- create_safe_filename(y_marker)
  safe_color <- create_safe_filename(color_by)
  fname <- paste0("biaxial_", safe_x, "_vs_", safe_y, "_", safe_color)
  if (!is.null(sample_label)) {
    fname <- paste0(fname, "_", create_safe_filename(sample_label))
  }
  save_plot_formats(p, file.path(output_dir, fname),
                    width = plot_width, height = plot_height,
                    formats = output_formats)

  cat("Biaxial plot saved to:", output_dir, "\n")
  invisible(p)
}
