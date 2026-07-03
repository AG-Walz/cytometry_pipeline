#!/usr/bin/env Rscript
# =============================================================================
# visualization/contour_plots.R — 2D Contour / Density Overlay Plots
# =============================================================================
# Contour density overlays for dimensionality reduction coordinates
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
source(file.path(shared_dir, "plot_utils.R"))

#' Create 2D contour density overlay plots
#'
#' Generates contour plots from 2D coordinates (e.g., UMAP, t-SNE, diffusion
#' map output). Supports coloring contours by group labels and two density
#' estimation methods: MASS::kde2d for filled contours, or ggplot2's built-in
#' geom_density_2d for line contours.
#'
#' Input from: plot_umap() / plot_umap_from_matrix() (visualization/umap.R) [coords],
#'             plot_tsne() / plot_tsne_from_matrix() (visualization/tsne.R) [coords],
#'             plot_diffusion_map() (visualization/diffusion_map.R) [coords],
#'             run_umap() (R_shared/stats_utils.R) [raw coords]
#' Output to:  standalone visualization output
#'
#' @param coords data.frame with x and y columns (e.g., UMAP or t-SNE coordinates)
#' @param group_labels Optional character/factor vector of group labels for colored
#'                      contours per group. Must have same length as nrow(coords).
#' @param x_col Name of the x coordinate column (default "x")
#' @param y_col Name of the y coordinate column (default "y")
#' @param method Density estimation method: "kde2d" for MASS::kde2d filled
#'               contours, "geom_density2d" for ggplot2 line contours
#' @param output_dir Output directory for plots
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @return Invisible ggplot object
plot_contour <- function(coords, group_labels = NULL,
                          x_col = "x", y_col = "y",
                          method = "kde2d",
                          output_dir,
                          output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                          plot_width = 10, plot_height = 8) {
  if (!all(c("x", "y") %in% colnames(coords))) {
    stop("coords must have 'x' and 'y' columns")
  }

  require(ggplot2)
  make_outdir(output_dir)

  # Validate coordinate columns
  if (!x_col %in% colnames(coords)) {
    stop("x_col '", x_col, "' not found in coords columns")
  }
  if (!y_col %in% colnames(coords)) {
    stop("y_col '", y_col, "' not found in coords columns")
  }

  plot_df <- data.frame(
    x = coords[[x_col]],
    y = coords[[y_col]]
  )

  has_groups <- !is.null(group_labels)
  if (has_groups) {
    plot_df$group <- as.factor(group_labels)
  }

  cat("Generating contour plot (method =", method, ") for",
      nrow(plot_df), "points...\n")

  if (method == "kde2d") {
    require(MASS)

    if (has_groups) {
      # Filled contour per group using stat_density_2d with filled contours
      p <- ggplot(plot_df, aes(x = x, y = y, color = group)) +
        stat_density_2d(aes(fill = group), geom = "polygon",
                        alpha = 0.15, linewidth = 0.5) +
        geom_point(size = 0.1, alpha = 0.2) +
        theme_minimal() +
        labs(title = "Contour Density by Group",
             x = x_col, y = y_col)
    } else {
      # Single group filled contour with viridis
      p <- ggplot(plot_df, aes(x = x, y = y)) +
        stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                        alpha = 0.7) +
        viridis::scale_fill_viridis(name = "Density") +
        geom_point(size = 0.1, alpha = 0.1) +
        theme_minimal() +
        labs(title = "Contour Density",
             x = x_col, y = y_col)
    }

  } else if (method == "geom_density2d") {
    if (has_groups) {
      p <- ggplot(plot_df, aes(x = x, y = y, color = group)) +
        geom_point(size = 0.1, alpha = 0.2) +
        geom_density_2d(linewidth = 0.5) +
        theme_minimal() +
        labs(title = "Contour Lines by Group",
             x = x_col, y = y_col)
    } else {
      p <- ggplot(plot_df, aes(x = x, y = y)) +
        geom_point(size = 0.1, alpha = 0.2, color = "grey40") +
        geom_density_2d(color = "steelblue", linewidth = 0.5) +
        theme_minimal() +
        labs(title = "Contour Lines",
             x = x_col, y = y_col)
    }

  } else {
    stop("Unknown method: '", method, "'. Must be 'kde2d' or 'geom_density2d'")
  }

  # Save plot
  fname <- paste0("contour_", method)
  if (has_groups) fname <- paste0(fname, "_grouped")
  save_plot_formats(p, file.path(output_dir, fname),
                    width = plot_width, height = plot_height,
                    formats = output_formats)

  cat("Contour plot saved to:", output_dir, "\n")
  invisible(p)
}
