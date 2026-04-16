#!/usr/bin/env Rscript
# =============================================================================
# visualization/density_plots.R — Density & Ridge Plots
# =============================================================================
# Merged from: R_scripts/06_visualization/Density_plots.R + Plot_raw.R
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Create density plots per marker across samples (flowFrame input)
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             downsample_flowframes() (preprocessing/downsampling.R),
#'             homogenize_flowframes() (harmonization/panel_homogenization.R),
#'             create_flowframes_from_matrices() (loading/create_flowframes.R)
#' Output to:  standalone visualization output
#'
#' @param flow_frames Named list of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param markers Character vector of marker names
#' @param output_dir Output directory
#' @param plot_type "density", "ridge", or "both"
#' @param group_by "sample", "batch", or "condition"
#' @param group_labels Named vector mapping samples to groups
#' @return Invisible NULL
plot_density <- function(flow_frames, sample_names, markers,
                          output_dir, plot_type = "density",
                          group_by = "sample", group_labels = NULL) {
  plot_type <- match.arg(plot_type, c("density", "ridge", "both"))
  require(ggplot2)
  make_outdir(output_dir)

  # Combine data
  expr_list <- list()
  for (i in seq_along(flow_frames)) {
    ff <- flow_frames[[i]]
    data <- explore_data(ff)
    expr <- as.data.frame(flowCore::exprs(ff))
    expr$sample <- sample_names[i]
    if (!is.null(group_labels)) {
      expr$group <- group_labels[sample_names[i]]
    }
    expr_list[[i]] <- expr
  }
  combined <- do.call(rbind, expr_list)

  color_col <- if (!is.null(group_labels)) "group" else "sample"

  for (marker in markers) {
    if (!marker %in% names(combined)) next
    cat("  Plotting:", marker, "\n")

    if (plot_type %in% c("density", "both")) {
      p <- ggplot(combined, aes(x = .data[[marker]],
                                color = .data[[color_col]])) +
        geom_density(alpha = 0.3) +
        theme_minimal() +
        labs(title = paste("Density:", marker), x = marker, y = "Density")
      save_plot_formats(p, file.path(output_dir, paste0("density_", marker)),
                        width = 10, height = 8)
    }

    if (plot_type %in% c("ridge", "both")) {
      if (requireNamespace("ggridges", quietly = TRUE)) {
        p2 <- ggplot(combined, aes(x = .data[[marker]],
                                    y = .data[[color_col]],
                                    fill = .data[[color_col]])) +
          ggridges::geom_density_ridges(alpha = 0.7) +
          theme_minimal() +
          labs(title = paste("Ridge:", marker), x = marker)
        save_plot_formats(p2, file.path(output_dir, paste0("ridge_", marker)),
                          width = 10, height = 8)
      }
    }
  }
  cat("Density plots saved to:", output_dir, "\n")
  invisible(NULL)
}

#' Create density plots from a pre-built expression matrix
#'
#' Matrix-based alternative to plot_density() that accepts a numeric matrix
#' directly, without requiring flowFrame objects. Enables CSV expression data
#' and any matrix-producing script to generate density visualizations.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             load_csv_expression() (loading/load_csv_expression.R)
#' Output to:  standalone visualization output
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param sample_labels Character vector of sample labels (one per row)
#' @param markers Character vector of marker names to plot (NULL = all columns)
#' @param output_dir Output directory
#' @param plot_type "density", "ridge", or "both"
#' @param group_labels Optional named vector mapping samples to groups
#' @return Invisible NULL
plot_density_from_matrix <- function(expr_matrix, sample_labels = NULL,
                                      markers = NULL, output_dir,
                                      plot_type = "density",
                                      group_labels = NULL) {
  plot_type <- match.arg(plot_type, c("density", "ridge", "both"))
  require(ggplot2)
  make_outdir(output_dir)

  if (is.null(sample_labels)) {
    sample_labels <- rep("all", nrow(expr_matrix))
  }
  if (is.null(markers)) {
    markers <- colnames(expr_matrix)
  }

  combined <- as.data.frame(expr_matrix)
  combined$sample <- sample_labels
  if (!is.null(group_labels)) {
    combined$group <- group_labels[sample_labels]
  }

  color_col <- if (!is.null(group_labels)) "group" else "sample"

  for (marker in markers) {
    if (!marker %in% names(combined)) next
    cat("  Plotting:", marker, "\n")

    if (plot_type %in% c("density", "both")) {
      p <- ggplot(combined, aes(x = .data[[marker]],
                                color = .data[[color_col]])) +
        geom_density(alpha = 0.3) +
        theme_minimal() +
        labs(title = paste("Density:", marker), x = marker, y = "Density")
      save_plot_formats(p, file.path(output_dir, paste0("density_", marker)),
                        width = 10, height = 8)
    }

    if (plot_type %in% c("ridge", "both")) {
      if (requireNamespace("ggridges", quietly = TRUE)) {
        p2 <- ggplot(combined, aes(x = .data[[marker]],
                                    y = .data[[color_col]],
                                    fill = .data[[color_col]])) +
          ggridges::geom_density_ridges(alpha = 0.7) +
          theme_minimal() +
          labs(title = paste("Ridge:", marker), x = marker)
        save_plot_formats(p2, file.path(output_dir, paste0("ridge_", marker)),
                          width = 10, height = 8)
      }
    }
  }
  cat("Density plots saved to:", output_dir, "\n")
  invisible(NULL)
}
