#!/usr/bin/env Rscript
# =============================================================================
# visualization/heatmaps.R — Expression Heatmaps
# =============================================================================
# Refactored from: R_scripts/06_visualization/Heatmap.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Create expression heatmap from a matrix
#'
#' Input from: combine_expression_matrices() (analysis/expression_matrices.R),
#'             create_expression_matrices() (analysis/expression_matrices.R)
#'
#' @param expr_matrix Numeric matrix (rows = clusters/samples, cols = markers)
#' @param output_path Output file path
#' @param title Heatmap title
#' @param cluster_rows Cluster rows
#' @param cluster_cols Cluster columns
#' @param scale "row", "column", or "none"
#' @param color_palette Color palette function
#' @param width Plot width
#' @param height Plot height
#' @return pheatmap object (invisibly)
plot_heatmap <- function(expr_matrix, output_path,
                          title = "Expression Heatmap",
                          cluster_rows = TRUE, cluster_cols = TRUE,
                          scale = "row",
                          color_palette = NULL,
                          width = 10, height = 8) {
  require(pheatmap)

  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(
      c("navy", "white", "firebrick3")
    )(100)
  }

  # Handle Inf/NaN
  expr_matrix[is.infinite(expr_matrix)] <- NA
  expr_matrix[is.nan(expr_matrix)] <- NA

  # Ensure output directory exists
  make_outdir(dirname(output_path))

  # Strip extension from output_path to get base name
  filename_base <- sub("\\.[^.]+$", "", output_path)

  draw_heatmap <- function() {
    pheatmap::pheatmap(
      expr_matrix,
      main = title,
      scale = scale,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      color = color_palette,
      border_color = NA,
      show_rownames = nrow(expr_matrix) <= 50,
      show_colnames = TRUE,
      fontsize_row = max(6, 12 - nrow(expr_matrix) / 10),
      fontsize_col = 8
    )
  }

  save_base_plot_formats(filename_base, draw_heatmap,
                         width = width, height = height)
  invisible(NULL)
}

#' Create median expression heatmap grouped by cluster
#'
#' @param expr_data Data frame with expression and cluster columns
#' @param cluster_col Name of cluster column
#' @param marker_cols Character vector of marker column names
#' @param output_path Output file path
#' @param title Heatmap title
plot_cluster_heatmap <- function(expr_data, cluster_col, marker_cols,
                                  output_path, title = "Cluster Expression") {
  # Calculate median expression per cluster
  median_expr <- aggregate(
    expr_data[, marker_cols],
    by = list(cluster = expr_data[[cluster_col]]),
    FUN = median, na.rm = TRUE
  )
  rownames(median_expr) <- median_expr$cluster
  median_expr$cluster <- NULL

  plot_heatmap(as.matrix(median_expr), output_path, title = title)
}
