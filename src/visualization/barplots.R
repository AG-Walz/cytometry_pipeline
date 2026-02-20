#!/usr/bin/env Rscript
# =============================================================================
# visualization/barplots.R — Composition Barplots
# =============================================================================
# Refactored from: R_scripts/06_visualization/Barplots.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Create composition barplots from cluster data
#'
#' Input from: tracking_matrix_to_long() (analysis/cluster_splitting.R)
#'
#' @param data Data frame with sample, cluster, and optionally condition columns
#' @param output_dir Output directory
#' @param x_col Column for x-axis (default "cluster")
#' @param fill_col Column for fill (default "sample")
#' @param position "stack", "dodge", or "fill"
#' @param title Plot title
plot_composition_barplot <- function(data, output_dir,
                                      x_col = "cluster",
                                      fill_col = "sample",
                                      position = "stack",
                                      title = "Cell Composition") {
  require(ggplot2)
  make_outdir(output_dir)

  p <- ggplot(data, aes(x = .data[[x_col]], fill = .data[[fill_col]])) +
    geom_bar(position = position, color = "black", linewidth = 0.2) +
    theme_bw() +
    labs(title = title, x = x_col, y = "Cell Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot_formats(p, file.path(output_dir, paste0("barplot_", x_col, "_", fill_col)),
                    width = 14, height = 8)
  cat("Barplot saved to:", output_dir, "\n")
  invisible(p)
}
