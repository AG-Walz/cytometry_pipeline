#!/usr/bin/env Rscript
# =============================================================================
# visualization/pca.R — PCA Visualization
# =============================================================================
# Refactored from: R_scripts/06_visualization/PCA.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Create PCA plots from expression data
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param sample_labels Character vector of sample labels
#' @param output_dir Output directory
#' @param n_components Number of PCs to compute
#' @param color_by "sample" or a character vector of group labels
#' @param title Plot title prefix
#' @return List with pca_result, plots
plot_pca <- function(expr_matrix, sample_labels = NULL,
                     output_dir, n_components = 2,
                     color_by = "sample", title = "PCA") {
  require(ggplot2)
  make_outdir(output_dir)

  pca <- run_pca(expr_matrix, n_components = n_components)
  pca_df <- as.data.frame(pca$coordinates)
  names(pca_df) <- paste0("PC", seq_len(n_components))

  if (!is.null(sample_labels)) {
    pca_df$sample <- sample_labels
  }
  if (is.character(color_by) && length(color_by) > 1) {
    pca_df$group <- color_by
    fill_col <- "group"
  } else {
    fill_col <- "sample"
  }

  var_expl <- round(pca$variance_explained * 100, 1)

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[fill_col]])) +
    geom_point(size = 1, alpha = 0.6) +
    theme_minimal() +
    labs(title = title,
         x = paste0("PC1 (", var_expl[1], "%)"),
         y = paste0("PC2 (", var_expl[2], "%)"))

  save_plot_formats(p, file.path(output_dir, "pca_plot"),
                    width = 8, height = 6)

  # Loadings plot
  loadings_df <- as.data.frame(pca$loadings)
  loadings_df$marker <- rownames(loadings_df)
  names(loadings_df)[1:2] <- c("PC1", "PC2")

  p_load <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = marker)) +
    geom_segment(aes(xend = 0, yend = 0), arrow = grid::arrow(length = grid::unit(0.2, "cm")),
                 color = "steelblue") +
    ggrepel::geom_text_repel(size = 3) +
    theme_minimal() +
    labs(title = paste(title, "- Loadings"))

  save_plot_formats(p_load, file.path(output_dir, "pca_loadings"),
                    width = 8, height = 8)

  cat("PCA plots saved to:", output_dir, "\n")
  list(pca_result = pca, scatter_plot = p, loadings_plot = p_load)
}
