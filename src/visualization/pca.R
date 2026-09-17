#!/usr/bin/env Rscript
# =============================================================================
# visualization/pca.R — PCA Visualization
# =============================================================================
# Refactored from: R_scripts/06_visualization/PCA.R
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
      # Rscript src/foo/bar.R carries no source frame; recover from --file=.
      if (is.null(of)) {
        a <- commandArgs(FALSE); m <- grep("^--file=", a)
        if (length(m)) of <- sub("^--file=", "", a[m[[1]]])
      }
      d <- if (!is.null(of)) dirname(normalizePath(of, mustWork = FALSE)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

# Guarded, NOT unconditional: nextflow/bin/*_cli.R and scripts/*.R set shared_dir
# themselves (the --src_dir convention), and clobbering it here would send every
# source() below to whichever tree this file happens to live in rather than the
# one the caller chose.
if (!exists("shared_dir")) shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Create PCA plots from expression data
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             load_csv_expression() (loading/load_csv_expression.R)
#' Output to:  standalone visualization output
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
                    width = 10, height = 8)

  # Loadings plot
  loadings_df <- as.data.frame(pca$loadings)
  loadings_df$marker <- rownames(loadings_df)
  names(loadings_df)[1:2] <- c("PC1", "PC2")

  p_load <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = marker)) +
    geom_segment(aes(xend = 0, yend = 0), arrow = grid::arrow(length = grid::unit(0.2, "cm")),
                 color = "steelblue") +
    theme_minimal() +
    labs(title = paste(title, "- Loadings"))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_load <- p_load + ggrepel::geom_text_repel(size = 3)
  } else {
    p_load <- p_load + ggplot2::geom_text(size = 3)
  }

  save_plot_formats(p_load, file.path(output_dir, "pca_loadings"),
                    width = 10, height = 10)

  cat("PCA plots saved to:", output_dir, "\n")
  list(pca_result = pca, scatter_plot = p, loadings_plot = p_load)
}
