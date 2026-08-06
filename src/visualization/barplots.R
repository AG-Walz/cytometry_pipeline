#!/usr/bin/env Rscript
# =============================================================================
# visualization/barplots.R — Composition Barplots
# =============================================================================
# Refactored from: R_scripts/06_visualization/Barplots.R
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
  if (!x_col %in% colnames(data)) {
    stop("x_col '", x_col, "' not found in data")
  }
  if (!fill_col %in% colnames(data)) {
    stop("fill_col '", fill_col, "' not found in data")
  }
  if (nrow(data) == 0) {
    warning("Empty data frame passed to plot_composition_barplot; returning NULL")
    return(invisible(NULL))
  }

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
