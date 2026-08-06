#!/usr/bin/env Rscript
# =============================================================================
# visualization/violin_plots.R — Violin, Box, and Sina Plots
# =============================================================================
# Expression distribution plots grouped by cluster or condition
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

#' Create violin, box, sina, or combined plots per marker
#'
#' Generates distribution plots of marker expression grouped by a categorical
#' variable (e.g., cluster, condition, sample). Supports violin plots, box
#' plots, sina strip plots (via ggforce), and combined violin + box overlays.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R) [as data.frame],
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             load_csv_expression() (loading/load_csv_expression.R),
#'             calculate_statistics() (R_shared/stats_utils.R)
#' Output to:  standalone visualization output
#'
#' @param expr_data data.table or data.frame with marker columns and a grouping column
#' @param markers Character vector of marker column names to plot
#' @param group_col Name of the grouping column (default "cluster")
#' @param plot_type One of "violin", "box", "sina", "violin_box" (combined)
#' @param output_dir Output directory for plots
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @return Invisible NULL
plot_violin <- function(expr_data, markers,
                         group_col = "cluster",
                         plot_type = "violin_box",
                         output_dir,
                         output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                         plot_width = 12, plot_height = 8) {
  require(ggplot2)
  make_outdir(output_dir)

  # Validate inputs
  if (!group_col %in% colnames(expr_data)) {
    stop("Grouping column '", group_col, "' not found in expr_data")
  }
  available_markers <- intersect(markers, colnames(expr_data))
  if (length(available_markers) == 0) {
    stop("None of the specified markers found in expr_data columns")
  }
  if (length(available_markers) < length(markers)) {
    missing <- setdiff(markers, colnames(expr_data))
    cat("Warning: markers not found in data:", paste(missing, collapse = ", "), "\n")
  }

  group_sizes <- table(expr_data[[group_col]])
  if (any(group_sizes < 2)) {
    small_groups <- names(group_sizes[group_sizes < 2])
    warning("Groups with fewer than 2 observations will be excluded from violin plots: ",
            paste(small_groups, collapse = ", "))
    expr_data <- expr_data[!expr_data[[group_col]] %in% small_groups, ]
  }

  expr_data[[group_col]] <- as.factor(expr_data[[group_col]])

  cat("Generating", plot_type, "plots for", length(available_markers),
      "markers grouped by", group_col, "...\n")

  for (marker in available_markers) {
    cat("  Plotting:", marker, "\n")

    p <- ggplot(expr_data, aes(x = .data[[group_col]], y = .data[[marker]],
                                fill = .data[[group_col]])) +
      theme_minimal() +
      labs(title = paste(marker, "by", group_col),
           x = group_col, y = marker) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    if (plot_type == "violin") {
      p <- p + geom_violin(alpha = 0.7, scale = "width")

    } else if (plot_type == "box") {
      p <- p + geom_boxplot(alpha = 0.7, outlier.size = 0.5)

    } else if (plot_type == "sina") {
      if (!requireNamespace("ggforce", quietly = TRUE)) {
        cat("    ggforce not available, falling back to jitter\n")
        p <- p + geom_jitter(width = 0.2, size = 0.3, alpha = 0.3)
      } else {
        p <- p + ggforce::geom_sina(size = 0.3, alpha = 0.3)
      }

    } else if (plot_type == "violin_box") {
      p <- p +
        geom_violin(alpha = 0.5, scale = "width") +
        geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.5)

    } else {
      stop("Unknown plot_type: '", plot_type, "'. ",
           "Must be one of: violin, box, sina, violin_box")
    }

    safe_marker <- create_safe_filename(marker)
    save_plot_formats(p,
                      file.path(output_dir, paste0(plot_type, "_", safe_marker)),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
  }

  cat("Violin/box plots saved to:", output_dir, "\n")
  invisible(NULL)
}
