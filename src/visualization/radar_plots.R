#!/usr/bin/env Rscript
# =============================================================================
# visualization/radar_plots.R — Radar / Spider Plots
# =============================================================================
# Cluster phenotype profiles visualized as radar charts using coord_polar
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

#' Create radar/spider plots of cluster marker profiles
#'
#' Generates radar (spider) charts showing the marker expression profile
#' of each cluster. When normalize = TRUE, each marker is scaled to 0-1
#' across all clusters. Produces one plot per cluster and one combined
#' overlay plot.
#'
#' Input from: calculate_statistics() (R_shared/stats_utils.R) [cluster medians],
#'             combine_expression_matrices() (analysis/expression_matrices.R),
#'             run_leiden_with_umap() (clustering/leiden.R) [aggregate per cluster]
#' Output to:  standalone visualization output
#'
#' @param cluster_medians data.frame with cluster rows and marker columns.
#'                         Row names (or a "cluster" column) identify clusters.
#' @param markers Character vector of marker column names to include
#' @param output_dir Output directory for plots
#' @param normalize Logical; if TRUE, scale each marker to 0-1 range
#' @param output_formats Named logical list for formats
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @return Invisible NULL
plot_radar <- function(cluster_medians, markers,
                        output_dir,
                        normalize = TRUE,
                        output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                        plot_width = 10, plot_height = 10) {
  require(ggplot2)
  make_outdir(output_dir)

  # Determine cluster identifiers
  if ("cluster" %in% colnames(cluster_medians)) {
    cluster_ids <- as.character(cluster_medians$cluster)
    data_cols <- cluster_medians[, markers, drop = FALSE]
  } else {
    cluster_ids <- rownames(cluster_medians)
    if (is.null(cluster_ids)) {
      cluster_ids <- paste0("Cluster_", seq_len(nrow(cluster_medians)))
    }
    data_cols <- cluster_medians[, markers, drop = FALSE]
  }

  # Validate markers
  available_markers <- intersect(markers, colnames(data_cols))
  if (length(available_markers) == 0) {
    stop("None of the specified markers found in cluster_medians columns")
  }
  data_cols <- data_cols[, available_markers, drop = FALSE]

  # Normalize each marker to 0-1 if requested
  if (normalize) {
    cat("Normalizing marker values to 0-1 range\n")
    for (col in available_markers) {
      vals <- data_cols[[col]]
      rng <- range(vals, na.rm = TRUE)
      if (rng[2] - rng[1] > 0) {
        data_cols[[col]] <- (vals - rng[1]) / (rng[2] - rng[1])
      } else {
        data_cols[[col]] <- 0
      }
    }
  }

  cat("Generating radar plots for", length(cluster_ids), "clusters x",
      length(available_markers), "markers...\n")

  # Helper: build radar plot for one or more clusters
  build_radar <- function(long_df, title_text) {
    p <- ggplot(long_df, aes(x = marker, y = value, group = cluster,
                              color = cluster, fill = cluster)) +
      geom_polygon(alpha = 0.1, linewidth = 0.8) +
      geom_point(size = 2) +
      coord_polar() +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 8),
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey85"),
        legend.position = "bottom"
      ) +
      labs(title = title_text)
    p
  }

  # Individual cluster plots
  for (i in seq_along(cluster_ids)) {
    cid <- cluster_ids[i]
    cat("  Plotting cluster:", cid, "\n")

    long_df <- data.frame(
      marker = factor(available_markers, levels = available_markers),
      value = as.numeric(data_cols[i, ]),
      cluster = cid
    )

    p <- build_radar(long_df, paste("Radar -", cid))
    safe_cid <- create_safe_filename(cid)
    save_plot_formats(p,
                      file.path(output_dir, paste0("radar_", safe_cid)),
                      width = plot_width, height = plot_height,
                      formats = output_formats)
  }

  # Combined overlay plot
  cat("  Plotting combined radar overlay\n")
  combined_long <- do.call(rbind, lapply(seq_along(cluster_ids), function(i) {
    data.frame(
      marker = factor(available_markers, levels = available_markers),
      value = as.numeric(data_cols[i, ]),
      cluster = cluster_ids[i]
    )
  }))

  p_combined <- build_radar(combined_long, "Radar - All Clusters")
  save_plot_formats(p_combined,
                    file.path(output_dir, "radar_combined"),
                    width = plot_width, height = plot_height,
                    formats = output_formats)

  cat("Radar plots saved to:", output_dir, "\n")
  invisible(NULL)
}
