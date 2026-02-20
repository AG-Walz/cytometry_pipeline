#!/usr/bin/env Rscript
# =============================================================================
# R_shared/dimred_dispatch.R — Dimensionality Reduction Dispatcher
# =============================================================================
# Routes to UMAP, t-SNE, or diffusion map via a single config parameter.
# All methods return list(coords, method, plots, extra).
#
# Usage:
#   source(file.path(shared_dir, "dimred_dispatch.R"))
#   result <- run_dimred(expr_matrix, method = "tsne", output_dir = "out/",
#                        sample_labels = labels, params = list(perplexity = 30))
# =============================================================================

# %||% polyfill for R < 4.0
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Dispatch dimensionality reduction to the specified method
#'
#' @param expr_matrix Numeric matrix (cells x markers)
#' @param method One of "umap", "tsne", "diffusion_map"
#' @param sample_labels Character vector (length == nrow), or NULL
#' @param cluster_labels Integer/character vector (length == nrow), or NULL
#' @param condition_labels Character vector (length == nrow), or NULL
#' @param output_dir Directory for plots. If NULL, compute coords only (no plotting).
#' @param params Named list of method-specific parameters
#' @return List with coords (numeric matrix), method (string), plots (list or NULL),
#'   extra (list — eigenvalues for diffusion_map, NULL otherwise)
run_dimred <- function(expr_matrix, method = "umap",
                       sample_labels = NULL, cluster_labels = NULL,
                       condition_labels = NULL, output_dir = NULL,
                       params = list()) {
  method <- tolower(method)
  valid_methods <- c("umap", "tsne", "diffusion_map")
  if (!method %in% valid_methods) {
    stop("Unknown dimred method '", method, "'. Valid: ",
         paste(valid_methods, collapse = ", "))
  }

  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  cat("Dispatching dimensionality reduction: method =", method, "\n")

  # Resolve visualization directory
  if (exists("shared_dir")) {
    viz_dir <- file.path(dirname(shared_dir), "visualization")
  } else {
    viz_dir <- file.path(dirname(getwd()), "src", "visualization")
  }

  # Common output format params
  output_formats <- params$output_formats %||% list(pdf = FALSE, png = TRUE, svg = FALSE)
  plot_width  <- params$plot_width %||% 10
  plot_height <- params$plot_height %||% 8
  seed <- params$seed %||% 42

  result <- switch(method,

    umap = {
      source(file.path(viz_dir, "umap.R"))
      umap_result <- plot_umap_from_matrix(
        expr_matrix     = expr_matrix,
        sample_labels   = sample_labels,
        cluster_labels  = cluster_labels,
        condition_labels = condition_labels,
        output_dir      = output_dir,
        n_neighbors     = params$n_neighbors %||% params$umap_n_neighbors %||% 15,
        min_dist        = params$umap_min_dist %||% params$min_dist %||% 0.1,
        umap_coords     = params$umap_coords %||% NULL,
        output_formats  = output_formats,
        plot_width      = plot_width,
        plot_height     = plot_height,
        seed            = seed
      )
      list(
        coords = umap_result$umap_coords,
        method = "umap",
        plots  = umap_result$plots,
        extra  = NULL
      )
    },

    tsne = {
      source(file.path(viz_dir, "tsne.R"))
      tsne_result <- plot_tsne_from_matrix(
        expr_matrix    = expr_matrix,
        sample_labels  = sample_labels,
        cluster_labels = cluster_labels,
        output_dir     = output_dir,
        perplexity     = params$tsne_perplexity %||% params$perplexity %||% 30,
        max_iter       = params$tsne_max_iter %||% params$max_iter %||% 1000,
        tsne_coords    = params$tsne_coords %||% NULL,
        output_formats = output_formats,
        plot_width     = plot_width,
        plot_height    = plot_height,
        seed           = seed
      )
      list(
        coords = tsne_result$tsne_coords,
        method = "tsne",
        plots  = tsne_result$plots,
        extra  = NULL
      )
    },

    diffusion_map = {
      source(file.path(viz_dir, "diffusion_map.R"))
      dm_result <- plot_diffusion_map(
        expr_matrix    = expr_matrix,
        sample_labels  = sample_labels,
        cluster_labels = cluster_labels,
        output_dir     = output_dir,
        n_components   = params$dm_n_components %||% params$n_components %||% 3,
        sigma          = params$dm_sigma %||% params$sigma %||% "local",
        seed           = seed,
        output_formats = output_formats,
        plot_width     = plot_width,
        plot_height    = plot_height
      )
      list(
        coords = dm_result$dm_coords,
        method = "diffusion_map",
        plots  = dm_result$plots,
        extra  = list(eigenvalues = dm_result$eigenvalues)
      )
    }
  )

  # Validate output
  stopifnot(is.matrix(result$coords) || is.data.frame(result$coords))
  stopifnot(nrow(result$coords) == nrow(expr_matrix))
  cat("Dimensionality reduction complete:", ncol(result$coords), "dimensions\n")

  result
}
