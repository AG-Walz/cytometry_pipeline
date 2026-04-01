# =============================================================================
# R_shared/pipeline_utils.R — Pipeline boilerplate utilities
# =============================================================================
# Eliminates repeated config loading, module sourcing, and param construction
# from pipeline scripts.
# =============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Initialize a pipeline: parse --config, source shared modules, load defaults
#'
#' @param pipeline_name  Short name for logging (e.g. "Pipeline A")
#' @param default_config Relative path from project root to the default config
#'                       (e.g. "config/pipeline_A_config.R")
#' @return Named list with base_dir, shared_dir, and the pipeline start_time
init_pipeline <- function(pipeline_name, default_config) {
  cat(sprintf("=== %s ===\n", pipeline_name))
  start_time <- Sys.time()

  base_dir   <- dirname(getwd())
  shared_dir <- file.path(base_dir, "src", "R_shared")

  # Parse --config from command line
  args <- commandArgs(trailingOnly = TRUE)
  config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
  if (!is.null(config_file)) {
    cat("Loading config:", config_file, "\n")
    source(config_file)
  } else {
    cat("No --config provided. Using default config.\n")
    source(file.path(base_dir, default_config))
  }

  # Source core shared modules
  source(file.path(shared_dir, "io_utils.R"))
  source(file.path(shared_dir, "memory_utils.R"))

  list(base_dir = base_dir, shared_dir = shared_dir, start_time = start_time)
}

#' Print pipeline completion message with elapsed time
#'
#' @param pipeline_name Short name (e.g. "Pipeline A")
#' @param start_time    The Sys.time() value from init_pipeline()
finish_pipeline <- function(pipeline_name, start_time) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("\n=== %s complete (%.1f minutes) ===\n", pipeline_name, elapsed))
}

#' Build the params list for run_clustering() from a user config
#'
#' Extracts only the params that the clustering dispatcher needs,
#' so pipeline scripts don't have to repeat this construction.
#'
#' @param params The user's params list
#' @return A list suitable for run_clustering(..., params = result)
clustering_params <- function(params) {
  list(
    n_neighbors    = params$n_neighbors    %||% 15,
    resolution     = params$resolution     %||% 0.5,
    n_clusters     = params$n_clusters     %||% 20,
    grid_xdim      = params$grid_xdim      %||% 10,
    grid_ydim      = params$grid_ydim      %||% 10,
    k              = params$k %||% params$n_neighbors %||% 30,
    dbscan_eps     = params$dbscan_eps,
    dbscan_min_pts = params$dbscan_min_pts %||% 5,
    seed           = params$seed           %||% 42
  )
}

#' Build the params list for run_dimred() from a user config
#'
#' @param params The user's params list
#' @return A list suitable for run_dimred(..., params = result)
dimred_params <- function(params) {
  list(
    n_neighbors     = params$n_neighbors     %||% 15,
    umap_min_dist   = params$umap_min_dist   %||% 0.1,
    tsne_perplexity = params$tsne_perplexity %||% 30,
    tsne_max_iter   = params$tsne_max_iter   %||% 1000,
    dm_n_components = params$dm_n_components %||% 3,
    dm_sigma        = params$dm_sigma        %||% "local",
    output_formats  = params$output_formats  %||% list(pdf = FALSE, png = TRUE, svg = FALSE),
    plot_width      = params$plot_width      %||% 10,
    plot_height     = params$plot_height     %||% 8,
    seed            = params$seed            %||% 42
  )
}

#' Build the params list for run_batch_correction_dispatch() from a user config
#'
#' @param params The user's params list
#' @return A list suitable for the batch correction dispatcher
batch_params <- function(params) {
  list(
    batch_do_pca           = params$batch_do_pca           %||% FALSE,
    combat_method          = params$combat_method          %||% "parametric",
    combat_covariates      = params$combat_covariates      %||% NULL,
    cytonorm_train_indices = params$cytonorm_train_indices,
    cytonorm_n_clusters    = params$cytonorm_n_clusters    %||% 10,
    cytonorm_goal          = params$cytonorm_goal          %||% "mean",
    seed                   = params$seed                   %||% 42
  )
}
