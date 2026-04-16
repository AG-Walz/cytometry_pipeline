#!/usr/bin/env Rscript
# =============================================================================
# R_shared/clustering_dispatch.R — Clustering Method Dispatcher
# =============================================================================
# Routes to any clustering method via a single config parameter.
# All methods return an integer vector of length == nrow(expr_matrix).
# DBSCAN output (list) is normalized: cluster_ids extracted, metadata stored as attr.
#
# Usage:
#   source(file.path(shared_dir, "clustering_dispatch.R"))
#   ids <- run_clustering(expr_matrix, method = "flowsom", params = list(n_clusters = 20))
# =============================================================================

# %||% polyfill for R < 4.0
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Dispatch clustering to the specified method
#'
#' @param expr_matrix Numeric matrix (cells x markers)
#' @param method One of "leiden", "flowsom", "phenograph", "dbscan"
#' @param params Named list of method-specific parameters
#' @return Integer vector of cluster IDs (length == nrow(expr_matrix))
run_clustering <- function(expr_matrix, method = "leiden", params = list()) {
  method <- tolower(method)
  valid_methods <- c("leiden", "flowsom", "phenograph", "dbscan")
  if (!method %in% valid_methods) {
    stop("Unknown clustering method '", method, "'. Valid: ",
         paste(valid_methods, collapse = ", "))
  }

  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  n_cells <- nrow(expr_matrix)
  cat("Dispatching clustering: method =", method, ", cells =", n_cells, "\n")

  # Resolve shared directories
  if (exists("shared_dir")) {
    cluster_dir <- file.path(dirname(shared_dir), "clustering")
  } else {
    cluster_dir <- file.path(dirname(getwd()), "src", "clustering")
  }

  cluster_ids <- switch(method,

    leiden = {
      tryCatch(source(file.path(cluster_dir, "leiden.R")),
               error = function(e) stop("Failed to load Leiden module: ", e$message,
                                        "\nExpected at: ", file.path(cluster_dir, "leiden.R")))
      run_leiden(
        expr_matrix,
        n_neighbors = params$n_neighbors %||% 30,
        resolution  = params$resolution %||% 0.5,
        seed        = params$seed %||% 42
      )
    },

    flowsom = {
      tryCatch(source(file.path(cluster_dir, "flowsom.R")),
               error = function(e) stop("Failed to load FlowSOM module: ", e$message,
                                        "\nExpected at: ", file.path(cluster_dir, "flowsom.R")))
      run_flowsom(
        expr_matrix,
        n_clusters = params$n_clusters %||% 20,
        grid_xdim  = params$grid_xdim %||% 10,
        grid_ydim  = params$grid_ydim %||% 10,
        seed       = params$seed %||% 42
      )
    },

    phenograph = {
      tryCatch(source(file.path(cluster_dir, "phenograph.R")),
               error = function(e) stop("Failed to load PhenoGraph module: ", e$message,
                                        "\nExpected at: ", file.path(cluster_dir, "phenograph.R")))
      run_phenograph(
        expr_matrix,
        k    = params$k %||% params$n_neighbors %||% 30,
        seed = params$seed %||% 42
      )
    },

    dbscan = {
      tryCatch(source(file.path(cluster_dir, "dbscan.R")),
               error = function(e) stop("Failed to load DBSCAN module: ", e$message,
                                        "\nExpected at: ", file.path(cluster_dir, "dbscan.R")))
      result <- run_dbscan(
        expr_matrix,
        eps      = params$dbscan_eps %||% params$eps %||% NULL,
        min_pts  = params$dbscan_min_pts %||% params$min_pts %||% 5,
        auto_eps = is.null(params$dbscan_eps %||% params$eps %||% NULL),
        seed     = params$seed %||% 42
      )
      # Normalize: extract cluster_ids from list, store metadata as attribute
      ids <- result$cluster_ids
      attr(ids, "dbscan_metadata") <- list(
        noise_indices  = result$noise_indices,
        n_clusters     = result$n_clusters,
        eps_used       = result$eps_used,
        noise_fraction = result$noise_fraction
      )
      ids
    }
  )

  # Ensure integer output (some methods like leiden return numeric)
  if (!is.integer(cluster_ids)) {
    meta_backup <- attr(cluster_ids, "dbscan_metadata")
    cluster_ids <- as.integer(cluster_ids)
    if (!is.null(meta_backup)) attr(cluster_ids, "dbscan_metadata") <- meta_backup
  }
  stopifnot(is.integer(cluster_ids))
  stopifnot(length(cluster_ids) == n_cells)
  n_unique <- length(unique(cluster_ids[cluster_ids > 0]))
  cat("Clustering complete:", n_unique, "clusters assigned\n")

  cluster_ids
}
