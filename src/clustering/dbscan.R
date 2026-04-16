#!/usr/bin/env Rscript
# =============================================================================
# clustering/dbscan.R — DBSCAN Density-Based Clustering
# =============================================================================
# DBSCAN clustering for rare population detection in flow cytometry data.
# Includes automatic epsilon estimation from k-distance plots.
# Noise points are assigned cluster ID 0.
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))

#' Run DBSCAN clustering on expression data
#'
#' Performs density-based spatial clustering (DBSCAN) on a numeric expression
#' matrix. If auto_eps is TRUE, the epsilon parameter is estimated from the
#' knee point of the k-distance plot. Noise points receive cluster ID 0.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             normalize_expression() (preprocessing/normalization.R)
#' Output to:  clusters_to_gating_csv() (R_shared/cluster_utils.R),
#'             validate_clusters() (clustering/cluster_validation.R),
#'             annotate_clusters() (clustering/cell_annotation.R)
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param eps Epsilon neighborhood radius. If NULL and auto_eps=TRUE, estimated
#'   automatically from the k-distance plot knee point.
#' @param min_pts Minimum points to form a dense region (default 5)
#' @param auto_eps Logical; if TRUE and eps is NULL, estimate eps from k-distance
#'   knee point using min_pts as k (default TRUE)
#' @param seed Random seed (used for subsampling during eps estimation)
#' @return List with components:
#'   \item{cluster_ids}{Integer vector of cluster assignments (0 = noise)}
#'   \item{noise_indices}{Integer vector of row indices classified as noise}
#'   \item{n_clusters}{Number of clusters found (excluding noise)}
#'   \item{eps_used}{The epsilon value used}
#'   \item{noise_fraction}{Fraction of cells classified as noise}
run_dbscan <- function(expr_matrix, eps = NULL, min_pts = 5,
                       auto_eps = TRUE, seed = 42) {
  require(dbscan)

  set.seed(seed)
  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  stopifnot(min_pts >= 1)

  n_cells <- nrow(expr_matrix)
  cat("DBSCAN clustering on", n_cells, "cells x", ncol(expr_matrix), "markers\n")

  # Estimate epsilon from k-distance plot if not provided
  if (is.null(eps) && auto_eps) {
    cat("Estimating eps from k-distance plot (k =", min_pts, ")...\n")
    eps <- estimate_eps(expr_matrix, k = min_pts, seed = seed)
    if (is.na(eps) || is.nan(eps) || is.infinite(eps)) {
      warning("Automatic eps estimation failed; using default eps = 0.5")
      eps <- 0.5
    }
    cat("Estimated eps =", round(eps, 4), "\n")
  } else if (is.null(eps)) {
    stop("eps must be provided when auto_eps = FALSE")
  }

  cat("Running DBSCAN (eps =", round(eps, 4), ", minPts =", min_pts, ")...\n")

  db_result <- tryCatch({
    dbscan::dbscan(expr_matrix, eps = eps, minPts = min_pts)
  }, error = function(e) {
    cat("DBSCAN failed:", e$message, "\n")
    cat("Retrying with scaled data...\n")
    scaled_matrix <- scale(expr_matrix)
    scaled_eps <- estimate_eps(scaled_matrix, k = min_pts, seed = seed)
    cat("Re-estimated eps on scaled data =", round(scaled_eps, 4), "\n")
    dbscan::dbscan(scaled_matrix, eps = scaled_eps, minPts = min_pts)
  })

  cluster_ids <- as.integer(db_result$cluster)
  noise_indices <- which(cluster_ids == 0)
  n_clusters <- max(cluster_ids)
  noise_fraction <- length(noise_indices) / n_cells

  cat("DBSCAN found", n_clusters, "clusters\n")
  cat("Noise points:", length(noise_indices),
      "(", round(noise_fraction * 100, 1), "%)\n")
  if (n_clusters > 0) {
    non_noise <- cluster_ids[cluster_ids > 0]
    cat("Cluster sizes:", paste(table(non_noise), collapse = ", "), "\n")
  }

  # Warn if too many noise points
  if (noise_fraction > 0.5) {
    cat("Warning: >50% noise — consider decreasing eps or min_pts\n")
  }
  if (n_clusters == 0) {
    cat("Warning: no clusters found — all points classified as noise.",
        " Try increasing eps or decreasing min_pts.\n")
  }

  list(
    cluster_ids = cluster_ids,
    noise_indices = noise_indices,
    n_clusters = n_clusters,
    eps_used = eps,
    noise_fraction = noise_fraction
  )
}

#' Estimate epsilon from k-distance plot knee point
#'
#' Computes the k-nearest-neighbor distances for each point, sorts them,
#' and finds the knee point using the maximum curvature method.
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param k Number of nearest neighbors (typically same as min_pts)
#' @param max_sample Maximum number of cells to sample for estimation
#' @param seed Random seed for subsampling
#' @return Numeric scalar: estimated epsilon value
estimate_eps <- function(expr_matrix, k = 5, max_sample = 10000, seed = 42) {
  require(dbscan)

  set.seed(seed)
  n_cells <- nrow(expr_matrix)

  # Subsample for efficiency if dataset is large
  if (n_cells > max_sample) {
    sample_idx <- sample(n_cells, max_sample)
    sample_matrix <- expr_matrix[sample_idx, , drop = FALSE]
  } else {
    sample_matrix <- expr_matrix
  }

  # Compute k-nearest-neighbor distances
  knn_dists <- tryCatch({
    dbscan::kNNdist(sample_matrix, k = k)
  }, error = function(e) {
    cat("kNNdist failed:", e$message, "\n")
    cat("Using manual distance computation...\n")
    dist_matrix <- as.matrix(dist(sample_matrix))
    apply(dist_matrix, 1, function(row) sort(row)[k + 1])
  })

  # Use the k-th neighbor distance (last column if matrix)
  if (is.matrix(knn_dists)) {
    k_dists <- knn_dists[, ncol(knn_dists)]
  } else {
    k_dists <- knn_dists
  }

  # Sort distances
  sorted_dists <- sort(k_dists)
  n <- length(sorted_dists)

  # Find knee point using maximum curvature
  # Line from first to last point; knee = max perpendicular distance
  x <- seq_len(n)
  y <- sorted_dists

  # Normalize to [0, 1] range for distance calculation
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / max(max(y) - min(y), 1e-10)

  # Line from (0, y_norm[1]) to (1, y_norm[n])
  line_vec <- c(1, y_norm[n] - y_norm[1])
  line_len <- sqrt(sum(line_vec^2))

  # Perpendicular distance of each point to the line
  perp_dist <- abs((x_norm - 0) * (y_norm[n] - y_norm[1]) -
                   (y_norm - y_norm[1]) * 1) / line_len

  knee_idx <- which.max(perp_dist)
  eps_estimate <- sorted_dists[knee_idx]

  eps_estimate
}
