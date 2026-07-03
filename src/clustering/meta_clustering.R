#!/usr/bin/env Rscript
# =============================================================================
# clustering/meta_clustering.R — Resolution Sweep & Consensus Clustering
# =============================================================================
# Provides resolution parameter optimization via sweep with internal
# validation metrics, plus consensus clustering from multiple independent
# clustering runs.
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

if (!exists("shared_dir")) {
  shared_dir <- cyto_shared_dir()
}
cluster_dir <- file.path(dirname(shared_dir), "clustering")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

# Source Leiden for run_leiden()
source(file.path(cluster_dir, "leiden.R"))

#' Run resolution sweep for clustering parameter optimization
#'
#' Tries multiple resolution values for a given clustering method, evaluates
#' each with an internal validation metric (silhouette or Calinski-Harabasz),
#' and returns the best clustering along with all scores.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             normalize_expression() (preprocessing/normalization.R)
#' Output to:  clusters_to_gating_csv() (clustering/leiden.R),
#'             validate_clusters() (clustering/cluster_validation.R)
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param method Clustering method: "leiden" (default)
#' @param resolutions Numeric vector of resolution values to try
#' @param metric Validation metric: "silhouette" or "calinski_harabasz"
#' @param n_neighbors Number of nearest neighbors for graph-based methods
#' @param seed Random seed
#' @return List with components:
#'   \item{best_clusters}{Integer vector of cluster assignments from best resolution}
#'   \item{best_param}{The resolution value that scored best}
#'   \item{scores_df}{Data frame with resolution, n_clusters, score columns}
#'   \item{all_clusterings}{Named list of cluster vectors (one per resolution)}
run_resolution_sweep <- function(expr_matrix, method = "leiden",
                                 resolutions = seq(0.2, 2.0, by = 0.2),
                                 metric = "silhouette",
                                 n_neighbors = 30, seed = 42) {
  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  stopifnot(method %in% c("leiden"))
  stopifnot(metric %in% c("silhouette", "calinski_harabasz"))
  stopifnot(length(resolutions) > 0)

  set.seed(seed)
  n_cells <- nrow(expr_matrix)

  cat("Resolution sweep:", length(resolutions), "values for", method,
      "(metric:", metric, ")\n")

  # Pre-compute distance matrix for metric evaluation (subsample if large)
  eval_sample_size <- min(5000, n_cells)
  if (n_cells > eval_sample_size) {
    eval_idx <- sample(n_cells, eval_sample_size)
    eval_matrix <- expr_matrix[eval_idx, , drop = FALSE]
    cat("Subsampling", eval_sample_size, "cells for metric evaluation\n")
  } else {
    eval_idx <- seq_len(n_cells)
    eval_matrix <- expr_matrix
  }
  eval_dist <- dist(eval_matrix)

  all_clusterings <- list()
  scores <- numeric(length(resolutions))
  n_clusters_vec <- integer(length(resolutions))

  for (i in seq_along(resolutions)) {
    res <- resolutions[i]
    cat("  Resolution", res, "...")

    clusters <- tryCatch({
      run_leiden(expr_matrix, n_neighbors = n_neighbors,
                 resolution = res, seed = seed)
    }, error = function(e) {
      cat(" FAILED:", e$message, "\n")
      return(NULL)
    })

    if (is.null(clusters)) {
      scores[i] <- NA
      n_clusters_vec[i] <- NA
      next
    }

    n_clust <- length(unique(clusters))
    n_clusters_vec[i] <- n_clust
    all_clusterings[[as.character(res)]] <- clusters

    # Evaluate metric on subsample
    eval_clusters <- clusters[eval_idx]

    if (n_clust < 2 || n_clust >= length(eval_clusters)) {
      cat(" k=", n_clust, "(invalid for metric)\n")
      scores[i] <- -Inf
      next
    }

    scores[i] <- tryCatch({
      if (metric == "silhouette") {
        compute_silhouette(eval_dist, eval_clusters)
      } else {
        compute_calinski_harabasz(eval_matrix, eval_clusters)
      }
    }, error = function(e) {
      cat(" metric error:", e$message)
      -Inf
    })

    cat(" k=", n_clust, " score=", round(scores[i], 4), "\n")
  }

  # Select best
  valid_scores <- !is.na(scores) & is.finite(scores)
  if (!any(valid_scores)) {
    cat("Warning: no valid scores found, returning first clustering\n")
    best_idx <- 1
  } else {
    best_idx <- which(valid_scores)[which.max(scores[valid_scores])]
  }

  best_res <- resolutions[best_idx]
  best_clusters <- all_clusterings[[as.character(best_res)]]

  scores_df <- data.frame(
    resolution = resolutions,
    n_clusters = n_clusters_vec,
    score = scores,
    metric = metric,
    stringsAsFactors = FALSE
  )

  cat("Best resolution:", best_res, "with", length(unique(best_clusters)),
      "clusters (score:", round(scores[best_idx], 4), ")\n")

  list(
    best_clusters = best_clusters,
    best_param = best_res,
    scores_df = scores_df,
    all_clusterings = all_clusterings
  )
}

#' Build consensus clusters from multiple clustering runs
#'
#' Takes a named list of cluster assignment vectors (one per method or
#' parameter setting), builds a co-association matrix, and derives
#' consensus cluster assignments via hierarchical clustering.
#'
#' Input from: run_resolution_sweep() (this file),
#'             run_leiden() (clustering/leiden.R),
#'             run_flowsom() (clustering/flowsom.R),
#'             run_phenograph() (clustering/phenograph.R)
#' Output to:  clusters_to_gating_csv() (clustering/leiden.R),
#'             validate_clusters() (clustering/cluster_validation.R)
#'
#' @param cluster_list Named list of integer vectors (cluster assignments).
#'   Each vector must have the same length (one entry per cell).
#' @param method Consensus method: "majority" (co-association + hclust) or
#'   "voting" (cell-level majority vote)
#' @param n_consensus_clusters Number of clusters to cut from the consensus
#'   dendrogram. If NULL, determined automatically from the co-association matrix.
#' @return Integer vector of consensus cluster assignments
consensus_clusters <- function(cluster_list, method = "majority",
                               n_consensus_clusters = NULL) {
  stopifnot(is.list(cluster_list), length(cluster_list) >= 2)
  stopifnot(method %in% c("majority", "voting"))

  lengths <- sapply(cluster_list, length)
  if (length(unique(lengths)) != 1) {
    stop("All cluster vectors must have the same length. Got: ",
         paste(unique(lengths), collapse = ", "))
  }
  n_cells <- lengths[1]
  n_runs <- length(cluster_list)

  cat("Building consensus from", n_runs, "clustering runs (",
      n_cells, "cells)\n")

  if (method == "voting") {
    # Simple majority vote per cell
    cat("Using majority voting...\n")
    # Build a matrix of cluster labels: cells x runs
    label_matrix <- do.call(cbind, cluster_list)

    consensus <- apply(label_matrix, 1, function(row) {
      tab <- table(row)
      as.integer(names(tab)[which.max(tab)])
    })

    # Relabel sequentially
    consensus <- as.integer(factor(consensus, levels = unique(consensus)))
    cat("Voting consensus:", length(unique(consensus)), "clusters\n")
    return(consensus)
  }

  # method == "majority": co-association matrix approach
  cat("Building co-association matrix...\n")

  # For large datasets, use sparse representation
  if (n_cells > 10000) {
    cat("Using block-wise computation for", n_cells, "cells...\n")
    consensus <- consensus_via_sampling(cluster_list, n_consensus_clusters)
    return(consensus)
  }

  # Co-association matrix: M[i,j] = fraction of runs where i and j are
  # in the same cluster
  coassoc <- matrix(0, nrow = n_cells, ncol = n_cells)

  for (run_name in names(cluster_list)) {
    clusters <- cluster_list[[run_name]]
    for (cl in unique(clusters)) {
      members <- which(clusters == cl)
      if (length(members) > 1) {
        pairs <- expand.grid(i = members, j = members)
        pairs <- pairs[pairs$i < pairs$j, ]
        for (r in seq_len(nrow(pairs))) {
          coassoc[pairs$i[r], pairs$j[r]] <- coassoc[pairs$i[r], pairs$j[r]] + 1
          coassoc[pairs$j[r], pairs$i[r]] <- coassoc[pairs$j[r], pairs$i[r]] + 1
        }
      }
    }
  }

  coassoc <- coassoc / n_runs
  diag(coassoc) <- 1

  # Convert to distance and run hierarchical clustering
  dissim <- as.dist(1 - coassoc)
  hc <- hclust(dissim, method = "average")

  # Determine number of clusters
  if (is.null(n_consensus_clusters)) {
    # Use the median number of clusters across runs
    n_per_run <- sapply(cluster_list, function(x) length(unique(x)))
    n_consensus_clusters <- round(median(n_per_run))
    cat("Auto-selected k =", n_consensus_clusters,
        "(median across runs)\n")
  }

  consensus <- as.integer(cutree(hc, k = n_consensus_clusters))

  cat("Consensus clustering:", length(unique(consensus)), "clusters\n")
  cat("Cluster sizes:", paste(table(consensus), collapse = ", "), "\n")

  consensus
}

#' Consensus clustering via landmark sampling (for large datasets)
#'
#' @param cluster_list Named list of cluster vectors
#' @param n_consensus_clusters Target number of clusters (NULL = auto)
#' @param n_landmarks Number of landmark cells to sample
#' @return Integer vector of consensus cluster assignments
consensus_via_sampling <- function(cluster_list, n_consensus_clusters = NULL,
                                   n_landmarks = 5000) {
  n_cells <- length(cluster_list[[1]])
  n_runs <- length(cluster_list)

  landmark_idx <- sample(n_cells, min(n_landmarks, n_cells))
  n_lm <- length(landmark_idx)

  cat("Computing co-association on", n_lm, "landmarks...\n")

  # Build landmark co-association matrix
  coassoc_lm <- matrix(0, nrow = n_lm, ncol = n_lm)

  for (clusters in cluster_list) {
    lm_clusters <- clusters[landmark_idx]
    for (cl in unique(lm_clusters)) {
      members <- which(lm_clusters == cl)
      if (length(members) > 1) {
        combs <- combn(members, 2)
        for (c in seq_len(ncol(combs))) {
          coassoc_lm[combs[1, c], combs[2, c]] <- coassoc_lm[combs[1, c], combs[2, c]] + 1
          coassoc_lm[combs[2, c], combs[1, c]] <- coassoc_lm[combs[2, c], combs[1, c]] + 1
        }
      }
    }
  }

  coassoc_lm <- coassoc_lm / n_runs
  diag(coassoc_lm) <- 1

  dissim_lm <- as.dist(1 - coassoc_lm)
  hc_lm <- hclust(dissim_lm, method = "average")

  if (is.null(n_consensus_clusters)) {
    n_per_run <- sapply(cluster_list, function(x) length(unique(x)))
    n_consensus_clusters <- round(median(n_per_run))
  }

  lm_labels <- cutree(hc_lm, k = n_consensus_clusters)

  # Assign remaining cells to nearest landmark's cluster
  cat("Assigning remaining cells to landmark clusters...\n")
  consensus <- integer(n_cells)
  consensus[landmark_idx] <- lm_labels

  non_lm_idx <- setdiff(seq_len(n_cells), landmark_idx)
  if (length(non_lm_idx) > 0) {
    # For each non-landmark cell, find the most co-associated landmark
    for (ci in non_lm_idx) {
      # Count co-occurrences with each landmark across runs
      cooccur <- integer(n_lm)
      for (clusters in cluster_list) {
        ci_cluster <- clusters[ci]
        lm_clusters <- clusters[landmark_idx]
        cooccur <- cooccur + as.integer(lm_clusters == ci_cluster)
      }
      nearest_lm <- which.max(cooccur)
      consensus[ci] <- lm_labels[nearest_lm]
    }
  }

  as.integer(consensus)
}

#' Compute mean silhouette width
#'
#' @param dist_obj dist object
#' @param cluster_ids Integer vector of cluster assignments
#' @return Numeric scalar: mean silhouette width
compute_silhouette <- function(dist_obj, cluster_ids) {
  require(cluster)
  sil <- cluster::silhouette(cluster_ids, dist_obj)
  mean(sil[, "sil_width"])
}

#' Compute Calinski-Harabasz index
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param cluster_ids Integer vector of cluster assignments
#' @return Numeric scalar: Calinski-Harabasz index
compute_calinski_harabasz <- function(expr_matrix, cluster_ids) {
  n <- nrow(expr_matrix)
  k <- length(unique(cluster_ids))

  overall_center <- colMeans(expr_matrix)

  # Between-cluster sum of squares
  bgss <- 0
  # Within-cluster sum of squares
  wgss <- 0

  for (cl in unique(cluster_ids)) {
    members <- which(cluster_ids == cl)
    n_cl <- length(members)
    cl_center <- colMeans(expr_matrix[members, , drop = FALSE])
    bgss <- bgss + n_cl * sum((cl_center - overall_center)^2)
    wgss <- wgss + sum(sweep(expr_matrix[members, , drop = FALSE], 2,
                              cl_center)^2)
  }

  if (wgss == 0) return(Inf)
  (bgss / (k - 1)) / (wgss / (n - k))
}
