#!/usr/bin/env Rscript
# =============================================================================
# clustering/cluster_validation.R — Cluster Quality Validation Metrics
# =============================================================================
# Computes internal and external validation metrics for evaluating
# cluster quality: silhouette width, Calinski-Harabasz index, and
# adjusted Rand index (when ground-truth labels are available).
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))

#' Validate clustering quality with internal and external metrics
#'
#' Computes a set of internal validation metrics on a clustering result.
#' When the dataset is large, subsamples for efficiency. If ground-truth
#' labels are provided, also computes the adjusted Rand index (ARI).
#'
#' Input from: run_leiden() (clustering/leiden.R),
#'             run_flowsom() (clustering/flowsom.R),
#'             run_phenograph() (clustering/phenograph.R),
#'             run_dbscan() (clustering/dbscan.R),
#'             consensus_clusters() (clustering/meta_clustering.R)
#' Output to:  run_resolution_sweep() (clustering/meta_clustering.R),
#'             analysis scripts for reporting
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param cluster_ids Integer vector of cluster assignments (length = nrow(expr_matrix))
#' @param true_labels Optional integer or character vector of ground-truth labels
#'   (same length as cluster_ids). When provided, ARI is computed.
#' @param metrics Character vector of metrics to compute. Options:
#'   "silhouette", "calinski_harabasz", "davies_bouldin"
#' @param sample_size Maximum number of cells for metric computation. If
#'   nrow(expr_matrix) > sample_size, a random subsample is drawn.
#' @param seed Random seed for subsampling
#' @return Data frame with columns: metric, value, description
validate_clusters <- function(expr_matrix, cluster_ids,
                              true_labels = NULL,
                              metrics = c("silhouette", "calinski_harabasz"),
                              sample_size = 5000, seed = 42) {
  set.seed(seed)
  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  stopifnot(length(cluster_ids) == nrow(expr_matrix))

  n_cells <- nrow(expr_matrix)
  n_clusters <- length(unique(cluster_ids))

  cat("Validating", n_clusters, "clusters on", n_cells, "cells\n")

  # Subsample if needed
  if (n_cells > sample_size) {
    cat("Subsampling", sample_size, "cells for efficiency\n")
    sample_idx <- sample(n_cells, sample_size)
    eval_matrix <- expr_matrix[sample_idx, , drop = FALSE]
    eval_clusters <- cluster_ids[sample_idx]
    if (!is.null(true_labels)) {
      eval_true <- true_labels[sample_idx]
    }
  } else {
    eval_matrix <- expr_matrix
    eval_clusters <- cluster_ids
    if (!is.null(true_labels)) {
      eval_true <- true_labels
    }
  }

  eval_n_clusters <- length(unique(eval_clusters))
  results <- list()

  # Check basic validity
  if (eval_n_clusters < 2) {
    warning("Fewer than 2 clusters in evaluation set (", eval_n_clusters,
            ") — validation metrics are not meaningful")
    cat("Warning: fewer than 2 clusters in evaluation set — metrics not meaningful\n")
    results[[length(results) + 1]] <- data.frame(
      metric = "n_clusters", value = eval_n_clusters,
      description = "Number of clusters (< 2, metrics not valid)",
      stringsAsFactors = FALSE
    )
    return(do.call(rbind, results))
  }

  if (eval_n_clusters >= nrow(eval_matrix)) {
    warning("As many clusters as cells (", eval_n_clusters,
            " clusters, ", nrow(eval_matrix), " cells) — validation metrics are not meaningful")
    cat("Warning: as many clusters as cells — metrics not meaningful\n")
    results[[length(results) + 1]] <- data.frame(
      metric = "n_clusters", value = eval_n_clusters,
      description = "Number of clusters (>= n_cells, metrics not valid)",
      stringsAsFactors = FALSE
    )
    return(do.call(rbind, results))
  }

  # Always report basic stats
  results[[length(results) + 1]] <- data.frame(
    metric = "n_clusters", value = n_clusters,
    description = "Number of clusters",
    stringsAsFactors = FALSE
  )
  results[[length(results) + 1]] <- data.frame(
    metric = "n_cells", value = n_cells,
    description = "Total number of cells",
    stringsAsFactors = FALSE
  )

  # Silhouette width
  if ("silhouette" %in% metrics) {
    cat("Computing silhouette width...\n")
    sil_value <- tryCatch({
      require(cluster)
      d <- dist(eval_matrix)
      sil <- cluster::silhouette(eval_clusters, d)
      mean(sil[, "sil_width"])
    }, error = function(e) {
      cat("  Silhouette failed:", e$message, "\n")
      NA_real_
    })
    results[[length(results) + 1]] <- data.frame(
      metric = "silhouette",
      value = round(sil_value, 6),
      description = "Mean silhouette width (-1 to 1; higher = better separation)",
      stringsAsFactors = FALSE
    )
    cat("  Silhouette:", round(sil_value, 4), "\n")
  }

  # Calinski-Harabasz index
  if ("calinski_harabasz" %in% metrics) {
    cat("Computing Calinski-Harabasz index...\n")
    ch_value <- tryCatch({
      compute_calinski_harabasz_internal(eval_matrix, eval_clusters)
    }, error = function(e) {
      cat("  Calinski-Harabasz failed:", e$message, "\n")
      NA_real_
    })
    results[[length(results) + 1]] <- data.frame(
      metric = "calinski_harabasz",
      value = round(ch_value, 4),
      description = "Calinski-Harabasz index (higher = better defined clusters)",
      stringsAsFactors = FALSE
    )
    cat("  Calinski-Harabasz:", round(ch_value, 4), "\n")
  }

  # Davies-Bouldin index
  if ("davies_bouldin" %in% metrics) {
    cat("Computing Davies-Bouldin index...\n")
    db_value <- tryCatch({
      compute_davies_bouldin(eval_matrix, eval_clusters)
    }, error = function(e) {
      cat("  Davies-Bouldin failed:", e$message, "\n")
      NA_real_
    })
    results[[length(results) + 1]] <- data.frame(
      metric = "davies_bouldin",
      value = round(db_value, 6),
      description = "Davies-Bouldin index (lower = better separation)",
      stringsAsFactors = FALSE
    )
    cat("  Davies-Bouldin:", round(db_value, 4), "\n")
  }

  # Adjusted Rand Index (external metric, requires true labels)
  if (!is.null(true_labels)) {
    cat("Computing Adjusted Rand Index...\n")
    ari_value <- tryCatch({
      require(mclust)
      mclust::adjustedRandIndex(eval_clusters, eval_true)
    }, error = function(e) {
      cat("  ARI failed:", e$message, "\n")
      cat("  Falling back to manual ARI computation...\n")
      tryCatch({
        compute_ari_manual(eval_clusters, eval_true)
      }, error = function(e2) {
        cat("  Manual ARI also failed:", e2$message, "\n")
        NA_real_
      })
    })
    results[[length(results) + 1]] <- data.frame(
      metric = "adjusted_rand_index",
      value = round(ari_value, 6),
      description = "Adjusted Rand Index vs true labels (-1 to 1; 1 = perfect match)",
      stringsAsFactors = FALSE
    )
    cat("  ARI:", round(ari_value, 4), "\n")
  }

  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL

  cat("Validation complete:", nrow(result_df), "metrics computed\n")

  result_df
}

#' Compute Calinski-Harabasz index (internal helper)
#'
#' @param expr_matrix Numeric matrix
#' @param cluster_ids Integer vector of cluster IDs
#' @return Numeric scalar: Calinski-Harabasz index
compute_calinski_harabasz_internal <- function(expr_matrix, cluster_ids) {
  n <- nrow(expr_matrix)
  k <- length(unique(cluster_ids))

  overall_center <- colMeans(expr_matrix)

  bgss <- 0  # between-group sum of squares
  wgss <- 0  # within-group sum of squares

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

#' Compute Davies-Bouldin index (internal helper)
#'
#' @param expr_matrix Numeric matrix
#' @param cluster_ids Integer vector of cluster IDs
#' @return Numeric scalar: Davies-Bouldin index
compute_davies_bouldin <- function(expr_matrix, cluster_ids) {
  unique_clusters <- unique(cluster_ids)
  k <- length(unique_clusters)

  # Compute cluster centroids and average within-cluster distances
  centroids <- list()
  avg_dists <- numeric(k)

  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    members <- which(cluster_ids == cl)
    cl_data <- expr_matrix[members, , drop = FALSE]
    centroids[[i]] <- colMeans(cl_data)
    # Average distance from each point to centroid
    diffs <- sweep(cl_data, 2, centroids[[i]])
    avg_dists[i] <- mean(sqrt(rowSums(diffs^2)))
  }

  # For each cluster, find the worst-case ratio
  db_ratios <- numeric(k)
  for (i in seq_len(k)) {
    max_ratio <- 0
    for (j in seq_len(k)) {
      if (i != j) {
        centroid_dist <- sqrt(sum((centroids[[i]] - centroids[[j]])^2))
        if (centroid_dist > 0) {
          ratio <- (avg_dists[i] + avg_dists[j]) / centroid_dist
          max_ratio <- max(max_ratio, ratio)
        }
      }
    }
    db_ratios[i] <- max_ratio
  }

  mean(db_ratios)
}

#' Compute Adjusted Rand Index manually (fallback when mclust unavailable)
#'
#' @param labels1 Integer/character vector of cluster assignments
#' @param labels2 Integer/character vector of cluster assignments
#' @return Numeric scalar: ARI value
compute_ari_manual <- function(labels1, labels2) {
  n <- length(labels1)
  stopifnot(length(labels2) == n)

  # Build contingency table
  cont <- table(labels1, labels2)

  # Sum of choose(n_ij, 2)
  sum_nij2 <- sum(choose(cont, 2))

  # Row and column marginal sums
  a <- rowSums(cont)
  b <- colSums(cont)
  sum_a2 <- sum(choose(a, 2))
  sum_b2 <- sum(choose(b, 2))

  n_choose_2 <- choose(n, 2)

  expected <- sum_a2 * sum_b2 / n_choose_2
  max_index <- 0.5 * (sum_a2 + sum_b2)

  if (max_index == expected) return(0)

  (sum_nij2 - expected) / (max_index - expected)
}
