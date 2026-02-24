#!/usr/bin/env Rscript
# =============================================================================
# clustering/phenograph.R — PhenoGraph Clustering via KNN + Louvain
# =============================================================================
# Implements PhenoGraph-style clustering: builds a K-nearest-neighbor
# graph using uwot, converts to an igraph object, and runs Louvain
# community detection to identify cell populations.
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Run PhenoGraph clustering on expression data
#'
#' Builds a shared nearest-neighbor (SNN) graph from a K-nearest-neighbor
#' search using uwot, then applies Louvain community detection via igraph
#' to identify cell clusters.
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             normalize_expression() (preprocessing/normalization.R)
#' Output to:  clusters_to_gating_csv() (clustering/leiden.R),
#'             validate_clusters() (clustering/cluster_validation.R),
#'             annotate_clusters() (clustering/cell_annotation.R)
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param k Number of nearest neighbors for graph construction
#' @param seed Random seed for reproducibility
#' @return Integer vector of cluster assignments (length = nrow(expr_matrix))
run_phenograph <- function(expr_matrix, k = 30, seed = 42) {
  require(igraph)
  require(uwot)

  set.seed(seed)
  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  stopifnot(k > 1)

  n_cells <- nrow(expr_matrix)
  if (k >= n_cells) {
    cat("Warning: k (", k, ") >= n_cells (", n_cells, "), reducing k to",
        n_cells - 1, "\n")
    k <- n_cells - 1
  }

  cat("Building KNN graph (k =", k, ") for", n_cells, "cells...\n")

  # Find k nearest neighbors using uwot
  nn <- tryCatch({
    uwot::umap(expr_matrix, n_neighbors = k,
               ret_nn = TRUE, n_components = 2,
               n_epochs = 0)
  }, error = function(e) {
    cat("uwot KNN failed:", e$message, "\n")
    cat("Falling back to FNN::get.knn...\n")
    require(FNN)
    knn_result <- FNN::get.knn(expr_matrix, k = k)
    list(nn = list(euclidean = list(
      idx = knn_result$nn.index,
      dist = knn_result$nn.dist
    )))
  })

  nn_idx <- nn$nn$euclidean$idx
  nn_dist <- nn$nn$euclidean$dist

  cat("Building shared nearest-neighbor graph...\n")

  # Build SNN graph: edge weight = number of shared neighbors / k
  # Use Jaccard-weighted edges between mutual neighbors
  edge_list <- vector("list", n_cells)
  edge_weights <- vector("list", n_cells)

  for (i in seq_len(n_cells)) {
    neighbors_i <- nn_idx[i, ]
    for (j_pos in seq_len(k)) {
      j <- neighbors_i[j_pos]
      if (j > i) {
        neighbors_j <- nn_idx[j, ]
        shared <- length(intersect(neighbors_i, neighbors_j))
        if (shared > 0) {
          jaccard <- shared / (2 * k - shared)
          edge_list[[i]] <- c(edge_list[[i]], i, j)
          edge_weights[[i]] <- c(edge_weights[[i]], jaccard)
        }
      }
    }
  }

  all_edges <- unlist(edge_list)
  all_weights <- unlist(edge_weights)

  if (length(all_edges) == 0) {
    cat("Warning: no shared neighbors found, falling back to direct KNN edges\n")
    # Use direct KNN edges with distance-based weights
    rows <- rep(seq_len(n_cells), each = k)
    cols <- as.vector(nn_idx)
    all_edges <- as.vector(rbind(rows, cols))
    all_weights <- 1 / (1 + as.vector(nn_dist))
  }

  cat("Constructing igraph object...\n")
  g <- igraph::make_undirected_graph(all_edges, n = n_cells)
  igraph::E(g)$weight <- all_weights

  # Simplify (remove multi-edges and self-loops, keeping max weight)
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "max"))

  cat("Running Louvain community detection...\n")
  communities <- tryCatch({
    igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  }, error = function(e) {
    cat("Weighted Louvain failed:", e$message, "\n")
    cat("Retrying without weights...\n")
    igraph::cluster_louvain(g)
  })

  cluster_ids <- as.integer(igraph::membership(communities))

  cat("PhenoGraph found", length(unique(cluster_ids)), "clusters\n")
  cat("Modularity:", round(igraph::modularity(communities), 4), "\n")
  cat("Cluster sizes:", paste(table(cluster_ids), collapse = ", "), "\n")

  cluster_ids
}
