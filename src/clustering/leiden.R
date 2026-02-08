#!/usr/bin/env Rscript
# =============================================================================
# clustering/leiden.R — Leiden Community Detection Clustering
# =============================================================================
# Merged from: R_scripts/04_gating/Leiden.R + Clustering.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Run Leiden clustering on expression data
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param n_neighbors Number of nearest neighbors for graph construction
#' @param resolution Leiden resolution parameter (higher = more clusters)
#' @param seed Random seed
#' @return Integer vector of cluster assignments
run_leiden <- function(expr_matrix, n_neighbors = 30,
                       resolution = 1.0, seed = 42) {
  require(leiden)

  set.seed(seed)
  cat("Building KNN graph (k =", n_neighbors, ")...\n")

  # Build adjacency matrix using uwot's nearest neighbors
  nn <- uwot::umap(expr_matrix, n_neighbors = n_neighbors,
                    ret_nn = TRUE, n_components = 2)

  # Create adjacency matrix from nearest neighbors
  n_cells <- nrow(expr_matrix)
  adj_matrix <- Matrix::sparseMatrix(
    i = rep(1:n_cells, each = n_neighbors),
    j = as.vector(nn$nn$euclidean$idx),
    x = 1,
    dims = c(n_cells, n_cells)
  )
  # Make symmetric
  adj_matrix <- adj_matrix + Matrix::t(adj_matrix)
  adj_matrix@x[adj_matrix@x > 1] <- 1

  cat("Running Leiden clustering (resolution =", resolution, ")...\n")
  clusters <- leiden::leiden(adj_matrix, resolution_parameter = resolution)

  cat("Found", length(unique(clusters)), "clusters\n")
  cat("Cluster sizes:", paste(table(clusters), collapse = ", "), "\n")

  clusters
}

#' Run Leiden clustering and save results with UMAP coordinates
#'
#' @param expr_matrix Numeric matrix
#' @param sample_labels Character vector of sample labels per cell
#' @param output_dir Output directory
#' @param n_neighbors KNN neighbors
#' @param resolution Leiden resolution
#' @param seed Random seed
#' @return Data frame with UMAP1, UMAP2, cluster, sample columns
run_leiden_with_umap <- function(expr_matrix, sample_labels = NULL,
                                  output_dir = NULL,
                                  n_neighbors = 30,
                                  resolution = 1.0,
                                  seed = 42) {
  set.seed(seed)

  # Run UMAP
  cat("Computing UMAP...\n")
  umap_coords <- run_umap(expr_matrix, n_neighbors = n_neighbors, seed = seed)

  # Run Leiden
  clusters <- run_leiden(expr_matrix, n_neighbors = n_neighbors,
                         resolution = resolution, seed = seed)

  # Build result data frame
  result <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    cluster = clusters
  )
  if (!is.null(sample_labels)) {
    result$sample <- sample_labels
  }

  # Save if output dir specified
  if (!is.null(output_dir)) {
    make_outdir(output_dir)
    write.csv(result, file.path(output_dir, "leiden_clusters.csv"),
              row.names = FALSE)
    cat("Results saved to:", output_dir, "\n")
  }

  result
}
