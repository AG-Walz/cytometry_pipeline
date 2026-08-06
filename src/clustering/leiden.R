#!/usr/bin/env Rscript
# =============================================================================
# clustering/leiden.R — Leiden Community Detection Clustering
# =============================================================================
# Merged from: R_scripts/04_gating/Leiden.R + Clustering.R
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
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' required for Leiden clustering. Install with: install.packages('Matrix')")
  }
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
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R)
#' Output to:  clusters_to_gating_csv() (R_shared/cluster_utils.R)
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

# clusters_to_gating_csv() lives in src/R_shared/cluster_utils.R
