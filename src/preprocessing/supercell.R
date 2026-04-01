#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/supercell.R — SuperCellCyto Integration (Metacell Compression)
# =============================================================================
# Compresses single-cell cytometry data into supercells (metacells) using
# SuperCellCyto. Supports global compression and per-cluster compression
# for finer-grained metacell construction.
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))

#' Compress cells into supercells (metacells) using SuperCellCyto
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R)
#' Output to:  run_leiden/run_flowsom/run_phenograph (clustering/),
#'             plot_umap_from_matrix/plot_tsne_from_matrix (visualization/),
#'             run_differential_abundance/run_differential_expression (analysis/)
#'
#' @param expr_matrix Numeric matrix (cells x markers)
#' @param sample_labels Character vector (one per row)
#' @param markers Character vector of marker names (colnames of expr_matrix to use)
#' @param gam Graining level: higher = fewer supercells (default 20 = ~1/20 compression)
#' @param k_knn Nearest neighbors for graph (default 30)
#' @param n_pc PCs for initial reduction (default 10)
#' @param do_pca Run PCA before supercell construction (default TRUE)
#' @param seed Random seed
#' @return list(supercell_matrix, supercell_labels, cell_map, supercell_sizes)
run_supercell <- function(expr_matrix, sample_labels, markers,
                          gam = 20, k_knn = 30, n_pc = 10,
                          do_pca = TRUE, seed = 42) {
  require(SuperCellCyto)
  require(data.table)

  # --- Validate inputs -------------------------------------------------------
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix or data.frame.")
  }
  if (nrow(expr_matrix) != length(sample_labels)) {
    stop("Number of rows in expr_matrix (", nrow(expr_matrix),
         ") must match length of sample_labels (", length(sample_labels), ").")
  }
  missing_markers <- setdiff(markers, colnames(expr_matrix))
  if (length(missing_markers) > 0) {
    stop("Markers not found in expr_matrix: ",
         paste(missing_markers, collapse = ", "))
  }

  n_cells <- nrow(expr_matrix)
  n_samples <- length(unique(sample_labels))
  cat("SuperCellCyto: compressing", n_cells, "cells from", n_samples,
      "samples (gam =", gam, ", k_knn =", k_knn, ")\n")

  # --- Build data.table with sample column -----------------------------------
  dt <- data.table::as.data.table(expr_matrix[, markers, drop = FALSE])
  dt[, sample := sample_labels]

  # --- Run SuperCellCyto -----------------------------------------------------
  result <- tryCatch({
    SuperCellCyto::runSuperCellCyto(
      dt       = dt,
      markers  = markers,
      sample_colname = "sample",
      gam      = gam,
      k_knn    = k_knn,
      n_pc     = n_pc,
      do_pca   = do_pca,
      seed     = seed
    )
  }, error = function(e) {
    stop("SuperCellCyto failed: ", e$message,
         "\nCheck that markers are numeric and gam/k_knn are appropriate ",
         "for the dataset size.")
  })

  # --- Extract supercell expression matrix -----------------------------------
  sc_expr_dt <- result$supercell_expression_matrix
  supercell_labels <- sc_expr_dt[["sample"]]

  # Drop the sample column to get numeric matrix
  sc_marker_cols <- intersect(markers, names(sc_expr_dt))
  supercell_matrix <- as.matrix(sc_expr_dt[, ..sc_marker_cols])

  # --- Extract cell map ------------------------------------------------------
  cell_map <- result$supercell_cell_map

  # --- Compute supercell sizes -----------------------------------------------
  supercell_sizes <- cell_map[, .N, by = .(SuperCellId)]
  setnames(supercell_sizes, "N", "size")

  # --- Log compression ratio -------------------------------------------------
  n_supercells <- nrow(supercell_matrix)
  compression_ratio <- n_cells / n_supercells
  cat("SuperCellCyto complete:", n_supercells, "supercells created\n")
  cat("Compression ratio:", round(compression_ratio, 1),
      ":1 (", n_cells, "->", n_supercells, ")\n")
  cat("Supercell size range:",
      min(supercell_sizes$size), "-", max(supercell_sizes$size),
      "(median:", median(supercell_sizes$size), ")\n")

  list(
    supercell_matrix = supercell_matrix,
    supercell_labels = supercell_labels,
    cell_map         = cell_map,
    supercell_sizes  = supercell_sizes
  )
}

#' Run SuperCellCyto per cluster (finer-grained compression)
#'
#' Applies supercell compression independently within each cluster, preserving
#' cluster boundaries. Useful when clusters have very different densities or
#' when downstream analysis requires cluster-aware metacells.
#'
#' Input from: run_leiden/run_flowsom/run_phenograph (clustering/),
#'             build_combined_expr_matrix() (R_shared/flowframe_utils.R)
#' Output to:  run_differential_abundance/run_differential_expression (analysis/),
#'             plot_umap_from_matrix (visualization/)
#'
#' @param expr_matrix Numeric matrix (cells x markers)
#' @param cluster_ids Integer vector of cluster assignments
#' @param sample_labels Character vector
#' @param markers Character vector
#' @param gam Graining level per cluster
#' @param seed Random seed
#' @return list(supercell_matrix, supercell_cluster_ids, supercell_labels, cell_map)
run_supercell_per_cluster <- function(expr_matrix, cluster_ids, sample_labels,
                                      markers, gam = 20, seed = 42) {
  require(data.table)

  # --- Validate inputs -------------------------------------------------------
  if (nrow(expr_matrix) != length(cluster_ids)) {
    stop("Number of rows in expr_matrix (", nrow(expr_matrix),
         ") must match length of cluster_ids (", length(cluster_ids), ").")
  }
  if (nrow(expr_matrix) != length(sample_labels)) {
    stop("Number of rows in expr_matrix (", nrow(expr_matrix),
         ") must match length of sample_labels (", length(sample_labels), ").")
  }

  unique_clusters <- sort(unique(cluster_ids))
  n_clusters <- length(unique_clusters)
  cat("Per-cluster SuperCellCyto:", n_clusters, "clusters, gam =", gam, "\n")

  # --- Process each cluster --------------------------------------------------
  sc_matrices <- list()
  sc_cluster_ids <- list()
  sc_labels <- list()
  sc_cell_maps <- list()

  for (cl in unique_clusters) {
    idx <- which(cluster_ids == cl)
    n_cells_cl <- length(idx)

    cat("  Cluster", cl, ":", n_cells_cl, "cells")

    # Skip supercell if cluster is too small for meaningful compression
    if (n_cells_cl < gam * 2) {
      cat(" (too small for compression, keeping as-is)\n")
      sc_matrices[[as.character(cl)]] <- expr_matrix[idx, markers, drop = FALSE]
      sc_cluster_ids[[as.character(cl)]] <- rep(cl, n_cells_cl)
      sc_labels[[as.character(cl)]] <- sample_labels[idx]

      # Build identity cell map
      identity_map <- data.table::data.table(
        CellId = seq_len(n_cells_cl),
        SuperCellId = seq_len(n_cells_cl),
        OriginalIndex = idx
      )
      sc_cell_maps[[as.character(cl)]] <- identity_map
      next
    }

    # Run supercell on this cluster's cells
    cl_result <- tryCatch({
      run_supercell(
        expr_matrix   = expr_matrix[idx, , drop = FALSE],
        sample_labels = sample_labels[idx],
        markers       = markers,
        gam           = gam,
        seed          = seed
      )
    }, error = function(e) {
      warning("SuperCellCyto failed for cluster ", cl, ": ", e$message,
              ". Keeping original cells.")
      NULL
    })

    if (is.null(cl_result)) {
      # Fallback: keep original cells
      sc_matrices[[as.character(cl)]] <- expr_matrix[idx, markers, drop = FALSE]
      sc_cluster_ids[[as.character(cl)]] <- rep(cl, n_cells_cl)
      sc_labels[[as.character(cl)]] <- sample_labels[idx]

      identity_map <- data.table::data.table(
        CellId = seq_len(n_cells_cl),
        SuperCellId = seq_len(n_cells_cl),
        OriginalIndex = idx
      )
      sc_cell_maps[[as.character(cl)]] <- identity_map
    } else {
      n_sc <- nrow(cl_result$supercell_matrix)
      cat("  -> ", n_sc, " supercells\n")
      sc_matrices[[as.character(cl)]] <- cl_result$supercell_matrix
      sc_cluster_ids[[as.character(cl)]] <- rep(cl, n_sc)
      sc_labels[[as.character(cl)]] <- cl_result$supercell_labels

      cl_cell_map <- cl_result$cell_map
      cl_cell_map[, OriginalIndex := idx[CellId]]
      sc_cell_maps[[as.character(cl)]] <- cl_cell_map
    }
  }

  # --- Combine results -------------------------------------------------------
  combined_matrix <- do.call(rbind, sc_matrices)
  combined_cluster_ids <- unlist(sc_cluster_ids, use.names = FALSE)
  combined_labels <- unlist(sc_labels, use.names = FALSE)
  combined_cell_map <- data.table::rbindlist(sc_cell_maps, idcol = "cluster")

  n_total_sc <- nrow(combined_matrix)
  cat("Per-cluster compression complete:", n_total_sc, "total supercells from",
      nrow(expr_matrix), "cells\n")
  cat("Overall compression ratio:",
      round(nrow(expr_matrix) / n_total_sc, 1), ":1\n")

  list(
    supercell_matrix      = combined_matrix,
    supercell_cluster_ids = combined_cluster_ids,
    supercell_labels      = combined_labels,
    cell_map              = combined_cell_map
  )
}
