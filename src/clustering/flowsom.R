#!/usr/bin/env Rscript
# =============================================================================
# clustering/flowsom.R — FlowSOM Self-Organizing Map Clustering
# =============================================================================
# FlowSOM meta-clustering for flow cytometry expression data.
# Builds a self-organizing map (SOM) grid and applies consensus
# clustering to identify meta-clusters.
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Run FlowSOM clustering on an expression matrix
#'
#' Builds a self-organizing map (SOM) on the expression data, then applies
#' ConsensusClusterPlus to determine meta-cluster assignments. Returns an
#' integer vector of cluster IDs (one per cell).
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             normalize_expression() (preprocessing/normalization.R)
#' Output to:  clusters_to_gating_csv() (R_shared/cluster_utils.R),
#'             validate_clusters() (clustering/cluster_validation.R),
#'             annotate_clusters() (clustering/cell_annotation.R)
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param n_clusters Target number of meta-clusters
#' @param grid_xdim SOM grid x dimension
#' @param grid_ydim SOM grid y dimension
#' @param seed Random seed for reproducibility
#' @return Integer vector of cluster assignments (length = nrow(expr_matrix))
run_flowsom <- function(expr_matrix, n_clusters = 20,
                        grid_xdim = 10, grid_ydim = 10,
                        seed = 42) {
  require(FlowSOM)
  require(flowCore)

  set.seed(seed)
  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  stopifnot(n_clusters > 1, grid_xdim > 0, grid_ydim > 0)

  cat("Building FlowSOM grid (", grid_xdim, "x", grid_ydim, ")...\n")

  # Convert matrix to flowFrame for FlowSOM input
  ff <- from_mtx_to_ffr(expr_matrix)

  # Build the SOM
  fsom <- tryCatch({
    FlowSOM::FlowSOM(ff,
                     colsToUse = seq_len(ncol(expr_matrix)),
                     xdim = grid_xdim, ydim = grid_ydim,
                     nClus = n_clusters,
                     seed = seed)
  }, error = function(e) {
    cat("FlowSOM with nClus failed:", e$message, "\n")
    cat("Retrying with default ConsensusClusterPlus settings...\n")
    fsom_raw <- FlowSOM::ReadInput(ff, transform = FALSE, scale = FALSE)
    fsom_raw <- FlowSOM::BuildSOM(fsom_raw, colsToUse = seq_len(ncol(expr_matrix)),
                                  xdim = grid_xdim, ydim = grid_ydim)
    fsom_raw <- FlowSOM::BuildMST(fsom_raw)
    meta <- FlowSOM::metaClustering_consensus(fsom_raw$map$codes,
                                               k = n_clusters, seed = seed)
    fsom_raw$metaclustering <- meta
    fsom_raw
  })

  # Extract meta-cluster assignments per cell
  cluster_ids <- tryCatch({
    as.integer(FlowSOM::GetMetaclusters(fsom))
  }, error = function(e) {
    cat("GetMetaclusters failed, extracting manually...\n")
    node_assignments <- fsom$map$mapping[, 1]
    meta <- fsom$metaclustering
    as.integer(meta[node_assignments])
  })

  cat("FlowSOM found", length(unique(cluster_ids)), "meta-clusters\n")
  cat("Cluster sizes:", paste(table(cluster_ids), collapse = ", "), "\n")

  cluster_ids
}

#' Run FlowSOM clustering on a list of flowFrames
#'
#' Extracts expression matrices from flowFrames for the specified markers,
#' concatenates them, runs FlowSOM, and returns cluster assignments with
#' sample-level tracking.
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             homogenize_flowframes() (harmonization/panel_homogenization.R)
#' Output to:  clusters_to_gating_csv() (R_shared/cluster_utils.R),
#'             run_differential_abundance() (analysis/differential_abundance.R)
#'
#' @param flow_frames List of flowFrame objects
#' @param markers_of_interest Character vector of marker names to use for clustering
#' @param n_clusters Target number of meta-clusters
#' @param seed Random seed
#' @return List with components:
#'   \item{cluster_ids}{Integer vector of cluster assignments}
#'   \item{flowsom_obj}{The FlowSOM object for downstream use}
#'   \item{sample_labels}{Character vector mapping each cell to its source sample}
run_flowsom_flowframes <- function(flow_frames, markers_of_interest,
                                   n_clusters = 20, seed = 42) {
  require(FlowSOM)
  require(flowCore)

  stopifnot(is.list(flow_frames), length(flow_frames) > 0)

  sample_names <- if (!is.null(names(flow_frames))) {
    names(flow_frames)
  } else {
    paste0("sample_", seq_along(flow_frames))
  }

  cat("Extracting expression from", length(flow_frames), "flowFrames...\n")

  # Use shared utility to build combined matrix
  combined <- build_combined_expr_matrix(flow_frames, sample_names,
                                          markers_of_interest)
  expr_matrix <- combined$combined_matrix
  sample_labels <- combined$sample_labels

  cat("Running FlowSOM on combined matrix:",
      nrow(expr_matrix), "cells x", ncol(expr_matrix), "markers\n")

  # Determine grid size based on data
  grid_dim <- max(5, min(15, round(sqrt(n_clusters * 2))))

  # Build FlowSOM
  ff <- from_mtx_to_ffr(expr_matrix)
  set.seed(seed)

  fsom <- tryCatch({
    FlowSOM::FlowSOM(ff,
                     colsToUse = seq_len(ncol(expr_matrix)),
                     xdim = grid_dim, ydim = grid_dim,
                     nClus = n_clusters,
                     seed = seed)
  }, error = function(e) {
    cat("FlowSOM failed:", e$message, "\n")
    cat("Falling back to manual SOM + consensus clustering...\n")
    fsom_raw <- FlowSOM::ReadInput(ff, transform = FALSE, scale = FALSE)
    fsom_raw <- FlowSOM::BuildSOM(fsom_raw, colsToUse = seq_len(ncol(expr_matrix)),
                                  xdim = grid_dim, ydim = grid_dim)
    fsom_raw <- FlowSOM::BuildMST(fsom_raw)
    meta <- FlowSOM::metaClustering_consensus(fsom_raw$map$codes,
                                               k = n_clusters, seed = seed)
    fsom_raw$metaclustering <- meta
    fsom_raw
  })

  cluster_ids <- tryCatch({
    as.integer(FlowSOM::GetMetaclusters(fsom))
  }, error = function(e) {
    node_assignments <- fsom$map$mapping[, 1]
    meta <- fsom$metaclustering
    as.integer(meta[node_assignments])
  })

  cat("FlowSOM complete:", length(unique(cluster_ids)), "meta-clusters from",
      length(flow_frames), "samples\n")

  list(
    cluster_ids = cluster_ids,
    flowsom_obj = fsom,
    sample_labels = sample_labels
  )
}
