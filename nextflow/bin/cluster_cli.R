#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_rds", type = "character", help = "expr_data.rds"),
  make_option("--output_rds", type = "character", default = "cluster_ids.rds"),
  make_option("--method", type = "character", default = "leiden"),
  make_option("--n_neighbors", type = "integer", default = 15),
  make_option("--resolution", type = "double", default = 0.5),
  make_option("--n_clusters", type = "integer", default = 20),
  make_option("--grid_xdim", type = "integer", default = 10),
  make_option("--grid_ydim", type = "integer", default = 10),
  make_option("--dbscan_eps", type = "double", default = NULL),
  make_option("--dbscan_min_pts", type = "integer", default = 5),
  make_option("--seed", type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$input_rds) && !file.exists(opt$input_rds)) {
  stop("Input file not found: ", opt$input_rds)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "clustering_dispatch.R"))

tryCatch({
  expr_data <- readRDS(opt$input_rds)

  cluster_ids <- run_clustering(
    expr_matrix = expr_data$combined_matrix,
    method      = opt$method,
    params      = list(
      n_neighbors    = opt$n_neighbors,
      resolution     = opt$resolution,
      n_clusters     = opt$n_clusters,
      grid_xdim      = opt$grid_xdim,
      grid_ydim      = opt$grid_ydim,
      k              = opt$n_neighbors,
      dbscan_eps     = opt$dbscan_eps,
      dbscan_min_pts = opt$dbscan_min_pts,
      seed           = opt$seed
    )
  )

  saveRDS(cluster_ids, opt$output_rds)
  cat("Clustered:", length(unique(cluster_ids)), "clusters from",
      length(cluster_ids), "cells\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] Clustering failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
