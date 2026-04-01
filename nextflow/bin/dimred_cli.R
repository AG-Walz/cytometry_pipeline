#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--expr_rds", type = "character", help = "expr_data.rds"),
  make_option("--cluster_rds", type = "character", help = "cluster_ids.rds"),
  make_option("--output_dir", type = "character", default = "dimred"),
  make_option("--method", type = "character", default = "umap"),
  make_option("--n_neighbors", type = "integer", default = 15),
  make_option("--umap_min_dist", type = "double", default = 0.1),
  make_option("--tsne_perplexity", type = "integer", default = 30),
  make_option("--tsne_max_iter", type = "integer", default = 1000),
  make_option("--output_formats", type = "character", default = "png"),
  make_option("--plot_width", type = "integer", default = 10),
  make_option("--plot_height", type = "integer", default = 8),
  make_option("--seed", type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "dimred_dispatch.R"))

expr_data   <- readRDS(opt$expr_rds)
cluster_ids <- readRDS(opt$cluster_rds)

# Parse output formats
fmt_parts <- strsplit(opt$output_formats, ",")[[1]]
output_formats <- list(
  pdf = "pdf" %in% fmt_parts,
  png = "png" %in% fmt_parts,
  svg = "svg" %in% fmt_parts
)

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

run_dimred(
  expr_matrix    = expr_data$combined_matrix,
  method         = opt$method,
  sample_labels  = expr_data$sample_labels,
  cluster_labels = cluster_ids,
  output_dir     = opt$output_dir,
  params         = list(
    n_neighbors     = opt$n_neighbors,
    umap_min_dist   = opt$umap_min_dist,
    tsne_perplexity = opt$tsne_perplexity,
    tsne_max_iter   = opt$tsne_max_iter,
    output_formats  = output_formats,
    plot_width      = opt$plot_width,
    plot_height     = opt$plot_height,
    seed            = opt$seed
  )
)
cat("Dimensionality reduction ->", opt$output_dir, "\n")
