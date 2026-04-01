#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--expr_rds", type = "character", help = "expr_data.rds"),
  make_option("--cluster_rds", type = "character", help = "cluster_ids.rds"),
  make_option("--output_dir", type = "character", default = "gating_csvs")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(src_dir, "clustering", "leiden.R"))

expr_data   <- readRDS(opt$expr_rds)
cluster_ids <- readRDS(opt$cluster_rds)

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

clusters_to_gating_csv(
  combined_matrix = expr_data$combined_matrix,
  sample_labels   = expr_data$sample_labels,
  cluster_ids     = cluster_ids,
  output_dir      = opt$output_dir
)
cat("Wrote gating CSVs to:", opt$output_dir, "\n")
