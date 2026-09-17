#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--expression_csv", type = "character", help = "Expression CSV (cells x markers)"),
  make_option("--output_dir", type = "character", default = "correlation"),
  make_option("--method", type = "character", default = "spearman"),
  make_option("--p_adjust_method", type = "character", default = "BH")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(src_dir, "analysis", "correlation_analysis.R"))

expr_df <- read.csv(opt$expression_csv, check.names = FALSE)
expr_matrix <- as.matrix(expr_df[, vapply(expr_df, is.numeric, logical(1)), drop = FALSE])

run_correlation_analysis(
  expr_matrix     = expr_matrix,
  method          = opt$method,
  output_dir      = opt$output_dir,
  cluster_labels  = NULL,
  p_adjust_method = opt$p_adjust_method
)
cat("Correlation analysis ->", opt$output_dir, "\n")
