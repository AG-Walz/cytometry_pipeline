#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_dir", type = "character", help = "Dir with per-cluster CSVs"),
  make_option("--conditions_csv", type = "character"),
  make_option("--output_dir", type = "character", default = "diff_expression"),
  make_option("--cluster_regex", type = "character",
              default = "Clustering annotated"),
  make_option("--markers", type = "character", default = NULL),
  make_option("--test_type", type = "character", default = "lmer"),
  make_option("--p_adjust_method", type = "character", default = "BH"),
  make_option("--seed", type = "integer", default = 123)
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(src_dir, "analysis", "differential_expression.R"))

markers <- if (!is.null(opt$markers)) {
  strsplit(opt$markers, ",")[[1]]
} else {
  NULL
}

run_differential_expression(
  input_dir       = opt$input_dir,
  conditions_csv  = opt$conditions_csv,
  output_dir      = opt$output_dir,
  cluster_regex   = opt$cluster_regex,
  markers         = markers,
  test_type       = opt$test_type,
  p_adjust_method = opt$p_adjust_method,
  seed            = opt$seed
)
cat("Differential expression results ->", opt$output_dir, "\n")
