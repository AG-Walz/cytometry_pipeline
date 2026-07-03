#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_csv", type = "character",
              help = "Expression CSV with condition and sample columns"),
  make_option("--output_dir", type = "character", default = "effect_sizes"),
  make_option("--condition_col", type = "character", default = "condition"),
  make_option("--value_cols", type = "character", default = NULL),
  make_option("--methods", type = "character", default = "cohens_d,cliffs_delta")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(src_dir, "analysis", "effect_sizes.R"))

value_cols <- if (!is.null(opt$value_cols)) {
  strsplit(opt$value_cols, ",")[[1]]
} else {
  NULL
}
methods <- strsplit(opt$methods, ",")[[1]]

run_effect_sizes(
  input_csv     = opt$input_csv,
  output_dir    = opt$output_dir,
  condition_col = opt$condition_col,
  value_cols    = value_cols,
  methods       = methods
)
cat("Effect sizes ->", opt$output_dir, "\n")
