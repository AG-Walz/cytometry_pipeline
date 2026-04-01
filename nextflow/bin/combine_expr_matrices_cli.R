#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_rds", type = "character",
              help = "expr_data_list.rds from load_csv_expression"),
  make_option("--markers", type = "character", default = NULL),
  make_option("--max_cells", type = "integer", default = 500),
  make_option("--output_rds", type = "character", default = "combined_matrix.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(src_dir, "analysis", "expression_matrices.R"))

expr_data_list <- readRDS(opt$input_rds)
markers <- if (!is.null(opt$markers)) strsplit(opt$markers, ",")[[1]] else NULL

combined <- combine_expression_matrices(
  expr_data_list = expr_data_list,
  markers        = markers,
  max_cells      = opt$max_cells
)

saveRDS(combined, opt$output_rds)
cat("Combined matrix:", nrow(combined), "x", ncol(combined), "\n")
