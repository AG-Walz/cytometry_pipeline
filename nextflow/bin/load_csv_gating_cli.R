#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_fcs_path", type = "character"),
  make_option("--input_csv_path", type = "character"),
  make_option("--output_dir", type = "character", default = "loaded_gating"),
  make_option("--gating_pattern", type = "character",
              default = "\\(Clustering annotated\\)$|^Clustering raw$"),
  make_option("--processing_mode", type = "character", default = "manual")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "fcs_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(src_dir, "loading", "load_csv_gating.R"))

gating_data <- load_csv_gating(
  input_fcs_path  = opt$input_fcs_path,
  input_csv_path  = opt$input_csv_path,
  output_path     = opt$output_dir,
  gating_pattern  = opt$gating_pattern,
  processing_mode = opt$processing_mode
)

saveRDS(gating_data, file.path(opt$output_dir, "gating_data.rds"))
cat("Loaded gating data ->", opt$output_dir, "\n")
