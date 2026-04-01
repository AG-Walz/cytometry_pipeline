#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_rds", type = "character", help = "Input flowFrame RDS"),
  make_option("--sample_name", type = "character"),
  make_option("--output_rds", type = "character", default = "downsampled.rds"),
  make_option("--target_cells", type = "integer", default = 10000),
  make_option("--method", type = "character", default = "Consistent"),
  make_option("--seed", type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(src_dir, "preprocessing", "downsampling.R"))

ff <- readRDS(opt$input_rds)
flow_frames <- list(ff)
names(flow_frames) <- opt$sample_name

set.seed(opt$seed)
downsampled <- downsample_flowframes(
  flow_frames  = flow_frames,
  sample_names = opt$sample_name,
  target_cells = opt$target_cells,
  method       = opt$method,
  seed         = opt$seed
)

saveRDS(downsampled[[1]], opt$output_rds)
cat("Downsampled:", opt$sample_name, "from", nrow(ff), "to",
    nrow(downsampled[[1]]), "cells\n")
