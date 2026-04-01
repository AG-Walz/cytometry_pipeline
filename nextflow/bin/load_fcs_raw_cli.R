#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--fcs_file", type = "character", help = "Path to single FCS file"),
  make_option("--sample_name", type = "character", help = "Sample name"),
  make_option("--output_rds", type = "character", default = "flowframe.rds"),
  make_option("--truncation", type = "logical", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "fcs_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

suppressPackageStartupMessages(library(flowCore))

ff <- read.FCS(opt$fcs_file, truncate_max_range = opt$truncation)

saveRDS(ff, opt$output_rds)
cat("Loaded:", opt$sample_name, "->", opt$output_rds,
    "(", nrow(ff), "events )\n")
