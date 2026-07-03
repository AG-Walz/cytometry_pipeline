#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_rds_dir", type = "character",
              help = "Directory containing per-sample RDS flowFrames"),
  make_option("--markers_of_interest", type = "character"),
  make_option("--output_rds", type = "character", default = "expr_data.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages(library(flowCore))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

markers <- strsplit(opt$markers_of_interest, ",")[[1]]

# Load all flowFrames
rds_files <- sort(list.files(opt$input_rds_dir, pattern = "\\.rds$",
                             full.names = TRUE))
flow_frames <- lapply(rds_files, readRDS)
sample_names <- tools::file_path_sans_ext(basename(rds_files))
names(flow_frames) <- sample_names

expr_data <- build_combined_expr_matrix(flow_frames, sample_names, markers)

saveRDS(expr_data, opt$output_rds)
cat("Built expression matrix:", nrow(expr_data$combined_matrix), "cells x",
    ncol(expr_data$combined_matrix), "markers\n")
