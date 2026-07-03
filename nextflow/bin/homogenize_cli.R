#!/usr/bin/env Rscript
library(optparse)
library(jsonlite)

option_list <- list(
  make_option("--input_rds_dir", type = "character",
              help = "Directory containing per-sample RDS files"),
  make_option("--markers_of_interest", type = "character",
              help = "Comma-separated marker names"),
  make_option("--markername_replacing", type = "character", default = "{}",
              help = "JSON string of marker name replacements"),
  make_option("--output_dir", type = "character", default = "homogenized")
)
opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages(library(flowCore))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(src_dir, "harmonization", "panel_homogenization.R"))

markers <- strsplit(opt$markers_of_interest, ",")[[1]]
replacing <- as.list(fromJSON(opt$markername_replacing))

# Load all flowFrames
rds_files <- sort(list.files(opt$input_rds_dir, pattern = "\\.rds$",
                             full.names = TRUE))
flow_frames <- lapply(rds_files, readRDS)
sample_names <- tools::file_path_sans_ext(basename(rds_files))
names(flow_frames) <- sample_names

# Homogenize
homogenized <- homogenize_flowframes(
  flow_frames         = flow_frames,
  sample_names        = sample_names,
  markers_of_interest = markers,
  markername_replacing = replacing,
  output_dir          = NULL
)

# Write one RDS per sample for re-scatter
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
for (nm in names(homogenized)) {
  saveRDS(homogenized[[nm]], file.path(opt$output_dir, paste0(nm, ".rds")))
}
cat("Homogenized", length(homogenized), "samples ->", opt$output_dir, "\n")
