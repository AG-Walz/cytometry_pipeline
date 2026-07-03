#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--flowframes_rds", type = "character",
              help = "RDS file with a list of flowFrame objects"),
  make_option("--experiment_info_csv", type = "character",
              help = "CSV with sample_id, group_id, patient_id columns"),
  make_option("--markers", type = "character", default = NULL),
  make_option("--output_dir", type = "character", default = "diffcyt"),
  make_option("--analysis_type", type = "character", default = "DA"),
  make_option("--method_DA", type = "character", default = "diffcyt-DA-edgeR"),
  make_option("--method_DS", type = "character", default = "diffcyt-DS-limma"),
  make_option("--n_clusters", type = "integer", default = 20),
  make_option("--seed", type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(src_dir, "analysis", "diffcyt_analysis.R"))

markers <- if (!is.null(opt$markers)) {
  strsplit(opt$markers, ",")[[1]]
} else {
  NULL
}

flow_frames     <- readRDS(opt$flowframes_rds)
experiment_info <- read.csv(opt$experiment_info_csv)

run_diffcyt_and_save(
  flow_frames         = flow_frames,
  experiment_info     = experiment_info,
  markers_of_interest = markers,
  output_dir          = opt$output_dir,
  analysis_type       = opt$analysis_type,
  method_DA           = opt$method_DA,
  method_DS           = opt$method_DS,
  n_clusters          = opt$n_clusters,
  seed                = opt$seed
)
cat("Diffcyt results ->", opt$output_dir, "\n")
