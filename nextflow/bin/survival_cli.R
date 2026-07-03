#!/usr/bin/env Rscript
library(optparse)
suppressPackageStartupMessages(library(data.table))

option_list <- list(
  make_option("--abundance_csv", type = "character",
              help = "CSV: sample rows x cluster columns (proportions)"),
  make_option("--clinical_csv", type = "character",
              help = "CSV with sample_id, time, event columns"),
  make_option("--output_dir", type = "character", default = "survival"),
  make_option("--method", type = "character", default = "both"),
  make_option("--split_method", type = "character", default = "median"),
  make_option("--p_adjust_method", type = "character", default = "BH")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(src_dir, "analysis", "survival_analysis.R"))

abundance <- as.data.frame(fread(opt$abundance_csv))

run_survival_analysis(
  abundance_matrix = abundance,
  clinical_csv     = opt$clinical_csv,
  output_dir       = opt$output_dir,
  method           = opt$method,
  split_method     = opt$split_method,
  p_adjust_method  = opt$p_adjust_method
)
cat("Survival analysis results ->", opt$output_dir, "\n")
