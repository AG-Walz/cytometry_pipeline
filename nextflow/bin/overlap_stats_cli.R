#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--label_lists_rds", type = "character",
              help = "label_lists.rds from compare_gatings"),
  make_option("--output_dir", type = "character", default = "overlap_stats")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(src_dir, "metadata_qc", "overlap_statistics.R"))

label_lists <- readRDS(opt$label_lists_rds)

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

kappa_results <- compute_fleiss_kappa(
  label_lists = label_lists,
  output_dir  = opt$output_dir
)
cat("Overlap statistics ->", opt$output_dir, "\n")
