#!/usr/bin/env Rscript
library(optparse)
library(jsonlite)

option_list <- list(
  make_option("--pop_strings_rds", type = "character",
              help = "pop_strings_list.rds"),
  make_option("--target_fcs", type = "character", default = NULL),
  make_option("--label_dict", type = "character", default = "{}",
              help = "JSON string for label dictionary"),
  make_option("--custom_name", type = "character", default = NULL),
  make_option("--output_dir", type = "character", default = "gating_comparison")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(src_dir, "metadata_qc", "compare_gatings.R"))

pop_strings_list <- readRDS(opt$pop_strings_rds)
label_dict <- as.list(fromJSON(opt$label_dict))
custom_name <- if (!is.null(opt$custom_name)) {
  strsplit(opt$custom_name, ",")[[1]]
} else {
  NULL
}

label_lists <- compare_gatings(
  pop_strings_list = pop_strings_list,
  target_fcs       = opt$target_fcs,
  label_dict       = label_dict,
  custom_name      = custom_name,
  output_dir       = opt$output_dir,
  scale_size       = TRUE
)

saveRDS(label_lists, file.path(opt$output_dir, "label_lists.rds"))
cat("Gating comparison ->", opt$output_dir, "\n")
