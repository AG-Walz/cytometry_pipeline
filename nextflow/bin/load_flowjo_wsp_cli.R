#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--wsp_paths", type = "character",
              help = "Comma-separated WSP file paths"),
  make_option("--wsp_names", type = "character",
              help = "Comma-separated names for workspaces"),
  make_option("--output_rds", type = "character", default = "pop_strings.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(src_dir, "loading", "load_flowjo_wsp.R"))

wsp_paths <- strsplit(opt$wsp_paths, ",")[[1]]
wsp_names <- strsplit(opt$wsp_names, ",")[[1]]

pop_strings_list <- load_multiple_workspaces(
  wsp_paths = wsp_paths,
  names     = wsp_names
)

saveRDS(pop_strings_list, opt$output_rds)
cat("Loaded", length(pop_strings_list), "workspaces\n")
