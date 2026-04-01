#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--data_dir", type = "character", help = "Dir with gating CSVs"),
  make_option("--conditions_csv", type = "character"),
  make_option("--output_dir", type = "character", default = "enrichment"),
  make_option("--cluster_regex", type = "character",
              default = "Clustering annotated"),
  make_option("--condition_levels", type = "character",
              default = "Control,Treatment"),
  make_option("--stat_method", type = "character", default = "pooled"),
  make_option("--p_adjust_method", type = "character", default = "BH"),
  make_option("--seed", type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(src_dir, "analysis", "cluster_enrichment.R"))

cond_levels <- strsplit(opt$condition_levels, ",")[[1]]

run_cluster_enrichment(
  data_dir         = opt$data_dir,
  conditions_csv   = opt$conditions_csv,
  output_dir       = opt$output_dir,
  cluster_regex    = opt$cluster_regex,
  condition_levels = cond_levels,
  stat_method      = opt$stat_method,
  p_adjust_method  = opt$p_adjust_method,
  seed             = opt$seed
)
cat("Enrichment results ->", opt$output_dir, "\n")
