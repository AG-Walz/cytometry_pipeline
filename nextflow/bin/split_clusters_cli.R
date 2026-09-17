#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_dir", type = "character", help = "Dir with gating CSVs"),
  make_option("--output_dir", type = "character", default = "split_clusters"),
  make_option("--cluster_regex", type = "character",
              default = "Clustering annotated"),
  make_option("--markers_of_interest", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$input_dir) && !dir.exists(opt$input_dir)) {
  stop("Input directory not found: ", opt$input_dir)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(src_dir, "analysis", "cluster_splitting.R"))

tryCatch({
  markers <- if (!is.null(opt$markers_of_interest)) {
    strsplit(opt$markers_of_interest, ",")[[1]]
  } else {
    NULL
  }

  split_result <- split_by_cluster(
    input_dir           = opt$input_dir,
    output_dir          = opt$output_dir,
    cluster_regex       = opt$cluster_regex,
    markers_of_interest = markers
  )

  # Save tracking matrix for downstream use
  if (!is.null(split_result$tracking_matrix)) {
    saveRDS(split_result, file.path(opt$output_dir, "split_result.rds"))
  }
  cat("Split clusters ->", opt$output_dir, "\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] Cluster splitting failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
