#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--split_result_rds", type = "character",
              help = "split_result.rds with tracking matrix"),
  make_option("--output_dir", type = "character", default = "barplots")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$split_result_rds) && !file.exists(opt$split_result_rds)) {
  stop("Input file not found: ", opt$split_result_rds)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(src_dir, "visualization", "barplots.R"))

tryCatch({
  split_result <- readRDS(opt$split_result_rds)

  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(split_result$tracking_matrix) &&
      nrow(split_result$tracking_matrix) > 0) {
    long_data <- tracking_matrix_to_long(split_result$tracking_matrix)
    plot_composition_barplot(
      data       = long_data,
      output_dir = opt$output_dir,
      x_col      = "cluster",
      fill_col   = "file",
      title      = "Cell Composition by Cluster"
    )
  }
  cat("Barplots ->", opt$output_dir, "\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] Barplot generation failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
