#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_folder", type = "character", help = "CSV folder"),
  make_option("--output_folder", type = "character", default = "marker_stats"),
  make_option("--markers", type = "character", default = NULL,
              help = "Comma-separated marker names"),
  make_option("--test_type", type = "character", default = "lmer")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$input_folder) && !dir.exists(opt$input_folder)) {
  stop("Input directory not found: ", opt$input_folder)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(src_dir, "analysis", "marker_statistics.R"))

tryCatch({
  markers <- if (!is.null(opt$markers)) {
    strsplit(opt$markers, ",")[[1]]
  } else {
    NULL
  }

  run_marker_statistics(
    input_folder  = opt$input_folder,
    output_folder = opt$output_folder,
    markers       = markers,
    test_type     = opt$test_type
  )
  cat("Marker statistics ->", opt$output_folder, "\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] Marker statistics failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
