#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_dir", type = "character", help = "Dir with gating CSVs"),
  make_option("--output_dir", type = "character", default = "abundance"),
  make_option("--cluster_regex", type = "character",
              default = "Clustering annotated"),
  make_option("--exclude_patterns", type = "character", default = ""),
  make_option("--metadata_file", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$input_dir) && !dir.exists(opt$input_dir)) {
  stop("Input directory not found: ", opt$input_dir)
}
if (!is.null(opt$metadata_file) && !file.exists(opt$metadata_file)) {
  stop("Metadata file not found: ", opt$metadata_file)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(src_dir, "analysis", "abundance_plotting.R"))

tryCatch({
  exclude <- if (nchar(opt$exclude_patterns) > 0) {
    strsplit(opt$exclude_patterns, ",")[[1]]
  } else {
    c()
  }

  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

  tracking_matrix <- build_tracking_matrix(
    input_dir        = opt$input_dir,
    cluster_regex    = opt$cluster_regex,
    exclude_patterns = exclude,
    output_dir       = opt$output_dir
  )
  plot_abundance(
    output_dir      = opt$output_dir,
    tracking_matrix = tracking_matrix,
    metadata_file   = opt$metadata_file
  )
  cat("Abundance plots ->", opt$output_dir, "\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] Abundance plotting failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
