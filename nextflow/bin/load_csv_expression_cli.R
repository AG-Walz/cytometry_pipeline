#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_dir", type = "character", help = "CSV folder"),
  make_option("--output_rds", type = "character", default = "expr_data_list.rds")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$input_dir) && !dir.exists(opt$input_dir)) {
  stop("Input directory not found: ", opt$input_dir)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(src_dir, "loading", "load_csv_expression.R"))

tryCatch({
  expr_data <- load_csv_expression(input_dir = opt$input_dir)

  saveRDS(expr_data, opt$output_rds)
  cat("Loaded", length(expr_data), "expression matrices ->", opt$output_rds, "\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] CSV expression loading failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
