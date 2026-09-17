#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--expr_rds", type = "character", help = "expr_data.rds or matrix RDS"),
  make_option("--output_path", type = "character", default = "heatmap.png"),
  make_option("--title", type = "character",
              default = "Expression Heatmap (sampled cells)"),
  make_option("--max_cells", type = "integer", default = 500),
  make_option("--seed", type = "integer", default = 42)
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$expr_rds) && !file.exists(opt$expr_rds)) {
  stop("Input file not found: ", opt$expr_rds)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(src_dir, "visualization", "heatmaps.R"))

tryCatch({
  data <- readRDS(opt$expr_rds)
  # Handle both raw matrix and expr_data list
  if (is.list(data) && !is.null(data$combined_matrix)) {
    mat <- data$combined_matrix
  } else {
    mat <- data
  }

  # Sample cells if needed
  if (nrow(mat) > opt$max_cells) {
    set.seed(opt$seed)
    mat <- mat[sample(nrow(mat), opt$max_cells), ]
  }

  plot_heatmap(
    expr_matrix = t(mat),
    output_path = opt$output_path,
    title       = opt$title
  )
  cat("Heatmap ->", opt$output_path, "\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] Heatmap generation failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
