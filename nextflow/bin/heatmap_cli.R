#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--expr_rds", type = "character", help = "expr_data.rds or matrix RDS"),
  make_option("--output_path", type = "character", default = "heatmap.png"),
  make_option("--title", type = "character",
              default = "Expression Heatmap (sampled cells)"),
  make_option("--max_cells", type = "integer", default = 500)
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(src_dir, "visualization", "heatmaps.R"))

data <- readRDS(opt$expr_rds)
# Handle raw matrix, expr_data list, or a named list of per-sample tables
if (is.list(data) && !is.null(data$combined_matrix)) {
  mat <- data$combined_matrix
} else if (is.list(data) && !is.data.frame(data)) {
  # named list of per-sample expression tables (e.g. from load_csv_expression)
  mat <- do.call(rbind, lapply(data, as.matrix))
} else {
  mat <- as.matrix(data)
}

# Sample cells if needed
if (nrow(mat) > opt$max_cells) {
  set.seed(42)
  mat <- mat[sample(nrow(mat), opt$max_cells), ]
}

plot_heatmap(
  expr_matrix = t(mat),
  output_path = opt$output_path,
  title       = opt$title
)
cat("Heatmap ->", opt$output_path, "\n")
