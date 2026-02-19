#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_B_csv_expression.R — CSV Expression -> Statistics Pipeline
# =============================================================================
# Flow: Load CSV expression -> Marker statistics -> Density/Heatmap
# =============================================================================

cat("=== Pipeline B: CSV Expression -> Statistics ===\n")
start_time <- Sys.time()

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
if (!is.null(config_file)) {
  source(config_file)
} else {
  source(file.path(dirname(getwd()), "config", "pipeline_B_config.R"))
}

base_dir <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))

# Step 1: Load CSV expression data
cat("\n--- Step 1: Loading CSV expression data ---\n")
source(file.path(base_dir, "src", "loading", "load_csv_expression.R"))
expr_data <- load_csv_expression(
  input_dir = params$csv_folder
)
cat("Loaded", length(expr_data), "expression matrices\n")

# Step 2: Marker statistics
cat("\n--- Step 2: Marker Statistics ---\n")
source(file.path(base_dir, "src", "analysis", "marker_statistics.R"))
stats_result <- run_marker_statistics(
  input_folder = params$csv_folder,
  output_folder = file.path(params$output_dir, "01_statistics"),
  markers      = params$markers,
  test_type    = params$test_type
)

# Step 3: Visualization
cat("\n--- Step 3: Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "heatmaps.R"))

# Build a combined expression matrix for heatmap
combined_list <- lapply(names(expr_data), function(nm) {
  dt <- expr_data[[nm]]
  marker_cols <- intersect(params$markers, names(dt))
  if (length(marker_cols) > 0) {
    as.matrix(dt[, ..marker_cols])
  } else {
    NULL
  }
})
combined_list <- Filter(Negate(is.null), combined_list)
if (length(combined_list) > 0) {
  combined_matrix <- do.call(rbind, combined_list)

  # Mean expression per marker for a summary heatmap
  viz_dir <- file.path(params$output_dir, "02_visualization")
  make_outdir(viz_dir)
  plot_heatmap(
    expr_matrix = t(combined_matrix[sample(nrow(combined_matrix), min(500, nrow(combined_matrix))), ]),
    output_path = file.path(viz_dir, "expression_heatmap.png"),
    title       = "Expression Heatmap (sampled cells)"
  )
}

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline B complete (%.1f minutes) ===\n", elapsed))
