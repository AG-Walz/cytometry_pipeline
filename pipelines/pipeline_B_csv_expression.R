#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_B_csv_expression.R — CSV Expression → Statistics Pipeline
# =============================================================================
# Flow: Load CSV expression → Marker statistics → Density/Heatmap
# =============================================================================

cat("=== Pipeline B: CSV Expression → Statistics ===\n")
start_time <- Sys.time()

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
if (!is.null(config_file)) {
  source(config_file)
} else {
  source(file.path(dirname(getwd()), "config", "Parameters_template.R"))
}

base_dir <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))

# Step 1: Load CSV expression data
cat("\n--- Step 1: Loading CSV expression data ---\n")
source(file.path(base_dir, "src", "loading", "load_csv_expression.R"))
expr_data <- load_csv_expression(
  csv_folder = params$csv_folder,
  markers = params$markers,
  output_dir = file.path(params$output_dir, "01_loaded")
)

# Step 2: Marker statistics
cat("\n--- Step 2: Marker Statistics ---\n")
source(file.path(base_dir, "src", "analysis", "marker_statistics.R"))
stats_result <- run_marker_statistics(
  data = expr_data$combined,
  markers = params$markers,
  group_col = params$condition_col,
  sample_col = params$sample_col,
  output_dir = file.path(params$output_dir, "02_statistics")
)

# Step 3: Visualization
cat("\n--- Step 3: Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "density_plots.R"))
source(file.path(base_dir, "src", "visualization", "heatmaps.R"))

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline B complete (%.1f minutes) ===\n", elapsed))
