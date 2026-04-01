#!/usr/bin/env Rscript
# =============================================================================
# Pipeline B: CSV Expression -> Statistics Pipeline
# =============================================================================
# Flow: Load CSV expression -> Marker statistics -> Density/Heatmap
# =============================================================================

# --- Init (config + shared modules) ---
base_dir   <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "pipeline_utils.R"))
env <- init_pipeline("Pipeline B: CSV Expression -> Statistics",
                     "config/pipeline_B_config.R")

# Validate required parameters
if (is.null(params$input_dir) && is.null(params$fcs_folder)) {
  stop("input_dir (or fcs_folder) is required in config")
}
if (is.null(params$output_dir)) {
  stop("output_dir is required in config")
}
# Validate output directory can be created
tryCatch(make_outdir(params$output_dir),
  error = function(e) stop("Cannot create output directory: ", params$output_dir, " (", e$message, ")")
)

# Step 1: Load CSV expression data
cat("\n--- Step 1: Loading CSV expression data ---\n")
source(file.path(env$base_dir, "src", "loading", "load_csv_expression.R"))
expr_data <- load_csv_expression(
  input_dir = params$csv_folder
)
cat("Loaded", length(expr_data), "expression matrices\n")

# Step 2: Marker statistics
cat("\n--- Step 2: Marker Statistics ---\n")
source(file.path(env$base_dir, "src", "analysis", "marker_statistics.R"))
stats_result <- run_marker_statistics(
  input_folder  = params$csv_folder,
  output_folder = file.path(params$output_dir, "01_statistics"),
  markers       = params$markers,
  test_type     = params$test_type
)

# Step 3: Visualization
cat("\n--- Step 3: Visualization ---\n")
source(file.path(env$base_dir, "src", "analysis", "expression_matrices.R"))
source(file.path(env$base_dir, "src", "visualization", "heatmaps.R"))

combined_matrix <- combine_expression_matrices(
  expr_data_list = expr_data,
  markers        = params$markers,
  max_cells      = params$max_cells_per_file %||% 500
)
if (!is.null(combined_matrix)) {
  viz_dir <- file.path(params$output_dir, "02_visualization")
  make_outdir(viz_dir)
  plot_heatmap(
    expr_matrix = t(combined_matrix),
    output_path = file.path(viz_dir, "expression_heatmap.png"),
    title       = "Expression Heatmap (sampled cells)"
  )
}

finish_pipeline("Pipeline B", env$start_time)
