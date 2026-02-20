#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_C_csv_gating.R — CSV Gating -> Cluster Analysis Pipeline
# =============================================================================
# Flow: Load gating CSV -> Split clusters -> Enrichment -> Barplots/Abundance
# =============================================================================

cat("=== Pipeline C: CSV Gating -> Cluster Analysis ===\n")
start_time <- Sys.time()

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
if (!is.null(config_file)) {
  source(config_file)
} else {
  source(file.path(dirname(getwd()), "config", "pipeline_C_config.R"))
}

base_dir <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))

# Step 1: Load CSV with gating data (also splits per population)
cat("\n--- Step 1: Loading CSV gating data ---\n")
source(file.path(base_dir, "src", "loading", "load_csv_gating.R"))
gating_data <- load_csv_gating(
  input_fcs_path  = params$csv_folder,
  input_csv_path  = params$csv_folder,
  output_path     = file.path(params$output_dir, "01_loaded"),
  gating_pattern  = params$gating_pattern,
  processing_mode = params$processing_mode
)

# Step 2: Cluster splitting from *original* input CSVs (which have gating columns)
cat("\n--- Step 2: Splitting by cluster ---\n")
source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
split_result <- split_by_cluster(
  input_dir     = params$csv_folder,
  output_dir    = file.path(params$output_dir, "02_clusters"),
  cluster_regex = params$cluster_regex
)

# Step 3: Statistical analysis (enrichment or GLMM)
cat("\n--- Step 3: Statistical Analysis ---\n")
if (params$analysis_method == "enrichment") {
  source(file.path(base_dir, "src", "analysis", "cluster_enrichment.R"))
  tryCatch({
    enrichment <- run_cluster_enrichment(
      data_dir         = params$csv_folder,
      conditions_csv   = params$conditions_csv,
      output_dir       = file.path(params$output_dir, "03_analysis"),
      cluster_regex    = params$cluster_regex,
      condition_levels = params$condition_levels,
      stat_method      = params$stat_method,
      p_adjust_method  = params$p_adjust_method,
      seed             = params$seed
    )
  }, error = function(e) {
    cat("  Enrichment analysis skipped:", e$message, "\n")
  })
} else {
  source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))
  tryCatch({
    abundance <- run_differential_abundance(
      data_dir         = params$csv_folder,
      conditions_csv   = params$conditions_csv,
      output_dir       = file.path(params$output_dir, "03_analysis"),
      cluster_regex    = params$cluster_regex,
      condition_levels = params$condition_levels,
      min_cluster_size = params$min_cluster_size,
      p_adjust_method  = params$p_adjust_method,
      seed             = params$seed
    )
  }, error = function(e) {
    cat("  Differential abundance skipped:", e$message, "\n")
  })
}

# Step 4: Visualization — composition barplot from tracking matrix
cat("\n--- Step 4: Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "barplots.R"))
if (!is.null(split_result$tracking_matrix) && nrow(split_result$tracking_matrix) > 0) {
  long_data <- tracking_matrix_to_long(split_result$tracking_matrix)
  plot_composition_barplot(
    data       = long_data,
    output_dir = file.path(params$output_dir, "04_visualization"),
    x_col      = "cluster",
    fill_col   = "file",
    title      = "Cell Composition by Cluster"
  )
}

# Step 5: Abundance analysis
cat("\n--- Step 5: Abundance Analysis ---\n")
source(file.path(base_dir, "src", "analysis", "abundance_plotting.R"))
tryCatch({
  tracking_matrix <- build_tracking_matrix(
    input_dir        = params$csv_folder,
    cluster_regex    = params$cluster_regex,
    exclude_patterns = params$exclude_patterns,
    output_dir       = file.path(params$output_dir, "05_abundance")
  )
  plot_abundance(
    output_dir      = file.path(params$output_dir, "05_abundance"),
    tracking_matrix = tracking_matrix,
    metadata_file   = params$cytolution_metadata_file
  )
}, error = function(e) {
  cat("  Abundance plotting skipped:", e$message, "\n")
})

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline C complete (%.1f minutes) ===\n", elapsed))
