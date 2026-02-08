#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_C_csv_gating.R — CSV Gating → Cluster Analysis Pipeline
# =============================================================================
# Flow: Load gating CSV → Split clusters → GLMM/Enrichment → UMAP/Barplots
# =============================================================================

cat("=== Pipeline C: CSV Gating → Cluster Analysis ===\n")
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

# Step 1: Load CSV with gating data
cat("\n--- Step 1: Loading CSV gating data ---\n")
source(file.path(base_dir, "src", "loading", "load_csv_gating.R"))
gating_data <- load_csv_gating(
  csv_folder = params$csv_folder,
  output_dir = file.path(params$output_dir, "01_loaded")
)

# Step 2: Cluster splitting
cat("\n--- Step 2: Splitting by cluster ---\n")
source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
split_result <- split_by_cluster(
  data = gating_data$combined,
  cluster_labels = gating_data$cluster_labels,
  output_dir = file.path(params$output_dir, "02_clusters")
)

# Step 3: Statistical analysis (GLMM or enrichment)
cat("\n--- Step 3: Statistical Analysis ---\n")
source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))
source(file.path(base_dir, "src", "analysis", "cluster_enrichment.R"))

if (exists("params") && !is.null(params$analysis_method) && params$analysis_method == "enrichment") {
  enrichment <- run_cluster_enrichment(
    tracking_matrix = split_result$tracking_matrix,
    condition_col = params$condition_col,
    output_dir = file.path(params$output_dir, "03_analysis")
  )
} else {
  abundance <- run_differential_abundance(
    tracking_matrix = split_result$tracking_matrix,
    condition_col = params$condition_col,
    sample_col = params$sample_col,
    output_dir = file.path(params$output_dir, "03_analysis")
  )
}

# Step 4: Visualization
cat("\n--- Step 4: Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "umap.R"))
source(file.path(base_dir, "src", "visualization", "barplots.R"))
source(file.path(base_dir, "src", "visualization", "heatmaps.R"))

# Step 5: Abundance analysis
cat("\n--- Step 5: Abundance Analysis ---\n")
source(file.path(base_dir, "src", "analysis", "abundance_plotting.R"))
tracking_matrix <- build_tracking_matrix(
  input_dir        = params$csv_folder,
  cluster_regex    = params$cluster_regex,
  exclude_patterns = params$exclude_patterns,
  output_dir       = file.path(params$output_dir, "05_abundance")
)
plot_abundance(
  output_dir       = file.path(params$output_dir, "05_abundance"),
  tracking_matrix  = tracking_matrix,
  metadata_file    = params$cytolution_metadata_file
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline C complete (%.1f minutes) ===\n", elapsed))
