#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_A_raw_fcs.R — Raw FCS → Full Analysis Pipeline
# =============================================================================
# Flow: Load FCS → QC → Transform → Homogenize → Downsample → Cluster → Analyze → Visualize
# =============================================================================

cat("=== Pipeline A: Raw FCS → Full Analysis ===\n")
start_time <- Sys.time()

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
if (!is.null(config_file)) {
  cat("Loading config:", config_file, "\n")
  source(config_file)
} else {
  cat("No --config provided. Using default parameters.\n")
  source(file.path(dirname(getwd()), "config", "Parameters_template.R"))
}

# --- Source modules ---
base_dir <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))

# Step 1: Load FCS files
cat("\n--- Step 1: Loading FCS files ---\n")
source(file.path(base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  fcs_folder = params$fcs_folder,
  truncate_max_range = params$truncate_max_range
)
cat("Loaded", length(flow_data$flow_frames), "samples\n")

# Step 2: Preprocessing (QC + Transform + Normalize)
cat("\n--- Step 2: Preprocessing ---\n")
source(file.path(base_dir, "src", "preprocessing", "preprocess_pipeline.R"))
preprocessed <- run_preprocessing(
  flow_frames = flow_data$flow_frames,
  sample_names = flow_data$sample_names,
  markers = params$markers,
  qc_method = params$qc_method,
  transform_method = params$transform_method,
  cofactor = params$cofactor,
  output_dir = file.path(params$output_dir, "01_preprocessing")
)

# Step 3: Panel homogenization
cat("\n--- Step 3: Panel Homogenization ---\n")
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames = preprocessed$flow_frames,
  output_dir = file.path(params$output_dir, "02_harmonization")
)

# Step 4: Downsampling
cat("\n--- Step 4: Downsampling ---\n")
source(file.path(base_dir, "src", "preprocessing", "downsampling.R"))
downsampled <- downsample_flowframes(
  flow_frames = homogenized$flow_frames,
  sample_names = flow_data$sample_names,
  method = params$downsample_method,
  n_cells = params$downsample_n,
  output_dir = file.path(params$output_dir, "03_downsampled")
)

# Step 5: Clustering
cat("\n--- Step 5: Clustering ---\n")
source(file.path(base_dir, "src", "clustering", "leiden.R"))
clusters <- run_leiden_with_umap(
  expr_matrix = downsampled$combined_matrix,
  markers = params$markers,
  n_neighbors = params$n_neighbors,
  resolution = params$resolution,
  output_dir = file.path(params$output_dir, "04_clustering")
)

# Step 6: Analysis
cat("\n--- Step 6: Statistical Analysis ---\n")
source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))
split_data <- split_by_cluster(
  data = downsampled$combined_data,
  cluster_labels = clusters$labels,
  output_dir = file.path(params$output_dir, "05_analysis")
)

# Step 7: Visualization
cat("\n--- Step 7: Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "umap.R"))
source(file.path(base_dir, "src", "visualization", "heatmaps.R"))
plot_umap(
  umap_coords = clusters$umap_coords,
  cluster_labels = clusters$labels,
  sample_labels = downsampled$sample_labels,
  expr_matrix = downsampled$combined_matrix,
  markers = params$markers,
  output_dir = file.path(params$output_dir, "06_visualization")
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline A complete (%.1f minutes) ===\n", elapsed))
