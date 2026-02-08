#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_E_preprocessed_umap.R — Preprocessed → UMAP/Clustering
# =============================================================================
# Flow: Load FCS (no QC) → Homogenize → Downsample → Cluster → UMAP
# =============================================================================

cat("=== Pipeline E: Preprocessed FCS → UMAP/Clustering ===\n")
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

# Step 1: Load preprocessed FCS (skip QC)
cat("\n--- Step 1: Loading preprocessed FCS files ---\n")
source(file.path(base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  fcs_folder = params$fcs_folder,
  truncate_max_range = FALSE
)
cat("Loaded", length(flow_data$flow_frames), "samples\n")

# Step 2: Homogenize panels
cat("\n--- Step 2: Panel Homogenization ---\n")
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames = flow_data$flow_frames,
  output_dir = file.path(params$output_dir, "01_harmonized")
)

# Step 3: Downsample
cat("\n--- Step 3: Downsampling ---\n")
source(file.path(base_dir, "src", "preprocessing", "downsampling.R"))
downsampled <- downsample_flowframes(
  flow_frames = homogenized$flow_frames,
  sample_names = flow_data$sample_names,
  method = params$downsample_method,
  n_cells = params$downsample_n,
  output_dir = file.path(params$output_dir, "02_downsampled")
)

# Step 4: Clustering + UMAP
cat("\n--- Step 4: Leiden Clustering + UMAP ---\n")
source(file.path(base_dir, "src", "clustering", "leiden.R"))
clusters <- run_leiden_with_umap(
  expr_matrix = downsampled$combined_matrix,
  markers = params$markers,
  n_neighbors = params$n_neighbors,
  resolution = params$resolution,
  output_dir = file.path(params$output_dir, "03_clustering")
)

# Step 5: UMAP visualization
cat("\n--- Step 5: UMAP Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "umap.R"))
plot_umap(
  umap_coords = clusters$umap_coords,
  cluster_labels = clusters$labels,
  sample_labels = downsampled$sample_labels,
  expr_matrix = downsampled$combined_matrix,
  markers = params$markers,
  output_dir = file.path(params$output_dir, "04_visualization")
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline E complete (%.1f minutes) ===\n", elapsed))
