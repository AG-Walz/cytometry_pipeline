#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_E_preprocessed_umap.R — Preprocessed -> UMAP/Clustering
# =============================================================================
# Flow: Load FCS (no QC) -> Homogenize -> Downsample -> Cluster -> UMAP
# =============================================================================

cat("=== Pipeline E: Preprocessed FCS -> UMAP/Clustering ===\n")
start_time <- Sys.time()

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
if (!is.null(config_file)) {
  source(config_file)
} else {
  source(file.path(dirname(getwd()), "config", "pipeline_E_config.R"))
}

base_dir <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

# Step 1: Load preprocessed FCS (skip QC)
cat("\n--- Step 1: Loading preprocessed FCS files ---\n")
source(file.path(base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  input_dir  = params$fcs_folder,
  truncation = FALSE
)
cat("Loaded", length(flow_data$flow_frames), "samples\n")

# Step 2: Homogenize panels
cat("\n--- Step 2: Panel Homogenization ---\n")
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames         = flow_data$flow_frames,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  markername_replacing = params$markername_replacing,
  output_dir          = file.path(params$output_dir, "01_harmonized")
)

# Step 3: Downsample
cat("\n--- Step 3: Downsampling ---\n")
source(file.path(base_dir, "src", "preprocessing", "downsampling.R"))
downsampled <- downsample_flowframes(
  flow_frames  = homogenized,
  sample_names = flow_data$sample_names,
  target_cells = params$downsample_n,
  method       = params$downsample_method,
  seed         = params$seed
)

# Step 4: Build combined expression matrix + Clustering + UMAP
cat("\n--- Step 4: Leiden Clustering + UMAP ---\n")
source(file.path(base_dir, "src", "clustering", "leiden.R"))

expr_list <- extract_exprmtx_list(downsampled, flow_data$sample_names, params$markers_of_interest)
combined_matrix <- do.call(rbind, expr_list)
sample_labels <- rep(names(expr_list), sapply(expr_list, nrow))

clusters <- run_leiden_with_umap(
  expr_matrix   = combined_matrix,
  sample_labels = sample_labels,
  output_dir    = file.path(params$output_dir, "02_clustering"),
  n_neighbors   = params$n_neighbors,
  resolution    = params$resolution,
  seed          = params$seed
)

# Step 5: UMAP visualization
cat("\n--- Step 5: UMAP Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "umap.R"))
plot_umap(
  flow_frames         = downsampled,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  output_dir          = file.path(params$output_dir, "03_visualization"),
  n_neighbors         = params$n_neighbors,
  output_formats      = params$output_formats,
  plot_width          = params$plot_width,
  plot_height         = params$plot_height,
  seed                = params$seed
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline E complete (%.1f minutes) ===\n", elapsed))
