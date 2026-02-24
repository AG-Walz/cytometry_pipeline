#!/usr/bin/env Rscript
# =============================================================================
# Pipeline E: Preprocessed FCS -> UMAP/Clustering
# =============================================================================
# Flow: Load FCS (no QC) -> Homogenize -> Downsample -> Cluster -> UMAP
# =============================================================================

# --- Init (config + shared modules) ---
base_dir   <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "pipeline_utils.R"))
env <- init_pipeline("Pipeline E: Preprocessed FCS -> UMAP/Clustering",
                     "config/pipeline_E_config.R")
source(file.path(env$shared_dir, "flowframe_utils.R"))

# Step 1: Load preprocessed FCS (skip QC)
cat("\n--- Step 1: Loading preprocessed FCS files ---\n")
source(file.path(env$base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  input_dir  = params$fcs_folder,
  truncation = FALSE
)
cat("Loaded", length(flow_data$flow_frames), "samples\n")

# Step 2: Homogenize panels
cat("\n--- Step 2: Panel Homogenization ---\n")
source(file.path(env$base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames         = flow_data$flow_frames,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  markername_replacing = params$markername_replacing,
  output_dir          = file.path(params$output_dir, "01_harmonized")
)

# Step 3: Downsample
cat("\n--- Step 3: Downsampling ---\n")
source(file.path(env$base_dir, "src", "preprocessing", "downsampling.R"))
downsampled <- downsample_flowframes(
  flow_frames  = homogenized,
  sample_names = flow_data$sample_names,
  target_cells = params$downsample_n,
  method       = params$downsample_method,
  seed         = params$seed
)

# Step 4: Clustering
cat("\n--- Step 4: Clustering ---\n")
source(file.path(env$base_dir, "src", "clustering", "leiden.R"))
source(file.path(env$shared_dir, "clustering_dispatch.R"))
expr_data <- build_combined_expr_matrix(
  downsampled, flow_data$sample_names, params$markers_of_interest
)
cluster_ids <- run_clustering(
  expr_matrix = expr_data$combined_matrix,
  method      = params$clustering_method %||% "leiden",
  params      = clustering_params(params)
)

# Step 5: Visualization
cat("\n--- Step 5: Visualization ---\n")
source(file.path(env$shared_dir, "dimred_dispatch.R"))
run_dimred(
  expr_matrix    = expr_data$combined_matrix,
  method         = params$dimred_method %||% "umap",
  sample_labels  = expr_data$sample_labels,
  cluster_labels = cluster_ids,
  output_dir     = file.path(params$output_dir, "03_visualization"),
  params         = dimred_params(params)
)

finish_pipeline("Pipeline E", env$start_time)
