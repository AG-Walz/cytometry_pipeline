#!/usr/bin/env Rscript
# =============================================================================
# Pipeline A: Raw FCS -> Full Analysis
# =============================================================================
# Flow: Load FCS -> QC -> Transform -> Homogenize -> Downsample -> Cluster
#       -> Analyze -> Visualize
# =============================================================================

# --- Init (config + shared modules) ---
base_dir   <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "pipeline_utils.R"))
env <- init_pipeline("Pipeline A: Raw FCS -> Full Analysis",
                     "config/pipeline_A_config.R")
source(file.path(env$shared_dir, "flowframe_utils.R"))
source(file.path(env$base_dir, "src", "loading", "create_flowframes.R"))

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

# Step 1: Load FCS files
cat("\n--- Step 1: Loading FCS files ---\n")
source(file.path(env$base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  input_dir  = params$fcs_folder,
  truncation = params$truncation
)
cat("Loaded", length(flow_data$flow_frames), "samples\n")

# Step 2: Preprocessing
cat("\n--- Step 2: Preprocessing ---\n")
source(file.path(env$base_dir, "src", "preprocessing", "preprocess_pipeline.R"))
preprocessed <- run_preprocessing(
  flow_frames = flow_data$flow_frames,
  params      = params,
  output_dir  = file.path(params$output_dir, "01_preprocessing")
)

# Step 3: Panel homogenization
cat("\n--- Step 3: Panel Homogenization ---\n")
source(file.path(env$base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames         = preprocessed,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  markername_replacing = params$markername_replacing,
  output_dir          = file.path(params$output_dir, "02_harmonization")
)

# Step 3b: Optional batch correction
if ((params$batch_method %||% "none") != "none") {
  cat("\n--- Step 3b: Batch Correction ---\n")
  source(file.path(env$shared_dir, "batch_correction_dispatch.R"))
  batch_result <- run_batch_correction_dispatch(
    data                = homogenized,
    method              = params$batch_method,
    batch_labels        = params$batch_names,
    markers_of_interest = params$markers_of_interest,
    sample_names        = flow_data$sample_names,
    params              = batch_params(params)
  )
  if (!is.null(batch_result$corrected_frames)) {
    homogenized <- batch_result$corrected_frames
  } else if (!is.null(batch_result$corrected_matrix)) {
    # Convert corrected matrix back to flowFrame list for Harmony/ComBat
    message("Converting batch-corrected matrix back to flowFrame list...")
    homogenized <- create_flowframes_from_matrices(
      expr_list = split(as.data.frame(batch_result$corrected_matrix),
                        batch_result$sample_labels),
      markers_of_interest = batch_result$markers
    )
  }
}

# Step 4: Downsampling
cat("\n--- Step 4: Downsampling ---\n")
source(file.path(env$base_dir, "src", "preprocessing", "downsampling.R"))
downsampled <- downsample_flowframes(
  flow_frames  = homogenized,
  sample_names = flow_data$sample_names,
  target_cells = params$downsample_n,
  method       = params$downsample_method,
  seed         = params$seed
)

# Step 5: Clustering
cat("\n--- Step 5: Clustering ---\n")
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

# Step 6: Statistical analysis
cat("\n--- Step 6: Statistical Analysis ---\n")
cluster_data_dir <- file.path(params$output_dir, "04_analysis", "cluster_data")
clusters_to_gating_csv(
  combined_matrix = expr_data$combined_matrix,
  sample_labels   = expr_data$sample_labels,
  cluster_ids     = cluster_ids,
  output_dir      = cluster_data_dir
)
source(file.path(env$base_dir, "src", "analysis", "cluster_splitting.R"))
split_data <- split_by_cluster(
  input_dir     = cluster_data_dir,
  output_dir    = file.path(params$output_dir, "04_analysis", "split_clusters"),
  cluster_regex = params$cluster_regex
)
source(file.path(env$base_dir, "src", "analysis", "differential_abundance.R"))
tryCatch({
  run_differential_abundance(
    data_dir         = cluster_data_dir,
    conditions_csv   = params$conditions_csv,
    output_dir       = file.path(params$output_dir, "04_analysis", "differential_abundance"),
    cluster_regex    = params$cluster_regex,
    condition_levels = params$condition_levels,
    min_cluster_size = params$min_cluster_size,
    p_adjust_method  = params$p_adjust_method,
    seed             = params$seed
  )
}, error = function(e) cat("  Differential abundance skipped:", e$message, "\n"))

# Step 7: Visualization
cat("\n--- Step 7: Visualization ---\n")
source(file.path(env$shared_dir, "dimred_dispatch.R"))
run_dimred(
  expr_matrix    = expr_data$combined_matrix,
  method         = params$dimred_method %||% "umap",
  sample_labels  = expr_data$sample_labels,
  cluster_labels = cluster_ids,
  output_dir     = file.path(params$output_dir, "05_visualization"),
  params         = dimred_params(params)
)

finish_pipeline("Pipeline A", env$start_time)
