#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_A_raw_fcs.R — Raw FCS -> Full Analysis Pipeline
# =============================================================================
# Flow: Load FCS -> QC -> Transform -> Homogenize -> Downsample -> Cluster -> Analyze -> Visualize
# =============================================================================

cat("=== Pipeline A: Raw FCS -> Full Analysis ===\n")
start_time <- Sys.time()

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
if (!is.null(config_file)) {
  cat("Loading config:", config_file, "\n")
  source(config_file)
} else {
  cat("No --config provided. Using default parameters.\n")
  source(file.path(dirname(getwd()), "config", "pipeline_A_config.R"))
}

# --- Source shared modules ---
base_dir <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

# Step 1: Load FCS files
cat("\n--- Step 1: Loading FCS files ---\n")
source(file.path(base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  input_dir  = params$fcs_folder,
  truncation = params$truncation
)
cat("Loaded", length(flow_data$flow_frames), "samples\n")

# Step 2: Preprocessing (QC + Transform + Normalize)
cat("\n--- Step 2: Preprocessing ---\n")
source(file.path(base_dir, "src", "preprocessing", "preprocess_pipeline.R"))
preprocessed <- run_preprocessing(
  flow_frames = flow_data$flow_frames,
  params      = params,
  output_dir  = file.path(params$output_dir, "01_preprocessing")
)

# Step 3: Panel homogenization
cat("\n--- Step 3: Panel Homogenization ---\n")
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames         = preprocessed,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  markername_replacing = params$markername_replacing,
  output_dir          = file.path(params$output_dir, "02_harmonization")
)

# Step 4: Downsampling
cat("\n--- Step 4: Downsampling ---\n")
source(file.path(base_dir, "src", "preprocessing", "downsampling.R"))
downsampled <- downsample_flowframes(
  flow_frames  = homogenized,
  sample_names = flow_data$sample_names,
  target_cells = params$downsample_n,
  method       = params$downsample_method,
  seed         = params$seed
)

# Step 5: Build combined expression matrix for clustering
cat("\n--- Step 5: Clustering ---\n")
source(file.path(base_dir, "src", "clustering", "leiden.R"))

# Build combined matrix + sample labels from downsampled flowFrames
expr_list <- extract_exprmtx_list(downsampled, flow_data$sample_names, params$markers_of_interest)
combined_matrix <- do.call(rbind, expr_list)
sample_labels <- rep(names(expr_list), sapply(expr_list, nrow))

clusters <- run_leiden_with_umap(
  expr_matrix   = combined_matrix,
  sample_labels = sample_labels,
  output_dir    = file.path(params$output_dir, "03_clustering"),
  n_neighbors   = params$n_neighbors,
  resolution    = params$resolution,
  seed          = params$seed
)

# Step 6: Save clustered data as CSVs for split_by_cluster (reads from directory)
cat("\n--- Step 6: Statistical Analysis ---\n")
cluster_data_dir <- file.path(params$output_dir, "04_analysis", "cluster_data")
make_outdir(cluster_data_dir)

clustered_df <- as.data.frame(combined_matrix)
clustered_df$sample <- sample_labels
clustered_df$cluster <- clusters$cluster

# Write per-sample CSVs with gating-style binary columns
for (sname in unique(sample_labels)) {
  sdata <- clustered_df[clustered_df$sample == sname, ]
  # Create binary cluster columns in "Clustering annotated" format
  for (cl in unique(clusters$cluster)) {
    col_name <- paste0("Cluster_", cl, " (Clustering annotated)")
    sdata[[col_name]] <- as.integer(sdata$cluster == cl)
  }
  sdata$cluster <- NULL
  sdata$sample <- NULL
  data.table::fwrite(sdata, file.path(cluster_data_dir, paste0(sname, ".csv")), sep = ";")
}

source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
split_data <- split_by_cluster(
  input_dir     = cluster_data_dir,
  output_dir    = file.path(params$output_dir, "04_analysis", "split_clusters"),
  cluster_regex = params$cluster_regex
)

# Run differential abundance
source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))
tryCatch({
  da_results <- run_differential_abundance(
    data_dir         = cluster_data_dir,
    conditions_csv   = params$conditions_csv,
    output_dir       = file.path(params$output_dir, "04_analysis", "differential_abundance"),
    cluster_regex    = params$cluster_regex,
    condition_levels = params$condition_levels,
    min_cluster_size = params$min_cluster_size,
    p_adjust_method  = params$p_adjust_method,
    seed             = params$seed
  )
}, error = function(e) {
  cat("  Differential abundance skipped:", e$message, "\n")
})

# Step 7: Visualization
cat("\n--- Step 7: Visualization ---\n")
source(file.path(base_dir, "src", "visualization", "umap.R"))
plot_umap(
  flow_frames         = downsampled,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  output_dir          = file.path(params$output_dir, "05_visualization"),
  n_neighbors         = params$n_neighbors,
  output_formats      = params$output_formats,
  plot_width          = params$plot_width,
  plot_height         = params$plot_height,
  seed                = params$seed
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline A complete (%.1f minutes) ===\n", elapsed))
