#!/usr/bin/env Rscript
# =============================================================================
# ICK Cluster Abundance Analysis
# =============================================================================
# Loads preprocessed CSVs with binary clustering columns, builds a sample-by-
# cluster tracking matrix (parent-level populations only), and generates
# abundance barplots.
#
# Usage:
#   cd pipelines
#   Rscript run_ICK_abundance.R
# =============================================================================

cat("=== ICK Cluster Abundance Analysis ===\n")
start_time <- Sys.time()

# --- Load config ---
base_dir <- dirname(getwd())
source(file.path(base_dir, "config", "study_params", "ICK_abundance.R"))

# --- Source shared utilities ---
source(file.path(base_dir, "src", "analysis", "abundance_plotting.R"))

# --- Step 1: Build tracking matrix ---
cat("\n--- Step 1: Building tracking matrix ---\n")
tracking_matrix <- build_tracking_matrix(
  input_dir        = input_dir,
  cluster_regex    = cluster_regex,
  exclude_patterns = exclude_patterns,
  output_dir       = output_dir
)

# --- Step 2: Abundance plots ---
cat("\n--- Step 2: Generating abundance plots ---\n")
plot_abundance(
  output_dir      = output_dir,
  tracking_matrix = tracking_matrix,
  metadata_file   = cytolution_metadata_file
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== ICK analysis complete (%.1f minutes) ===\n", elapsed))
cat("Output directory:", output_dir, "\n")
