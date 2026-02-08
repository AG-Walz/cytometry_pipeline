# =============================================================================
# Standalone config: src/analysis/marker_statistics.R
# =============================================================================
# Computes per-marker summary statistics (mean, SD, variance components, ICC)
# within each cluster using linear mixed-effects models.
#
# Usage:
#   source("config/standalone/marker_statistics_config.R")
#   source("src/analysis/marker_statistics.R")
#   result <- run_marker_statistics(
#     input_folder   = input_folder,
#     output_folder  = output_folder,
#     filename_regex = filename_regex,
#     markers        = markers,
#     test_type      = test_type
#   )
#
# Or from CLI:
#   Rscript src/analysis/marker_statistics.R <input_folder> <output_folder>
# =============================================================================

# ========================= EDIT THESE ========================================

input_folder  <- "/path/to/cluster/csvs"     # Per-cluster CSV files
output_folder <- "/path/to/output/stats"      # Results saved here

# ========================= OPTIONAL ==========================================

# Regex to extract (population, sample) from filenames
# Default expects files like "ClusterName_SampleID.csv"
filename_regex <- "^(.+)_(.+)\\.csv$"

# Markers to analyze (NULL = use all numeric columns)
markers <- NULL

# Statistical model: "lmer" (linear mixed-effects) or "basic" (simple summary)
test_type <- "lmer"
