# =============================================================================
# Standalone config: src/analysis/differential_abundance.R
# =============================================================================
# Tests differential abundance of cell clusters between two conditions using
# negative binomial GLMM with sample-level random effects.
#
# Usage:
#   source("config/standalone/differential_abundance_config.R")
#   source("src/analysis/differential_abundance.R")
#   result <- run_differential_abundance(
#     data_dir            = data_dir,
#     conditions_csv      = conditions_csv,
#     output_dir          = output_dir,
#     cluster_regex       = cluster_regex,
#     condition_levels    = condition_levels,
#     markers_of_interest = markers_of_interest,
#     min_cluster_size    = min_cluster_size,
#     p_adjust_method     = p_adjust_method,
#     seed                = seed
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

data_dir       <- "/path/to/cluster/csvs"    # Directory with per-cluster CSV files
conditions_csv <- "/path/to/conditions.csv"   # CSV: columns "sample" and "condition"
output_dir     <- "/path/to/output/da"        # Results saved here

# Two-level condition for comparison (reference level first)
condition_levels <- c("Control", "Treatment")

# ========================= OPTIONAL ==========================================

cluster_regex       <- "Clustering annotated"
markers_of_interest <- NULL    # NULL = use all markers
min_cluster_size    <- 10      # Drop small clusters
p_adjust_method     <- "BH"    # "BH", "bonferroni", "holm"
max_cells_per_file  <- 100000  # Cap per file for memory
seed                <- 42
