# =============================================================================
# Standalone config: analysis/cluster_enrichment.R
# =============================================================================
# Tests whether clusters are enriched (over/under-represented) in one condition
# vs. another using Fisher's exact test or mixed-effects logistic regression.
#
# Usage:
#   source("config/standalone/cluster_enrichment_config.R")
#   source("analysis/cluster_enrichment.R")
#   result <- run_cluster_enrichment(
#     data_dir         = data_dir,
#     conditions_csv   = conditions_csv,
#     output_dir       = output_dir,
#     cluster_regex    = cluster_regex,
#     condition_levels = condition_levels,
#     stat_method      = stat_method,
#     min_cluster_frac = min_cluster_frac,
#     p_adjust_method  = p_adjust_method,
#     seed             = seed
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

data_dir       <- "/path/to/cluster/csvs"    # Directory with per-cluster CSV files
conditions_csv <- "/path/to/conditions.csv"   # CSV: columns "sample" and "condition"
output_dir     <- "/path/to/output/enrichment"

condition_levels <- c("Control", "Treatment")

# ========================= OPTIONAL ==========================================

cluster_regex    <- "Clustering annotated"
stat_method      <- "mixed_effects"  # "mixed_effects" or "pooled" (Fisher's exact)
min_cluster_frac <- 0.001            # Minimum fraction of total cells to include cluster
p_adjust_method  <- "BH"
resume           <- FALSE            # Resume from checkpoint if previous run was interrupted
seed             <- 42
