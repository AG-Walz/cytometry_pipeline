# =============================================================================
# Standalone config: metadata_qc/overlap_statistics.R
# =============================================================================
# Computes Fleiss' kappa for inter-rater agreement across multiple gating
# strategies. Measures how well different analysts agree on population
# assignments for the same cells.
#
# Usage:
#   source("config/standalone/overlap_statistics_config.R")
#   source("metadata_qc/overlap_statistics.R")
#   result <- compute_fleiss_kappa(
#     label_lists   = label_lists,     # Named list of label vectors (from compare_gatings)
#     output_dir    = output_dir,
#     subset_specs  = subset_specs,
#     subset_labels = subset_labels,
#     row_range     = row_range
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

# label_lists is typically the output of compare_gatings()
# Each element is a character vector of population labels for each cell.
output_dir <- "/path/to/output/overlap_stats"

# ========================= OPTIONAL ==========================================

# Subset analysis: compute kappa for specific populations only
subset_specs  <- NULL   # List of vectors specifying population subsets
subset_labels <- NULL   # Labels for each subset
row_range     <- NULL   # c(start, end) to restrict to specific cell rows
