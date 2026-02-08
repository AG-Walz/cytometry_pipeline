# =============================================================================
# Standalone config: analysis/expression_matrices.R
# =============================================================================
# Generates grouped expression matrices from per-cluster CSV data. Useful for
# creating input for heatmaps or downstream analysis tools.
#
# Usage:
#   source("config/standalone/expression_matrices_config.R")
#   source("analysis/expression_matrices.R")
#   result <- create_expression_matrices(
#     input_dir        = input_dir,
#     output_dir       = output_dir,
#     markers          = markers,
#     group_by_pattern = group_by_pattern
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

input_dir  <- "/path/to/cluster/csvs"         # Per-cluster CSV files
output_dir <- "/path/to/output/expr_matrices"  # Grouped matrices saved here

# ========================= OPTIONAL ==========================================

# Markers to include (NULL = all numeric columns)
markers <- NULL

# Regex pattern for grouping files (NULL = no grouping)
group_by_pattern <- NULL
