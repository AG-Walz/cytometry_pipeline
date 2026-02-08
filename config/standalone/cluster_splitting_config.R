# =============================================================================
# Standalone config: analysis/cluster_splitting.R
# =============================================================================
# Splits combined expression CSVs into per-cluster files and generates a
# sample-by-cluster tracking matrix (cell counts).
#
# Usage:
#   source("config/standalone/cluster_splitting_config.R")
#   source("analysis/cluster_splitting.R")
#   result <- split_by_cluster(
#     input_dir             = input_dir,
#     output_dir            = output_dir,
#     cluster_regex         = cluster_regex,
#     markers_of_interest   = markers_of_interest,
#     min_cells_per_cluster = min_cells_per_cluster,
#     chunk_size            = chunk_size
#   )
#
# Or from CLI:
#   Rscript analysis/cluster_splitting.R <input_dir> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

input_dir  <- "/path/to/expression/csvs"   # Directory with combined expression CSVs
output_dir <- "/path/to/output/clusters"    # Per-cluster CSVs saved here

# ========================= OPTIONAL ==========================================

# Regex to identify cluster/gating columns in your CSVs
cluster_regex <- "Clustering annotated"

# Markers to keep in the per-cluster files (NULL = keep all numeric columns)
markers_of_interest <- NULL

# Drop clusters with fewer cells than this
min_cells_per_cluster <- 10

# Rows processed at a time (lower = less memory)
chunk_size <- 10000
