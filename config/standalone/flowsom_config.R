# =============================================================================
# Standalone config: src/clustering/flowsom.R
# =============================================================================
# Runs FlowSOM self-organizing map clustering on expression data (FCS or CSV).
# Builds a grid, then metaclusters into k groups using consensus clustering.
#
# Usage:
#   source("config/standalone/flowsom_config.R")
#   source("src/clustering/flowsom.R")
#   result <- run_flowsom(
#     input_dir           = input_dir,
#     output_dir          = output_dir,
#     markers_of_interest = markers_of_interest,
#     n_clusters          = n_clusters,
#     grid_xdim           = grid_xdim,
#     grid_ydim           = grid_ydim,
#     seed                = seed
#   )
#
# Or from CLI:
#   Rscript src/clustering/flowsom.R <input_dir> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

input_dir           <- "/path/to/expression/data"   # Directory with FCS or expression CSV files
output_dir          <- "/path/to/output/flowsom"
markers_of_interest <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

# Number of metaclusters
n_clusters <- 20

# SOM grid dimensions (total nodes = xdim * ydim)
grid_xdim <- 10
grid_ydim <- 10

seed <- 42
