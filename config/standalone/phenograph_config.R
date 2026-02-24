# =============================================================================
# Standalone config: src/clustering/phenograph.R
# =============================================================================
# Runs PhenoGraph clustering on expression data. Builds a k-nearest-neighbor
# graph and applies Louvain community detection to identify cell populations.
#
# Usage:
#   source("config/standalone/phenograph_config.R")
#   source("src/clustering/phenograph.R")
#   result <- run_phenograph(
#     input_dir           = input_dir,
#     output_dir          = output_dir,
#     markers_of_interest = markers_of_interest,
#     k                   = k,
#     seed                = seed
#   )
#
# Or from CLI:
#   Rscript src/clustering/phenograph.R <input_dir> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

input_dir           <- "/path/to/expression/data"   # Directory with FCS or expression CSV files
output_dir          <- "/path/to/output/phenograph"
markers_of_interest <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

# Number of nearest neighbors for the kNN graph
k <- 30

seed <- 42
