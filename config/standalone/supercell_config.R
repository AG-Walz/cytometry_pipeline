# =============================================================================
# Standalone config: src/preprocessing/supercell.R
# =============================================================================
# Aggregates single cells into supercells (metacells) using the SuperCell
# algorithm. Reduces dataset size while preserving biological heterogeneity,
# enabling faster downstream clustering and visualization.
#
# Usage:
#   source("config/standalone/supercell_config.R")
#   source("src/preprocessing/supercell.R")
#   result <- run_supercell(
#     input_dir           = input_dir,
#     output_dir          = output_dir,
#     markers_of_interest = markers_of_interest,
#     gam                 = gam,
#     k_knn               = k_knn,
#     n_pc                = n_pc,
#     do_pca              = do_pca,
#     seed                = seed
#   )
#
# Or from CLI:
#   Rscript src/preprocessing/supercell.R <input_dir> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

input_dir           <- "/path/to/expression/data"      # Directory with FCS or expression CSV files
output_dir          <- "/path/to/output/supercell"
markers_of_interest <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

# Graining level: average number of cells per supercell (higher = more compression)
gam <- 20

# Number of nearest neighbors for the kNN graph
k_knn <- 30

# Number of principal components for dimensionality reduction
n_pc <- 10

# Whether to run PCA before building the kNN graph
do_pca <- TRUE

seed <- 42
