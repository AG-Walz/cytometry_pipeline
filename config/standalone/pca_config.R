# =============================================================================
# Standalone config: src/visualization/pca.R
# =============================================================================
# Performs PCA and generates scatter plots of principal components plus
# loadings plots showing which markers drive each component.
#
# Usage:
#   source("config/standalone/pca_config.R")
#   source("src/visualization/pca.R")
#   result <- plot_pca(
#     expr_matrix   = expr_matrix,     # numeric matrix (cells x markers)
#     sample_labels = sample_labels,   # character vector of sample IDs
#     output_dir    = output_dir,
#     n_components  = n_components,
#     color_by      = color_by,
#     title         = title
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

output_dir <- "/path/to/output/pca"

# ========================= OPTIONAL ==========================================

n_components <- 2          # Number of PCs to plot
color_by     <- "sample"   # "sample", "cluster", or a column name
title        <- "PCA"
