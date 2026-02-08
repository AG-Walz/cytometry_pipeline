# =============================================================================
# Standalone config: src/visualization/heatmaps.R
# =============================================================================
# Generates clustered heatmaps from expression matrices using pheatmap.
# Two modes: raw expression heatmap, or median-expression per cluster.
#
# Usage (expression heatmap):
#   source("config/standalone/heatmap_config.R")
#   source("src/visualization/heatmaps.R")
#   plot_heatmap(
#     expr_matrix  = expr_matrix,   # numeric matrix (cells x markers)
#     output_path  = output_path,
#     title        = title,
#     scale        = scale
#   )
#
# Usage (cluster median heatmap):
#   plot_cluster_heatmap(
#     expr_data    = expr_data,      # data.frame with cluster + marker columns
#     cluster_col  = "cluster",
#     marker_cols  = markers,
#     output_path  = output_path
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

output_path <- "/path/to/output/heatmap.png"

# For plot_cluster_heatmap:
markers <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

title         <- "Expression Heatmap"
cluster_rows  <- TRUE
cluster_cols  <- TRUE
scale         <- "row"    # "row", "column", or "none"
color_palette <- NULL     # NULL = viridis default
width         <- 10
height        <- 8
