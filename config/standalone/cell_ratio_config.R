# =============================================================================
# Standalone config: src/analysis/cell_ratio.R
# =============================================================================
# Calculates ratios between specified cluster populations per sample
# (e.g., CD4/CD8 ratio). Outputs a table with per-sample ratios and log2 ratios.
#
# Usage:
#   source("config/standalone/cell_ratio_config.R")
#   source("src/analysis/cell_ratio.R")
#   result <- calculate_cell_ratios(
#     abundance_csv        = abundance_csv,
#     numerator_clusters   = numerator_clusters,
#     denominator_clusters = denominator_clusters,
#     output_path          = output_path
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

abundance_csv <- "/path/to/tracking_matrix.csv"  # From cluster_splitting output

# Cluster names (as they appear in the tracking matrix) for numerator and denominator.
# Multiple clusters per side are summed before computing the ratio.
numerator_clusters   <- c("CD4+ T cells")
denominator_clusters <- c("CD8+ T cells")

# ========================= OPTIONAL ==========================================

# Output path for the ratio table (NULL = return without saving)
output_path <- "/path/to/output/cd4_cd8_ratios.csv"
