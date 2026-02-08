# =============================================================================
# Standalone config: src/analysis/abundance_plotting.R
# =============================================================================
# Creates barplots from a cluster-by-sample tracking matrix showing cell
# distribution across clusters. Supports stacked, grouped, and proportional.
#
# Usage:
#   source("config/standalone/abundance_plotting_config.R")
#   source("src/analysis/abundance_plotting.R")
#   plot_abundance(
#     tracking_matrix_path = tracking_matrix_path,
#     conditions_csv       = conditions_csv,
#     output_dir           = output_dir,
#     plot_type            = plot_type
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

tracking_matrix_path <- "/path/to/tracking_matrix.csv"  # From cluster_splitting output
output_dir           <- "/path/to/output/plots"

# ========================= OPTIONAL ==========================================

# CSV with condition assignments (columns: sample, condition). NULL = no condition grouping.
conditions_csv <- NULL

# Plot layout: "stacked", "grouped", or "proportional"
plot_type <- "stacked"
