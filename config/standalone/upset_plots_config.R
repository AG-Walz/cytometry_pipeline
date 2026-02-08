# =============================================================================
# Standalone config: src/visualization/upset_plots.R
# =============================================================================
# Creates UpSet plots for visualizing set intersections (e.g., shared
# populations across gating strategies or shared markers across panels).
#
# Usage (from data.frame with binary set columns):
#   source("config/standalone/upset_plots_config.R")
#   source("src/visualization/upset_plots.R")
#   plot_upset(
#     data        = data,           # data.frame with binary 0/1 set columns
#     set_columns = set_columns,    # names of binary columns
#     output_path = output_path,
#     title       = title,
#     min_size    = min_size
#   )
#
# Usage (from named list of character vectors):
#   plot_upset_from_list(
#     sets        = sets,           # named list of character vectors
#     output_path = output_path,
#     title       = title
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

output_path <- "/path/to/output/upset_plot.png"

# ========================= OPTIONAL ==========================================

title    <- "UpSet Plot"
min_size <- 1    # Minimum intersection size to show
