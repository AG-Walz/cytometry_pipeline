# =============================================================================
# Standalone config: src/visualization/barplots.R
# =============================================================================
# Creates composition barplots showing cell distribution across clusters
# and samples. Supports stacked, grouped, and proportional layouts.
#
# Usage:
#   source("config/standalone/barplots_config.R")
#   source("src/visualization/barplots.R")
#   plot_composition_barplot(
#     data       = data,        # data.frame with cluster and sample columns
#     output_dir = output_dir,
#     x_col      = x_col,
#     fill_col   = fill_col,
#     position   = position,
#     title      = title
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

output_dir <- "/path/to/output/barplots"

# ========================= OPTIONAL ==========================================

x_col    <- "cluster"             # Column for x-axis categories
fill_col <- "sample"              # Column for fill/color
position <- "stack"               # "stack", "dodge", or "fill" (proportional)
title    <- "Cell Composition"
