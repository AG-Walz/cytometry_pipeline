# =============================================================================
# Standalone config: src/visualization/venn_diagrams.R
# =============================================================================
# Creates Venn diagrams (non-proportional) and area-proportional Euler
# diagrams from named lists of sets. Commonly used with Pipeline D output
# to compare gating strategies.
#
# Usage:
#   source("config/standalone/venn_diagrams_config.R")
#   source("src/visualization/venn_diagrams.R")
#
#   # Area-proportional Euler diagram
#   plot_euler(
#     sets       = sets,           # named list of character vectors
#     output_dir = output_dir,
#     filename   = filename
#   )
#
#   # Traditional Venn diagram
#   plot_venn(
#     sets       = sets,
#     filename   = filename,
#     output_dir = output_dir
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

output_dir <- "/path/to/output/venn"
filename   <- "gating_comparison"

# ========================= OPTIONAL ==========================================

custom_name <- NULL    # Custom title for the plot
scale_size  <- FALSE   # Scale circle areas proportionally
