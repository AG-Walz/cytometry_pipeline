# =============================================================================
# Standalone config: src/visualization/density_plots.R
# =============================================================================
# Creates per-marker density or ridge plots, grouped by sample or condition.
#
# Usage:
#   source("config/standalone/density_plots_config.R")
#   source("src/loading/load_fcs_raw.R")
#   source("src/visualization/density_plots.R")
#   flow_data <- load_fcs_raw(input_dir = fcs_folder)
#   plot_density(
#     flow_frames  = flow_data$flow_frames,
#     sample_names = flow_data$sample_names,
#     markers      = markers,
#     output_dir   = output_dir,
#     plot_type    = plot_type,
#     group_by     = group_by
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

fcs_folder <- "/path/to/fcs/files"
output_dir <- "/path/to/output/density"
markers    <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

plot_type    <- "density"   # "density" or "ridge"
group_by     <- "sample"    # "sample" or "condition"
group_labels <- NULL        # Named vector mapping sample names to group labels
