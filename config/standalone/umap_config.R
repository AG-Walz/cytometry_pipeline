# =============================================================================
# Standalone config: src/visualization/umap.R
# =============================================================================
# Computes UMAP embeddings and generates per-sample, per-marker, and
# per-cluster overlay plots.
#
# Usage:
#   source("config/standalone/umap_config.R")
#   source("src/loading/load_fcs_raw.R")
#   source("src/visualization/umap.R")
#   flow_data <- load_fcs_raw(input_dir = fcs_folder)
#   result <- plot_umap(
#     flow_frames           = flow_data$flow_frames,
#     sample_names          = flow_data$sample_names,
#     markers_of_interest   = markers_of_interest,
#     output_dir            = output_dir,
#     n_neighbors           = n_neighbors,
#     min_dist              = min_dist,
#     marker_selection_mode = marker_selection_mode,
#     output_formats        = output_formats,
#     seed                  = seed
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

fcs_folder          <- "/path/to/fcs/files"      # Or use pre-loaded flow_frames
output_dir          <- "/path/to/output/umap"
markers_of_interest <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

n_neighbors           <- 15       # UMAP neighbors (higher = more global structure)
min_dist              <- 0.1      # UMAP min distance (lower = tighter clusters)
marker_selection_mode <- "all"    # "all", "gui", or "manual"
output_formats        <- list(pdf = FALSE, png = TRUE, svg = FALSE)
plot_width            <- 8
plot_height           <- 6
seed                  <- 42
