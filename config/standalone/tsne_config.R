# =============================================================================
# Standalone config: src/visualization/tsne.R
# =============================================================================
# Computes t-SNE embeddings and generates per-sample, per-marker, and
# per-cluster overlay plots. Accepts FCS files or a pre-built expression matrix.
#
# Usage:
#   source("config/standalone/tsne_config.R")
#   source("src/loading/load_fcs_raw.R")
#   source("src/visualization/tsne.R")
#   flow_data <- load_fcs_raw(input_dir = fcs_folder)
#   result <- plot_tsne(
#     flow_frames         = flow_data$flow_frames,
#     sample_names        = flow_data$sample_names,
#     markers_of_interest = markers_of_interest,
#     output_dir          = output_dir,
#     perplexity          = perplexity,
#     max_iter            = max_iter,
#     output_formats      = output_formats,
#     seed                = seed
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

fcs_folder          <- "/path/to/fcs/files"        # Or use pre-loaded input_matrix
output_dir          <- "/path/to/output/tsne"
markers_of_interest <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

# Pre-built expression matrix (if provided, fcs_folder is ignored)
input_matrix <- NULL

# t-SNE perplexity (higher = more global structure; typical range 5-50)
perplexity <- 30

# Maximum number of iterations for optimization
max_iter <- 1000

output_formats <- list(pdf = FALSE, png = TRUE, svg = FALSE)
seed           <- 42
