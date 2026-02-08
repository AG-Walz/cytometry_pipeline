# =============================================================================
# config/pipeline_E_config.R — Pipeline E: Preprocessed FCS → UMAP/Clustering
# =============================================================================
# Flow: Load FCS (skip QC) → Homogenize → Downsample → Leiden → UMAP
#
# Use this pipeline when you have already preprocessed FCS files (QC'd,
# transformed, normalized) and want to quickly run unsupervised clustering
# and UMAP visualization. Skips quality control and transformation steps.
#
# Usage:
#   Rscript pipelines/pipeline_E_preprocessed_umap.R --config config/pipeline_E_config.R
# Or in R:
#   source("config/pipeline_E_config.R")
#   source("pipelines/pipeline_E_preprocessed_umap.R")
# =============================================================================

params <- list(

  # ========================= EDIT THESE ======================================

  # -- Input / Output --
  fcs_folder = "/path/to/preprocessed/fcs/files",  # Preprocessed FCS files
  output_dir = "/path/to/output",                   # All results saved here

  # -- Markers --
  markers = c(
    "CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"
  ),

  # ========================= OPTIONAL (sensible defaults) ====================

  # -- Panel Homogenization --
  markername_replacing = list(
    # `PD1 (CD279)` = "PD1",
    # `Fas (CD95)` = "CD95"
  ),

  # -- Downsampling --
  downsample_method = "Consistent",  # "Consistent" (fixed n per sample) or "Maximal" (proportional)
  downsample_n      = 10000,         # Target cells per sample (Consistent) or total (Maximal)

  # -- Clustering (Leiden) --
  n_neighbors = 30,     # KNN graph neighbors (higher = coarser clusters)
  resolution  = 0.8,    # Leiden resolution (higher = more clusters)

  # -- Output --
  output_formats = list(pdf = TRUE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # -- Reproducibility --
  seed = 42
)
