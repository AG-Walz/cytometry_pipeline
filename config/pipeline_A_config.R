# =============================================================================
# config/pipeline_A_config.R — Pipeline A: Raw FCS → Full Analysis
# =============================================================================
# Flow: Load FCS → QC → Transform → Homogenize → Downsample → Cluster → Analyze → Visualize
#
# Usage:
#   Rscript pipelines/pipeline_A_raw_fcs.R --config config/pipeline_A_config.R
# Or in R:
#   source("config/pipeline_A_config.R")
#   source("pipelines/pipeline_A_raw_fcs.R")
# =============================================================================

params <- list(

  # ========================= EDIT THESE ======================================

  # -- Input / Output --
  fcs_folder  = "/path/to/your/fcs/files",       # Directory containing .fcs files
  output_dir  = "/path/to/output",                # All results saved here

  # -- Markers --
  # List the markers (channel names) to include in the analysis.
  # Must match the marker names in your FCS files. Use markername_replacing
  # below if your FCS files use different names than you want in the output.
  markers = c(
    "CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"
  ),

  # -- Conditions --
  # Two-level condition for differential abundance testing (e.g., "Control" vs "Treatment").
  conditions_csv  = "/path/to/conditions.csv",     # CSV with columns: sample, condition
  condition_levels = c("Control", "Treatment"),

  # ========================= OPTIONAL (sensible defaults) ====================

  # -- Loading --
  truncate_max_range = FALSE,   # Truncate extreme values when loading FCS

  # -- Quality Control --
  qc_method        = "PeacoQC",  # "PeacoQC", "flowAI", or "None"

  # -- Transformation --
  transform_method = "arcsinh",  # "arcsinh", "logicle", or "biexponential"
  cofactor         = 5,          # Cofactor for arcsinh (5 for mass, 150 for fluorescence)

  # -- Panel Homogenization --
  # Rename markers across files for consistency (FCS name -> display name)
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

  # -- Analysis --
  cluster_regex    = "Clustering annotated",  # Regex to identify cluster columns
  min_cluster_size = 10,                      # Drop clusters with fewer cells
  p_adjust_method  = "BH",                    # p-value correction: "BH", "bonferroni", "holm"
  stat_method      = "mixed_effects",         # "mixed_effects" or "pooled"

  # -- Output --
  output_formats = list(pdf = TRUE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # -- Reproducibility --
  seed = 42
)
