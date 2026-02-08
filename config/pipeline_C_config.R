# =============================================================================
# config/pipeline_C_config.R — Pipeline C: CSV Gating → Cluster Analysis
# =============================================================================
# Flow: Load gating CSV → Split clusters → GLMM/Enrichment → UMAP/Barplots
#
# Use this pipeline when you have CSV files with binary gating columns
# (e.g., exported from Cytolution or manual gating) alongside expression data.
# The pipeline splits data by cluster, runs differential abundance or
# enrichment analysis, and produces UMAP and barplot visualizations.
#
# Usage:
#   Rscript pipelines/pipeline_C_csv_gating.R --config config/pipeline_C_config.R
# Or in R:
#   source("config/pipeline_C_config.R")
#   source("pipelines/pipeline_C_csv_gating.R")
# =============================================================================

params <- list(

  # ========================= EDIT THESE ======================================

  # -- Input / Output --
  csv_folder = "/path/to/csv/gating/files",  # Directory with gating CSV files
  output_dir = "/path/to/output",             # All results saved here

  # -- Conditions --
  conditions_csv   = "/path/to/conditions.csv",  # CSV with columns: sample, condition
  condition_col    = "condition",                 # Column name for condition labels
  sample_col       = "sample",                    # Column name for sample identifiers
  condition_levels = c("Control", "Treatment"),   # Two-level condition for DA testing

  # ========================= OPTIONAL (sensible defaults) ====================

  # -- Analysis --
  # Choose which statistical method to use:
  #   "glmm"       — Negative binomial GLMM for differential abundance (default)
  #   "enrichment" — Fisher's exact test or mixed-effects logistic regression
  analysis_method = "glmm",

  # -- Cluster Splitting --
  cluster_regex    = "Clustering annotated",  # Regex to identify cluster/gating columns
  min_cluster_size = 10,                      # Drop clusters with fewer cells

  # -- Statistics --
  stat_method     = "mixed_effects",  # "mixed_effects" or "pooled"
  p_adjust_method = "BH",             # p-value correction: "BH", "bonferroni", "holm"

  # -- Markers (optional) --
  # If NULL, all non-gating numeric columns are used as markers.
  markers = NULL,

  # -- Output --
  output_formats = list(pdf = TRUE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # -- Reproducibility --
  seed = 42
)
