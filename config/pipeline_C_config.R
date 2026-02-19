# =============================================================================
# config/pipeline_C_config.R — Pipeline C: CSV Gating -> Cluster Analysis
# =============================================================================
# Flow: Load gating CSV -> Split clusters -> Enrichment -> Barplots/Abundance
#
# Usage:
#   cd pipelines && Rscript pipeline_C_csv_gating.R --config ../config/pipeline_C_config.R
# =============================================================================

# Compute project root relative to this config file
.config_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) {
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) return(normalizePath(f))
    }
    stop("Cannot determine config file path")
  }
)
base_dir <- dirname(dirname(.config_path))

params <- list(

  # ========================= Input / Output ==================================
  csv_folder     = file.path(base_dir, "test_data", "csv_gating"),
  output_dir     = file.path(base_dir, "output", "pipeline_C_results"),

  # -- Conditions --
  conditions_csv   = file.path(base_dir, "test_data", "csv_gating_conditions.csv"),
  condition_levels = c("Control", "Treatment"),

  # ========================= Analysis ========================================
  analysis_method  = "enrichment",
  stat_method      = "pooled",
  cluster_regex    = "Clustering annotated",
  min_cluster_size = 5,
  p_adjust_method  = "BH",

  # -- Gating loading --
  gating_pattern   = "\\(Clustering annotated\\)$|^Clustering raw$",
  processing_mode  = "manual",

  # -- Abundance --
  exclude_patterns         = c(),
  cytolution_metadata_file = NULL,

  # ========================= Output ==========================================
  output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # ========================= Reproducibility =================================
  seed = 42
)
