# =============================================================================
# config/pipeline_A_config.R — Pipeline A: Raw FCS -> Full Analysis
# =============================================================================
# Flow: Load FCS -> QC -> Transform -> Homogenize -> Downsample -> Cluster -> Analyze -> Visualize
#
# Usage:
#   cd pipelines && Rscript pipeline_A_raw_fcs.R --config ../config/pipeline_A_config.R
# =============================================================================

# Compute project root relative to this config file
.config_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) {
    # Fallback: search the call stack for source()'d file path
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) return(normalizePath(f))
    }
    stop("Cannot determine config file path")
  }
)
base_dir <- dirname(dirname(.config_path))  # config/ -> project root

params <- list(

  # ========================= Input / Output ==================================
  fcs_folder  = file.path(base_dir, "test_data", "fcs_samples"),
  output_dir  = file.path(base_dir, "output", "pipeline_A_results"),

  # -- Markers --
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),

  # -- Conditions --
  conditions_csv   = file.path(base_dir, "test_data", "fcs_samples", "conditions.csv"),
  condition_levels = c("Control", "Treatment"),

  # ========================= Loading =========================================
  truncation = FALSE,

  # ========================= Preprocessing ===================================
  pp_order       = "QT",
  qc_type        = "None",
  transform_type = "arcsinh",
  offset_list    = rep(0, 6),   # one per marker (default: no offset)
  width_list     = rep(5, 6),   # cofactor per marker (5 = mass cytometry default)
  quant_norm     = FALSE,
  scaling        = FALSE,
  save_pp_fcs    = FALSE,
  save_pp_csv    = FALSE,

  # ========================= Homogenization ==================================
  markername_replacing = list(),

  # ========================= Downsampling ====================================
  downsample_method = "Consistent",
  downsample_n      = 100,

  # ========================= Clustering (Leiden) =============================
  n_neighbors = 15,
  resolution  = 0.5,

  # ========================= Analysis ========================================
  cluster_regex    = "Clustering annotated",
  min_cluster_size = 5,
  p_adjust_method  = "BH",
  stat_method      = "pooled",

  # ========================= Output ==========================================
  output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # ========================= Reproducibility =================================
  seed = 42
)
