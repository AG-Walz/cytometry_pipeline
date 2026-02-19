# =============================================================================
# config/pipeline_B_config.R — Pipeline B: CSV Expression -> Statistics
# =============================================================================
# Flow: Load CSV expression matrices -> Marker statistics -> Density/Heatmap
#
# Usage:
#   cd pipelines && Rscript pipeline_B_csv_expression.R --config ../config/pipeline_B_config.R
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
  csv_folder = file.path(base_dir, "test_data", "csv_expression"),
  output_dir = file.path(base_dir, "output", "pipeline_B_results"),

  # -- Markers --
  markers = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),

  # -- Statistics --
  test_type = "lmer",

  # ========================= Output ==========================================
  output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # ========================= Reproducibility =================================
  seed = 42
)
