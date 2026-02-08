# =============================================================================
# config/pipeline_B_config.R — Pipeline B: CSV Expression → Statistics
# =============================================================================
# Flow: Load CSV expression matrices → Marker statistics → Density/Heatmap
#
# Use this pipeline when you already have expression data exported as CSV files
# (e.g., from FlowJo, Cytolution, or another tool) and want to run marker-level
# statistics and visualizations without any FCS processing.
#
# Usage:
#   Rscript pipelines/pipeline_B_csv_expression.R --config config/pipeline_B_config.R
# Or in R:
#   source("config/pipeline_B_config.R")
#   source("pipelines/pipeline_B_csv_expression.R")
# =============================================================================

params <- list(

  # ========================= EDIT THESE ======================================

  # -- Input / Output --
  csv_folder = "/path/to/csv/expression/files",  # Directory with expression CSV files
  output_dir = "/path/to/output",                 # All results saved here

  # -- Markers --
  # Column names in your CSVs that correspond to markers.
  # Only these columns will be used for statistics and plots.
  markers = c(
    "CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"
  ),

  # -- Grouping --
  # Column names in your CSVs used for grouping in statistical tests and plots.
  condition_col = "condition",   # Column with group/condition labels
  sample_col    = "sample",      # Column with sample identifiers

  # ========================= OPTIONAL (sensible defaults) ====================

  # -- Output --
  output_formats = list(pdf = TRUE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # -- Reproducibility --
  seed = 42
)
