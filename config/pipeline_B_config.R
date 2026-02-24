# =============================================================================
# Pipeline B: CSV Expression -> Quick Statistics
# =============================================================================
# Edit the values below, then run:
#   cd pipelines
#   Rscript pipeline_B_csv_expression.R --config ../config/pipeline_B_config.R
# =============================================================================

source(file.path(dirname(normalizePath(sys.frame(1)$ofile %||%
  Filter(Negate(is.null), lapply(rev(seq_len(sys.nframe())),
    function(i) sys.frame(i)$ofile))[[1]])), "defaults.R"))

params <- pipeline_defaults(

  # === EDIT THESE ============================================================
  input_dir           = "../test_data/csv_expression",
  output_dir          = "../output/pipeline_B_results",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),

  # === STATISTICS ============================================================
  test_type           = "lmer"          # "lmer", "wilcox", "kruskal"

  # All other parameters use sensible defaults. See config/defaults.R.
)
