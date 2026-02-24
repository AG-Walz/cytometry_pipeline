# =============================================================================
# Pipeline C: CSV Gating -> Cluster Abundance Analysis
# =============================================================================
# Edit the values below, then run:
#   cd pipelines
#   Rscript pipeline_C_csv_gating.R --config ../config/pipeline_C_config.R
# =============================================================================

source(file.path(dirname(normalizePath(sys.frame(1)$ofile %||%
  Filter(Negate(is.null), lapply(rev(seq_len(sys.nframe())),
    function(i) sys.frame(i)$ofile))[[1]])), "defaults.R"))

params <- pipeline_defaults(

  # === EDIT THESE ============================================================
  input_dir           = "../test_data/csv_gating",
  output_dir          = "../output/pipeline_C_results",

  # === CONDITIONS ============================================================
  conditions_csv      = "../test_data/csv_gating_conditions.csv",
  condition_levels    = c("Control", "Treatment"),

  # === ANALYSIS METHOD =======================================================
  analysis_method     = "enrichment"    # "enrichment" (Fisher) or "abundance" (GLMM)

  # All other parameters use sensible defaults. See config/defaults.R.
)
