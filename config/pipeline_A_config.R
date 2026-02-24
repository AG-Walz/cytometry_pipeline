# =============================================================================
# Pipeline A: Raw FCS -> Full Analysis
# =============================================================================
# Edit the values below, then run:
#   cd pipelines
#   Rscript pipeline_A_raw_fcs.R --config ../config/pipeline_A_config.R
# =============================================================================

# Load defaults (provides sensible values for everything not listed here)
source(file.path(dirname(normalizePath(sys.frame(1)$ofile %||%
  Filter(Negate(is.null), lapply(rev(seq_len(sys.nframe())),
    function(i) sys.frame(i)$ofile))[[1]])), "defaults.R"))

params <- pipeline_defaults(

  # === EDIT THESE ============================================================
  input_dir           = "../test_data/fcs_samples",
  output_dir          = "../output/pipeline_A_results",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),

  # === CONDITIONS (for differential abundance) ===============================
  conditions_csv      = "../test_data/fcs_samples/conditions.csv",
  condition_levels    = c("Control", "Treatment"),

  # === METHODS (change these to swap algorithms) =============================
  clustering_method   = "leiden",       # "leiden", "flowsom", "phenograph", "dbscan"
  batch_method        = "none",         # "harmony", "cytonorm", "combat", "none"
  dimred_method       = "umap",         # "umap", "tsne", "diffusion_map"
  transform_type      = "arcsinh"       # "arcsinh", "logicle", "biexponential"

  # All other parameters (clustering resolution, downsample count, output
  # formats, etc.) use sensible defaults. See config/defaults.R to override.
)
