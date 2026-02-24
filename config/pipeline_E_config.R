# =============================================================================
# Pipeline E: Preprocessed FCS -> Fast UMAP / Clustering
# =============================================================================
# Edit the values below, then run:
#   cd pipelines
#   Rscript pipeline_E_preprocessed_umap.R --config ../config/pipeline_E_config.R
# =============================================================================

source(file.path(dirname(normalizePath(sys.frame(1)$ofile %||%
  Filter(Negate(is.null), lapply(rev(seq_len(sys.nframe())),
    function(i) sys.frame(i)$ofile))[[1]])), "defaults.R"))

params <- pipeline_defaults(

  # === EDIT THESE ============================================================
  input_dir           = "../test_data/fcs_samples",
  output_dir          = "../output/pipeline_E_results",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),

  # === METHODS (change these to swap algorithms) =============================
  clustering_method   = "leiden",       # "leiden", "flowsom", "phenograph", "dbscan"
  dimred_method       = "umap"          # "umap", "tsne", "diffusion_map"

  # All other parameters use sensible defaults. See config/defaults.R.
)
