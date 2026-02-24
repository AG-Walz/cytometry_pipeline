# =============================================================================
# config/Parameters_template.R — Full Parameter Template
# =============================================================================
# Copy this file for your study and edit the values below.
# Only the "EDIT THESE" section is required — everything else has defaults.
#
# For a quick start, see the pipeline-specific configs instead:
#   pipeline_A_config.R, pipeline_B_config.R, etc.
#
# For the complete list of available defaults, see config/defaults.R.
# =============================================================================

source(file.path(dirname(normalizePath(sys.frame(1)$ofile %||%
  Filter(Negate(is.null), lapply(rev(seq_len(sys.nframe())),
    function(i) sys.frame(i)$ofile))[[1]])), "defaults.R"))

params <- pipeline_defaults(

  # ========================= EDIT THESE (required) ============================
  input_dir           = "path/to/your/data",
  output_dir          = "path/to/results",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),

  # ========================= Conditions =======================================
  conditions_csv      = NULL,                # Path to conditions CSV (sample -> group)
  condition_levels    = c("Control", "Treatment"),
  sample_names        = NULL,                # NULL = auto-detect from filenames
  batch_names         = NULL,                # NULL = no batch correction
  patient_names       = NULL,

  # ========================= Methods ==========================================
  clustering_method   = "leiden",            # "leiden", "flowsom", "phenograph", "dbscan"
  batch_method        = "none",              # "harmony", "cytonorm", "combat", "none"
  dimred_method       = "umap",              # "umap", "tsne", "diffusion_map"
  transform_type      = "arcsinh",           # "arcsinh", "logicle", "biexponential"
  analysis_method     = "enrichment",        # "enrichment" or "abundance"

  # ========================= Marker Renaming ==================================
  markername_replacing = list(
    # `PD1 (CD279)` = "PD1",
    # `Fas (CD95)`  = "CD95"
  ),

  # ========================= Tuning (optional) ================================
  # Uncomment and change only if you need non-default values:
  # n_neighbors       = 15,
  # resolution        = 0.5,
  # downsample_n      = 10000,
  # downsample_method = "Consistent",
  # qc_type           = "None",
  # seed              = 42

  seed = 42
)
