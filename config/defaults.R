# =============================================================================
# config/defaults.R — Sensible defaults for all pipeline parameters
# =============================================================================
# Users only need to set input_dir, output_dir, and markers in their config.
# Everything else here provides a reasonable starting point.
#
# Usage (in a pipeline config):
#   source("config/defaults.R")
#   params <- pipeline_defaults(
#     input_dir  = "path/to/my/data",
#     output_dir = "path/to/results",
#     markers_of_interest = c("CD3", "CD4", "CD8")
#   )
# =============================================================================

#' Build a complete params list by merging user overrides with defaults
#'
#' @param ... Named arguments that override the defaults (input_dir, output_dir, etc.)
#' @return A complete params list with all pipeline parameters
pipeline_defaults <- function(...) {
  user <- list(...)

  defaults <- list(

    # === Paths (no defaults — user MUST provide these) =======================
    input_dir           = NULL,
    output_dir          = NULL,

    # === Markers (no default — user MUST provide) ============================
    markers_of_interest = NULL,

    # === Methods (swappable via config) ======================================
    clustering_method   = "leiden",          # "leiden", "flowsom", "phenograph", "dbscan"
    batch_method        = "none",            # "harmony", "cytonorm", "combat", "none"
    dimred_method       = "umap",            # "umap", "tsne", "diffusion_map"

    # === Conditions / Metadata ===============================================
    conditions_csv      = NULL,
    condition_levels    = c("Control", "Treatment"),
    sample_names        = NULL,
    batch_names         = NULL,
    patient_names       = NULL,

    # === Loading =============================================================
    truncation          = FALSE,
    min_limit           = NULL,
    check_panels        = TRUE,

    # === Preprocessing =======================================================
    pp_order            = "QT",              # "QT" (QC then Transform) or "TQ"
    qc_type             = "None",            # "None", "flowAI", "PeacoQC"
    transform_type      = "arcsinh",         # "arcsinh", "logicle", "biexponential", "log"
    offset_list         = NULL,
    width_list          = NULL,
    quant_norm          = FALSE,
    scaling             = FALSE,
    save_pp_fcs         = FALSE,
    save_pp_csv         = FALSE,

    # === Homogenization ======================================================
    markername_replacing = list(),

    # === Downsampling ========================================================
    downsample_method   = "Consistent",      # "Consistent" or "Maximal"
    downsample_n        = 10000,

    # === Clustering parameters ===============================================
    n_neighbors         = 15,                # Leiden / UMAP
    resolution          = 0.5,               # Leiden
    n_clusters          = 20,                # FlowSOM meta-clusters
    grid_xdim           = 10,                # FlowSOM SOM grid
    grid_ydim           = 10,
    dbscan_eps          = NULL,              # DBSCAN (NULL = auto-estimate)
    dbscan_min_pts      = 5,

    # === Batch correction parameters =========================================
    batch_do_pca            = FALSE,
    cytonorm_train_indices  = NULL,
    cytonorm_n_clusters     = 10,
    combat_method           = "parametric",  # "parametric" or "nonparametric"

    # === Dimensionality reduction parameters =================================
    umap_min_dist       = 0.1,
    tsne_perplexity     = 30,
    tsne_max_iter       = 1000,
    dm_n_components     = 3,
    dm_sigma            = "local",

    # === Analysis ============================================================
    analysis_method     = "enrichment",      # "enrichment" or "abundance"
    cluster_regex       = "Clustering annotated",
    gating_pattern      = "\\(Clustering annotated\\)$|^Clustering raw$",
    processing_mode     = "manual",
    min_cluster_size    = 5,
    p_adjust_method     = "BH",
    alpha_level         = 0.05,
    stat_method         = "pooled",          # "pooled" or "mixed_effects"
    test_type           = "lmer",
    exclude_patterns    = c(),
    cytolution_metadata_file = NULL,

    # === Output ==============================================================
    output_formats      = list(pdf = FALSE, png = TRUE, svg = FALSE),
    plot_width          = 10,
    plot_height         = 8,

    # === Memory ==============================================================
    chunk_size          = 50000,
    max_memory_gb       = 4,
    max_cells_per_file  = 500000,

    # === Reproducibility =====================================================
    seed                = 42
  )

  # Merge: user values override defaults
  for (nm in names(user)) {
    defaults[[nm]] <- user[[nm]]
  }

  # Convenience aliases used by different pipelines
  if (is.null(defaults$fcs_folder))  defaults$fcs_folder  <- defaults$input_dir
  if (is.null(defaults$csv_folder))  defaults$csv_folder  <- defaults$input_dir
  if (is.null(defaults$markers))     defaults$markers     <- defaults$markers_of_interest

  # Derived values
  if (!is.null(defaults$sample_names) && !is.null(defaults$batch_names)) {
    defaults$sample_batch_dict <- stats::setNames(defaults$batch_names, defaults$sample_names)
  }

  defaults
}
