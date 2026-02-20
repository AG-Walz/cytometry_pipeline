# =============================================================================
# config/Parameters_template.R — Pipeline Parameter Template
# =============================================================================
# Copy this file and modify for your study. Pass as argument to pipeline scripts.
# =============================================================================

params <- list(
  # ========================= General =========================================
  input_dir        = "../test_data/",
  project_name     = "MyStudy",
  sample_names     = c("Sample1", "Sample2", "Sample3"),
  batch_names      = c("1", "1", "2"),
  patient_names    = c("P1", "P2", "P3"),

  # ========================= Markers =========================================
  markers_of_interest = c(
    "CD4", "CD8", "CD19", "CD56", "HLA-DR"
    # Add your markers here
  ),

  # Marker renaming dictionary (FCS name -> display name)
  markername_replacing = list(
    # `PD1 (CD279)` = "PD1",
    # `Fas (CD95)` = "CD95"
  ),

  # ========================= Loading =========================================
  truncation                   = FALSE,
  min_limit                    = NULL,
  check_panels                 = TRUE,

  # ========================= Plotting ========================================
  plot_density_raw_per_file    = FALSE,
  plot_density_per_samplenames = FALSE,
  plot_density_per_batch       = FALSE,
  plot_umap_raw_per_file       = FALSE,

  # ========================= Downsampling ====================================
  target_ncell_subsampling = 10000,
  subsampling_method       = "Maximal",   # "Consistent" or "Maximal"
  save_downsampling_fcs    = FALSE,

  # ========================= Homogenization ==================================
  homogenize      = "before",   # "before", "after", or "None"
  extract_ro_csv  = FALSE,
  extract_ro_fcs  = FALSE,

  # ========================= Preprocessing ===================================
  transform_type  = "arcsinh",   # "arcsinh", "logicle", or "biexponential"
  qc_type         = "PeacoQC",  # "None", "flowAI", "PeacoQC"
  pp_order        = "QT",       # "QT" or "TQ"
  quant_norm      = FALSE,
  scaling         = FALSE,
  save_pp_fcs     = TRUE,
  save_pp_csv     = TRUE,
  normalize       = FALSE,
  norm_covar      = "batch",

  # ========================= Analysis ========================================
  cluster_regex              = "Clustering annotated",
  condition_levels           = c("Control", "Treatment"),
  min_cluster_size           = 10,
  p_adjust_method            = "BH",
  alpha_level                = 0.05,
  stat_method                = "mixed_effects",  # "pooled" or "mixed_effects"

  # ========================= Memory Management ===============================
  chunk_size       = 50000,
  max_memory_gb    = 4,
  max_cells_per_file = 500000,
  sample_large_files = TRUE,

  # ========================= Output Formats ==================================
  output_formats = list(pdf = TRUE, png = TRUE, svg = TRUE),
  plot_width     = 14,
  plot_height    = 10,
  png_dpi        = 200,

  # ========================= Abundance Analysis =============================
  exclude_patterns = c("^AI_Cluster_", "^Unlabeled_\\d+", "^Clustering raw$"),
  cytolution_metadata_file = NULL,  # Path to metadata CSV for grouped abundance plots

  # ========================= Method Selection =================================
  clustering_method   = "leiden",       # "leiden", "flowsom", "phenograph", "dbscan"
  batch_method        = "none",         # "harmony", "cytonorm", "combat", "none"
  dimred_method       = "umap",         # "umap", "tsne", "diffusion_map"

  # Clustering-specific
  n_clusters          = 20,             # FlowSOM meta-clusters
  grid_xdim           = 10,             # FlowSOM SOM grid
  grid_ydim           = 10,
  dbscan_eps          = NULL,           # DBSCAN eps (NULL = auto)
  dbscan_min_pts      = 5,

  # Batch correction-specific
  batch_do_pca        = FALSE,          # Harmony: run PCA before correction
  cytonorm_train_indices = NULL,        # CytoNorm: NULL = use all frames
  cytonorm_n_clusters = 10,
  combat_method       = "parametric",   # ComBat: "parametric" or "nonparametric"

  # DimRed-specific
  umap_min_dist       = 0.1,
  tsne_perplexity     = 30,
  tsne_max_iter       = 1000,
  dm_n_components     = 3,
  dm_sigma            = "local",

  # ========================= Seed ============================================
  seed = 42
)

# Derived values
params$sample_batch_dict <- stats::setNames(params$batch_names, params$sample_names)
