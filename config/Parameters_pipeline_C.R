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
  output_formats = list(pdf = TRUE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,
  png_dpi        = 300,

  # ========================= Seed ============================================
  seed = 42
)

# Derived values
params$sample_batch_dict <- stats::setNames(params$batch_names, params$sample_names)
