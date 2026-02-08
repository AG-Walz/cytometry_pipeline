# =============================================================================
# Example Configuration File for Cytometry Pipeline
# =============================================================================
# Copy this file to config/study_params/my_study.R and edit as needed.
# See config/Parameters_template.R for the full list of parameters.
# =============================================================================

params <- list()

# --- Paths ---
params$input_dir     <- "test_data/examples/"
params$project_name  <- "Example_Study"
params$output_dir    <- "output/example/"

# --- Sample metadata ---
params$sample_names  <- c("sample_001", "sample_002", "sample_003",
                          "sample_004", "sample_005", "sample_006")
params$batch_names   <- c("Batch1", "Batch1", "Batch1",
                          "Batch2", "Batch2", "Batch2")

# --- Markers ---
params$markers_of_interest <- c("CD4", "CD8", "CD19", "CD56", "HLA-DR",
                                 "CD45", "CD3", "CD14", "CD16", "CD25")

# --- Preprocessing ---
params$qc_type        <- "PeacoQC"
params$transform_type <- "arcsinh"
params$pp_order       <- "QT"
params$quant_norm     <- FALSE
params$scaling        <- FALSE
params$save_pp_fcs    <- TRUE
params$save_pp_csv    <- TRUE

# --- Downsampling ---
params$target_ncell_subsampling <- 10000
params$subsampling_method       <- "Maximal"

# --- Clustering ---
params$n_neighbors <- 15
params$resolution  <- 1.0

# --- Analysis ---
params$condition_levels <- c("Control", "Treatment")
params$stat_method      <- "mixed_effects"
params$p_adjust_method  <- "BH"
params$alpha_level      <- 0.05
params$min_cluster_size <- 10

# --- Output ---
params$output_formats <- list(pdf = TRUE, png = TRUE, svg = FALSE)
params$plot_width     <- 10
params$plot_height    <- 8
params$png_dpi        <- 300

# --- Memory ---
params$chunk_size         <- 50000
params$max_memory_gb      <- 4
params$max_cells_per_file <- 500000

# --- Reproducibility ---
params$seed <- 42

# --- Derived (auto-generated) ---
params$sample_batch_dict <- setNames(params$batch_names, params$sample_names)
