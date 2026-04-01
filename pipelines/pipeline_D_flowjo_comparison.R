#!/usr/bin/env Rscript
# =============================================================================
# Pipeline D: FlowJo Workspace -> Gating Comparison
# =============================================================================
# Flow: Load WSP -> Compare gatings -> Overlap statistics -> Venn/UpSet
# =============================================================================

# --- Init (config + shared modules) ---
base_dir   <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "pipeline_utils.R"))
env <- init_pipeline("Pipeline D: FlowJo Workspace -> Gating Comparison",
                     "config/pipeline_D_config.R")
source(file.path(env$shared_dir, "cluster_utils.R"))

# Validate required parameters
if (is.null(params$input_dir) && is.null(params$fcs_folder)) {
  stop("input_dir (or fcs_folder) is required in config")
}
if (is.null(params$output_dir)) {
  stop("output_dir is required in config")
}
# Validate output directory can be created
tryCatch(make_outdir(params$output_dir),
  error = function(e) stop("Cannot create output directory: ", params$output_dir, " (", e$message, ")")
)

# Check for placeholder paths in config
if (any(grepl("/path/to/", unlist(params)))) {
  stop("Config contains placeholder paths ('/path/to/...'). ",
       "Edit your config file with real paths before running Pipeline D.")
}

# Step 1: Load FlowJo workspaces
cat("\n--- Step 1: Loading FlowJo workspaces ---\n")
source(file.path(env$base_dir, "src", "loading", "load_flowjo_wsp.R"))
pop_strings_list <- load_multiple_workspaces(
  wsp_paths = wsp_paths,
  names     = custom_name
)

# Step 2: Compare gatings
cat("\n--- Step 2: Comparing gating strategies ---\n")
source(file.path(env$base_dir, "src", "metadata_qc", "compare_gatings.R"))
label_lists <- compare_gatings(
  pop_strings_list = pop_strings_list,
  target_fcs       = target_fcs,
  label_dict       = label_dict,
  custom_name      = custom_name,
  output_dir       = file.path(params$output_dir, "01_gating_comparison"),
  idx_labels       = if (exists("idx_labels", inherits = FALSE)) idx_labels else NULL,
  scale_size       = TRUE
)

# Step 3: Overlap statistics
cat("\n--- Step 3: Computing overlap statistics ---\n")
source(file.path(env$base_dir, "src", "metadata_qc", "overlap_statistics.R"))
kappa_results <- compute_fleiss_kappa(
  label_lists = label_lists,
  output_dir  = file.path(params$output_dir, "02_overlap_stats")
)

finish_pipeline("Pipeline D", env$start_time)
