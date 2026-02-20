#!/usr/bin/env Rscript
# =============================================================================
# pipelines/pipeline_D_flowjo_comparison.R — FlowJo Gating Comparison Pipeline
# =============================================================================
# Flow: Load WSP → Compare gatings → Overlap statistics → Venn/UpSet
# =============================================================================

cat("=== Pipeline D: FlowJo Workspace → Gating Comparison ===\n")
start_time <- Sys.time()

# --- Configuration ---
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0 && args[1] == "--config") args[2] else NULL
if (!is.null(config_file)) {
  source(config_file)
} else {
  stop("Pipeline D requires a study-specific config. Use: --config config/study_params/Compare_gatings_*.R")
}

base_dir <- dirname(getwd())
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

# Step 1: Load FlowJo workspaces
cat("\n--- Step 1: Loading FlowJo workspaces ---\n")
source(file.path(base_dir, "src", "loading", "load_flowjo_wsp.R"))

pop_strings_list <- load_multiple_workspaces(
  wsp_paths = wsp_paths,
  names     = custom_name
)

# Step 2: Compare gatings
cat("\n--- Step 2: Comparing gating strategies ---\n")
source(file.path(base_dir, "src", "metadata_qc", "compare_gatings.R"))
label_lists <- compare_gatings(
  pop_strings_list = pop_strings_list,
  target_fcs = target_fcs,
  label_dict = label_dict,
  custom_name = custom_name,
  output_dir = file.path(params$output_dir, "01_gating_comparison"),
  idx_labels = if (exists("idx_labels", inherits = FALSE)) idx_labels else NULL,
  scale_size = TRUE
)

# Step 3: Overlap statistics
cat("\n--- Step 3: Computing overlap statistics ---\n")
source(file.path(base_dir, "src", "metadata_qc", "overlap_statistics.R"))
kappa_results <- compute_fleiss_kappa(
  label_lists = label_lists,
  output_dir = file.path(params$output_dir, "02_overlap_stats")
)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== Pipeline D complete (%.1f minutes) ===\n", elapsed))
