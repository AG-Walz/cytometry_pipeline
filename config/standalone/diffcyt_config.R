# =============================================================================
# Standalone config: src/analysis/diffcyt_analysis.R
# =============================================================================
# Runs differential analysis using the diffcyt framework. Supports differential
# abundance (DA) and differential state (DS) testing with edgeR or limma backends.
#
# Usage:
#   source("config/standalone/diffcyt_config.R")
#   source("src/analysis/diffcyt_analysis.R")
#   result <- run_diffcyt_analysis(
#     fcs_folder          = fcs_folder,
#     experiment_info_csv = experiment_info_csv,
#     output_dir          = output_dir,
#     markers_of_interest = markers_of_interest,
#     analysis_type       = analysis_type,
#     method_DA           = method_DA,
#     n_clusters          = n_clusters,
#     seed                = seed
#   )
#
# Or from CLI:
#   Rscript src/analysis/diffcyt_analysis.R <fcs_folder> <experiment_info_csv> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

fcs_folder          <- "/path/to/fcs/files"
experiment_info_csv <- "/path/to/experiment_info.csv"  # CSV: sample_id, group_id, patient_id
output_dir          <- "/path/to/output/diffcyt"
markers_of_interest <- c("CD3", "CD4", "CD8", "CD19", "CD56")

# ========================= OPTIONAL ==========================================

# Analysis type: "DA" (differential abundance) or "DS" (differential state)
analysis_type <- "DA"

# Method for DA testing: "diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"
method_DA <- "diffcyt-DA-edgeR"

# Number of clusters for FlowSOM step within diffcyt
n_clusters <- 20

seed <- 42
