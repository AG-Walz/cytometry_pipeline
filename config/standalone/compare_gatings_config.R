# =============================================================================
# Standalone config: metadata_qc/compare_gatings.R
# =============================================================================
# Compares gating strategies from multiple FlowJo workspaces for a single
# FCS file. Produces Venn/Euler diagrams per population showing agreement
# between analysts or strategies.
#
# Typically used after loading workspaces with load_flowjo_wsp.R.
#
# Usage:
#   source("config/standalone/compare_gatings_config.R")
#   source("loading/load_flowjo_wsp.R")
#   source("metadata_qc/compare_gatings.R")
#
#   # Load each workspace
#   pop_strings_list <- list()
#   for (i in seq_along(wsp_paths)) {
#     wsp <- load_flowjo_wsp(wsp_path = wsp_paths[i])
#     pop_strings_list[[custom_name[i]]] <- get_population_strings_list(wsp)
#   }
#
#   # Compare
#   label_lists <- compare_gatings(
#     pop_strings_list = pop_strings_list,
#     target_fcs       = target_fcs,
#     label_dict       = label_dict,
#     custom_name      = custom_name,
#     output_dir       = output_dir
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

wsp_paths <- c(
  "/path/to/workspace_analyst1.wsp",
  "/path/to/workspace_analyst2.wsp"
)

custom_name <- c("Analyst1", "Analyst2")
target_fcs  <- "Sample001_Unmixed.fcs"
output_dir  <- "/path/to/output/compare_gatings"

# Label dictionary: map FlowJo gate paths → simplified population names.
# Set to "Removed" to exclude from comparison.
label_dict <- c(
  "root"                                     = "Removed",
  "/Singlets"                                = "Removed",
  "/Singlets/Live/Lymphocytes/CD3+"          = "T cells",
  "/Singlets/Live/Lymphocytes/CD3+/CD4+"     = "CD4+ T cells",
  "/Singlets/Live/Lymphocytes/CD3+/CD8+"     = "CD8+ T cells",
  "/Singlets/Live/Lymphocytes/CD3-/CD19+"    = "B cells",
  "/Singlets/Live/Lymphocytes/CD3-/CD56+"    = "NK cells"
)

# ========================= OPTIONAL ==========================================

idx_labels <- NULL     # Per-workspace FCS indices (NULL if identical)
scale_size <- TRUE     # Scale Euler diagram areas proportionally
