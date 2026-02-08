# =============================================================================
# Standalone config: src/loading/load_flowjo_wsp.R
# =============================================================================
# Loads a FlowJo workspace (.wsp) and extracts gating sets, population
# assignments, and population statistics.
#
# Usage:
#   source("config/standalone/load_flowjo_wsp_config.R")
#   source("src/loading/load_flowjo_wsp.R")
#   wsp_data <- load_flowjo_wsp(
#     wsp_path = wsp_path,
#     fcs_dir  = fcs_dir,
#     name     = name
#   )
#   # Explore results:
#   get_population_stats(wsp_data$gating_set)
#   pop_strings <- get_population_strings_list(wsp_data)
#   summarize_population_counts(pop_strings)
# =============================================================================

# ========================= EDIT THESE ========================================

wsp_path <- "/path/to/workspace.wsp"

# ========================= OPTIONAL ==========================================

fcs_dir <- NULL   # Directory with FCS files (NULL = same directory as WSP)
name    <- NULL   # Custom name for this workspace (NULL = derived from filename)
