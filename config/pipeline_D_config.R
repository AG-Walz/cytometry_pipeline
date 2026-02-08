# =============================================================================
# config/pipeline_D_config.R — Pipeline D: FlowJo Gating Comparison
# =============================================================================
# Flow: Load WSP → Compare gating strategies → Overlap statistics → Venn/UpSet
#
# Use this pipeline to compare gating strategies from multiple FlowJo workspaces.
# For example, compare how different analysts gated the same samples, or compare
# broad vs. narrow gating strategies. Produces Venn/Euler diagrams per population
# and Fleiss' kappa inter-rater agreement scores.
#
# Usage:
#   Rscript pipelines/pipeline_D_flowjo_comparison.R --config config/pipeline_D_config.R
# Or in R:
#   source("config/pipeline_D_config.R")
#   source("pipelines/pipeline_D_flowjo_comparison.R")
#
# NOTE: Pipeline D requires a config file (no default template works).
# =============================================================================

# ========================= EDIT THESE ========================================

# -- FlowJo Workspace Paths --
# List all .wsp files to compare. Each workspace should contain gating for the
# same set of FCS files.
wsp_paths <- c(
  "/path/to/workspace_analyst1.wsp",
  "/path/to/workspace_analyst2.wsp"
)

# -- Workspace Labels --
# Short names identifying each workspace (same order as wsp_paths).
# These appear as labels in Venn/Euler diagrams and overlap tables.
custom_name <- c("Analyst1", "Analyst2")

# -- Target FCS --
# The specific FCS filename to compare across workspaces.
# All workspaces must contain gating for this file.
target_fcs <- "Sample001_Unmixed.fcs"

# -- Label Dictionary --
# Maps FlowJo gate paths to simplified population labels.
# Gates mapped to "Removed" are excluded from the comparison.
# Gates mapped to the same label across different hierarchy paths are unified.
#
# To find your gate paths, load a workspace and inspect:
#   wsp <- load_flowjo_wsp("/path/to/workspace.wsp")
#   get_population_stats(wsp$gating_set)
label_dict <- c(
  # Exclude pre-gating hierarchy
  "root"                                              = "Removed",
  "/Singlets"                                         = "Removed",
  "/Singlets/Live"                                    = "Removed",
  "/Singlets/Live/Lymphocytes"                        = "Removed",

  # Map gate paths to population labels
  # Analyst 1 paths
  "/Singlets/Live/Lymphocytes/CD3+"                   = "T cells",
  "/Singlets/Live/Lymphocytes/CD3+/CD4+"              = "CD4+ T cells",
  "/Singlets/Live/Lymphocytes/CD3+/CD8+"              = "CD8+ T cells",
  "/Singlets/Live/Lymphocytes/CD3-/CD19+"             = "B cells",
  "/Singlets/Live/Lymphocytes/CD3-/CD56+"             = "NK cells",

  # Analyst 2 paths (may differ in naming)
  "/Single Cells/Live/Lymphos/CD3 positive"           = "T cells",
  "/Single Cells/Live/Lymphos/CD3 positive/CD4 pos"   = "CD4+ T cells",
  "/Single Cells/Live/Lymphos/CD3 positive/CD8 pos"   = "CD8+ T cells",
  "/Single Cells/Live/Lymphos/CD3 neg/CD19 pos"       = "B cells",
  "/Single Cells/Live/Lymphos/CD3 neg/CD56 pos"       = "NK cells"
)

# -- Index Labels (optional) --
# If your workspaces use different FCS indexing, specify per-workspace indices.
# Set to NULL if all workspaces index identically.
idx_labels <- NULL

# -- General params --
params <- list(
  output_dir = "/path/to/output"
)
