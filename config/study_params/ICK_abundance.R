# =============================================================================
# Study config: ICK — Cluster Abundance Analysis
# =============================================================================
# Input: Preprocessed CSVs with expression data + binary clustering columns
# Goal:  Build a sample-by-cluster tracking matrix and produce abundance plots
#
# Usage:
#   cd pipelines
#   Rscript run_ICK_abundance.R
# =============================================================================

# ========================= PATHS ============================================

input_dir  <- "/mnt/volume/nacho/data/ICK"
output_dir <- "/mnt/volume/nacho/data/ICK/abundance_output"

# ========================= CLUSTERING COLUMNS ===============================

# Regex to identify clustering annotation columns
cluster_regex <- "\\(Clustering annotated\\)$"

# Regex patterns to EXCLUDE from the tracking matrix (fine-grained sub-clusters)
exclude_patterns <- c("^AI_Cluster_", "^Unlabeled_\\d+", "^Clustering raw$")

# ========================= MARKERS ==========================================

# Expression markers to keep (NULL = auto-detect all non-scatter/non-gating cols)
markers_of_interest <- c(
  "CD366 | PE.Fire810.A",
  "CD45RO | BUV395.A",
  "CD16 | BUV496.A",
  "CD62L | BUV563.A",
  "CD137 | BUV615.A",
  "CD154 | BUV661.A",
  "CD56 | BUV737.A",
  "CD8 | BUV805.A",
  "LOX1 | BV421.A",
  "CD127 | BV480.A",
  "CD3 | BV510.A",
  "CD11B | BV570.A",
  "TIGIT | BV605.A",
  "CD28 | BV650.A",
  "TCRGD | BV711.A",
  "CD14 | BV785.A",
  "CD57 | FITC.A",
  "HLADR | Spark.Blue.550.A",
  "CD66B | PerCP.Cy5.5.A",
  "CD279 | RB705.A",
  "CD69 | RB744.A",
  "KLRG1 | PerCP.Fire.780.A",
  "CD4 | Spark.YG.581.A",
  "CD152 | PE.CF594.A",
  "CD223 | PE.Cy5.A",
  "CD25 | PE.Fire.700.A",
  "CD137L | PE.Cy7.A",
  "CD226 | ALEXA.FLUOR.647.A",
  "CD19 | Spark.NIR.685.A",
  "CD278 | Spark.Red.718.A",
  "CD134 | APC.FIRE.750.A",
  "CD274 | APC.FIRE.810.A"
)

# ========================= PLOT OPTIONS =====================================

# Plot types to generate: "stacked", "grouped", "proportional"
plot_types <- c("stacked", "proportional")
