# =============================================================================
# Standalone config: src/preprocessing/compensation.R
# =============================================================================
# Applies fluorescence compensation to FCS files using either an embedded
# spillover matrix or an externally supplied compensation matrix CSV.
#
# Usage:
#   source("config/standalone/compensation_config.R")
#   source("src/preprocessing/compensation.R")
#   result <- run_compensation(
#     fcs_folder     = fcs_folder,
#     output_dir     = output_dir,
#     comp_matrix_csv = comp_matrix_csv,
#     auto_extract   = auto_extract
#   )
#
# Or from CLI:
#   Rscript src/preprocessing/compensation.R <fcs_folder> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

fcs_folder <- "/path/to/fcs/files"       # Directory containing raw FCS files
output_dir <- "/path/to/output/compensated"

# ========================= OPTIONAL ==========================================

# Path to external compensation matrix CSV (NULL = use embedded spillover)
comp_matrix_csv <- NULL

# Auto-extract spillover matrix from FCS $SPILL or $SPILLOVER keywords
auto_extract <- TRUE
