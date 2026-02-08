# =============================================================================
# Standalone config: src/loading/load_fcs_raw.R
# =============================================================================
# Loads raw FCS files from a directory into flowFrame objects. Optionally
# saves expression matrices as CSV and truncates extreme values.
#
# Usage:
#   source("config/standalone/load_fcs_config.R")
#   source("src/loading/load_fcs_raw.R")
#   flow_data <- load_fcs_raw(
#     input_dir    = input_dir,
#     sample_names = sample_names,
#     truncation   = truncation,
#     extract_csv  = extract_csv,
#     output_dir   = output_dir
#   )
#
# Or from CLI:
#   Rscript src/loading/load_fcs_raw.R /path/to/fcs [/path/to/output]
# =============================================================================

# ========================= EDIT THESE ========================================

input_dir <- "/path/to/fcs/files"

# ========================= OPTIONAL ==========================================

sample_names <- NULL     # NULL = auto-detect from filenames
truncation   <- FALSE    # Truncate extreme fluorescence values
min_limit    <- NULL     # Minimum cell count to keep a sample (NULL = keep all)
extract_csv  <- FALSE    # Also save expression matrices as CSV
output_dir   <- NULL     # Where to save CSVs (required if extract_csv = TRUE)
