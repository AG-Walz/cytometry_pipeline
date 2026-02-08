# =============================================================================
# Standalone config: src/metadata_qc/check_metadata.R
# =============================================================================
# Validates and visualizes study metadata. Detects missing values, temporal
# patterns, and metadata inconsistencies. Generates barplot summaries.
#
# Usage:
#   source("config/standalone/check_metadata_config.R")
#   source("src/metadata_qc/check_metadata.R")
#   check_metadata(
#     metadata_file = metadata_file,
#     output_dir    = output_dir,
#     sep           = sep,
#     skip_columns  = skip_columns,
#     placeholders  = placeholders
#   )
#
#   # Optional: temporal heatmap
#   plot_temporal_heatmap(
#     metadata_file = metadata_file,
#     output_dir    = output_dir,
#     date_column   = date_column
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

metadata_file <- "/path/to/metadata.csv"
output_dir    <- "/path/to/output/metadata_qc"

# ========================= OPTIONAL ==========================================

sep          <- ","                                      # CSV separator
skip_columns <- NULL                                     # Columns to exclude from checks
only_cols    <- NULL                                     # Only check these columns (NULL = all)
sort_column  <- NULL                                     # Sort metadata by this column
sort_desc    <- FALSE                                    # Sort descending
n_splits     <- 3                                        # Number of temporal splits for heatmap
placeholders <- c("?", "na", "n/a", "null", "-")        # Values treated as missing
date_column  <- "Datum"                                  # Column name for acquisition date
