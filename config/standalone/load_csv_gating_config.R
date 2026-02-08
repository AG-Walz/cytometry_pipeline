# =============================================================================
# Standalone config: src/loading/load_csv_gating.R
# =============================================================================
# Loads CSV files containing binary gating columns alongside expression data.
# Auto-classifies columns (index, scatter, marker, gating, dim-reduction,
# autofluorescence) and optionally matches with FCS files for validation.
#
# Usage:
#   source("config/standalone/load_csv_gating_config.R")
#   source("src/loading/load_csv_gating.R")
#   result <- load_csv_gating(
#     input_fcs_path  = input_fcs_path,
#     input_csv_path  = input_csv_path,
#     output_path     = output_path,
#     gating_pattern  = gating_pattern,
#     csv_label_dict  = csv_label_dict,
#     processing_mode = processing_mode,
#     chunk_size      = chunk_size
#   )
# =============================================================================

# ========================= EDIT THESE ========================================

input_fcs_path <- "/path/to/fcs/files"    # FCS directory (for row count validation)
input_csv_path <- "/path/to/csv/files"    # CSV directory with gating columns

# ========================= OPTIONAL ==========================================

output_path     <- NULL   # NULL = input_fcs_path/Expression_output
gating_pattern  <- "\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$"
csv_label_dict  <- NULL   # Named vector for renaming gating columns
gating_tree     <- NULL   # Hierarchical gating structure
processing_mode <- "automatic"   # "automatic" or "manual"
output_format   <- "wide"        # "wide" or "long"
chunk_size      <- 100000        # Rows per chunk for memory efficiency
