# =============================================================================
# Standalone config: metadata_qc/check_acquisition.R
# =============================================================================
# Extracts acquisition dates and metadata from FCS file headers. Detects
# likely duplicate files by matching date + file size combinations.
#
# Usage:
#   source("config/standalone/check_acquisition_config.R")
#   source("metadata_qc/check_acquisition.R")
#   result <- check_acquisition_dates(
#     parent_folder = parent_folder,
#     output_dir    = output_dir
#   )
#
#   # Or inspect a single FCS file:
#   get_fcs_metadata("/path/to/file.fcs")
#   explore_fcs_keywords("/path/to/file.fcs")
# =============================================================================

# ========================= EDIT THESE ========================================

parent_folder <- "/path/to/fcs/folder"         # Directory (recursive) with FCS files
output_dir    <- "/path/to/output/acquisition"  # CSV output with dates and metadata
