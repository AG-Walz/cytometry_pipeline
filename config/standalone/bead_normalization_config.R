# =============================================================================
# Standalone config: src/preprocessing/bead_normalization.R
# =============================================================================
# Normalizes mass/spectral cytometry data using bead standards to correct
# for signal drift over time. Optionally removes bead events from the output.
#
# Usage:
#   source("config/standalone/bead_normalization_config.R")
#   source("src/preprocessing/bead_normalization.R")
#   result <- run_bead_normalization(
#     fcs_file      = fcs_file,
#     output_dir    = output_dir,
#     bead_type     = bead_type,
#     bead_channels = bead_channels,
#     remove_beads  = remove_beads,
#     baseline      = baseline,
#     window_size   = window_size
#   )
#
# Or from CLI:
#   Rscript src/preprocessing/bead_normalization.R <fcs_file> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

fcs_file   <- "/path/to/sample.fcs"              # FCS file containing bead + cell events
output_dir <- "/path/to/output/normalized"

# ========================= OPTIONAL ==========================================

# Bead type: "four_element", "six_element", or "custom"
bead_type <- "four_element"

# Bead channels (NULL = auto-detect from bead_type)
bead_channels <- NULL

# Remove bead events from the output FCS
remove_beads <- TRUE

# Baseline for normalization: "median" or "mean"
baseline <- "median"

# Sliding window size (number of events) for smoothing bead signal
window_size <- 500
