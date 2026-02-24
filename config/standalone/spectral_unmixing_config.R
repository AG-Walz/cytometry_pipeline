# =============================================================================
# Standalone config: src/preprocessing/spectral_unmixing.R
# =============================================================================
# Performs spectral unmixing of a sample FCS file using single-stained control
# references and an unstained control for autofluorescence estimation.
#
# Usage:
#   source("config/standalone/spectral_unmixing_config.R")
#   source("src/preprocessing/spectral_unmixing.R")
#   result <- run_spectral_unmixing(
#     sample_fcs      = sample_fcs,
#     controls_folder = controls_folder,
#     unstained_fcs   = unstained_fcs,
#     output_dir      = output_dir,
#     channels        = channels
#   )
#
# Or from CLI:
#   Rscript src/preprocessing/spectral_unmixing.R <sample_fcs> <controls_folder> <unstained_fcs> <output_dir>
# =============================================================================

# ========================= EDIT THESE ========================================

sample_fcs      <- "/path/to/sample.fcs"           # FCS file to unmix
controls_folder <- "/path/to/single-stain/controls" # Directory with single-stained control FCS files
unstained_fcs   <- "/path/to/unstained.fcs"         # Unstained control for autofluorescence
output_dir      <- "/path/to/output/unmixed"

# ========================= OPTIONAL ==========================================

# Channels to unmix (NULL = auto-detect from controls)
channels <- NULL
