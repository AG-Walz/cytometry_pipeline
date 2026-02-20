# =============================================================================
# config/pipeline_E_config.R — Pipeline E: Preprocessed FCS -> UMAP/Clustering
# =============================================================================
# Flow: Load FCS (skip QC) -> Homogenize -> Downsample -> Leiden -> UMAP
#
# Usage:
#   cd pipelines && Rscript pipeline_E_preprocessed_umap.R --config ../config/pipeline_E_config.R
# =============================================================================

# Compute project root relative to this config file
.config_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) {
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) return(normalizePath(f))
    }
    stop("Cannot determine config file path")
  }
)
base_dir <- dirname(dirname(.config_path))

params <- list(

  # ========================= Input / Output ==================================
  fcs_folder = file.path(base_dir, "test_data", "fcs_samples"),
  output_dir = file.path(base_dir, "output", "pipeline_E_results"),

  # -- Markers --
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),

  # ========================= Homogenization ==================================
  markername_replacing = list(),

  # ========================= Downsampling ====================================
  downsample_method = "Consistent",
  downsample_n      = 100,

  # ========================= Method Selection =================================
  clustering_method = "leiden",         # "leiden", "flowsom", "phenograph", "dbscan"
  batch_method      = "none",           # "harmony", "cytonorm", "combat", "none"
  dimred_method     = "umap",           # "umap", "tsne", "diffusion_map"

  # ========================= Clustering (Leiden) =============================
  n_neighbors = 15,
  resolution  = 0.5,

  # ========================= Output ==========================================
  output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
  plot_width     = 10,
  plot_height    = 8,

  # ========================= Reproducibility =================================
  seed = 42
)
