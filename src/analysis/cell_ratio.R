#!/usr/bin/env Rscript
# =============================================================================
# analysis/cell_ratio.R — Cell Type Ratio Calculations
# =============================================================================
# Refactored from: R_scripts/05_analysis/Cell_ratio.R
# =============================================================================

# ---- cytometry src bootstrap (Nextflow-safe); see src/R_shared/paths.R ----
if (!exists(".cyto_paths_loaded", inherits = TRUE)) {
  local({
    root <- Sys.getenv("CYTOMETRY_SRC", "")
    if (!nzchar(root) || !dir.exists(file.path(root, "R_shared"))) {
      of <- NULL
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        f <- frames[[i]]$ofile
        if (!is.null(f) && nzchar(f)) { of <- f; break }
      }
      d <- if (!is.null(of)) dirname(normalizePath(of)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))

#' Calculate cell type ratios from abundance data
#'
#' @param abundance_csv Path to CSV with cluster counts per sample
#' @param numerator_clusters Character vector of cluster names for numerator
#' @param denominator_clusters Character vector of cluster names for denominator
#' @param output_path Output CSV path
#' @return Data frame with ratios per sample
calculate_cell_ratios <- function(abundance_csv,
                                   numerator_clusters,
                                   denominator_clusters,
                                   output_path = NULL) {
  require(data.table)

  abundance <- fread(abundance_csv)
  id_col <- names(abundance)[1]

  # Calculate numerator and denominator sums
  num_cols <- intersect(numerator_clusters, names(abundance))
  den_cols <- intersect(denominator_clusters, names(abundance))

  if (length(num_cols) == 0) stop("No numerator clusters found")
  if (length(den_cols) == 0) stop("No denominator clusters found")

  result <- data.frame(
    sample = abundance[[id_col]],
    numerator_sum = rowSums(abundance[, ..num_cols]),
    denominator_sum = rowSums(abundance[, ..den_cols])
  )
  result$ratio <- result$numerator_sum / pmax(result$denominator_sum, 1)
  result$log2_ratio <- log2(result$ratio + 1e-6)

  cat("Ratio:", paste(num_cols, collapse = "+"), "/",
      paste(den_cols, collapse = "+"), "\n")
  cat("Mean ratio:", round(mean(result$ratio, na.rm = TRUE), 4), "\n")

  if (!is.null(output_path)) {
    fwrite(result, output_path, sep = ";")
    cat("Saved to:", output_path, "\n")
  }

  result
}
