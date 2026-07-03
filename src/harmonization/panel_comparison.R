#!/usr/bin/env Rscript
# =============================================================================
# harmonization/panel_comparison.R — Compare Panels Across Samples
# =============================================================================
# Merged from: R_scripts/03_normalization/Panels.R + 00_utilities/Compare.R
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
source(file.path(shared_dir, "flowframe_utils.R"))

#' Compare panels across a list of flowFrames and report differences
#'
#' @param flow_frames Named list of flowFrame objects
#' @return List with common_markers, unique_markers, panel_groups
compare_panels_detailed <- function(flow_frames) {
  cat("Comparing panels across", length(flow_frames), "samples...\n")

  # Extract markers from each sample
  all_markers <- lapply(flow_frames, function(ff) {
    data <- explore_data(ff)
    data$markers[data$mcols]
  })

  # Find common and unique markers
  common <- Reduce(intersect, all_markers)
  all_unique <- unique(unlist(all_markers))
  sample_names <- names(flow_frames)

  cat("Common markers:", length(common), "\n")
  cat("Total unique markers:", length(all_unique), "\n")

  # Build presence matrix
  presence <- matrix(FALSE, nrow = length(all_unique), ncol = length(flow_frames),
                     dimnames = list(all_unique, sample_names))
  for (i in seq_along(all_markers)) {
    presence[all_markers[[i]], i] <- TRUE
  }

  # Report differences
  not_common <- setdiff(all_unique, common)
  if (length(not_common) > 0) {
    cat("\nMarkers NOT present in all samples:\n")
    for (m in not_common) {
      present_in <- sample_names[presence[m, ]]
      cat(sprintf("  %s: present in %d/%d (%s)\n",
                  m, length(present_in), length(flow_frames),
                  paste(present_in, collapse = ", ")))
    }
  } else {
    cat("All samples have identical panels.\n")
  }

  # Group samples by panel
  panel_groups <- compare_panels(flow_frames)

  list(
    common_markers = common,
    all_markers = all_unique,
    not_common = not_common,
    presence_matrix = presence,
    panel_groups = panel_groups,
    per_sample_markers = all_markers
  )
}
