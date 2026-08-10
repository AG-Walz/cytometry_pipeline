#!/usr/bin/env Rscript
# =============================================================================
# harmonization/panel_comparison.R — Compare Panels Across Samples
# =============================================================================
# Merged from: R_scripts/03_normalization/Panels.R + 00_utilities/Compare.R
# =============================================================================

# Sourced either standalone (from src/<subproject>/) or by a CLI wrapper that
# has already set `shared_dir` and loaded the shared utilities. Only the
# standalone case has to go looking for them, and the getwd()-relative guess is
# what makes that case fragile — so it is not made when it is not needed.
if (!exists("explore_data", mode = "function")) {
  .panel_shared <- Find(dir.exists, c(
    if (exists("shared_dir")) shared_dir,
    if (exists("cyto_shared_dir", mode = "function")) cyto_shared_dir(),
    file.path(getwd(), "src", "R_shared"),
    file.path(dirname(getwd()), "src", "R_shared")))
  if (!is.null(.panel_shared)) {
    source(file.path(.panel_shared, "io_utils.R"))
    source(file.path(.panel_shared, "flowframe_utils.R"))
  }
}

#' Compare marker panels given as plain character vectors.
#'
#' The pure core of compare_panels_detailed(): everything that does not need a
#' flowFrame. A "panel" here is one named source's marker names — a sample, a
#' FlowJo workspace, a batch — so the same comparison serves any of them.
#' Returns data, prints nothing; the caller decides how to say it.
#'
#' @param panels Named list of character vectors (names = source labels).
#' @return List with common, all_markers, not_common, presence (logical matrix,
#'   markers x sources), and per_source counts.
compare_marker_sets <- function(panels) {
  panels <- lapply(panels, function(x) unique(as.character(x[!is.na(x)])))
  sources <- names(panels)
  if (is.null(sources)) sources <- as.character(seq_along(panels))

  all_unique <- unique(unlist(panels, use.names = FALSE))
  common <- if (length(panels)) Reduce(intersect, panels) else character(0)
  # First source's order, so the result is deterministic and reads the way the
  # first panel does rather than alphabetically.
  if (length(panels)) common <- panels[[1]][panels[[1]] %in% common]

  presence <- matrix(FALSE, nrow = length(all_unique), ncol = length(panels),
                     dimnames = list(all_unique, sources))
  for (i in seq_along(panels)) presence[panels[[i]], i] <- TRUE

  list(
    common      = common,
    all_markers = all_unique,
    not_common  = setdiff(all_unique, common),
    presence    = presence,
    per_source  = vapply(panels, length, integer(1))
  )
}

#' Compare panels across a list of flowFrames and report differences
#'
#' @param flow_frames Named list of flowFrame objects
#' @return List with common_markers, unique_markers, panel_groups
compare_panels_detailed <- function(flow_frames) {
  cat("Comparing panels across", length(flow_frames), "samples...\n")

  # Extract markers from each sample, then hand the comparison itself to the
  # pure core so this and the workspace path cannot drift apart.
  all_markers <- lapply(flow_frames, function(ff) {
    data <- explore_data(ff)
    data$markers[data$mcols]
  })
  cmp <- compare_marker_sets(all_markers)
  sample_names <- names(flow_frames)

  cat("Common markers:", length(cmp$common), "\n")
  cat("Total unique markers:", length(cmp$all_markers), "\n")

  # Report differences
  if (length(cmp$not_common) > 0) {
    cat("\nMarkers NOT present in all samples:\n")
    for (m in cmp$not_common) {
      present_in <- sample_names[cmp$presence[m, ]]
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
    common_markers = cmp$common,
    all_markers = cmp$all_markers,
    not_common = cmp$not_common,
    presence_matrix = cmp$presence,
    panel_groups = panel_groups,
    per_sample_markers = all_markers
  )
}
