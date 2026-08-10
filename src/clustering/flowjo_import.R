#!/usr/bin/env Rscript
# =============================================================================
# clustering/flowjo_import.R — Import Gates from FlowJo Workspaces
# =============================================================================
# Merged from: R_scripts/04_gating/FlowJo_wsp.R + flowWorkspace.R
# =============================================================================

# Respect a shared_dir the caller already set; only guess from the working
# directory when there is none (see the same note in loading/load_flowjo_wsp.R).
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
      # Rscript src/foo/bar.R carries no source frame; recover from --file=.
      if (is.null(of)) {
        a <- commandArgs(FALSE); m <- grep("^--file=", a)
        if (length(m)) of <- sub("^--file=", "", a[m[[1]]])
      }
      d <- if (!is.null(of)) dirname(normalizePath(of, mustWork = FALSE)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

# Guarded, NOT unconditional: nextflow/bin/*_cli.R and scripts/*.R set shared_dir
# themselves (the --src_dir convention), and clobbering it here would send every
# source() below to whichever tree this file happens to live in rather than the
# one the caller chose.
if (!exists("shared_dir")) shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

# All CytoML workspace calls live in src/loading/flowjo_workspace.R.
if (!exists("wsp_open", mode = "function")) {
  source(file.path(dirname(shared_dir), "loading", "flowjo_workspace.R"))
}

#' Import cluster/gate assignments from a FlowJo workspace
#'
#' @param wsp_path Path to .wsp file
#' @param fcs_dir Directory with FCS files
#' @param target_populations Character vector of population paths to extract
#'   (NULL = extract all)
#' @return List with gating_set and assignments data.frame
import_flowjo_gates <- function(wsp_path, fcs_dir = dirname(wsp_path),
                                 target_populations = NULL) {
  require(flowWorkspace)
  require(CytoML)

  stopifnot(file.exists(wsp_path))

  # Parse and import (via loading/flowjo_workspace.R)
  ws <- wsp_open(wsp_path)
  gs <- wsp_to_gatingset(ws, group = 1, path = fcs_dir)

  sample_names <- flowWorkspace::sampleNames(gs)
  all_pops <- flowWorkspace::gs_get_pop_paths(gs, path = "full")

  if (!is.null(target_populations)) {
    pops_to_extract <- intersect(all_pops, target_populations)
  } else {
    pops_to_extract <- all_pops
  }

  cat("Extracting", length(pops_to_extract), "populations from",
      length(sample_names), "samples\n")

  # Build assignment data frame
  assignments_list <- list()
  for (sample in sample_names) {
    gh <- gs[[sample]]
    n_events <- nrow(flowCore::exprs(flowWorkspace::gh_pop_get_data(gh)))
    sample_df <- data.frame(cell_id = seq_len(n_events), sample = sample)

    for (pop in pops_to_extract) {
      indices <- flowWorkspace::gh_pop_get_indices(gh, pop)
      pop_clean <- clean_cluster_name(basename(pop))
      sample_df[[pop_clean]] <- as.integer(indices)
    }
    assignments_list[[sample]] <- sample_df
  }

  assignments <- do.call(rbind, assignments_list)
  cat("Total:", nrow(assignments), "cells with",
      length(pops_to_extract), "population columns\n")

  list(
    gating_set = gs,
    assignments = assignments,
    populations = pops_to_extract,
    sample_names = sample_names
  )
}
