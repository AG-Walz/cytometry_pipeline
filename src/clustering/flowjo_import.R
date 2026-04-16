#!/usr/bin/env Rscript
# =============================================================================
# clustering/flowjo_import.R — Import Gates from FlowJo Workspaces
# =============================================================================
# Merged from: R_scripts/04_gating/FlowJo_wsp.R + flowWorkspace.R
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

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

  # Parse and import
  ws <- CytoML::open_flowjo_xml(wsp_path)
  gs <- CytoML::flowjo_to_gatingset(ws, name = 1, path = fcs_dir)

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
