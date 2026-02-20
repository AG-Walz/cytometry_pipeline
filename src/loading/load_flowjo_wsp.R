#!/usr/bin/env Rscript
# =============================================================================
# loading/load_flowjo_wsp.R — Load FlowJo Workspace Files
# =============================================================================
# Refactored from: R_scripts/01_loading/Load_wsp.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Load a FlowJo workspace (.wsp) file
#'
#' @param wsp_path Path to the .wsp file
#' @param fcs_dir Directory containing the FCS files referenced by the workspace
#' @param name Label for this workspace
#' @return List with gating_set, sample_names, population_names
load_flowjo_wsp <- function(wsp_path, fcs_dir = dirname(wsp_path), name = NULL) {
  require(flowWorkspace)
  require(CytoML)

  stopifnot(file.exists(wsp_path))
  if (is.null(name)) name <- tools::file_path_sans_ext(basename(wsp_path))

  cat("Loading workspace:", basename(wsp_path), "\n")
  cat("FCS directory:", fcs_dir, "\n")

  # Parse workspace
  ws <- CytoML::open_flowjo_xml(wsp_path)
  cat("Workspace groups:", paste(flowWorkspace::fj_ws_get_sample_groups(ws)$groupName,
                                  collapse = ", "), "\n")

  # Import as GatingSet
  gs <- CytoML::flowjo_to_gatingset(ws, name = 1, path = fcs_dir)

  sample_names <- flowWorkspace::sampleNames(gs)
  cat("Loaded", length(sample_names), "samples\n")

  # Get population tree
  pop_names <- flowWorkspace::gs_get_pop_paths(gs, path = "full")
  cat("Populations:", length(pop_names), "\n")
  for (p in pop_names) cat("  ", p, "\n")

  list(
    gating_set = gs,
    sample_names = sample_names,
    population_names = pop_names,
    workspace_name = name
  )
}

#' Get population statistics from a GatingSet
#'
#' @param gating_set A GatingSet object
#' @param type NULL for counts, "freq" for frequencies
#' @return Data frame with Population, Sample, Value columns
get_population_stats <- function(gating_set, type = NULL) {
  sample_names <- flowWorkspace::sampleNames(gating_set)
  stats_list <- list()
  for (sample in sample_names) {
    if (is.null(type)) {
      stats <- tibble::tibble(flowWorkspace::gh_pop_get_stats(gating_set[[sample]]))
      stats_list[[sample]] <- data.frame(
        Population = stats$pop, Sample = sample, Value = stats$count
      )
    } else {
      stats <- tibble::tibble(flowWorkspace::gh_pop_get_stats(gating_set[[sample]], type = type))
      stats_list[[sample]] <- data.frame(
        Population = stats$pop, Sample = sample, Value = stats$percent
      )
    }
  }
  do.call(rbind, stats_list)
}

#' Extract population string assignments for each cell
#'
#' Output to: load_multiple_workspaces() (this file),
#'            compare_gatings() (metadata_qc/compare_gatings.R)
#'
#' @param gating_set A GatingSet object
#' @param sample_names Character vector of sample names
#' @return Named list of character vectors (one per sample)
get_population_strings_list <- function(gating_set, sample_names) {
  pop_strings <- list()
  for (i in seq_along(sample_names)) {
    gh <- gating_set[[i]]
    pop_names <- flowWorkspace::gs_get_pop_paths(gh, path = "full")
    cytoframe <- flowWorkspace::gh_pop_get_data(gh)
    n_events <- nrow(flowCore::exprs(cytoframe))
    strings <- rep(NA, n_events)
    for (pop_name in pop_names) {
      indices <- flowWorkspace::gh_pop_get_indices(gh, pop_name)
      strings[indices] <- pop_name
    }
    pop_strings[[sample_names[i]]] <- strings
  }
  pop_strings
}

#' Summarize population counts across samples
#'
#' @param population_strings_list Output from get_population_strings_list
#' @return Data frame with Sample, Population, Count columns
summarize_population_counts <- function(population_strings_list) {
  pop_data <- do.call(rbind, lapply(names(population_strings_list), function(nm) {
    data.frame(Population = population_strings_list[[nm]], Sample = nm)
  }))
  pop_data %>%
    dplyr::group_by(Sample, Population) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop")
}

#' Save population strings list as RDS
#'
#' @param population_strings_list Named list of population strings
#' @param out_dir Output directory
save_population_strings <- function(population_strings_list, out_dir) {
  out_file <- file.path(out_dir, "population_strings_list.rds")
  saveRDS(population_strings_list, file = out_file)
  cat("Population strings saved to:", out_file, "\n")
}

#' Load population strings from multiple FlowJo workspaces
#'
#' Iterates over workspace paths, loads each with load_flowjo_wsp(), and
#' extracts population string assignments using get_population_strings_list().
#'
#' Input from: .wsp files on disk
#' Output to:  compare_gatings() (metadata_qc/compare_gatings.R),
#'             compute_fleiss_kappa() (metadata_qc/overlap_statistics.R)
#'
#' @param wsp_paths Character vector of paths to .wsp files
#' @param names Character vector of labels for each workspace (same length as wsp_paths)
#' @param fcs_dirs Character vector of FCS directories (default: dirname of each wsp_path).
#'   If length 1, recycled for all workspaces.
#' @return Named list where each element is the population_strings_list for one workspace
load_multiple_workspaces <- function(wsp_paths, names, fcs_dirs = NULL) {
  stopifnot(length(wsp_paths) == length(names))
  if (is.null(fcs_dirs)) fcs_dirs <- dirname(wsp_paths)
  if (length(fcs_dirs) == 1) fcs_dirs <- rep(fcs_dirs, length(wsp_paths))

  pop_strings_list <- list()
  for (i in seq_along(wsp_paths)) {
    cat("Loading workspace", i, "of", length(wsp_paths), ":", wsp_paths[i], "\n")
    wsp_data <- load_flowjo_wsp(wsp_path = wsp_paths[i], fcs_dir = fcs_dirs[i])
    pop_strings_list[[names[i]]] <- get_population_strings_list(
      wsp_data$gating_set, wsp_data$sample_names
    )
  }

  cat("Loaded", length(pop_strings_list), "workspaces:",
      paste(names, collapse = ", "), "\n")
  pop_strings_list
}
