#!/usr/bin/env Rscript
# =============================================================================
# loading/load_fcs_raw.R — Load Raw FCS Files
# =============================================================================
# Refactored from: R_scripts/01_loading/Load_raw.R + Load_fcs.R
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

# Source shared utilities
shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "fcs_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Load raw FCS files from a directory
#'
#' @param input_dir Path to directory containing FCS files
#' @param sample_names Optional character vector of sample names (derived from filenames if NULL)
#' @param truncation Logical: truncate values to max range (default FALSE)
#' @param min_limit Optional lower truncation limit
#' @param extract_csv Logical: save raw expression matrices as CSV
#' @param output_dir Output directory for CSVs (required if extract_csv = TRUE)
#' @return List with flow_frames, sample_names, input_files
load_fcs_raw <- function(input_dir,
                         sample_names = NULL,
                         truncation = FALSE,
                         min_limit = NULL,
                         extract_csv = FALSE,
                         output_dir = NULL) {

  # Validate input
  stopifnot(dir.exists(input_dir))

  # Find FCS files
  input_fcs_files <- list.files(path = input_dir, pattern = "\\.fcs$",
                                full.names = TRUE)
  if (length(input_fcs_files) == 0) {
    stop("No FCS files found in: ", input_dir)
  }
  cat("Found", length(input_fcs_files), "FCS files in", input_dir, "\n")

  # Derive sample names if not provided
  if (is.null(sample_names)) {
    sample_names <- gsub("\\.fcs$", "", basename(input_fcs_files))
  }

  if (length(sample_names) != length(input_fcs_files)) {
    stop("Number of sample names (", length(sample_names),
         ") does not match number of FCS files (", length(input_fcs_files), ")")
  }

  # Create output dir if needed
  if (extract_csv) {
    if (is.null(output_dir)) output_dir <- file.path(input_dir, "raw_csv")
    dir_create(output_dir)
  }

  # Load each FCS file
  flow_frames <- list()
  for (i in seq_along(input_fcs_files)) {
    cat("Loading", basename(input_fcs_files[i]), "...")
    flow_frames[[i]] <- flowCore::read.FCS(
      input_fcs_files[i],
      transformation = FALSE,
      truncate_max_range = truncation,
      min.limit = min_limit
    )
    cat(" Sample", sample_names[i], "has", nrow(flow_frames[[i]]),
        "cells and", ncol(flow_frames[[i]]), "markers\n")

    if (extract_csv) {
      csv_path <- file.path(output_dir, paste0(sample_names[i], "_raw.csv"))
      write.csv(flowCore::exprs(flow_frames[[i]]), csv_path)
      cat("  Saved CSV:", csv_path, "\n")
    }
  }

  names(flow_frames) <- sample_names

  list(
    flow_frames = flow_frames,
    sample_names = sample_names,
    input_files = input_fcs_files
  )
}
