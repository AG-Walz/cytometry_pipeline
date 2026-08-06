#!/usr/bin/env Rscript
# =============================================================================
# metadata_qc/check_acquisition.R — FCS Acquisition Date Checker
# =============================================================================
# Refactored from: R_scripts/00_utilities/Check_AqDate.R
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

#' Extract FCS keyword from a single file
#'
#' @param fcs_path Path to FCS file
#' @param keyword_name FCS keyword to extract (e.g. "$DATE")
#' @return The keyword value, or "Unknown"/"Error"
get_fcs_keyword <- function(fcs_path, keyword_name) {
  require(flowCore)
  tryCatch({
    ff <- flowCore::read.FCS(fcs_path, transformation = FALSE,
                             which.lines = 1, truncate_max_range = FALSE)
    # flowCore::keyword(ff, name) returns a named list; unwrap to a scalar
    # character so downstream code can treat the result like a single value.
    val <- flowCore::keyword(ff, keyword_name)
    if (is.list(val)) val <- if (length(val) >= 1) val[[1]] else NULL
    if (is.null(val) || length(val) == 0 ||
        (length(val) == 1 && (is.na(val) || val == ""))) {
      "Unknown"
    } else {
      as.character(val)
    }
  }, error = function(e) {
    message("Error reading ", fcs_path, ": ", e$message)
    "Error"
  })
}

#' Explore available keywords in an FCS file
#'
#' @param fcs_path Path to FCS file
#' @return Named list of keywords (invisible)
explore_fcs_keywords <- function(fcs_path) {
  require(flowCore)
  tryCatch({
    ff <- flowCore::read.FCS(fcs_path, transformation = FALSE,
                             which.lines = 1, truncate_max_range = FALSE)
    keywords <- flowCore::keyword(ff)
    cat("File:", basename(fcs_path), "\n")
    if ("$DATATYPE" %in% names(keywords)) {
      cat("  $DATATYPE:", keywords[["$DATATYPE"]], "\n")
    }
    bit_kw <- names(keywords)[grepl("P[0-9]+B", names(keywords))]
    if (length(bit_kw) > 0) {
      cat("  Parameter bit depths:\n")
      for (kw in bit_kw) cat("    ", kw, ":", keywords[[kw]], "\n")
    }
    invisible(keywords)
  }, error = function(e) {
    message("Error exploring ", fcs_path, ": ", e$message)
    invisible(NULL)
  })
}

#' Get acquisition date and file size from an FCS file
#'
#' @param fcs_path Path to FCS file
#' @return Named character vector with acquisition_date and file_size
get_fcs_metadata <- function(fcs_path) {
  require(flowCore)
  tryCatch({
    ff <- flowCore::read.FCS(fcs_path, transformation = FALSE,
                             which.lines = 1, truncate_max_range = FALSE)
    # keyword() returns a named list, so extract the scalar value.
    acq_date <- flowCore::keyword(ff, "$DATE")
    if (is.list(acq_date)) acq_date <- if (length(acq_date) >= 1) acq_date[[1]] else NULL
    if (is.null(acq_date) || length(acq_date) == 0 ||
        (length(acq_date) == 1 && (is.na(acq_date) || acq_date == ""))) {
      acq_date <- "Unknown"
    }
    file_sz <- file.size(fcs_path)
    if (is.na(file_sz)) file_sz <- "Unknown" else file_sz <- as.character(file_sz)
    c(acquisition_date = as.character(acq_date), file_size = file_sz)
  }, error = function(e) {
    message("Error reading ", fcs_path, ": ", e$message)
    c(acquisition_date = "Error", file_size = "Error")
  })
}

#' Check acquisition dates for all FCS files in a folder
#'
#' Recursively scans for .fcs files, extracts acquisition dates and file sizes,
#' assigns unique IDs for files with matching date+size (likely duplicates).
#'
#' @param parent_folder Root folder to scan recursively
#' @param output_dir Output directory for CSV results
#' @return Data frame with file_relative_path, acquisition_date, file_size, unique_file_id
check_acquisition_dates <- function(parent_folder, output_dir) {
  require(data.table)
  make_outdir(output_dir)

  if (!dir.exists(parent_folder)) stop("Folder does not exist: ", parent_folder)

  all_fcs <- list.files(parent_folder, pattern = "\\.fcs$",
                        full.names = TRUE, recursive = TRUE)
  if (length(all_fcs) == 0) stop("No FCS files found in: ", parent_folder)
  cat("Found", length(all_fcs), "FCS files\n")

  # Relative paths
  parent_norm <- normalizePath(parent_folder, winslash = "/")
  if (!endsWith(parent_norm, "/")) {
    parent_norm <- paste0(parent_norm, "/")
  }
  rel_paths <- substring(normalizePath(all_fcs, winslash = "/"), nchar(parent_norm) + 1)

  # Extract metadata with per-file error handling
  cat("Extracting metadata...\n")
  results_list <- lapply(all_fcs, function(f) {
    tryCatch({
      cat("  Processing:", basename(f), "\n")
      get_fcs_metadata(f)
    }, error = function(e) {
      warning("Failed to read FCS file: ", f, " (", e$message, ")")
      c(acquisition_date = "Error", file_size = "Error")
    })
  })
  meta_list <- do.call(rbind, results_list)
  colnames(meta_list) <- c("acquisition_date", "file_size")

  result <- data.frame(
    file_relative_path = rel_paths,
    acquisition_date = as.character(meta_list[, "acquisition_date"]),
    file_size = as.character(meta_list[, "file_size"]),
    stringsAsFactors = FALSE
  )

  # Assign unique IDs for likely duplicates (same date + size)
  paste_vals <- paste(result$acquisition_date, result$file_size, sep = "_")
  paste_vals[is.na(paste_vals)] <- "Unknown_Unknown"
  result$unique_file_id <- as.integer(as.factor(paste_vals))
  result <- result[order(result$unique_file_id, result$acquisition_date), ]

  out_file <- file.path(output_dir, "acquisition_dates.csv")
  data.table::fwrite(result, file = out_file)
  cat("Acquisition date report saved to:", out_file, "\n")
  result
}
