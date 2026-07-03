#!/usr/bin/env Rscript
# =============================================================================
# metadata_qc/overlap_statistics.R — Inter-rater Agreement (Fleiss' Kappa)
# =============================================================================
# Refactored from: R_scripts/05_analysis/Overlap_statistics.R
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

#' Compute Fleiss' kappa for label list subsets
#'
#' Given a data frame where each column represents a rater's labels,
#' computes Fleiss' kappa for all columns and optionally for pairwise subsets.
#'
#' @param label_lists Data frame of label assignments (cols = raters)
#' @param output_dir Output directory for log file
#' @param subset_specs Optional named list of column index vectors for subsets.
#'        If NULL, only the full dataset is analyzed.
#' @param subset_labels Optional character vector of labels for each subset
#' @param row_range Optional integer vector of row indices to subset
#' @return List of irr::kappam.fleiss results
compute_fleiss_kappa <- function(label_lists, output_dir,
                                  subset_specs = NULL,
                                  subset_labels = NULL,
                                  row_range = NULL) {
  require(irr)
  require(data.table)

  # Validate label_lists input
  if (!is.data.frame(label_lists)) {
    stop("label_lists must be a data frame with one column per rater")
  }
  if (ncol(label_lists) < 2) {
    stop("At least 2 raters (columns) required for Fleiss' kappa")
  }

  make_outdir(output_dir)

  if (!is.null(row_range)) {
    label_lists <- label_lists[row_range, , drop = FALSE]
  }

  if (is.null(subset_specs)) {
    subset_specs <- list(label_lists)
    subset_labels <- c("Full data")
  } else {
    subset_specs <- lapply(subset_specs, function(cols) label_lists[, cols, drop = FALSE])
    if (is.null(subset_labels)) {
      subset_labels <- paste0("Subset_", seq_along(subset_specs))
    }
  }

  log_file <- file.path(output_dir, "fleiss_kappa_log.txt")
  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  results <- list()
  for (i in seq_along(subset_specs)) {
    step_start <- Sys.time()
    cat("====================================\n", file = log_con, append = TRUE)
    cat(paste0("Analysis ", i, ": ", subset_labels[i], "\n"), file = log_con, append = TRUE)
    cat("Analysing:", subset_labels[i], "\n")

    results[[i]] <- irr::kappam.fleiss(ratings = subset_specs[[i]],
                                       exact = FALSE, detail = TRUE)
    capture.output(results[[i]], file = log_con, append = TRUE)

    elapsed <- as.numeric(difftime(Sys.time(), step_start, units = "mins"))
    cat(sprintf("Elapsed: %.2f minutes\n", elapsed), file = log_con, append = TRUE)
  }

  names(results) <- subset_labels
  cat("Fleiss' kappa log saved to:", log_file, "\n")
  results
}
