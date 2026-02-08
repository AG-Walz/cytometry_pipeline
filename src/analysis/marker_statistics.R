#!/usr/bin/env Rscript
# =============================================================================
# analysis/marker_statistics.R — Per-Marker Statistical Tests
# =============================================================================
# Merged from: R_scripts/05_analysis/Stats.R + Wilcox.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Run per-marker statistical analysis on expression CSVs
#'
#' @param input_folder Directory with expression CSVs
#' @param output_folder Output directory
#' @param filename_regex Regex to parse sample/population from filenames
#' @param markers Character vector of marker columns (NULL = interactive)
#' @param test_type "lmer" (mixed-effects), "wilcox", or "both"
#' @param n_workers Number of parallel workers
#' @return Data frame with per-marker per-population statistics
run_marker_statistics <- function(input_folder,
                                   output_folder,
                                   filename_regex = "^(.+)_(.+)\\.csv$",
                                   markers = NULL,
                                   test_type = "lmer",
                                   n_workers = NULL) {
  require(data.table)
  require(lme4)
  require(lmerTest)

  make_outdir(output_folder)

  if (is.null(n_workers)) {
    n_workers <- max(1, parallel::detectCores() - 1)
  }

  # Setup parallel
  future::plan(future::multisession, workers = n_workers)

  # Discover files
  files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
  stopifnot(length(files) > 0)
  cat("Found", length(files), "CSV files\n")

  # Parse metadata from filenames
  file_meta <- data.table(
    path = files,
    filename = basename(files)
  )
  parsed <- stringr::str_match(file_meta$filename, filename_regex)
  file_meta$sample_id <- parsed[, 2]
  file_meta$population_id <- parsed[, 3]

  # Column selection
  sample_cols <- names(fread(files[1], nrows = 0))
  if (is.null(markers)) {
    markers <- select_columns_gui(sample_cols, "Select Markers for Statistics")
    if (is.null(markers) || length(markers) == 0) markers <- sample_cols
  }

  cat("Analyzing", length(markers), "markers across",
      length(unique(file_meta$population_id)), "populations\n")

  # Analyze per population
  populations <- unique(file_meta$population_id)
  all_results <- list()

  for (pop in populations) {
    cat("\n--- Population:", pop, "---\n")
    pop_files <- file_meta[population_id == pop]

    # Load data for this population
    pop_data_list <- lapply(seq_len(nrow(pop_files)), function(i) {
      dt <- fread(pop_files$path[i], select = markers)
      dt$sample_id <- pop_files$sample_id[i]
      dt
    })
    pop_data <- rbindlist(pop_data_list)
    cat("  Loaded", nrow(pop_data), "cells from", nrow(pop_files), "files\n")

    # Per-marker analysis
    for (marker in markers) {
      if (!marker %in% names(pop_data)) next

      marker_stats <- tryCatch({
        if (test_type %in% c("lmer", "both")) {
          # Mixed-effects model
          model <- lmerTest::lmer(
            as.formula(paste0("`", marker, "` ~ 1 + (1|sample_id)")),
            data = pop_data
          )
          fixef_val <- lme4::fixef(model)
          var_comp <- as.data.frame(lme4::VarCorr(model))

          data.frame(
            population = pop, marker = marker,
            mean = fixef_val[1],
            between_sample_var = var_comp$vcov[1],
            within_sample_var = var_comp$vcov[2],
            icc = var_comp$vcov[1] / sum(var_comp$vcov),
            n_cells = nrow(pop_data),
            n_samples = length(unique(pop_data$sample_id)),
            method = "lmer"
          )
        }
      }, error = function(e) {
        # Fallback to simple summary
        data.frame(
          population = pop, marker = marker,
          mean = mean(pop_data[[marker]], na.rm = TRUE),
          between_sample_var = NA, within_sample_var = NA, icc = NA,
          n_cells = nrow(pop_data),
          n_samples = length(unique(pop_data$sample_id)),
          method = "simple_mean"
        )
      })
      all_results[[paste0(pop, "_", marker)]] <- marker_stats
    }
  }

  results_df <- do.call(rbind, all_results)
  fwrite(results_df, file.path(output_folder, "marker_statistics.csv"), sep = ";")
  cat("\nResults saved to:", output_folder, "\n")

  results_df
}
