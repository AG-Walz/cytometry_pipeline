#!/usr/bin/env Rscript
# =============================================================================
# analysis/marker_statistics.R — Per-Marker Statistical Tests
# =============================================================================
# Merged from: R_scripts/05_analysis/Stats.R + Wilcox.R
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
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Run per-marker statistical analysis on expression CSVs
#'
#' @param input_folder Directory with expression CSVs
#' @param output_folder Output directory
#' @param filename_regex Regex to parse sample/population from filenames
#' @param markers Character vector of marker columns (NULL = interactive)
#' @param test_type "lmer" (mixed-effects), "wilcox", "kruskal", or "both" (lmer + wilcox)
#' @param conditions_csv Optional path to CSV with sample_id and condition columns
#'   for group comparisons (required for "wilcox", "kruskal", and "both")
#' @param n_workers Number of parallel workers
#' @return Data frame with per-marker per-population statistics
run_marker_statistics <- function(input_folder,
                                   output_folder,
                                   filename_regex = "^(.+)_(.+)\\.csv$",
                                   markers = NULL,
                                   test_type = "lmer",
                                   conditions_csv = NULL,
                                   n_workers = NULL) {
  require(data.table)
  require(lme4)
  require(lmerTest)
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required for run_marker_statistics(). ",
         "Install with: install.packages('stringr')")
  }

  make_outdir(output_folder)

  if (is.null(n_workers)) {
    n_workers <- max(1, parallel::detectCores() - 1)
  }

  # Setup parallel
  future::plan(future::multisession, workers = n_workers)

  # Load conditions mapping if provided
  conditions_map <- NULL
  if (!is.null(conditions_csv)) {
    conditions_map <- fread(conditions_csv)
    if (!all(c("sample_id", "condition") %in% names(conditions_map))) {
      stop("conditions_csv must have 'sample_id' and 'condition' columns.")
    }
    cat("Loaded condition mapping:", nlevels(factor(conditions_map$condition)),
        "groups\n")
  } else if (test_type %in% c("wilcox", "kruskal", "both")) {
    stop("conditions_csv is required for test_type '", test_type, "'.")
  }

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

    # Merge condition labels if available
    if (!is.null(conditions_map)) {
      pop_data <- merge(pop_data, conditions_map, by = "sample_id", all.x = TRUE)
      n_conditions <- length(unique(pop_data$condition[!is.na(pop_data$condition)]))
      cat("  Conditions:", paste(unique(na.omit(pop_data$condition)), collapse = ", "), "\n")
    }

    # Per-marker analysis
    for (marker in markers) {
      if (!marker %in% names(pop_data)) next

      result_key <- paste0(pop, "_", marker)

      # --- lmer branch ---
      if (test_type %in% c("lmer", "both")) {
        marker_stats <- tryCatch({
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
            pvalue = NA, statistic = NA,
            n_cells = nrow(pop_data),
            n_samples = length(unique(pop_data$sample_id)),
            method = "lmer"
          )
        }, error = function(e) {
          data.frame(
            population = pop, marker = marker,
            mean = mean(pop_data[[marker]], na.rm = TRUE),
            between_sample_var = NA, within_sample_var = NA, icc = NA,
            pvalue = NA, statistic = NA,
            n_cells = nrow(pop_data),
            n_samples = length(unique(pop_data$sample_id)),
            method = "simple_mean"
          )
        })
        all_results[[paste0(result_key, "_lmer")]] <- marker_stats
      }

      # --- wilcox branch ---
      if (test_type %in% c("wilcox", "both")) {
        marker_stats <- tryCatch({
          # Pseudobulk: median per sample
          sample_medians <- pop_data[, .(median_val = median(get(marker), na.rm = TRUE)),
                                     by = .(sample_id, condition)]
          conditions <- unique(sample_medians$condition)
          if (length(conditions) != 2) {
            stop("Wilcoxon requires exactly 2 conditions, found ", length(conditions))
          }

          g1 <- sample_medians[condition == conditions[1]]$median_val
          g2 <- sample_medians[condition == conditions[2]]$median_val
          wt <- wilcox.test(g1, g2)

          data.frame(
            population = pop, marker = marker,
            mean = mean(pop_data[[marker]], na.rm = TRUE),
            between_sample_var = NA, within_sample_var = NA, icc = NA,
            pvalue = wt$p.value, statistic = wt$statistic,
            n_cells = nrow(pop_data),
            n_samples = length(unique(pop_data$sample_id)),
            method = "wilcox"
          )
        }, error = function(e) {
          warning("Wilcoxon test failed for ", marker, " in ", pop, ": ", e$message)
          data.frame(
            population = pop, marker = marker,
            mean = mean(pop_data[[marker]], na.rm = TRUE),
            between_sample_var = NA, within_sample_var = NA, icc = NA,
            pvalue = NA, statistic = NA,
            n_cells = nrow(pop_data),
            n_samples = length(unique(pop_data$sample_id)),
            method = "wilcox_failed"
          )
        })
        all_results[[paste0(result_key, "_wilcox")]] <- marker_stats
      }

      # --- kruskal branch ---
      if (test_type == "kruskal") {
        marker_stats <- tryCatch({
          # Pseudobulk: median per sample
          sample_medians <- pop_data[, .(median_val = median(get(marker), na.rm = TRUE)),
                                     by = .(sample_id, condition)]
          kt <- kruskal.test(median_val ~ condition, data = sample_medians)

          data.frame(
            population = pop, marker = marker,
            mean = mean(pop_data[[marker]], na.rm = TRUE),
            between_sample_var = NA, within_sample_var = NA, icc = NA,
            pvalue = kt$p.value, statistic = kt$statistic,
            n_cells = nrow(pop_data),
            n_samples = length(unique(pop_data$sample_id)),
            method = "kruskal"
          )
        }, error = function(e) {
          warning("Kruskal-Wallis test failed for ", marker, " in ", pop, ": ", e$message)
          data.frame(
            population = pop, marker = marker,
            mean = mean(pop_data[[marker]], na.rm = TRUE),
            between_sample_var = NA, within_sample_var = NA, icc = NA,
            pvalue = NA, statistic = NA,
            n_cells = nrow(pop_data),
            n_samples = length(unique(pop_data$sample_id)),
            method = "kruskal_failed"
          )
        })
        all_results[[paste0(result_key, "_kruskal")]] <- marker_stats
      }
    }
  }

  results_df <- do.call(rbind, all_results)
  fwrite(results_df, file.path(output_folder, "marker_statistics.csv"), sep = ";")
  cat("\nResults saved to:", output_folder, "\n")

  results_df
}
