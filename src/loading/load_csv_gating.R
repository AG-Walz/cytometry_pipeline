#!/usr/bin/env Rscript
# =============================================================================
# loading/load_csv_gating.R — Load CSV Files with Gating Data
# =============================================================================
# Refactored from: R_scripts/01_loading/Load_csv_gating.R
# Core merge/pivot logic retained; GUI, I/O, cluster utils moved to R_shared
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
source(file.path(shared_dir, "fcs_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))

#' Load CSV gating data and match with FCS files
#'
#' @param input_fcs_path Directory with FCS files
#' @param input_csv_path Directory with CSV files (defaults to input_fcs_path)
#' @param output_path Output directory for processed files
#' @param gating_pattern Regex to identify gating columns
#' @param csv_label_dict Named vector for label translation
#' @param gating_tree Named numeric vector of gating hierarchy depths
#' @param processing_mode "automatic" (GUI), "manual" (dict-based), or "expression_only"
#' @param output_format "wide" (per sample/gate CSVs) or "long" (single merged CSV)
#' @param chunk_size Rows per processing chunk
#' @return List with processed data paths and summary
load_csv_gating <- function(input_fcs_path,
                            input_csv_path = input_fcs_path,
                            output_path = file.path(input_fcs_path, "Expression_output"),
                            gating_pattern = "\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$",
                            csv_label_dict = NULL,
                            gating_tree = NULL,
                            processing_mode = "automatic",
                            output_format = "wide",
                            chunk_size = 100000) {
  require(data.table)
  make_outdir(output_path)

  # Discover files
  csv_files <- list.files(input_csv_path, pattern = "\\.csv$", full.names = TRUE)
  fcs_files <- list.files(input_fcs_path, pattern = "\\.fcs$", full.names = TRUE)

  if (length(csv_files) == 0) stop("No CSV files found in: ", input_csv_path)
  if (length(fcs_files) == 0) stop("No FCS files found in: ", input_fcs_path)

  cat("Found", length(csv_files), "CSV files and", length(fcs_files), "FCS files\n")

  # Create CSV-FCS mapping
  csv_to_fcs <- create_csv_fcs_mapping(csv_files, fcs_files)

  # Filter to files with gating columns
  valid_csvs <- filter_valid_csv_files(csv_files, pattern = gating_pattern)
  if (length(valid_csvs) == 0) {
    message("No CSV files with gating columns found. Processing as expression-only.")
    processing_mode <- "expression_only"
    valid_csvs <- csv_files
  }

  # Interactive column selection if automatic mode
  if (processing_mode == "automatic") {
    sample_csv <- valid_csvs[1]
    sep <- detect_csv_sep(sample_csv)
    header <- names(fread(sample_csv, sep = sep, nrows = 0))

    gating_cols <- grep(gating_pattern, header, value = TRUE)
    cat("\nDetected gating columns:\n")
    for (col in gating_cols) cat("  -", col, "\n")

    # Offer gating column selection
    selected_cols <- select_columns_gui(gating_cols, title = "Select Gating Columns")
    if (!is.null(selected_cols) && length(selected_cols) > 0) {
      gating_pattern <- paste0("^(", paste(gsub("([\\(\\)])", "\\\\\\1",
                                                 selected_cols), collapse = "|"), ")$")
    }

    # Classify all columns and let user select markers
    col_classes <- classify_csv_columns(header, gating_pattern)

    cat("\n--- Column classification ---\n")
    cat("  Index:             ", length(col_classes$index), "\n")
    cat("  Scatter/Time:      ", length(col_classes$scatter), "\n")
    cat("  Gating/Clustering: ", length(col_classes$gating), "\n")
    cat("  Dim. reduction:    ", length(col_classes$dimreduction), "\n")
    cat("  Autofluorescence:  ", length(col_classes$autofluorescence), "\n")
    cat("  Expression markers:", length(col_classes$markers), "\n")

    if (length(col_classes$scatter) > 0) {
      cat("\nExcluded scatter/time columns:\n")
      for (col in col_classes$scatter) cat("  -", col, "\n")
    }
    if (length(col_classes$autofluorescence) > 0) {
      cat("\nExcluded autofluorescence columns:\n")
      for (col in col_classes$autofluorescence) cat("  -", col, "\n")
    }
    if (length(col_classes$dimreduction) > 0) {
      cat("\nExcluded dimensionality reduction columns:\n")
      for (col in col_classes$dimreduction) cat("  -", col, "\n")
    }

    # Let user confirm/adjust marker selection
    selected_markers <- select_columns_gui(col_classes$markers,
                                            title = "Select Expression Markers")
    if (is.null(selected_markers) || length(selected_markers) == 0) {
      selected_markers <- col_classes$markers
    }
  }

  # In non-automatic modes, classify columns without GUI
  if (processing_mode != "automatic") {
    sample_csv <- valid_csvs[1]
    sep <- detect_csv_sep(sample_csv)
    sample_header <- names(fread(sample_csv, sep = sep, nrows = 0))
    col_classes <- classify_csv_columns(sample_header, gating_pattern)
    selected_markers <- col_classes$markers
  }

  # Process each CSV file
  results <- list()
  for (csv_file in valid_csvs) {
    result <- tryCatch({
      process_single_csv_gating(
        csv_file = csv_file,
        gating_pattern = gating_pattern,
        expr_cols = selected_markers,
        csv_label_dict = csv_label_dict,
        gating_tree = gating_tree,
        csv_to_fcs = csv_to_fcs,
        output_path = output_path,
        output_format = output_format,
        chunk_size = chunk_size
      )
    }, error = function(e) {
      warning("Error processing ", basename(csv_file), ": ", e$message)
      NULL
    })
    if (!is.null(result)) results[[basename(csv_file)]] <- result
    check_memory()
  }

  cat("\nProcessed", length(results), "/", length(valid_csvs), "files successfully\n")

  list(
    output_dir = output_path,
    results = results,
    n_processed = length(results),
    n_total = length(valid_csvs)
  )
}

#' Process a single CSV file with gating data
#' @keywords internal
process_single_csv_gating <- function(csv_file, gating_pattern, expr_cols = NULL,
                                       csv_label_dict, gating_tree, csv_to_fcs,
                                       output_path, output_format, chunk_size) {
  require(data.table)
  file_size_mb <- get_file_size_mb(csv_file)
  cat("Processing:", basename(csv_file), sprintf("(%.1f MB)...", file_size_mb))

  # Read CSV
  sep <- detect_csv_sep(csv_file)
  csv_data <- fread(csv_file, sep = sep, header = TRUE)

  # Identify gating and expression columns
  gating_cols <- grep(gating_pattern, names(csv_data), value = TRUE)
  if (is.null(expr_cols)) {
    expr_cols <- setdiff(names(csv_data), c(gating_cols, "index", "Index"))
  } else {
    # Keep only expr_cols that actually exist in this file
    expr_cols <- intersect(expr_cols, names(csv_data))
  }

  cat(sprintf(" %d gating cols, %d expression cols, %d rows",
              length(gating_cols), length(expr_cols), nrow(csv_data)))

  # Get matching FCS event count for validation
  fcs_file <- csv_to_fcs[[csv_file]]
  if (!is.null(fcs_file)) {
    fcs_events <- get_fcs_event_count(fcs_file)
    if (nrow(csv_data) != fcs_events) {
      cat(sprintf(" [WARNING: CSV %d != FCS %d rows]", nrow(csv_data), fcs_events))
    }
  }

  if (output_format == "wide") {
    # One file per population
    if (length(gating_cols) > 0) {
      for (col in gating_cols) {
        pop_cells <- csv_data[csv_data[[col]] == 1, ..expr_cols]
        pop_name <- clean_cluster_name(col)
        safe_name <- create_safe_filename(pop_name)
        out_file <- file.path(output_path,
                              paste0(gsub("\\.csv$", "", basename(csv_file)),
                                     "_", safe_name, ".csv"))
        fwrite(pop_cells, out_file, sep = ";")
      }
    }
  } else {
    # Long format: single file with population column
    if (length(gating_cols) > 0) {
      label_lists <- labels_from_binary_df(csv_data[, ..gating_cols])
      if (!is.null(gating_tree)) {
        csv_data$population <- flatten_label_list(label_lists, gating_tree)
      } else {
        csv_data$population <- sapply(label_lists, function(x) {
          if (length(x) == 0) "Unassigned" else x[1]
        })
      }
      if (!is.null(csv_label_dict)) {
        csv_data$population <- translate_labels(csv_data$population, csv_label_dict)
      }
    } else {
      csv_data$population <- "All"
    }
    out_file <- file.path(output_path, paste0("long_", basename(csv_file)))
    fwrite(csv_data[, c(..expr_cols, "population")], out_file, sep = ";")
  }

  cat(" -> done\n")
  list(file = csv_file, n_cells = nrow(csv_data), n_gating_cols = length(gating_cols))
}

# Wrapper for creating safe filenames (uses R_shared)
create_safe_filename <- function(name) {
  gsub("[^A-Za-z0-9_.-]", "_", gsub("_+", "_", name))
}
