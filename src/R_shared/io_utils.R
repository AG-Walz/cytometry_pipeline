# =============================================================================
# R_shared/io_utils.R — I/O Utility Functions
# =============================================================================
# Deduplicated from: Functions.R, data_loading_functions.R, analysis scripts
# =============================================================================

#' Check if a path is a directory
#' @param path Character string path to check
#' @return Logical TRUE if path exists and is a directory
is_directory <- function(path) {
  if (!file.exists(path)) return(FALSE)
  file.info(path)$isdir
}

#' Create directories from a list of paths
#' @param dir_path_list Character vector of directory paths
dir_create <- function(dir_path_list) {
  for (path in dir_path_list) {
    path <- normalizePath(path, mustWork = FALSE)
    if (!is_directory(path)) {
      dir.create(path, showWarnings = FALSE, recursive = TRUE)
    }
  }
}

#' Detect CSV separator (semicolon vs comma)
#' @param file Path to CSV file
#' @return Character ";" or ","
detect_csv_sep <- function(file) {
  first_line <- readLines(file, n = 1)
  if (grepl(";", first_line)) ";" else ","
}

#' Get file size in MB
#' @param file_path Path to file
#' @return Numeric file size in MB
get_file_size_mb <- function(file_path) {
  file.info(file_path)$size / 1024^2
}

#' Save data.table result and clean up memory
#' @param result data.table to save
#' @param output_file Output path
#' @param cleanup_vars Optional character vector of variable names to remove from parent frame
#' @return output_file path (invisible)
save_and_cleanup <- function(result, output_file, cleanup_vars = NULL) {
  require(data.table)
  data.table::fwrite(result, output_file, sep = ";")
  if (!is.null(cleanup_vars)) {
    for (var in cleanup_vars) {
      if (exists(var, envir = parent.frame())) {
        rm(list = var, envir = parent.frame())
      }
    }
  }
  gc()
  invisible(output_file)
}

#' Format seconds into human-readable time strings
#' @param x Numeric vector of seconds
#' @return Character vector of formatted strings
time_formatter <- function(x) {
  sapply(x, function(val) {
    if (val < 60) sprintf("%.0fs", val)
    else if (val < 3600) sprintf("%.0fmin", val / 60)
    else sprintf("%.1fh", val / 3600)
  })
}

#' List CSV files with metadata
#' @param dir_path Directory to scan
#' @return tibble with file_path, file_name, file_size_mb columns
list_csvs <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  if (length(files) == 0) stop("No .csv files found in: ", dir_path)
  result <- tibble::tibble(
    file_path = files,
    file_name = tools::file_path_sans_ext(basename(files)),
    file_size_mb = round(file.info(files)$size / 1024^2, 2)
  )
  dplyr::arrange(result, file_name)
}

#' Ensure packages are installed and loaded
#' @param req_packages Character vector of package names
ensure_and_load_packages <- function(req_packages) {
  for (pkg in req_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#' Remove columns matching a pattern from all CSVs in a folder
#' @param folder_path Path to folder with CSV files
#' @param colname_pattern Regex pattern to match column names for removal
#' @param output_folder Output folder (default: folder_path/filtered)
#' @return List of results per file (invisible)
remove_columns_from_csvs <- function(folder_path, colname_pattern, output_folder = NULL) {
  require(data.table)
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the specified folder")

  output_folder <- output_folder %||% file.path(folder_path, "filtered")
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

  results <- purrr::map(csv_files, ~ {
    cat("Processing:", basename(.x), "...")
    data <- data.table::fread(.x, sep = ";", header = TRUE)
    matching_cols <- grep(colname_pattern, names(data), value = TRUE)
    if (length(matching_cols) > 0) {
      data_filtered <- data[, !names(data) %like% colname_pattern, with = FALSE]
      cat(" removed", length(matching_cols), "columns")
    } else {
      data_filtered <- data
      cat(" no matches found")
    }
    output_file <- file.path(output_folder, basename(.x))
    data.table::fwrite(data_filtered, output_file, sep = ";")
    cat(" ->", ncol(data_filtered), "cols saved\n")
    list(original_cols = ncol(data), filtered_cols = ncol(data_filtered),
         removed_cols = matching_cols)
  }) %>% stats::setNames(basename(csv_files))

  cat("Summary: Processed", length(csv_files), "files -", output_folder, "\n")
  invisible(results)
}

#' Classify Cytolution CSV columns by type
#'
#' Auto-classifies column names from Cytolution CSV exports into categories
#' based on known naming patterns.
#'
#' @param col_names Character vector of column names
#' @param gating_pattern Regex pattern identifying gating/clustering columns
#' @return Named list with character vectors: index, scatter, markers, gating,
#'   dimreduction, autofluorescence
classify_csv_columns <- function(col_names,
                                  gating_pattern = "\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$") {
  # Index columns
  index_cols <- grep("^[Ii]ndex$", col_names, value = TRUE)

  # Scatter/time columns: SSC, FSC, TIME (with separator variants)
  scatter_cols <- grep("^(SSC|FSC|TIME)[\\.|\\s\\|]", col_names, value = TRUE, perl = TRUE)

  # Gating/clustering columns
  gating_cols <- grep(gating_pattern, col_names, value = TRUE)

  # Dimensionality reduction columns (Cytolution UMAP coords)
  dimreduction_cols <- grep("^cyto_2d_", col_names, value = TRUE)

  # Autofluorescence channels: AF.P<number> pattern or where marker name == channel name
  # e.g. "AF.P1 | AF.P1.A" or columns where both sides of " | " share the same root
  af_cols <- grep("^AF[\\.-]", col_names, value = TRUE, perl = TRUE)

  # Also catch columns where marker == channel (e.g. "SSC.H | SSC.H" already caught,

  # but some channels lack a marker: "BUV496.A | BUV496.A")
  pipe_cols <- grep(" \\| ", col_names, value = TRUE)
  for (col in pipe_cols) {
    parts <- trimws(strsplit(col, " \\| ")[[1]])
    if (length(parts) == 2) {
      # Strip trailing .A/.H/.W from BOTH sides so "BUV496.A | BUV496.A" is
      # detected as autofluorescence (marker root matches channel root).
      marker_root  <- gsub("\\.[AHW]$", "", parts[1])
      channel_root <- gsub("\\.[AHW]$", "", parts[2])
      if (marker_root == channel_root && !col %in% c(scatter_cols, gating_cols, af_cols)) {
        af_cols <- c(af_cols, col)
      }
    }
  }

  # Already classified columns
  classified <- c(index_cols, scatter_cols, gating_cols, dimreduction_cols, af_cols)

  # Markers: everything not yet classified

  marker_cols <- setdiff(col_names, classified)

  list(
    index         = index_cols,
    scatter       = scatter_cols,
    markers       = marker_cols,
    gating        = gating_cols,
    dimreduction  = dimreduction_cols,
    autofluorescence = af_cols
  )
}

#' Stop with message and exit (for scripts)
#' @param ... Message parts passed to message()
die <- function(...) {
  message(...)
  quit(status = 1)
}

#' Create output directory if it doesn't exist
#' @param path Directory path
make_outdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Machine-readable progress, for a GUI that is watching the log.
#
# A long run's only honest progress signal is the work itself: how many
# workspaces have been built, how many samples extracted, how many plots
# written. Only the analysis knows those numbers, so it declares them and
# whoever is watching renders them. The alternative — a front-end counting
# lines of prose it recognises — is a parser for a format nobody agreed to, and
# it breaks the first time someone rewords a message.
#
# One line, prefixed so it can be told apart from anything a human is meant to
# read, and so a front-end can hide it from the visible log:
#
#   ::progress:: phase=extract done=57 total=311 elapsed=812 label=MTB023_204275
#
# Emitted to stdout in the same stream as everything else, because the runner
# already line-buffers and streams that, and a side-channel file would need its
# own polling and its own staleness problem.
#
# `elapsed` is seconds since this process started. It is here because the GUI
# re-reads the WHOLE log on every poll to recover the state of a detached run —
# it keeps nothing in memory that survives a browser reload — so how long a
# phase took has to be readable from the log itself. Without it, closing the tab
# during a two-hour run loses every timing behind it and with them any estimate
# of what is left. It goes BEFORE `label=` because a label may contain anything
# at all, including "elapsed=", so `label` has to stay the last field.
#
# Silent when `total` is unknown or zero: a progress bar that cannot say what
# fraction is done is worse than none.
# -----------------------------------------------------------------------------
PROGRESS_PREFIX <- "::progress::"

#' Seconds since the R process started, for the `elapsed` field.
#'
#' proc.time()["elapsed"] is measured from the start of THIS process, so it
#' needs no shared epoch with the reader and cannot be skewed by the clock
#' changing under a long run.
progress_elapsed <- function() {
  as.numeric(proc.time()[["elapsed"]])
}

emit_progress <- function(phase, done, total, label = "") {
  total <- suppressWarnings(as.integer(total))
  done  <- suppressWarnings(as.integer(done))
  if (is.na(total) || is.na(done) || total <= 0L) return(invisible(NULL))
  # A label can be any sample name, and those contain spaces, '=' and worse.
  # Newlines are the only character that could break the one-line contract.
  label <- gsub("[\r\n]", " ", as.character(label))
  cat(sprintf("%s phase=%s done=%d total=%d elapsed=%.1f label=%s\n",
              PROGRESS_PREFIX, phase, done, total, progress_elapsed(), label))
  flush(stdout())
  invisible(NULL)
}
