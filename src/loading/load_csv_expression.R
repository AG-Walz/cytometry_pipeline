#!/usr/bin/env Rscript
# =============================================================================
# loading/load_csv_expression.R — Load CSV Expression Matrices
# =============================================================================
# Refactored from: R_scripts/01_loading/Load_exprs_csv.R
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))

#' Load expression CSV files from a directory
#'
#' @param input_dir Directory containing CSV expression files
#' @param sep CSV separator (auto-detected if NULL)
#' @param max_memory_gb Max memory before chunked reading
#' @return Named list of data.tables (one per CSV file)
load_csv_expression <- function(input_dir, sep = NULL, max_memory_gb = 4) {
  require(data.table)
  stopifnot(dir.exists(input_dir))

  csv_info <- list_csvs(input_dir)
  cat("Found", nrow(csv_info), "CSV files\n")

  results <- list()
  for (i in seq_len(nrow(csv_info))) {
    file_path <- csv_info$file_path[i]
    file_name <- csv_info$file_name[i]
    cat("Loading:", file_name, sprintf("(%.1f MB)...", csv_info$file_size_mb[i]))

    file_sep <- if (is.null(sep)) detect_csv_sep(file_path) else sep
    # Use fread directly for small files (respects sep); chunked reader for large files
    file_size_gb <- csv_info$file_size_mb[i] / 1000
    dt <- if (file_size_gb > max_memory_gb) {
      read_csv_chunked(file_path, max_memory_gb = max_memory_gb)
    } else {
      data.table::fread(file_path, sep = file_sep)
    }
    results[[file_name]] <- dt
    cat(sprintf(" %d rows x %d cols\n", nrow(dt), ncol(dt)))
  }

  cat("Loaded", length(results), "expression matrices\n")
  results
}

# --- CLI execution ---
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    cat("Usage: Rscript load_csv_expression.R <input_dir>\n")
    quit(status = 1)
  }
  result <- load_csv_expression(args[1])
  cat("Summary:\n")
  for (nm in names(result)) {
    cat(sprintf("  %s: %d x %d\n", nm, nrow(result[[nm]]), ncol(result[[nm]])))
  }
}
