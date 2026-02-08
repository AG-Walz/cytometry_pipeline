# =============================================================================
# R_shared/fcs_utils.R — FCS File Handling Utilities
# =============================================================================
# Deduplicated from: Functions.R, data_loading_functions.R
# =============================================================================

#' Read the text segment of an FCS file
#' @param con Connection object to FCS file
#' @param offsets Named vector with textstart, textend, FCSversion
#' @param emptyValue Include empty values in output
#' @param cpp Use C++ parser if available
#' @return Named character vector of text segment key-value pairs
readFCStext <- function(con, offsets, emptyValue = TRUE, cpp = TRUE, ...) {
  seek(con, offsets["textstart"])
  txt <- iconv(
    rawToChar(readBin(con, "raw", offsets["textend"] - offsets["textstart"] + 1)),
    "", "latin1", sub = "byte"
  )
  txt <- trimws(txt, "right")

  if (offsets["FCSversion"] <= 2) {
    delimiter <- substr(txt, 1, 1)
    sp <- strsplit(substr(txt, 2, nchar(txt)), split = delimiter, fixed = TRUE)[[1]]
    rv <- c(offsets["FCSversion"], sp[seq(2, length(sp), by = 2)])
    names(rv) <- gsub("^ *| *$", "", c("FCSversion", sp[seq(1, length(sp) - 1, by = 2)]))
  } else {
    rv <- if (cpp) {
      flowCore::fcsTextParse(txt, emptyValue = emptyValue)
    } else {
      flowCore::fcs_text_parse(txt, emptyValue = emptyValue)
    }
    rv <- c(offsets["FCSversion"], rv)
    names(rv) <- gsub("^ *| *$", "", names(rv))
  }
  rv
}

#' Retrieve parameters from an FCS text segment
#' @param x Named character vector (FCS text segment)
#' @param pnam Character vector of parameter names to retrieve
#' @param strict If TRUE, stop on missing parameters
#' @return Character vector of parameter values
readFCSgetPar <- function(x, pnam, strict = TRUE) {
  stopifnot(is.character(x), is.character(pnam))
  i <- match(pnam, names(x))
  if (any(is.na(i)) && strict) {
    stop(paste("Parameter(s)", pnam[is.na(i)], "not contained in 'x'\n"))
  }
  if (!strict) {
    i[!is.na(i)] <- x[i[!is.na(i)]]
    names(i) <- pnam
    return(i)
  }
  x[i]
}

#' Check $TOT keyword vs actual event count
#' @param flowframe A flowFrame object
#' @param index Optional sample index for messages
check_tot <- function(flowframe, index = "") {
  tot_keyword <- as.numeric(flowCore::keyword(flowframe, "$TOT")[[1]])
  n_events <- nrow(flowframe)
  if (tot_keyword == n_events) {
    message("$TOT matches event count (", n_events, ") for sample ", index)
  } else {
    warning("$TOT (", tot_keyword, ") != event count (", n_events,
            ") for sample ", index)
  }
}

#' Create CSV to FCS file mapping based on filename similarity
#' @param csv_files Character vector of CSV file paths
#' @param fcs_files Character vector of FCS file paths
#' @return Named character vector (names = csv paths, values = fcs paths)
create_csv_fcs_mapping <- function(csv_files, fcs_files) {
  csv_to_fcs <- character(length(csv_files))
  names(csv_to_fcs) <- csv_files
  fcs_bases_norm <- gsub("[-_]", "", gsub("\\.fcs$", "", basename(fcs_files)))

  for (i in seq_along(csv_files)) {
    csv_base_norm <- gsub("[-_]", "", gsub("\\.csv$", "", basename(csv_files[i])))
    match_idx <- which(sapply(fcs_bases_norm, function(fcs_norm) {
      grepl(fcs_norm, csv_base_norm, fixed = TRUE)
    }))
    if (length(match_idx) != 1) {
      stop(sprintf("Could not uniquely match FCS file for CSV '%s' (found %d matches)",
                   basename(csv_files[i]), length(match_idx)))
    }
    csv_to_fcs[csv_files[i]] <- fcs_files[match_idx]
  }
  csv_to_fcs
}

#' Get FCS event count without loading full file
#' @param fcs_file Path to FCS file
#' @return Integer event count
get_fcs_event_count <- function(fcs_file) {
  tryCatch({
    ff <- flowCore::read.FCS(fcs_file, truncate_max_range = FALSE,
                             transformation = FALSE, which.lines = 1)
    as.numeric(flowCore::keyword(ff, "$TOT")[[1]])
  }, error = function(e) {
    ff <- flowCore::read.FCS(fcs_file, truncate_max_range = FALSE)
    nrow(ff)
  })
}

#' Check if CSV has manual gating columns
#' @param csv_file Path to CSV file
#' @param pattern Regex pattern for gating columns
#' @return Logical
has_manual_gating <- function(csv_file, pattern = "Manual gating") {
  tryCatch({
    separator <- detect_csv_sep(csv_file)
    header <- data.table::fread(csv_file, sep = separator, nrows = 0)
    any(grepl(pattern, names(header)))
  }, error = function(e) FALSE)
}

#' Filter CSV files to those containing manual gating columns
#' @param csv_files Character vector of CSV file paths
#' @param pattern Regex pattern for gating columns
#' @return Filtered character vector
filter_valid_csv_files <- function(csv_files, pattern = "Manual gating") {
  valid_files <- csv_files[sapply(csv_files, has_manual_gating, pattern = pattern)]
  cat("Found", length(valid_files), "valid CSV files with pattern:", pattern, "\n")
  valid_files
}

#' Read CSV and filter to essential columns
#' @param csv_file Path to CSV file
#' @param pattern Regex pattern for columns to keep
#' @return data.table with filtered columns, or NULL
read_csv_filtered <- function(csv_file, pattern = "Manual gating") {
  require(data.table)
  separator <- detect_csv_sep(csv_file)
  csv_data <- data.table::fread(csv_file, sep = separator, header = TRUE)
  manual_cols <- grep(pattern, names(csv_data), value = TRUE)
  if (length(manual_cols) == 0) return(NULL)
  essential_cols <- c("index", manual_cols)
  available_cols <- intersect(essential_cols, names(csv_data))
  csv_data[, ..available_cols]
}

#' Load a list of FCS files into flowFrames
#' @param fcs_file_list Character vector of FCS file paths
#' @return List of flowFrame objects
load_fcs_list <- function(fcs_file_list) {
  purrr::map2(fcs_file_list, seq_along(fcs_file_list), ~ {
    cat("Loading", basename(.x), "...")
    ff <- flowWorkspace::load_cytoframe_from_fcs(.x, transformation = FALSE,
                                                 truncate_max_range = FALSE)
    cat(" Sample", .y, "has", nrow(ff), "cells and", ncol(ff), "markers\n")
    ff
  })
}

#' Load all FCS files from a folder
#' @param input_file_path Path to folder containing FCS files
#' @return List of flowFrame objects
load_fcs_folder <- function(input_file_path) {
  cat("Reading FCS-files from", input_file_path, "\n")
  input_fcs_files <- list.files(path = input_file_path,
                                pattern = "\\.fcs$", full.names = TRUE)
  sample_names <- input_fcs_files %>%
    basename() %>%
    gsub("\\.fcs$", "", .) %>%
    sapply(function(x) strsplit(x, "_")[[1]][1:2] %>% paste(collapse = "_"))
  cat("Found", length(input_fcs_files), "FCS files in", input_file_path, "\n")
  cat("Sample names:", paste(sample_names, collapse = ", "), "\n")
  load_fcs_list(input_fcs_files)
}

# Note: detect_csv_sep is in io_utils.R. Source io_utils.R before this file.
