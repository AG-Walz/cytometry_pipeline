# =============================================================================
# R_shared/memory_utils.R — Memory Management Utilities
# =============================================================================
# Deduplicated from: Fine_clustering_cytolution.R, GLMM.R, Functions.R
# =============================================================================

#' Check system memory usage and warn if high
#' @param warn_threshold Fraction (0-1) at which to warn (default 0.85)
#' @return Logical TRUE if memory is OK, FALSE if under pressure
check_memory <- function(warn_threshold = 0.85) {
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -m", intern = TRUE)
    mem_line <- mem_info[2]
    mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
    used_pct <- mem_values[2] / mem_values[1]
    if (used_pct > warn_threshold) {
      message("Warning: Memory usage at ", round(used_pct * 100, 1),
              "%. Running garbage collection...")
      gc()
      return(FALSE)
    }
  }
  TRUE
}

#' Print current memory usage summary
#' @param warn_threshold Fraction at which to trigger gc()
check_memory_usage <- function(warn_threshold = 0.8) {
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -m", intern = TRUE)
    mem_line <- mem_info[2]
    mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
    total_mb <- mem_values[1]
    used_mb <- mem_values[2]
    used_pct <- used_mb / total_mb
    cat(sprintf("Memory: %d / %d MB (%.1f%%)\n", used_mb, total_mb, used_pct * 100))
    if (used_pct > warn_threshold) {
      cat("Warning: Memory usage at", round(used_pct * 100, 1), "%\n")
      gc()
    }
  }
}

#' Memory-efficient CSV reading with chunking for large files
#' @param file_path Path to CSV file
#' @param chunk_size Number of rows per chunk
#' @param max_memory_gb Max memory threshold in GB before triggering chunked reading
#' @return data.table with file contents
read_csv_chunked <- function(file_path, chunk_size = 50000, max_memory_gb = 4) {
  require(data.table)
  file_size_mb <- file.info(file_path)$size / (1024^2)

  if (file_size_mb > max_memory_gb * 1000) {
    message("Large file detected (", round(file_size_mb, 1),
            " MB). Using chunked reading...")

    header <- data.table::fread(file_path, nrows = 0)

    # Count total rows
    total_rows <- tryCatch({
      wc_result <- system2("wc", args = c("-l", shQuote(file_path)),
                           stdout = TRUE, stderr = FALSE)
      if (length(wc_result) > 0 && !is.na(wc_result)) {
        as.numeric(stringr::str_extract(wc_result, "^\\d+"))
      } else stop("wc failed")
    }, error = function(e) {
      message("wc command failed, using R fallback...")
      tryCatch(
        length(readLines(file_path)),
        error = function(e2) round(file.size(file_path) / 1000)
      )
    })

    chunks <- list()
    n_chunks <- ceiling(total_rows / chunk_size)

    if (requireNamespace("progress", quietly = TRUE)) {
      pb <- progress::progress_bar$new(
        format = "Reading [:bar] :percent :current/:total chunks",
        total = n_chunks, clear = FALSE
      )
    } else {
      pb <- NULL
    }

    for (i in seq_len(n_chunks)) {
      skip_rows <- (i - 1) * chunk_size
      chunk <- data.table::fread(file_path, skip = skip_rows + 1,
                                 nrows = chunk_size, header = FALSE)
      if (nrow(chunk) > 0) {
        names(chunk) <- names(header)
        chunks[[i]] <- chunk
      }
      if (!is.null(pb)) pb$tick()
    }

    return(data.table::rbindlist(chunks, fill = TRUE))
  } else {
    return(data.table::fread(file_path))
  }
}

#' Check available system memory in GB
#' @return Numeric available memory in GB, or NA if cannot determine
get_available_memory_gb <- function() {
  if (Sys.info()["sysname"] == "Linux") {
    tryCatch({
      mem_info <- system("free -m", intern = TRUE)
      mem_line <- mem_info[2]
      mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
      available_mb <- mem_values[6]  # "available" column
      available_mb / 1024
    }, error = function(e) NA_real_)
  } else {
    NA_real_
  }
}
