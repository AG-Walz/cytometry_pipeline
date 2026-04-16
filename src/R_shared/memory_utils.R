# =============================================================================
# R_shared/memory_utils.R — Memory Management Utilities
# =============================================================================
# Deduplicated from: Fine_clustering_cytolution.R, GLMM.R, Functions.R
# Cross-platform: Linux, macOS, Windows
# =============================================================================

# %||% polyfill for R < 4.0
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

.get_os <- function() {
  Sys.info()[["sysname"]]
}

#' Query system memory on any OS (total_mb, used_mb, available_mb)
#' @return Named numeric vector with total_mb, used_mb, available_mb; or NULL on failure
.query_system_memory <- function() {
  os <- .get_os()
  tryCatch({
    if (os == "Linux") {
      mem_info <- system("free -m", intern = TRUE)
      mem_line <- mem_info[2]
      mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
      c(total_mb = mem_values[1], used_mb = mem_values[2], available_mb = mem_values[6])
    } else if (os == "Windows") {
      # Try PowerShell Get-CimInstance first (Win11+; wmic was deprecated and
      # may not be present). Fall back to wmic on older Windows.
      parse_kv <- function(x) {
        val_line <- grep("=", x, value = TRUE)
        if (length(val_line) == 0) return(NA_real_)
        as.numeric(sub(".*=", "", trimws(val_line[1])))
      }
      ps_query <- function() {
        ps_cmd <- paste(
          "$o = Get-CimInstance Win32_OperatingSystem;",
          "Write-Output ('TotalVisibleMemorySize=' + $o.TotalVisibleMemorySize);",
          "Write-Output ('FreePhysicalMemory=' + $o.FreePhysicalMemory)"
        )
        out <- tryCatch(
          system2("powershell", args = c("-NoProfile", "-Command", ps_cmd),
                  stdout = TRUE, stderr = FALSE),
          error = function(e) NULL
        )
        if (is.null(out) || length(out) == 0) return(NULL)
        total <- parse_kv(grep("^TotalVisibleMemorySize=", out, value = TRUE))
        free  <- parse_kv(grep("^FreePhysicalMemory=", out, value = TRUE))
        if (is.na(total) || is.na(free)) return(NULL)
        list(total_kb = total, free_kb = free)
      }
      wmic_query <- function() {
        out_t <- tryCatch(
          system2("wmic", args = c("OS", "get", "TotalVisibleMemorySize", "/value"),
                  stdout = TRUE, stderr = FALSE),
          error = function(e) NULL
        )
        out_f <- tryCatch(
          system2("wmic", args = c("OS", "get", "FreePhysicalMemory", "/value"),
                  stdout = TRUE, stderr = FALSE),
          error = function(e) NULL
        )
        if (is.null(out_t) || is.null(out_f)) return(NULL)
        total <- parse_kv(out_t)
        free  <- parse_kv(out_f)
        if (is.na(total) || is.na(free)) return(NULL)
        list(total_kb = total, free_kb = free)
      }
      kb <- ps_query() %||% wmic_query()
      if (is.null(kb)) return(NULL)
      total_mb <- kb$total_kb / 1024
      free_mb  <- kb$free_kb / 1024
      c(total_mb = total_mb, used_mb = total_mb - free_mb, available_mb = free_mb)
    } else if (os == "Darwin") {
      # macOS: use vm_stat and sysctl
      total_bytes <- as.numeric(system2("sysctl", args = c("-n", "hw.memsize"),
                                        stdout = TRUE, stderr = FALSE))
      total_mb <- total_bytes / (1024^2)
      vm <- system("vm_stat", intern = TRUE)
      page_size <- 16384  # default; try to parse
      ps_line <- grep("page size", vm, value = TRUE)
      if (length(ps_line) > 0) {
        ps_val <- as.numeric(gsub("[^0-9]", "", ps_line[1]))
        if (!is.na(ps_val) && ps_val > 0) page_size <- ps_val
      }
      free_line <- grep("Pages free", vm, value = TRUE)
      free_pages <- if (length(free_line) > 0) {
        as.numeric(gsub("[^0-9]", "", free_line[1]))
      } else 0
      inactive_line <- grep("Pages inactive", vm, value = TRUE)
      inactive_pages <- if (length(inactive_line) > 0) {
        as.numeric(gsub("[^0-9]", "", inactive_line[1]))
      } else 0
      available_mb <- (free_pages + inactive_pages) * page_size / (1024^2)
      c(total_mb = total_mb, used_mb = total_mb - available_mb, available_mb = available_mb)
    } else {
      NULL
    }
  }, error = function(e) NULL)
}

#' Check system memory usage and warn if high
#' @param warn_threshold Fraction (0-1) at which to warn (default 0.85)
#' @return Logical TRUE if memory is OK, FALSE if under pressure
check_memory <- function(warn_threshold = 0.85) {
  mem <- .query_system_memory()
  if (is.null(mem)) return(TRUE)  # cannot determine; assume OK
  used_pct <- mem["used_mb"] / mem["total_mb"]
  if (used_pct > warn_threshold) {
    message("Warning: Memory usage at ", round(used_pct * 100, 1),
            "%. Running garbage collection...")
    gc()
    return(FALSE)
  }
  TRUE
}

#' Print current memory usage summary
#' @param warn_threshold Fraction at which to trigger gc()
check_memory_usage <- function(warn_threshold = 0.8) {
  mem <- .query_system_memory()
  if (is.null(mem)) {
    cat("Memory: unable to query on this OS\n")
    return(invisible(NULL))
  }
  total_mb <- mem["total_mb"]
  used_mb  <- mem["used_mb"]
  used_pct <- used_mb / total_mb
  cat(sprintf("Memory: %d / %d MB (%.1f%%)\n", used_mb, total_mb, used_pct * 100))
  if (used_pct > warn_threshold) {
    cat("Warning: Memory usage at", round(used_pct * 100, 1), "%\n")
    gc()
  }
}

#' Count lines in a file (cross-platform)
#' @param file_path Path to file
#' @return Integer line count
.count_lines <- function(file_path) {
  os <- .get_os()
  # Try fast OS-native approach first
  result <- tryCatch({
    if (os == "Windows") {
      # Use find /c /v "" which counts lines on Windows
      cmd_out <- system2("cmd", args = c("/c", "find", "/c", "/v", '""',
                                         shQuote(normalizePath(file_path, winslash = "\\"))),
                         stdout = TRUE, stderr = FALSE)
      as.numeric(gsub("[^0-9]", "", tail(cmd_out, 1)))
    } else {
      wc_result <- system2("wc", args = c("-l", shQuote(file_path)),
                           stdout = TRUE, stderr = FALSE)
      as.numeric(stringr::str_extract(wc_result, "^\\d+"))
    }
  }, error = function(e) NA_real_)
  # Fallback: pure R line counting
  if (is.na(result)) {
    tryCatch(
      length(readLines(file_path, warn = FALSE)),
      error = function(e) round(file.size(file_path) / 1000)
    )
  } else {
    result
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
    total_rows <- .count_lines(file_path)

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
  mem <- .query_system_memory()
  if (is.null(mem)) return(NA_real_)
  mem["available_mb"] / 1024
}
