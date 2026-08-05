# =============================================================================
# R_shared/memory_utils.R — Memory Management Utilities
# =============================================================================
# Deduplicated from: Fine_clustering_cytolution.R, GLMM.R, Functions.R
# Cross-platform: Linux, macOS, Windows
# =============================================================================

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
      # Use wmic for total and free physical memory (returns KB)
      total_kb <- system2("wmic", args = c("OS", "get", "TotalVisibleMemorySize", "/value"),
                          stdout = TRUE, stderr = FALSE)
      free_kb  <- system2("wmic", args = c("OS", "get", "FreePhysicalMemory", "/value"),
                          stdout = TRUE, stderr = FALSE)
      parse_wmic <- function(x) {
        val_line <- grep("=", x, value = TRUE)
        if (length(val_line) == 0) return(NA_real_)
        as.numeric(sub(".*=", "", trimws(val_line[1])))
      }
      total_mb <- parse_wmic(total_kb) / 1024
      free_mb  <- parse_wmic(free_kb) / 1024
      if (is.na(total_mb) || is.na(free_mb)) return(NULL)
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

# =============================================================================
# CPU: how many workers a run may actually use
# =============================================================================
# detectCores() reports the HOST's CPUs, which is the wrong number twice over on
# an analysis server: it ignores a cgroup CPU limit (so a container allotted two
# cores happily forks sixteen and thrashes), and it ignores the fact that other
# people are using the machine. Both matter here — highmem2 is shared, and the
# app is served from a container by PBIAppLauncher.
# =============================================================================

#' The CPU quota this process is actually allowed, or NA when unlimited.
#'
#' cgroup v2 writes "$MAX $PERIOD" (or "max $PERIOD") to cpu.max; v1 splits it
#' across cpu.cfs_quota_us and cpu.cfs_period_us with -1 meaning unlimited.
#' Both are read because hosts in the wild still run either.
.cgroup_cpu_quota <- function() {
  v2 <- "/sys/fs/cgroup/cpu.max"
  if (file.exists(v2)) {
    parts <- tryCatch(strsplit(trimws(readLines(v2, warn = FALSE))[1], "\\s+")[[1]],
                      error = function(e) NULL)
    if (length(parts) == 2L && parts[1] != "max") {
      q <- suppressWarnings(as.numeric(parts[1]))
      p <- suppressWarnings(as.numeric(parts[2]))
      if (!is.na(q) && !is.na(p) && p > 0) return(q / p)
    }
    return(NA_real_)
  }
  q_file <- "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
  p_file <- "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
  if (file.exists(q_file) && file.exists(p_file)) {
    q <- suppressWarnings(as.numeric(readLines(q_file, warn = FALSE)[1]))
    p <- suppressWarnings(as.numeric(readLines(p_file, warn = FALSE)[1]))
    if (!is.na(q) && q > 0 && !is.na(p) && p > 0) return(q / p)
  }
  NA_real_
}

#' How many workers to use.
#'
#' @param requested 0 or NULL = decide automatically; N = use N.
#' @return An integer >= 1.
#'
#' Automatic leaves one core free. On a shared host that is the difference
#' between a run that finishes and a run everyone else notices; on a laptop it
#' is the difference between working through it and watching a frozen cursor.
#' An explicit N is honoured but still capped at what the cgroup allows, because
#' asking for more workers than the quota does not buy throughput — it buys
#' context switching.
resolve_cores <- function(requested = 0L) {
  detected <- tryCatch(parallel::detectCores(logical = TRUE),
                       error = function(e) NA_integer_)
  if (is.na(detected) || detected < 1L) detected <- 1L
  quota <- .cgroup_cpu_quota()
  allowed <- if (is.na(quota)) detected else max(1L, min(detected, floor(quota)))

  requested <- suppressWarnings(as.integer(requested %||% 0L))
  if (is.na(requested) || requested <= 0L) {
    return(as.integer(max(1L, allowed - 1L)))
  }
  as.integer(max(1L, min(requested, allowed)))
}

#' lapply, over `n_cores` forked workers.
#'
#' n_cores = 1 runs plain lapply — same code path as before this existed, which
#' is what makes "did parallelism break it?" answerable by setting one flag.
#'
#' mclapply returns a try-error OBJECT rather than raising, so a failure in a
#' worker would otherwise be discovered later as a corrupt result. Every element
#' is checked and the first failure is re-raised with its original message.
#'
#' Forking only: the payload here is large shared data (a coordinate frame and an
#' expression matrix) that copy-on-write shares for free, where a socket cluster
#' would serialise it to every worker. That does mean Linux and macOS only —
#' this repo's stated target — and Windows falls back to serial.
parallel_lapply <- function(X, FUN, n_cores = 1L, ...) {
  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores <= 1L || .Platform$OS.type != "unix" ||
      length(X) <= 1L) {
    return(lapply(X, FUN, ...))
  }
  res <- parallel::mclapply(X, FUN, ..., mc.cores = min(n_cores, length(X)),
                            mc.preschedule = TRUE)

  # An R-level error comes back as a try-error object and is re-raised with its
  # own message.
  failed <- vapply(res, function(r) inherits(r, "try-error"), logical(1))
  if (any(failed)) {
    stop("a parallel worker failed: ",
         conditionMessage(attr(res[[which(failed)[1]]], "condition")),
         call. = FALSE)
  }

  # A worker killed by a SIGNAL does not. mclapply returns NULL for every
  # element that worker was carrying, with nothing to distinguish it from a
  # function that legitimately returned NULL — and the commonest way to be
  # killed by a signal here is the OOM killer, on exactly the large run where
  # forking many workers is most tempting. Silently losing a third of the plots
  # to that is the worst outcome available, so the count is checked: mclapply
  # always returns one element per input, so a short result or a NULL where the
  # caller expects otherwise means a worker died.
  if (length(res) != length(X)) {
    stop("a parallel worker was killed (", length(res), " of ", length(X),
         " results came back). The usual cause is the out-of-memory killer: ",
         "each of the ", min(n_cores, length(X)), " workers needs its own ",
         "working memory. Re-run with fewer cores.", call. = FALSE)
  }
  n_null <- sum(vapply(res, is.null, logical(1)))
  if (n_null > 0L && !isTRUE(attr(FUN, "returns_null"))) {
    stop(n_null, " of ", length(X), " parallel task(s) returned nothing, which ",
         "is what a worker killed by a signal looks like. The usual cause is ",
         "the out-of-memory killer. Re-run with fewer cores.", call. = FALSE)
  }
  res
}
