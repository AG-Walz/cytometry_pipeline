# =============================================================================
# R_shared/run_record.R — what a run was asked to do, written before it starts
# =============================================================================
# run_info.json (written at the END of a UMAP run) says what a run found. It
# could not answer "what exactly did that run use?" for the runs that matter
# most when debugging:
#
#   * a run that crashed or was killed never got to write it;
#   * it recorded the argv as typed, so every option left at its default was
#     missing -- and defaults change (2026-11-03: fit cells, output formats,
#     the UMAP start), so replaying an old argv quietly gives a different run;
#   * a re-run into the same folder overwrote it.
#
# So every run now writes, before any work:
#
#   _run/settings/<YYYYmmdd-HHMMSS>/settings.json   every option's effective
#                                                    value (given or default),
#                                                    code version, packages,
#                                                    host, inputs, outcome
#   _run/settings/<YYYYmmdd-HHMMSS>/command.sh      the exact command, with
#                                                    EVERY option spelt out
#   _run/settings/<...>/uncommitted.diff            when the checkout had local
#                                                    changes
#   _run/settings.json, _run/command.sh             copies of the latest
#
# and settings.json gets the outcome (completed / failed + the error) when the
# run ends. `bash _run/command.sh` repeats the run with the same settings into
# the same folder; edit --output_dir in it to keep both.
# =============================================================================

RUN_RECORD_DIR <- "_run"

# Options that describe the invocation, not the analysis; never replayed.
.RUN_RECORD_SKIP <- c("help")

#' Every option of an optparse option_list with its effective value.
#'
#' @param opt Parsed options (parse_args()).
#' @param option_list The option_list given to OptionParser().
#' @param argv The arguments as given (commandArgs(trailingOnly = TRUE)),
#'   used only to mark which options were typed rather than defaulted.
#' @return A list of list(flag, dest, value, given, is_flag), in option order.
run_record_options <- function(opt, option_list, argv = character(0)) {
  given_flags <- sub("=.*$", "", argv[startsWith(argv, "--")])
  out <- list()
  for (o in option_list) {
    flag <- methods::slot(o, "long_flag")
    dest <- methods::slot(o, "dest")
    if (!nzchar(flag) || dest %in% .RUN_RECORD_SKIP) next
    action <- methods::slot(o, "action")
    out[[length(out) + 1L]] <- list(
      flag = flag, dest = dest,
      value = opt[[dest]],
      given = flag %in% given_flags,
      is_flag = action %in% c("store_true", "store_false"),
      action = action)
  }
  out
}

#' The command that repeats a run: every option written out, defaults too.
#'
#' A store_true flag appears when TRUE; a store_false flag when FALSE. An
#' option whose value is NULL (no default, not given) is left out, which is
#' what leaving it out meant the first time.
#'
#' @param script Path of the R script that was run.
#' @param options Output of run_record_options().
#' @param rscript The R front end.
#' @return One string, shell-quoted, continued over lines with backslashes.
run_record_command <- function(script, options, rscript = "Rscript") {
  parts <- c(shQuote(rscript), shQuote(normalizePath(script, mustWork = FALSE)))
  for (o in options) {
    v <- o$value
    if (is.null(v) || (length(v) == 1L && is.na(v))) next
    if (isTRUE(o$is_flag)) {
      on <- if (identical(o$action, "store_false")) identical(v, FALSE)
            else identical(v, TRUE)
      if (on) parts <- c(parts, o$flag)
      next
    }
    if (length(v) > 1L) v <- paste(v, collapse = ",")
    if (is.numeric(v)) v <- format(v, scientific = FALSE, trim = TRUE,
                                   digits = 15)
    # --flag=value, never "--flag value": a value that starts with "-" (an
    # automatic -1) would otherwise be read back as another flag.
    parts <- c(parts, paste0(o$flag, "=", shQuote(as.character(v))))
  }
  paste(parts, collapse = " \\\n  ")
}

# Size, mtime and md5 of a file; size and mtime only above `md5_limit` bytes
# or when md5 = FALSE (an FCS folder of hundreds of files is fingerprinted by
# its listing, not by reading every byte).
.run_record_file <- function(path, md5 = TRUE, md5_limit = 2e8) {
  info <- file.info(path)
  rec <- list(path = normalizePath(path, mustWork = FALSE),
              bytes = unname(info$size),
              modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"))
  if (md5 && isTRUE(info$size <= md5_limit)) {
    rec$md5 <- unname(tools::md5sum(path))
  }
  rec
}

#' Fingerprint the inputs: each file given, and each folder's listing.
#'
#' @param inputs Named list of paths (files or folders), NULLs ignored.
#' @param dir A run settings directory; folder listings are written there as
#'   <name>_files.tsv (name, bytes, modified), so a later "which samples were
#'   there?" has an answer.
.run_record_inputs <- function(inputs, dir) {
  out <- list()
  for (name in names(inputs)) {
    p <- inputs[[name]]
    if (is.null(p) || !nzchar(p) || !file.exists(p)) next
    if (dir.exists(p)) {
      files <- list.files(p, recursive = TRUE, full.names = TRUE)
      info <- file.info(files)
      listing <- data.frame(file = sub(paste0("^", normalizePath(p), "/?"), "",
                                       normalizePath(files)),
                            bytes = info$size,
                            modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S"),
                            stringsAsFactors = FALSE)
      tsv <- file.path(dir, paste0(name, "_files.tsv"))
      utils::write.table(listing, tsv, sep = "\t", row.names = FALSE,
                         quote = FALSE)
      # Small files (workspaces, metadata) are hashed as well.
      small <- files[!is.na(info$size) & info$size <= 2e8 &
                     grepl("\\.(wsp|csv|tsv|txt|json|ya?ml)$", files,
                           ignore.case = TRUE)]
      out[[name]] <- list(path = normalizePath(p), folder = TRUE,
                          files = length(files),
                          bytes = sum(info$size, na.rm = TRUE),
                          listing = basename(tsv),
                          hashed = lapply(small, .run_record_file))
    } else {
      out[[name]] <- .run_record_file(p)
    }
  }
  out
}

# git state of the checkout holding src_dir: commit, branch, local changes.
.run_record_git <- function(repo, dir) {
  git <- function(args) tryCatch(
    suppressWarnings(system2("git", shQuote(c("-C", repo, args)), stdout = TRUE,
                             stderr = FALSE)),
    error = function(e) structure(character(0), status = 1L))
  commit <- git(c("rev-parse", "HEAD"))
  if (length(commit) != 1L || !is.null(attr(commit, "status"))) {
    return(list(repo = repo, commit = NA_character_))
  }
  status <- git(c("status", "--porcelain"))
  rec <- list(repo = repo, commit = commit,
              branch = paste(git(c("rev-parse", "--abbrev-ref", "HEAD")),
                             collapse = ""),
              describe = paste(git(c("log", "-1", "--format=%h %cI %s")),
                               collapse = ""),
              dirty = length(status) > 0L)
  if (rec$dirty) {
    diff <- git(c("diff", "HEAD"))
    writeLines(c(paste("#", status), diff), file.path(dir, "uncommitted.diff"))
    rec$uncommitted <- "uncommitted.diff"
  }
  rec
}

#' Write the settings record of a run that is about to start.
#'
#' @param output_dir The run's output folder.
#' @param script Path of the R script being run.
#' @param opt,option_list Parsed options and the parser's option list.
#' @param inputs Named list of input paths to fingerprint (files or folders).
#' @param repo The repository checkout the code comes from.
#' @param packages Package names whose versions are recorded.
#' @return The record (a list with $dir), for finish_run_record(). NULL when it
#'   could not be written -- a missing record must not stop the analysis.
write_run_record <- function(output_dir, script, opt, option_list,
                             inputs = list(), repo = getwd(),
                             packages = character(0)) {
  tryCatch({
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      message("Settings record not written: install jsonlite.")
      return(invisible(NULL))
    }
    started <- Sys.time()
    stamp <- format(started, "%Y%m%d-%H%M%S")
    run_dir <- file.path(output_dir, RUN_RECORD_DIR)
    dir <- file.path(run_dir, "settings", stamp)
    # Two runs started in the same second (a retry) keep both records.
    k <- 1L
    while (dir.exists(dir)) {
      k <- k + 1L
      dir <- file.path(run_dir, "settings", paste0(stamp, "-", k))
    }
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    argv <- commandArgs(trailingOnly = TRUE)
    options <- run_record_options(opt, option_list, argv)
    command <- run_record_command(script, options)
    env_keep <- c("TMPDIR", "OMP_NUM_THREADS", "CONDA_DEFAULT_ENV",
                  "CONDA_PREFIX", "R_LIBS_USER", "UMAP_APP_TMPDIR")
    env <- Sys.getenv(env_keep, unset = NA)
    env <- as.list(env[!is.na(env) & nzchar(env)])
    pkg <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                                error = function(e) NA_character_)

    record <- list(
      status = "running",
      started = format(started, "%Y-%m-%dT%H:%M:%S%z"),
      script = normalizePath(script, mustWork = FALSE),
      command = "command.sh",
      argv_as_given = I(argv),
      options = stats::setNames(lapply(options, function(o) o$value),
                                vapply(options, function(o) o$dest, "")),
      options_given = I(vapply(Filter(function(o) o$given, options),
                               function(o) o$dest, "")),
      git = .run_record_git(repo, dir),
      R = R.version.string,
      packages = lapply(stats::setNames(nm = packages), pkg),
      host = list(name = Sys.info()[["nodename"]], user = Sys.info()[["user"]],
                  cores = parallel::detectCores(), cwd = getwd()),
      environment = env,
      inputs = .run_record_inputs(inputs, dir))

    header <- c(
      "#!/usr/bin/env bash",
      paste0("# Repeats the run started ", record$started, " on ",
             record$host$name, " by ", record$host$user, "."),
      paste0("# Code: ", if (is.null(record$git$describe)) "(not a git checkout)"
                         else record$git$describe,
             if (isTRUE(record$git$dirty)) "  -- WITH LOCAL CHANGES, see uncommitted.diff" else ""),
      if (!is.na(record$git$commit))
        paste0("# Same code: git -C ", shQuote(repo), " checkout ",
               record$git$commit),
      "# Every option is written out, defaults included, so a later change of",
      "# a default does not change what this repeats.",
      "set -euo pipefail",
      paste("cd", shQuote(getwd())),
      if (!is.null(env$TMPDIR)) paste0("export TMPDIR=", shQuote(env$TMPDIR)),
      "")
    writeLines(c(header, command), file.path(dir, "command.sh"))
    Sys.chmod(file.path(dir, "command.sh"), "755")
    record$dir <- normalizePath(dir)
    .write_run_record_json(record)
    file.copy(file.path(dir, "command.sh"), file.path(run_dir, "command.sh"),
              overwrite = TRUE)
    cat("Settings recorded: ", file.path(run_dir, "settings.json"),
        " (repeat with: bash ", file.path(run_dir, "command.sh"), ")\n",
        sep = "")
    invisible(record)
  }, error = function(e) {
    message("Settings record not written: ", conditionMessage(e))
    invisible(NULL)
  })
}

#' Record how a run ended. Safe to call with record = NULL.
#'
#' @param record What write_run_record() returned.
#' @param status "completed" or "failed".
#' @param error The condition, when it failed.
finish_run_record <- function(record, status, error = NULL) {
  if (is.null(record)) return(invisible(NULL))
  tryCatch({
    finished <- Sys.time()
    record$status <- status
    record$finished <- format(finished, "%Y-%m-%dT%H:%M:%S%z")
    record$minutes <- round(as.numeric(difftime(
      finished, as.POSIXct(record$started, format = "%Y-%m-%dT%H:%M:%S%z"),
      units = "mins")), 1)
    if (!is.null(error)) record$error <- conditionMessage(error)
    .write_run_record_json(record)
  }, error = function(e) NULL)
  invisible(record)
}

.write_run_record_json <- function(record) {
  json <- jsonlite::toJSON(record[setdiff(names(record), "dir")],
                           auto_unbox = TRUE, pretty = TRUE, null = "null",
                           na = "null", digits = NA)
  writeLines(json, file.path(record$dir, "settings.json"))
  writeLines(json, file.path(dirname(dirname(record$dir)), "settings.json"))
}
