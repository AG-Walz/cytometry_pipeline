#!/usr/bin/env Rscript
# =============================================================================
# metadata_qc/harmonize_marker_names.R : one spelling per marker across a cohort
# =============================================================================
# A cohort's FCS files name the same antibody differently ($PnS = INFy / IFN /
# IFNg, TNFa / TNF, ALD / Aqua LiveDead), and many carry no $PnS at all. Every
# downstream step then patches the names with its own regexes. This module
# rewrites $PnS in the files once, reversibly, so the rest of the repo can
# address channels by one canonical marker name.
#
# Base R only. Sources nothing itself: the caller provides
#   src/R_shared/fcs_keywords.R                 fcs_keywords, fcs_inspect,
#                                               fcs_write_keywords, fcs_text_get
#   src/harmonization/panel_homogenization.R    marker_alias_key,
#                                               suggest_marker_aliases,
#                                               apply_marker_aliases,
#                                               parse_marker_aliases
#   src/harmonization/drift/drift_panel.R       fcs_parameter_table,
#                                               assign_canonical
#
# Safety posture mirrors sync_fcs_names.R: dry run by default, in_place is
# explicit, symlinks are refused, event bytes are never decoded, and every
# replaced value is recorded in a CSV ledger that revert_marker_plan() replays.
# Only $PnS is written: $PnN and the structural parameter keywords stay as they
# are, so gates that address channels by $PnN are unaffected.
#
# Function definitions only, no top-level side effects.
# =============================================================================

# Channels whose $PnS is never touched: scatter and time carry no antibody.
# FSC/SSC may carry a prefix ("SP SSC-A" on a Symphony), so a word boundary
# rather than a line start.
MARKER_SKIP_CHANNEL_REGEX <- "(^|[ _-])(FSC|SSC)([-_ ]|$)|^TIME$"

MARKER_LEDGER_ABSENT <- "<ABSENT>"
MARKER_LEDGER_DEFAULT <- "marker_names_ledger.csv"

.marker_rel <- function(files, fcs_dir) {
  if (is.null(fcs_dir)) return(basename(files))
  root <- normalizePath(fcs_dir, mustWork = FALSE)
  # Resolve the folder only: normalising the file would follow a symlink and
  # give it its target's name.
  full <- file.path(normalizePath(dirname(files), mustWork = FALSE), basename(files))
  pre <- paste0(root, .Platform$file.sep)
  ifelse(startsWith(full, pre), substring(full, nchar(pre) + 1L), basename(files))
}

# -----------------------------------------------------------------------------
# Census
# -----------------------------------------------------------------------------

#' One row per file x parameter: filename, index, channel ($PnN), marker ($PnS
#' or NA when absent or blank). Read-only.
#'
#' @param files Paths to FCS files.
#' @param fcs_dir When given, `filename` is the path relative to it, so files in
#'   subfolders stay distinguishable; otherwise the basename.
marker_census <- function(files, fcs_dir = NULL) {
  rel <- .marker_rel(files, fcs_dir)
  rows <- lapply(seq_along(files), function(i) {
    kw <- tryCatch(fcs_keywords(files[i]), error = function(e) character(0))
    if (!length(kw)) {
      return(data.frame(filename = rel[i], index = NA_integer_,
                        channel = NA_character_, marker = NA_character_,
                        stringsAsFactors = FALSE))
    }
    pt <- fcs_parameter_table(kw)
    pt$marker <- trimws(pt$marker)
    pt$marker[!nzchar(pt$marker)] <- NA_character_
    data.frame(filename = rel[i], index = pt$index, channel = pt$channel,
               marker = pt$marker, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(filename = character(0), index = integer(0),
                      channel = character(0), marker = character(0),
                      stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  out
}

#' Channel x marker spelling x number of files: the cohort's naming at a glance.
census_summary <- function(census) {
  c2 <- census[!is.na(census$index), , drop = FALSE]
  if (!nrow(c2)) {
    return(data.frame(channel = character(0), marker = character(0),
                      n_files = integer(0), stringsAsFactors = FALSE))
  }
  key <- paste(c2$channel, ifelse(is.na(c2$marker), "", c2$marker), sep = "\r")
  parts <- split(c2, key)
  out <- do.call(rbind, lapply(parts, function(d) data.frame(
    channel = d$channel[1], marker = d$marker[1],
    n_files = length(unique(d$filename)), stringsAsFactors = FALSE)))
  out <- out[order(out$channel, -out$n_files, out$marker, na.last = TRUE), ]
  rownames(out) <- NULL
  out
}

# -----------------------------------------------------------------------------
# Planning
# -----------------------------------------------------------------------------

#' Decide the target $PnS for every file and parameter.
#'
#' Precedence, highest first:
#'   alias      an explicit map, as written -> to use (named character vector or
#'              the output of parse_marker_aliases()). Applied to $PnS only.
#'   canonical  the `markers` block of a drift config: named list canonical ->
#'              regex patterns, matched case-insensitively against $PnS first
#'              and then $PnN through assign_canonical(). Fills blank $PnS from
#'              the channel name.
#'   spelling   suggest_marker_aliases() groups with "spelling" confidence
#'              (TNFa / TNF), voted across the cohort. "transposition" groups
#'              (INFy / IFN) are never applied; they come back as suggestions.
#' Scatter and Time channels are never changed.
#'
#' @return list(plan = data.frame(filename, index, channel, marker_old,
#'   marker_new, source), suggestions = data.frame(key, names, suggested,
#'   sources, confidence)) where source is alias / canonical / spelling /
#'   unchanged.
plan_marker_names <- function(census, canonical = NULL, aliases = NULL) {
  plan <- census[!is.na(census$index), , drop = FALSE]
  names(plan)[names(plan) == "marker"] <- "marker_old"
  plan$marker_new <- plan$marker_old
  plan$source <- "unchanged"
  skip <- grepl(MARKER_SKIP_CHANNEL_REGEX, plan$channel, ignore.case = TRUE)
  open <- !skip

  if (!is.null(aliases) && length(aliases)) {
    aliases <- unlist(aliases)
    hit <- open & !is.na(plan$marker_old) & plan$marker_old %in% names(aliases)
    plan$marker_new[hit] <- apply_marker_aliases(plan$marker_old[hit], aliases)
    plan$source[hit] <- "alias"
    open <- open & !hit
  }

  if (!is.null(canonical) && length(canonical)) {
    for (f in unique(plan$filename)) {
      rows <- which(plan$filename == f)
      ptab <- data.frame(index = plan$index[rows], channel = plan$channel[rows],
                         marker = ifelse(is.na(plan$marker_old[rows]), "",
                                         plan$marker_old[rows]),
                         stringsAsFactors = FALSE)
      canon <- assign_canonical(ptab, canonical)$canonical
      hit <- open[rows] & !is.na(canon)
      plan$marker_new[rows][hit] <- canon[hit]
      plan$source[rows][hit] <- "canonical"
      open[rows][hit] <- FALSE
    }
  }

  # Vote over the names as they stand after the explicit rules, so a canonical
  # spelling counts as a source and wins its group.
  named <- !is.na(plan$marker_new)
  panels <- split(plan$marker_new[named], plan$filename[named])
  sugg <- suggest_marker_aliases(panels)
  spelling <- sugg[sugg$confidence == "spelling", , drop = FALSE]
  if (nrow(spelling)) {
    from <- character(0); to <- character(0)
    for (i in seq_len(nrow(spelling))) {
      grp <- strsplit(spelling$names[i], ", ", fixed = TRUE)[[1]]
      grp <- setdiff(grp, spelling$suggested[i])
      from <- c(from, grp); to <- c(to, rep(spelling$suggested[i], length(grp)))
    }
    hit <- open & !is.na(plan$marker_new) & plan$marker_new %in% from
    plan$marker_new[hit] <- unname(to[match(plan$marker_new[hit], from)])
    plan$source[hit] <- "spelling"
  }

  same <- !is.na(plan$marker_new) & !is.na(plan$marker_old) &
    plan$marker_new == plan$marker_old
  plan$source[same] <- "unchanged"
  plan$marker_new[plan$source == "unchanged"] <- plan$marker_old[plan$source == "unchanged"]
  rownames(plan) <- NULL

  suggestions <- sugg[sugg$confidence != "spelling", , drop = FALSE]
  rownames(suggestions) <- NULL
  list(plan = plan[, c("filename", "index", "channel", "marker_old",
                       "marker_new", "source")],
       suggestions = suggestions)
}

#' The plan grouped by (channel, marker_old -> marker_new): one row per
#' distinct rename with the number of files it touches.
summarise_marker_plan <- function(plan) {
  ch <- plan[plan$source != "unchanged", , drop = FALSE]
  if (!nrow(ch)) {
    return(data.frame(channel = character(0), marker_old = character(0),
                      marker_new = character(0), source = character(0),
                      n_files = integer(0), stringsAsFactors = FALSE))
  }
  key <- paste(ch$channel, ifelse(is.na(ch$marker_old), "", ch$marker_old),
               ch$marker_new, sep = "\r")
  out <- do.call(rbind, lapply(split(ch, key), function(d) data.frame(
    channel = d$channel[1], marker_old = d$marker_old[1],
    marker_new = d$marker_new[1], source = d$source[1],
    n_files = length(unique(d$filename)), stringsAsFactors = FALSE)))
  out <- out[order(out$channel, -out$n_files), ]
  rownames(out) <- NULL
  out
}

.blank <- function(x) ifelse(is.na(x), "<blank>", x)

#' The census summary, plan summary and suggestions as text.
format_marker_report <- function(census, plan = NULL, suggestions = NULL) {
  L <- character(0)
  n_files <- length(unique(census$filename))
  cs <- census_summary(census)
  L <- c(L, sprintf("Scanned %d FCS file(s), %d parameter(s) in total.",
                    n_files, sum(!is.na(census$index))))
  n_blank <- sum(!is.na(census$index) & is.na(census$marker))
  L <- c(L, sprintf("  parameters without $PnS: %d", n_blank), "",
         "Marker spellings per channel ($PnN: $PnS x files):")
  for (ch in unique(cs$channel)) {
    d <- cs[cs$channel == ch, ]
    L <- c(L, sprintf("  %-22s %s", ch,
                      paste(sprintf("%s x%d", .blank(d$marker), d$n_files),
                            collapse = ", ")))
  }
  if (!is.null(plan)) {
    sp <- summarise_marker_plan(plan)
    L <- c(L, "", sprintf("Planned $PnS changes: %d parameter(s) in %d file(s).",
                          sum(plan$source != "unchanged"),
                          length(unique(plan$filename[plan$source != "unchanged"]))))
    for (i in seq_len(nrow(sp))) {
      L <- c(L, sprintf("  %-22s %-18s -> %-14s %-10s %d file(s)", sp$channel[i],
                        .blank(sp$marker_old[i]), sp$marker_new[i],
                        sp$source[i], sp$n_files[i]))
    }
    left <- plan[is.na(plan$marker_new) &
                   !grepl(MARKER_SKIP_CHANNEL_REGEX, plan$channel, ignore.case = TRUE), ]
    if (nrow(left)) {
      tab <- table(left$channel)
      L <- c(L, sprintf("  still without $PnS afterwards: %s",
                        paste(sprintf("%s x%d", names(tab), as.integer(tab)),
                              collapse = ", ")))
    }
  }
  if (!is.null(suggestions) && nrow(suggestions)) {
    L <- c(L, "", "Low-confidence groups, NOT applied. Confirm with --aliases if they are the same marker:")
    for (i in seq_len(nrow(suggestions))) {
      L <- c(L, sprintf("  %-14s %s  (suggest %s)", suggestions$confidence[i],
                        suggestions$names[i], suggestions$suggested[i]))
    }
  }
  paste(L, collapse = "\n")
}

# -----------------------------------------------------------------------------
# Ledger
# -----------------------------------------------------------------------------

.marker_ledger_append <- function(ledger_path, run_id, stage, path, keyword = "",
                                  old_value = "", new_value = "", status = "") {
  rows <- data.frame(run_id = run_id, stage = stage, path = path,
                     keyword = keyword, old_value = old_value,
                     new_value = new_value, status = status,
                     stringsAsFactors = FALSE)
  new <- !file.exists(ledger_path)
  utils::write.table(rows, ledger_path, sep = ",", row.names = FALSE,
                     col.names = new, append = !new, qmethod = "double",
                     quote = TRUE, fileEncoding = "UTF-8")
}

.marker_ledger_read <- function(ledger_path) {
  if (!file.exists(ledger_path)) {
    stop("no ledger to revert: ", ledger_path, call. = FALSE)
  }
  utils::read.table(ledger_path, sep = ",", header = TRUE, quote = "\"",
                    stringsAsFactors = FALSE, colClasses = "character",
                    fileEncoding = "UTF-8")
}

# -----------------------------------------------------------------------------
# Apply and revert
# -----------------------------------------------------------------------------

#' Write the planned $PnS values, one file at a time.
#'
#' @param plan The `plan` element of plan_marker_names().
#' @param fcs_dir Folder the plan's filenames are relative to.
#' @param out_dir Write patched copies here (subfolders mirrored) and leave the
#'   originals alone. Unchanged files are copied verbatim so the folder is a
#'   complete cohort. Mutually exclusive with in_place.
#' @param in_place Patch the originals where they sit. Must be said explicitly.
#' @param dry_run TRUE (the default) plans every write and touches nothing.
#' @param ledger_path Undo record, one row per replaced keyword. Default: next
#'   to the written data.
#' @return data.frame(filename, status, strategy, n_changes, message),
#'   invisibly, with the ledger path and run id as attributes.
apply_marker_plan <- function(plan, fcs_dir, out_dir = NULL, in_place = FALSE,
                              dry_run = TRUE, ledger_path = NULL,
                              verbose = TRUE) {
  if (!dir.exists(fcs_dir)) stop("folder does not exist: ", fcs_dir, call. = FALSE)
  fcs_dir <- normalizePath(fcs_dir, mustWork = TRUE)
  if (in_place && !is.null(out_dir)) {
    stop("in_place and out_dir are mutually exclusive.", call. = FALSE)
  }
  if (!dry_run && !in_place && is.null(out_dir)) {
    stop("refusing to guess where to write: pass in_place = TRUE to patch the ",
         "originals, or out_dir to write copies.", call. = FALSE)
  }
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_dir <- normalizePath(out_dir, mustWork = TRUE)
    if (identical(out_dir, fcs_dir) ||
        startsWith(out_dir, paste0(fcs_dir, .Platform$file.sep)) ||
        startsWith(fcs_dir, paste0(out_dir, .Platform$file.sep))) {
      stop("out_dir must not be the same as, inside, or the parent of fcs_dir",
           call. = FALSE)
    }
  }
  if (is.null(ledger_path)) {
    ledger_path <- file.path(if (is.null(out_dir)) fcs_dir else out_dir,
                             MARKER_LEDGER_DEFAULT)
  }
  run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%OS3"), "-", Sys.getpid())

  files <- unique(plan$filename)
  res <- data.frame(filename = files, status = NA_character_,
                    strategy = NA_character_, n_changes = 0L, message = "",
                    stringsAsFactors = FALSE)
  for (i in seq_along(files)) {
    p <- file.path(fcs_dir, files[i])
    dest <- if (is.null(out_dir)) NULL else file.path(out_dir, files[i])
    if (!is.null(dest) && !dir.exists(dirname(dest))) {
      dir.create(dirname(dest), recursive = TRUE)
    }
    rows <- plan[plan$filename == files[i] & plan$source != "unchanged", , drop = FALSE]
    set <- if (nrow(rows)) stats::setNames(rows$marker_new, paste0("$P", rows$index, "S")) else character(0)
    res$n_changes[i] <- length(set)

    link <- Sys.readlink(p)
    if (!is.na(link) && nzchar(link)) {
      res$status[i] <- "refused"; res$message[i] <- "refused: symlink"
      next
    }
    if (!length(set)) {
      res$status[i] <- "unchanged"; res$strategy[i] <- "none"
      if (!dry_run && !is.null(dest)) {
        if (!file.copy(p, dest, overwrite = FALSE)) {
          res$status[i] <- "failed"; res$message[i] <- "copy failed"
        } else res$message[i] <- "copied verbatim"
      }
      next
    }

    if (!dry_run) {
      ins <- tryCatch(fcs_inspect(p), error = function(e) NULL)
      for (k in names(set)) {
        old <- if (is.null(ins) || is.null(ins$text)) NA_character_ else
          fcs_text_get(ins$text, k)
        .marker_ledger_append(ledger_path, run_id, "before", p, k,
                              if (is.na(old)) MARKER_LEDGER_ABSENT else old,
                              set[[k]])
      }
    }
    r <- tryCatch(
      fcs_write_keywords(p, set = set, output_path = dest, dry_run = dry_run,
                         allow_parameter_labels = TRUE),
      error = function(e) list(status = "failed", strategy = NA_character_,
                               message = conditionMessage(e)))
    res$status[i] <- r$status
    res$strategy[i] <- r$strategy
    res$message[i] <- r$message
    if (!dry_run) {
      .marker_ledger_append(ledger_path, run_id, "after", p,
                            status = r$status)
    }
  }

  if (verbose) {
    n_ok <- sum(res$status == "ok"); n_same <- sum(res$status == "unchanged")
    n_bad <- nrow(res) - n_ok - n_same
    cat(sprintf("%d file(s) %s, %d already harmonised, %d not touched.\n",
                n_ok, if (dry_run) "would be written" else "written", n_same, n_bad))
    if (n_bad) {
      bad <- res[!(res$status %in% c("ok", "unchanged")), ]
      cat(sprintf("  %-9s %s  %s\n", bad$status, bad$filename, bad$message), sep = "")
    }
    if (dry_run) cat("Nothing was written. Re-run with apply.\n")
    else cat("Undo record: ", ledger_path, "\n", sep = "")
  }
  attr(res, "ledger_path") <- ledger_path
  attr(res, "run_id") <- run_id
  invisible(res)
}

#' Undo one run from its ledger: the most recent not yet reverted, or `run_id`.
#'
#' Restores every $PnS the run replaced and removes the ones it created. Only
#' files whose "after" row says the write succeeded are touched.
revert_marker_plan <- function(ledger_path, run_id = NULL, verbose = TRUE) {
  led <- .marker_ledger_read(ledger_path)
  if (is.null(run_id)) {
    done <- unique(led$run_id[led$stage == "reverted"])
    cand <- setdiff(unique(led$run_id[led$stage == "after"]), done)
    if (!length(cand)) {
      stop("every run in ", ledger_path, " has already been reverted", call. = FALSE)
    }
    run_id <- sort(cand, decreasing = TRUE)[1]
  }
  led <- led[led$run_id == run_id, , drop = FALSE]
  if (verbose) cat("Reverting run ", run_id, "\n", sep = "")

  paths <- unique(led$path[led$stage == "after" & led$status == "ok"])
  results <- list()
  for (p in paths) {
    before <- led[led$stage == "before" & led$path == p & nzchar(led$keyword), ]
    if (!nrow(before)) next
    before <- before[!duplicated(before$keyword, fromLast = TRUE), , drop = FALSE]
    restore <- stats::setNames(before$old_value, before$keyword)
    unset <- names(restore)[restore == MARKER_LEDGER_ABSENT]
    restore <- restore[restore != MARKER_LEDGER_ABSENT]
    if (!file.exists(p)) {
      results[[p]] <- list(status = "missing", message = "file is gone")
      next
    }
    r <- fcs_write_keywords(p, set = restore, unset = unset,
                            allow_parameter_labels = TRUE)
    results[[p]] <- r
    if (verbose) cat(sprintf("  %-9s %s  %s\n", r$status, basename(p), r$message))
    if (identical(r$status, "ok")) {
      .marker_ledger_append(ledger_path, run_id, "reverted", p, status = "ok")
    }
  }
  invisible(results)
}

# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------

#' Harmonise the $PnS marker names of every FCS file under a folder.
#'
#' @param fcs_dir Folder, searched recursively and case-insensitively.
#' @param canonical `markers` block: named list canonical -> regex patterns.
#' @param aliases Named character vector, as written -> to use.
#' @param probe_only Census only: print the spelling table and stop.
#' @param apply FALSE (default) previews and writes nothing.
#' @param in_place,out_dir Where to write; one must be given with apply.
#' @param revert Undo the last run from the ledger instead.
#' @param ledger Ledger file name, relative to the written data, or a path.
#' @return invisible list(census, plan, suggestions, results, report).
harmonize_marker_names <- function(fcs_dir, canonical = NULL, aliases = NULL,
                                   probe_only = FALSE, apply = FALSE,
                                   in_place = FALSE, out_dir = NULL,
                                   revert = FALSE, ledger = MARKER_LEDGER_DEFAULT,
                                   recursive = TRUE, verbose = TRUE) {
  if (!dir.exists(fcs_dir)) stop("Folder does not exist: ", fcs_dir, call. = FALSE)
  fcs_dir <- normalizePath(fcs_dir, mustWork = TRUE)
  ledger_dir <- if (is.null(out_dir)) fcs_dir else out_dir
  ledger_path <- if (basename(ledger) == ledger) file.path(ledger_dir, ledger) else ledger

  if (revert) {
    res <- revert_marker_plan(ledger_path, verbose = verbose)
    return(invisible(list(results = res)))
  }
  if (apply && !in_place && is.null(out_dir)) {
    stop("refusing to guess where to write: pass in_place = TRUE to patch the ",
         "originals, or out_dir to write copies.", call. = FALSE)
  }

  files <- list.files(fcs_dir, pattern = "\\.fcs$", recursive = recursive,
                      full.names = TRUE, ignore.case = TRUE)
  files <- files[!dir.exists(files)]
  if (!is.null(out_dir) && dir.exists(out_dir)) {
    od <- normalizePath(out_dir)
    files <- files[!startsWith(normalizePath(files), paste0(od, .Platform$file.sep))]
  }
  if (!length(files)) stop("No .fcs files found under ", fcs_dir, call. = FALSE)
  if (verbose) cat("Reading ", length(files), " FCS file(s) under ", fcs_dir, " ...\n", sep = "")

  census <- marker_census(sort(files), fcs_dir)
  out <- list(census = census, plan = NULL, suggestions = NULL, results = NULL,
              report = "", ledger = ledger_path)
  if (probe_only) {
    out$report <- format_marker_report(census)
    if (verbose) cat(out$report, "\n", sep = "")
    return(invisible(out))
  }

  pl <- plan_marker_names(census, canonical = canonical, aliases = aliases)
  out$plan <- pl$plan
  out$suggestions <- pl$suggestions
  out$report <- format_marker_report(census, pl$plan, pl$suggestions)
  if (verbose) {
    cat(out$report, "\n\n", if (apply) "Applying" else "PREVIEW, nothing will be written",
        ":\n", sep = "")
  }
  out$results <- apply_marker_plan(pl$plan, fcs_dir, out_dir = out_dir,
                                   in_place = in_place, dry_run = !apply,
                                   ledger_path = ledger_path, verbose = verbose)
  n_bad <- sum(!(out$results$status %in% c("ok", "unchanged")))
  if (apply && n_bad > 0) {
    stop(n_bad, " of ", nrow(out$results), " file(s) could not be written. ",
         "See the report above.", call. = FALSE)
  }
  invisible(out)
}
