#!/usr/bin/env Rscript
# =============================================================================
# metadata_qc/sync_fcs_names.R — make the name inside an FCS file match the
# name on disk
# =============================================================================
# The problem this solves: a plate loader or a Diva tube template acquires a
# whole run under one file name. The analyst renames the files afterwards to
# something meaningful, but the name is also recorded INSIDE each file, in the
# $FIL keyword, and that is what FlowJo displays by default (Preferences ->
# Workspace -> "How should data files be named" -> "FCS Default"). So the
# workspace shows n identically-labelled rows. FlowJo cannot repair it — it
# never writes back to raw data — and re-exporting re-encodes every event.
#
# This module walks a folder of ALREADY-RENAMED files and writes each file's
# own name into its metadata, using src/R_shared/fcs_keywords.R so the event
# data is never decoded. Base R only: no flowCore, no Bioconductor, no optparse.
#
# DEPENDENCY: source R_shared/fcs_keywords.R first. io_utils.R is optional
# (die/make_outdir); local fallbacks are used when it is absent.
#
# THREE THINGS THIS DELIBERATELY DOES NOT DO
#   1. It does not write by default. The default mode is --probe, a read-only
#      census, because "$FIL collides" is an inference about FlowJo's behaviour
#      and the census is how you check it before mutating a cohort.
#   2. It does not rewrite anything but $FIL unless asked. GUID, $SRC, TUBE
#      NAME and FILENAME are all opt-in — see the keyword notes on
#      plan_name_sync().
#   3. It does not stamp provenance into the files. A "renamed by / renamed on"
#      pair costs ~120 bytes, which forces the slow rewrite path on every
#      zero-slack file, re-injects the colliding value into a field a LIMS will
#      index, and makes a second run non-idempotent. The ledger CSV is the
#      provenance record, and it is what --revert replays.
#
# This file holds only function definitions (no top-level side effects).
# =============================================================================

# --- local fallbacks, so the module stands alone ------------------------------
if (!exists("make_outdir", mode = "function")) {
  make_outdir <- function(path) if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# Keywords whose value is a name and which some site somewhere displays. Not
# rewritten by default — this is the menu for --also, and the list the census
# scans when it looks for other copies of the stale name.
FCS_NAME_LIKE_KEYWORDS <- c("$FIL", "FILENAME", "GUID", "$SRC", "TUBE NAME",
                            "$TUBE NAME", "$SMNO", "WELL ID", "PLATE ID",
                            "PLATE NAME", "$PLATEID", "$WELLID", "EXPORT NAME")

# Where the name the file was ACQUIRED under is preserved when $FIL is rewritten.
#
# No "$" prefix on purpose: FCS 3.1 3.2.12 reserves that namespace for keywords
# the standard itself defines, so "$ORIGINALFIL" would be a non-conforming
# invention that a strict reader is entitled to reject. This is the same
# convention FlowJo uses for its own export provenance (ORIGINALGUID).
FCS_ORIGINAL_NAME_KEYWORD <- "ORIGINALFIL"

# The keywords a cohort is usually compared on. Ordering for the keyword table,
# and the subset shown before "show every keyword" is switched on.
FCS_IDENTITY_KEYWORDS <- c("$FIL", FCS_ORIGINAL_NAME_KEYWORD, "GUID", "$SRC",
                           "TUBE NAME", "$SMNO", "WELL ID", "FILENAME",
                           "$TOT", "$PAR", "$DATE", "$BTIM", "$CYT", "$CYTSN")

# -----------------------------------------------------------------------------
# Discovery
# -----------------------------------------------------------------------------

#' Every .fcs file under `dir`.
#'
#' Case-insensitive, because Windows-sourced cohorts routinely carry .FCS.
#' Symlinks are listed but the writer refuses them, so they are reported rather
#' than silently skipped.
scan_fcs_dir <- function(dir, recursive = TRUE, exclude = character()) {
  if (!dir.exists(dir)) stop("folder does not exist: ", dir, call. = FALSE)
  files <- list.files(dir, pattern = "\\.fcs$", recursive = recursive,
                      full.names = TRUE, ignore.case = TRUE)
  files <- files[!dir.exists(files)]
  if (length(exclude) && length(files)) {
    ex <- normalizePath(exclude, mustWork = FALSE)
    keep <- vapply(normalizePath(files, mustWork = FALSE), function(f)
      !any(startsWith(f, paste0(ex, .Platform$file.sep))) && !(f %in% ex),
      logical(1))
    files <- files[keep]
  }
  sort(files)
}

# -----------------------------------------------------------------------------
# Census
# -----------------------------------------------------------------------------

#' Read every file's layout and identity keywords. Writes nothing.
#'
#' @return data.frame, one row per file.
probe_fcs_dir <- function(files, fcs_dir = NULL, allow = character()) {
  n <- length(files)
  blank <- function(v) rep(v, n)
  out <- data.frame(
    path = files, name = basename(files),
    rel = if (is.null(fcs_dir)) basename(files) else
      sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", fcs_dir),
                 .Platform$file.sep, "?"), "", files),
    tot = blank(NA_character_),
    fil = blank(NA_character_), filename_kw = blank(NA_character_),
    guid = blank(NA_character_), src = blank(NA_character_),
    tube = blank(NA_character_), originalfil = blank(NA_character_),
    strategy = blank(NA_character_),
    refusals = blank(""), warnings = blank(""),
    stringsAsFactors = FALSE)

  for (i in seq_len(n)) {
    ins <- tryCatch(fcs_inspect(files[i], allow = allow),
                    error = function(e) list(refusals = "unreadable_file",
                                             warnings = character(0)))
    out$refusals[i] <- paste(ins$refusals, collapse = ",")
    out$warnings[i] <- paste(ins$warnings, collapse = ",")
    if (is.null(ins$header)) next
    if (!is.null(ins$text)) {
      out$fil[i] <- fcs_text_get(ins$text, "$FIL")
      out$filename_kw[i] <- fcs_text_get(ins$text, "FILENAME")
      out$guid[i] <- fcs_text_get(ins$text, "GUID")
      out$src[i] <- fcs_text_get(ins$text, "$SRC")
      out$tube[i] <- fcs_text_get(ins$text, "TUBE NAME")
      out$originalfil[i] <- fcs_text_get(ins$text, FCS_ORIGINAL_NAME_KEYWORD)
      out$tot[i] <- fcs_text_get(ins$text, "$TOT")
    }
  }
  out
}

#' Groups of files sharing a value. The diagnosis half of the census.
#' @return named list, value -> file names, only for values seen more than once.
collision_groups <- function(values, names_) {
  keep <- !is.na(values) & nzchar(values)
  if (!any(keep)) return(list())
  v <- values[keep]; nm <- names_[keep]
  tab <- split(nm, v)
  tab[vapply(tab, length, integer(1)) > 1L]
}

#' Every keyword of every file, as a long-format CSV.
#'
#' The TEXT segment is a delimiter-packed blob that nothing outside a parser can
#' read, so "what does this cohort actually say about itself" is a question you
#' cannot answer by looking at the files. This answers it for the whole folder
#' at once.
#'
#' Long format (one row per file/keyword) rather than a wide grid because files
#' do not have to carry the same keywords -- a wide table would need a union of
#' columns and a rule for the holes, and the caller that wants a grid can pivot
#' far more easily than one that wants tidy rows can unpivot.
#'
#' `keyword_rank` lets a consumer sort identity keywords to the top without
#' having to know which those are: 1..n for FCS_IDENTITY_KEYWORDS in order,
#' then 999 for everything else. `name_like` says whether a value shared by two
#' files is a PROBLEM: a cohort sharing $TOT or $DATE is ordinary, a cohort
#' sharing $FIL or GUID is the thing this tool exists to fix. Both are emitted
#' so a consumer never has to hard-code either list.
#'
#' Read-only. Writes exactly one file: `path`.
#'
#' @return the data.frame it wrote, invisibly.
write_keyword_table <- function(files, path, fcs_dir = NULL, allow = character()) {
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    f <- files[i]
    rel <- if (is.null(fcs_dir)) basename(f) else
      sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", fcs_dir),
                 .Platform$file.sep, "?"), "", f)
    kw <- tryCatch(fcs_keywords(f), error = function(e) character(0))
    if (!length(kw)) {
      rows[[i]] <- data.frame(file = rel, keyword = "(unreadable)", value = "",
                              keyword_rank = 999L, name_like = FALSE,
                              stringsAsFactors = FALSE)
      next
    }
    keys <- names(kw)
    upper <- toupper(trimws(keys))
    rank <- match(upper, toupper(FCS_IDENTITY_KEYWORDS))
    rank[is.na(rank)] <- 999L
    rows[[i]] <- data.frame(file = rel, keyword = keys, value = unname(kw),
                            keyword_rank = as.integer(rank),
                            # ORIGINALFIL is deliberately NOT in here. It is
                            # SUPPOSED to be shared -- it records the one name
                            # the whole cohort used to have -- so flagging it
                            # would report the fix as if it were the fault.
                            name_like = upper %in% toupper(FCS_NAME_LIKE_KEYWORDS),
                            stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(file = character(0), keyword = character(0),
                      value = character(0), keyword_rank = integer(0),
                      name_like = logical(0), stringsAsFactors = FALSE)
  }
  make_outdir(dirname(path))
  utils::write.csv(out, path, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(out)
}

# -----------------------------------------------------------------------------
# Planning
# -----------------------------------------------------------------------------

.expand_template <- function(tmpl, path, fcs_dir, old_name) {
  nm <- basename(path)
  rel <- if (is.null(fcs_dir)) nm else
    sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", fcs_dir),
               .Platform$file.sep, "?"), "", path)
  out <- tmpl
  out <- gsub("%name%", nm, out, fixed = TRUE)
  out <- gsub("%stem%", sub("\\.[Ff][Cc][Ss]$", "", nm), out, fixed = TRUE)
  out <- gsub("%oldname%", if (is.na(old_name)) "" else old_name, out, fixed = TRUE)
  out <- gsub("%relpath%", rel, out, fixed = TRUE)
  out <- gsub("%dir%", basename(dirname(path)), out, fixed = TRUE)
  out
}

#' Decide, per file, exactly which keywords get which values.
#'
#' Keyword policy, and why it is this narrow:
#'   $FIL      always. It is FlowJo's display name under the shipped preference,
#'             and the whole point of the exercise. Created when absent.
#'   FILENAME  opt-in. A flowCore::read.FCS artefact (IO.R:360) recording the
#'             absolute path it read; FlowJo does not use it and it is not in
#'             FlowJo's keyword list. Writing a bare basename there would
#'             quietly redefine the field for no benefit.
#'   GUID      "never" by default. CytoML's sample matcher never reads it and
#'             nothing in this repo derives a sample name from it. "auto" only
#'             rewrites it when its current value is already just an echo of the
#'             old file name — which is the flowCore round-trip case, and the
#'             case in this repo's own fixtures. Never mint a fresh UUID: that
#'             destroys the only cross-tool linkage the file has, and GUID is
#'             not what is broken.
#'   ORIGINALFIL  written WHENEVER $FIL is being replaced and the keyword is not
#'             already there, so the name the file was acquired under survives
#'             the rename and can be read back. Write-once: a second run leaves
#'             it alone, which is both what "original" should mean and what
#'             keeps repeat runs idempotent. No "$" prefix — that namespace
#'             belongs to the standard (3.2.12). Disable with write_original.
#'   others    opt-in via --also, for sites whose FlowJo is set to "Keyword
#'             Combination" over $SRC / TUBE NAME.
#'
#' @return data.frame with one row per file and a `sets` list-column.
plan_name_sync <- function(probe, fcs_dir = NULL, also = character(),
                           set = character(), unset = character(),
                           guid = c("never", "auto", "always"),
                           strip_extension = FALSE,
                           write_original = TRUE) {
  guid <- match.arg(guid)
  n <- nrow(probe)
  plan <- probe
  plan$new_fil <- NA_character_
  plan$sets <- vector("list", n)
  plan$unsets <- vector("list", n)

  for (i in seq_len(n)) {
    nm <- probe$name[i]
    value <- if (strip_extension) sub("\\.[Ff][Cc][Ss]$", "", nm) else nm
    sets <- c("$FIL" = value)

    for (k in also) {
      k <- trimws(k)
      if (!nzchar(k)) next
      sets[[k]] <- value
    }

    old_guid <- probe$guid[i]
    write_guid <- switch(guid,
      never = FALSE,
      always = TRUE,
      auto = {
        # "Name-shaped" means the GUID is already just a copy of a file name,
        # so replacing it loses nothing. A real vendor GUID is left alone.
        cand <- c(probe$fil[i], probe$filename_kw[i], probe$name[i])
        cand <- c(cand, basename(cand[!is.na(cand)]),
                  sub("\\.[Ff][Cc][Ss]$", "", cand[!is.na(cand)]))
        !is.na(old_guid) && nzchar(old_guid) &&
          tolower(old_guid) %in% tolower(cand[!is.na(cand)])
      })
    if (isTRUE(write_guid)) sets[["GUID"]] <- value

    # Preserve the acquisition name, once. Write-only-if-absent is what makes
    # "original" keep meaning original: a second run must not overwrite it with
    # the name the FIRST run wrote, and skipping it keeps repeat runs no-ops.
    # There is nothing to preserve when the file never carried a $FIL, so
    # nothing is written in that case rather than recording the new name as if
    # it were the old one.
    if (write_original && is.na(probe$originalfil[i]) &&
        !is.na(probe$fil[i]) && nzchar(probe$fil[i]) &&
        !identical(probe$fil[i], value)) {
      sets[[FCS_ORIGINAL_NAME_KEYWORD]] <- probe$fil[i]
    }

    for (nmk in names(set)) {
      sets[[nmk]] <- .expand_template(set[[nmk]], probe$path[i], fcs_dir,
                                      probe$fil[i])
    }

    plan$new_fil[i] <- value
    plan$sets[[i]] <- sets
    plan$unsets[[i]] <- unset
  }
  plan
}

# Case- and accent-insensitive enough to catch "A.fcs" vs "a.fcs" on a
# case-insensitive filesystem. iconv to ASCII is approximate (stringi is not
# available here), so a non-ASCII name is reported for a human to look at
# rather than silently folded.
.fold_name <- function(x) {
  y <- toupper(x)
  z <- iconv(y, "", "ASCII//TRANSLIT", sub = "?")
  ifelse(is.na(z), y, z)
}

#' Names that would still collide AFTER the sync. Fatal, by design.
#'
#' Two files in different subfolders can share a basename, and syncing both
#' would recreate exactly the duplicate the user came here to remove. There is
#' no --force for this: the only fix is to rename one of the files, which the
#' user is better placed to do than we are.
detect_collisions <- function(plan) {
  exact <- collision_groups(plan$new_fil, plan$rel)
  folded <- collision_groups(.fold_name(plan$new_fil), plan$rel)
  near <- folded[!(names(folded) %in% .fold_name(names(exact)))]
  nonascii <- plan$rel[grepl("[^ -~]", plan$new_fil)]
  list(exact = exact, near = near, non_ascii = nonascii)
}

# -----------------------------------------------------------------------------
# Ledger — the undo record
# -----------------------------------------------------------------------------

FCS_LEDGER_ABSENT <- "<ABSENT>"

.ledger_write <- function(rows, ledger_path) {
  new <- !file.exists(ledger_path)
  utils::write.table(rows, ledger_path, sep = ";", row.names = FALSE,
                     col.names = new, append = !new, qmethod = "double",
                     quote = TRUE, fileEncoding = "UTF-8")
}

#' Append rows for one file and flush them to disk immediately.
#'
#' Flushed BEFORE the file is opened for writing, so a crash mid-write still
#' leaves a record of what the previous value was.
append_ledger <- function(ledger_path, run_id, stage, path, strategy = "",
                          status = "", changes = list()) {
  if (!length(changes)) {
    changes <- list(list(keyword = "", old = "", new = ""))
  }
  rows <- do.call(rbind, lapply(changes, function(ch) data.frame(
    run_id = run_id, stage = stage, path = path,
    keyword = ch$keyword, old_value = ch$old, new_value = ch$new,
    strategy = strategy, status = status, stringsAsFactors = FALSE)))
  .ledger_write(rows, ledger_path)
  invisible(rows)
}

#' Undo a previous run from its ledger.
#'
#' Reverts ONE run: by default the most recent one that has not been reverted
#' already. That matters because a ledger accumulates every run against a
#' folder, and replaying all of them at once would restore each keyword to
#' whatever the oldest entry happened to say — which for a folder synced twice
#' is not any state the files were ever in. Calling --revert repeatedly walks
#' back through the history one run at a time.
#'
#' Only files whose "after" row says the write succeeded are touched, and a
#' keyword that did not exist before is removed rather than blanked.
revert_from_ledger <- function(ledger_path, run_id = NULL, verbose = TRUE) {
  if (!file.exists(ledger_path)) {
    stop("no ledger to revert: ", ledger_path, call. = FALSE)
  }
  led <- utils::read.table(ledger_path, sep = ";", header = TRUE,
                           quote = "\"", stringsAsFactors = FALSE,
                           colClasses = "character", fileEncoding = "UTF-8")
  if (is.null(run_id)) {
    already <- unique(led$run_id[led$stage == "reverted"])
    candidates <- setdiff(unique(led$run_id[led$stage == "after"]), already)
    if (!length(candidates)) {
      stop("every run in ", ledger_path, " has already been reverted",
           call. = FALSE)
    }
    run_id <- sort(candidates, decreasing = TRUE)[1]
  }
  led <- led[led$run_id == run_id, , drop = FALSE]
  if (verbose) cat("Reverting run ", run_id, "\n", sep = "")

  done <- unique(led$path[led$stage == "after" & led$status == "ok"])
  results <- list()
  for (p in done) {
    before <- led[led$stage == "before" & led$path == p & nzchar(led$keyword), ]
    if (!nrow(before)) next
    before <- before[!duplicated(before$keyword, fromLast = TRUE), , drop = FALSE]
    restore <- before$old_value
    names(restore) <- before$keyword
    unset <- names(restore)[restore == FCS_LEDGER_ABSENT]
    restore <- restore[restore != FCS_LEDGER_ABSENT]
    if (!file.exists(p)) {
      results[[p]] <- list(path = p, status = "missing", strategy = NA_character_,
                           changed = character(0), message = "file is gone")
      next
    }
    r <- fcs_write_keywords(p, set = restore, unset = unset)
    results[[p]] <- r
    if (verbose) {
      cat(sprintf("  %-9s %s  %s\n", r$status, basename(p), r$message))
    }
  }
  # Record the revert, so a second --revert steps to the run before this one
  # rather than undoing the same thing twice.
  if (length(results)) {
    ok_paths <- names(results)[vapply(results, function(r)
      identical(r$status, "ok"), logical(1))]
    if (length(ok_paths)) {
      for (p in ok_paths) {
        append_ledger(ledger_path, run_id, "reverted", p, status = "ok")
      }
    }
  }
  invisible(results)
}

# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------

.fmt_group <- function(groups, label, max_show = 8L) {
  if (!length(groups)) return(character(0))
  out <- sprintf("  %s: %d group(s)", label, length(groups))
  for (nm in names(groups)[seq_len(min(length(groups), max_show))] ) {
    out <- c(out, sprintf("    \"%s\"  <- %s", nm,
                          paste(groups[[nm]], collapse = ", ")))
  }
  if (length(groups) > max_show) {
    out <- c(out, sprintf("    ... and %d more", length(groups) - max_show))
  }
  out
}

#' The census, as text.
format_probe_report <- function(probe, plan = NULL, collisions = NULL,
                                workspaces = character()) {
  L <- character(0)
  n <- nrow(probe)
  L <- c(L, sprintf("Scanned %d FCS file(s).", n))
  if (n == 0L) return(paste(L, collapse = "\n"))

  bad <- probe[nzchar(probe$refusals), , drop = FALSE]
  L <- c(L, sprintf("  readable and safe to edit: %d", n - nrow(bad)))
  if (nrow(bad)) {
    L <- c(L, sprintf("  refused: %d", nrow(bad)))
    tab <- sort(table(unlist(strsplit(bad$refusals, ","))), decreasing = TRUE)
    for (nm in names(tab)) L <- c(L, sprintf("    %-32s %d", nm, tab[[nm]]))
  }
  warn <- probe[nzchar(probe$warnings), , drop = FALSE]
  if (nrow(warn)) {
    tab <- sort(table(unlist(strsplit(warn$warnings, ","))), decreasing = TRUE)
    for (nm in names(tab)) L <- c(L, sprintf("  note: %-30s %d", nm, tab[[nm]]))
  }
  if (!is.null(plan) && "strategy" %in% names(plan)) {
    st <- table(plan$strategy[!is.na(plan$strategy)])
    if (length(st)) {
      L <- c(L, "  planned write strategy:")
      for (nm in names(st)) L <- c(L, sprintf("    %-32s %d", nm, st[[nm]]))
    }
  }

  L <- c(L, "", "Values shared by more than one file (before the sync):")
  found <- FALSE
  for (col in list(c("fil", "$FIL"), c("guid", "GUID"), c("src", "$SRC"),
                   c("tube", "TUBE NAME"), c("name", "file name on disk"))) {
    g <- collision_groups(probe[[col[1]]], probe$rel)
    if (length(g)) { found <- TRUE; L <- c(L, .fmt_group(g, col[2])) }
  }
  if (!found) L <- c(L, "  none — no identity keyword is duplicated.")

  # The check that matters more than the naming. CytoML resolves a workspace
  # sample to a file by $FIL when the on-disk name does not match, then
  # tie-breaks on event count; this repo passes greedy_match = TRUE
  # (flowjo_workspace.R:819), which takes the first $TOT-correct hit instead of
  # erroring on ambiguity. Two files sharing BOTH values are being told apart
  # by luck, and may already have been matched to the wrong gates.
  # A control byte the TEXT segment cannot contain (3.2.15 bars ASCII 0-31
  # outside a value), so the join is unambiguous and the split below is exact.
  sep <- rawToChar(as.raw(1))
  key <- paste(probe$fil, probe$tot, sep = sep)
  key[is.na(probe$fil) | is.na(probe$tot)] <- NA_character_
  lucky <- collision_groups(key, probe$rel)
  L <- c(L, "")
  if (length(lucky)) {
    L <- c(L, sprintf(paste0("WARNING: %d group(s) of files share BOTH $FIL and ",
                             "$TOT."), length(lucky)))
    L <- c(L, paste0("  CytoML cannot tell these apart, and with greedy_match ",
                     "it picks the first."))
    L <- c(L, paste0("  Any gating already run against these files may have ",
                     "used the wrong sample."))
    for (nm in names(lucky)[seq_len(min(length(lucky), 8L))]) {
      parts <- strsplit(nm, sep, fixed = TRUE)[[1]]
      L <- c(L, sprintf("    $FIL \"%s\", $TOT %s  <- %s", parts[1], parts[2],
                        paste(lucky[[nm]], collapse = ", ")))
    }
  } else {
    L <- c(L, "($FIL, $TOT) is unique for every file — nothing is being matched by luck.")
  }

  if (!is.null(collisions)) {
    if (length(collisions$exact)) {
      L <- c(L, "", "FATAL: these files would STILL collide after the sync:")
      L <- c(L, .fmt_group(collisions$exact, "new $FIL"))
      L <- c(L, "  Rename one of each pair on disk first.")
    }
    if (length(collisions$near)) {
      L <- c(L, "", "Near-collisions (differ only by case or accent):")
      L <- c(L, .fmt_group(collisions$near, "folded name"))
    }
    if (length(collisions$non_ascii)) {
      L <- c(L, "", paste0("Non-ASCII file names, worth a look by eye: ",
                           paste(collisions$non_ascii, collapse = ", ")))
    }
  }

  if (length(workspaces)) {
    L <- c(L, "", strrep("-", 74),
      sprintf("%d FlowJo workspace(s) sit under this folder:", length(workspaces)),
      paste0("  ", paste(workspaces, collapse = "\n  ")),
      "",
      "Renaming the files on disk already broke those workspaces: this repo's",
      "resolve_fcs_names() (src/loading/flowjo_workspace.R:192-234) pairs a",
      "workspace sample to a file BY ON-DISK BASENAME. Syncing $FIL also removes",
      "CytoML's fallback of matching a sample to a file whose $FIL still carries",
      "the old label. If you intend to keep using those workspaces rather than",
      "re-gating, re-open each one in FlowJo, re-point its samples and re-save",
      "AFTER this run, so FlowJo refreshes its own frozen copy of $FIL.",
      strrep("-", 74))
  }
  paste(L, collapse = "\n")
}

#' The per-file change list, as text.
format_sync_report <- function(results) {
  if (!length(results)) return("")
  L <- character(0)
  for (r in results) {
    L <- c(L, sprintf("  %-9s %-9s %s", r$status,
                      if (is.na(r$strategy)) "" else r$strategy,
                      basename(r$path)))
    if (length(r$changed)) {
      for (nm in names(r$changed)) {
        L <- c(L, sprintf("              %s = %s", nm, r$changed[[nm]]))
      }
    }
    if (nzchar(r$message) && r$status != "ok") {
      L <- c(L, sprintf("              %s", r$message))
    }
  }
  paste(L, collapse = "\n")
}

# -----------------------------------------------------------------------------
# The entry point
# -----------------------------------------------------------------------------

#' Sync the embedded name of every FCS file in a folder to its file name.
#'
#' @param fcs_dir Folder, searched recursively by default.
#' @param apply FALSE (the default) reports and writes nothing.
#' @param in_place Patch the originals. Mutually exclusive with output_dir.
#' @param output_dir Write patched copies here instead, leaving originals alone.
#' @param also Extra keywords to give the same value, e.g. c("$SRC","TUBE NAME").
#' @param set Named character vector of KEY = template. Templates understand
#'   %name%, %stem%, %oldname%, %relpath% and %dir%.
#' @param unset Keywords to remove.
#' @param guid "never" (default), "auto" or "always".
#' @param strip_extension Write $FIL without the .fcs suffix.
#' @param probe_only Census only.
#' @param revert Undo a previous run from its ledger.
#' @param allow Refusal codes to downgrade.
#' @param ledger Ledger file name, relative to the output/ledger directory.
#' @param keyword_table Path to write every keyword of every file to, as a
#'   long-format CSV. Read-only, and done before any write, so it also works
#'   with probe_only.
#' @param write_original Preserve the replaced $FIL in ORIGINALFIL (write-once).
#' @return invisible list(probe, plan, collisions, results, report).
sync_fcs_names <- function(fcs_dir, apply = FALSE, in_place = TRUE,
                           output_dir = NULL, also = character(),
                           set = character(), unset = character(),
                           guid = c("never", "auto", "always"),
                           strip_extension = FALSE, recursive = TRUE,
                           probe_only = FALSE,
                           revert = FALSE, allow = character(),
                           ledger = "fcs_name_sync_ledger.csv",
                           keyword_table = NULL, write_original = TRUE,
                           verbose = TRUE) {
  guid <- match.arg(guid)
  if (!dir.exists(fcs_dir)) {
    stop("Folder does not exist: ", fcs_dir, call. = FALSE)
  }
  fcs_dir <- normalizePath(fcs_dir, mustWork = TRUE)

  ledger_dir <- if (is.null(output_dir)) fcs_dir else output_dir
  ledger_path <- if (basename(ledger) == ledger) file.path(ledger_dir, ledger)
                 else ledger

  if (revert) {
    if (verbose) cat("Reverting from ", ledger_path, "\n", sep = "")
    res <- revert_from_ledger(ledger_path, verbose = verbose)
    return(invisible(list(results = res, report = "")))
  }

  # --- destination sanity, before anything is read --------------------------
  if (!is.null(output_dir)) {
    make_outdir(output_dir)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    # file.copy(x, x) truncates, and a destination nested inside the source
    # would be re-scanned on the next run. Both are fatal, not warnings.
    if (identical(output_dir, fcs_dir) ||
        startsWith(output_dir, paste0(fcs_dir, .Platform$file.sep)) ||
        startsWith(fcs_dir, paste0(output_dir, .Platform$file.sep))) {
      stop("--output_dir must not be the same as, inside, or the parent of ",
           "--fcs_dir (got ", output_dir, ")", call. = FALSE)
    }
    in_place <- FALSE
  }

  files <- scan_fcs_dir(fcs_dir, recursive = recursive,
                        exclude = if (is.null(output_dir)) character() else output_dir)
  if (!length(files)) {
    stop("No .fcs files found under ", fcs_dir, call. = FALSE)
  }

  # Leftovers from a run that died mid-write. Naming them is more useful than
  # cleaning them up silently, because one of them may be the only copy.
  orphans <- list.files(c(fcs_dir, if (!is.null(output_dir)) output_dir),
                        pattern = "\\.fcs(prefix|tmp[0-9]*|old[0-9]*)$",
                        recursive = recursive, full.names = TRUE)
  if (length(orphans)) {
    stop("a previous run left these behind — inspect them before re-running:\n  ",
         paste(orphans, collapse = "\n  "), call. = FALSE)
  }

  if (verbose) cat("Reading ", length(files), " FCS file(s) under ", fcs_dir,
                   " ...\n", sep = "")
  probe <- probe_fcs_dir(files, fcs_dir = fcs_dir, allow = allow)
  plan <- plan_name_sync(probe, fcs_dir = fcs_dir, also = also, set = set,
                         unset = unset, guid = guid,
                         strip_extension = strip_extension,
                         write_original = write_original)

  if (!is.null(keyword_table) && nzchar(keyword_table)) {
    kwt <- write_keyword_table(files, keyword_table, fcs_dir = fcs_dir,
                               allow = allow)
    if (verbose) {
      cat("Keyword table: ", nrow(kwt), " rows -> ", keyword_table, "\n", sep = "")
    }
  }

  collisions <- detect_collisions(plan)
  workspaces <- list.files(fcs_dir, pattern = "\\.(wsp|jo)$", recursive = TRUE,
                           ignore.case = TRUE)
  plan$strategy <- NA_character_

  out <- list(probe = probe, plan = plan, collisions = collisions,
              workspaces = workspaces, results = list(), report = "",
              ledger = ledger_path)

  # A census stops here: no per-file simulation, which is what keeps it cheap
  # enough to run just for the keyword table.
  if (probe_only) {
    out$report <- format_probe_report(probe, plan, collisions, workspaces)
    if (verbose) cat(out$report, "\n", sep = "")
    return(invisible(out))
  }

  if (apply && length(collisions$exact)) {
    cat(format_probe_report(probe, plan, collisions, workspaces), "\n",
        sep = "", file = stderr())
    stop(length(collisions$exact), " group(s) of files would still share a name ",
         "after the sync. Nothing was written. See the report above.",
         call. = FALSE)
  }

  run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%OS3"), "-", Sys.getpid())
  results <- vector("list", nrow(plan))

  # One pass: the strategies come out of the results, and the report is
  # printed after the pass rather than before it.
  for (i in seq_len(nrow(plan))) {
    p <- plan$path[i]
    dest <- if (is.null(output_dir)) NULL else
      file.path(output_dir, plan$rel[i])
    if (!is.null(dest)) make_outdir(dirname(dest))

    if (apply) {
      # Flushed before the file is opened, so a crash still leaves the record.
      ins <- tryCatch(fcs_inspect(p, allow = allow), error = function(e) NULL)
      changes <- lapply(names(plan$sets[[i]]), function(k) {
        old <- if (is.null(ins) || is.null(ins$text)) NA_character_ else
          fcs_text_get(ins$text, k)
        list(keyword = k, new = plan$sets[[i]][[k]],
             old = if (is.na(old)) FCS_LEDGER_ABSENT else old)
      })
      append_ledger(ledger_path, run_id, "before", p, changes = changes)
    }

    r <- tryCatch(
      fcs_write_keywords(p, set = plan$sets[[i]], unset = plan$unsets[[i]],
                         output_path = dest, allow = allow, dry_run = !apply),
      error = function(e) list(path = p, status = "failed", strategy = NA_character_,
                               changed = character(0), refusals = character(0),
                               message = conditionMessage(e)))
    r$path <- p
    results[[i]] <- r
    plan$strategy[i] <- if (identical(r$status, "ok")) r$strategy else r$status

    if (apply) {
      append_ledger(ledger_path, run_id, "after", p, strategy = r$strategy,
                    status = r$status)
    }
  }

  out$plan <- plan
  out$report <- format_probe_report(probe, plan, collisions, workspaces)
  if (verbose) {
    cat(out$report, "\n", sep = "")
    cat("\n", if (apply) "Applying" else "PREVIEW — nothing will be written",
        ":\n", sep = "")
    for (i in seq_len(nrow(plan))) {
      cat(sprintf("  %-9s %-9s %s\n", results[[i]]$status,
                  if (is.na(results[[i]]$strategy)) "" else results[[i]]$strategy,
                  plan$rel[i]))
    }
  }
  out$results <- results

  status <- vapply(results, function(r) r$status, character(1))
  n_ok <- sum(status == "ok")
  n_same <- sum(status == "unchanged")
  n_bad <- sum(!(status %in% c("ok", "unchanged")))
  if (verbose) {
    cat("\n", sprintf("%d %s, %d already in sync, %d not touched.",
                      n_ok, if (apply) "written" else "would be written",
                      n_same, n_bad), "\n", sep = "")
    if (apply) cat("Undo record: ", ledger_path, "\n", sep = "")
    if (!apply) cat("Nothing was written. Re-run with --apply.\n")
  }

  # Gated on `apply`: a preview that reports unsyncable files has done its job,
  # so only a write fails.
  if (n_bad > 0) {
    # The detail goes to stderr BEFORE stop(), because R truncates the message
    # it prints from stop() at warning.length and would swallow most of a
    # cohort-sized report.
    detail <- format_sync_report(results[!(status %in% c("ok", "unchanged"))])
    cat("\nFiles that were not touched:\n", detail, "\n", sep = "", file = stderr())
    if (apply) {
      stop(n_bad, " of ", length(results), " file(s) could not be synced. ",
           "See the report above.", call. = FALSE)
    }
  }
  invisible(out)
}
