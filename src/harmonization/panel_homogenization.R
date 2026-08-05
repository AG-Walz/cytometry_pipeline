#!/usr/bin/env Rscript
# =============================================================================
# harmonization/panel_homogenization.R — Panel Homogenization
# =============================================================================
# Combined from: R_scripts/03_normalization/Homogenize.R + Unify.R
# =============================================================================

# Sourced either standalone (from src/<subproject>/) or by a CLI wrapper that
# has already set `shared_dir` and loaded the shared utilities. See the same
# guard in panel_comparison.R.
#
# The flowFrame helpers below need R_shared; the marker-alias layer does not.
# So a shared_dir that cannot be found is not fatal here — it only means the
# flowFrame functions are unusable, which is the honest outcome, rather than
# every caller of the pure functions dying on a source() that guessed wrong.
if (!exists("explore_data", mode = "function")) {
  .panel_shared <- Find(dir.exists, c(
    if (exists("shared_dir")) shared_dir,
    file.path(getwd(), "src", "R_shared"),
    file.path(dirname(getwd()), "src", "R_shared")))
  if (!is.null(.panel_shared)) {
    source(file.path(.panel_shared, "io_utils.R"))
    source(file.path(.panel_shared, "flowframe_utils.R"))
  }
}

# -----------------------------------------------------------------------------
# Marker ALIASES: the same antibody written differently in different panels.
#
# One study's workspaces will happily disagree about spelling — INFy / IFN,
# TNFa / TNF, CD107a / CD107 — and an intersection over those names comes back
# EMPTY, which reads as "these panels have nothing in common" when in fact they
# are the same six markers typed by different people on different days.
#
# Nothing here renames anything on its own. A rename decides which column's
# numbers are plotted under which label, and getting it wrong is a wrong figure,
# not a wrong message — so the alias map is always the caller's, and the only
# thing offered automatically is a SUGGESTION to look at.
#
# The vocabulary matches replace_markernames() below: a named character vector,
# names = the name as written, values = the name to use. That way a map built
# here can be handed straight to the flowFrame path without translation.
# -----------------------------------------------------------------------------

#' Comparison key for two marker names that are "the same, written differently".
#'
#' Case, punctuation and whitespace are dropped, and a single trailing Greek
#' letter written in Latin (a/b/g/y for alpha/beta/gamma) is stripped: TNFa and
#' TNF collapse, IFNg and IFN collapse. That last rule is the only inference
#' made, and it is deliberately narrow — CD4 and CD40 must NOT collapse, so only
#' a trailing letter after a digit-or-letter stem is removed, never a digit.
marker_alias_key <- function(x) {
  key <- tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))
  sub("(?<=[a-z0-9])[abgy]$", "", key, perl = TRUE)
}

#' Grouping key one step looser than marker_alias_key(): the letters, sorted.
#'
#' This catches a TRANSPOSITION — "INFy" for "IFNy" — which is the single most
#' common way a marker name is mistyped and is invisible to any prefix or suffix
#' rule. Sorting the characters cannot merge CD4 with CD40 (different letters),
#' but it will merge genuinely unrelated names that happen to be anagrams, so it
#' is reported at a lower confidence and never applied on its own.
marker_anagram_key <- function(x) {
  vapply(marker_alias_key(x),
         function(k) paste(sort(strsplit(k, "")[[1]]), collapse = ""),
         character(1), USE.NAMES = FALSE)
}

#' Suggest which marker names across panels are probably the same marker.
#'
#' Two rules, reported separately, neither ever applied on its own:
#'
#'   "spelling"      — identical after case, punctuation and a trailing Latin
#'                     Greek letter are dropped. TNFa / TNF. Safe enough that
#'                     the only reason not to take it is a panel that really
#'                     does carry both.
#'   "transposition" — identical after the letters are sorted. INFy / IFN. A
#'                     weaker claim: it is right nearly always and wrong
#'                     spectacularly when it is wrong, so it is a prompt to
#'                     look, not an answer.
#'
#' @param panels Named list of character vectors, one per source.
#' @return A data.frame with one row per suggested group: `key`, `names` (the
#'   spellings seen, comma-separated), `suggested` (the spelling to keep — the
#'   one used by the most sources, ties broken by the longest name, then
#'   alphabetically, so the answer does not depend on argument order),
#'   `sources` (how many sources are covered) and `confidence`. Only groups with
#'   MORE THAN ONE spelling are returned: an agreed name needs no suggestion.
suggest_marker_aliases <- function(panels) {
  empty <- data.frame(key = character(0), names = character(0),
                      suggested = character(0), sources = integer(0),
                      confidence = character(0), stringsAsFactors = FALSE)
  seen <- unlist(lapply(panels, function(x) unique(as.character(x))),
                 use.names = FALSE)
  seen <- seen[!is.na(seen) & nzchar(seen)]
  if (!length(seen)) return(empty)

  counts    <- table(seen)
  spellings <- unique(seen)

  group_by <- function(keys, confidence) {
    rows <- lapply(split(spellings, keys), function(grp) {
      if (length(grp) < 2L) return(NULL)
      n <- as.integer(counts[grp])
      # Most widely used spelling wins; a tie goes to the longer name, which is
      # nearly always the more explicit one (TNFa over TNF), then alphabetically
      # so the result is stable.
      ord <- order(-n, -nchar(grp), grp)
      data.frame(key = keys[match(grp[1], spellings)],
                 names = paste(grp[ord], collapse = ", "),
                 suggested = grp[ord][1],
                 sources = sum(n),
                 confidence = confidence,
                 stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) empty else do.call(rbind, rows)
  }

  strong <- group_by(marker_alias_key(spellings), "spelling")
  loose  <- group_by(marker_anagram_key(spellings), "transposition")
  # A group the strong rule already found needs no weaker restatement: the
  # anagram key merges everything the alias key merges.
  if (nrow(loose) && nrow(strong)) {
    covered <- unlist(strsplit(strong$names, ", ", fixed = TRUE))
    keep <- vapply(strsplit(loose$names, ", ", fixed = TRUE),
                   function(g) !all(g %in% covered), logical(1))
    loose <- loose[keep, , drop = FALSE]
  }
  out <- rbind(strong, loose)
  if (!nrow(out)) return(empty)
  out[order(out$confidence, out$key), , drop = FALSE]
}

#' Apply an alias map to marker names.
#'
#' Unmapped names are returned UNCHANGED — unlike translate_vector(), which
#' turns anything it does not know into "Unknown". That is right for a strict
#' dictionary lookup and wrong here: an alias map covers the handful of names
#' that disagree, and everything else is already correct.
#'
#' @param names Character vector of marker names.
#' @param aliases Named character vector: names = as written, values = to use.
apply_marker_aliases <- function(names, aliases) {
  if (is.null(aliases) || !length(aliases)) return(as.character(names))
  aliases <- unlist(aliases)
  hit <- match(as.character(names), names(aliases))
  out <- as.character(names)
  out[!is.na(hit)] <- unname(aliases[hit[!is.na(hit)]])
  out
}

#' Parse an alias map from the CLI form: "INFy=IFN;TNFa=TNF", or a two-column
#' file. Mirrors parse_population_map() in visualization/umap_from_flowjo.R so
#' the two options are typed the same way.
#'
#' @return A named character vector, or NULL when `spec` is empty.
parse_marker_aliases <- function(spec) {
  if (is.null(spec) || !nzchar(trimws(spec))) return(NULL)
  spec <- trimws(spec)

  if (file.exists(spec)) {
    tbl <- utils::read.table(spec, sep = if (grepl("\\.tsv$", spec)) "\t" else ",",
                             header = FALSE, stringsAsFactors = FALSE,
                             comment.char = "#", quote = "\"'")
    if (ncol(tbl) < 2L) {
      stop("marker alias file needs two columns (as written, to use): ", spec)
    }
    from <- trimws(as.character(tbl[[1]]))
    to   <- trimws(as.character(tbl[[2]]))
    # A header row is common and harmless to drop by name.
    keep <- !(tolower(from) %in% c("from", "marker", "as_written", "alias"))
    from <- from[keep]; to <- to[keep]
  } else {
    pairs <- strsplit(spec, ";", fixed = TRUE)[[1]]
    pairs <- trimws(pairs[nzchar(trimws(pairs))])
    split_one <- function(p) {
      i <- regexpr("=", p, fixed = TRUE)
      if (i < 1L) stop("marker alias '", p, "' is not of the form written=used")
      c(trimws(substr(p, 1L, i - 1L)), trimws(substr(p, i + 1L, nchar(p))))
    }
    parts <- vapply(pairs, split_one, character(2), USE.NAMES = FALSE)
    from <- parts[1, ]; to <- parts[2, ]
  }

  if (!length(from)) return(NULL)
  if (any(!nzchar(from)) || any(!nzchar(to))) {
    stop("a marker alias has an empty side; both names are required.")
  }
  if (anyDuplicated(from)) {
    stop("marker alias names the same spelling twice: ",
         paste(unique(from[duplicated(from)]), collapse = ", "))
  }
  stats::setNames(to, from)
}

# -----------------------------------------------------------------------------
# Recovering marker names from the CHANNEL they were measured on.
#
# A depressingly ordinary state for a real cohort: half the FCS files carry no
# $PnS at all. The detector names ($PnN — "PE-A", "APC-Cy7-A") are always there,
# because the instrument writes them; the marker names are typed by a human into
# the acquisition template and are simply absent when nobody filled the template
# in. Those files then have NO biological columns, contribute nothing to a panel
# intersection, and take the whole run down with them.
#
# The recovery is: within one study, a detector carries one antibody. So the
# files that DO name their markers define channel -> marker, and the files that
# do not can borrow it. That is an inference, and it is stated as one — but it is
# a narrow one, and the alternative is discarding every unlabelled file.
#
# It is only ever applied where $PnS is EMPTY. A file that names its own markers
# keeps them, always, even when it disagrees with the rest of the cohort: that
# disagreement is a fact about the data and belongs in a report, not under a
# silent overwrite.
# -----------------------------------------------------------------------------

#' Build a channel -> marker dictionary from panels that name their markers.
#'
#' @param panels List of list(channels=, markers=) — one per sample.
#' @param aliases Optional alias map applied to marker names before agreeing, so
#'   TNFa and TNF do not count as a conflict.
#' @return list(dictionary = named character (names = channel, values = marker),
#'   conflicts = data.frame(channel, names, chosen, counts) for channels that
#'   were named inconsistently). A conflicted channel still gets an entry — the
#'   most common name wins — because refusing to guess at all would strand every
#'   unlabelled file over one disagreement; the conflict is reported so the
#'   choice is visible and can be overridden with an alias.
channel_marker_dictionary <- function(panels, aliases = NULL) {
  chan <- character(0)
  mark <- character(0)
  for (p in panels) {
    ch <- as.character(p$channels)
    mk <- apply_marker_aliases(as.character(p$markers), aliases)
    keep <- !is.na(mk) & nzchar(mk) & !is.na(ch) & nzchar(ch)
    chan <- c(chan, ch[keep])
    mark <- c(mark, mk[keep])
  }
  empty <- data.frame(channel = character(0), names = character(0),
                      chosen = character(0), counts = character(0),
                      stringsAsFactors = FALSE)
  if (!length(chan)) return(list(dictionary = character(0), conflicts = empty))

  by_channel <- split(mark, chan)
  chosen <- vapply(by_channel, function(m) {
    tab <- sort(table(m), decreasing = TRUE)
    names(tab)[1]
  }, character(1))

  conflicted <- names(by_channel)[vapply(by_channel,
                                         function(m) length(unique(m)) > 1L,
                                         logical(1))]
  conflicts <- if (length(conflicted)) {
    do.call(rbind, lapply(conflicted, function(ch) {
      tab <- sort(table(by_channel[[ch]]), decreasing = TRUE)
      data.frame(channel = ch,
                 names = paste(names(tab), collapse = ", "),
                 chosen = names(tab)[1],
                 counts = paste(sprintf("%s x%d", names(tab), as.integer(tab)),
                                collapse = ", "),
                 stringsAsFactors = FALSE)
    }))
  } else empty

  list(dictionary = chosen, conflicts = conflicts)
}

#' Fill in marker names that the file did not carry, from a channel dictionary.
#'
#' @param channels Character vector of channel names.
#' @param markers Character vector of marker names, "" or NA where absent.
#' @param dictionary Named character vector from channel_marker_dictionary().
#' @return list(markers = the filled vector, recovered = channels that were
#'   filled). Channels the dictionary does not know stay empty.
recover_marker_names <- function(channels, markers, dictionary) {
  markers <- as.character(markers)
  channels <- as.character(channels)
  if (is.null(dictionary) || !length(dictionary)) {
    return(list(markers = markers, recovered = character(0)))
  }
  gap <- is.na(markers) | !nzchar(markers)
  hit <- gap & channels %in% names(dictionary)
  markers[hit] <- unname(dictionary[channels[hit]])
  list(markers = markers, recovered = channels[hit])
}

#' One-line summary of a channel dictionary, for the run log.
format_channel_dictionary <- function(dict) {
  if (!length(dict$dictionary)) return("")
  lines <- c(sprintf("Channel -> marker dictionary built from the files that name their markers (%d channel(s)):",
                     length(dict$dictionary)),
             sprintf("    %-24s %s", names(dict$dictionary),
                     unname(dict$dictionary)))
  if (nrow(dict$conflicts)) {
    lines <- c(lines,
      sprintf("  %d channel(s) are named INCONSISTENTLY across the cohort. The most common",
              nrow(dict$conflicts)),
      "  name was used; override with a marker alias if it is the wrong one:",
      sprintf("    %-24s %s  (using %s)", dict$conflicts$channel,
              dict$conflicts$counts, dict$conflicts$chosen))
  }
  paste(lines, collapse = "\n")
}

#' Translate elements of a vector using a dictionary
#' @param vector Character vector
#' @param dictionary Named list (names = original, values = replacement)
#' @return Translated vector
translate_vector <- function(vector, dictionary) {
  sapply(vector, function(x) {
    if (x %in% names(dictionary)) dictionary[[x]] else "Unknown"
  })
}

#' Replace marker names in a flowFrame
#' @param flowframe A flowFrame object
#' @param replace_map Named list (names = old markers, values = new markers)
#' @return flowFrame with updated marker names
replace_markernames <- function(flowframe, replace_map) {
  marker_list <- flowCore::markernames(flowframe)
  marker_pos <- explore_data(flowframe)$mcols
  for (i in seq_along(marker_list)) {
    if (marker_list[i] %in% names(replace_map)) {
      new_marker <- replace_map[[marker_list[i]]]
      cat("Marker", marker_list[i], "->", new_marker, "\n")
      names(new_marker) <- flowCore::featureNames(flowframe)[marker_pos[i]]
      flowCore::markernames(flowframe)[i] <- new_marker
    }
  }
  flowframe
}

#' Reorder expression matrix to match target marker order
#' @param flowframe A flowFrame object
#' @param marker_name_list Character vector of target marker order
#' @return Reordered expression matrix
reorder_exprmtx <- function(flowframe, marker_name_list) {
  data <- explore_data(flowframe)
  result <- matrix(NA, nrow = nrow(flowCore::exprs(flowframe)),
                   ncol = length(marker_name_list))
  colnames(result) <- marker_name_list
  missing_markers <- c()
  for (j in seq_along(marker_name_list)) {
    col_idx <- which(data$markers == marker_name_list[j])
    if (length(col_idx) > 0) {
      result[, j] <- flowCore::exprs(flowframe)[, col_idx[1]]
    } else {
      missing_markers <- c(missing_markers, marker_name_list[j])
    }
  }
  if (length(missing_markers) > 0) {
    warning("Markers filled with NA (not found in flowFrame): ",
            paste(missing_markers, collapse = ", "))
  }
  result
}

#' Homogenize a list of flowFrames to a common panel
#'
#' @param flow_frames Named list of flowFrame objects
#' @param sample_names Character vector of sample names
#' @param markers_of_interest Target marker list
#' @param markername_replacing Named list for marker renaming
#' @param save_csv Save reordered CSVs
#' @param save_fcs Save reordered FCS files
#' @param output_dir Output directory
#' @return Named list of homogenized flowFrame objects
homogenize_flowframes <- function(flow_frames, sample_names,
                                   markers_of_interest,
                                   markername_replacing = list(),
                                   save_csv = FALSE,
                                   save_fcs = FALSE,
                                   output_dir = NULL) {
  cat("Homogenizing", length(flow_frames), "samples to",
      length(markers_of_interest), "markers\n")

  if (!is.null(output_dir)) make_outdir(output_dir)

  hg_frames <- list()
  cell_ids <- generate_cellID_list(flow_frames)

  for (i in seq_along(flow_frames)) {
    sname <- sample_names[i]
    cat("  Homogenizing", sname, "...")

    ff <- flow_frames[[i]]

    # Rename markers if needed
    if (length(markername_replacing) > 0) {
      # Check for collisions: different old names mapping to the same new name
      new_names <- unlist(markername_replacing)
      if (any(duplicated(new_names))) {
        dupes <- unique(new_names[duplicated(new_names)])
        stop("markername_replacing collision: multiple markers map to the same name: ",
             paste(dupes, collapse = ", "))
      }
      ff <- replace_markernames(ff, markername_replacing)
    }

    # Reorder to target panel
    mtx <- reorder_exprmtx(ff, markers_of_interest)

    # Convert back to flowFrame
    hg_ff <- from_mtx_to_ffr(mtx, markers_of_interest)
    hg_frames[[sname]] <- hg_ff

    cat(" ", nrow(mtx), "cells x", ncol(mtx), "markers\n")

    # Save if requested
    if (save_csv && !is.null(output_dir)) {
      write.csv(mtx, file.path(output_dir, paste0("hg_", sname, ".csv")))
    }
    if (save_fcs && !is.null(output_dir)) {
      flowCore::write.FCS(hg_ff, file.path(output_dir, paste0("hg_", sname, ".fcs")))
    }
  }

  hg_frames
}

# Helper: generate cell IDs for a list of flowFrames
generate_cellID_list <- function(flow_frames) {
  total <- 0
  ids <- list()
  for (i in seq_along(flow_frames)) {
    n <- nrow(flow_frames[[i]])
    ids[[i]] <- (total + 1):(total + n)
    total <- total + n
  }
  ids
}
