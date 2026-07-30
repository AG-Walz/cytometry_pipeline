#!/usr/bin/env Rscript
# =============================================================================
# loading/sample_annotations.R — Per-sample annotations (patient, timepoint, …)
# =============================================================================
# Attaches arbitrary categorical annotations to samples so an analysis can group
# by them. Nothing here is specific to "patient" or "timepoint": a column is a
# column, so response, cohort arm, batch and treatment all work through the same
# path and no new code is needed to add one.
#
# Two sources, in priority order:
#
#   1. A metadata CSV keyed on sample. This is authoritative and follows the
#      repo's existing `sample_id` convention (see metadata_qc/check_metadata.R),
#      so a clinical table you already maintain can be used as-is.
#   2. A regular expression over the sample name, whose NAMED capture groups
#      become the columns. Nothing to write when the names already encode it —
#      "MTB048_TP1.fcs" with "(?<patient>MTB[0-9]+)_(?<timepoint>TP[0-9]+)".
#
# Matching a CSV row to a sample is deliberately tolerant, because a FlowJo
# sample name ("MTB048_TP1.fcs_0"), the FCS file name ("MTB048_TP1.fcs") and
# whatever someone typed into a spreadsheet ("MTB048_TP1") are rarely identical.
#
# A sample that cannot be annotated gets "Unassigned" rather than being dropped
# — the same placeholder check_metadata.R uses — so a partial metadata table
# still produces a usable plot and the gap is visible in the legend instead of
# silently shrinking the data.
#
# This file holds only function definitions (no top-level side effects) and
# needs no cytometry packages, so it is testable with base R alone.
# =============================================================================

#' Placeholder for a sample with no annotation. Matches check_metadata.R.
.ANNOT_UNASSIGNED <- "Unassigned"

# -----------------------------------------------------------------------------
# Normalised join key: lower-cased base name without the .fcs extension or the
# cytoset "_0" suffix a GatingSet appends. Windows separators are handled
# because base R's basename() only splits on "/".
# -----------------------------------------------------------------------------
annotation_key <- function(x) {
  x <- sub(".*[\\\\/]", "", as.character(x))   # base name, either separator
  x <- sub("\\.fcs_[0-9]+$", "", x, ignore.case = TRUE)
  x <- sub("\\.fcs$", "", x, ignore.case = TRUE)
  tolower(trimws(x))
}

# -----------------------------------------------------------------------------
# Read a metadata CSV. The separator is auto-detected (semicolon-separated
# exports are common in this project), and every column is read as text: an
# annotation is a label, and a timepoint like "01" must not become the number 1.
# -----------------------------------------------------------------------------
read_sample_metadata <- function(path, sep = NULL) {
  if (!file.exists(path)) stop("sample_metadata file not found: ", path)
  if (is.null(sep)) {
    first <- readLines(path, n = 1, warn = FALSE)
    sep <- if (length(first) && grepl(";", first)) ";" else ","
  }
  md <- utils::read.table(path, sep = sep, header = TRUE, quote = "\"",
                          stringsAsFactors = FALSE, colClasses = "character",
                          check.names = FALSE, comment.char = "")
  if (nrow(md) == 0) stop("sample_metadata file has no rows: ", path)
  md
}

# -----------------------------------------------------------------------------
# Which column identifies the sample. `sample_id` by convention; otherwise the
# first column, which is what a hand-made table almost always uses.
# -----------------------------------------------------------------------------
metadata_id_column <- function(md, id_column = NULL) {
  if (!is.null(id_column) && nzchar(id_column)) {
    if (!id_column %in% colnames(md)) {
      stop("sample_id_column '", id_column, "' is not in the metadata (saw: ",
           paste(colnames(md), collapse = ", "), ").")
    }
    return(id_column)
  }
  for (cand in c("sample_id", "sample", "sampleID", "sample_name", "filename",
                 "file")) {
    if (cand %in% colnames(md)) return(cand)
  }
  colnames(md)[1]
}

# -----------------------------------------------------------------------------
# Annotations parsed out of the sample names themselves.
#
# `pattern` must use NAMED capture groups — (?<patient>...) — because the group
# names become the column names. A name the pattern does not match yields NA for
# every column, which the caller turns into "Unassigned".
# -----------------------------------------------------------------------------
annotations_from_pattern <- function(sample_names, pattern) {
  m      <- regexpr(pattern, sample_names, perl = TRUE)
  starts <- attr(m, "capture.start")
  lens   <- attr(m, "capture.length")
  nms    <- attr(m, "capture.names")

  if (is.null(starts) || is.null(nms) || !any(nzchar(nms))) {
    stop("sample_pattern has no named capture groups, so there is nothing to ",
         "name the columns after. Use (?<name>...), e.g.\n  ",
         "\"(?<patient>MTB[0-9]+)_(?<timepoint>TP[0-9]+)\"")
  }

  keep <- which(nzchar(nms))
  out <- lapply(keep, function(j) {
    v <- substring(sample_names, starts[, j], starts[, j] + lens[, j] - 1L)
    v[m == -1L | lens[, j] == 0L] <- NA_character_
    v
  })
  names(out) <- nms[keep]
  as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
}

# -----------------------------------------------------------------------------
# Resolve annotations for every sample in the run.
#
# `samples` is a data.frame with at least `display` (the run-wide label) and
# `raw` (the GatingSet sample name); `workspace` is used as a fallback key so a
# metadata table keyed on workspace still works.
#
# Returns a data.frame with one row per sample, in the SAME ORDER as `samples`,
# holding only the requested annotation columns. Every value is a character; a
# missing one is .ANNOT_UNASSIGNED.
# -----------------------------------------------------------------------------
resolve_sample_annotations <- function(samples, metadata_file = NULL,
                                       pattern = NULL, annotate = NULL,
                                       id_column = NULL) {
  n <- nrow(samples)
  empty <- data.frame(row.names = seq_len(n))

  if (is.null(metadata_file) && is.null(pattern)) return(empty)

  if (!is.null(metadata_file)) {
    md   <- read_sample_metadata(metadata_file)
    idc  <- metadata_id_column(md, id_column)
    keys <- annotation_key(md[[idc]])
    if (anyDuplicated(keys)) {
      dup <- unique(keys[duplicated(keys)])
      warning("sample_metadata has duplicate sample ids (first row wins): ",
              paste(utils::head(dup, 5), collapse = ", "))
    }

    # Try each identity a sample is known by, most specific first, so a table
    # keyed on the FCS file name and one keyed on the display label both work.
    cand <- list(samples$raw, samples$display, samples$clean,
                 samples$workspace)
    cand <- cand[!vapply(cand, is.null, logical(1))]
    row_of <- rep(NA_integer_, n)
    for (k in cand) {
      todo <- is.na(row_of)
      if (!any(todo)) break
      row_of[todo] <- match(annotation_key(k[todo]), keys)
    }

    if (any(is.na(row_of))) {
      miss <- samples$display[is.na(row_of)]
      warning(sum(is.na(row_of)), " of ", n, " sample(s) are not in ",
              "sample_metadata and will be '", .ANNOT_UNASSIGNED, "': ",
              paste(utils::head(miss, 5), collapse = ", "),
              if (length(miss) > 5) ", ..." else "")
    }

    value_cols <- setdiff(colnames(md), idc)
    ann <- as.data.frame(
      lapply(value_cols, function(cl) md[[cl]][row_of]),
      stringsAsFactors = FALSE, optional = TRUE)
    names(ann) <- value_cols

  } else {
    ann <- annotations_from_pattern(samples$display, pattern)
    unmatched <- rowSums(!is.na(ann)) == 0L
    if (any(unmatched)) {
      warning(sum(unmatched), " of ", n, " sample name(s) do not match ",
              "sample_pattern and will be '", .ANNOT_UNASSIGNED, "': ",
              paste(utils::head(samples$display[unmatched], 5),
                    collapse = ", "))
    }
  }

  # Restrict to the requested columns, in the order asked for.
  if (!is.null(annotate) && length(annotate)) {
    missing_cols <- setdiff(annotate, colnames(ann))
    if (length(missing_cols)) {
      stop("annotate column(s) not available: ",
           paste(missing_cols, collapse = ", "),
           "\nAvailable: ",
           if (ncol(ann)) paste(colnames(ann), collapse = ", ") else "(none)")
    }
    ann <- ann[, annotate, drop = FALSE]
  }

  if (ncol(ann) == 0) return(empty)

  # Blanks and the usual placeholders read as missing, matching the convention
  # in metadata_qc/check_metadata.R.
  placeholders <- c("?", "na", "n/a", "null", "-", "")
  for (cl in colnames(ann)) {
    v <- trimws(as.character(ann[[cl]]))
    v[is.na(v) | tolower(v) %in% placeholders] <- .ANNOT_UNASSIGNED
    ann[[cl]] <- v
  }
  rownames(ann) <- NULL
  ann
}
