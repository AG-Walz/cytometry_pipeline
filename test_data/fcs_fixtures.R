# =============================================================================
# test_data/fcs_fixtures.R — build FCS files byte by byte, including broken ones
# =============================================================================
# A keyword rewriter has to be trusted with irreplaceable data, so its tests
# cannot lean on flowCore to produce their inputs: flowCore writes one dialect
# (form-feed-free, zero slack, no $FIL, no ANALYSIS) and that is precisely the
# dialect least likely to expose a bug. These builders emit the header, the
# TEXT segment, the DATA bytes and the CRC field directly, with knobs for every
# layout the FCS standard permits and several it does not.
#
# Needs nothing but base R. Sourced by test_fcs_keywords.R and
# test_sync_fcs_names.R.
# =============================================================================

# Deterministic event bytes, so "the data segment is unchanged" is a real
# assertion rather than a comparison of two blocks of zeros.
fcs_fixture_data <- function(n, seed = 1L) {
  as.raw(((seq_len(n) * 7L + seed * 13L) - 1L) %% 256L)
}

# Double any delimiter byte appearing inside a keyword or value (3.2.7).
.fx_escape <- function(s, delim) {
  gsub(delim, paste0(delim, delim), s, fixed = TRUE)
}

.fx_text <- function(pairs, delim) {
  parts <- character(0)
  for (nm in names(pairs)) {
    parts <- c(parts, .fx_escape(nm, delim), .fx_escape(as.character(pairs[[nm]]), delim))
  }
  paste0(delim, paste(parts, collapse = delim), delim)
}

.fx_field8 <- function(v) formatC(sprintf("%.0f", v), width = 8, flag = " ")

#' Write a synthetic FCS file.
#'
#' @param path Destination.
#' @param keywords Named list of extra keywords, merged over the required set
#'   (same name replaces, new name is appended).
#' @param delimiter One-character TEXT delimiter.
#' @param slack Space bytes inserted between the end of TEXT and the next
#'   segment — the room an in-place patch can grow into.
#' @param analysis Raw vector for an ANALYSIS segment, or NULL.
#' @param header_zero,header_sentinel Segment names ("data", "analysis") whose
#'   HEADER offsets are written as 0 / 99999999 instead of the real value, so
#'   the truth lives only in the TEXT keywords.
#' @param offset_keywords FALSE drops $BEGINDATA/$ENDDATA entirely.
#' @param raw_text Verbatim TEXT segment (already delimited). Bypasses the
#'   keyword machinery — for deliberately malformed segments.
#' @return `path`, invisibly.
fcs_fixture <- function(path,
                        keywords = list(),
                        delimiter = "|",
                        version = "FCS3.0",
                        par = 2L, tot = 5L,
                        datatype = "F", pnb = 32L,
                        data_raw = NULL,
                        slack = 0L,
                        analysis = NULL,
                        analysis_before_data = FALSE,
                        trailer = "00000000",
                        header_zero = character(),
                        header_sentinel = character(),
                        text_after_data = FALSE,
                        force_textend = NULL,
                        offset_keywords = TRUE,
                        other_bytes = 0L,
                        raw_text = NULL) {

  if (is.null(data_raw)) data_raw <- fcs_fixture_data(tot * par * pnb / 8)
  n_data <- length(data_raw)
  n_ana  <- if (is.null(analysis)) 0L else length(analysis)

  base <- list()
  if (offset_keywords) {
    base[["$BEGINANALYSIS"]] <- "0"
    base[["$BEGINDATA"]] <- "0"
  }
  base[["$BEGINSTEXT"]] <- "0"
  base[["$BYTEORD"]] <- "4,3,2,1"
  base[["$DATATYPE"]] <- datatype
  if (offset_keywords) {
    base[["$ENDANALYSIS"]] <- "0"
    base[["$ENDDATA"]] <- "0"
  }
  base[["$ENDSTEXT"]] <- "0"
  base[["$MODE"]] <- "L"
  base[["$NEXTDATA"]] <- "0"
  for (i in seq_len(par)) {
    base[[paste0("$P", i, "B")]] <- as.character(pnb)
    base[[paste0("$P", i, "N")]] <- paste0("FL", i)
    base[[paste0("$P", i, "R")]] <- "1024"
  }
  base[["$PAR"]] <- as.character(par)
  base[["$TOT"]] <- as.character(tot)
  for (nm in names(keywords)) base[[nm]] <- keywords[[nm]]

  textstart <- 58L + other_bytes

  # $BEGINDATA's digits are part of the length that decides $BEGINDATA, so the
  # builder has the same fixed point the writer does. Three rounds is plenty;
  # ten is paranoia.
  layout <- NULL
  for (round in 1:10) {
    text <- if (!is.null(raw_text)) raw_text else .fx_text(base, delimiter)
    tlen <- nchar(text, type = "bytes")

    if (text_after_data) {
      datastart <- textstart
      dataend   <- datastart + n_data - 1L
      anastart  <- if (n_ana) dataend + 1L else 0L
      anaend    <- if (n_ana) anastart + n_ana - 1L else 0L
      ts <- (if (n_ana) anaend else dataend) + 1L + slack
      te <- ts + tlen - 1L
      new_layout <- list(textstart = ts, textend = te, datastart = datastart,
                         dataend = dataend, anastart = anastart, anaend = anaend,
                         text = text)
    } else {
      te <- textstart + tlen - 1L
      first <- te + 1L + slack
      if (n_ana && analysis_before_data) {
        anastart <- first; anaend <- anastart + n_ana - 1L
        datastart <- anaend + 1L; dataend <- datastart + n_data - 1L
      } else {
        datastart <- first; dataend <- datastart + n_data - 1L
        anastart <- if (n_ana) dataend + 1L else 0L
        anaend   <- if (n_ana) anastart + n_ana - 1L else 0L
      }
      new_layout <- list(textstart = textstart, textend = te,
                         datastart = datastart, dataend = dataend,
                         anastart = anastart, anaend = anaend, text = text)
    }

    if (!is.null(raw_text)) { layout <- new_layout; break }
    if (offset_keywords) {
      base[["$BEGINDATA"]] <- as.character(new_layout$datastart)
      base[["$ENDDATA"]]   <- as.character(new_layout$dataend)
      base[["$BEGINANALYSIS"]] <- as.character(new_layout$anastart)
      base[["$ENDANALYSIS"]]   <- as.character(new_layout$anaend)
    }
    if (!is.null(layout) && identical(layout[names(layout) != "text"],
                                      new_layout[names(new_layout) != "text"])) {
      layout <- new_layout
      break
    }
    layout <- new_layout
  }

  textend <- if (is.null(force_textend)) layout$textend else force_textend

  hdr <- charToRaw(paste0(
    formatC(version, width = -6),
    "    ",
    .fx_field8(layout$textstart),
    .fx_field8(textend),
    .fx_field8(if ("data" %in% header_zero) 0
               else if ("data" %in% header_sentinel) 99999999
               else layout$datastart),
    .fx_field8(if ("data" %in% header_zero) 0
               else if ("data" %in% header_sentinel) 99999999
               else layout$dataend),
    .fx_field8(if ("analysis" %in% header_zero) 0
               else if ("analysis" %in% header_sentinel) 99999999
               else layout$anastart),
    .fx_field8(if ("analysis" %in% header_zero) 0
               else if ("analysis" %in% header_sentinel) 99999999
               else layout$anaend)))
  stopifnot(length(hdr) == 58L)

  last_end <- max(layout$textend, layout$dataend, layout$anaend)
  total <- last_end + 1L + (if (is.null(trailer)) 0L else 8L)
  buf <- rep(as.raw(0x20), total)
  buf[1:58] <- hdr
  put <- function(start, bytes) {
    if (length(bytes)) buf[(start + 1L):(start + length(bytes))] <<- bytes
  }
  put(layout$textstart, charToRaw(layout$text))
  put(layout$datastart, data_raw)
  if (n_ana) put(layout$anastart, analysis)
  if (!is.null(trailer)) put(last_end + 1L, charToRaw(trailer))

  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(buf, con)
  invisible(path)
}

#' md5 of a file — the "we did not touch it" witness used all over the suites.
fx_md5 <- function(path) unname(tools::md5sum(path))

#' Overwrite bytes at an offset, in place.
#'
#' Uses rw = "write" deliberately: an "r+b" connection keeps separate read and
#' write positions and seek() moves the READ one by default, so the obvious
#' version of this helper silently writes at offset 0.
fx_poke <- function(path, offset, bytes) {
  if (is.character(bytes)) bytes <- charToRaw(bytes)
  con <- file(path, "r+b")
  on.exit(close(con), add = TRUE)
  seek(con, offset, origin = "start", rw = "write")
  writeBin(bytes, con)
  invisible(path)
}

#' Read an arbitrary byte range of a file, for spot assertions on layout.
fx_bytes <- function(path, from, n) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, from)
  readBin(con, "raw", n = as.integer(n))
}
