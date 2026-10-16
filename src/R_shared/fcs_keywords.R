# =============================================================================
# R_shared/fcs_keywords.R — read and REWRITE FCS keywords, without touching events
# =============================================================================
# Everything here is base R. No flowCore, no Bioconductor, no optparse — which
# is deliberate: this is the one layer that has to run on a plain R install
# (an analyst's laptop, a Windows box, a CI container) because its whole job is
# to repair files before the cytometry stack will accept them. It also means
# the test suite runs anywhere. `fcs_utils.R` is the flowCore-backed reader and
# stays that way; this file does not depend on it and it does not depend on this.
#
# WHY IT EXISTS
# An FCS file records the name it was acquired under in its TEXT segment
# ($FIL). Rename the file on disk and that keyword still says the old thing, so
# a folder of files acquired under one plate template all claim the same name
# and FlowJo shows n identical rows. FlowJo cannot fix it — it never writes back
# to raw data — and re-exporting re-encodes every event.
#
# THE IDEA
# The DATA segment is an opaque byte range addressed by ASCII offsets in the
# 58-byte HEADER (and duplicated in TEXT keywords). So it never has to be
# decoded, only located and, at worst, moved. Three write strategies, chosen per
# file by how much room sits between the start of TEXT and whatever follows it
# (`room = body_start - textstart`):
#
#   neutral   new TEXT is no longer than the old TEXT. Overwrite [textstart,
#             body_start) with the new TEXT plus space fill. DATA never moves;
#             $BEGINDATA/$ENDDATA keep their exact original bytes; file size is
#             unchanged. Only `textend` in the header moves, and only if the
#             length actually changed.
#   patch     new TEXT is longer, but still fits in `room` (the slack some
#             writers leave before DATA). Identical write, identical guarantees.
#   rewrite   new TEXT does not fit. The whole file is re-emitted: patched
#             header, new TEXT, then [body_start, EOF) copied VERBATIM in
#             chunks. Bytes are copied, never reconstructed.
#
# Note the consequence that makes this tractable: the notorious circular
# dependency — $BEGINDATA's value determines TEXT's length which determines
# $BEGINDATA — only exists on the `rewrite` path. On neutral/patch nothing
# moves, so the offset keywords are simply left alone.
#
# CONVENTIONS THAT ARE NOT NEGOTIABLE, AND WHY
#   * Offsets are DOUBLES from parse to write. R integers stop at 2^31 and a
#     spectral run passes that routinely, so `as.integer`/`strtoi` on an offset
#     is a silent NA. Formatting is sprintf("%.0f") / sprintf("%8.0f") only —
#     paste()/format() would hand you "2e+09".
#   * Byte offsets are 0-BASED and segment ends are INCLUSIVE, exactly as the
#     spec defines them (FCS 3.1 2.2.4: "The last byte of a segment is
#     understood as the end of a segment"). seek() is 0-based too, so the two
#     line up and there is nothing to convert.
#   * Keys and values are stored as RAW, never character. The wild is full of
#     latin-1 bytes (µ, °, accented $PnS/$CYT) and a character round-trip
#     mangles them; length arithmetic must be in bytes regardless.
#
# Function definitions only — no top-level side effects.
# =============================================================================

# -----------------------------------------------------------------------------
# Small byte helpers
# -----------------------------------------------------------------------------

# The 8-byte header offset fields are right-justified ASCII integers, and the
# spec allows an absent ANALYSIS segment to be written as zeros OR left blank:
# "these last two byte offsets can be set to zero (right justified) or left
# blank (filled with space characters)". So blank means 0. Anything that is not
# a run of digits is NA, which the validator turns into a refusal rather than a
# guess.
#
# Byte-level on purpose: this is also pointed at the OTHER-offset region, whose
# bytes may be arbitrary binary. rawToChar() on those produces a string that is
# not valid in the session encoding, and trimws()/grepl() on it throw rather
# than return FALSE.
.fcs_num_field <- function(bytes) {
  b <- as.integer(bytes)
  b <- b[b != 32L]
  if (!length(b)) return(0)
  if (any(b < 48L | b > 57L)) return(NA_real_)
  as.numeric(rawToChar(as.raw(b)))
}

# Render an offset into an 8-byte header field. Returns NULL when it does not
# fit, which is the caller's signal to fall back to the "0 in the header, truth
# in TEXT" convention rather than truncating.
.fcs_field8 <- function(value) {
  if (!is.finite(value) || value < 0 || value > 99999999) return(NULL)
  s <- sprintf("%8.0f", value)
  if (nchar(s, type = "bytes") != 8L) return(NULL)
  charToRaw(s)
}

# Keyword names are ASCII 32-126 per the spec and case-insensitive per 3.2.10.
# A key with bytes outside that range is not something we should be indexing,
# so it comes back NA and becomes a refusal.
.fcs_key_string <- function(bytes) {
  if (length(bytes) == 0L) return(NA_character_)
  n <- as.integer(bytes)
  if (any(n < 32L | n > 126L)) return(NA_character_)
  toupper(trimws(rawToChar(bytes)))
}

# Values are for display and for numeric keywords only. latin-1 with sub="byte"
# is what the rest of this repo does (fcs_utils.R:19-22) and it never errors.
.fcs_value_string <- function(bytes) {
  if (length(bytes) == 0L) return("")
  iconv(rawToChar(bytes), "", "latin1", sub = "byte")
}

.fcs_nbytes <- function(x) nchar(x, type = "bytes")

# Position a read/write connection for WRITING, and prove it landed.
#
# This is not defensive padding. A connection opened "r+b" keeps SEPARATE read
# and write positions, and seek()'s `rw` argument defaults to the read one — so
# the obvious `seek(con, offset); writeBin(...)` silently writes wherever the
# write pointer happened to be. It is invisible in the common case, because
# after writing the 58-byte header the write pointer sits at 58, which is
# exactly where TEXT starts in most files. It corrupts any file with an OTHER
# segment. Hence: always rw = "write", always assert.
.fcs_seek_write <- function(con, where) {
  seek(con, where, origin = "start", rw = "write")
  at <- seek(con, NA, rw = "write")
  if (!isTRUE(at == where)) {
    stop(sprintf("could not position the file at %.0f (landed at %.0f)",
                 where, at))
  }
  invisible(where)
}

# -----------------------------------------------------------------------------
# HEADER
# -----------------------------------------------------------------------------

#' Read the 58-byte FCS HEADER.
#'
#' @param path Path to an FCS file.
#' @return list with version, fcs_major, the six offsets (doubles, 0-based,
#'   ends inclusive), the raw OTHER-offset region, file_size and the raw
#'   header bytes. `bad_fields` names any offset field that was not numeric.
fcs_read_header <- function(path) {
  size <- file.info(path)$size
  if (is.na(size)) stop("FCS file not found: ", path)
  size <- as.numeric(size)
  if (size < 58) stop("file is shorter than an FCS header (58 bytes): ", path)

  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  hdr <- readBin(con, "raw", n = 58L)

  version <- rawToChar(hdr[1:6])
  fcs_major <- suppressWarnings(as.numeric(substr(version, 4, 4)))
  if (is.na(fcs_major)) fcs_major <- NA_real_

  fields <- list(textstart = 11:18, textend = 19:26,
                 datastart = 27:34, dataend  = 35:42,
                 anastart  = 43:50, anaend   = 51:58)
  vals <- vapply(fields, function(ix) .fcs_num_field(hdr[ix]), numeric(1))
  bad  <- names(vals)[is.na(vals)]

  # A TEXT segment that does not start at 58 means user-defined OTHER segments
  # were written between the header and TEXT. We keep their bytes so we can
  # copy them verbatim, and refuse to shift them (see fcs_validate).
  # Capped: in a TEXT-after-DATA layout this "gap" is the entire data segment,
  # and reading a 4 GB block to check whether it is padding would be absurd.
  # Oversize means we cannot vouch for it, which the validator treats the same
  # as a real OTHER segment.
  other_raw <- raw(0)
  other_oversize <- FALSE
  ts <- vals[["textstart"]]
  if (!is.na(ts) && ts > 58 && ts <= size) {
    if (ts - 58 > 8192) {
      other_oversize <- TRUE
    } else {
      seek(con, 58)
      other_raw <- readBin(con, "raw", n = as.integer(ts - 58))
    }
  }

  list(path = path, version = version, fcs_major = fcs_major,
       textstart = vals[["textstart"]], textend = vals[["textend"]],
       datastart = vals[["datastart"]], dataend = vals[["dataend"]],
       anastart  = vals[["anastart"]],  anaend  = vals[["anaend"]],
       bad_fields = bad, other_raw = other_raw,
       other_oversize = other_oversize,
       file_size = size, raw = hdr)
}

# -----------------------------------------------------------------------------
# TEXT: parse and serialize
# -----------------------------------------------------------------------------

#' Parse an FCS TEXT segment from its raw bytes.
#'
#' The delimiter is whatever byte the segment starts with, and a literal
#' delimiter inside a keyword or value is written twice (3.2.7: "If the
#' delimiter appears in a keyword or keyword value, it must be immediately
#' followed by a second delimiter"). So a run of k consecutive delimiters
#' contributes floor(k/2) literal bytes to the current token and ends that
#' token only when k is odd.
#'
#' This is why strsplit(fixed = TRUE) — what fcs_utils.R:27 uses on the FCS2.0
#' path — cannot be reused here: it splits on every delimiter and silently
#' shears an escaped value into two tokens, which for a WRITER would scramble
#' the key/value alternation of the whole segment.
#'
#' @param txt_raw Raw vector: the whole TEXT segment, delimiter to delimiter.
#' @return An "fcs_text" object. `malformed` is a code when the segment could
#'   not be read as key/value pairs, otherwise NA.
fcs_parse_text <- function(txt_raw) {
  mk <- function(delim, keys, values, malformed = NA_character_) {
    keys_upper <- vapply(keys, .fcs_key_string, character(1), USE.NAMES = FALSE)
    structure(list(delimiter = delim, keys_raw = keys, values_raw = values,
                   keys_upper = keys_upper, n_pairs = length(keys),
                   malformed = malformed), class = "fcs_text")
  }

  if (length(txt_raw) < 2L) {
    return(mk(if (length(txt_raw)) txt_raw[1] else as.raw(0),
              list(), list(), "malformed_text_unpaired"))
  }
  delim <- txt_raw[1]

  body <- txt_raw[-1]
  runs <- rle(body == delim)

  tokens <- vector("list", 64L)
  ntok <- 0L
  pieces <- vector("list", 8L)
  npc <- 0L
  idx <- 1L

  for (i in seq_along(runs$lengths)) {
    len <- runs$lengths[i]
    if (!runs$values[i]) {
      npc <- npc + 1L
      if (npc > length(pieces)) length(pieces) <- 2L * npc
      pieces[[npc]] <- body[idx:(idx + len - 1L)]
    } else {
      half <- len %/% 2L
      if (half > 0L) {
        npc <- npc + 1L
        if (npc > length(pieces)) length(pieces) <- 2L * npc
        pieces[[npc]] <- rep(delim, half)
      }
      if (len %% 2L == 1L) {
        ntok <- ntok + 1L
        if (ntok > length(tokens)) length(tokens) <- 2L * ntok
        tokens[[ntok]] <- if (npc == 0L) raw(0) else
          unlist(pieces[seq_len(npc)], use.names = FALSE)
        pieces <- vector("list", 8L)
        npc <- 0L
      }
    }
    idx <- idx + len
  }

  # Anything still accumulating means the segment did not end on a delimiter.
  if (npc > 0L) return(mk(delim, list(), list(), "malformed_text_unterminated"))
  tokens <- tokens[seq_len(ntok)]
  if (ntok == 0L || ntok %% 2L != 0L) {
    return(mk(delim, list(), list(), "malformed_text_unpaired"))
  }

  keys   <- tokens[seq(1L, ntok, by = 2L)]
  values <- tokens[seq(2L, ntok, by = 2L)]
  # 3.2.9: "Keywords and keyword values must have lengths greater than zero."
  if (any(vapply(keys, length, integer(1)) == 0L) ||
      any(vapply(values, length, integer(1)) == 0L)) {
    return(mk(delim, keys, values, "empty_keyword_or_value"))
  }
  mk(delim, keys, values)
}

#' Serialize an "fcs_text" back to raw bytes. Exact inverse of fcs_parse_text().
fcs_serialize_text <- function(txt) {
  delim <- txt$delimiter
  esc <- function(tok) {
    if (!any(tok == delim)) return(tok)
    out <- vector("list", length(tok))
    for (i in seq_along(tok)) {
      out[[i]] <- if (tok[i] == delim) c(delim, delim) else tok[i]
    }
    unlist(out, use.names = FALSE)
  }
  n <- txt$n_pairs
  parts <- vector("list", 2L * n + 1L)
  parts[[1]] <- delim
  for (i in seq_len(n)) {
    parts[[2L * i]]      <- c(esc(txt$keys_raw[[i]]), delim)
    parts[[2L * i + 1L]] <- c(esc(txt$values_raw[[i]]), delim)
  }
  unlist(parts, use.names = FALSE)
}

#' Index of a keyword, case-insensitively (3.2.10: "an FCS file reader must
#' ignore keyword case"). Returns integer(0) when absent.
fcs_text_index <- function(txt, key) {
  which(txt$keys_upper == toupper(trimws(key)))
}

#' Value of a keyword as a display string, or `default` when absent.
fcs_text_get <- function(txt, key, default = NA_character_) {
  i <- fcs_text_index(txt, key)
  if (length(i) == 0L) return(default)
  .fcs_value_string(txt$values_raw[[i[1]]])
}

#' Numeric value of a keyword. Leading zeros are legal (3.2.17); space padding
#' is not, but we tolerate it on read because writers do it anyway.
fcs_text_num <- function(txt, key, default = NA_real_) {
  v <- fcs_text_get(txt, key, NA_character_)
  if (is.na(v)) return(default)
  v <- trimws(v)
  if (!grepl("^[0-9]+$", v)) return(NA_real_)
  as.numeric(v)
}

#' Set / remove keywords, returning a new "fcs_text".
#'
#' A keyword that exists is replaced in place (its on-disk key bytes are kept,
#' so an unchanged file still re-serializes bit-for-bit); one that does not is
#' appended. Matching is case-insensitive.
fcs_text_set <- function(txt, set = character(), unset = character()) {
  keys   <- txt$keys_raw
  values <- txt$values_raw
  upper  <- txt$keys_upper

  for (k in unset) {
    i <- which(upper == toupper(trimws(k)))
    if (length(i)) {
      keys   <- keys[-i]
      values <- values[-i]
      upper  <- upper[-i]
    }
  }
  for (nm in names(set)) {
    val_raw <- charToRaw(enc2utf8(as.character(set[[nm]])))
    i <- which(upper == toupper(trimws(nm)))
    if (length(i)) {
      values[[i[1]]] <- val_raw
    } else {
      keys[[length(keys) + 1L]]     <- charToRaw(enc2utf8(nm))
      values[[length(values) + 1L]] <- val_raw
      upper <- c(upper, toupper(trimws(nm)))
    }
  }
  structure(list(delimiter = txt$delimiter, keys_raw = keys,
                 values_raw = values, keys_upper = upper,
                 n_pairs = length(keys), malformed = txt$malformed),
            class = "fcs_text")
}

# -----------------------------------------------------------------------------
# Segment layout
# -----------------------------------------------------------------------------

# Resolve one segment's begin/end from the header field and the TEXT keyword.
#
# The spec's own rule (3.1, and repeated per keyword): when a segment fits in
# the first 99,999,999 bytes the offsets appear in BOTH places and agree; when
# it does not, the header carries 0 and "the relevant offset information is
# stored in the keyword value only". 99999999 is a fourth convention that is not
# in the spec but that flowCore's checkOffset() hard-codes, so we honour it too.
#
# style is recorded because the rewrite path must reproduce the ORIGINAL
# writer's convention rather than re-deriving one — flowCore-written files are
# already rejected by BD FACSDiva 6 (flowCore #147), so "flowCore reads it back"
# is not evidence that an instrument vendor will.
.fcs_resolve_segment <- function(hdr_begin, hdr_end, kw_begin, kw_end) {
  sentinel <- function(v) !is.na(v) && (v == 0 || v == 99999999)
  style <- if (is.na(hdr_begin)) "value"
           else if (hdr_begin == 0) "zero"
           else if (hdr_begin == 99999999) "sentinel"
           else "value"

  conflict <- FALSE
  if (!is.na(hdr_begin) && !sentinel(hdr_begin) && !is.na(kw_begin) &&
      hdr_begin != kw_begin) conflict <- TRUE
  if (!is.na(hdr_end) && !sentinel(hdr_end) && !is.na(kw_end) &&
      hdr_end != kw_end) conflict <- TRUE

  begin <- if (!is.na(hdr_begin) && !sentinel(hdr_begin)) hdr_begin else kw_begin
  end   <- if (!is.na(hdr_end)   && !sentinel(hdr_end))   hdr_end   else kw_end

  # Header zeroed and no keyword: for ANALYSIS that legitimately means "absent"
  # (the spec says to check the keywords to decide); the caller distinguishes.
  if (is.na(begin) && sentinel(hdr_begin)) begin <- 0
  if (is.na(end)   && sentinel(hdr_end))   end   <- 0

  list(begin = begin, end = end, style = style, conflict = conflict,
       from_keyword = is.na(hdr_begin) || sentinel(hdr_begin))
}

#' Work out where every segment lives, and therefore how much room TEXT has.
#'
#' The layout model is one number: `body_start`, the first byte of the first
#' segment that begins after TEXT. Everything in [body_start, EOF) is an opaque
#' tail that either stays exactly where it is or shifts as one block; nothing
#' in this file ever looks inside it.
fcs_resolve_offsets <- function(hdr, txt) {
  data <- .fcs_resolve_segment(hdr$datastart, hdr$dataend,
                               fcs_text_num(txt, "$BEGINDATA"),
                               fcs_text_num(txt, "$ENDDATA"))
  ana  <- .fcs_resolve_segment(hdr$anastart, hdr$anaend,
                               fcs_text_num(txt, "$BEGINANALYSIS"),
                               fcs_text_num(txt, "$ENDANALYSIS"))
  # Supplemental TEXT has no header field at all — it exists only as keywords.
  stext_begin <- fcs_text_num(txt, "$BEGINSTEXT", 0)
  stext_end   <- fcs_text_num(txt, "$ENDSTEXT", 0)
  if (is.na(stext_begin)) stext_begin <- 0
  if (is.na(stext_end))   stext_end   <- 0

  textstart <- hdr$textstart
  textend   <- hdr$textend

  starts <- c(data$begin, ana$begin, stext_begin)
  starts <- starts[!is.na(starts) & starts > 0 & starts > textend]
  body_start <- if (length(starts)) min(starts) else hdr$file_size

  ends <- c(textend, data$end, ana$end, stext_end)
  ends <- ends[!is.na(ends) & ends > 0]
  last_seg_end <- if (length(ends)) max(ends) else textend

  trailer_present <- isTRUE(hdr$file_size == last_seg_end + 9)

  list(textstart = textstart, textend = textend,
       text_len = textend - textstart + 1,
       datastart = data$begin, dataend = data$end,
       data_style = data$style, data_from_keyword = data$from_keyword,
       data_conflict = data$conflict,
       anastart = ana$begin, anaend = ana$end,
       ana_style = ana$style, ana_from_keyword = ana$from_keyword,
       ana_conflict = ana$conflict,
       stextstart = stext_begin, stextend = stext_end,
       body_start = body_start,
       room = body_start - textstart,
       gap = body_start - (textend + 1),
       last_seg_end = last_seg_end,
       trailer_present = trailer_present,
       file_size = hdr$file_size,
       nextdata = fcs_text_num(txt, "$NEXTDATA", 0))
}

# -----------------------------------------------------------------------------
# Validation — the refuse-to-touch list
# -----------------------------------------------------------------------------

# Structural, hardware and audit keywords. Rewriting any of these would make
# the file describe data it does not contain, so no caller can reach them.
FCS_PROTECTED_KEYWORDS <- c(
  "$TOT", "$PAR", "$DATATYPE", "$BYTEORD", "$MODE", "$NEXTDATA",
  "$BEGINDATA", "$ENDDATA", "$BEGINANALYSIS", "$ENDANALYSIS",
  "$BEGINSTEXT", "$ENDSTEXT",
  "$CYT", "$CYTSN", "$DATE", "$BTIM", "$ETIM", "$TIMESTEP",
  "ORIGINALGUID"
)

# $PnB/$PnN/$PnS/$PnE/$PnR/$PnG and friends, for any n.
.fcs_is_parameter_keyword <- function(key) {
  grepl("^\\$P[0-9]+[A-Z]+$", toupper(trimws(key)))
}

#' Everything that makes a file unsafe to edit, as a character vector of codes.
#'
#' A file that yields any code is left byte-identical. The codes are the user's
#' interface for the "why did you skip this one" question, so they are stable
#' strings, and `allow` can downgrade the handful that are judgement calls.
fcs_validate <- function(hdr, txt, seg, allow = character()) {
  codes <- character(0)
  add <- function(x) if (!(x %in% allow)) codes <<- c(codes, x)

  if (!grepl("^FCS3\\.[012]$", hdr$version)) add("unsupported_version")
  if (length(hdr$bad_fields)) add("bad_header_field")

  if (!is.na(txt$malformed)) add(txt$malformed)
  if (txt$n_pairs > 0L) {
    if (any(is.na(txt$keys_upper))) add("non_ascii_keyword")
    dup <- txt$keys_upper[!is.na(txt$keys_upper)]
    if (any(duplicated(dup))) add("duplicate_keyword")
  }

  d <- as.integer(txt$delimiter)
  if (length(d) == 1L && !is.na(d)) {
    # 2.2.15 allows ASCII 1-126, but a space or digit delimiter makes the
    # padding and numeric conventions ambiguous and every reader disagrees.
    if (d == 32L || (d >= 48L && d <= 57L) || d == 0L || d > 126L) {
      add("bad_delimiter")
    }
  }

  if (isTRUE(seg$data_conflict) || isTRUE(seg$ana_conflict)) {
    add("header_text_offset_conflict")
  }
  if (is.na(seg$datastart) || is.na(seg$dataend)) add("offsets_unknown")
  if (!is.na(seg$nextdata) && seg$nextdata != 0) add("nextdata_nonzero")
  if (!is.na(seg$stextstart) && seg$stextstart > 0) add("supplemental_text")

  if (!is.na(seg$datastart) && !is.na(seg$textstart) &&
      seg$datastart > 0 && seg$datastart < seg$textstart) add("text_after_data")
  if (!is.na(seg$textend) && seg$textend > 99999999) add("text_exceeds_limit")
  if (!is.na(seg$room) && seg$room < 1) add("no_room_for_text")

  # OTHER segments live between the header and TEXT and carry their own
  # absolute offsets. We never shift them, so they are tolerated only when they
  # are padding rather than a real segment table.
  if (!is.na(hdr$textstart) && hdr$textstart > 58) {
    other <- hdr$other_raw
    ints  <- as.integer(other)
    known <- !isTRUE(hdr$other_oversize) && length(other) > 0L

    # (a) Nothing but spaces and/or NULs. An OTHER offset is a right-justified
    # ASCII integer, so a region with no digits in it cannot encode one: every
    # 8-byte window reads as 0, which the spec defines as "no segment". This is
    # its own case rather than a subset of (b) because the region's LENGTH need
    # not be a multiple of 8 — BD Diva starts TEXT at 256 and space-pads bytes
    # 58..255, which is 198 bytes, and the pair-wise check below would reject a
    # cohort of perfectly ordinary files on the strength of that remainder
    # alone.
    blank_pad <- known && all(ints == 32L | ints == 0L)

    # (b) A well-formed table of 8-byte offset fields that are all zero, i.e.
    # the writer reserved the space and declared no segments.
    zero_pairs <- known && length(other) %% 8L == 0L &&
      all(vapply(seq_len(length(other) %/% 8L), function(i) {
        v <- .fcs_num_field(other[((i - 1L) * 8L + 1L):(i * 8L)])
        !is.na(v) && v == 0
      }, logical(1)))

    if (!isTRUE(blank_pad || zero_pairs)) add("other_segments")
  }

  if (!is.na(seg$last_seg_end) && hdr$file_size < seg$last_seg_end + 1) {
    add("truncated_file")
  } else if (!is.na(seg$last_seg_end) &&
             hdr$file_size > seg$last_seg_end + 9) {
    # More bytes than "last segment plus a CRC field" means something is there
    # that we cannot account for — most often a second data set whose
    # $NEXTDATA was never written.
    add("unexplained_trailing_bytes")
  }

  codes
}

# The 8-byte CRC field sits immediately after the last segment. A real CRC
# covers the HEADER and the TEXT segment, so any edit invalidates it — but
# essentially nobody computes one and the opt-out is eight ASCII zeros.
.fcs_read_trailer <- function(path, seg) {
  if (!isTRUE(seg$trailer_present)) return(NA_character_)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, seg$last_seg_end + 1)
  rawToChar(readBin(con, "raw", n = 8L))
}

# -----------------------------------------------------------------------------
# Inspection — one read, everything a caller needs to decide
# -----------------------------------------------------------------------------

#' Read a file's header, TEXT and layout and say whether it is safe to edit.
#'
#' @param path Path to an FCS file.
#' @param allow Refusal codes to downgrade to warnings.
#' @return list(ok, refusals, warnings, header, text, seg, trailer, ...).
#'   Never throws for a malformed file — that is what `refusals` is for.
fcs_inspect <- function(path, allow = character()) {
  out <- list(path = path, ok = FALSE, refusals = character(0),
              warnings = character(0))

  hdr <- tryCatch(fcs_read_header(path), error = function(e) e)
  if (inherits(hdr, "error")) {
    out$refusals <- "unreadable_file"
    out$message <- conditionMessage(hdr)
    return(out)
  }

  info <- file.info(path)
  if (isTRUE(info$isdir)) {
    out$refusals <- "not_a_file"
    return(out)
  }
  # Note: base R's file.info() does not expose st_nlink, so a hard link cannot
  # be detected without shelling out. An in-place patch of a hard-linked file
  # edits every name pointing at that inode, which is arguably what a hard link
  # means. Symlinks we can see, and we do not follow them.
  if (!is.na(Sys.readlink(path)) && nzchar(Sys.readlink(path))) {
    out$refusals <- c(out$refusals, "symlink")
  }

  ts <- hdr$textstart
  te <- hdr$textend

  # One vendor bug is worth absorbing rather than refusing: writers that treat
  # textend as EXCLUSIVE, so it equals datastart and the segment reads one byte
  # long. fcsparser carries the same workaround. Checked before the sanity
  # tests below, because such a file is not otherwise malformed.
  if (!is.na(ts) && !is.na(te) && !is.na(hdr$datastart) &&
      hdr$datastart > 0 && te == hdr$datastart && te > ts) {
    te <- te - 1
    hdr$textend <- te
    out$warnings <- c(out$warnings, "textend_off_by_one")
  }

  if (is.na(ts) || is.na(te) || ts < 58 || te < ts || te >= hdr$file_size) {
    out$refusals <- c(out$refusals, "bad_text_offsets")
    out$header <- hdr
    return(out)
  }

  con <- file(path, "rb")
  txt_raw <- tryCatch({
    seek(con, ts)
    readBin(con, "raw", n = as.integer(te - ts + 1))
  }, finally = close(con))

  txt <- fcs_parse_text(txt_raw)
  seg <- fcs_resolve_offsets(hdr, txt)

  out$header <- hdr
  out$text <- txt
  out$text_raw <- txt_raw
  out$seg <- seg
  out$refusals <- c(out$refusals, fcs_validate(hdr, txt, seg, allow))

  # The precondition witness. Everything downstream compares the file we write
  # against the structure we parsed; if the parse is not lossless for THIS
  # file's bytes then all of those comparisons are self-referential.
  if (is.na(txt$malformed) &&
      !identical(fcs_serialize_text(txt), txt_raw)) {
    out$refusals <- c(out$refusals, "text_roundtrip_failed")
  }

  trailer <- tryCatch(.fcs_read_trailer(path, seg), error = function(e) NA_character_)
  out$trailer <- trailer
  if (!is.na(trailer) && trailer != "00000000") {
    out$refusals <- c(out$refusals, "crc_present")
  }

  # Advisory only. $PnB may be "*" (variable width) and $MODE may not be list
  # mode, in which case the arithmetic does not apply — and an NA reaching an
  # if() would abort the whole folder rather than skip one check.
  out$data_length_ok <- NA
  mode <- fcs_text_get(txt, "$MODE", NA_character_)
  tot  <- fcs_text_num(txt, "$TOT")
  par  <- fcs_text_num(txt, "$PAR")
  if (!is.na(mode) && toupper(trimws(mode)) == "L" &&
      !is.na(tot) && !is.na(par) && par >= 1 &&
      !is.na(seg$datastart) && !is.na(seg$dataend)) {
    pnb <- vapply(seq_len(par), function(i) fcs_text_num(txt, paste0("$P", i, "B")),
                  numeric(1))
    if (all(is.finite(pnb))) {
      expected <- tot * sum(pnb) / 8
      actual <- seg$dataend - seg$datastart + 1
      out$data_length_ok <- isTRUE(expected == actual)
      out$data_length_expected <- expected
      out$data_length_actual <- actual
      if (!out$data_length_ok) {
        out$warnings <- c(out$warnings, "data_length_mismatch")
      }
    }
  }

  out$refusals <- setdiff(unique(out$refusals), allow)
  out$ok <- length(out$refusals) == 0L
  out
}

#' All keywords of a file as a named character vector. For reports only —
#' the writer never round-trips through character.
fcs_keywords <- function(path) {
  ins <- fcs_inspect(path, allow = character())
  if (is.null(ins$text) || ins$text$n_pairs == 0L) return(character(0))
  v <- vapply(ins$text$values_raw, .fcs_value_string, character(1))
  names(v) <- vapply(ins$text$keys_raw, function(k)
    .fcs_value_string(k), character(1))
  v
}

#' One keyword of a file, case-insensitively.
fcs_keyword <- function(path, key, default = NA_character_) {
  ins <- fcs_inspect(path, allow = character())
  if (is.null(ins$text)) return(default)
  fcs_text_get(ins$text, key, default)
}

# -----------------------------------------------------------------------------
# Layout planning
# -----------------------------------------------------------------------------

#' Choose a write strategy from three numbers. Pure — no I/O, no file needed.
#'
#' @param text_len_old Length in bytes of the existing TEXT segment.
#' @param room Bytes available from textstart up to the next segment.
#' @param natural_len Length in bytes the new TEXT serializes to.
#' @return list(strategy, target_len, fill, grow).
#'   `fill` is the number of space bytes written after the new TEXT to erase
#'   whatever used to be there; `grow` is how far DATA has to move.
fcs_plan_layout <- function(text_len_old, room, natural_len) {
  if (natural_len <= room) {
    strategy <- if (natural_len <= text_len_old) "neutral" else "patch"
    return(list(strategy = strategy, target_len = natural_len,
                fill = room - natural_len, grow = 0))
  }
  list(strategy = "rewrite", target_len = natural_len,
       fill = 0, grow = natural_len - room)
}

# Serialize the candidate TEXT with the six offset keywords carrying `shift`.
#
# This is where the circular dependency lives: the digits of $BEGINDATA are
# part of the length that determines $BEGINDATA. It only bites on the rewrite
# path, and it resolves by iteration to a fixed point below — the difference
# from flowCore's loop (issue #129) being that ours re-derives from scratch each
# round, so it converges when a value gets SHORTER too.
.fcs_shift_offsets <- function(txt, seg, shift) {
  set <- character(0)
  # The guard is on the ORIGINAL value, not the shifted one. A zero offset means
  # "this segment does not exist" (the spec has readers check $BEGINANALYSIS to
  # decide whether an ANALYSIS segment is present at all), and 0 + shift is a
  # perfectly plausible-looking offset — so shifting it would invent a segment
  # pointing into the middle of somebody's events.
  add <- function(key, value) {
    if (is.na(value) || value <= 0) return()
    if (!length(fcs_text_index(txt, key))) return()
    set[[key]] <<- sprintf("%.0f", value + shift)
  }
  add("$BEGINDATA", seg$datastart)
  add("$ENDDATA", seg$dataend)
  add("$BEGINANALYSIS", seg$anastart)
  add("$ENDANALYSIS", seg$anaend)
  if (length(set) == 0L) return(txt)
  fcs_text_set(txt, set = set)
}

# -----------------------------------------------------------------------------
# Verification helpers
# -----------------------------------------------------------------------------

#' Compare two byte ranges in two files without loading either into memory.
#'
#' Used to prove the event data survived. It reads the WHOLE range, not a
#' sample: a 4 KB spot check cannot see an overrun at byte 2e9, which is exactly
#' the failure this is here to catch.
fcs_range_identical <- function(path_a, from_a, path_b, from_b, n,
                                chunk_bytes = 8 * 1024^2) {
  if (n <= 0) return(TRUE)
  ca <- file(path_a, "rb"); cb <- file(path_b, "rb")
  on.exit({ close(ca); close(cb) }, add = TRUE)
  seek(ca, from_a)
  seek(cb, from_b)
  left <- n
  while (left > 0) {
    k <- as.integer(min(chunk_bytes, left))
    a <- readBin(ca, "raw", n = k)
    b <- readBin(cb, "raw", n = k)
    if (length(a) != k || length(b) != k || !identical(a, b)) return(FALSE)
    left <- left - k
  }
  TRUE
}

# Ordered (key, value) byte sequence — the structural witness. Comparing this
# as raw catches a value that decoded and re-encoded differently, which a
# character comparison would not.
.fcs_pair_signature <- function(txt) {
  lapply(seq_len(txt$n_pairs), function(i)
    list(k = txt$keys_raw[[i]], v = txt$values_raw[[i]]))
}

# -----------------------------------------------------------------------------
# The writer
# -----------------------------------------------------------------------------

# Patch the six offset fields of a 58-byte header, preserving each field's
# original convention: a field written as 0 stays 0, a field written as
# 99999999 stays 99999999, and only a field that carried a real value is
# updated. The spec zeroes offsets PER SEGMENT, so when one member of a pair
# overflows the 8 digits both members go to zero together.
.fcs_patch_header <- function(hdr_raw, textend, seg, shift) {
  out <- hdr_raw
  put <- function(ix, value) {
    f <- .fcs_field8(value)
    if (is.null(f)) return(FALSE)
    out[ix] <<- f
    TRUE
  }
  zero <- function(ix) out[ix] <<- charToRaw("       0")

  if (!put(19:26, textend)) return(NULL)   # TEXT must fit; validator enforced it

  pair <- function(ix_b, ix_e, begin, end, style) {
    if (identical(style, "zero")) { zero(ix_b); zero(ix_e); return(TRUE) }
    if (identical(style, "sentinel")) return(TRUE)   # leave 99999999 in place
    if (is.na(begin) || is.na(end) || begin <= 0) {
      zero(ix_b); zero(ix_e); return(TRUE)
    }
    if (begin + shift > 99999999 || end + shift > 99999999) {
      zero(ix_b); zero(ix_e)
      return(NA)   # caller must confirm the TEXT keywords carry the truth
    }
    put(ix_b, begin + shift) && put(ix_e, end + shift)
  }

  ok_data <- pair(27:34, 35:42, seg$datastart, seg$dataend, seg$data_style)
  ok_ana  <- pair(43:50, 51:58, seg$anastart, seg$anaend, seg$ana_style)
  list(raw = out, data_overflow = is.na(ok_data), ana_overflow = is.na(ok_ana))
}

#' Rewrite keywords in an FCS file without decoding the event data.
#'
#' @param path Source file.
#' @param set Named character vector of keyword = value to write.
#' @param unset Character vector of keywords to remove.
#' @param output_path Where to write. NULL means patch `path` itself.
#' @param allow Refusal codes to downgrade.
#' @param dry_run Plan and report, write nothing.
#' @param allow_parameter_labels Let $PnS (the marker label of parameter n)
#'   through the parameter-keyword guard. $PnS is a free-text display name with
#'   no structural role, unlike $PnB/$PnN/$PnE/$PnR, which stay protected. Off
#'   by default; only the marker-name harmoniser turns it on.
#' @param .fault Test-only fault injection; NULL in every real call.
#' @return list(status = ok|unchanged|refused|failed, strategy, changed, ...).
fcs_write_keywords <- function(path, set = character(), unset = character(),
                               output_path = NULL, allow = character(),
                               dry_run = FALSE, chunk_bytes = 8 * 1024^2,
                               allow_parameter_labels = FALSE,
                               .fault = NULL) {
  res <- list(path = path, output_path = output_path, status = "failed",
              strategy = NA_character_, changed = character(0),
              refusals = character(0), warnings = character(0),
              text_bytes_old = NA_real_, text_bytes_new = NA_real_,
              room = NA_real_, grow = 0, fill = NA_real_,
              data_identical = NA, message = "")

  # Collapse a repeated keyword to its last value before anything looks at it.
  # Otherwise the write applies both (last wins) while the arithmetic witness
  # counts the delta twice, and the two disagree by exactly the difference —
  # a verification failure on a write that was actually fine.
  if (length(set) && !is.null(names(set))) {
    keep <- !duplicated(toupper(trimws(names(set))), fromLast = TRUE)
    set <- set[keep]
  }
  unset <- unique(unset)

  ins <- fcs_inspect(path, allow = allow)
  res$warnings <- ins$warnings
  if (!ins$ok) {
    res$status <- "refused"
    res$refusals <- ins$refusals
    res$message <- paste("refused:", paste(ins$refusals, collapse = ", "))
    return(res)
  }

  hdr <- ins$header; txt <- ins$text; seg <- ins$seg
  res$room <- seg$room
  res$text_bytes_old <- seg$text_len

  # --- what is actually changing -------------------------------------------
  protected <- toupper(trimws(c(names(set), unset)))
  is_param <- vapply(protected, .fcs_is_parameter_keyword, logical(1))
  if (allow_parameter_labels) is_param <- is_param & !grepl("^\\$P[0-9]+S$", protected)
  bad <- protected[protected %in% FCS_PROTECTED_KEYWORDS | is_param]
  if (length(bad)) {
    res$status <- "refused"
    res$refusals <- "protected_keyword"
    res$message <- paste("refusing to write protected keyword(s):",
                         paste(unique(bad), collapse = ", "))
    return(res)
  }

  delim_chr <- rawToChar(txt$delimiter)
  for (nm in names(set)) {
    v <- as.character(set[[nm]])
    if (grepl(delim_chr, v, fixed = TRUE)) {
      res$status <- "refused"
      res$refusals <- "new_value_contains_delimiter"
      res$message <- sprintf(
        paste0("value for %s contains this file's delimiter '%s'. ",
               "Rename the file instead — writing an escaped delimiter would ",
               "break readers that split on it."), nm, delim_chr)
      return(res)
    }
    b <- as.integer(charToRaw(enc2utf8(v)))
    if (length(b) == 0L || any(b < 32L) || any(b == 127L)) {
      res$status <- "refused"
      res$refusals <- "new_value_unprintable"
      res$message <- sprintf("value for %s is empty or has control bytes", nm)
      return(res)
    }
  }

  changed <- character(0)
  for (nm in names(set)) {
    old <- fcs_text_get(txt, nm, NA_character_)
    new <- as.character(set[[nm]])
    if (is.na(old) || !identical(old, new)) changed[[nm]] <- new
  }
  removed <- unset[vapply(unset, function(k)
    length(fcs_text_index(txt, k)) > 0L, logical(1))]

  if (length(changed) == 0L && length(removed) == 0L) {
    res$status <- "unchanged"
    res$strategy <- "none"
    res$text_bytes_new <- seg$text_len
    res$data_identical <- TRUE
    res$message <- "already in sync"
    return(res)
  }
  res$changed <- changed

  # --- plan -----------------------------------------------------------------
  new_txt <- fcs_text_set(txt, set = set, unset = unset)
  natural <- length(fcs_serialize_text(new_txt))
  plan <- fcs_plan_layout(seg$text_len, seg$room, natural)

  shift <- 0
  if (plan$strategy == "rewrite") {
    # Fixed point: the shift changes the digit count of the offset keywords,
    # which changes the length, which changes the shift. Bounded and monotone.
    converged <- FALSE
    for (i in 1:16) {
      cand <- .fcs_shift_offsets(new_txt, seg, shift)
      len <- length(fcs_serialize_text(cand))
      new_body <- seg$textstart + len + seg$gap
      new_shift <- new_body - seg$body_start
      if (new_shift == shift) { converged <- TRUE; break }
      shift <- new_shift
    }
    if (!converged) {
      res$status <- "failed"
      res$message <- "offset layout did not converge"
      return(res)
    }
    new_txt <- .fcs_shift_offsets(new_txt, seg, shift)
    natural <- length(fcs_serialize_text(new_txt))
    plan <- list(strategy = "rewrite", target_len = natural, fill = 0,
                 grow = shift)
  }

  res$strategy <- plan$strategy
  res$text_bytes_new <- natural
  res$fill <- plan$fill
  res$grow <- plan$grow

  # Arithmetic witness, computed without consulting the serializer: the byte
  # delta of the segment must equal the sum of the per-keyword byte deltas.
  delta <- 0
  for (nm in names(set)) {
    i <- fcs_text_index(txt, nm)
    new_bytes <- .fcs_nbytes(enc2utf8(as.character(set[[nm]])))
    if (length(i)) {
      delta <- delta + new_bytes - length(txt$values_raw[[i[1]]])
    } else {
      delta <- delta + new_bytes + .fcs_nbytes(enc2utf8(nm)) + 2  # + 2 delimiters
    }
  }
  for (k in removed) {
    i <- fcs_text_index(txt, k)
    delta <- delta - (length(txt$keys_raw[[i[1]]]) +
                        length(txt$values_raw[[i[1]]]) + 2)
  }
  res$expected_delta <- delta

  if (dry_run) {
    res$status <- "ok"
    res$data_identical <- TRUE
    res$message <- sprintf("DRY RUN: %s, TEXT %.0f -> %.0f bytes",
                           plan$strategy, seg$text_len, natural)
    return(res)
  }

  new_text_raw <- fcs_serialize_text(new_txt)
  in_place <- is.null(output_path) || identical(normalizePath(output_path,
                                                             mustWork = FALSE),
                                                normalizePath(path,
                                                              mustWork = FALSE))
  dest <- if (in_place) path else output_path

  out <- tryCatch(
    if (plan$strategy == "rewrite" || !in_place) {
      .fcs_write_rewrite(path, dest, hdr, seg, new_text_raw, plan, shift,
                         in_place, chunk_bytes, .fault)
    } else {
      .fcs_write_patch(path, hdr, seg, new_text_raw, plan, .fault)
    },
    error = function(e) e)

  if (inherits(out, "error")) {
    res$status <- "failed"
    res$message <- conditionMessage(out)
    return(res)
  }
  res$output_path <- if (!is.null(out$final_dest)) out$final_dest else out$dest
  res$warnings <- c(res$warnings, out$warnings)

  # --- verify ---------------------------------------------------------------
  v <- .fcs_verify(path, out$dest, ins, new_txt, plan, shift, in_place,
                   chunk_bytes, expected_delta = delta)
  res$data_identical <- v$data_identical
  if (!v$ok) {
    if (!is.null(out$restore)) out$restore()
    res$status <- "failed"
    res$message <- paste("verification failed:", v$message,
                         if (!is.null(out$restore)) "(original restored)" else "")
    return(res)
  }
  if (!is.null(out$commit)) out$commit()

  res$status <- "ok"
  res$message <- sprintf("%s: TEXT %.0f -> %.0f bytes%s", plan$strategy,
                         seg$text_len, natural,
                         if (plan$grow > 0) sprintf(", DATA moved +%.0f",
                                                    plan$grow) else "")
  res
}

# --- in-place patch (neutral / patch) ----------------------------------------
#
# One contiguous write of exactly `room` bytes at textstart: the new TEXT
# followed by spaces. The spaces matter — 2.2.11 requires unused space between
# segments to be filled with ASCII 32, and writing them is what guarantees no
# stale fragment of the old TEXT is left where a lenient parser could reach it.
#
# The crash-recovery record is a sidecar holding just [0, body_start) — a
# kilobyte or so, not a copy of the events — written and closed BEFORE the file
# is opened for writing. R exposes no fsync, so a power cut mid-write can still
# tear the segment; this is what makes that recoverable.
.fcs_write_patch <- function(path, hdr, seg, new_text_raw, plan, .fault = NULL) {
  room <- seg$room
  stopifnot(length(new_text_raw) <= room)
  textend_new <- seg$textstart + length(new_text_raw) - 1

  patched <- .fcs_patch_header(hdr$raw, textend_new, seg, 0)
  if (is.null(patched)) stop("textend does not fit the 8-byte header field")

  sidecar <- paste0(path, ".fcsprefix")
  con <- file(path, "rb")
  prefix <- tryCatch({
    seek(con, 0)
    readBin(con, "raw", n = as.integer(seg$body_start))
  }, finally = close(con))
  saveRDS(list(path = path, body_start = seg$body_start,
               file_size = seg$file_size, prefix = prefix), sidecar)

  buf <- c(new_text_raw, rep(as.raw(0x20), room - length(new_text_raw)))
  if (identical(.fault, "corrupt_prefix")) buf[1] <- as.raw(0x00)

  con <- file(path, "r+b")
  tryCatch({
    .fcs_seek_write(con, 0)
    writeBin(patched$raw, con)
    .fcs_seek_write(con, seg$textstart)
    writeBin(buf, con)
  }, finally = close(con))

  restore <- function() {
    rec <- readRDS(sidecar)
    c2 <- file(rec$path, "r+b")
    tryCatch({
      .fcs_seek_write(c2, 0)
      writeBin(rec$prefix, c2)
    }, finally = close(c2))
    unlink(sidecar)
  }
  list(dest = path, warnings = character(0), restore = restore,
       commit = function() unlink(sidecar))
}

# --- full rewrite -------------------------------------------------------------
#
# Used when TEXT no longer fits, and also whenever the caller asked for a copy
# into an output directory. The tail is copied with readBin/writeBin in chunks
# so a 4 GB file never lands in memory.
.fcs_write_rewrite <- function(path, dest, hdr, seg, new_text_raw, plan, shift,
                               in_place, chunk_bytes, .fault = NULL) {
  warnings <- character(0)
  textend_new <- seg$textstart + length(new_text_raw) - 1
  patched <- .fcs_patch_header(hdr$raw, textend_new, seg, shift)
  if (is.null(patched)) stop("textend does not fit the 8-byte header field")

  # If an offset outgrew its 8 digits the header field goes to zero, which is
  # only legal if the TEXT keyword carries the real value.
  if (isTRUE(patched$data_overflow)) {
    if (!length(fcs_text_index(fcs_parse_text(new_text_raw), "$BEGINDATA"))) {
      stop("DATA offset exceeds 8 digits and $BEGINDATA is absent")
    }
    warnings <- c(warnings, "data_offsets_moved_to_keywords")
  }

  tmp <- paste0(dest, ".fcstmp", Sys.getpid())
  tail_bytes <- seg$file_size - seg$body_start

  # The tail has to land at exactly body_start + grow, so the space fill is
  # whatever is left over between the end of the new TEXT and that point. For a
  # true rewrite that works out to the file's original gap; for a copy that did
  # not need to grow it is the leftover room. One formula covers both.
  pad <- seg$body_start + plan$grow - (seg$textstart + length(new_text_raw))
  if (pad < 0) stop("internal: new TEXT overruns the data segment")

  ci <- file(path, "rb")
  co <- file(tmp, "wb")
  tryCatch({
    writeBin(patched$raw, co)
    if (seg$textstart > 58) writeBin(hdr$other_raw, co)
    writeBin(new_text_raw, co)
    if (pad > 0) writeBin(rep(as.raw(0x20), pad), co)
    seek(ci, seg$body_start)
    left <- tail_bytes
    while (left > 0) {
      k <- as.integer(min(chunk_bytes, left))
      chunk <- readBin(ci, "raw", n = k)
      if (length(chunk) != k) stop("short read copying the data segment")
      if (identical(.fault, "corrupt_tail") && left == tail_bytes) {
        chunk[1] <- as.raw(if (chunk[1] == as.raw(0)) 1 else 0)
      }
      writeBin(chunk, co)
      left <- left - k
    }
  }, finally = { close(ci); close(co) })

  mode <- file.info(path)$mode
  if (!is.na(mode)) suppressWarnings(Sys.chmod(tmp, mode))

  old <- paste0(dest, ".fcsold", Sys.getpid())
  commit <- function() {
    if (file.exists(dest)) {
      if (!file.rename(dest, old)) stop("could not move the existing file aside")
    }
    if (!file.rename(tmp, dest)) {
      if (file.exists(old)) file.rename(old, dest)
      stop("could not put the new file in place")
    }
    unlink(old)
  }
  # Verification runs against the temp file, so "restore" is just cleanup.
  restore <- function() unlink(tmp)

  list(dest = tmp, final_dest = dest, warnings = warnings,
       restore = restore,
       commit = function() { commit(); invisible(TRUE) })
}

# --- the five witnesses -------------------------------------------------------
.fcs_verify <- function(orig_path, cand_path, ins, new_txt, plan, shift,
                        in_place, chunk_bytes, expected_delta = NULL) {
  fail <- function(msg) list(ok = FALSE, message = msg, data_identical = FALSE)
  seg <- ins$seg

  cand <- fcs_inspect(cand_path, allow = "unexplained_trailing_bytes")
  if (is.null(cand$text) || !is.na(cand$text$malformed)) {
    return(fail("the written TEXT segment does not re-parse"))
  }

  # (3) structural: the ordered (key, value) byte sequence must be what we planned
  if (!identical(.fcs_pair_signature(cand$text), .fcs_pair_signature(new_txt))) {
    return(fail("keyword sequence differs from the plan"))
  }

  # (2) arithmetic: the length must match the plan, and — where no offset
  # keyword was allowed to change width — the independent per-keyword sum too.
  if (!identical(as.numeric(cand$seg$text_len), as.numeric(plan$target_len))) {
    return(fail(sprintf("TEXT length is %.0f, planned %.0f",
                        cand$seg$text_len, plan$target_len)))
  }
  if (plan$grow == 0 && !is.null(expected_delta)) {
    if (!isTRUE(cand$seg$text_len == seg$text_len + expected_delta)) {
      return(fail(sprintf(
        "TEXT grew by %.0f bytes, but the keywords account for %.0f",
        cand$seg$text_len - seg$text_len, expected_delta)))
    }
  }

  # (4) range: every byte of [body_start, EOF) survived, and the file grew by
  # exactly `grow`. The whole range, not a sample.
  tail_bytes <- seg$file_size - seg$body_start
  same <- fcs_range_identical(orig_path, seg$body_start,
                              cand_path, seg$body_start + plan$grow,
                              tail_bytes, chunk_bytes)
  if (!same) return(fail("the data segment changed"))
  expected_size <- seg$file_size + plan$grow
  if (!isTRUE(cand$header$file_size == expected_size)) {
    return(list(ok = FALSE, data_identical = TRUE,
                message = sprintf("file is %.0f bytes, expected %.0f",
                                  cand$header$file_size, expected_size)))
  }

  # (5) independent: the data segment still has its original extent, and the
  # header and TEXT still agree wherever both carry a value.
  if (!is.na(cand$seg$datastart) && !is.na(seg$datastart)) {
    old_len <- seg$dataend - seg$datastart + 1
    new_len <- cand$seg$dataend - cand$seg$datastart + 1
    if (!isTRUE(old_len == new_len)) {
      return(list(ok = FALSE, data_identical = TRUE,
                  message = "the recorded data segment length changed"))
    }
    if (!isTRUE(cand$seg$datastart == seg$datastart + plan$grow)) {
      return(list(ok = FALSE, data_identical = TRUE,
                  message = "the recorded data offset is not where it was moved to"))
    }
  }
  if (isTRUE(cand$seg$data_conflict) || isTRUE(cand$seg$ana_conflict)) {
    return(list(ok = FALSE, data_identical = TRUE,
                message = "header and TEXT offsets disagree after the write"))
  }
  if (identical(ins$data_length_ok, TRUE) &&
      identical(cand$data_length_ok, FALSE)) {
    return(list(ok = FALSE, data_identical = TRUE,
                message = "$TOT x $PnB no longer matches the data segment"))
  }

  list(ok = TRUE, message = "", data_identical = TRUE)
}
