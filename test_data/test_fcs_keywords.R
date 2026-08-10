#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_fcs_keywords.R — the FCS keyword read/write layer
# =============================================================================
# Covers src/R_shared/fcs_keywords.R: header parsing, the TEXT parser and its
# doubled-delimiter escape, segment resolution across the header/keyword split,
# the three write strategies, every refusal code, and the post-write witnesses.
#
# These need NEITHER flowCore NOR CytoML NOR optparse — the fixtures are built
# byte by byte by test_data/fcs_fixtures.R — so they run on any machine with
# base R. That is the point: the layer whose job is to repair files that the
# cytometry stack rejects must be testable without the cytometry stack.
#
# Run from the project root:
#   Rscript test_data/test_fcs_keywords.R
# =============================================================================

cat("=== FCS keyword layer tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_fcs_keywords.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
source(file.path(base_dir, "test_data", "fcs_fixtures.R"))

# --- tiny test harness -------------------------------------------------------
.pass <- 0L
.fail <- 0L

ok <- function(label, cond) {
  if (isTRUE(cond)) {
    .pass <<- .pass + 1L
    cat("  PASS  ", label, "\n", sep = "")
  } else {
    .fail <<- .fail + 1L
    cat("  FAIL  ", label, "\n", sep = "")
  }
}

eq <- function(label, actual, expected) {
  same <- isTRUE(all.equal(actual, expected))
  if (!same) {
    cat("  FAIL  ", label, "\n",
        "        expected: ", paste(format(expected), collapse = " | "), "\n",
        "        actual:   ", paste(format(actual), collapse = " | "), "\n",
        sep = "")
    .fail <<- .fail + 1L
  } else {
    .pass <<- .pass + 1L
    cat("  PASS  ", label, "\n", sep = "")
  }
}

throws <- function(label, expr, pattern = NULL) {
  msg <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(msg)) {
    ok(label, FALSE)
  } else if (!is.null(pattern) && !grepl(pattern, msg)) {
    cat("  FAIL  ", label, "\n        message did not match '", pattern,
        "':\n        ", msg, "\n", sep = "")
    .fail <<- .fail + 1L
  } else {
    ok(label, TRUE)
  }
}

root <- file.path(tempdir(), paste0("fcskw_", Sys.getpid()))
unlink(root, recursive = TRUE)
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
fx <- function(name) file.path(root, name)

# Refuses, and leaves the file exactly as it found it.
refuses <- function(label, path, code, set = c("$FIL" = "new.fcs")) {
  before <- fx_md5(path)
  ins <- fcs_inspect(path)
  got_inspect <- code %in% ins$refusals
  r <- fcs_write_keywords(path, set = set)
  unchanged <- identical(fx_md5(path), before)
  if (!got_inspect || r$status != "refused" || !unchanged) {
    cat("  FAIL  ", label, "\n",
        "        refusals: ", paste(ins$refusals, collapse = ", "), "\n",
        "        status:   ", r$status, "  untouched: ", unchanged, "\n", sep = "")
    .fail <<- .fail + 1L
  } else {
    ok(label, TRUE)
  }
}

# --- the real fixtures --------------------------------------------------------
cat("\n-- committed fixtures ---------------------------------------------\n")
real <- file.path(base_dir, "test_data", "fcs_samples", "sample_001.fcs")
if (file.exists(real)) {
  h <- fcs_read_header(real)
  eq("real: version", h$version, "FCS3.0")
  eq("real: header offsets", c(h$textstart, h$textend, h$datastart, h$dataend,
                               h$anastart, h$anaend),
     c(58, 801, 802, 2001, 0, 0))
  ok("real: offsets are doubles",
     all(vapply(list(h$textstart, h$dataend), typeof, character(1)) == "double"))

  ins <- fcs_inspect(real)
  # The delimiter is '|' (0x7C). An early note in this work said form feed
  # (0x0C); a parser that hard-coded that would read the whole segment as one
  # keyword, so this assertion is here to stay.
  eq("real: delimiter is | not form feed", ins$text$delimiter, as.raw(0x7c))
  eq("real: keyword pairs", ins$text$n_pairs, 58L)
  eq("real: room", ins$seg$room, 744)
  eq("real: TEXT length", ins$seg$text_len, 744)
  ok("real: no $FIL", is.na(fcs_text_get(ins$text, "$FIL")))
  eq("real: GUID", fcs_text_get(ins$text, "GUID"), "file6385e7193252e")
  ok("real: TEXT round-trips byte-for-byte",
     identical(fcs_serialize_text(ins$text), ins$text_raw))
  eq("real: nothing to refuse", ins$refusals, character(0))
  eq("real: CRC field is the opt-out", ins$trailer, "00000000")
  ok("real: $TOT x $PnB matches the data segment", isTRUE(ins$data_length_ok))
  eq("real: fcs_keyword reads $TOT", fcs_keyword(real, "$TOT"), "50")
  eq("real: keyword lookup is case-insensitive",
     fcs_keyword(real, "$tot"), "50")

  # --- acceptance: the real files, renamed, gaining a $FIL they never had ----
  acc <- fx("acc"); dir.create(acc)
  for (i in 1:3) {
    src <- file.path(base_dir, "test_data", "fcs_samples",
                     sprintf("sample_00%d.fcs", i))
    file.copy(src, file.path(acc, sprintf("Patient%s.fcs", LETTERS[i])))
  }
  target <- file.path(acc, "PatientA.fcs")
  before <- fcs_inspect(target)
  data0 <- fx_bytes(target, before$seg$datastart,
                    before$seg$dataend - before$seg$datastart + 1)
  r <- fcs_write_keywords(target, set = c("$FIL" = "PatientA.fcs"))
  eq("acceptance: status", r$status, "ok")
  eq("acceptance: zero slack forces a rewrite", r$strategy, "rewrite")
  after <- fcs_inspect(target)
  eq("acceptance: $FIL written", fcs_text_get(after$text, "$FIL"), "PatientA.fcs")
  ok("acceptance: 1200 data bytes are byte-identical",
     identical(fx_bytes(target, after$seg$datastart, 1200), data0))
  eq("acceptance: data segment length unchanged",
     after$seg$dataend - after$seg$datastart + 1, 1200)
  eq("acceptance: DATA moved by exactly the growth",
     after$seg$datastart - before$seg$datastart, r$grow)
  eq("acceptance: header and $BEGINDATA agree",
     fcs_text_num(after$text, "$BEGINDATA"), after$header$datastart)
  eq("acceptance: header and $ENDDATA agree",
     fcs_text_num(after$text, "$ENDDATA"), after$header$dataend)
  eq("acceptance: file grew by the same amount",
     after$header$file_size - before$header$file_size, r$grow)
  eq("acceptance: CRC field still trails the data", after$trailer, "00000000")
  eq("acceptance: still clean", after$refusals, character(0))
  # Everything except the two offset keywords that HAVE to move is bit-identical
  # — same bytes, same order, same position. $BEGINANALYSIS/$ENDANALYSIS are in
  # that list on purpose: they are 0 here, meaning "no ANALYSIS segment", and a
  # shifted zero would invent a segment pointing into the events.
  moved <- c("$BEGINDATA", "$ENDDATA")
  sig <- function(t, n) lapply(seq_len(n), function(i) {
    if (t$keys_upper[i] %in% moved) return(list(k = t$keys_raw[[i]], v = "moved"))
    list(k = t$keys_raw[[i]], v = t$values_raw[[i]])
  })
  ok("acceptance: every other keyword survives bit-identical, in order",
     identical(sig(after$text, before$text$n_pairs),
               sig(before$text, before$text$n_pairs)))
  eq("acceptance: an absent ANALYSIS segment stays absent",
     c(fcs_text_get(after$text, "$BEGINANALYSIS"),
       fcs_text_get(after$text, "$ENDANALYSIS")), c("0", "0"))
  eq("acceptance: exactly one keyword added",
     after$text$n_pairs - before$text$n_pairs, 1L)
} else {
  cat("  SKIP  committed fixtures not present\n")
}

# --- delimiters and escaping --------------------------------------------------
cat("\n-- delimiters and escaping ----------------------------------------\n")
for (d in list(c("\f", "form feed"), c("/", "slash"), c("!", "bang"))) {
  p <- fx(paste0("delim_", utf8ToInt(d[1]), ".fcs"))
  fcs_fixture(p, delimiter = d[1], keywords = list("$FIL" = "old.fcs"))
  ins <- fcs_inspect(p)
  eq(paste0("delimiter read from the file: ", d[2]),
     ins$text$delimiter, charToRaw(d[1]))
  ok(paste0("round-trips: ", d[2]),
     identical(fcs_serialize_text(ins$text), ins$text_raw))
  r <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"))
  eq(paste0("writes through: ", d[2]), fcs_keyword(p, "$FIL"), "new.fcs")
}

# 3.2.7's own example: $SYS/RSX-11//M/ is the value RSX-11/M
p <- fx("escape.fcs")
fcs_fixture(p, delimiter = "/", keywords = list("$SYS" = "RSX-11/M",
                                                "$FIL" = "old.fcs"))
ins <- fcs_inspect(p)
eq("doubled delimiter decodes to one literal", fcs_text_get(ins$text, "$SYS"),
   "RSX-11/M")
ok("escaped value round-trips byte-exact",
   identical(fcs_serialize_text(ins$text), ins$text_raw))
r <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"))
eq("an escaped neighbour survives a write", fcs_keyword(p, "$SYS"), "RSX-11/M")
eq("...and the file is still clean", fcs_inspect(p)$refusals, character(0))

p <- fx("escape_tail.fcs")
fcs_fixture(p, delimiter = "/", keywords = list("$OP" = "ends with/",
                                                "$FIL" = "old.fcs"))
ins <- fcs_inspect(p)
eq("a value ending in the delimiter decodes", fcs_text_get(ins$text, "$OP"),
   "ends with/")
ok("...and round-trips", identical(fcs_serialize_text(ins$text), ins$text_raw))

# --- non-ASCII bytes ----------------------------------------------------------
cat("\n-- byte fidelity ---------------------------------------------------\n")
p <- fx("latin1.fcs")
mu <- rawToChar(as.raw(c(0x43, 0x44, 0x33, 0x20, 0xB5, 0x67)))  # "CD3 <b5>g"
Encoding(mu) <- "latin1"
fcs_fixture(p, keywords = list("$P1S" = mu, "$FIL" = "old.fcs"))
raw_before <- readBin(p, "raw", n = file.size(p))
ins <- fcs_inspect(p)
ok("latin-1 bytes survive the parser",
   identical(fcs_serialize_text(ins$text), ins$text_raw))
i1 <- fcs_text_index(ins$text, "$P1S")
ok("the 0xB5 byte is still 0xB5",
   any(ins$text$values_raw[[i1]] == as.raw(0xB5)))
r <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"))
eq("a write next to non-ASCII succeeds", r$status, "ok")
after <- fcs_inspect(p)
ok("the 0xB5 byte is untouched after the write",
   identical(after$text$values_raw[[fcs_text_index(after$text, "$P1S")]],
             ins$text$values_raw[[i1]]))

# --- strategy selection -------------------------------------------------------
cat("\n-- strategy selection ----------------------------------------------\n")
eq("planner: shrink is neutral",
   fcs_plan_layout(500, 500, 480)$strategy, "neutral")
eq("planner: growth into slack is a patch",
   fcs_plan_layout(500, 628, 520)$strategy, "patch")
eq("planner: growth past the room is a rewrite",
   fcs_plan_layout(500, 504, 520)$strategy, "rewrite")
eq("planner: the rewrite grows by room-shortfall",
   fcs_plan_layout(500, 504, 520)$grow, 16)
eq("planner: fill erases the rest of the room",
   fcs_plan_layout(500, 500, 480)$fill, 20)

# neutral — DATA never moves, the file never changes size
p <- fx("neutral.fcs")
fcs_fixture(p, keywords = list("$FIL" = "a_very_long_original_name.fcs"))
i0 <- fcs_inspect(p)
size0 <- i0$header$file_size
hdr0  <- fx_bytes(p, 0, 58)
data0 <- fx_bytes(p, i0$seg$datastart, 40)
bd0   <- fcs_text_get(i0$text, "$BEGINDATA")
r <- fcs_write_keywords(p, set = c("$FIL" = "s.fcs"))
eq("neutral: strategy", r$strategy, "neutral")
i1 <- fcs_inspect(p)
eq("neutral: file size unchanged", i1$header$file_size, size0)
eq("neutral: datastart unchanged", i1$seg$datastart, i0$seg$datastart)
ok("neutral: data bytes unchanged", identical(fx_bytes(p, i1$seg$datastart, 40), data0))
eq("neutral: $BEGINDATA never touched", fcs_text_get(i1$text, "$BEGINDATA"), bd0)
ok("neutral: only textend moved in the header",
   identical(fx_bytes(p, 0, 18), hdr0[1:18]) &&
     identical(fx_bytes(p, 26, 32), hdr0[27:58]) &&
     !identical(fx_bytes(p, 18, 8), hdr0[19:26]))
eq("neutral: TEXT still ends on its delimiter",
   fx_bytes(p, i1$seg$textend, 1), i1$text$delimiter)
ok("neutral: the freed bytes are spaces, not stale TEXT",
   all(fx_bytes(p, i1$seg$textend + 1,
                i1$seg$body_start - i1$seg$textend - 1) == as.raw(0x20)))

# patch — grows into slack, DATA still never moves
p <- fx("patch.fcs")
fcs_fixture(p, keywords = list("$FIL" = "s.fcs"), slack = 128)
i0 <- fcs_inspect(p)
data0 <- fx_bytes(p, i0$seg$datastart, 40)
r <- fcs_write_keywords(p, set = c("$FIL" = "a_considerably_longer_name.fcs"))
eq("patch: strategy", r$strategy, "patch")
i1 <- fcs_inspect(p)
eq("patch: datastart unchanged", i1$seg$datastart, i0$seg$datastart)
eq("patch: file size unchanged", i1$header$file_size, i0$header$file_size)
ok("patch: data bytes unchanged", identical(fx_bytes(p, i1$seg$datastart, 40), data0))
eq("patch: textend follows the longer TEXT",
   i1$seg$textend, i1$seg$textstart + i1$seg$text_len - 1)
ok("patch: leftover slack is spaces",
   all(fx_bytes(p, i1$seg$textend + 1,
                i1$seg$body_start - i1$seg$textend - 1) == as.raw(0x20)))

# An in-place patch on a file whose TEXT does NOT start at byte 58.
#
# Regression test with teeth: an "r+b" connection keeps separate read and write
# positions, and seek() moves the read one by default. With textstart == 58 the
# resulting misplaced write lands in the right spot by accident — every other
# in-place test here would pass with the bug present. This one would not.
p <- fx("offset_patch.fcs")
fcs_fixture(p, keywords = list("$FIL" = "a_long_original_name.fcs"),
            other_bytes = 16L, slack = 64)
i0 <- fcs_inspect(p)
eq("shifted TEXT: starts past the header", i0$seg$textstart, 74)
data0 <- fx_bytes(p, i0$seg$datastart, 40)
r <- fcs_write_keywords(p, set = c("$FIL" = "short.fcs"))
eq("shifted TEXT: patched in place", r$strategy, "neutral")
i1 <- fcs_inspect(p)
eq("shifted TEXT: the new value landed at textstart, not at byte 58",
   fcs_text_get(i1$text, "$FIL"), "short.fcs")
eq("shifted TEXT: version string intact", i1$header$version, "FCS3.0")
eq("shifted TEXT: OTHER padding intact", i1$seg$textstart, 74)
ok("shifted TEXT: data untouched", identical(fx_bytes(p, i1$seg$datastart, 40), data0))
eq("shifted TEXT: file is still clean", i1$refusals, character(0))

# room, not textend, is the threshold
p <- fx("just_short.fcs")
fcs_fixture(p, keywords = list("$FIL" = "s.fcs"), slack = 4)
r <- fcs_write_keywords(p, set = c("$FIL" = "sxxxxxxxxxxxxxxxxxxxx.fcs"),
                        dry_run = TRUE)
eq("4 bytes of slack cannot absorb 20 bytes of growth", r$strategy, "rewrite")
p <- fx("just_enough.fcs")
fcs_fixture(p, keywords = list("$FIL" = "s.fcs"), slack = 64)
r <- fcs_write_keywords(p, set = c("$FIL" = "sxxxxxxxxxxxxxxxxxxxx.fcs"),
                        dry_run = TRUE)
eq("64 bytes of slack absorbs the same growth", r$strategy, "patch")

# --- segment resolution -------------------------------------------------------
cat("\n-- segment resolution ----------------------------------------------\n")
# ANALYSIS between TEXT and DATA, hidden from the header
p <- fx("ana.fcs")
ana <- as.raw(rev(seq_len(64)))
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), analysis = ana,
            analysis_before_data = TRUE, header_zero = "analysis")
i0 <- fcs_inspect(p)
eq("ANALYSIS resolved from its keywords", i0$seg$anastart,
   fcs_text_num(i0$text, "$BEGINANALYSIS"))
eq("body_start is the ANALYSIS segment, not DATA",
   i0$seg$body_start, i0$seg$anastart)
ana0 <- fx_bytes(p, i0$seg$anastart, 64)
data0 <- fx_bytes(p, i0$seg$datastart, 40)
r <- fcs_write_keywords(p, set = c("$FIL" = "a_much_longer_replacement.fcs"))
eq("ANALYSIS file: write succeeds", r$status, "ok")
i1 <- fcs_inspect(p)
ok("ANALYSIS bytes are byte-identical", identical(fx_bytes(p, i1$seg$anastart, 64), ana0))
ok("DATA bytes are byte-identical", identical(fx_bytes(p, i1$seg$datastart, 40), data0))
eq("both ANALYSIS keywords shifted", i1$seg$anastart - i0$seg$anastart, r$grow)
eq("$ENDANALYSIS shifted too",
   fcs_text_num(i1$text, "$ENDANALYSIS") - fcs_text_num(i0$text, "$ENDANALYSIS"),
   r$grow)
eq("the header ANALYSIS field stays zeroed", i1$header$anastart, 0)

# ANALYSIS after DATA
p <- fx("ana_after.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), analysis = ana)
i0 <- fcs_inspect(p)
ana0 <- fx_bytes(p, i0$seg$anastart, 64)
r <- fcs_write_keywords(p, set = c("$FIL" = "a_much_longer_replacement.fcs"))
i1 <- fcs_inspect(p)
eq("ANALYSIS-after-DATA: write succeeds", r$status, "ok")
ok("ANALYSIS-after-DATA: bytes identical",
   identical(fx_bytes(p, i1$seg$anastart, 64), ana0))
eq("ANALYSIS-after-DATA: header field shifted",
   i1$header$anastart - i0$header$anastart, r$grow)

# header zeroed / sentinel, truth in TEXT
for (style in c("zero", "sentinel")) {
  p <- fx(paste0("hdr_", style, ".fcs"))
  fcs_fixture(p, keywords = list("$FIL" = "old.fcs"),
              header_zero = if (style == "zero") "data" else character(),
              header_sentinel = if (style == "sentinel") "data" else character())
  i0 <- fcs_inspect(p)
  eq(paste0(style, ": DATA resolved from $BEGINDATA"), i0$seg$datastart,
     fcs_text_num(i0$text, "$BEGINDATA"))
  data0 <- fx_bytes(p, i0$seg$datastart, 40)
  r <- fcs_write_keywords(p, set = c("$FIL" = "a_much_longer_replacement.fcs"))
  eq(paste0(style, ": write succeeds"), r$status, "ok")
  i1 <- fcs_inspect(p)
  eq(paste0(style, ": header convention preserved, not repaired"),
     i1$header$datastart, if (style == "zero") 0 else 99999999)
  eq(paste0(style, ": $BEGINDATA carries the shift"),
     i1$seg$datastart - i0$seg$datastart, r$grow)
  ok(paste0(style, ": data bytes identical"),
     identical(fx_bytes(p, i1$seg$datastart, 40), data0))
}

# --- refusals -----------------------------------------------------------------
cat("\n-- refusals --------------------------------------------------------\n")
p <- fx("conflict.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"))
i0 <- fcs_inspect(p)
# Patch the HEADER's datastart field directly, so it disagrees with the
# $BEGINDATA the builder wrote. Same width, so nothing else shifts.
fx_poke(p, 26, sprintf("%8.0f", i0$seg$datastart + 98))
refuses("header and $BEGINDATA disagree", p, "header_text_offset_conflict")

p <- fx("nextdata.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs", "$NEXTDATA" = "4096"))
refuses("$NEXTDATA is non-zero (a second data set)", p, "nextdata_nonzero")

p <- fx("stext.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs", "$BEGINSTEXT" = "9000",
                               "$ENDSTEXT" = "9100"))
refuses("supplemental TEXT is present", p, "supplemental_text")

p <- fx("textafter.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), text_after_data = TRUE)
refuses("TEXT sits after DATA", p, "text_after_data")

p <- fx("spacedelim.fcs")
fcs_fixture(p, delimiter = " ", keywords = list("$FIL" = "old.fcs"))
refuses("a space delimiter", p, "bad_delimiter")

p <- fx("digitdelim.fcs")
fcs_fixture(p, delimiter = "7", keywords = list("$FIL" = "old.fcs"))
refuses("a digit delimiter", p, "bad_delimiter")

p <- fx("dupkey.fcs")
fcs_fixture(p, keywords = list("$FIL" = "one.fcs", "$fil" = "two.fcs"))
refuses("the same keyword twice in different case", p, "duplicate_keyword")

p <- fx("unpaired.fcs")
fcs_fixture(p, raw_text = "|$TOT|5|$PAR|1|$FIL|")
refuses("an odd number of TEXT tokens", p, "malformed_text_unpaired")

p <- fx("unterminated.fcs")
fcs_fixture(p, raw_text = "|$TOT|5|$FIL|old.fcs")
refuses("TEXT that does not end on a delimiter", p, "malformed_text_unterminated")

p <- fx("emptykey.fcs")
fcs_fixture(p, raw_text = "||A|1|B|")
refuses("a zero-length keyword", p, "empty_keyword_or_value")

p <- fx("crc.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), trailer = "A3F19C22")
refuses("a real CRC we would invalidate", p, "crc_present")
r <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"), allow = "crc_present")
eq("--allow crc_present lets it through", r$status, "ok")

p <- fx("trailing.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"))
con <- file(p, "ab"); writeBin(charToRaw("FCS3.0    extra junk"), con); close(con)
refuses("bytes past the CRC field (an unannounced second data set)", p,
        "unexplained_trailing_bytes")

p <- fx("other.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), other_bytes = 16L)
fx_poke(p, 58, "     900    1000")
refuses("a real OTHER segment we would have to shift", p, "other_segments")

p <- fx("other_pad.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), other_bytes = 16L)
ins <- fcs_inspect(p)
eq("all-zero OTHER padding is tolerated", ins$refusals, character(0))
eq("...and TEXT is found where the header says", ins$seg$textstart, 74)
r <- fcs_write_keywords(p, set = c("$FIL" = "a_much_longer_replacement.fcs"))
eq("...and a rewrite preserves it", r$status, "ok")
eq("...with TEXT still starting at 74", fcs_inspect(p)$seg$textstart, 74)

# BD FACSDiva starts TEXT at byte 256 and space-pads 58..255. That gap is 198
# bytes — NOT a multiple of 8 — so a check that only understands 8-byte offset
# pairs refuses an entire Diva cohort on the strength of the remainder. A region
# with no digits in it cannot encode an offset, whatever its length.
p <- fx("diva_pad.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), other_bytes = 198L)
ins <- fcs_inspect(p)
eq("a 198-byte space-padded gap is padding, not an OTHER segment",
   ins$refusals, character(0))
eq("...and TEXT starts where Diva puts it", ins$seg$textstart, 256)
eq("...and the gap length is deliberately not a multiple of 8", 198L %% 8L, 6L)
r <- fcs_write_keywords(p, set = c("$FIL" = "a_much_longer_replacement.fcs"))
eq("...and a rewrite through it succeeds", r$status, "ok")
eq("...with TEXT still starting at 256", fcs_inspect(p)$seg$textstart, 256)
eq("...and the new value readable", fcs_keywords(p)[["$FIL"]],
   "a_much_longer_replacement.fcs")

# The refusal must still fire on a gap that carries real offsets, whatever its
# length — otherwise the fix above would have disarmed the check entirely.
p <- fx("other_odd.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"), other_bytes = 22L)
fx_poke(p, 58, "     900    1000")
refuses("a real OTHER segment in a non-multiple-of-8 gap", p, "other_segments")

p <- fx("v2.fcs")
fcs_fixture(p, version = "FCS2.0", keywords = list("$FIL" = "old.fcs"))
refuses("FCS2.0, where $BEGINDATA is not defined", p, "unsupported_version")

p <- fx("symlink.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"))
link <- fx("link.fcs")
if (file.symlink(p, link)) {
  refuses("a symlink (we do not follow it)", link, "symlink")
} else {
  cat("  SKIP  symlinks unavailable here\n")
}

# --- tolerated quirks ---------------------------------------------------------
cat("\n-- tolerated quirks ------------------------------------------------\n")
p <- fx("offbyone.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"))
i0 <- fcs_inspect(p)
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"),
            force_textend = i0$seg$datastart)
ins <- fcs_inspect(p)
ok("textend written as exclusive is absorbed",
   "textend_off_by_one" %in% ins$warnings)
eq("...and the real TEXT span is used", ins$seg$textend, i0$seg$textend)
eq("...and it still parses", fcs_text_get(ins$text, "$FIL"), "old.fcs")

p <- fx("modec.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs", "$MODE" = "C"))
ins <- fcs_inspect(p)
ok("$MODE C skips the length invariant rather than erroring",
   is.na(ins$data_length_ok) && length(ins$refusals) == 0L)
eq("...and still writes", fcs_write_keywords(p, set = c("$FIL" = "new.fcs"))$status,
   "ok")

p <- fx("varwidth.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs", "$P1B" = "*"))
ins <- fcs_inspect(p)
ok("$PnB = * skips the length invariant", is.na(ins$data_length_ok))
eq("...and still writes", fcs_write_keywords(p, set = c("$FIL" = "new.fcs"))$status,
   "ok")

# --- what we refuse to write --------------------------------------------------
cat("\n-- what we refuse to write -----------------------------------------\n")
p <- fx("guard.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"))
before <- fx_md5(p)
r <- fcs_write_keywords(p, set = c("$FIL" = "has|a|pipe.fcs"))
eq("a new value containing the delimiter", r$refusals, "new_value_contains_delimiter")
ok("...tells the user to rename the file instead", grepl("Rename the file", r$message))
r <- fcs_write_keywords(p, set = c("$FIL" = "has\nnewline.fcs"))
eq("a new value containing a control byte", r$refusals, "new_value_unprintable")
for (k in c("$TOT", "$PAR", "$P1B", "$DATE", "ORIGINALGUID", "$BEGINDATA")) {
  r <- fcs_write_keywords(p, set = structure(list("x"), names = k))
  eq(paste0("protected keyword ", k), r$refusals, "protected_keyword")
}
ok("none of that touched the file", identical(fx_md5(p), before))

# --- dry run, idempotency, copies ---------------------------------------------
cat("\n-- dry run, idempotency, copies ------------------------------------\n")
p <- fx("dry.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"))
before <- fx_md5(p); mt <- file.mtime(p)
r <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"), dry_run = TRUE)
eq("dry run reports ok", r$status, "ok")
ok("dry run says so", grepl("DRY RUN", r$message))
ok("dry run leaves the bytes alone", identical(fx_md5(p), before))
ok("dry run leaves the mtime alone", identical(file.mtime(p), mt))

r1 <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"))
eq("first apply writes", r1$status, "ok")
after1 <- fx_md5(p)
r2 <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"))
eq("second apply is a no-op", r2$status, "unchanged")
ok("...and wrote nothing", identical(fx_md5(p), after1))

p <- fx("copysrc.fcs")
fcs_fixture(p, keywords = list("$FIL" = "old.fcs"))
before <- fx_md5(p)
dst <- fx("copydst.fcs")
r <- fcs_write_keywords(p, set = c("$FIL" = "copied.fcs"), output_path = dst)
eq("copy: status", r$status, "ok")
eq("copy: reports the real destination, not the temp file", r$output_path, dst)
ok("copy: source untouched", identical(fx_md5(p), before))
eq("copy: destination carries the new value", fcs_keyword(dst, "$FIL"), "copied.fcs")
eq("copy: destination is clean", fcs_inspect(dst)$refusals, character(0))
ok("copy: no temp files left behind",
   length(list.files(root, pattern = "fcstmp|fcsold|fcsprefix")) == 0L)

# --- failure paths ------------------------------------------------------------
cat("\n-- failure paths ---------------------------------------------------\n")
p <- fx("fault_rw.fcs")
fcs_fixture(p, keywords = list())            # no $FIL -> rewrite
before <- fx_md5(p)
r <- fcs_write_keywords(p, set = c("$FIL" = "new.fcs"), .fault = "corrupt_tail")
eq("rewrite with a corrupted tail fails", r$status, "failed")
ok("...naming the data segment", grepl("data segment", r$message))
ok("...leaving the original byte-identical", identical(fx_md5(p), before))
ok("...and no temp file behind",
   length(list.files(root, pattern = "fcstmp|fcsold")) == 0L)

p <- fx("fault_patch.fcs")
fcs_fixture(p, keywords = list("$FIL" = "a_long_original_name.fcs"))
before <- fx_md5(p)
r <- fcs_write_keywords(p, set = c("$FIL" = "s.fcs"), .fault = "corrupt_prefix")
eq("in-place patch with a corrupted write fails", r$status, "failed")
ok("...says the original was restored", grepl("restored", r$message))
ok("...and the original really is back", identical(fx_md5(p), before))
ok("...with the recovery sidecar cleaned up",
   length(list.files(root, pattern = "fcsprefix")) == 0L)

# --- large-offset arithmetic (no huge file needed) ----------------------------
cat("\n-- large offsets ---------------------------------------------------\n")
big <- 4.5e9
eq("offsets past 2^31 stay doubles", typeof(big + 1), "double")
eq("sprintf renders them without scientific notation", sprintf("%.0f", big),
   "4500000000")
ok("an 8-byte header field refuses to hold them", is.null(.fcs_field8(big)))
ok("...and holds a value that does fit", !is.null(.fcs_field8(99999999)))
ok("...rejecting one digit more", is.null(.fcs_field8(100000000)))
pl <- fcs_plan_layout(big, big, big + 100)
eq("the planner works at that scale", pl$strategy, "rewrite")
eq("...and reports the growth as a double", pl$grow, 100)
eq("as.integer would have destroyed it", suppressWarnings(as.integer(big)),
   NA_integer_)
src <- readLines(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
ok("the layer never calls strtoi (32-bit only)", !any(grepl("strtoi\\(", src)))
ok("the layer never calls readChar (encoding round-trip)",
   !any(grepl("readChar\\(", src)))

cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")
if (.fail > 0) quit(status = 1)
