# src/harmonization/drift/drift_metadata.R — tube table from file names and FCS keywords
# Base R only, on top of R_shared/fcs_keywords.R, so it runs anywhere and is
# testable with the byte-level fixtures. Produces one row per FCS file with the
# biological fields from the file name and the technical fields from TEXT.
#
# Which fields a file name carries, and how instruments are labelled, come from
# the cohort config (`tubes` and `instrument_labels` blocks), so no cohort
# naming lives here.

#' Parse one file name with the named-group regex of the config `tubes` block.
#' Unmatched names give a row of NAs. `is_neg` is TRUE when the `cocktail`
#' field matches `neg_regex`; such tubes get cocktail "neg".
parse_filename <- function(filename, spec) {
  stem <- sub("\\.fcs$", "", basename(filename), ignore.case = TRUE)
  m <- regexpr(spec$filename_regex, stem, perl = TRUE)
  fields <- attr(m, "capture.names")
  if (is.null(fields) || any(!nzchar(fields)))
    stop("tubes$filename_regex must name every capture group, e.g. (?<patient>\\d+)", call. = FALSE)
  vals <- if (m > 0) {
    st <- attr(m, "capture.start")[1, ]; ln <- attr(m, "capture.length")[1, ]
    substring(stem, st, st + ln - 1L)
  } else rep(NA_character_, length(fields))
  out <- as.list(setNames(vals, fields))
  for (f in intersect(fields, spec$integer_fields)) out[[f]] <- as.integer(out[[f]])
  out <- c(list(filename = basename(filename)), out)
  df <- as.data.frame(out, stringsAsFactors = FALSE)
  flag_neg(df, spec$neg_regex)
}

#' Add `is_neg` from the cocktail column (unless the table already has one).
flag_neg <- function(df, neg_regex) {
  if (is.null(df$is_neg)) {
    if (is.null(df$cocktail)) stop("tube metadata needs a `cocktail` field or an `is_neg` column", call. = FALSE)
    df$is_neg <- grepl(neg_regex, df$cocktail, ignore.case = TRUE)
  }
  df$is_neg <- as.logical(df$is_neg)
  df$is_neg[is.na(df$cocktail) & is.na(df$is_neg)] <- NA
  df$cocktail[df$is_neg %in% TRUE] <- "neg"
  df
}

#' Instrument label from $CYT: the first `instrument_labels` key found in it
#' (case-insensitive) wins; otherwise the cleaned $CYT string.
instrument_label <- function(cyt, labels) {
  cyt <- ifelse(is.na(cyt), "", cyt)
  out <- ifelse(nzchar(cyt), gsub("[^A-Za-z0-9]", "", cyt), NA_character_)
  for (key in rev(names(labels))) out <- ifelse(grepl(key, cyt, ignore.case = TRUE), labels[[key]], out)
  out
}

#' $DATE is dd-MMM-yyyy in Diva files; parse locale-independently.
parse_fcs_date <- function(x) {
  mon <- c(JAN = 1, FEB = 2, MAR = 3, APR = 4, MAY = 5, JUN = 6, JUL = 7,
           AUG = 8, SEP = 9, OCT = 10, NOV = 11, DEC = 12)
  m <- regmatches(x, regexec("(\\d{1,2})-([A-Za-z]{3})-(\\d{4})", x))
  out <- rep(as.Date(NA), length(x))
  for (i in seq_along(m)) {
    if (length(m[[i]]) == 4L) {
      mm <- mon[toupper(m[[i]][3])]
      if (!is.na(mm)) out[i] <- as.Date(sprintf("%s-%02d-%02d", m[[i]][4], mm,
                                                 as.integer(m[[i]][2])))
    }
  }
  out
}

#' Technical keywords of one file. Missing keywords become NA, never an error.
fcs_technical_keywords <- function(path, instrument_labels = list()) {
  kw <- tryCatch(fcs_keywords(path), error = function(e) character(0))
  get <- function(k) {
    hit <- which(toupper(names(kw)) == toupper(k))
    if (length(hit)) kw[[hit[1]]] else NA_character_
  }
  data.frame(
    acq_date     = as.character(parse_fcs_date(get("$DATE"))),
    instrument   = instrument_label(get("$CYT"), instrument_labels),
    cytnum       = get("CYTNUM"),
    software     = sub("^.*Version\\s*", "", get("CREATOR")),
    n_par        = suppressWarnings(as.integer(get("$PAR"))),
    n_events     = suppressWarnings(as.numeric(get("$TOT"))),
    operator     = get("$OP"),
    experiment   = get("EXPERIMENT NAME"),
    tube_name    = get("TUBE NAME"),
    btim         = get("$BTIM"),
    has_spill    = !is.na(get("SPILL")) || !is.na(get("$SPILLOVER")),
    stringsAsFactors = FALSE)
}

#' The tube table: one row per .fcs under fcs_dir (recursive). Biological fields come
#' from `sample_sheet` (data frame with a `filename` column) when given, otherwise
#' from the file name via `cfg$tubes`; technical fields from the FCS TEXT segment.
build_tube_table <- function(fcs_dir, cfg, sample_sheet = NULL) {
  files <- list.files(fcs_dir, pattern = "\\.fcs$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("no .fcs files under ", fcs_dir, call. = FALSE)
  bio <- if (is.null(sample_sheet)) {
    do.call(rbind, lapply(files, parse_filename, spec = cfg$tubes))
  } else {
    if (is.null(sample_sheet$filename)) stop("sample sheet needs a `filename` column", call. = FALSE)
    sample_sheet$filename <- basename(sample_sheet$filename)
    if (anyDuplicated(sample_sheet$filename)) stop("sample sheet has duplicated file names", call. = FALSE)
    ss <- sample_sheet[match(basename(files), sample_sheet$filename), , drop = FALSE]
    ss$filename <- basename(files)
    rownames(ss) <- NULL
    flag_neg(ss, cfg$tubes$neg_regex)
  }
  need <- c("patient", "cocktail", "is_neg")
  if (length(miss <- setdiff(need, names(bio))))
    stop("tube metadata lacks required field(s): ", paste(miss, collapse = ", "), call. = FALSE)
  tech <- do.call(rbind, lapply(files, fcs_technical_keywords, instrument_labels = cfg$instrument_labels))
  tt <- cbind(bio, path = files, tech, stringsAsFactors = FALSE)
  tt$acq_day <- tt$acq_date                       # the technical "day" factor
  tt$patient_day <- paste(tt$patient, tt$acq_day, sep = "_")
  ord <- if (!is.null(tt$week)) order(tt$patient, tt$week, tt$cocktail) else order(tt$patient, tt$cocktail)
  tt[ord, ]
}

#' Which patients make day drift and instrument shift identifiable.
bridge_report <- function(tt) {
  per_pat <- split(tt, tt$patient)
  data.frame(
    patient       = names(per_pat),
    n_tubes       = vapply(per_pat, nrow, integer(1)),
    n_negs        = vapply(per_pat, function(d) sum(d$is_neg %in% TRUE), integer(1)),
    n_days        = vapply(per_pat, function(d) length(unique(d$acq_day)), integer(1)),
    instruments   = vapply(per_pat, function(d)
      paste(sort(unique(d$instrument)), collapse = "+"), character(1)),
    day_bridge    = vapply(per_pat, function(d) length(unique(d$acq_day)) > 1L, logical(1)),
    instr_bridge  = vapply(per_pat, function(d) length(unique(d$instrument)) > 1L, logical(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}
