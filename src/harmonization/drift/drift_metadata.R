# =============================================================================
# src/harmonization/drift/drift_metadata.R — tube table from file names and FCS keywords
# =============================================================================
# Base R only, on top of R_shared/fcs_keywords.R, so it runs anywhere and is
# testable with the byte-level fixtures. Produces one row per FCS file with the
# biological fields from the file name and the technical fields from TEXT.
#
# The file-name field called "batch" is PRIOR VACCINATION STATUS (2 = previously
# vaccinated). It is kept under the honest name `prior_vacc` and is never a
# batch variable.
# =============================================================================

MTB_NAME_RE <- "^MTB(\\d+)_(\\d+)_V(\\d+)W(\\d+)_(.+)$"

#' Parse MTB<patient>_<prior_vacc>_V<vaccine>W<week>_<cocktail>.fcs
parse_mtb_filename <- function(filename) {
  stem <- sub("\\.fcs$", "", basename(filename), ignore.case = TRUE)
  m <- regmatches(stem, regexec(MTB_NAME_RE, stem))[[1]]
  if (length(m) == 0L) {
    return(data.frame(filename = basename(filename), patient = NA_integer_,
                      prior_vacc = NA_integer_, vaccine = NA_integer_,
                      week = NA_integer_, cocktail = NA_character_,
                      is_neg = NA, stringsAsFactors = FALSE))
  }
  cocktail <- m[6]
  is_neg <- grepl("neg$", cocktail, ignore.case = TRUE)
  data.frame(filename = basename(filename), patient = as.integer(m[2]),
             prior_vacc = as.integer(m[3]), vaccine = as.integer(m[4]),
             week = as.integer(m[5]), cocktail = if (is_neg) "neg" else cocktail,
             is_neg = is_neg, stringsAsFactors = FALSE)
}

#' Instrument label from $CYT, collapsing Diva's model strings.
instrument_label <- function(cyt) {
  cyt <- ifelse(is.na(cyt), "", cyt)
  ifelse(grepl("Symphony", cyt, ignore.case = TRUE), "FACSymphonyA1",
  ifelse(grepl("Canto", cyt, ignore.case = TRUE), "FACSCantoII",
  ifelse(nzchar(cyt), gsub("[^A-Za-z0-9]", "", cyt), NA_character_)))
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
fcs_technical_keywords <- function(path) {
  kw <- tryCatch(fcs_keywords(path), error = function(e) character(0))
  get <- function(k) {
    hit <- which(toupper(names(kw)) == toupper(k))
    if (length(hit)) kw[[hit[1]]] else NA_character_
  }
  data.frame(
    acq_date     = as.character(parse_fcs_date(get("$DATE"))),
    instrument   = instrument_label(get("$CYT")),
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

#' The tube table: one row per .fcs under fcs_dir (recursive).
build_tube_table <- function(fcs_dir) {
  files <- list.files(fcs_dir, pattern = "\\.fcs$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("no .fcs files under ", fcs_dir, call. = FALSE)
  rows <- lapply(files, function(f) {
    cbind(parse_mtb_filename(f), path = f, fcs_technical_keywords(f),
          stringsAsFactors = FALSE)
  })
  tt <- do.call(rbind, rows)
  tt$acq_day <- tt$acq_date                       # the technical "day" factor
  tt$patient_day <- paste(tt$patient, tt$acq_day, sep = "_")
  tt[order(tt$patient, tt$week, tt$cocktail), ]
}

#' Which patients make day drift and instrument shift identifiable.
bridge_report <- function(tt) {
  per_pat <- split(tt, tt$patient)
  data.frame(
    patient       = as.integer(names(per_pat)),
    n_tubes       = vapply(per_pat, nrow, integer(1)),
    n_negs        = vapply(per_pat, function(d) sum(d$is_neg %in% TRUE), integer(1)),
    n_days        = vapply(per_pat, function(d) length(unique(d$acq_day)), integer(1)),
    instruments   = vapply(per_pat, function(d)
      paste(sort(unique(d$instrument)), collapse = "+"), character(1)),
    day_bridge    = vapply(per_pat, function(d) length(unique(d$acq_day)) > 1L, logical(1)),
    instr_bridge  = vapply(per_pat, function(d) length(unique(d$instrument)) > 1L, logical(1)),
    row.names = NULL, stringsAsFactors = FALSE)
}
