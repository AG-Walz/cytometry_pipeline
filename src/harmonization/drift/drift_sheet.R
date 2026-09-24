# src/harmonization/drift/drift_sheet.R — the sample sheet: the single metadata input
# Base R only. The sheet (docs/plans/drift_usability_plan.md, section 2) is pre-filled
# by the pipeline and completed by the user. Its user-facing column names are mapped
# here onto the internal names every later step uses (patient, week, cocktail, is_neg,
# one column per protected covariate), so nothing downstream knows about the sheet.

#' Fixed sheet columns, in file order. Any number of `group_*` columns (protected
#' biological covariates) may follow `control_for`.
SHEET_COLUMNS <- c("filename", "acq_date", "instrument", "n_events", "software",
                   "subject", "visit", "weeks_since_first_vaccination", "stimulation",
                   "control_for", "exclude", "note")
SHEET_HEADER_COLUMNS <- c("acq_date", "instrument", "n_events", "software")
SHEET_UNSTIM <- "unstimulated"

#' Supported `control_for` keys and the tubes.csv columns an unstimulated tube is
#' matched on (the `controls.match_keys` of the cohort config).
CONTROL_FOR_KEYS <- list("subject+visit" = c("patient", "week"),
                         "subject+acq_date" = c("patient", "acq_day"),
                         "subject" = "patient")

sheet_group_names <- function(sheet) sub("^group_", "", grep("^group_", names(sheet), value = TRUE))

#' Pre-filled sheet for every .fcs under `fcs_dir`: technical columns from the FCS
#' header, biological columns guessed from the file names when `cfg$tubes$filename_regex`
#' exists (patient -> subject, week -> weeks_since_first_vaccination, cocktail ->
#' stimulation with the neg rule, prior_vacc -> group_prior_peptide_vaccination,
#' vaccine + week -> visit "V3W6"), blank otherwise. Ordered by subject and visit.
init_sample_sheet <- function(fcs_dir, cfg = list()) {
  files <- list.files(fcs_dir, pattern = "\\.fcs$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("no .fcs files under ", fcs_dir, call. = FALSE)
  tech <- do.call(rbind, lapply(files, fcs_technical_keywords, instrument_labels = cfg$instrument_labels))
  n <- length(files)
  blank <- rep(NA_character_, n)
  sheet <- data.frame(filename = basename(files), tech[, SHEET_HEADER_COLUMNS],
                      subject = blank, visit = blank, weeks_since_first_vaccination = rep(NA_real_, n),
                      stimulation = blank, control_for = rep(names(CONTROL_FOR_KEYS)[1], n),
                      stringsAsFactors = FALSE)
  if (!is.null(cfg$tubes$filename_regex)) {
    bio <- do.call(rbind, lapply(files, parse_filename, spec = cfg$tubes))
    sheet$subject <- bio$patient
    if (!is.null(bio$week)) {
      sheet$weeks_since_first_vaccination <- bio$week
      sheet$visit <- if (!is.null(bio$vaccine)) paste0("V", bio$vaccine, "W", bio$week) else paste0("W", bio$week)
      sheet$visit[is.na(bio$week)] <- NA_character_
    }
    sheet$stimulation <- ifelse(bio$is_neg %in% TRUE, SHEET_UNSTIM, bio$cocktail)
    if (!is.null(bio$prior_vacc)) sheet$group_prior_peptide_vaccination <- bio$prior_vacc
  }
  sheet$exclude <- FALSE
  sheet$note <- blank
  sheet <- sheet[order(sheet$subject, sheet$weeks_since_first_vaccination, sheet$visit, sheet$stimulation), ]
  rownames(sheet) <- NULL
  sheet
}

#' Read and validate a sample sheet; return the internal table used by
#' build_tube_table(): filename, patient, visit, week, cocktail, is_neg, control_for,
#' one column per group_* (prefix stripped), exclude, note. The header-derived columns
#' are kept in attr(, "header") for check_sheet_against_files() only.
read_sample_sheet <- function(path) {
  if (!file.exists(path)) stop("sample sheet not found: ", path, call. = FALSE)
  raw <- read.csv(path, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE,
                  na.strings = character(0))
  for (j in seq_along(raw)) raw[[j]] <- trimws(raw[[j]])
  need <- setdiff(SHEET_COLUMNS, SHEET_HEADER_COLUMNS)
  if (length(miss <- setdiff(need, names(raw))))
    stop("sample sheet lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  raw$filename <- basename(raw$filename)
  if (any(!nzchar(raw$filename))) stop("sample sheet has rows without a filename", call. = FALSE)
  if (anyDuplicated(raw$filename))
    stop("sample sheet has duplicated file names: ",
         paste(unique(raw$filename[duplicated(raw$filename)]), collapse = ", "), call. = FALSE)

  ex <- toupper(raw$exclude)
  yes <- c("TRUE", "T", "1", "YES"); no <- c("FALSE", "F", "0", "NO", "")
  if (any(bad <- !ex %in% c(yes, no)))
    stop("exclude must be TRUE or FALSE; offending rows: ", paste(raw$filename[bad], collapse = ", "), call. = FALSE)
  exclude <- ex %in% yes

  active <- !exclude
  for (col in c("subject", "stimulation")) {
    if (any(bad <- active & !nzchar(raw[[col]])))
      stop("`", col, "` is empty for: ", paste(raw$filename[bad], collapse = ", "), call. = FALSE)
  }
  wk <- raw$weeks_since_first_vaccination
  week <- suppressWarnings(as.numeric(wk))
  if (any(bad <- nzchar(wk) & is.na(week)))
    stop("weeks_since_first_vaccination must be numeric; offending rows: ",
         paste(raw$filename[bad], collapse = ", "), call. = FALSE)
  if (all(is.na(week) | week == round(week))) week <- as.integer(week)

  cf <- raw$control_for
  cf[!nzchar(cf)] <- names(CONTROL_FOR_KEYS)[1]
  if (any(bad <- !cf %in% names(CONTROL_FOR_KEYS)))
    stop("control_for must be one of ", paste(names(CONTROL_FOR_KEYS), collapse = ", "),
         "; offending rows: ", paste(raw$filename[bad], collapse = ", "), call. = FALSE)

  as_value <- function(x) { x[!nzchar(x)] <- NA; type.convert(x, as.is = TRUE) }
  is_neg <- tolower(raw$stimulation) == SHEET_UNSTIM
  out <- data.frame(filename = raw$filename, patient = as_value(raw$subject), visit = as_value(raw$visit),
                    week = week, cocktail = ifelse(is_neg, "neg", raw$stimulation), is_neg = is_neg,
                    control_for = cf, stringsAsFactors = FALSE)
  out$cocktail[!nzchar(raw$stimulation)] <- NA_character_
  groups <- grep("^group_", names(raw), value = TRUE)
  for (g in groups) out[[sub("^group_", "", g)]] <- as_value(raw[[g]])
  out$exclude <- exclude
  out$note <- as_value(raw$note)
  attr(out, "header") <- raw[, intersect(c("filename", SHEET_HEADER_COLUMNS), names(raw)), drop = FALSE]
  attr(out, "groups") <- sub("^group_", "", groups)
  out
}

#' Warn when acq_date or instrument in the sheet differ from the FCS header (`tt`,
#' the tube table). Returns the mismatching rows invisibly.
check_sheet_against_files <- function(sheet, tt) {
  hdr <- attr(sheet, "header")
  cols <- intersect(c("acq_date", "instrument"), names(hdr))
  if (is.null(hdr) || !length(cols)) return(invisible(NULL))
  i <- match(hdr$filename, tt$filename)
  bad <- rep(FALSE, nrow(hdr))
  for (col in cols) {
    s <- hdr[[col]]; h <- as.character(tt[[col]][i])
    bad <- bad | (!is.na(i) & nzchar(s) & !is.na(h) & s != h) | (!is.na(i) & nzchar(s) & is.na(h))
  }
  mm <- hdr[bad, , drop = FALSE]
  if (nrow(mm)) warning("sheet and FCS header disagree on acq_date/instrument for: ",
                        paste(mm$filename, collapse = ", "), call. = FALSE)
  invisible(mm)
}

#' The `design` block implied by a sheet (either form): the group_* columns are the
#' protected fixed covariates, everything else is the standard layout.
sheet_design <- function(sheet) {
  groups <- attr(sheet, "groups")
  if (is.null(groups)) groups <- sheet_group_names(sheet)
  list(unit = "patient", technical = c("acq_day", "instrument"), protected_fixed = groups,
       time = "week", stimulation = "cocktail")
}

#' `controls.match_keys` for a `control_for` key (one value, or the sheet's values,
#' which must agree).
control_match_keys <- function(control_for) {
  control_for <- unique(control_for[!is.na(control_for) & nzchar(control_for)])
  if (!length(control_for)) control_for <- names(CONTROL_FOR_KEYS)[1]
  if (length(control_for) > 1L)
    stop("control_for must be the same for every row of the sheet; found: ",
         paste(control_for, collapse = ", "), call. = FALSE)
  if (!control_for %in% names(CONTROL_FOR_KEYS))
    stop("unsupported control_for key: ", control_for, call. = FALSE)
  CONTROL_FOR_KEYS[[control_for]]
}

#' Fill the `{{cohort}}`, `{{protected_fixed}}` and `{{match_keys}}` placeholders of
#' `template_path` and write the cohort config to `out_path`.
write_cohort_config <- function(template_path, out_path, cohort, design, controls_match_keys) {
  txt <- readLines(template_path, warn = FALSE)
  flow <- function(x) paste0("[", paste(x, collapse = ", "), "]")
  txt <- gsub("{{cohort}}", cohort, txt, fixed = TRUE)
  txt <- gsub("{{protected_fixed}}", flow(unlist(design$protected_fixed)), txt, fixed = TRUE)
  txt <- gsub("{{match_keys}}", flow(controls_match_keys), txt, fixed = TRUE)
  left <- grep("\\{\\{[a-z_]+\\}\\}", txt, value = TRUE)
  if (length(left)) stop("unfilled placeholder(s) in template: ", paste(left, collapse = "; "), call. = FALSE)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(txt, out_path)
  invisible(out_path)
}
