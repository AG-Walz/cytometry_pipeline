#!/usr/bin/env Rscript
# test_data/test_drift_sheet.R — unit tests for src/harmonization/drift/drift_sheet.R
# Base R plus test_data/fcs_fixtures.R.
#   Rscript test_data/test_drift_sheet.R
cat("=== drift sample sheet tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_sheet.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
for (m in list.files(file.path(base_dir, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)
source(file.path(base_dir, "test_data", "fcs_fixtures.R"))

n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") }
}
fails <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  !is.null(e) && grepl(pattern, e)
}
tmp <- tempfile("sheet_"); dir.create(tmp)
fcs <- file.path(tmp, "fcs"); dir.create(fcs)

cfg <- list(
  tubes = list(filename_regex = "^MTB(?<patient>\\d+)_(?<prior_vacc>\\d+)_V(?<vaccine>\\d+)W(?<week>\\d+)_(?<cocktail>.+)$",
               integer_fields = c("patient", "prior_vacc", "vaccine", "week"), neg_regex = "neg$"),
  instrument_labels = list(Symphony = "FACSymphonyA1", Canto = "FACSCantoII"))
mk <- function(name, date, cyt, tot = 5L)
  fcs_fixture(file.path(fcs, name), par = 2L, tot = tot, keywords = list("$DATE" = date, "$CYT" = cyt,
              "CREATOR" = "BD FACSDiva Software Version 6.1.2", "$P1N" = "FSC-A", "$P2N" = "PE-A"))
mk("MTB012_2_V4W70_215008.fcs", "04-JUN-2025", "FACSCantoII", 7L)
mk("MTB012_2_V4W70_1_neg.fcs", "04-JUN-2025", "FACSCantoII")
mk("MTB011_1_V3W6_215008.fcs", "07-MAR-2024", "FACSymphonyA1")
mk("MTB011_1_V3W6_neg.fcs", "07-MAR-2024", "FACSymphonyA1")
mk("odd_name.fcs", "07-MAR-2024", "FACSymphonyA1")

# --- init with the MTB regex ---------------------------------------------------------
cat("init_sample_sheet with regex\n")
s <- init_sample_sheet(fcs, cfg)
check("fixed columns plus one group column", identical(setdiff(names(s), "group_prior_peptide_vaccination"), SHEET_COLUMNS) &&
        "group_prior_peptide_vaccination" %in% names(s))
check("five rows ordered by subject then weeks", nrow(s) == 5 && identical(s$subject[1:4], c(11L, 11L, 12L, 12L)) &&
        is.na(s$subject[5]))
r <- s[s$filename == "MTB012_2_V4W70_215008.fcs", ]
check("technical columns from the header", r$acq_date == "2025-06-04" && r$instrument == "FACSCantoII" && r$n_events == 7 && r$software == "6.1.2")
check("guesses: subject, weeks, visit, stimulation, group", r$subject == 12 && r$weeks_since_first_vaccination == 70 &&
        r$visit == "V4W70" && r$stimulation == "215008" && r$group_prior_peptide_vaccination == 2)
check("neg suffix -> unstimulated", s$stimulation[s$filename == "MTB012_2_V4W70_1_neg.fcs"] == SHEET_UNSTIM)
check("default control_for, exclude FALSE, note blank", all(s$control_for == "subject+visit") && !any(s$exclude) && all(is.na(s$note)))
odd <- s[s$filename == "odd_name.fcs", ]
check("unparsable name -> blank biology, header still filled", is.na(odd$subject) && is.na(odd$stimulation) && is.na(odd$visit) &&
        odd$instrument == "FACSymphonyA1")

cat("init_sample_sheet without regex\n")
s0 <- init_sample_sheet(fcs, list())
check("no group column, biology blank", identical(names(s0), SHEET_COLUMNS) && all(is.na(s0$subject)) && all(is.na(s0$stimulation)) &&
        all(is.na(s0$weeks_since_first_vaccination)))
check("instrument falls back to cleaned $CYT", all(s0$instrument %in% c("FACSCantoII", "FACSymphonyA1")))

# --- read_sample_sheet -----------------------------------------------------------------
cat("read_sample_sheet\n")
write_sheet <- function(df, name = "sheet.csv") { p <- file.path(tmp, name); write.csv(df, p, row.names = FALSE, na = ""); p }
s2 <- s
s2$exclude[s2$filename == "odd_name.fcs"] <- TRUE
s2$note[s2$filename == "odd_name.fcs"] <- "not part of the study"
p <- write_sheet(s2)
n <- read_sample_sheet(p)
check("internal columns", identical(names(n), c("filename", "patient", "visit", "week", "cocktail", "is_neg", "control_for",
                                                "prior_peptide_vaccination", "exclude", "note")))
check("types: integer patient and week, logical is_neg/exclude", is.integer(n$patient) && is.integer(n$week) && is.logical(n$is_neg) && is.logical(n$exclude))
check("unstimulated -> cocktail neg, is_neg", all(n$cocktail[n$is_neg] == "neg") && sum(n$is_neg) == 2)
check("group prefix stripped, values kept", identical(n$prior_peptide_vaccination[n$patient %in% 12], c(2L, 2L)))
check("excluded row may have blank subject", n$exclude[n$filename == "odd_name.fcs"] && is.na(n$patient[n$filename == "odd_name.fcs"]) &&
        n$note[n$filename == "odd_name.fcs"] == "not part of the study")
check("header kept as attribute, groups recorded", identical(names(attr(n, "header")), c("filename", SHEET_HEADER_COLUMNS)) &&
        identical(attr(n, "groups"), "prior_peptide_vaccination"))
s3 <- s2[s2$filename != "odd_name.fcs", ]; s3$control_for <- ""; s3$exclude <- ""
n3 <- read_sample_sheet(write_sheet(s3))
check("blank control_for and exclude take the defaults", all(n3$control_for == "subject+visit") && !any(n3$exclude))

cat("read_sample_sheet validation\n")
check("missing column", fails(read_sample_sheet(write_sheet(s2[, names(s2) != "stimulation"])), "lacks column.*stimulation"))
check("duplicated filename", fails(read_sample_sheet(write_sheet(rbind(s2, s2[1, ]))), "duplicated file names: MTB011_1_V3W6_215008.fcs"))
bad <- s2; bad$subject[1] <- NA
check("blank subject on active row", fails(read_sample_sheet(write_sheet(bad)), "`subject` is empty for: MTB011_1_V3W6_215008.fcs"))
bad <- s2; bad$stimulation[2] <- NA
check("blank stimulation on active row", fails(read_sample_sheet(write_sheet(bad)), "`stimulation` is empty"))
bad <- s2; bad$weeks_since_first_vaccination <- as.character(bad$weeks_since_first_vaccination); bad$weeks_since_first_vaccination[1] <- "six"
check("non-numeric weeks", fails(read_sample_sheet(write_sheet(bad)), "weeks_since_first_vaccination must be numeric"))
bad <- s2; bad$control_for[1] <- "visit"
check("unknown control_for", fails(read_sample_sheet(write_sheet(bad)), "control_for must be one of"))
bad <- s2; bad$exclude <- as.character(bad$exclude); bad$exclude[1] <- "maybe"
check("non-logical exclude", fails(read_sample_sheet(write_sheet(bad)), "exclude must be TRUE or FALSE"))
check("missing file", fails(read_sample_sheet(file.path(tmp, "nope.csv")), "not found"))

# --- header check, design, match keys, config ------------------------------------------
cat("check_sheet_against_files / sheet_design / write_cohort_config\n")
tt <- build_tube_table(fcs, cfg, n)
check("no warning when header agrees", is.null(withCallingHandlers(check_sheet_against_files(n, tt), warning = function(w) stop(w))) ||
        nrow(check_sheet_against_files(n, tt)) == 0)
s4 <- s2; s4$acq_date[1] <- "1999-01-01"; s4$instrument[2] <- "Other"
n4 <- read_sample_sheet(write_sheet(s4))
w <- tryCatch(check_sheet_against_files(n4, tt), warning = function(w) conditionMessage(w))
check("warns on acq_date/instrument mismatch, names both files", grepl("disagree", w) && grepl(s4$filename[1], w) && grepl(s4$filename[2], w))
d <- sheet_design(n)
check("design from the normalised sheet", identical(d$protected_fixed, "prior_peptide_vaccination") && d$time == "week" &&
        d$unit == "patient" && identical(d$technical, c("acq_day", "instrument")) && d$stimulation == "cocktail")
check("design from the init sheet", identical(sheet_design(s)$protected_fixed, "prior_peptide_vaccination") &&
        identical(sheet_design(s0)$protected_fixed, character(0)))
check("control_match_keys", identical(control_match_keys(s$control_for), c("patient", "week")) &&
        identical(control_match_keys("subject+acq_date"), c("patient", "acq_day")) &&
        identical(control_match_keys("subject"), "patient") &&
        fails(control_match_keys(c("subject", "subject+visit")), "same for every row"))
tpl <- file.path(base_dir, "config", "drift_template.yaml")
out_cfg <- file.path(tmp, "drift_test.yaml")
write_cohort_config(tpl, out_cfg, "test cohort", d, control_match_keys(n$control_for))
txt <- readLines(out_cfg)
check("placeholders filled", any(txt == 'cohort: "test cohort"') && any(grepl("^  protected_fixed: \\[prior_peptide_vaccination\\]", txt)) &&
        any(grepl("^  match_keys: \\[patient, week\\]", txt)) && !any(grepl("\\{\\{", txt)))
write_cohort_config(tpl, out_cfg, "c", sheet_design(s0), "patient")
check("empty protected_fixed -> []", any(grepl("^  protected_fixed: \\[\\]", readLines(out_cfg))))
tpl_bad <- file.path(tmp, "bad.yaml"); writeLines(c("cohort: {{cohort}}", "x: {{unknown}}"), tpl_bad)
check("unknown placeholder is an error", fails(write_cohort_config(tpl_bad, out_cfg, "c", d, "patient"), "unfilled placeholder"))
if (requireNamespace("yaml", quietly = TRUE)) {
  write_cohort_config(tpl, out_cfg, "test cohort", d, control_match_keys(n$control_for))
  y <- yaml::read_yaml(out_cfg)
  check("generated config parses as YAML", y$cohort == "test cohort" && identical(y$design$protected_fixed, "prior_peptide_vaccination") &&
          identical(y$controls$match_keys, c("patient", "week")))
}

# --- round trip: init -> edit -> read -> build_tube_table == file-name path -----------------
cat("round trip\n")
unlink(file.path(fcs, "odd_name.fcs"))
s5 <- init_sample_sheet(fcs, cfg)
n5 <- read_sample_sheet(write_sheet(s5))
tt_sheet <- build_tube_table(fcs, cfg, n5)
tt_name <- build_tube_table(fcs, cfg)
shared <- setdiff(intersect(names(tt_sheet), names(tt_name)), "path")
check("shared columns equal", isTRUE(all.equal(tt_sheet[shared], tt_name[shared], check.attributes = FALSE)))
check("group column equals prior_vacc", identical(tt_sheet$prior_peptide_vaccination, tt_name$prior_vacc))
check("extra sheet columns carried", all(c("visit", "control_for", "exclude", "note") %in% names(tt_sheet)))
check("edited subject overrides the file name", {
  s6 <- s5; s6$subject[s6$filename == "MTB011_1_V3W6_neg.fcs"] <- 99L
  t6 <- build_tube_table(fcs, cfg, read_sample_sheet(write_sheet(s6)))
  t6$patient[t6$filename == "MTB011_1_V3W6_neg.fcs"] == 99 && t6$patient_day[t6$filename == "MTB011_1_V3W6_neg.fcs"] == "99_2024-03-07"
})

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
