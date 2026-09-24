#!/usr/bin/env Rscript
# Step 0: pre-fill the sample sheet from the FCS headers (and the file names when the config has a
# filename_regex) and write a cohort config from config/drift_template.yaml. The user completes the
# sheet (subject, visit, weeks, stimulation, group_* columns, exclusions) and the EDIT fields of the
# config before running step 1 with --sample_sheet.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "00_init_sheet.R --fcs_dir DIR --out DIR [--config CFG] [--cohort NAME]"
a <- parse_args(list(fcs_dir = "", out = "", config = "none", cohort = "none"), usage = usage)
cfg <- if (a$config == "none") list(tubes = list(), instrument_labels = list()) else read_drift_config(a$config)
cohort <- if (a$cohort == "none") basename(normalizePath(a$fcs_dir)) else a$cohort

sheet <- init_sample_sheet(a$fcs_dir, cfg)
invisible(ensure_dir(a$out))
sheet_path <- file.path(a$out, "sample_sheet.csv")
write.csv(sheet, sheet_path, row.names = FALSE, na = "")
cat("wrote", sheet_path, "(", nrow(sheet), "rows )\n")

bio <- c("subject", "visit", "weeks_since_first_vaccination", "stimulation", grep("^group_", names(sheet), value = TRUE))
blank <- vapply(sheet[bio], function(x) sum(is.na(x) | !nzchar(as.character(x))), integer(1))
full <- sum(rowSums(is.na(sheet[bio]) | sheet[bio] == "") == 0)
cat(sprintf("%d of %d rows have every biological field guessed", full, nrow(sheet)))
if (any(blank > 0)) cat("; blank cells per column: ", paste(names(blank)[blank > 0], blank[blank > 0], sep = "=", collapse = ", "))
cat("\n")

cfg_path <- file.path(a$out, paste0("drift_", cohort, ".yaml"))
write_cohort_config(file.path(REPO, "config", "drift_template.yaml"), cfg_path, cohort,
                    sheet_design(sheet), control_match_keys(sheet$control_for))
cat("wrote", cfg_path, "\n")
cat("\nNext: complete the sheet (subject, visit, weeks_since_first_vaccination, stimulation with\n",
    "'unstimulated' for controls, group_* covariates, exclude/note), then check the EDIT fields of\n",
    "the config (instrument_labels, markers, gate names, reference_instrument) and run\n",
    "  01_extract_metadata.R --fcs_dir", a$fcs_dir, "--out OUT --config", cfg_path, "--sample_sheet", sheet_path, "\n")
