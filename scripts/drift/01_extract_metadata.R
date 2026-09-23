#!/usr/bin/env Rscript
# Step 1: tube table (biology from a sample sheet or the file names, FCS technical keywords) and bridge report.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "01_extract_metadata.R --fcs_dir DIR --out DIR --config FILE [--sample_sheet CSV]"
a <- parse_args(list(fcs_dir = "", out = "", config = "", sample_sheet = "none"), usage = usage)
cfg <- read_drift_config(a$config)

sheet <- if (a$sample_sheet == "none") NULL else read.csv(a$sample_sheet, stringsAsFactors = FALSE)
tt <- build_tube_table(a$fcs_dir, cfg, sheet)
bad <- tt[is.na(tt$patient), "filename"]
if (length(bad)) cat("WARNING: files without metadata (unparsable name or not in the sample sheet):", paste(bad, collapse = ", "), "\n")

write_table(tt, file.path(a$out, "tubes.csv"))
write_table(bridge_report(tt), file.path(a$out, "bridge_report.csv"))

cat("\nTubes per instrument x acquisition day:\n")
print(table(tt$instrument, tt$acq_day, useNA = "ifany"))
if (!is.null(tt$prior_vacc)) {
  cat("\nPatients with inconsistent prior_vacc label:\n")
  inc <- tapply(tt$prior_vacc, tt$patient, function(x) length(unique(x)) > 1L)
  print(names(inc)[inc])
}
if (any(!tt$has_spill)) cat("\nWARNING: files without SPILL keyword:", sum(!tt$has_spill), "\n")
