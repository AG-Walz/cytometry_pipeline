#!/usr/bin/env Rscript
# Step 1: tube table (file-name biology + FCS technical keywords) and bridge report.
#   Rscript scripts/mtb_drift/01_extract_metadata.R --fcs_dir <dir> --out <results_dir>
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "01_extract_metadata.R --fcs_dir DIR --out DIR"
a <- parse_args(list(fcs_dir = "", out = ""), usage = usage)

tt <- build_tube_table(a$fcs_dir)
bad <- tt[is.na(tt$patient), "filename"]
if (length(bad)) cat("WARNING: unparsable file names:", paste(bad, collapse = ", "), "\n")

write_table(tt, file.path(a$out, "tubes.csv"))
write_table(bridge_report(tt), file.path(a$out, "bridge_report.csv"))

cat("\nTubes per instrument x acquisition day:\n")
print(table(tt$instrument, tt$acq_day, useNA = "ifany"))
cat("\nPatients with inconsistent prior_vacc label:\n")
inc <- tapply(tt$prior_vacc, tt$patient, function(x) length(unique(x)) > 1L)
print(names(inc)[inc])
if (any(!tt$has_spill)) cat("\nWARNING: files without SPILL keyword:", sum(!tt$has_spill), "\n")
