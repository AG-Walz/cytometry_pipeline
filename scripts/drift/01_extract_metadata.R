#!/usr/bin/env Rscript
# Step 1: tube table (biology from a sample sheet or the file names, FCS technical keywords) and bridge report.
# With --sample_sheet the group_* columns of the sheet must be the config's design.protected_fixed
# (00_init_sheet.R writes a matching config); rows marked exclude stay in tubes.csv and are listed in
# excluded_manual.csv, which steps 4 and 5 honour through excluded_tubes().
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "01_extract_metadata.R --fcs_dir DIR --out DIR --config FILE [--sample_sheet CSV]"
a <- parse_args(list(fcs_dir = "", out = "", config = "", sample_sheet = "none"), usage = usage)
cfg <- read_drift_config(a$config)

sheet <- if (a$sample_sheet == "none") NULL else read_sample_sheet(a$sample_sheet)
if (!is.null(sheet)) {
  groups <- sheet_design(sheet)$protected_fixed
  pf <- as.character(unlist(cfg$design$protected_fixed))
  if (!setequal(groups, pf))
    stop("the sheet's group_* columns [", paste(groups, collapse = ", "), "] differ from design.protected_fixed [",
         paste(pf, collapse = ", "), "] in ", a$config,
         "; regenerate the config with 00_init_sheet.R or edit design.protected_fixed", call. = FALSE)
}
tt <- build_tube_table(a$fcs_dir, cfg, sheet)
bad <- tt[is.na(tt$patient), "filename"]
if (length(bad)) cat("WARNING: files without metadata (unparsable name or not in the sample sheet):", paste(bad, collapse = ", "), "\n")
if (!is.null(sheet)) {
  check_sheet_against_files(sheet, tt)
  missing <- setdiff(sheet$filename, tt$filename)
  if (length(missing)) cat("WARNING: sheet rows without an FCS file:", paste(missing, collapse = ", "), "\n")
  write_table(sheet[sheet$exclude %in% TRUE, c("filename", "note")], file.path(a$out, "excluded_manual.csv"))
}

write_table(tt, file.path(a$out, "tubes.csv"))
write_table(bridge_report(tt), file.path(a$out, "bridge_report.csv"))

cat("\nTubes per instrument (acquisition days, first to last):\n")
for (ins in sort(unique(tt$instrument))) { d <- tt[tt$instrument %in% ins, ]
  cat(sprintf("  %-16s %3d tubes on %2d day(s), %s to %s\n", ins, nrow(d), length(unique(d$acq_day)), min(d$acq_day, na.rm = TRUE), max(d$acq_day, na.rm = TRUE))) }
cat("  (full table: tubes.csv, columns instrument and acq_day)\n")
for (col in intersect(unlist(cfg$design$protected_fixed), names(tt))) {
  cat(sprintf("\nPatients with inconsistent %s label:\n", col))
  inc <- tapply(tt[[col]], tt$patient, function(x) length(unique(x)) > 1L)
  print(names(inc)[inc])
}
if (any(!tt$has_spill)) cat("\nWARNING: files without SPILL keyword:", sum(!tt$has_spill), "\n")
