#!/usr/bin/env Rscript
# Step 2: map every file's channels to canonical markers; flag missing markers.
#   Rscript scripts/mtb_drift/02_panel_reconcile.R --out <results_dir>   (needs tubes.csv from step 1)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "02_panel_reconcile.R --out DIR [--config FILE]"
a <- parse_args(list(out = "", config = file.path(REPO, "config", "mtb_drift.yaml")), usage = usage)
cfg <- read_drift_config(a$config)

tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
pm <- build_panel_map(tt, cfg$markers)
required <- c(cfg$technical_markers, cfg$cytokine_markers)
gaps <- panel_gaps(pm$map, required)

write_table(pm$map, file.path(a$out, "panel_map.csv"))
write_table(pm$signature, file.path(a$out, "panel_signatures.csv"))
write_table(gaps, file.path(a$out, "panel_gaps.csv"))

cat("\nDistinct panel signatures:", length(unique(pm$signature$panel_signature)), "\n")
print(table(pm$signature$panel_signature))
if (nrow(gaps)) { cat("\nFiles missing required markers:\n"); print(gaps) } else cat("\nAll files carry every required marker.\n")
unmapped <- pm$map[is.na(pm$map$canonical) & !grepl("Time|-H$|-W$", pm$map$channel), ]
if (nrow(unmapped)) { cat("\nChannels not mapped to any canonical marker (check config$markers):\n")
  print(unique(unmapped[, c("channel", "marker")])) }
