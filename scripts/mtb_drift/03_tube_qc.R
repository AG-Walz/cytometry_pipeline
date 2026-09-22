#!/usr/bin/env Rscript
# Step 3: within-tube flow-stability QC (PeacoQC, flowAI fallback). Writes qc.csv and
# excluded_tubes.csv. A tube is excluded when it has fewer than --min_events events, when QC
# fails, or when PeacoQC flags more than --max_removed of its events. The QC flags are NOT applied
# to the data (frequencies come from the FlowJo gates on the full file); they only decide exclusion,
# so the default threshold is deliberately lenient: PeacoQC flags 30 to 50 percent of events on some
# perfectly gateable tubes here, and excluding a neg tube costs a whole patient-week.
#   Rscript scripts/mtb_drift/03_tube_qc.R --out DIR [--max_removed 0.6] [--min_events 1000] [--workers 8]
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "03_tube_qc.R --out DIR [--config FILE] [--max_removed 0.6] [--min_events 1000] [--workers 8]"
a <- parse_args(list(out = "", config = file.path(REPO, "config", "mtb_drift.yaml"),
                     max_removed = 0.6, min_events = 1000, workers = 8), usage = usage)
cfg <- read_drift_config(a$config)
need_packages(c("flowCore"))
use_peacoqc <- requireNamespace("PeacoQC", quietly = TRUE)
if (!use_peacoqc) need_packages("flowAI")

tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
pm <- read.csv(file.path(a$out, "panel_map.csv"), stringsAsFactors = FALSE)

qc_one <- function(i) {
  f <- tt$path[i]
  ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
  sp <- flowCore::spillover(ff); sp <- sp[[which(!vapply(sp, is.null, logical(1)))[1]]]
  if (!is.null(sp)) ff <- flowCore::compensate(ff, sp)
  chans <- pm$channel[pm$filename == tt$filename[i] & pm$canonical %in%
                        c(cfg$technical_markers, cfg$cytokine_markers) &
                        !pm$canonical %in% c("FSC_A", "SSC_A")]
  tl <- flowCore::transformList(chans, function(x) asinh(x / cfg$arcsinh_cofactor))
  fft <- flowCore::transform(ff, tl)
  n0 <- nrow(fft)
  res <- tryCatch({
    if (use_peacoqc) {
      q <- PeacoQC::PeacoQC(fft, channels = chans, plot = FALSE, save_fcs = FALSE,
                            report = FALSE, output_directory = NULL)
      list(kept = nrow(q$FinalFF), method = "PeacoQC")
    } else {
      q <- flowAI::flow_auto_qc(fft, output = 1, html_report = FALSE, mini_report = FALSE,
                                fcs_QC = FALSE, folder_results = FALSE)
      list(kept = nrow(q), method = "flowAI")
    }
  }, error = function(e) list(kept = NA_integer_, method = paste("failed:", conditionMessage(e))))
  data.frame(filename = tt$filename[i], n_events = n0, n_kept = res$kept,
             frac_removed = 1 - res$kept / n0, qc_method = res$method,
             stringsAsFactors = FALSE)
}

rows <- if (a$workers > 1 && requireNamespace("parallel", quietly = TRUE))
  parallel::mclapply(seq_len(nrow(tt)), qc_one, mc.cores = a$workers) else lapply(seq_len(nrow(tt)), qc_one)
qc <- do.call(rbind, rows)
qc$flag <- ifelse(qc$n_events < a$min_events, "too few events",
           ifelse(is.na(qc$frac_removed), "QC failed",
           ifelse(qc$frac_removed > a$max_removed, "unstable acquisition",
           ifelse(qc$frac_removed > 0.3, "check: >30% flagged by PeacoQC", ""))))
excl <- qc[qc$flag %in% c("too few events", "QC failed", "unstable acquisition"), ]
write_table(qc, file.path(a$out, "qc.csv"))
write_table(excl, file.path(a$out, "excluded_tubes.csv"))
cat(sprintf("\n%d of %d tubes excluded; %d more flagged for a look:\n", nrow(excl), nrow(qc), sum(startsWith(qc$flag, "check"))))
print(qc[nzchar(qc$flag), c("filename", "n_events", "frac_removed", "flag")], row.names = FALSE)
