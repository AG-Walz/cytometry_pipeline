#!/usr/bin/env Rscript
# Step 14: figures and table comparing the data representations (raw / arcsinh MFI, %, arcsinh % and logit
# frequencies, each plain and stimulated minus unstimulated; the responder call is not a setting): which setting keeps the
# most biological and the least batch variance? Reads variance_partition_brms.csv (lme4 when brms is missing).
# Writes partition_grid.csv (per readout), partition_grid_summary.csv (medians per setting) and figures 01 to 16
# under figures/ as PNG and PDF (names say what each shows).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "partition_plots.R"))
usage <- "14_partition_figures.R --out DIR --config FILE"
a <- parse_args(list(out = "", config = ""), usage = usage)
need_packages("ggplot2")
vp <- read_partition(a$out)
if (is.null(vp)) stop("no variance_partition_*.csv in ", a$out, "; run step 12 first", call. = FALSE)
cat("engine:", vp$engine[1], "\n")
write_table(partition_grid(vp), file.path(a$out, "partition_grid.csv"))
s <- partition_grid(vp, summary = TRUE); write_table(s, file.path(a$out, "partition_grid_summary.csv"))
print(s, digits = 3, row.names = FALSE)
fig <- ensure_dir(file.path(a$out, "figures"))
unlink(list.files(fig, "\\.(png|pdf)$", full.names = TRUE))   # every figure is rebuilt; renamed ones must not linger
save <- function(p, name, w, h) for (ext in c("png", "pdf")) ggplot2::ggsave(file.path(fig, paste0(name, ".", ext)), p, width = w, height = h, dpi = 150)
n_k <- function(k) length(unique(vp$readout[vapply(vp$representation, function(r) PARTITION_REPS[[r]]$kind, "") == k]))
fq <- read.csv(file.path(a$out, "gated_frequencies.csv"), stringsAsFactors = FALSE); tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
if (!is.null(cfg_tm <- read_drift_config(a$config)$technical_markers)) PARTITION_CONTROL_MARKERS <- unlist(cfg_tm)
save(plot_explain_workflow(), "01_analysis_workflow", 15, 3.4)
save(plot_explain_example(fq, tt, "CD4:IFNg"), "02_CD4_IFNg_per_tube_by_acquisition_day", 12, 9)
save(plot_explain_scorecard(vp), "03_combinations_ranked_biology_vs_batch", 12, 9)
save(plot_combinations(vp), "04_combinations_kept_and_dropped", 15, 5.5)
save(plot_partition_paired(vp, "freq"), "05_frequency_biology_vs_batch_by_source", 15, 6)
save(plot_partition_paired(vp, "mfi"), "06_mfi_biology_vs_batch_by_source", 12, 8)
save(plot_partition_full(vp, "freq"), "07_frequency_all_sources_of_variance", 13, 6)
save(plot_partition_full(vp, "mfi"), "08_mfi_all_sources_of_variance", 13, 6)
save(plot_partition_paired(vp, "freq", by_readout = TRUE), "09_frequency_per_readout_biology_vs_batch", 15, 1.9 * n_k("freq") + 2)
save(plot_partition_paired(vp, "mfi", by_readout = TRUE), "10_mfi_per_readout_biology_vs_batch_with_negative_controls", 12, 1.9 * n_k("mfi") + 2)
# acquisition audit (step 11b)
rdo <- function(f) { p <- file.path(a$out, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) }
fo <- rdo("factor_overlap.csv"); if (!is.null(fo)) { M <- as.matrix(fo[, -1]); rownames(M) <- fo[[1]]; save(plot_factor_overlap(M), "11_sources_overlap", 10, 8) }
acq <- rdo("acquisition_settings.csv"); if (!is.null(acq) && any(grepl("^V_", names(acq)))) save(plot_voltages(acq, tt), "12_detector_voltages_over_time", 13, 6)
gp <- rdo("gate_positions.csv"); if (!is.null(gp)) save(plot_gate_positions(gp, tt), "13_cytokine_gate_positions_over_time", 13, 6)
ea <- plot_engine_agreement(a$out)
if (!is.null(ea)) { save(ea$plot, "14_engines_brms_vs_lme4_variance_shares", 14, 11); write_table(ea$table, file.path(a$out, "partition_engine_agreement.csv"))
  print(ea$table, digits = 3, row.names = FALSE) }
es <- plot_engine_stan(a$out)
if (!is.null(es)) { save(es$plot, "15_engines_brms_vs_stan_correction_effects", 12, 9); print(es$table, row.names = FALSE) }
bh <- plot_bayesflow_history(a$out)
if (!is.null(bh)) save(bh$plot, paste0("16_bayesflow_vs_brms_run_", bh$date), 12, 7)
cat("figures in", fig, "\n")
