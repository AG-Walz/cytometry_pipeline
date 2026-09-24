#!/usr/bin/env Rscript
# Pre-flight design check (docs/plans/drift_usability_plan.md, section 3): are the biological
# clock (`week`) and the technical clock (`acq_day`) separable for this cohort? Reads tubes.csv
# from step 1, prints the bridge report, the confounding index and a simulation-based power
# check (lme4), writes separability.csv and separability_summary.csv and ends with a verdict
# line. With --wsp_dir and CytoML installed it also checks that every tube is referenced by a
# workspace and that the configured gates exist there. Runs in minutes, before any model fit.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "check_design.R --out DIR --config FILE [--n_rep 10] [--wsp_dir DIR --fcs_dir DIR]"
a <- parse_args(list(out = "", config = "", n_rep = 10, wsp_dir = NA_character_, fcs_dir = NA_character_), usage = usage)
cfg <- read_drift_config(a$config)
tubes_csv <- file.path(a$out, "tubes.csv")
if (!file.exists(tubes_csv)) stop("no tubes.csv in ", a$out, ": run scripts/drift/01_tube_table.R first", call. = FALSE)
tt <- read.csv(tubes_csv, stringsAsFactors = FALSE)
check_design(cfg, tt)
tt <- tt[!tt$filename %in% excluded_tubes(a$out), ]
t0 <- Sys.time()

cat(sprintf("=== Design check: %s (%d tubes, %d subjects, %d acquisition days) ===\n\n",
            cfg$cohort, nrow(tt), length(unique(tt$patient)), length(unique(tt$acq_day))))

# --- 1. bridge report ------------------------------------------------------------------------
br <- bridge_report(tt)
cat("Bridge report\n")
cat(sprintf("  %d subject(s) span more than one day: %s\n", sum(br$day_bridge),
            if (any(br$day_bridge)) paste(br$patient[br$day_bridge], collapse = ", ") else "none"))
cat(sprintf("  %d subject(s) span more than one instrument: %s\n", sum(br$instr_bridge),
            if (any(br$instr_bridge)) paste(br$patient[br$instr_bridge], collapse = ", ") else "none"))

# --- 2. confounding index --------------------------------------------------------------------
ci <- confounding_index(tt)
cat("\nConfounding index (within-subject share of week variance explained by acq_day)\n")
cat(sprintf("  index %s | %d multi-week patient-day(s) | %d subject(s) on >1 day | %d on >1 instrument\n",
            format(round(ci$index, 4)), ci$n_multi_week_days, ci$n_subjects_multi_day, ci$n_subjects_multi_instrument))

# --- 3. marker consistency -------------------------------------------------------------------
pm_csv <- file.path(a$out, "panel_map.csv")
if (file.exists(pm_csv)) {
  pm <- read.csv(pm_csv, stringsAsFactors = FALSE)
  fl <- pm[!is.na(pm$canonical) & !pm$canonical %in% c("FSC_A", "SSC_A", "Time") & pm$filename %in% tt$filename, ]
  blank <- unique(fl$filename[is.na(fl$marker) | !nzchar(fl$marker)])
  noncanon <- unique(fl$filename[!is.na(fl$marker) & nzchar(fl$marker) & fl$marker != fl$canonical])
  cat("\nMarker names ($PnS) across", length(unique(fl$filename)), "files\n")
  cat(sprintf("  %d file(s) with blank $PnS on a fluorescence channel, %d with non-canonical spellings\n",
              length(blank), length(noncanon)))
  if (length(blank) + length(noncanon))
    cat("  the panel map resolves them by channel name; to write one spelling into the files run\n",
        "   Rscript scripts/harmonize_marker_names.R --canonical", a$config, "(see its header)\n")
} else cat("\nMarker names: no panel_map.csv in", a$out, "(step 2 not run), skipped\n")

# --- 4. workspaces: every tube referenced, configured gates present --------------------------
ws_note <- NA_character_
if (!is.na(a$wsp_dir)) {
  if (!requireNamespace("CytoML", quietly = TRUE) || !requireNamespace("flowWorkspace", quietly = TRUE)) {
    ws_note <- "CytoML not installed, workspace check skipped"
  } else if (is.na(a$fcs_dir)) {
    ws_note <- "--fcs_dir is needed with --wsp_dir, workspace check skipped"
  } else {
    sig_csv <- file.path(a$out, "panel_signatures.csv")
    sigs <- if (file.exists(sig_csv)) read.csv(sig_csv, stringsAsFactors = FALSE) else NULL
    wanted <- c(neg = cfg$neg_population, unlist(cfg$parent_populations))
    seen <- character(0); miss <- list(); empty <- list()
    wsps <- list_workspaces(a$wsp_dir)
    cat(sprintf("\nWorkspaces: %d under %s (reading the gates, this takes a while)\n", length(wsps), a$wsp_dir))
    for (w in wsps) for (gs in open_workspace(w, a$fcs_dir, sigs, quiet = TRUE)) for (s in flowWorkspace::sampleNames(gs)) {
      gh <- gs[[s]]; fn <- gh_filename(gh); seen <- c(seen, fn)
      for (k in names(wanted)) {
        node <- find_population(gh, wanted[[k]])
        if (is.na(node)) { miss[[length(miss) + 1L]] <- data.frame(filename = fn, population = wanted[[k]], workspace = basename(w)); next }
        n <- tryCatch(flowWorkspace::gh_pop_get_count(gh, node), error = function(e) NA_integer_)
        if (isTRUE(n == 0L)) empty[[length(empty) + 1L]] <- data.frame(filename = fn, population = node, workspace = basename(w))
      }
    }
    unref <- setdiff(tt$filename, seen)
    cat(sprintf("\n  %d of %d tubes referenced by a workspace", nrow(tt) - length(unref), nrow(tt)))
    if (length(unref)) cat(":\n  not referenced:", paste(unref, collapse = ", "), "\n") else cat("\n")
    if (length(miss)) { m <- do.call(rbind, miss); cat("  configured gates missing:\n"); print(m, row.names = FALSE) }
    else cat(sprintf("  gates %s present for every referenced tube\n", paste(wanted, collapse = ", ")))
    if (length(empty)) { e <- do.call(rbind, empty); cat("  gates present but empty (0 events):\n"); print(e, row.names = FALSE) }
    ws_note <- sprintf("%d unreferenced tube(s), %d missing gate(s), %d empty gate(s)", length(unref), length(miss), length(empty))
  }
  cat(" ", ws_note, "\n")
} else cat("\nWorkspaces: --wsp_dir not given, gate check skipped\n")

# --- 5. design power by simulation -----------------------------------------------------------
power <- NULL
if (requireNamespace("lme4", quietly = TRUE)) {
  cat(sprintf("\nDesign power: %d reps x 3 scenarios of one synthetic endpoint on the real design (lme4 glmer)\n", a$n_rep))
  power <- design_power(tt, cfg, n_rep = a$n_rep, seed = if (is.null(cfg$seed)) 1 else cfg$seed)
  write_table(power, file.path(a$out, "separability.csv"))
  ps <- power_summary(power)
  with(power[power$scenario == "a_both", ], cat(sprintf(
    "  a  both effects:   week slope %.2f (true %.2f), day SD %.2f (true %.2f), instrument %.2f (true %.2f); leakage week %.0f%%, day %.0f%%\n",
    mean(week_slope_hat, na.rm = TRUE), week_slope_true[1], mean(day_sd_hat, na.rm = TRUE), day_sd_true[1],
    mean(instr_hat, na.rm = TRUE), instr_true[1], 100 * ps$leakage_week, 100 * ps$leakage_day)))
  with(power[power$scenario == "b_week_only", ], cat(sprintf(
    "  b  week only:      day SD %.3f (true 0); false day detection (SD > 0.1 logit) %.0f%%\n",
    mean(day_sd_hat, na.rm = TRUE), 100 * ps$false_day_rate)))
  with(power[power$scenario == "c_day_only", ], cat(sprintf(
    "  c  day shift only: week slope %.3f logit/year (true 0), bias %.3f\n", mean(week_slope_hat, na.rm = TRUE), ps$week_bias)))
  cat(sprintf("  %d of %d fits succeeded, %d without lme4 warnings (a singular fit of the stimulation term is expected: the simulation has no cocktail effect)\n",
              ps$n_fit, nrow(power), sum(power$converged)))
} else cat("\nDesign power: lme4 not installed, simulation skipped (verdict from the index alone)\n")

# --- verdict -----------------------------------------------------------------------------------
v <- separability_verdict(ci, power)
ps <- if (is.null(power)) list(leakage_week = NA, leakage_day = NA, false_day_rate = NA, week_bias = NA, n_fit = 0L) else power_summary(power)
write_table(data.frame(cohort = cfg$cohort, n_tubes = nrow(tt), n_subjects = length(unique(tt$patient)), n_days = length(unique(tt$acq_day)),
  index = ci$index, n_multi_week_days = ci$n_multi_week_days, n_subjects_multi_day = ci$n_subjects_multi_day,
  n_subjects_multi_instrument = ci$n_subjects_multi_instrument, leakage_week = ps$leakage_week, leakage_day = ps$leakage_day,
  false_day_rate = ps$false_day_rate, week_bias = ps$week_bias, n_fit = ps$n_fit, workspace_check = ws_note, verdict = v$label,
  stringsAsFactors = FALSE), file.path(a$out, "separability_summary.csv"))
cat("\n", strwrap(v$explanation, 100, prefix = "  "), sep = "\n")
cat(sprintf("\n(%.0f s)\nVERDICT: %s\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), v$label))
