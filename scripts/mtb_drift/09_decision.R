#!/usr/bin/env Rscript
# Step 9: combine posteriors with the decision rule (plan section 3). For every layer A statistic
# and layer B response, and for the instrument contrast and each acquisition day:
#   rule 1: P(|effect| > replicate SD) >= replicate_noise_prob  -> "detected"
#   rule 2: removing the effect flips positivity calls (layer B)  -> "correct"
# Also compares every fitter (brms, stan, bayesflow_v1, bayesflow_v2) against brms.
#   Rscript scripts/mtb_drift/09_decision.R --out DIR [--primary brms]
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "09_decision.R --out DIR [--config FILE] [--primary brms]"
a <- parse_args(list(out = "", config = file.path(REPO, "config", "mtb_drift.yaml"), primary = "brms"), usage = usage)
cfg <- read_drift_config(a$config)
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)

eff_files <- list.files(a$out, "^effects_.*\\.csv$", full.names = TRUE)
if (!length(eff_files)) stop("no effects_*.csv in ", a$out, "; run step 6 first")
eff <- do.call(rbind, lapply(eff_files, read.csv, stringsAsFactors = FALSE))
models <- unique(eff$model)
cat("posteriors found from:", paste(models, collapse = ", "), "\n")
if (!a$primary %in% models) stop("primary model ", a$primary, " has no effects file")
prim <- eff[eff$model == a$primary, ]

# --- replicate noise -----------------------------------------------------------------------
fa_all <- layer_a_frame(read.csv(file.path(a$out, "neg_summaries.csv")), read.csv(file.path(a$out, "neg_emd.csv")),
                        tt, cfg$technical_markers)
rep_a <- layer_a_replicate_sd(fa_all)                       # NA when step 4 ran on negs only
fb <- layer_b_frame(read.csv(file.path(a$out, "gated_frequencies.csv")))
rep_b <- do.call(rbind, lapply(split(fb, fb$response), function(d) data.frame(
  response = d$response[1], replicate_sd = replicate_sd(qlogis(pmin(pmax(d$pct / 100, 1e-6), 1 - 1e-6)), d$patient_day))))

# --- decisions -----------------------------------------------------------------------------
rows <- list()
for (key in unique(paste(prim$layer, prim$response, prim$term, prim$level, sep = "|"))) {
  d <- prim[paste(prim$layer, prim$response, prim$term, prim$level, sep = "|") == key, ]
  layer <- d$layer[1]; resp <- d$response[1]; term <- d$term[1]; lev <- d$level[1]
  if (layer == "A") {
    mk <- sub(":.*$", "", resp); st <- sub("^.*:", "", resp)
    rs <- rep_a$replicate_sd[rep_a$marker == mk & rep_a$statistic == st]
    if (!length(rs) || is.na(rs)) {  # fallback: posterior residual SD from the fit summary
      fs <- read.csv(file.path(a$out, "fits_brms_summary.csv")); rs <- fs$sigma[fs$layer == "A" & fs$response == resp]
      if (!length(rs)) rs <- NA }
    flips <- NA
  } else {
    rs <- rep_b$replicate_sd[rep_b$response == resp]
    db <- fb[fb$response == resp, ]
    which <- if (term == "instrument") db$instrument == lev else db$acq_day == lev
    flips <- count_positivity_flips(db, mean(d$value), which, cfg$positivity$min_subtracted_pct, cfg$positivity$min_fold_over_neg)
  }
  dec <- drift_decision(d$value, rs, flips, cfg$replicate_noise_prob)
  n_tubes <- if (term == "instrument") sum(tt$instrument == lev) else sum(tt$acq_day == lev)
  rows[[length(rows) + 1L]] <- cbind(layer = layer, response = resp, term = term, level = lev, n_tubes = n_tubes, dec)
}
decisions <- do.call(rbind, rows)
decisions <- decisions[order(decisions$layer, decisions$term, -abs(decisions$effect_mean)), ]
write_table(decisions, file.path(a$out, "decisions.csv"))

# --- headline -----------------------------------------------------------------------------
cat("\n=== Instrument contrast ===\n")
print(decisions[decisions$term == "instrument", c("layer", "response", "effect_mean", "ci_low", "ci_high", "p_exceeds_noise", "positivity_flips", "verdict")], row.names = FALSE)
cat("\n=== Days detected above replicate noise ===\n")
dd <- decisions[decisions$term == "day" & grepl("^detected", decisions$verdict), ]
if (nrow(dd)) print(dd[, c("layer", "response", "level", "n_tubes", "effect_mean", "positivity_flips", "verdict")], row.names = FALSE) else cat("none\n")

# --- gate drift diagnostic (rank 3): do CD4/CD8 gate bounds move with day or instrument? --
gb_path <- file.path(a$out, "gate_bounds.csv")
if (file.exists(gb_path) && file.size(gb_path) > 5) {
  gb <- merge(read.csv(gb_path), tt[, c("filename", "patient", "acq_day", "instrument")], by = "filename")
  gd <- do.call(rbind, lapply(split(gb, list(gb$parent, gb$marker), drop = TRUE), function(g) {
    if (nrow(g) < 6) return(NULL)
    m <- lm(lower ~ instrument + patient, data = transform(g, patient = factor(patient)))
    co <- summary(m)$coefficients; ic <- grep("^instrument", rownames(co))
    data.frame(parent = g$parent[1], marker = g$marker[1], n = nrow(g),
               instrument_shift_of_lower_bound = if (length(ic)) co[ic, 1] else NA,
               p = if (length(ic)) co[ic, 4] else NA, sd_lower_across_days = sd(tapply(g$lower, g$acq_day, mean)))
  }))
  write_table(gd, file.path(a$out, "gate_drift.csv"))
}

# --- fitter agreement ---------------------------------------------------------------------
cmp <- do.call(rbind, lapply(setdiff(models, a$primary), function(m) compare_effects(eff, a$primary, m)))
if (!is.null(cmp)) {
  write_table(cmp, file.path(a$out, "fitter_agreement.csv"))
  cat("\n=== Agreement with", a$primary, "(posterior mean inside its 95% interval) ===\n")
  print(aggregate(agree ~ model_b + layer + term, cmp, function(x) sprintf("%d/%d", sum(x), length(x))), row.names = FALSE)
}
# BayesFlow acceptance checks (plan section 6.3)
for (v in c("v1", "v2")) {
  rf <- file.path(a$out, sprintf("bf_%s_recovery.csv", v)); sf <- file.path(a$out, sprintf("bf_%s_sbc.csv", v))
  if (file.exists(rf)) {
    rec <- read.csv(rf); sbc <- read.csv(sf)
    ag <- if (!is.null(cmp)) mean(cmp$agree[cmp$model_b == paste0("bayesflow_", v)]) else NA
    ppc <- if (file.exists(file.path(a$out, "bf_v2_ppc.csv")) && v == "v2") mean(tolower(read.csv(file.path(a$out, "bf_v2_ppc.csv"))$inside_90) == "true") else NA
    cat(sprintf("\nBayesFlow %s: min R2 %.2f (pass >= 0.9), min SBC p %.3f (pass > 0.05), agreement with %s %.2f (pass = 1.00)%s\n",
                v, min(rec$r2), min(sbc$sbc_p), a$primary, ag, if (is.na(ppc)) "" else sprintf(", PPC coverage %.2f (pass >= 0.85)", ppc)))
    pass <- min(rec$r2) >= 0.9 && min(sbc$sbc_p) > 0.05 && isTRUE(ag >= 0.999) && (is.na(ppc) || ppc >= 0.85)
    cat(sprintf("  -> %s\n", if (pass) "PASSES: may drive the decision" else "FAILS: decision falls back to the yardstick"))
  }
}
