# src/harmonization/drift/drift_models.R — model frames and posterior handling
# Base R. Builds the two model frames (layer A: neg summaries; layer B: gated
# counts) from the step 4/5 tables, and turns posterior draws of the day and
# instrument effects into the decision table. brms, Stan, lme4 and BayesFlow
# all consume/produce the same shapes, so they are comparable in step 9.

#' Layer A frame: one row per neg tube x marker x statistic.
#' Statistics used: median, spread (q90-q10), valley and frac_pos (bimodal
#' markers only), emd. `summaries`, `emd`, `tubes` are the CSVs of steps 4 and 1.
layer_a_frame <- function(summaries, emd, tubes, technical_markers) {
  s <- summaries[summaries$marker %in% technical_markers, ]
  long <- rbind(
    data.frame(filename = s$filename, marker = s$marker, statistic = "median", value = s$median),
    data.frame(filename = s$filename, marker = s$marker, statistic = "spread", value = s$spread))
  if ("valley" %in% names(s)) {
    v <- s[!is.na(s$valley), ]
    long <- rbind(long,
      data.frame(filename = v$filename, marker = v$marker, statistic = "valley", value = v$valley),
      data.frame(filename = v$filename, marker = v$marker, statistic = "frac_pos", value = v$frac_pos))
  }
  e <- emd[emd$marker %in% technical_markers, ]
  long <- rbind(long, data.frame(filename = e$filename, marker = e$marker, statistic = "emd", value = e$emd))
  fr <- merge(long, tubes[, c("filename", "patient", "prior_vacc", "week", "is_neg", "acq_day",
                              "instrument", "patient_day", "cocktail")], by = "filename")
  fr$patient <- factor(fr$patient); fr$acq_day <- factor(fr$acq_day)
  fr$instrument <- factor(fr$instrument)
  fr[order(fr$marker, fr$statistic, fr$patient, fr$acq_day), ]
}

#' Layer B frame: cocktail tubes only, one row per tube x parent x gate.
layer_b_frame <- function(freqs) {
  fb <- freqs[!freqs$is_neg & !is.na(freqs$neg_pct), ]
  fb$patient <- factor(fb$patient); fb$acq_day <- factor(fb$acq_day)
  fb$instrument <- factor(fb$instrument); fb$cocktail <- factor(fb$cocktail)
  fb$prior_vacc <- factor(fb$prior_vacc)
  fb$week_c <- (fb$week - mean(fb$week)) / 52   # years since cohort-mean week, keeps slopes O(1)
  fb$response <- paste(fb$parent, fb$gate, sep = ":")
  fb
}

#' The two frames from the step 4/5 tables in `out`. With `ref` the instrument
#' factor is relevelled so that the fixed effect is the contrast against it.
read_layer_a <- function(out, cfg, tubes, ref = NULL) {
  fa <- layer_a_frame(read.csv(file.path(out, "neg_summaries.csv")),
                      read.csv(file.path(out, "neg_emd.csv")), tubes, cfg$technical_markers)
  if (!is.null(ref)) fa$instrument <- relevel(fa$instrument, ref = ref)
  fa
}
read_layer_b <- function(out, ref = NULL) {
  fb <- layer_b_frame(read.csv(file.path(out, "gated_frequencies.csv")))
  if (!is.null(ref)) fb$instrument <- relevel(fb$instrument, ref = ref)
  fb
}

#' Replicate SD per marker x statistic from ALL tubes of the same patient-day
#' (technical markers should not depend on the peptide). Needs step 4 run with
#' --all_tubes; with negs only, returns NA and step 9 falls back to the
#' posterior residual SD.
layer_a_replicate_sd <- function(frame_all_tubes) {
  sp <- split(frame_all_tubes, list(frame_all_tubes$marker, frame_all_tubes$statistic), drop = TRUE)
  do.call(rbind, lapply(sp, function(d) data.frame(
    marker = d$marker[1], statistic = d$statistic[1],
    replicate_sd = replicate_sd(d$value, d$patient_day), stringsAsFactors = FALSE)))
}

#' Count positivity calls that change when an effect on the logit scale is
#' removed from the tubes it applies to. `fb` is a layer B frame for ONE
#' response; `effect` is a scalar (posterior mean) applied to rows in `which`.
count_positivity_flips <- function(fb, effect, which, min_subtracted_pct, min_fold_over_neg) {
  before <- positivity_call(fb$pct, fb$neg_pct, min_subtracted_pct, min_fold_over_neg)$positive
  p <- fb$pct / 100; p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  adj <- p
  adj[which] <- plogis(qlogis(p[which]) - effect)
  # the neg tube of the same day/instrument is shifted the same way
  pn <- fb$neg_pct / 100; pn <- pmin(pmax(pn, 1e-6), 1 - 1e-6)
  adjn <- pn; adjn[which] <- plogis(qlogis(pn[which]) - effect)
  after <- positivity_call(100 * adj, 100 * adjn, min_subtracted_pct, min_fold_over_neg)$positive
  sum(before != after, na.rm = TRUE)
}

#' Standard shape for posterior draws of the effects, shared by every fitter:
#' data.frame(model, layer, response, term, level, draw, value)
#' term is "instrument" (level = non-reference instrument) or "day" (level = acq_day).
effects_long <- function(model, layer, response, instrument_draws, day_draws, instrument_level) {
  out <- data.frame(model = model, layer = layer, response = response, term = "instrument",
                    level = instrument_level, draw = seq_along(instrument_draws),
                    value = instrument_draws, stringsAsFactors = FALSE)
  for (d in colnames(day_draws)) {
    out <- rbind(out, data.frame(model = model, layer = layer, response = response, term = "day",
                                 level = d, draw = seq_len(nrow(day_draws)), value = day_draws[, d],
                                 stringsAsFactors = FALSE))
  }
  out
}

#' All effects_*.csv written by the step 6 fitters, or NULL when there are none.
read_effects <- function(out) {
  files <- list.files(out, "^effects_.*\\.csv$", full.names = TRUE)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

#' Variance components and the instrument fixed effect of one lme4 fit as a
#' summary row. `resid_grp` is the VarCorr group holding the residual variance
#' ("Residual" for lmer, the observation-level term for glmer).
lme_variance_row <- function(m, layer, response, n, resid_grp) {
  vc <- as.data.frame(lme4::VarCorr(m)); v <- setNames(vc$vcov, vc$grp)
  fe <- lme4::fixef(m); se <- sqrt(diag(as.matrix(vcov(m))))
  ic <- grep("^instrument", names(fe))
  data.frame(layer = layer, response = response, n = n,
    var_patient = v["patient"], var_day = v["acq_day"], var_resid = v[resid_grp],
    icc_day = v["acq_day"] / sum(v), drift_to_biology = sqrt(v["acq_day"] / v["patient"]),
    instrument_effect = if (length(ic)) fe[ic] else NA, instrument_se = if (length(ic)) se[ic] else NA,
    stringsAsFactors = FALSE, row.names = NULL)
}

#' Number of divergent transitions of a brms fit.
n_divergent <- function(fit) {
  np <- brms::nuts_params(fit)
  sum(np$Value[np$Parameter == "divergent__"])
}

#' Summary row of one brms fit: group SDs and sigma (multiplied by `scale` to
#' undo a standardised response; sigma is NA for families without one) plus
#' convergence diagnostics.
brms_fit_row <- function(fit, layer, response, n, scale = 1) {
  s <- as.data.frame(brms::VarCorr(fit)$acq_day$sd); sp <- as.data.frame(brms::VarCorr(fit)$patient$sd)
  spec <- summary(fit)$spec_pars
  data.frame(layer = layer, response = response, n = n,
    tau_day = s$Estimate * scale, tau_patient = sp$Estimate[1] * scale,
    sigma = if ("sigma" %in% rownames(spec)) spec["sigma", "Estimate"] * scale else NA,
    max_rhat = max(brms::rhat(fit), na.rm = TRUE), n_divergent = n_divergent(fit),
    stringsAsFactors = FALSE)
}

#' Agreement of two fitters on the same effects: does the posterior mean of B
#' fall inside the 95 percent interval of A, per response x term x level.
compare_effects <- function(eff, model_a, model_b) {
  key <- c("layer", "response", "term", "level")
  agg <- function(m) {
    d <- eff[eff$model == m, ]
    sp <- split(d, d[key], drop = TRUE)
    do.call(rbind, lapply(sp, function(x) cbind(x[1, key], mean = mean(x$value),
      lo = unname(quantile(x$value, 0.025)), hi = unname(quantile(x$value, 0.975)))))
  }
  a <- agg(model_a); b <- agg(model_b)
  names(b)[5:7] <- c("mean_b", "lo_b", "hi_b")
  m <- merge(a, b, by = key)
  m$model_a <- model_a; m$model_b <- model_b
  m$agree <- m$mean_b >= m$lo & m$mean_b <= m$hi
  m
}

#' Agreement counts "agreed/total" per model x layer x term.
agreement_table <- function(cmp) {
  aggregate(agree ~ model_b + layer + term, cmp, function(x) sprintf("%d/%d", sum(x), length(x)))
}
