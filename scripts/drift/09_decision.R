#!/usr/bin/env Rscript
# Step 9: combine posteriors with the decision rule (plan section 3, revised 2026-10-13 after
# docs/plans/drift_methods_review.md). For every layer A statistic, layer BG (neg background) response
# and layer B response, for every technical contrast (instrument, operator, ...: the fixed technical
# factors of the design) and each acquisition day:
#   rule 1, three states against a relevance margin (relevance_margin(): count CV for counts, the
#     replicate SD otherwise): "detected" when P(|effect| > margin) >= replicate_noise_prob,
#     "negligible" when P(|effect| < margin) >= replicate_noise_prob, else "inconclusive"
#   rule 2 (layers B and BG counts): calls flipped when the effect is removed, over posterior draws,
#     against the flips that counting noise alone produces on the same tubes -> "correct" only above it
# Also compares every fitter (brms, stan, bayesflow_v1, bayesflow_v2) against brms.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "09_decision.R --out DIR --config FILE [--primary brms]"
a <- parse_args(list(out = "", config = "", primary = "brms"), usage = usage)
cfg <- read_drift_config(a$config)
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)

eff <- read_effects(a$out)
if (is.null(eff)) stop("no effects_*.csv in ", a$out, "; run step 6 first")
models <- unique(eff$model)
cat("posteriors found from:", paste(models, collapse = ", "), "\n")
if (!a$primary %in% models) stop("primary model ", a$primary, " has no effects file")
prim <- eff[eff$model == a$primary, ]

# replicate noise. Technical markers (layer A): SD across the tubes of the same sample (patient, day,
# visit), the peptide being irrelevant to them. Counts (layers B and BG): the tube-level residual SD of
# the lme4 fit (observation-level term, after patient, week, pool and the technical terms), because
# tubes of one sample carry different patient-specific pools and so are not replicates; the
# same-sample empirical-logit SD is the fallback when step 6c has not run.
# empirical logit with a half-count continuity correction: neg gates often hold 0 events, and a
# clipped logit of 0 percent (about -14) would swamp the replicate SD with the clipping constant
elogit <- function(k, n) qlogis((k + 0.5) / (n + 1))
rep_by <- function(d, value) do.call(rbind, lapply(split(d, d$response), function(x) data.frame(
  response = x$response[1], replicate_sd = replicate_sd(value(x), same_sample(x)), stringsAsFactors = FALSE)))
lme_path <- file.path(a$out, "lme4_variance_components.csv")
lme_vc <- if (file.exists(lme_path)) read.csv(lme_path, stringsAsFactors = FALSE) else NULL
resid_sd <- function(layer, resp) {
  if (is.null(lme_vc)) return(NA)
  v <- lme_vc$var_resid[lme_vc$layer == layer & lme_vc$response == resp]
  if (length(v) && is.finite(v[1])) sqrt(v[1]) else NA
}
rep_a <- layer_a_replicate_sd(read_layer_a(a$out, cfg, tt))   # NA when step 4 ran on negs only
fb <- read_layer_b(a$out, cfg = cfg)
day_col <- design_day_term(cfg)
tube_level <- function(df, term) as.character(df[[if (term == "day") day_col else term]])
rep_b <- rep_by(fb, function(x) elogit(x$count, x$parent_count))
# layer BG: negs of the same patient-day (several exist where two visits were acquired together)
fq <- read.csv(file.path(a$out, "gated_frequencies.csv"), stringsAsFactors = FALSE)
ng <- fq[fq$is_neg & !is.na(fq$parent_count) & fq$parent_count > 0, ]
ng$response <- paste(ng$parent, ng$gate, sep = ":")
rep_bg_count <- rep_by(ng, function(x) elogit(x$count, x$parent_count))
ns <- merge(read.csv(file.path(a$out, "neg_summaries.csv"), stringsAsFactors = FALSE),
            tt[tt$is_neg, c("filename", "patient_day", "acq_day", "instrument")], by = "filename")
ns <- ns[ns$marker %in% cfg$cytokine_markers, ]; ns$response <- paste0(ns$marker, ":median")
rep_bg_int <- if (nrow(ns)) rep_by(ns, function(x) x$median) else NULL
# fallback for a missing replicate SD: posterior residual SD from the brms fit summary
sigma_fallback <- function(layer, resp) {
  fs <- file.path(a$out, "fits_brms_summary.csv"); if (!file.exists(fs)) return(NA)
  fs <- read.csv(fs); rs <- fs$sigma[fs$layer == layer & fs$response == resp]
  if (length(rs)) rs[1] else NA
}
lookup <- function(tab, resp) { rs <- if (is.null(tab)) numeric(0) else tab$replicate_sd[tab$response == resp]; if (length(rs)) rs[1] else NA }

# decisions
pos <- cfg$positivity
rows <- list()
prim_idx <- split(seq_len(nrow(prim)), paste(prim$layer, prim$response, prim$term, prim$level, sep = "|"))
for (key in names(prim_idx)) {
  d <- prim[prim_idx[[key]], ]
  layer <- d$layer[1]; resp <- d$response[1]; term <- d$term[1]; lev <- d$level[1]
  db <- fb[fb$response == resp, ]
  if (layer == "A") {
    mk <- sub(":.*$", "", resp); st <- sub("^.*:", "", resp)
    rs <- rep_a$replicate_sd[rep_a$marker == mk & rep_a$statistic == st]; src <- "same-sample tubes"
    if (!length(rs) || is.na(rs)) { rs <- sigma_fallback("A", resp); src <- "brms residual" }
    flips <- NA; counts <- FALSE
  } else if (layer == "BG" && sub(":.*$", "", resp) %in% names(cfg$parent_populations)) {   # neg counts
    rs <- resid_sd("BG", resp); src <- "lme4 tube residual"
    if (is.na(rs)) { rs <- lookup(rep_bg_count, resp); src <- "same-sample negs" }
    # stimulated rows whose matched neg sits on the affected day / instrument
    which_neg <- matched_neg_level(db, term, tt) %in% lev
    flips <- flip_summary(d$value, function(e) count_background_flips(db, e, which_neg, pos$min_subtracted_pct, pos$min_fold_over_neg),
                          function() counting_noise_flips(db, which_neg, pos$min_subtracted_pct, pos$min_fold_over_neg))
    counts <- TRUE
  } else if (layer == "BG") {   # neg cytokine median intensity
    rs <- lookup(rep_bg_int, resp); src <- "same-sample negs"
    if (is.na(rs)) { rs <- sigma_fallback("BG", resp); src <- "brms residual" }
    flips <- NA; counts <- FALSE
  } else {
    rs <- resid_sd("B", resp); src <- "lme4 tube residual"
    if (is.na(rs)) { rs <- lookup(rep_b, resp); src <- "same-sample tubes" }
    which <- tube_level(db, term) %in% lev
    flips <- flip_summary(d$value, function(e) count_positivity_flips(db, e, which, pos$min_subtracted_pct, pos$min_fold_over_neg),
                          function() counting_noise_flips(db, which, pos$min_subtracted_pct, pos$min_fold_over_neg))
    counts <- TRUE
  }
  mg <- relevance_margin(cfg, counts, rs)
  dec <- drift_decision(d$value, rs, flips, cfg$replicate_noise_prob, margin = mg$margin)
  sel <- if (layer == "BG") tt$is_neg else rep(TRUE, nrow(tt))
  n_tubes <- sum(sel & tube_level(tt, term) %in% lev)
  rows[[length(rows) + 1L]] <- cbind(layer = layer, response = resp, term = term, level = lev, n_tubes = n_tubes, dec,
                                     noise_source = src, margin_source = mg$source)
}
decisions <- do.call(rbind, rows)
decisions <- decisions[order(decisions$layer, decisions$term, -abs(decisions$effect_mean)), ]
write_table(decisions, file.path(a$out, "decisions.csv"))

# headline
cat("\n=== Verdicts per layer and term (rule 1: three states against the relevance margin) ===\n")
print(table(paste(decisions$layer, decisions$term), sub(",.*$", "", decisions$verdict)))
cat("\n=== Technical contrasts (instrument and other technical factors) ===\n")
print(decisions[decisions$term != "day", c("layer", "response", "term", "level", "effect_mean", "ci_low", "ci_high", "p_exceeds_noise", "positivity_flips", "verdict")], row.names = FALSE)
cat("\n=== Days detected above the relevance margin ===\n")
dd <- decisions[decisions$term == "day" & grepl("^detected", decisions$verdict), ]
if (nrow(dd)) print(dd[, c("layer", "response", "level", "n_tubes", "effect_mean", "positivity_flips", "verdict")], row.names = FALSE) else cat("none\n")
cat("\n=== Background (neg) shifts ===\n")
bg <- decisions[decisions$layer == "BG" & (decisions$term != "day" | grepl("^detected", decisions$verdict)), ]
if (nrow(bg)) { print(bg[, c("response", "term", "level", "n_tubes", "effect_mean", "p_exceeds_noise", "positivity_flips", "verdict")], row.names = FALSE)
} else cat(if (any(decisions$layer == "BG")) "no background effect above the relevance margin\n" else "no layer BG effects (step 6 background fits not run)\n")

# control matching: how the negs were matched, and the stimulated rows that needed a fallback
if ("neg_match" %in% names(fq)) {
  st <- fq[!fq$is_neg, ]
  cnt <- as.data.frame(table(parent = st$parent, neg_match = st$neg_match), responseName = "n", stringsAsFactors = FALSE)
  cnt <- cnt[cnt$n > 0, ]; cnt$patient <- NA; cnt$week <- NA; cnt$gate <- NA
  fbk <- unique(st[st$neg_match != "exact", c("patient", "week", "parent", "gate", "neg_match")])
  fbk$n <- rep(NA, nrow(fbk))
  write_table(rbind(cnt, fbk[, names(cnt)]), file.path(a$out, "controls_summary.csv"))
}

# gate drift diagnostic (rank 3): do the CD4/CD8 gates and the cytokine gates move with day or instrument?
# One row per parent and gate, on the gate's own dimension (`gate_marker`; files written before
# 2026-10-13 carry parent gates only, whose own dimension is the parent's marker: the CD4+ gate on CD4).
gb_path <- file.path(a$out, "gate_bounds.csv")
if (file.exists(gb_path) && file.size(gb_path) > 5) {
  gb <- merge(read.csv(gb_path, stringsAsFactors = FALSE), tt[, c("filename", "patient", "acq_day", "instrument")], by = "filename")
  if (!"gate" %in% names(gb)) { gb$gate <- gb$parent; gb$gate_marker <- gb$parent }
  gb <- gb[!is.na(gb$marker) & !is.na(gb$gate_marker) & gb$marker == gb$gate_marker, ]
  # a rectangle gate has one finite side per dimension (CD4+ is "above x"), so use whichever bound is finite
  gb$position <- ifelse(is.na(gb$lower), gb$upper, gb$lower)
  gb <- gb[!is.na(gb$position), ]
  gd <- do.call(rbind, lapply(split(gb, list(gb$parent, gb$gate_marker), drop = TRUE), function(g) {
    if (nrow(g) < 6) return(NULL)
    g$patient <- factor(g$patient); two <- length(unique(g$instrument)) > 1
    if (two) g$instrument <- relevel(factor(g$instrument), ref = cfg$reference_instrument)
    m <- if (two) lm(position ~ instrument + patient, data = g) else lm(position ~ patient, data = g)
    co <- summary(m)$coefficients; ic <- grep("^instrument", rownames(co))
    data.frame(parent = g$parent[1], gate = g$gate_marker[1], kind = if (g$gate_marker[1] == g$parent[1]) "parent" else "cytokine",
               n = nrow(g), instrument_shift_of_gate = if (length(ic)) co[ic[1], 1] else NA,
               p = if (length(ic)) co[ic[1], 4] else NA,
               sd_gate_across_days = sd(tapply(g$position, g$acq_day, mean)),
               sd_gate_within_day = sqrt(mean(tapply(g$position, g$acq_day, var), na.rm = TRUE)))
  }))
  write_table(gd, file.path(a$out, "gate_drift.csv"))
  cat("\n=== Gate positions (arcsinh): shift on the non-reference instrument and day-to-day movement ===\n")
  print(gd, digits = 3, row.names = FALSE)
}

# fitter agreement
cmp <- do.call(rbind, lapply(setdiff(models, a$primary), function(m) compare_effects(eff, a$primary, m)))
if (!is.null(cmp)) {
  write_table(cmp, file.path(a$out, "fitter_agreement.csv"))
  cat("\n=== Agreement with", a$primary, "(posterior mean inside its 95% interval) ===\n")
  print(agreement_table(cmp), row.names = FALSE)
}
# BayesFlow acceptance checks (plan section 6.3). With one neg tube per patient-day most day
# effects are only weakly identified by ANY method, so recovery R2 is reported but not gated on;
# what is gated is that the posteriors are calibrated (SBC, coverage) and agree with the yardstick.
for (v in c("v1", "v2")) {
  rf <- file.path(a$out, sprintf("bf_%s_recovery.csv", v)); sf <- file.path(a$out, sprintf("bf_%s_sbc.csv", v))
  if (file.exists(rf)) {
    rec <- read.csv(rf); sbc <- read.csv(sf)
    ag <- if (!is.null(cmp) && any(cmp$model_b == paste0("bayesflow_", v))) mean(cmp$agree[cmp$model_b == paste0("bayesflow_", v)]) else NA
    ppc <- if (file.exists(file.path(a$out, "bf_v2_ppc.csv")) && v == "v2") mean(tolower(read.csv(file.path(a$out, "bf_v2_ppc.csv"))$inside_90) == "true") else NA
    sbc_ok <- mean(sbc$sbc_p > 0.01); cov_ok <- mean(rec$coverage95 >= 0.90 & rec$coverage95 <= 0.99)
    cat(sprintf("\nBayesFlow %s: SBC p > 0.01 for %.0f%% of parameters (pass >= 90%%); 95%% coverage in [0.90, 0.99] for %.0f%% (pass >= 90%%); agreement with %s %.2f (pass >= 0.95)%s; median recovery R2 %.2f (information only)\n",
                v, 100 * sbc_ok, 100 * cov_ok, a$primary, ag, if (is.na(ppc)) "" else sprintf("; PPC coverage %.2f (pass >= 0.85)", ppc), median(rec$r2)))
    pass <- sbc_ok >= 0.9 && cov_ok >= 0.9 && isTRUE(ag >= 0.95) && (is.na(ppc) || ppc >= 0.85)
    cat(sprintf("  -> %s\n", if (pass) "PASSES: may drive the decision" else "FAILS: decision falls back to the yardstick"))
    write_table(data.frame(variant = v, sbc_frac_ok = sbc_ok, coverage_frac_ok = cov_ok, agreement = ag, ppc = ppc,
                           median_r2 = median(rec$r2), pass = pass), file.path(a$out, sprintf("bf_%s_acceptance.csv", v)))
  }
}
