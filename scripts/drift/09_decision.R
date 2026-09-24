#!/usr/bin/env Rscript
# Step 9: combine posteriors with the decision rule (plan section 3). For every layer A statistic,
# layer BG (neg background) response and layer B response, for every instrument contrast and each
# acquisition day:
#   rule 1: P(|effect| > replicate SD) >= replicate_noise_prob  -> "detected"
#   rule 2: removing the effect flips positivity calls (layers B and BG counts)  -> "correct"
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

# replicate noise
# empirical logit with a half-count continuity correction: neg gates often hold 0 events, and a
# clipped logit of 0 percent (about -14) would swamp the replicate SD with the clipping constant
elogit <- function(k, n) qlogis((k + 0.5) / (n + 1))
rep_by <- function(d, value) do.call(rbind, lapply(split(d, d$response), function(x) data.frame(
  response = x$response[1], replicate_sd = replicate_sd(value(x), x$patient_day), stringsAsFactors = FALSE)))
rep_a <- layer_a_replicate_sd(read_layer_a(a$out, cfg, tt))   # NA when step 4 ran on negs only
fb <- read_layer_b(a$out)
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
rows <- list()
prim_idx <- split(seq_len(nrow(prim)), paste(prim$layer, prim$response, prim$term, prim$level, sep = "|"))
for (key in names(prim_idx)) {
  d <- prim[prim_idx[[key]], ]
  layer <- d$layer[1]; resp <- d$response[1]; term <- d$term[1]; lev <- d$level[1]
  db <- fb[fb$response == resp, ]
  if (layer == "A") {
    mk <- sub(":.*$", "", resp); st <- sub("^.*:", "", resp)
    rs <- rep_a$replicate_sd[rep_a$marker == mk & rep_a$statistic == st]
    if (!length(rs) || is.na(rs)) rs <- sigma_fallback("A", resp)
    flips <- NA
  } else if (layer == "BG" && sub(":.*$", "", resp) %in% names(cfg$parent_populations)) {   # neg counts
    rs <- lookup(rep_bg_count, resp)
    # stimulated rows whose matched neg sits on the affected day / instrument
    which_neg <- if (term == "instrument") db$neg_instrument %in% lev else db$neg_acq_day %in% lev
    flips <- count_background_flips(db, mean(d$value), which_neg, cfg$positivity$min_subtracted_pct, cfg$positivity$min_fold_over_neg)
  } else if (layer == "BG") {   # neg cytokine median intensity
    rs <- lookup(rep_bg_int, resp)
    if (is.na(rs)) rs <- sigma_fallback("BG", resp)
    flips <- NA
  } else {
    rs <- lookup(rep_b, resp)
    which <- if (term == "instrument") db$instrument == lev else db$acq_day == lev
    flips <- count_positivity_flips(db, mean(d$value), which, cfg$positivity$min_subtracted_pct, cfg$positivity$min_fold_over_neg)
  }
  dec <- drift_decision(d$value, rs, flips, cfg$replicate_noise_prob)
  sel <- if (layer == "BG") tt$is_neg else rep(TRUE, nrow(tt))
  n_tubes <- if (term == "instrument") sum(sel & tt$instrument == lev) else sum(sel & tt$acq_day == lev)
  rows[[length(rows) + 1L]] <- cbind(layer = layer, response = resp, term = term, level = lev, n_tubes = n_tubes, dec)
}
decisions <- do.call(rbind, rows)
decisions <- decisions[order(decisions$layer, decisions$term, -abs(decisions$effect_mean)), ]
write_table(decisions, file.path(a$out, "decisions.csv"))

# headline
cat("\n=== Instrument contrast ===\n")
print(decisions[decisions$term == "instrument", c("layer", "response", "effect_mean", "ci_low", "ci_high", "p_exceeds_noise", "positivity_flips", "verdict")], row.names = FALSE)
cat("\n=== Days detected above replicate noise ===\n")
dd <- decisions[decisions$term == "day" & grepl("^detected", decisions$verdict), ]
if (nrow(dd)) print(dd[, c("layer", "response", "level", "n_tubes", "effect_mean", "positivity_flips", "verdict")], row.names = FALSE) else cat("none\n")
cat("\n=== Background (neg) shifts ===\n")
bg <- decisions[decisions$layer == "BG" & (decisions$term == "instrument" | grepl("^detected", decisions$verdict)), ]
if (nrow(bg)) { print(bg[, c("response", "term", "level", "n_tubes", "effect_mean", "p_exceeds_noise", "positivity_flips", "verdict")], row.names = FALSE)
} else cat(if (any(decisions$layer == "BG")) "no background effect above replicate noise\n" else "no layer BG effects (step 6 background fits not run)\n")

# control matching: how the negs were matched, and the stimulated rows that needed a fallback
if ("neg_match" %in% names(fq)) {
  st <- fq[!fq$is_neg, ]
  cnt <- as.data.frame(table(parent = st$parent, neg_match = st$neg_match), responseName = "n", stringsAsFactors = FALSE)
  cnt <- cnt[cnt$n > 0, ]; cnt$patient <- NA; cnt$week <- NA; cnt$gate <- NA
  fbk <- unique(st[st$neg_match != "exact", c("patient", "week", "parent", "gate", "neg_match")])
  fbk$n <- rep(NA, nrow(fbk))
  write_table(rbind(cnt, fbk[, names(cnt)]), file.path(a$out, "controls_summary.csv"))
}

# gate drift diagnostic (rank 3): do CD4/CD8 gate bounds move with day or instrument?
gb_path <- file.path(a$out, "gate_bounds.csv")
if (file.exists(gb_path) && file.size(gb_path) > 5) {
  gb <- merge(read.csv(gb_path), tt[, c("filename", "patient", "acq_day", "instrument")], by = "filename")
  # a rectangle gate has one finite side per dimension (CD4+ is "above x", CD8+ "above y"),
  # so use whichever bound is finite as the gate position
  gb$position <- ifelse(is.na(gb$lower), gb$upper, gb$lower)
  gb <- gb[!is.na(gb$position) & !is.na(gb$marker), ]
  gd <- do.call(rbind, lapply(split(gb, list(gb$parent, gb$marker), drop = TRUE), function(g) {
    if (nrow(g) < 6) return(NULL)
    g$patient <- factor(g$patient); two <- length(unique(g$instrument)) > 1
    m <- if (two) lm(position ~ instrument + patient, data = g) else lm(position ~ patient, data = g)
    co <- summary(m)$coefficients; ic <- grep("^instrument", rownames(co))
    data.frame(parent = g$parent[1], marker = g$marker[1], n = nrow(g),
               instrument_shift_of_gate = if (length(ic)) co[ic, 1] else NA,
               p = if (length(ic)) co[ic, 4] else NA,
               sd_gate_across_days = sd(tapply(g$position, g$acq_day, mean)),
               sd_gate_within_day = sqrt(mean(tapply(g$position, g$acq_day, var), na.rm = TRUE)))
  }))
  write_table(gd, file.path(a$out, "gate_drift.csv"))
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
