#!/usr/bin/env Rscript
# Step 15: sensitivity analyses of the drift verdict (docs/plans/drift_methods_review.md, sections 3 and 4).
# lme4 only, minutes. Each part writes one table and is skipped when its input is missing.
#   calendar  layer A and B refitted with a smooth of the acquisition date (natural spline, --calendar_df)
#             next to the day intercept: does gradual drift (laser ageing, tandem dyes, lots) explain the
#             day variance? -> sensitivity_calendar.csv
#   fixed     layer B on fixed_gate_frequencies.csv (step 11: one threshold per channel for every tube)
#             against the manual FlowJo gates: technical effects and variance components when nobody
#             moves a gate -> sensitivity_fixed_gates.csv
#   calls     probability that each lab-rule call changes under counting noise alone, and the calls
#             kept at several minimum parent counts -> call_stability.csv, call_stability_summary.csv,
#             min_parent_sensitivity.csv
#   ungated   layer A on neg summaries taken without the hand-set live gate (step 4 run with
#             neg_population "root" into --ungated_dir) against the gated ones -> sensitivity_ungated.csv
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "15_sensitivity.R --out DIR --config FILE [--calendar_df 4] [--ungated_dir DIR] [--min_parent 0,1000,2500,5000,10000] [--n_noise 500]"
a <- parse_args(list(out = "", config = "", calendar_df = 4, ungated_dir = "none", min_parent = "0,1000,2500,5000,10000", n_noise = 500),
                usage = usage)
cfg <- read_drift_config(a$config)
need_packages("lme4")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
check_design(cfg, tt)
ref <- cfg$reference_instrument; day_term <- design_day_term(cfg); unit <- cfg$design$unit
pos <- cfg$positivity
ctl <- lme4::glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
quiet <- function(expr) tryCatch(suppressMessages(suppressWarnings(expr)), error = function(e) { cat("   skipped:", conditionMessage(e), "\n"); NULL })
vc_of <- function(m, grp) { if (is.null(m)) return(NA_real_); v <- as.data.frame(lme4::VarCorr(m))
  x <- v$vcov[v$grp == grp & is.na(v$var2) & (is.na(v$var1) | v$var1 == "(Intercept)")]; if (length(x)) x[1] else NA_real_ }
fe_of <- function(m, pattern) { if (is.null(m)) return(c(NA, NA)); fe <- lme4::fixef(m); se <- sqrt(diag(as.matrix(vcov(m))))
  i <- grep(pattern, names(fe)); if (length(i)) c(unname(fe[i[1]]), unname(se[i[1]])) else c(NA, NA) }
inst_pat <- "^instrument"

# ---- calendar: smooth of the acquisition date next to the day intercept -------------------------------
cal_cols <- function(d) {
  x <- as.numeric(as.Date(as.character(d[[day_term]])))
  b <- splines::ns(x, df = a$calendar_df)
  for (j in seq_len(ncol(b))) d[[paste0("cal", j)]] <- as.numeric(scale(b[, j]))
  d
}
add_cal <- function(f) update(f, as.formula(paste(". ~ . +", paste0("cal", seq_len(a$calendar_df), collapse = " + "))))
cal_row <- function(layer, resp, m0, m1, resid) data.frame(layer = layer, response = resp,
  var_day = vc_of(m0, day_term), var_day_with_calendar = vc_of(m1, day_term),
  instrument = fe_of(m0, inst_pat)[1], instrument_with_calendar = fe_of(m1, inst_pat)[1],
  aic = if (is.null(m0)) NA else AIC(m0), aic_with_calendar = if (is.null(m1)) NA else AIC(m1), stringsAsFactors = FALSE)
rows <- list()
fa <- layer_a_rows(read_layer_a(a$out, cfg, tt, ref), cfg)
for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
  d <- cal_cols(fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]); d$value <- as.numeric(scale(d$value))
  f <- design_formula(cfg, "technical", "lme4", "value")
  m0 <- quiet(lme4::lmer(f, data = d, REML = FALSE)); m1 <- quiet(lme4::lmer(add_cal(f), data = d, REML = FALSE))
  rows[[length(rows) + 1L]] <- cal_row("A", key, m0, m1)
}
fb <- read_layer_b(a$out, ref, cfg)
for (resp in unique(fb$response)) {
  d <- cal_cols(fb[fb$response == resp, ]); f <- design_formula(cfg, "endpoint", "lme4")
  m0 <- quiet(lme4::glmer(f, data = d, family = binomial, control = ctl)); m1 <- quiet(lme4::glmer(add_cal(f), data = d, family = binomial, control = ctl))
  rows[[length(rows) + 1L]] <- cal_row("B", resp, m0, m1)
}
cal <- do.call(rbind, rows); cal$delta_aic <- cal$aic_with_calendar - cal$aic
write_table(cal, file.path(a$out, "sensitivity_calendar.csv"))
cat("\n=== Calendar-time smooth (", a$calendar_df, "df) next to the day intercept; layer A responses standardised ===\n", sep = "")
print(cal[, c("layer", "response", "var_day", "var_day_with_calendar", "delta_aic")], digits = 3, row.names = FALSE)

# ---- fixed: one threshold per channel instead of the manual gates --------------------------------------
fx_path <- file.path(a$out, "fixed_gate_frequencies.csv")
if (file.exists(fx_path)) {
  fx <- layer_b_frame(with_tube_columns(read.csv(fx_path, stringsAsFactors = FALSE), a$out, cfg), cfg, ref)
  fx$response <- sub(" fixed$", "", fx$response)
  tech <- design_technical_factors(cfg)
  rows <- list()
  for (resp in intersect(unique(fb$response), unique(fx$response))) for (src in c("manual", "fixed")) {
    d <- if (src == "manual") fb[fb$response == resp, ] else fx[fx$response == resp, ]
    m <- quiet(lme4::glmer(design_formula(cfg, "endpoint", "lme4"), data = d, family = binomial, control = ctl))
    r <- data.frame(response = resp, gates = src, n = nrow(d), var_patient = vc_of(m, unit), var_day = vc_of(m, day_term),
                    var_visit = vc_of(m, "visit"), var_tube = vc_of(m, "filename"), positive_calls = sum(d$positive %in% TRUE),
                    stringsAsFactors = FALSE)
    for (tf in tech) { e <- fe_of(m, paste0("^", tf)); r[[paste0(tf, "_effect")]] <- e[1]; r[[paste0(tf, "_se")]] <- e[2] }
    rows[[length(rows) + 1L]] <- r
  }
  fg <- do.call(rbind, rows)
  write_table(fg, file.path(a$out, "sensitivity_fixed_gates.csv"))
  cat("\n=== Manual FlowJo gates against one fixed threshold per channel (layer B, lme4, logit scale) ===\n")
  print(fg, digits = 3, row.names = FALSE)
} else cat("\nno fixed_gate_frequencies.csv (step 11): fixed-gate sensitivity skipped\n")

# ---- calls: counting-noise stability and minimum parent count -----------------------------------------
fq <- read.csv(file.path(a$out, "gated_frequencies.csv"), stringsAsFactors = FALSE)
st <- fq[!fq$is_neg & !is.na(fq$neg_pct) & !is.na(fq$parent_count) & fq$parent_count > 0, ]
st$response <- paste(st$parent, st$gate, sep = ":")
st$positive <- positivity_call(st$pct, st$neg_pct, pos$min_subtracted_pct, pos$min_fold_over_neg)$positive
st$p_flip <- call_flip_probability(st, pos$min_subtracted_pct, pos$min_fold_over_neg, n = a$n_noise, seed = cfg$seed)
write_table(st[, c("filename", "response", "patient", "week", "acq_day", "parent_count", "count", "neg_count", "neg_parent_count",
                   "pct", "neg_pct", "positive", "p_flip")], file.path(a$out, "call_stability.csv"))
cs <- do.call(rbind, lapply(split(st, st$response), function(d) data.frame(response = d$response[1], calls = nrow(d),
  positive = sum(d$positive), expected_flips = sum(d$p_flip), unstable_10 = sum(d$p_flip > 0.10), unstable_25 = sum(d$p_flip > 0.25),
  median_parent = median(d$parent_count), neg_le5_events = sum(d$neg_count <= 5), stringsAsFactors = FALSE)))
write_table(cs, file.path(a$out, "call_stability_summary.csv"))
cat("\n=== Lab-rule calls under counting noise alone (binomial resampling,", a$n_noise, "times) ===\n")
print(cs, digits = 3, row.names = FALSE)
mins <- as.numeric(strsplit(a$min_parent, ",")[[1]])
mp <- do.call(rbind, lapply(mins, function(m) { d <- st[st$parent_count >= m, ]
  data.frame(min_parent_events = m, calls = nrow(d), share_of_calls = nrow(d) / nrow(st), positive = sum(d$positive),
             expected_flips = sum(d$p_flip), unstable_10 = sum(d$p_flip > 0.10), stringsAsFactors = FALSE) }))
write_table(mp, file.path(a$out, "min_parent_sensitivity.csv"))
cat("\n=== Calls kept at a minimum parent event count (all responses) ===\n"); print(mp, digits = 3, row.names = FALSE)

# ---- ungated: technical markers without the hand-set live gate ----------------------------------------
if (a$ungated_dir != "none" && file.exists(file.path(a$ungated_dir, "neg_summaries.csv"))) {
  fu <- layer_a_rows(read_layer_a(a$ungated_dir, cfg, tt, ref), cfg)
  rows <- list()
  for (src in c("gated", "ungated")) {
    d0 <- if (src == "gated") fa else fu
    for (key in unique(paste(d0$marker, d0$statistic, sep = ":"))) {
      d <- d0[paste(d0$marker, d0$statistic, sep = ":") == key, ]
      m <- quiet(lme4::lmer(design_formula(cfg, "technical", "lme4", "value"), data = d, REML = TRUE))
      e <- fe_of(m, inst_pat)
      rows[[length(rows) + 1L]] <- data.frame(response = key, population = src, n = nrow(d), instrument = e[1], instrument_se = e[2],
        var_patient = vc_of(m, unit), var_day = vc_of(m, day_term), var_resid = if (is.null(m)) NA else sigma(m)^2, stringsAsFactors = FALSE)
    }
  }
  ug <- do.call(rbind, rows)
  write_table(ug, file.path(a$out, "sensitivity_ungated.csv"))
  cat("\n=== Technical markers in the hand-gated live cells against all events (lme4, raw units) ===\n")
  print(ug[order(ug$response, ug$population), ], digits = 3, row.names = FALSE)
} else cat("\nno --ungated_dir with neg_summaries.csv: ungated sensitivity skipped\n")
