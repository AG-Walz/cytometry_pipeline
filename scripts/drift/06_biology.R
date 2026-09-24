#!/usr/bin/env Rscript
# Step 6d: the biological effect. For every layer B endpoint the fixed week slope (change per year
# since first vaccination, logit scale) is estimated with glmer, once with the full design and once
# without the technical terms (instrument, acq_day), then with each acquisition day left out in
# turn. If the slope does not move, technical drift is not what the response is measuring.
# Writes week_slope_sensitivity.csv, week_slope_lodo.csv and week_slope_subjects.csv. Minutes to run.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_biology.R --out DIR --config FILE"
a <- parse_args(list(out = "", config = ""), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("lme4")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
check_design(cfg, tt)
ref <- cfg$reference_instrument
unit <- cfg$design$unit
day_term <- setdiff(unlist(cfg$design$technical), "instrument")[1]
if (is.na(day_term)) stop("the design names no day term among `technical`; nothing to leave out", call. = FALSE)

fit_glmer <- function(d, f, what) tryCatch(
  lme4::glmer(f, data = d, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa")),
  error = function(e) { cat("   ", what, "failed:", conditionMessage(e), "\n"); NULL })
# instrument is dropped from the formula when the subset holds a single instrument
lodo_formula <- function(d) {
  c2 <- cfg
  if ("instrument" %in% unlist(c2$design$technical) && nlevels(droplevels(d$instrument)) < 2)
    c2$design$technical <- setdiff(unlist(c2$design$technical), "instrument")
  biology_formula(c2)
}

f_full <- biology_formula(cfg); f_nt <- biology_formula(cfg, technical = FALSE)
cat("full:         ", deparse1(f_full), "\nno_technical: ", deparse1(f_nt), "\n")
fb <- read_layer_b(a$out, ref, cfg)
sens <- list(); lodo <- list(); subj <- list()
for (resp in unique(fb$response)) {
  d <- fb[fb$response == resp, ]
  cat(resp, ":", nrow(d), "tubes\n")
  m_full <- fit_glmer(d, f_full, "full")
  sens[[length(sens) + 1L]] <- list(response = resp, model = "full", slope = slope_of(m_full), n = nrow(d))
  m_nt <- fit_glmer(d, f_nt, "no_technical")
  sens[[length(sens) + 1L]] <- list(response = resp, model = "no_technical", slope = slope_of(m_nt), n = nrow(d))
  if (!is.null(m_full)) subj[[length(subj) + 1L]] <- subject_slopes(coef(m_full)[[unit]], resp)
  for (day in sort(unique(as.character(d[[day_term]])))) {
    keep <- d[[day_term]] != day
    d2 <- droplevels(d[keep, ])
    if ("instrument" %in% names(d2) && ref %in% levels(d2$instrument)) d2$instrument <- relevel(d2$instrument, ref = ref)
    m <- fit_glmer(d2, lodo_formula(d2), paste("without", day))
    lodo[[length(lodo) + 1L]] <- list(response = resp, dropped_day = day, slope = slope_of(m), n_tubes_dropped = sum(!keep))
  }
}
sens_df <- sensitivity_table(sens)
write_table(sens_df, file.path(a$out, "week_slope_sensitivity.csv"))
write_table(lodo_table(lodo), file.path(a$out, "week_slope_lodo.csv"))
empty_subj <- data.frame(response = character(0), patient = character(0), slope = numeric(0))
write_table(if (length(subj)) do.call(rbind, subj) else empty_subj, file.path(a$out, "week_slope_subjects.csv"))
cat("\nweek slope, logit per year (SE):\n")
print(sensitivity_pairs(sens_df)[, c("response", "week_slope", "se", "week_slope_nt", "se_nt", "shift_se")], digits = 3, row.names = FALSE)
