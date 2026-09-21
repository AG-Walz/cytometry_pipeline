#!/usr/bin/env Rscript
# Step 6c: frequentist cross-check. Variance components and ICC of day within patient for every
# layer A statistic (lmer), binomial GLMM for layer B (glmer). No priors, seconds to run.
#   Rscript scripts/mtb_drift/06_lme4.R --out DIR
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_lme4.R --out DIR [--config FILE]"
a <- parse_args(list(out = "", config = file.path(REPO, "config", "mtb_drift.yaml")), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("lme4")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
ref <- cfg$reference_instrument
rows <- list()

fa <- layer_a_frame(read.csv(file.path(a$out, "neg_summaries.csv")),
                    read.csv(file.path(a$out, "neg_emd.csv")), tt, cfg$technical_markers)
fa <- fa[fa$is_neg, ]; fa$instrument <- relevel(fa$instrument, ref = ref)
for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
  d <- fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]
  m <- tryCatch(lme4::lmer(value ~ instrument + (1 | patient) + (1 | acq_day), data = d, REML = TRUE),
                error = function(e) NULL)
  if (is.null(m)) next
  vc <- as.data.frame(lme4::VarCorr(m)); v <- setNames(vc$vcov, vc$grp)
  fe <- lme4::fixef(m); se <- sqrt(diag(as.matrix(vcov(m))))
  ic <- grep("^instrument", names(fe))
  rows[[length(rows) + 1L]] <- data.frame(layer = "A", response = key, n = nrow(d),
    var_patient = v["patient"], var_day = v["acq_day"], var_resid = v["Residual"],
    icc_day = v["acq_day"] / sum(v), drift_to_biology = sqrt(v["acq_day"] / v["patient"]),
    instrument_effect = if (length(ic)) fe[ic] else NA, instrument_se = if (length(ic)) se[ic] else NA,
    stringsAsFactors = FALSE, row.names = NULL)
}
fb <- layer_b_frame(read.csv(file.path(a$out, "gated_frequencies.csv")))
fb$instrument <- relevel(fb$instrument, ref = ref)
for (resp in unique(fb$response)) {
  d <- fb[fb$response == resp, ]
  m <- tryCatch(lme4::glmer(cbind(count, parent_count - count) ~ instrument + prior_vacc + week_c +
                              (1 | patient) + (1 | acq_day) + (1 | cocktail) + (1 | filename),   # obs-level RE = overdispersion
                            data = d, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa")),
                error = function(e) NULL)
  if (is.null(m)) next
  vc <- as.data.frame(lme4::VarCorr(m)); v <- setNames(vc$vcov, vc$grp)
  fe <- lme4::fixef(m); se <- sqrt(diag(as.matrix(vcov(m)))); ic <- grep("^instrument", names(fe))
  rows[[length(rows) + 1L]] <- data.frame(layer = "B", response = resp, n = nrow(d),
    var_patient = v["patient"], var_day = v["acq_day"], var_resid = v["filename"],
    icc_day = v["acq_day"] / sum(v), drift_to_biology = sqrt(v["acq_day"] / v["patient"]),
    instrument_effect = if (length(ic)) fe[ic] else NA, instrument_se = if (length(ic)) se[ic] else NA,
    stringsAsFactors = FALSE, row.names = NULL)
}
write_table(do.call(rbind, rows), file.path(a$out, "lme4_variance_components.csv"))
