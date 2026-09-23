#!/usr/bin/env Rscript
# Step 6c: frequentist cross-check. Variance components and ICC of day within patient for every
# layer A statistic (lmer), binomial GLMM for layer B (glmer). No priors, seconds to run.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_lme4.R --out DIR --config FILE"
a <- parse_args(list(out = "", config = ""), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("lme4")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
ref <- cfg$reference_instrument
rows <- list()

fa <- read_layer_a(a$out, cfg, tt, ref); fa <- fa[fa$is_neg, ]
for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
  d <- fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]
  m <- tryCatch(lme4::lmer(value ~ instrument + (1 | patient) + (1 | acq_day), data = d, REML = TRUE),
                error = function(e) NULL)
  if (is.null(m)) next
  rows[[length(rows) + 1L]] <- lme_variance_row(m, "A", key, nrow(d), resid_grp = "Residual")
}
fb <- read_layer_b(a$out, ref)
for (resp in unique(fb$response)) {
  d <- fb[fb$response == resp, ]
  m <- tryCatch(lme4::glmer(cbind(count, parent_count - count) ~ instrument + prior_vacc + week_c +
                              (1 | patient) + (1 | acq_day) + (1 | cocktail) + (1 | filename),   # obs-level RE = overdispersion
                            data = d, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa")),
                error = function(e) NULL)
  if (is.null(m)) next
  rows[[length(rows) + 1L]] <- lme_variance_row(m, "B", resp, nrow(d), resid_grp = "filename")
}
write_table(do.call(rbind, rows), file.path(a$out, "lme4_variance_components.csv"))
