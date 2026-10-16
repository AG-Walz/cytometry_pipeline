#!/usr/bin/env Rscript
# Step 6c: frequentist cross-check. Variance components and ICC of day within patient for every
# layer A statistic (lmer), layer BG response (glmer counts, lmer intensity) and layer B response
# (glmer), with the formulas of design_formula(). Also the layer A assumption check: every
# technical response is refitted with each protected covariate and the time slope added, and
# tested against the technical model by likelihood ratio (assumption_check.csv). Seconds to run.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_lme4.R --out DIR --config FILE"
a <- parse_args(list(out = "", config = ""), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("lme4")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
check_design(cfg, tt)
ref <- cfg$reference_instrument
unit <- cfg$design$unit; day_term <- setdiff(unlist(cfg$design$technical), "instrument")[1]
rows <- list(); checks <- list()
try_fit <- function(expr) tryCatch(expr, error = function(e) { cat("   skipped:", conditionMessage(e), "\n"); NULL })
fit_lmer <- function(d, class) try_fit(lme4::lmer(design_formula(cfg, class, "lme4", "value"), data = d, REML = TRUE))
fit_glmer <- function(d, class) try_fit(lme4::glmer(design_formula(cfg, class, "lme4"), data = d, family = binomial,
                                                    control = lme4::glmerControl(optimizer = "bobyqa")))
add_row <- function(m, layer, response, n, resid_grp)
  if (!is.null(m)) rows[[length(rows) + 1L]] <<- lme_variance_row(m, layer, response, n, resid_grp, unit, day_term)

# Layer A, with the assumption check: a technical response must not depend on the biology
covariates <- c(unlist(cfg$design$protected_fixed), if (!is.null(cfg$design$time)) "week_c")
fa <- read_layer_a(a$out, cfg, tt, ref); fa <- fa[fa$is_neg, ]
for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
  d <- fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]
  m0 <- fit_lmer(d, "technical")
  if (is.null(m0)) next
  add_row(m0, "A", key, nrow(d), "Residual")
  for (cv in covariates) {
    m1 <- try_fit(update(m0, as.formula(paste(". ~ . +", cv))))
    if (is.null(m1)) next
    p <- anova(m0, m1, refit = TRUE)[["Pr(>Chisq)"]][2]   # both refitted with ML
    fe <- lme4::fixef(m1); se <- sqrt(diag(as.matrix(vcov(m1))))
    ic <- grep(paste0("^", cv), names(fe))
    checks[[length(checks) + 1L]] <- data.frame(layer = "A", response = key, covariate = names(fe)[ic],
                                                estimate = unname(fe[ic]), se = unname(se[ic]), p_lrt = p,
                                                stringsAsFactors = FALSE)
  }
}
# Layer BG: cytokine background in the negs
bc <- read_layer_bg_counts(a$out, cfg, ref)
for (resp in unique(bc$response)) {
  d <- bc[bc$response == resp, ]
  add_row(fit_glmer(d, "background_counts"), "BG", resp, nrow(d), "filename")
}
bi <- read_layer_bg_intensity(a$out, cfg, tt, ref)
for (resp in unique(bi$response)) {
  d <- bi[bi$response == resp, ]
  add_row(fit_lmer(d, "background_intensity"), "BG", resp, nrow(d), "Residual")
}
# Layer B
fb <- read_layer_b(a$out, ref, cfg)
for (resp in unique(fb$response)) {
  d <- fb[fb$response == resp, ]
  add_row(fit_glmer(d, "endpoint"), "B", resp, nrow(d), "filename")
}
write_table(do.call(rbind, rows), file.path(a$out, "lme4_variance_components.csv"))
empty_check <- data.frame(layer = character(0), response = character(0), covariate = character(0),
                          estimate = numeric(0), se = numeric(0), p_lrt = numeric(0))
write_table(if (length(checks)) do.call(rbind, checks) else empty_check, file.path(a$out, "assumption_check.csv"))
