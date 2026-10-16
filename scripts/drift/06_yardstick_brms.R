#!/usr/bin/env Rscript
# Step 6a: Bayesian hierarchical yardstick with brms. Fits every layer A statistic (technical),
# every layer BG response (cytokine background in negs: gate counts and median intensity) and
# every layer B response (endpoint). Formulas come from design_formula() and the config `design`
# block. Writes posterior draws of the day and instrument effects in the shared long format
# (effects_brms_<layer>.csv), the protected covariates (covariates_brms.csv) and per-fit summaries.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_yardstick_brms.R --out DIR --config FILE [--chains 4] [--iter 2000] [--cores 4] [--layer A,BG,B]"
a <- parse_args(list(out = "", config = "",
                     chains = 4, iter = 2000, cores = 4, layer = "A,BG,B"), usage = usage)
cfg <- read_drift_config(a$config)
need_packages(c("brms"))
backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
options(mc.cores = a$cores, brms.backend = backend)
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
check_design(cfg, tt)
ref <- cfg$reference_instrument
layers <- if (a$layer == "AB") c("A", "B") else trimws(strsplit(a$layer, ",")[[1]])   # "AB" kept for old calls
unit <- cfg$design$unit; day_term <- setdiff(unlist(cfg$design$technical), "instrument")[1]
eff <- list(); covs <- list(); fits_summary <- list()
save_fit <- function(fit, layer, response)
  saveRDS(fit, file.path(ensure_dir(file.path(a$out, "fits_brms")), paste0(layer, "_", gsub("[^A-Za-z0-9]", "_", response), ".rds")))

gauss_priors <- c(brms::prior(normal(0, 1), class = "b"), brms::prior(normal(0, 1), class = "sd"),
                  brms::prior(normal(0, 1), class = "sigma"))
count_priors <- c(brms::prior(normal(-5, 2.5), class = "Intercept"), brms::prior(normal(0, 1), class = "b"),
                  brms::prior(normal(0, 1), class = "sd"))

# Instrument contrasts (every non-reference level) and day effects, plus the protected covariates,
# all multiplied by `scale` to undo a standardised response.
extract <- function(fit, layer, response, frame, class, scale = 1) {
  dr <- as.data.frame(fit)
  ft <- fixed_terms_table(frame, design_fixed_vars(cfg, class)); ft$coef <- paste0("b_", ft$coef)
  fx <- split_fixed_draws(dr[, intersect(ft$coef, names(dr)), drop = FALSE], ft)
  day_cols <- grep(sprintf("^r_%s\\[", day_term), names(dr), value = TRUE)
  days <- as.matrix(dr[, day_cols, drop = FALSE])
  colnames(days) <- sub(sprintf("^r_%s\\[(.*),Intercept\\]$", day_term), "\\1", day_cols)
  e <- effects_long("brms", layer, response, fx$instrument * scale, days * scale)
  cv <- covariates_long("brms", layer, response, lapply(fx$covariates, function(m) m * scale))
  list(effects = e, covariates = cv)
}
keep <- function(fit, layer, response, frame, class, scale = 1) {
  x <- extract(fit, layer, response, frame, class, scale)
  eff[[length(eff) + 1L]] <<- x$effects; covs[[length(covs) + 1L]] <<- x$covariates
  fits_summary[[length(fits_summary) + 1L]] <<- brms_fit_row(fit, layer, response, nrow(frame), scale, unit, day_term)
  save_fit(fit, layer, response)
}
# Gaussian fit on the standardised value so priors are comparable across statistics.
fit_gaussian <- function(d, class) {
  d$y <- as.numeric(scale(d$value))
  brms::brm(design_formula(cfg, class, "brms"), data = d, family = gaussian(), prior = gauss_priors,
            chains = a$chains, iter = a$iter, cores = a$cores, refresh = 0, silent = 2,
            control = list(adapt_delta = 0.95))
}
fit_counts <- function(d, class, iter) {
  brms::brm(design_formula(cfg, class, "brms"), data = d, family = brms::beta_binomial(), prior = count_priors,
            chains = a$chains, iter = iter, cores = a$cores, refresh = 0, silent = 2,
            control = list(adapt_delta = 0.99, max_treedepth = 12))   # beta-binomial with unit slopes needs it
}

if ("A" %in% layers) {
  fa <- read_layer_a(a$out, cfg, tt, ref); fa <- fa[fa$is_neg, ]
  for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
    d <- fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]
    cat("\n[A]", key, "n =", nrow(d), "\n")
    keep(fit_gaussian(d, "technical"), "A", key, d, "technical", scale = sd(d$value))
  }
}

if ("BG" %in% layers) {
  bc <- read_layer_bg_counts(a$out, cfg, ref)
  for (resp in unique(bc$response)) {
    d <- bc[bc$response == resp, ]
    cat("\n[BG counts]", resp, "n =", nrow(d), "\n")
    keep(fit_counts(d, "background_counts", a$iter), "BG", resp, d, "background_counts")
  }
  bi <- read_layer_bg_intensity(a$out, cfg, tt, ref)
  for (resp in unique(bi$response)) {
    d <- bi[bi$response == resp, ]
    cat("\n[BG intensity]", resp, "n =", nrow(d), "\n")
    keep(fit_gaussian(d, "background_intensity"), "BG", resp, d, "background_intensity", scale = sd(d$value))
  }
}

if ("B" %in% layers) {
  fb <- read_layer_b(a$out, ref, cfg)
  for (resp in unique(fb$response)) {
    d <- fb[fb$response == resp, ]
    cat("\n[B]", resp, "n =", nrow(d), "\n")
    keep(fit_counts(d, "endpoint", max(a$iter, 4000)), "B", resp, d, "endpoint")
  }
}

# one effects file per layer so a layer can be refitted alone without losing the others
for (ly in unique(vapply(eff, function(e) e$layer[1], character(1))))
  write_table(do.call(rbind, eff[vapply(eff, function(e) e$layer[1] == ly, logical(1))]),
              file.path(a$out, sprintf("effects_brms_%s.csv", ly)))
# summary tables hold all layers: rows of the layers fitted now replace the old ones
merge_layers <- function(new, path) {
  if (file.exists(path)) { old <- read.csv(path, stringsAsFactors = FALSE); old <- old[!old$layer %in% new$layer, ]; new <- rbind(old, new) }
  write_table(new, path)
}
merge_layers(do.call(rbind, fits_summary), file.path(a$out, "fits_brms_summary.csv"))
merge_layers(do.call(rbind, covs), file.path(a$out, "covariates_brms.csv"))
