#!/usr/bin/env Rscript
# Step 6b: the same models written by hand in Stan (scripts/drift/stan/), fitted with cmdstanr.
# A software cross-check of 06_yardstick_brms.R. The fixed effects (instrument dummies and the
# protected covariates from the config `design` block) enter through one design matrix X; the
# gaussian program serves layer A and BG intensity, the beta-binomial program layer BG counts
# (no stimulation term) and layer B. Writes effects_stan_<layer>.csv and covariates_stan.csv.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_yardstick_stan.R --out DIR --config FILE [--chains 4] [--iter 2000] [--cores 4] [--layer A,BG,B]"
a <- parse_args(list(out = "", config = "",
                     chains = 4, iter = 2000, cores = 4, layer = "A,BG,B"), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("cmdstanr")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
check_design(cfg, tt)
ref <- cfg$reference_instrument
layers <- if (a$layer == "AB") c("A", "B") else trimws(strsplit(a$layer, ",")[[1]])   # "AB" kept for old calls
stan_dir <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stan")
unit <- cfg$design$unit; day_term <- setdiff(unlist(cfg$design$technical), "instrument")[1]
eff <- list(); covs <- list()

# Fixed-effects design matrix of one response frame (no intercept) and the table mapping its
# columns back to term and level. Unused factor levels are dropped first so K is exact.
design_matrix <- function(d, class) {
  vars <- design_fixed_vars(cfg, class)
  for (v in vars) if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  ft <- fixed_terms_table(d, vars)
  X <- if (length(vars)) model.matrix(as.formula(paste("~", paste(vars, collapse = " + "))), d)[, -1, drop = FALSE]
       else matrix(0, nrow(d), 0)
  list(X = X[, ft$coef, drop = FALSE], ft = ft)
}
# Posterior draws: beta split into instrument and covariate columns, c_d named by day.
run <- function(model, data, ft, days) {
  fit <- model$sample(data = data, chains = a$chains, parallel_chains = a$cores,
                      iter_warmup = a$iter / 2, iter_sampling = a$iter / 2, refresh = 0,
                      adapt_delta = 0.95, show_messages = FALSE)
  # strip the posterior::draws_matrix class so the columns rbind as plain numerics
  as_plain <- function(x, nm) matrix(as.numeric(x), nrow = nrow(x), dimnames = list(NULL, nm))
  c_d <- as_plain(fit$draws("c_d", format = "matrix"), days)
  beta <- if (nrow(ft)) as_plain(fit$draws("beta", format = "matrix"), ft$coef) else matrix(0, nrow(c_d), 0)
  fx <- split_fixed_draws(beta, ft)
  # the beta-binomial program carries the population-level time slope outside beta
  if ("mu_week" %in% fit$metadata()$stan_variables && !is.null(cfg$design$time))
    fx$covariates$week_c <- as.numeric(fit$draws("mu_week", format = "matrix"))
  list(instrument = fx$instrument, covariates = fx$covariates, c_d = c_d,
       n_divergent = sum(fit$sampler_diagnostics(format = "matrix")[, "divergent__"]))
}
keep <- function(r, layer, response, scale = 1) {
  eff[[length(eff) + 1L]] <<- effects_long("stan", layer, response, r$instrument * scale, r$c_d * scale)
  covs[[length(covs) + 1L]] <<- covariates_long("stan", layer, response, lapply(r$covariates, function(m) m * scale))
  if (r$n_divergent > 0) cat("  ", r$n_divergent, "divergent transitions\n")
}
prep <- function(d) {   # drop unused unit and day levels so P and D match the frame
  d[[unit]] <- droplevels(d[[unit]]); d[[day_term]] <- droplevels(d[[day_term]]); d
}
fit_gaussian <- function(model, d, class) {
  dm <- design_matrix(d, class)
  run(model, list(N = nrow(d), P = nlevels(d[[unit]]), D = nlevels(d[[day_term]]), K = ncol(dm$X),
                  unit = as.integer(d[[unit]]), day = as.integer(d[[day_term]]), X = dm$X,
                  y = as.numeric(scale(d$value))),
      dm$ft, levels(d[[day_term]]))
}
fit_counts <- function(model, d, class, stim = NULL) {
  dm <- design_matrix(d, class)
  has_stim <- !is.null(stim) && stim %in% names(d)
  s <- if (has_stim) droplevels(d[[stim]]) else factor(rep(1L, nrow(d)))
  run(model, list(N = nrow(d), P = nlevels(d[[unit]]), D = nlevels(d[[day_term]]), C = nlevels(s), K = ncol(dm$X),
                  has_stim = as.integer(has_stim), unit = as.integer(d[[unit]]), day = as.integer(d[[day_term]]),
                  stim = as.integer(s), X = dm$X, week_c = if ("week_c" %in% names(d)) d$week_c else rep(0, nrow(d)),
                  n_parent = as.integer(d$parent_count), k = as.integer(d$count)),
      dm$ft, levels(d[[day_term]]))
}

if (any(c("A", "BG") %in% layers)) mN <- cmdstanr::cmdstan_model(file.path(stan_dir, "design_normal.stan"))
if (any(c("BG", "B") %in% layers)) mBB <- cmdstanr::cmdstan_model(file.path(stan_dir, "design_betabinom.stan"))

if ("A" %in% layers) {
  fa <- read_layer_a(a$out, cfg, tt, ref); fa <- fa[fa$is_neg, ]
  for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
    d <- prep(fa[paste(fa$marker, fa$statistic, sep = ":") == key, ])
    cat("\n[A]", key, "\n")
    keep(fit_gaussian(mN, d, "technical"), "A", key, scale = sd(d$value))
  }
}
if ("BG" %in% layers) {
  bc <- read_layer_bg_counts(a$out, cfg, ref)
  for (resp in unique(bc$response)) {
    d <- prep(bc[bc$response == resp, ])
    cat("\n[BG counts]", resp, "\n")
    keep(fit_counts(mBB, d, "background_counts"), "BG", resp)
  }
  bi <- read_layer_bg_intensity(a$out, cfg, tt, ref)
  for (resp in unique(bi$response)) {
    d <- prep(bi[bi$response == resp, ])
    cat("\n[BG intensity]", resp, "\n")
    keep(fit_gaussian(mN, d, "background_intensity"), "BG", resp, scale = sd(d$value))
  }
}
if ("B" %in% layers) {
  fb <- read_layer_b(a$out, ref, cfg)
  for (resp in unique(fb$response)) {
    d <- prep(fb[fb$response == resp, ])
    cat("\n[B]", resp, "\n")
    keep(fit_counts(mBB, d, "endpoint", stim = cfg$design$stimulation), "B", resp)
  }
}
for (ly in unique(vapply(eff, function(e) e$layer[1], character(1))))
  write_table(do.call(rbind, eff[vapply(eff, function(e) e$layer[1] == ly, logical(1))]),
              file.path(a$out, sprintf("effects_stan_%s.csv", ly)))
cv_path <- file.path(a$out, "covariates_stan.csv"); cv_new <- do.call(rbind, covs)
if (file.exists(cv_path)) { old <- read.csv(cv_path, stringsAsFactors = FALSE); old <- old[!old$layer %in% cv_new$layer, ]; cv_new <- rbind(old, cv_new) }
write_table(cv_new, cv_path)
