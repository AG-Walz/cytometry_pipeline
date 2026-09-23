#!/usr/bin/env Rscript
# Step 6b: the same two models written by hand in Stan (scripts/drift/stan/), fitted with cmdstanr.
# A software cross-check of 06_yardstick_brms.R; writes effects_stan.csv in the shared long format.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_yardstick_stan.R --out DIR --config FILE [--chains 4] [--iter 2000] [--cores 4] [--layer A|B|AB]"
a <- parse_args(list(out = "", config = "",
                     chains = 4, iter = 2000, cores = 4, layer = "AB"), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("cmdstanr")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
ref <- cfg$reference_instrument
stan_dir <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stan")
eff <- list()

# Posterior draws of the instrument contrast and of the day effects (columns named by `days`).
run <- function(model, data, days) {
  fit <- model$sample(data = data, chains = a$chains, parallel_chains = a$cores,
                      iter_warmup = a$iter / 2, iter_sampling = a$iter / 2, refresh = 0,
                      adapt_delta = 0.95, show_messages = FALSE)
  # strip the posterior::draws_matrix class: a subset of it keeps the class and
  # data.frame() then names the column after the parameter, breaking rbind
  c_d <- fit$draws("c_d", format = "matrix")
  c_d <- matrix(as.numeric(c_d), nrow = nrow(c_d), dimnames = list(NULL, days))
  list(b_instr = as.numeric(fit$draws("b_instr", format = "matrix")), c_d = c_d)
}

if (grepl("A", a$layer)) {
  mA <- cmdstanr::cmdstan_model(file.path(stan_dir, "layer_a_normal.stan"))
  fa <- read_layer_a(a$out, cfg, tt); fa <- fa[fa$is_neg, ]
  for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
    d <- fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]
    d$patient <- droplevels(d$patient); d$acq_day <- droplevels(d$acq_day)
    sc <- sd(d$value); other <- setdiff(levels(d$instrument), ref)
    cat("\n[A]", key, "\n")
    r <- run(mA, list(N = nrow(d), P = nlevels(d$patient), D = nlevels(d$acq_day),
                      patient = as.integer(d$patient), day = as.integer(d$acq_day),
                      instr = as.numeric(d$instrument != ref), y = as.numeric(scale(d$value))),
             days = levels(d$acq_day))
    eff[[length(eff) + 1L]] <- effects_long("stan", "A", key, r$b_instr * sc, r$c_d * sc,
                                            if (length(other)) other[1] else NA)
  }
}
if (grepl("B", a$layer)) {
  mB <- cmdstanr::cmdstan_model(file.path(stan_dir, "layer_b_betabinom.stan"))
  fb <- read_layer_b(a$out)
  for (resp in unique(fb$response)) {
    d <- fb[fb$response == resp, ]
    d$patient <- droplevels(d$patient); d$acq_day <- droplevels(d$acq_day); d$cocktail <- droplevels(d$cocktail)
    other <- setdiff(levels(d$instrument), ref)
    cat("\n[B]", resp, "\n")
    r <- run(mB, list(N = nrow(d), P = nlevels(d$patient), D = nlevels(d$acq_day), C = nlevels(d$cocktail),
                      patient = as.integer(d$patient), day = as.integer(d$acq_day), cocktail = as.integer(d$cocktail),
                      instr = as.numeric(d$instrument != ref), prior_vacc = as.numeric(d$prior_vacc == "2"),
                      week_c = d$week_c, n_parent = as.integer(d$parent_count), k = as.integer(d$count)),
             days = levels(d$acq_day))
    eff[[length(eff) + 1L]] <- effects_long("stan", "B", resp, r$b_instr, r$c_d,
                                            if (length(other)) other[1] else NA)
  }
}
for (ly in unique(vapply(eff, function(e) e$layer[1], character(1))))
  write_table(do.call(rbind, eff[vapply(eff, function(e) e$layer[1] == ly, logical(1))]),
              file.path(a$out, sprintf("effects_stan_%s.csv", ly)))
