#!/usr/bin/env Rscript
# Step 6b: the same two models written by hand in Stan (scripts/mtb_drift/stan/), fitted with cmdstanr.
# A software cross-check of 06_yardstick_brms.R; writes effects_stan.csv in the shared long format.
#   Rscript scripts/mtb_drift/06_yardstick_stan.R --out DIR [--chains 4] [--iter 2000]
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_yardstick_stan.R --out DIR [--config FILE] [--chains 4] [--iter 2000] [--cores 4] [--layer A|B|AB]"
a <- parse_args(list(out = "", config = file.path(REPO, "config", "mtb_drift.yaml"),
                     chains = 4, iter = 2000, cores = 4, layer = "AB"), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("cmdstanr")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
ref <- cfg$reference_instrument
stan_dir <- file.path(REPO, "scripts", "mtb_drift", "stan")
eff <- list()

run <- function(model, data) {
  fit <- model$sample(data = data, chains = a$chains, parallel_chains = a$cores,
                      iter_warmup = a$iter / 2, iter_sampling = a$iter / 2, refresh = 0,
                      adapt_delta = 0.95, show_messages = FALSE)
  # strip the posterior::draws_matrix class: a subset of it keeps the class and
  # data.frame() then names the column after the parameter, breaking rbind
  plain <- function(d) { m <- matrix(as.numeric(d), nrow = nrow(d)); colnames(m) <- colnames(d); m }
  list(b_instr = as.numeric(fit$draws("b_instr", format = "matrix")),
       c_d = plain(fit$draws("c_d", format = "matrix")), fit = fit)
}
day_matrix <- function(cd, levels) { colnames(cd) <- levels; cd }

if (grepl("A", a$layer)) {
  mA <- cmdstanr::cmdstan_model(file.path(stan_dir, "layer_a_normal.stan"))
  fa <- layer_a_frame(read.csv(file.path(a$out, "neg_summaries.csv")),
                      read.csv(file.path(a$out, "neg_emd.csv")), tt, cfg$technical_markers)
  fa <- fa[fa$is_neg, ]
  for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
    d <- fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]
    d$patient <- droplevels(d$patient); d$acq_day <- droplevels(d$acq_day)
    sc <- sd(d$value); other <- setdiff(levels(d$instrument), ref)
    cat("\n[A]", key, "\n")
    r <- run(mA, list(N = nrow(d), P = nlevels(d$patient), D = nlevels(d$acq_day),
                      patient = as.integer(d$patient), day = as.integer(d$acq_day),
                      instr = as.numeric(d$instrument != ref), y = as.numeric(scale(d$value))))
    eff[[length(eff) + 1L]] <- effects_long("stan", "A", key, r$b_instr * sc,
                                            day_matrix(r$c_d, levels(d$acq_day)) * sc,
                                            if (length(other)) other[1] else NA)
  }
}
if (grepl("B", a$layer)) {
  mB <- cmdstanr::cmdstan_model(file.path(stan_dir, "layer_b_betabinom.stan"))
  fb <- layer_b_frame(read.csv(file.path(a$out, "gated_frequencies.csv")))
  for (resp in unique(fb$response)) {
    d <- fb[fb$response == resp, ]
    d$patient <- droplevels(d$patient); d$acq_day <- droplevels(d$acq_day); d$cocktail <- droplevels(d$cocktail)
    other <- setdiff(levels(d$instrument), ref)
    cat("\n[B]", resp, "\n")
    r <- run(mB, list(N = nrow(d), P = nlevels(d$patient), D = nlevels(d$acq_day), C = nlevels(d$cocktail),
                      patient = as.integer(d$patient), day = as.integer(d$acq_day), cocktail = as.integer(d$cocktail),
                      instr = as.numeric(d$instrument != ref), prior_vacc = as.numeric(d$prior_vacc == "2"),
                      week_c = d$week_c, n_parent = as.integer(d$parent_count), k = as.integer(d$count)))
    eff[[length(eff) + 1L]] <- effects_long("stan", "B", resp, r$b_instr, day_matrix(r$c_d, levels(d$acq_day)),
                                            if (length(other)) other[1] else NA)
  }
}
for (ly in unique(vapply(eff, function(e) e$layer[1], character(1))))
  write_table(do.call(rbind, eff[vapply(eff, function(e) e$layer[1] == ly, logical(1))]),
              file.path(a$out, sprintf("effects_stan_%s.csv", ly)))
