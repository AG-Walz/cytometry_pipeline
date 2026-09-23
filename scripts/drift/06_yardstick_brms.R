#!/usr/bin/env Rscript
# Step 6a: Bayesian hierarchical yardstick with brms. Fits every layer A statistic and every
# layer B response; writes posterior draws of the day and instrument effects in the shared
# long format (effects_brms.csv) plus per-model summaries.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "06_yardstick_brms.R --out DIR --config FILE [--chains 4] [--iter 2000] [--cores 4] [--layer A|B|AB]"
a <- parse_args(list(out = "", config = "",
                     chains = 4, iter = 2000, cores = 4, layer = "AB"), usage = usage)
cfg <- read_drift_config(a$config)
need_packages(c("brms"))
backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
options(mc.cores = a$cores, brms.backend = backend)
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
ref <- cfg$reference_instrument
eff <- list(); fits_summary <- list()
save_fit <- function(fit, layer, response)
  saveRDS(fit, file.path(ensure_dir(file.path(a$out, "fits_brms")), paste0(layer, "_", gsub("[^A-Za-z0-9]", "_", response), ".rds")))

extract_effects <- function(fit, layer, response, frame) {
  dr <- as.data.frame(fit)
  other <- setdiff(levels(frame$instrument), ref)
  instr_col <- grep(paste0("^b_instrument", other), names(dr), value = TRUE)
  inst <- if (length(instr_col)) dr[[instr_col[1]]] else rep(0, nrow(dr))
  day_cols <- grep("^r_acq_day\\[", names(dr), value = TRUE)
  days <- as.matrix(dr[, day_cols, drop = FALSE])
  colnames(days) <- sub("^r_acq_day\\[(.*),Intercept\\]$", "\\1", day_cols)
  effects_long("brms", layer, response, inst, days, if (length(other)) other[1] else NA)
}

if (grepl("A", a$layer)) {
  fa <- read_layer_a(a$out, cfg, tt, ref); fa <- fa[fa$is_neg, ]
  for (key in unique(paste(fa$marker, fa$statistic, sep = ":"))) {
    d <- fa[paste(fa$marker, fa$statistic, sep = ":") == key, ]
    d$y <- as.numeric(scale(d$value))          # standardised so priors are comparable across statistics
    sc <- sd(d$value)
    cat("\n[A]", key, "n =", nrow(d), "\n")
    fit <- brms::brm(y ~ 1 + instrument + (1 | patient) + (1 | acq_day), data = d, family = gaussian(),
                     prior = c(brms::prior(normal(0, 1), class = "b"), brms::prior(normal(0, 1), class = "sd"),
                               brms::prior(normal(0, 1), class = "sigma")),
                     chains = a$chains, iter = a$iter, cores = a$cores, refresh = 0, silent = 2,
                     control = list(adapt_delta = 0.95))
    e <- extract_effects(fit, "A", key, d); e$value <- e$value * sc   # back to the data scale
    eff[[length(eff) + 1L]] <- e
    fits_summary[[length(fits_summary) + 1L]] <- brms_fit_row(fit, "A", key, nrow(d), scale = sc)
    save_fit(fit, "A", key)
  }
}

if (grepl("B", a$layer)) {
  fb <- read_layer_b(a$out, ref)
  for (resp in unique(fb$response)) {
    d <- fb[fb$response == resp, ]
    cat("\n[B]", resp, "n =", nrow(d), "\n")
    fit <- brms::brm(count | trials(parent_count) ~ 1 + instrument + prior_vacc + (1 + week_c | patient) +
                       (1 | acq_day) + (1 | cocktail), data = d, family = brms::beta_binomial(),
                     prior = c(brms::prior(normal(-5, 2.5), class = "Intercept"), brms::prior(normal(0, 1), class = "b"),
                               brms::prior(normal(0, 1), class = "sd")),
                     chains = a$chains, iter = max(a$iter, 4000), cores = a$cores, refresh = 0, silent = 2,
                     control = list(adapt_delta = 0.99, max_treedepth = 12))   # beta-binomial with patient slopes needs it
    eff[[length(eff) + 1L]] <- extract_effects(fit, "B", resp, d)
    fits_summary[[length(fits_summary) + 1L]] <- brms_fit_row(fit, "B", resp, nrow(d))
    save_fit(fit, "B", resp)
  }
}
# one effects file per layer so a layer can be refitted alone without losing the other
for (ly in unique(vapply(eff, function(e) e$layer[1], character(1))))
  write_table(do.call(rbind, eff[vapply(eff, function(e) e$layer[1] == ly, logical(1))]),
              file.path(a$out, sprintf("effects_brms_%s.csv", ly)))
fs_path <- file.path(a$out, "fits_brms_summary.csv"); fs_new <- do.call(rbind, fits_summary)
if (file.exists(fs_path)) { old <- read.csv(fs_path); old <- old[!old$layer %in% fs_new$layer, ]; fs_new <- rbind(old, fs_new) }
write_table(fs_new, fs_path)
