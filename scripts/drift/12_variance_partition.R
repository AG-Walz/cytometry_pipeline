#!/usr/bin/env Rscript
# Step 12: where do the differences between tubes come from, per data representation?
# (docs/plans/drift_three_questions_plan.md, parts C and D). For every readout (MFI of a marker in a
# population from step 11, cytokine gate in a parent from step 5) and every representation in
# PARTITION_REPS, one model with the unwanted factors (technical factors, acquisition day) and the
# wanted ones (patient, time, protected covariates, stimulated, pool) is fitted and its variance split
# into shares (src/harmonization/drift/drift_partition.R).
#   --engine lme4  point estimates in minutes, plus a refit without each non-instrument technical
#                  factor (variant "without_<factor>") to show how variance moves between the nested
#                  technical factors
#   --engine brms  posterior shares with 95 percent intervals; one file per fit under partition_brms/
#                  so an interrupted run resumes; --workers fits run in parallel, 4 chains each
# Writes variance_partition_<engine>.csv and, for brms, partition_fits_brms.csv (R-hat, divergences).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "12_variance_partition.R --out DIR --config FILE [--engine lme4|brms] [--reps all|name,...] [--workers 6] [--chains 4] [--iter 2000]"
a <- parse_args(list(out = "", config = "", engine = "lme4", reps = "all", workers = 6, chains = 4, iter = 2000), usage = usage)
cfg <- read_drift_config(a$config)
need_packages(if (a$engine == "brms") c("brms", "lme4") else "lme4")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
check_design(cfg, tt)
tt <- tt[!tt$filename %in% excluded_tubes(a$out), ]
rd <- function(f) { p <- file.path(a$out, f); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL }
mfi <- rd("mfi.csv"); fq <- rd("gated_frequencies.csv")
quality <- rd("sample_quality.csv")   # acquisition audit (step 11b); visit-level viability and live cells
if (is.null(quality)) cat("no sample_quality.csv (step 11b not run): no sample-quality term\n")
if (is.null(mfi)) cat("no mfi.csv (step 11 not run): MFI representations skipped\n")
reps <- if (a$reps == "all") names(PARTITION_REPS) else trimws(strsplit(a$reps, ",")[[1]])
ref <- cfg$reference_instrument

frames <- list()
for (r in reps) { d <- partition_frame(r, mfi, fq, tt, cfg, ref, quality); if (!is.null(d) && nrow(d)) frames[[r]] <- d }
tasks <- do.call(rbind, lapply(names(frames), function(r) data.frame(rep = r, readout = unique(frames[[r]]$readout), stringsAsFactors = FALSE)))
cat(nrow(tasks), "fits:", paste(sprintf("%s (%d)", names(table(tasks$rep)), table(tasks$rep)), collapse = ", "), "\n")

one_frame <- function(r, ro) {
  d <- frames[[r]]; d <- d[d$readout == ro, ]; d[] <- lapply(d, function(x) if (is.factor(x)) droplevels(x) else x)
  if ("neg_offset" %in% names(d)) d$neg_offset <- d$neg_offset - mean(d$neg_offset)   # centring moves only the intercept
  d
}

if (a$engine == "lme4") {
  fit_lme4 <- function(d, terms, r) {
    f <- partition_formula(terms, r, "lme4"); fam <- PARTITION_REPS[[r]]$family
    if (fam == "gaussian") { d$y <- as.numeric(scale(d$y)); lme4::lmer(f, data = d, REML = TRUE) }
    else lme4::glmer(f, data = d, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  }
  out <- list()
  for (i in seq_len(nrow(tasks))) {
    r <- tasks$rep[i]; ro <- tasks$readout[i]; d <- one_frame(r, ro)
    terms <- partition_terms(cfg, d, r)
    drops <- setdiff(intersect(names(terms$groups)[terms$groups == "unwanted"], design_technical_factors(cfg)), "instrument")
    for (v in c("full", paste0("without_", drops))) {
      tm <- terms
      if (v != "full") { x <- sub("^without_", "", v); tm$fixed <- setdiff(tm$fixed, x); tm$groups <- tm$groups[names(tm$groups) != x] }
      fit <- tryCatch(suppressMessages(suppressWarnings(fit_lme4(d, tm, r))), error = function(e) { cat("  ", r, ro, v, "failed:", conditionMessage(e), "\n"); NULL })
      if (is.null(fit)) next
      s <- summarise_shares(lme4_shares(fit, tm, r), scale_var = if (PARTITION_REPS[[r]]$family == "gaussian") var(d$y) else 1)
      s <- split_counting_noise(s, counting_noise_var(d, r, cfg))
      out[[length(out) + 1L]] <- cbind(readout = ro, representation = r, variant = v, engine = "lme4", n = nrow(d),
                                       singular = lme4::isSingular(fit), s, stringsAsFactors = FALSE)
    }
    if (i %% 10 == 0) cat(" ", i, "of", nrow(tasks), "\n")
  }
  res <- do.call(rbind, out)
  write_table(res, file.path(a$out, "variance_partition_lme4.csv"))
} else {
  options(brms.backend = if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan")
  dir_fits <- ensure_dir(file.path(a$out, "partition_brms"))
  priors <- list(
    gaussian = c(brms::prior(normal(0, 1), class = "b"), brms::prior(normal(0, 1), class = "sd"), brms::prior(normal(0, 1), class = "sigma")),
    binomial = c(brms::prior(normal(-5, 2.5), class = "Intercept"), brms::prior(normal(0, 1), class = "b"), brms::prior(normal(0, 1), class = "sd")),
    bernoulli = c(brms::prior(normal(-1, 2), class = "Intercept"), brms::prior(normal(0, 1), class = "b"), brms::prior(normal(0, 1), class = "sd")))
  run_task <- function(i) {
    r <- tasks$rep[i]; ro <- tasks$readout[i]
    stub <- file.path(dir_fits, paste0(r, "__", gsub("[^A-Za-z0-9]", "_", ro)))
    if (file.exists(paste0(stub, ".csv"))) return(invisible(stub))
    d <- one_frame(r, ro); terms <- partition_terms(cfg, d, r); fam <- PARTITION_REPS[[r]]$family
    sv <- if (fam == "gaussian") var(d$y) else 1
    if (fam == "gaussian") d$y <- as.numeric(scale(d$y))
    t0 <- Sys.time()
    fit <- tryCatch(brms::brm(partition_formula(terms, r, "brms"), data = d,
      family = switch(fam, gaussian = gaussian(), binomial = binomial(), bernoulli = brms::bernoulli()),   # counts: binomial + (1 | filename)
      prior = priors[[fam]], chains = a$chains, cores = a$chains, iter = a$iter, refresh = 0, silent = 2,
      control = list(adapt_delta = if (fam == "gaussian") 0.95 else 0.99, max_treedepth = 12)),
      error = function(e) { cat("  ", r, ro, "failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(fit)) return(invisible(NULL))
    s <- summarise_shares(brms_shares(fit, terms, r, d), scale_var = sv)
    s <- split_counting_noise(s, counting_noise_var(d, r, cfg))
    np <- brms::nuts_params(fit)
    write.csv(cbind(readout = ro, representation = r, variant = "full", engine = "brms", n = nrow(d), singular = NA, s),
              paste0(stub, ".csv"), row.names = FALSE)
    write.csv(data.frame(readout = ro, representation = r, n = nrow(d), max_rhat = max(brms::rhat(fit), na.rm = TRUE),
                         n_divergent = sum(np$Value[np$Parameter == "divergent__"]),
                         seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))),
              paste0(stub, ".diag"), row.names = FALSE)
    cat(sprintf("  done %s %s (%.0f s)\n", r, ro, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    invisible(stub)
  }
  invisible(parallel::mclapply(seq_len(nrow(tasks)), run_task, mc.cores = a$workers, mc.preschedule = FALSE))
  rd_all <- function(pat) { fs <- list.files(dir_fits, pat, full.names = TRUE); if (length(fs)) do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE)) }
  write_table(rd_all("\\.csv$"), file.path(a$out, "variance_partition_brms.csv"))
  write_table(rd_all("\\.diag$"), file.path(a$out, "partition_fits_brms.csv"))
}
