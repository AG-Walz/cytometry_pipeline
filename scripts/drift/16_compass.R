#!/usr/bin/env Rscript
# Step 16: Boolean cytokine subsets, COMPASS and the two-sided subset model
# (docs/plans/compass_plan_2026-10-14.md). Three stages, run in order:
#   --stage counts   per tube and parent, the 2^M Boolean subsets of the cytokine gates from the FlowJo
#                    hand gates (per-cell gate membership). Writes boolean_counts.csv. Stops unless the
#                    subsets sum to the parent and every marker's marginal equals its count in
#                    gated_frequencies.csv.
#   --stage compass  COMPASS (Lin et al. 2015) per parent, one individual per stimulated tube against its
#                    matched neg, one subset list shared by both parents, several seeds. The response
#                    map and the scores come from the mean model (mean_gamma averaged over the seeds),
#                    the intervals from the seeds' draws pooled. --variant no_low_info refits without
#                    tubes whose own or matched-neg parent is below low_info_parent; --reuse_fits
#                    recomputes the tables from saved fits instead of sampling again. Writes
#                    compass_gamma<v>.csv, compass_scores<v>.csv, compass_mcmc<v>.csv, compass_fits/, figures.
#   --stage model    two-sided subset model (brms): stimulated and neg counts of every kept subset of both
#                    parents in one binomial model, the stimulation effect varying by subset and time.
#                    --time visit|weeks, --variant full|no_low_info|no_fallback_neg,
#                    --synthetic none|effect|null (simulated counts on the real design), --patient_stim
#                    (sensitivity: each patient also gets its own stimulation effect per subset). Writes
#                    compass_model/<name>.rds, subset_effects_<name>.csv, hypotheses_<name>.csv.
#   --stage contrasts exploratory late-vs-reference contrasts from a saved fit (same --time, --variant).
#   --stage summary  descriptive tables (compass_pools.csv per patient x pool, compass_by_visit.csv) and figures.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- paste("16_compass.R --stage counts|compass|model|contrasts|summary --out DIR --config FILE",
               "[--wsp_dir DIR --fcs_dir DIR] [--time visit|weeks] [--variant full|no_low_info|no_fallback_neg]",
               "[--synthetic none|effect|null] [--patient_stim] [--reuse_fits] [--cores 4]")
a <- parse_args(list(stage = "", out = "", config = "", wsp_dir = "none", fcs_dir = "none", time = "visit",
                     variant = "full", synthetic = "none", patient_stim = FALSE, reuse_fits = FALSE, cores = 4), usage = usage)
# outputs of a variant carry its name; the full run keeps the plain file names
sfx <- if (a$variant == "full") "" else paste0("_", a$variant)
cfg <- read_drift_config(a$config)
cc <- compass_cfg(cfg)
mk <- unlist(cc$markers); nul <- null_subset(mk)
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------------- counts
if (a$stage == "counts") {
  need_packages(c("CytoML", "flowWorkspace"))
  if (a$wsp_dir == "none" || a$fcs_dir == "none") stop("--stage counts needs --wsp_dir and --fcs_dir", call. = FALSE)
  sigs <- read.csv(file.path(a$out, "panel_signatures.csv"), stringsAsFactors = FALSE)
  excl <- excluded_tubes(a$out)
  rows <- list(); seen <- character(0)
  for (w in list_workspaces(a$wsp_dir)) {
    for (gs in open_workspace(w, a$fcs_dir, sigs)) for (s in flowWorkspace::sampleNames(gs)) {
      gh <- gs[[s]]; fn <- gh_filename(gh)
      if (!fn %in% tt$filename || fn %in% excl || fn %in% seen) next
      seen <- c(seen, fn)
      for (pp in cc$parents) {
        node <- find_population(gh, cfg$parent_populations[[pp]])
        if (is.na(node)) stop("parent ", pp, " not found in ", fn, call. = FALSE)
        kids <- flowWorkspace::gh_pop_get_children(gh, node, path = "full")
        leaf <- sub("^.*/", "", kids)
        kids <- kids[grepl(cfg$cytokine_gate_regex, leaf, ignore.case = TRUE)]
        canon <- canon_cytokine_gate(sub("^.*/", "", kids), cfg)
        if (!setequal(intersect(canon, mk), mk) || anyDuplicated(canon[canon %in% mk]))
          stop(fn, " ", pp, ": cytokine gates ", paste(sub("^.*/", "", kids), collapse = ", "),
               " do not map one-to-one onto ", paste(mk, collapse = ", "), call. = FALSE)
        par <- flowWorkspace::gh_pop_get_indices(gh, node)
        hits <- lapply(mk, function(m) flowWorkspace::gh_pop_get_indices(gh, kids[match(m, canon)]))
        names(hits) <- mk
        bc <- boolean_counts(par, hits)
        rows[[length(rows) + 1L]] <- data.frame(filename = fn, parent = pp, subset = names(bc), degree = subset_degree(names(bc)),
                                                count = as.integer(bc), parent_count = sum(par), stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  # hard checks: subsets sum to the parent; marginals equal the single-gate counts of step 5
  sums <- aggregate(count ~ filename + parent + parent_count, out, sum)
  if (any(sums$count != sums$parent_count)) stop("Boolean subsets do not sum to the parent in ", sum(sums$count != sums$parent_count), " tube-parents", call. = FALSE)
  fq <- read.csv(file.path(a$out, "gated_frequencies.csv"), stringsAsFactors = FALSE)
  marg <- do.call(rbind, lapply(split(out, paste(out$filename, out$parent)), function(d) {
    v <- boolean_marginals(setNames(d$count, d$subset), mk)
    data.frame(filename = d$filename[1], parent = d$parent[1], gate = mk, boolean = unname(v), parent_boolean = d$parent_count[1])
  }))
  cmp <- merge(marg, fq[, c("filename", "parent", "gate", "count", "parent_count")], by = c("filename", "parent", "gate"), all = TRUE)
  bad <- cmp[is.na(cmp$boolean) | is.na(cmp$count) | cmp$boolean != cmp$count | cmp$parent_boolean != cmp$parent_count, ]
  write_table(cmp, file.path(a$out, "boolean_marginal_check.csv"))
  if (nrow(bad)) { print(utils::head(bad, 20)); stop(nrow(bad), " marker marginals differ from gated_frequencies.csv", call. = FALSE) }
  cat("check passed: subsets sum to the parent and all", nrow(cmp), "marker marginals equal gated_frequencies.csv\n")
  write_table(out, file.path(a$out, "boolean_counts.csv"))
  quit(status = 0)
}

bc <- read.csv(file.path(a$out, "boolean_counts.csv"), stringsAsFactors = FALSE)
fq <- read.csv(file.path(a$out, "gated_frequencies.csv"), stringsAsFactors = FALSE)
subsets_all <- compass_subset_names(mk)

#' Wide count matrix (filenames x all subsets) of one parent.
wide_counts <- function(pp, files) {
  d <- bc[bc$parent == pp & bc$filename %in% files, ]
  m <- matrix(0, length(files), length(subsets_all), dimnames = list(files, subsets_all))
  m[cbind(match(d$filename, files), match(d$subset, subsets_all))] <- d$count
  m
}

# the matched neg of every stimulated tube (step 5 policy; pooled fallbacks are comma-joined)
pairs <- unique(fq[!fq$is_neg & !is.na(fq$neg_filename), c("filename", "parent", "neg_filename", "neg_match")])
stim_files <- sort(unique(pairs$filename))

# subsets kept: the COMPASS filter on the stimulated counts of each parent, united over parents
keep <- sort(unique(unlist(lapply(cc$parents, function(pp)
  compass_filter(wide_counts(pp, stim_files), nul, cc$category_filter$min_count, cc$category_filter$min_tubes)))))
keep <- subsets_all[subsets_all %in% keep]   # COMPASS order

# ---------------------------------------------------------------------------------------- compass
if (a$stage == "compass") {
  need_packages("COMPASS")
  fdir <- ensure_dir(file.path(a$out, "compass_fits"))
  cat("kept subsets (", length(keep), "of", length(subsets_all) - 1, "non-null):\n ", paste(keep, collapse = "\n  "), "\n")
  gam_rows <- list(); sc_rows <- list(); mc_rows <- list()
  for (pp in cc$parents) {
    pr <- pairs[pairs$parent == pp, ]; pr <- pr[order(pr$filename), ]
    n_s_all <- wide_counts(pp, pr$filename)
    negs <- strsplit(pr$neg_filename, ",")
    n_u_all <- t(vapply(negs, function(f) colSums(wide_counts(pp, f)), numeric(length(subsets_all))))
    dimnames(n_u_all) <- dimnames(n_s_all)
    if (a$variant == "no_low_info") {
      ok <- pmin(rowSums(n_s_all), rowSums(n_u_all)) >= cc$low_info_parent
      cat(pp, ": dropping", sum(!ok), "of", length(ok), "tubes below", cc$low_info_parent, "parent events\n")
      pr <- pr[ok, ]; n_s_all <- n_s_all[ok, , drop = FALSE]; n_u_all <- n_u_all[ok, , drop = FALSE]
    }
    n_s <- collapse_to_kept(n_s_all, keep, nul, rowSums(n_s_all))
    n_u <- collapse_to_kept(n_u_all, keep, nul, rowSums(n_u_all))
    meta <- data.frame(filename = pr$filename, stringsAsFactors = FALSE)
    fits <- list()
    for (sd in cc$seeds) {
      ff <- file.path(fdir, sprintf("compass_%s%s_seed%d.rds", pp, sfx, sd))
      if (a$reuse_fits) {
        if (!file.exists(ff)) stop("--reuse_fits: no saved fit at ", ff, call. = FALSE)
        cat("reusing", basename(ff), "\n")
        f <- readRDS(ff)
        if (!identical(dim(f$fit$mean_gamma), c(nrow(n_s), ncol(n_s)))) stop("saved fit ", basename(ff), " does not match the current data", call. = FALSE)
      } else {
        cat("\nCOMPASS", pp, if (nzchar(sfx)) a$variant, "seed", sd, "\n")
        f <- COMPASS::SimpleCOMPASS(n_s, n_u, meta, "filename", iterations = cc$iterations, replications = cc$replications,
                                    verbose = FALSE, seed = sd)
        saveRDS(f, ff)
      }
      fits[[as.character(sd)]] <- f
    }
    f1 <- fits[[1]]
    deg <- subset_degree(colnames(n_s))
    # mean model: mean_gamma averaged over the seeds, and the seeds' draws pooled for the intervals.
    # The response map and the scores below come from it, not from one chain.
    mg <- mean_gamma_across_seeds(fits); dimnames(mg) <- list(pr$filename, colnames(n_s))
    fmean <- f1; fmean$fit$mean_gamma <- mg
    sc <- compass_score_draws(pooled_gamma_draws(fits), deg, length(mk))
    fs_seed <- sapply(fits, function(f) COMPASS::FunctionalityScore(f))
    pfs_seed <- sapply(fits, function(f) COMPASS::PolyfunctionalityScore(f))
    mg1 <- f1$fit$mean_gamma                       # seed 1 alone, kept to report what the mean model changed
    gam_rows[[pp]] <- data.frame(filename = rep(pr$filename, ncol(mg) - 1), parent = pp,
                                 subset = rep(colnames(mg)[-ncol(mg)], each = nrow(mg)), degree = rep(deg[-ncol(mg)], each = nrow(mg)),
                                 mean_gamma = as.vector(mg[, -ncol(mg)]), mean_gamma_seed1 = as.vector(mg1[, -ncol(mg1)]),
                                 stringsAsFactors = FALSE)
    small <- pmin(rowSums(n_s_all), rowSums(n_u_all)) < cc$low_info_parent
    sc_rows[[pp]] <- data.frame(filename = pr$filename, parent = pp, neg_filename = pr$neg_filename, neg_match = pr$neg_match,
                                n_stim_sharing_neg = neg_sharing(pr$neg_filename),
                                parent_count = rowSums(n_s_all), neg_parent_count = rowSums(n_u_all), low_info = small,
                                draw_summary(sc$FS, "FS"), draw_summary(sc$PFS, "PFS"),
                                FS_meangamma = COMPASS::FunctionalityScore(fmean), PFS_meangamma = COMPASS::PolyfunctionalityScore(fmean),
                                n_responding = rowSums(mg[, -ncol(mg)] > 0.5),
                                any_IFNg = apply(mg[, grepl("(^|&)IFNg", colnames(mg)) & colnames(mg) != nul, drop = FALSE] > 0.5, 1, any),
                                stringsAsFactors = FALSE)
    diag <- tryCatch(COMPASS::COMPASSMCMCDiagnosis(fits), error = function(e) NULL)
    # Gelman-Rubin R-hat of the Dirichlet hyperparameters alpha_s and alpha_u across the seeds
    rhat <- if (is.null(diag)) NULL else rbind(diag$diag$alpha_s$psrf, diag$diag$alpha_u$psrf)
    cors <- utils::combn(seq_along(fits), 2, function(ix) c(stats::cor(fs_seed[, ix[1]], fs_seed[, ix[2]]),
                                                             stats::cor(pfs_seed[, ix[1]], pfs_seed[, ix[2]])))
    mc_rows[[pp]] <- data.frame(parent = pp, variant = a$variant, n_tubes = nrow(n_s), n_subsets = ncol(n_s) - 1,
                                seeds = paste(cc$seeds, collapse = " "),
                                min_cor_FS = min(cors[1, ]), min_cor_PFS = min(cors[2, ]),
                                max_rhat_alpha = if (is.null(rhat)) NA_real_ else max(rhat[, 1], na.rm = TRUE),
                                max_abs_meangamma_mean_vs_seed1 = max(abs(mg - mg1)),
                                n_calls_flipped_mean_vs_seed1 = sum((mg[, -ncol(mg)] > 0.5) != (mg1[, -ncol(mg1)] > 0.5)),
                                n_calls = length(mg[, -ncol(mg)]), stringsAsFactors = FALSE)
  }
  gam <- do.call(rbind, gam_rows); sc <- do.call(rbind, sc_rows); mc <- do.call(rbind, mc_rows)
  sc <- merge(sc, tt[, intersect(c("filename", "patient", "week", "cocktail", "prior_vacc", "instrument", "operator", "acq_day", cc$visit_column), names(tt))], by = "filename")
  sc$dose <- visit_dose(as.character(sc[[cc$visit_column]]), cc$dose_regex)
  write_table(gam, file.path(a$out, paste0("compass_gamma", sfx, ".csv")))
  write_table(sc, file.path(a$out, paste0("compass_scores", sfx, ".csv")))
  write_table(mc, file.path(a$out, paste0("compass_mcmc", sfx, ".csv")))
  print(mc, row.names = FALSE)

  # agreement with step 13 MIMOSA IFNg calls: "any IFNg-containing subset responds" vs MIMOSA per parent
  cf <- file.path(a$out, "calls_compare.csv")
  if (file.exists(cf)) {
    calls <- read.csv(cf, stringsAsFactors = FALSE)
    ag <- do.call(rbind, lapply(cc$parents, function(pp) {
      m <- merge(sc[sc$parent == pp, c("filename", "any_IFNg")], calls[calls$response == paste0(pp, ":IFNg"), c("filename", "mimosa", "lab")], by = "filename")
      m$compass <- m$any_IFNg
      out <- calls_agreement(data.frame(response = paste0(pp, ":IFNg"), lab = m$lab, mimosa = m$mimosa, fisher = m$compass))
      names(out) <- sub("fisher", "compass", names(out)); out
    }))
    write_table(ag, file.path(a$out, paste0("compass_calls_agreement", sfx, ".csv")))
    print(ag, row.names = FALSE, digits = 2)
  }

  quit(status = 0)
}

# ---------------------------------------------------------------------------------------- model
if (a$stage == "model") {
  need_packages("brms")
  name <- paste0(paste(a$time, a$variant, sep = "_"), if (a$synthetic != "none") paste0("_syn_", a$synthetic), if (a$patient_stim) "_pstim")
  mdir <- ensure_dir(file.path(a$out, "compass_model"))
  d <- subset_model_frame(bc, tt, keep, cc)
  d <- d[d$parent %in% cc$parents, ]
  if (a$variant == "no_low_info") {
    # drop per parent: a tube whose own or matched-neg parent is small, and small neg tubes
    lowp <- unique(rbind(bc[bc$parent_count < cc$low_info_parent, c("filename", "parent")],
                         do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
                           f <- strsplit(pairs$neg_filename[i], ",")[[1]]
                           if (any(bc$parent_count[bc$filename %in% f & bc$parent == pairs$parent[i]] < cc$low_info_parent))
                             data.frame(filename = pairs$filename[i], parent = pairs$parent[i]) else NULL }))))
    d <- d[!paste(d$filename, d$parent) %in% paste(lowp$filename, lowp$parent), ]
  }
  if (a$variant == "no_fallback_neg") d <- d[!d$filename %in% pairs$filename[pairs$neg_match != "exact"], ]
  if (a$synthetic != "none") {
    # baseline per subset from the observed neg tubes; the generating effects are written next to the fit
    neg <- d[d$stim == 0, ]
    base <- tapply(seq_len(nrow(neg)), neg$ps, function(i) qlogis((sum(neg$count[i]) + 0.5) / (sum(neg$parent_count[i]) + 1)))
    poly <- paste0("CD4:", paste(c("IFNg", "IL2", "TNFa", "!CD107a"), collapse = "&"))
    if (!poly %in% names(base)) stop("synthetic target subset ", poly, " is not among the kept subsets", call. = FALSE)
    ramp <- c(V0 = 0, V1 = 1, V2 = 1.5, V3 = 2.5, V4 = 3)
    effect <- switch(a$synthetic,
      effect = function(ps, v) ifelse(ps == poly, ramp[v], ifelse(ps == "CD4:IFNg&!IL2&!TNFa&!CD107a", ifelse(v == "V0", 0, 1), 0)),
      null = function(ps, v) ifelse(grepl("^CD[48]:IFNg&!IL2&!TNFa&!CD107a$", ps) & v != "V0", 1, 0))
    d <- simulate_subset_counts(d, base, effect, sd_visit = 0.4, sd_day = if (a$synthetic == "null") 0.5 else 0.2,
                                sd_pool = 0.3, sd_obs = 0.3, seed = cc$seed)
    truth <- unique(data.frame(ps = d$ps, visit = as.character(d$visit), stringsAsFactors = FALSE))
    truth$true_log_or <- effect(truth$ps, truth$visit)
    write_table(truth, file.path(mdir, paste0("truth_", name, ".csv")))
  }
  d <- design_frame(d, modifyList(cfg, list(design = modifyList(cfg$design, list(time = NULL, time_df = 1)))),
                    ref = cfg$reference_instrument)
  d$prior_vacc <- factor(d$prior_vacc)
  tech <- intersect(design_technical_factors(cfg), names(d))
  tech <- tech[vapply(tech, function(v) nlevels(droplevels(factor(d[[v]]))) > 1, logical(1))]
  if (a$time == "visit") {
    tslope <- "stim:visit"; tvars <- "visit"
  } else {
    wb <- week_basis(d$week, cc$weeks_df); d <- cbind(d, wb)
    wn <- names(wb); tslope <- paste0("stim:", wn, collapse = " + "); tvars <- paste0("(", paste(wn, collapse = " + "), ")")
  }
  d$parent <- factor(d$parent)
  # fixed effects per parent, so the subsets of CD4 and of CD8 are each pooled towards their own
  # mean (one shared mean pulled null CD8 subsets towards the strong CD4 effects in the smoke test);
  # each subset deviates from its parent's stimulation effect and time course (random slopes by ps)
  rhs <- c(sprintf("parent * stim * %s", tvars), "parent * stim * prior_vacc", tech,
           sprintf("(1 + stim + %s | ps)", tslope), if (a$patient_stim) "(1 + stim | ps:patient)" else "(1 | ps:patient)",
           "(1 | ps:visit_id)", "(0 + stim | ps:cocktail)",
           sprintf("(1 | %s)", design_day_term(cfg)))
  # beta-binomial: overdispersion per subset through phi, in place of an observation-level random
  # effect (9,000 latent parameters that made sampling about 10 times slower in the smoke test)
  fml <- brms::bf(stats::as.formula(paste("count | trials(parent_count) ~", paste(rhs, collapse = " + "))), phi ~ 1 + (1 | ps))
  cat("model", name, ":", nrow(d), "rows,", length(unique(d$ps)), "subsets\n"); print(fml)
  pri <- c(brms::prior(normal(0, 1.5), class = "b"), brms::prior(normal(0, 1), class = "sd"),
           brms::prior(lkj(2), class = "cor"), brms::prior(student_t(3, -6, 2.5), class = "Intercept"),
           brms::prior(normal(4, 2), class = "Intercept", dpar = "phi"), brms::prior(normal(0, 1), class = "sd", dpar = "phi"))
  backend <- if (requireNamespace("cmdstanr", quietly = TRUE) && !is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL))) "cmdstanr" else "rstan"
  t0 <- Sys.time()
  fit <- brms::brm(fml, data = d, family = brms::beta_binomial(), prior = pri, chains = cc$chains, iter = cc$iter, warmup = cc$warmup,
                   cores = min(a$cores, cc$chains), threads = if (backend == "cmdstanr") brms::threading(cc$threads) else NULL,
                   backend = backend, seed = cc$seed, control = list(adapt_delta = cc$adapt_delta), refresh = 100, silent = 1)
  cat("sampling took", format(Sys.time() - t0), "\n")
  saveRDS(list(fit = fit, data = d, keep = keep, basis = if (a$time == "weeks") attr(wb, "basis")), file.path(mdir, paste0(name, ".rds")))

  # population-level predictions per subset at each time, typical sample (all random effects but the
  # subset's own set to zero), technical factors at their reference levels, prior_vacc at its first level
  if (a$time == "visit") {
    times <- data.frame(time = levels(d$visit), visit = factor(levels(d$visit), levels = levels(d$visit)))
  } else {
    wk <- c(V0 = 0, tapply(d$week[d$stim == 1], d$visit[d$stim == 1], stats::median))
    wk <- wk[!is.na(wk)]; wk <- wk[!duplicated(names(wk))]
    times <- data.frame(time = names(wk), week = unname(wk)); times <- cbind(times, week_basis_at(attr(wb, "basis"), times$week))
    late_wk <- stats::median(d$week[d$stim == 1 & d$dose %in% cc$late_doses])
    times <- rbind(times, cbind(data.frame(time = "late", week = late_wk), week_basis_at(attr(wb, "basis"), late_wk)))
  }
  pss <- sort(unique(d$ps))
  nd <- merge(merge(data.frame(ps = pss, stringsAsFactors = FALSE), data.frame(stim = 0:1)), times)
  nd$parent <- factor(sub(":.*$", "", nd$ps), levels = levels(d$parent))
  nd$parent_count <- 1L
  nd$prior_vacc <- factor(levels(d$prior_vacc)[1], levels = levels(d$prior_vacc))
  for (v in tech) nd[[v]] <- factor(levels(d[[v]])[1], levels = levels(d[[v]]))
  pfun <- function(x) brms::posterior_linpred(fit, newdata = x, re_formula = stats::as.formula(sprintf("~ (1 + stim + %s | ps)", tslope)),
                                              transform = TRUE, allow_new_levels = FALSE)
  P <- pfun(nd)   # draws x rows: probability of a cell in the subset
  # subset-level effects: log odds ratio stimulated vs neg, per subset and time
  eff <- do.call(rbind, lapply(split(seq_len(nrow(nd)), paste(nd$ps, nd$time)), function(ix) {
    i1 <- ix[nd$stim[ix] == 1]; i0 <- ix[nd$stim[ix] == 0]
    lor <- qlogis(P[, i1]) - qlogis(P[, i0]); dif <- P[, i1] - P[, i0]
    s <- post_summary(lor); s2 <- post_summary(dif)
    data.frame(ps = nd$ps[i1], time = nd$time[i1], log_or = s$mean, log_or_lo = s$lo, log_or_hi = s$hi, p_up = s$p_gt0,
               diff_pct = 100 * s2$mean, diff_pct_lo = 100 * s2$lo, diff_pct_hi = 100 * s2$hi,
               neg_pct = 100 * mean(P[, i0]), stim_pct = 100 * mean(P[, i1]), stringsAsFactors = FALSE)
  }))
  eff$parent <- sub(":.*$", "", eff$ps); eff$subset <- sub("^[^:]*:", "", eff$ps); eff$degree <- subset_degree(eff$subset)
  eff$call <- ifelse(eff$p_up >= 0.95, "up", ifelse(eff$p_up <= 0.05, "down", "unclear"))
  write_table(eff, file.path(a$out, paste0("subset_effects_", name, ".csv")))

  # hypotheses per parent: H1 response at late times, H2 polyfunctional share, H3 subsets up
  late <- if (a$time == "visit") paste0("V", cc$late_doses) else "late"
  late <- intersect(late, times$time); ref <- if (a$time == "visit") paste0("V", cc$reference_dose) else "V1"
  hd <- list()
  for (pp in cc$parents) for (tm in times$time) {
    cols <- nd$ps %in% pss[startsWith(pss, paste0(pp, ":"))] & nd$time == tm
    i1 <- which(cols & nd$stim == 1); i0 <- which(cols & nd$stim == 0)
    i0 <- i0[match(nd$ps[i1], nd$ps[i0])]
    hd[[paste(pp, tm)]] <- cbind(parent = pp, time = tm, draw = seq_len(nrow(P)),
                                 hypothesis_draws(P[, i1, drop = FALSE], P[, i0, drop = FALSE], subset_degree(sub("^[^:]*:", "", nd$ps[i1]))))
  }
  hd <- do.call(rbind, hd)
  late_mean <- function(pp, col) { x <- sapply(late, function(tm) hd[hd$parent == pp & hd$time == tm, col]); if (is.matrix(x)) rowMeans(x) else mean(x) }
  at <- function(pp, tm, col) hd[hd$parent == pp & hd$time == tm, col]
  H <- list()
  for (pp in cc$parents) {
    for (tm in times$time) for (col in c("response", "share_poly", "n_up", "n_down"))
      H[[length(H) + 1L]] <- cbind(data.frame(hypothesis = "level", parent = pp, time = tm, quantity = col), post_summary(at(pp, tm, col)))
    H[[length(H) + 1L]] <- cbind(data.frame(hypothesis = "H1", parent = pp, time = paste(late, collapse = "+"), quantity = "response_late"),
                                 post_summary(late_mean(pp, "response")))
    H[[length(H) + 1L]] <- cbind(data.frame(hypothesis = "H2", parent = pp, time = paste0(paste(late, collapse = "+"), " - ", ref), quantity = "share_poly_change"),
                                 post_summary(late_mean(pp, "share_poly") - at(pp, ref, "share_poly")))
    H[[length(H) + 1L]] <- cbind(data.frame(hypothesis = "H3", parent = pp, time = paste0(paste(late, collapse = "+"), " - ", ref), quantity = "n_up_change"),
                                 post_summary(late_mean(pp, "n_up") - at(pp, ref, "n_up")))
  }
  both <- late_mean(cc$parents[1], "response") > 0 & late_mean(cc$parents[2], "response") > 0
  H[[length(H) + 1L]] <- data.frame(hypothesis = "H1", parent = paste(cc$parents, collapse = "&"), time = paste(late, collapse = "+"),
                                    quantity = "P(both respond late)", mean = mean(both), lo = NA, hi = NA, p_gt0 = mean(both))
  H <- do.call(rbind, H); H$model <- name
  write_table(H, file.path(a$out, paste0("hypotheses_", name, ".csv")))
  # diagnostics
  rh <- brms::rhat(fit); ess <- brms::neff_ratio(fit)
  diag <- data.frame(model = name, rows = nrow(d), subsets = length(pss), max_rhat = max(rh, na.rm = TRUE),
                     min_neff_ratio = min(ess, na.rm = TRUE), divergent = n_divergent(fit),
                     minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  write_table(diag, file.path(mdir, paste0("diagnostics_", name, ".csv")))
  print(diag, row.names = FALSE)
  print(H[H$hypothesis != "level", ], row.names = FALSE, digits = 3)
  quit(status = 0)
}
# ---------------------------------------------------------------------------------------- contrasts
# Exploratory, chosen after H1 to H3 were seen (H3 proved uninformative on synthetic data): late minus
# reference-visit change of the net response, of the mean log odds ratio of the subsets with >= 2 markers
# and with 1 marker, of their difference, and of every subset, from a saved fit (--time, --variant).
if (a$stage == "contrasts") {
  need_packages("brms")
  name <- paste0(paste(a$time, a$variant, sep = "_"), if (a$patient_stim) "_pstim")
  x <- readRDS(file.path(a$out, "compass_model", paste0(name, ".rds"))); fit <- x$fit; d <- x$data
  pss <- sort(unique(d$ps)); ref <- paste0("V", cc$reference_dose)
  first_level <- function(v) factor(levels(d[[v]])[1], levels = levels(d[[v]]))
  if (a$time == "visit") {
    times <- data.frame(time = levels(d$visit), visit = factor(levels(d$visit), levels = levels(d$visit)))
    late <- intersect(paste0("V", cc$late_doses), times$time); tsl <- "stim:visit"
  } else {
    wk <- c(stats::median(d$week[d$stim == 1 & d$dose == cc$reference_dose]), stats::median(d$week[d$stim == 1 & d$dose %in% cc$late_doses]))
    times <- cbind(data.frame(time = c(ref, "late"), week = wk), week_basis_at(x$basis, wk)); late <- "late"
    tsl <- paste0("stim:", grep("^w[0-9]+$", names(times), value = TRUE), collapse = " + ")
  }
  nd <- merge(merge(data.frame(ps = pss, stringsAsFactors = FALSE), data.frame(stim = 0:1)), times)
  nd$parent <- factor(sub(":.*$", "", nd$ps), levels = levels(d$parent)); nd$parent_count <- 1L
  for (v in intersect(c("prior_vacc", design_technical_factors(cfg)), names(d))) if (is.factor(d[[v]])) nd[[v]] <- first_level(v)
  P <- brms::posterior_linpred(fit, newdata = nd, re_formula = stats::as.formula(sprintf("~ (1 + stim + %s | ps)", tsl)), transform = TRUE)
  col <- function(ps, tm, st) which(nd$ps == ps & nd$time == tm & nd$stim == st)
  lor <- function(ps, tm) qlogis(P[, col(ps, tm, 1)]) - qlogis(P[, col(ps, tm, 0)])
  dif <- function(ps, tm) P[, col(ps, tm, 1)] - P[, col(ps, tm, 0)]
  change <- function(f, ps) Reduce(`+`, lapply(late, function(tm) f(ps, tm))) / length(late) - f(ps, ref)
  mean_change <- function(f, s) Reduce(`+`, lapply(s, function(p) change(f, p))) / length(s)
  out <- list()
  add_row <- function(pp, q, x) out[[length(out) + 1L]] <<- cbind(data.frame(model = name, parent = pp, quantity = q), post_summary(x))
  for (pp in cc$parents) {
    s <- pss[startsWith(pss, paste0(pp, ":"))]; deg <- subset_degree(sub("^[^:]*:", "", s))
    poly <- mean_change(lor, s[deg >= 2]); single <- mean_change(lor, s[deg == 1])
    add_row(pp, "net response, late - ref (fraction of parent)", Reduce(`+`, lapply(s, function(p) change(dif, p))))
    add_row(pp, "mean log OR of >=2-marker subsets, late - ref", poly)
    add_row(pp, "mean log OR of single-marker subsets, late - ref", single)
    add_row(pp, "(>=2-marker minus single-marker) log OR, late - ref", poly - single)
    for (p in s) add_row(pp, paste("log OR, late - ref:", sub("^[^:]*:", "", p)), change(lor, p))
  }
  out <- do.call(rbind, out); out$late <- paste(late, collapse = "+"); out$ref <- ref
  write_table(out, file.path(a$out, paste0("exploratory_late_vs_ref_", name, ".csv")))
  print(out[!grepl("^log OR, late", out$quantity), ], digits = 3, row.names = FALSE)
  quit(status = 0)
}

# ---------------------------------------------------------------------------------------- summary
# Descriptive tables and figures from the COMPASS output (after --stage compass) and, when present,
# from the subset model fits (after --stage model): per patient x pool, per visit, figures.
if (a$stage == "summary") {
  figdir <- ensure_dir(file.path(a$out, "figures"))
  # COMPASS pairs every individual with its own control and treats individuals as independent. Here an
  # individual is a stimulated tube, so both assumptions are strained: one neg serves several peptide
  # pools of a visit, and one patient contributes many tubes. The figures say so.
  share_note <- function(d) {
    paste0("COMPASS treats each tube as an independent individual paired with its own control.\n",
           "Neither holds exactly here: ", sum(d$n_stim_sharing_neg > 1), " of ", nrow(d),
           " tubes share their unstimulated control with another stimulated tube\n(up to ",
           max(d$n_stim_sharing_neg), " tubes per control, ", length(unique(d$neg_filename)), " controls in all), and the ",
           nrow(d), " tubes come from only ", length(unique(d$patient)), " patients.\n",
           "Tubes sharing one control respond or fail together, so the rows are correlated in blocks.",
           if (any(d$low_info)) paste0("\n", sum(d$low_info), " tubes have fewer than 1,000 parent events in the tube or its control.") else "")
  }
  g <- read.csv(file.path(a$out, paste0("compass_gamma", sfx, ".csv")), stringsAsFactors = FALSE)
  sc <- read.csv(file.path(a$out, paste0("compass_scores", sfx, ".csv")), stringsAsFactors = FALSE)
  lab <- function(x) sub("&$", "", gsub("![A-Za-z0-9]+&?", "", x))   # "IFNg&!IL2&TNFa&!CD107a" -> "IFNg&TNFa"
  g$label <- lab(g$subset); g$resp <- g$mean_gamma > 0.5
  g <- merge(g, sc[, c("filename", "parent", "patient", "cocktail", "dose", "low_info")], by = c("filename", "parent"))
  # per patient x peptide pool: every tube of that pool, over visits
  pools <- do.call(rbind, lapply(split(g, paste(g$patient, g$cocktail, g$parent)), function(d) {
    s1 <- sc[sc$patient == d$patient[1] & sc$cocktail == d$cocktail[1] & sc$parent == d$parent[1], ]
    r <- unique(d$label[d$resp])
    data.frame(patient = d$patient[1], cocktail = d$cocktail[1], parent = d$parent[1], n_tubes = nrow(s1),
               doses = paste(sort(unique(s1$dose)), collapse = " "), FS_mean = mean(s1$FS), FS_max = max(s1$FS),
               PFS_mean = mean(s1$PFS), tubes_any_response = sum(s1$n_responding > 0), low_info_tubes = sum(s1$low_info),
               subsets_responding = paste(r[order(-lengths(strsplit(r, "&")), r)], collapse = "; "), stringsAsFactors = FALSE)
  }))
  pools <- pools[order(pools$parent, pools$patient, pools$cocktail), ]
  write_table(pools, file.path(a$out, paste0("compass_pools", sfx, ".csv")))
  # per parent x visit: tubes, patients, scores, responding subsets and their degree
  sc$any <- sc$n_responding > 0
  byv <- do.call(rbind, lapply(split(sc, paste(sc$parent, sc$dose)), function(d) {
    gg <- g[g$parent == d$parent[1] & g$dose == d$dose[1] & g$resp, ]
    pm <- aggregate(cbind(FS, PFS) ~ patient, d, mean)
    data.frame(parent = d$parent[1], dose = d$dose[1], n_tubes = nrow(d), n_patients = length(unique(d$patient)),
               FS_mean = mean(d$FS), PFS_mean = mean(d$PFS), FS_patient_mean = mean(pm$FS), PFS_patient_mean = mean(pm$PFS),
               pct_tubes_any_response = 100 * mean(d$any), mean_n_responding = mean(d$n_responding),
               share_poly_of_responding = if (nrow(gg)) mean(gg$degree >= 2) else NA_real_, stringsAsFactors = FALSE)
  }))
  byv <- byv[order(byv$parent, byv$dose), ]
  write_table(byv, file.path(a$out, paste0("compass_by_visit", sfx, ".csv")))
  print(byv, row.names = FALSE, digits = 3)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # PFS by visit: faint individual data (each patient x pool followed across visits, every tube a dot,
    # hollow when low-information), cohort mean of patient means per visit with a 95 % interval from
    # resampling patients in front. Means of patient means, so a patient with many pools counts once.
    ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e6e3"; faint <- "#c4c7cc"; dot <- "#8f949c"; accent <- "#2a78d6"
    sc$pool <- paste(sc$patient, sc$cocktail)
    pm <- aggregate(PFS ~ parent + patient + dose, sc, mean)
    # pool lines run through the pool's mean per visit: some pools were tested at several visits with the
    # same dose number (e.g. patient 23 at four V3 visits), which would otherwise draw vertical streaks
    plm <- aggregate(PFS ~ parent + patient + pool + dose, sc, mean)
    set.seed(cc$seed)
    cm <- do.call(rbind, lapply(split(pm, paste(pm$parent, pm$dose)), function(d) {
      boot <- replicate(2000, mean(sample(d$PFS, replace = TRUE)))
      data.frame(parent = d$parent[1], dose = d$dose[1], PFS = mean(d$PFS), lo = unname(stats::quantile(boot, 0.025)),
                 hi = unname(stats::quantile(boot, 0.975)), n_pat = nrow(d))
    }))
    nt <- aggregate(filename ~ dose, sc[sc$parent == cc$parents[1], ], length)
    np <- aggregate(patient ~ dose, sc[sc$parent == cc$parents[1], ], function(x) length(unique(x)))
    xlab <- setNames(sprintf("V%d\n%d patients\n%d tubes", nt$dose, np$patient, nt$filename), nt$dose)
    theme_pfs <- ggplot2::theme_minimal(base_size = 10.5) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_line(colour = grid_col, linewidth = 0.3),
                     axis.text = ggplot2::element_text(colour = ink2), axis.title = ggplot2::element_text(colour = ink2),
                     strip.text = ggplot2::element_text(colour = ink, face = "bold", size = 11, hjust = 0),
                     plot.title = ggplot2::element_text(colour = ink, face = "bold", size = 12.5),
                     plot.subtitle = ggplot2::element_text(colour = ink2, size = 10), plot.caption = ggplot2::element_text(colour = ink2, size = 8.5, hjust = 0),
                     plot.title.position = "plot", plot.caption.position = "plot", legend.position = "none",
                     plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    p <- ggplot2::ggplot(sc, ggplot2::aes(dose, PFS)) +
      ggplot2::geom_line(data = plm, ggplot2::aes(group = pool), colour = faint, linewidth = 0.3, alpha = 0.8) +
      ggplot2::geom_point(ggplot2::aes(shape = low_info), colour = dot, size = 1.1, alpha = 0.55,
                          position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 1)) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
      ggplot2::geom_ribbon(data = cm, ggplot2::aes(ymin = lo, ymax = hi), fill = accent, alpha = 0.14) +
      ggplot2::geom_line(data = cm, colour = accent, linewidth = 0.9) +
      ggplot2::geom_point(data = cm, colour = accent, fill = accent, size = 2.8, shape = 21, stroke = 0.6) +
      ggplot2::facet_wrap(~ parent, nrow = 1) +
      ggplot2::scale_x_continuous(breaks = nt$dose, labels = xlab[as.character(nt$dose)], expand = ggplot2::expansion(add = 0.25)) +
      ggplot2::scale_y_continuous(limits = c(0, 0.82), breaks = seq(0, 0.8, 0.2), expand = ggplot2::expansion(add = c(0.01, 0))) +
      ggplot2::labs(x = NULL, y = "COMPASS polyfunctionality score",
                    title = "CD4 responses are polyfunctional from the first dose, with no clear rise afterwards; CD8 stays near zero",
                    subtitle = "Blue: cohort mean per visit (mean of patient means) with 95 % interval from resampling patients",
                    caption = paste("Grey lines follow one patient's peptide pool across visits (mean per visit); grey dots are single stimulated tubes",
                                    "(hollow: fewer than 1,000 parent events in the tube or its neg).",
                                    "\nV0 is before vaccination. V3 and V4 come mostly from patients not seen earlier, so the line compares different people.")) +
      theme_pfs
    ggplot2::ggsave(file.path(figdir, paste0("16_pfs_by_visit", sfx, ".png")), p, width = 9.5, height = 5.4, dpi = 170)
    # one panel per patient: every peptide pool its own line, the patient's mean per visit in blue
    for (pp in cc$parents) {
      d <- sc[sc$parent == pp, ]
      first <- tapply(d$dose, d$patient, min); ord <- names(sort(first + as.numeric(names(first)) / 1e4))
      # star: prior peptide vaccination in any of the patient's tubes (patient 16's file names disagree)
      pv <- tapply(d$prior_vacc == 2, d$patient, any)
      lab_of <- function(pt) paste0("Patient ", pt, ifelse(pv[as.character(pt)], " *", ""))
      d$panel <- factor(lab_of(d$patient), levels = lab_of(ord))
      dm <- aggregate(PFS ~ panel + dose, d, mean)
      dl <- aggregate(PFS ~ panel + pool + dose, d, mean)
      q <- ggplot2::ggplot(d, ggplot2::aes(dose, PFS)) +
        ggplot2::geom_line(data = dl, ggplot2::aes(group = pool), colour = faint, linewidth = 0.4) +
        ggplot2::geom_point(ggplot2::aes(shape = low_info), colour = dot, size = 1.2, alpha = 0.8) +
        ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
        ggplot2::geom_line(data = dm, colour = accent, linewidth = 0.8) +
        ggplot2::geom_point(data = dm, colour = accent, size = 2) +
        ggplot2::facet_wrap(~ panel, ncol = 5) +
        ggplot2::scale_x_continuous(breaks = 0:4, labels = paste0("V", 0:4), limits = c(-0.2, 4.2)) +
        ggplot2::scale_y_continuous(limits = c(0, 0.82), breaks = c(0, 0.4, 0.8)) +
        ggplot2::labs(x = NULL, y = "COMPASS polyfunctionality score",
                      title = sprintf("%s T cells: each patient's peptide pools respond differently", pp),
                      subtitle = "Grey: one line per peptide pool across visits (dots: single tubes, hollow = fewer than 1,000 parent events). Blue: patient mean per visit.\nSeveral dots at one visit for one pool: that pool was tested at more than one visit with the same dose number.",
                      caption = "* previously peptide-vaccinated in at least one file name (patient 16's files disagree: MTB016_2 at V1/V2, MTB16_1 at V4). Panels ordered by the patient's first visit.") +
        theme_pfs + ggplot2::theme(strip.text = ggplot2::element_text(colour = ink, face = "bold", size = 9.5, hjust = 0),
                                   panel.spacing = ggplot2::unit(0.8, "lines"))
      ggplot2::ggsave(file.path(figdir, sprintf("16_pfs_by_patient_%s%s.png", pp, sfx)), q, width = 11, height = 9, dpi = 160)
    }
    comp <- aggregate(resp ~ parent + dose + degree, g, sum)
    p2 <- ggplot2::ggplot(comp, ggplot2::aes(factor(dose), resp, fill = factor(degree))) + ggplot2::geom_col(position = "fill") +
      ggplot2::scale_fill_manual(values = c(`1` = "#c9d3ee", `2` = "#8ea0d8", `3` = "#4b63b8", `4` = "#1d2f73"), name = "markers per cell") +
      ggplot2::facet_wrap(~ parent) + ggplot2::labs(x = "vaccine dose (visit)", y = "share of responding tube x subset calls",
      title = "Which kinds of subsets respond (COMPASS P > 0.5), by number of markers") + ggplot2::theme_minimal(base_size = 10)
    ggplot2::ggsave(file.path(figdir, paste0("16_responding_degree_by_visit", sfx, ".png")), p2, width = 8, height = 4.5, dpi = 150)
  }
  # subset model: stimulation effect (log odds ratio stimulated vs neg) per subset and time, and the hypotheses
  ef <- file.path(a$out, "subset_effects_visit_full.csv")
  if (file.exists(ef) && requireNamespace("ggplot2", quietly = TRUE)) {
    e <- read.csv(ef, stringsAsFactors = FALSE)
    e$label <- lab(e$subset); e$label <- factor(e$label, levels = rev(unique(e$label[order(-e$degree, e$label)])))
    e$txt <- ifelse(e$call == "up", "+", ifelse(e$call == "down", "-", ""))
    p3 <- ggplot2::ggplot(e, ggplot2::aes(time, label, fill = pmax(pmin(log_or, 3), -3))) + ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = txt), size = 4) +
      ggplot2::scale_fill_gradient2(low = "#b2593a", mid = "#f7f7f7", high = "#1d3b8f", limits = c(-3, 3), name = "log OR\nstim vs neg") +
      ggplot2::facet_wrap(~ parent) + ggplot2::labs(x = "visit (vaccine dose)", y = NULL,
      title = "Stimulation effect per Boolean subset and visit (subset model)", subtitle = "+ / -: P(up) >= 0.95 or <= 0.05") +
      ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank())
    ggplot2::ggsave(file.path(figdir, "16_subset_effects_visit.png"), p3, width = 8, height = 6, dpi = 150)
  }
  # the collected file matches the glob itself, so exclude it or every rerun doubles its rows
  hf <- setdiff(Sys.glob(file.path(a$out, "hypotheses_*.csv")), file.path(a$out, "hypotheses_all_models.csv"))
  if (length(hf)) {
    H <- do.call(rbind, lapply(hf, read.csv, stringsAsFactors = FALSE))
    write_table(H[H$hypothesis != "level", ], file.path(a$out, "hypotheses_all_models.csv"))
  }
  # heatmap per parent: tubes (rows, ordered by patient and dose) x subsets (columns, by degree)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- merge(read.csv(file.path(a$out, paste0("compass_gamma", sfx, ".csv")), stringsAsFactors = FALSE), sc[, c("filename", "parent", "patient", "dose", "prior_vacc", "instrument", "low_info")], by = c("filename", "parent"))
    for (pp in cc$parents) {
      d <- gg[gg$parent == pp, ]
      ord <- unique(d[order(d$patient, d$dose, d$filename), "filename"])
      d$row <- factor(d$filename, levels = rev(ord))
      d$subset <- factor(d$subset, levels = subsets_all)
      d$label <- sub("&$", "", gsub("![A-Za-z0-9]+&?", "", as.character(d$subset)))
      d$label <- factor(d$label, levels = unique(d$label[order(d$subset)]))
      p <- ggplot2::ggplot(d, ggplot2::aes(label, row, fill = mean_gamma)) + ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient(low = "#f4f6fb", high = "#1d3b8f", limits = c(0, 1), name = "P(response)") +
        ggplot2::facet_grid(paste0("V", dose) ~ ., scales = "free_y", space = "free_y") +
        ggplot2::labs(x = "Boolean subset (positive markers)", y = "stimulated tube (patient order)",
                      title = paste(pp, "T cells: COMPASS posterior probability of response per subset"),
                      subtitle = sprintf("Mean model over %d seeds. Each row is one stimulated tube, not one patient.", length(cc$seeds)),
                      caption = share_note(sc[sc$parent == pp, ])) +
        ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 50, hjust = 1), axis.text.y = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(angle = 0),
                       plot.caption = ggplot2::element_text(hjust = 0, size = 7.5, colour = "#52514e"),
                       plot.caption.position = "plot", plot.margin = ggplot2::margin(5, 5, 5, 40))
      ggplot2::ggsave(file.path(figdir, sprintf("16_compass_heatmap_%s%s.png", pp, sfx)), p, width = 8, height = 11.6, dpi = 150)
    }
  }
  quit(status = 0)
}
stop("unknown --stage ", a$stage, "\n", usage, call. = FALSE)
