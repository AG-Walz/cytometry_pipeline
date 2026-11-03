# src/harmonization/drift/drift_stats.R — pure statistics helpers, base R
# Distribution summaries for neg tubes (layer A), 1-D Earth Mover's Distance,
# bimodal valley finder, background subtraction and the positivity rule
# (layer B), control matching, background flips and the decision logic that
# combines model posteriors with the rule. No cytometry packages, so every
# function here is unit-tested locally.

#' arcsinh with cofactor; scatter is left linear by the caller.
asinh_cof <- function(x, cofactor) asinh(x / cofactor)

#' 1-D EMD between two samples = integral |F_a - F_b| = mean |Q_a - Q_b| on a
#' common probability grid. Exact for equal-weight samples, fast, no package.
emd_1d <- function(a, b, n_grid = 512L) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  p <- (seq_len(n_grid) - 0.5) / n_grid
  mean(abs(quantile(a, p, names = FALSE, type = 7) -
           quantile(b, p, names = FALSE, type = 7)))
}

#' Position of the density minimum between the two highest peaks.
#' Returns NA when the marker is not bimodal in this tube.
bimodal_valley <- function(x, bw = "nrd0", n = 512L, min_peak_frac = 0.10, max_valley_frac = 0.8) {
  x <- x[is.finite(x)]
  if (length(x) < 200L) return(NA_real_)
  d <- density(x, bw = bw, n = n)
  y <- d$y
  is_peak <- c(FALSE, diff(sign(diff(y))) < 0, FALSE)
  peaks <- which(is_peak & y >= min_peak_frac * max(y))
  if (length(peaks) < 2L) return(NA_real_)
  top2 <- sort(peaks[order(y[peaks], decreasing = TRUE)][1:2])
  seg <- top2[1]:top2[2]
  v <- seg[which.min(y[seg])]
  # a real valley is clearly lower than the smaller of the two peaks; density
  # wiggles on a unimodal marker are not
  if (y[v] > max_valley_frac * min(y[top2])) return(NA_real_)
  d$x[v]
}

#' Summary row for one marker vector of one tube.
marker_summary <- function(x, probs = c(0.10, 0.50, 0.90), valley = FALSE) {
  x <- x[is.finite(x)]
  q <- quantile(x, probs, names = FALSE)
  # valley / frac_pos are always present (NA for unimodal markers) so rows rbind
  out <- data.frame(n = length(x), median = median(x), mad = mad(x),
                    q10 = q[1], q90 = q[3], spread = q[3] - q[1],
                    valley = NA_real_, frac_pos = NA_real_,
                    stringsAsFactors = FALSE)
  if (valley) {
    v <- bimodal_valley(x)
    out$valley <- v
    out$frac_pos <- if (is.na(v)) NA_real_ else mean(x > v)
  }
  out
}

#' Background subtraction and positivity, vectorised over tubes.
#' `pct` and `neg_pct` are percentages of the same parent population.
positivity_call <- function(pct, neg_pct, min_subtracted_pct = 0.05,
                            min_fold_over_neg = 3) {
  sub <- pct - neg_pct
  fold <- ifelse(neg_pct > 0, pct / neg_pct, Inf)
  tol <- 1e-9   # 0.06 - 0.01 is 0.04999... in floating point
  data.frame(subtracted_pct = sub, fold_over_neg = fold,
             positive = (sub >= min_subtracted_pct - tol) & (fold >= min_fold_over_neg - tol))
}

#' Replicate noise: SD of a summary across the cocktail tubes of one
#' patient-day (same cells, same day, different peptides). For technical
#' markers the peptide should not matter, so this is pure assay noise.
replicate_sd <- function(value, patient_day) {
  s <- tapply(value, patient_day, function(v) if (length(v) > 1L) sd(v) else NA_real_)
  sqrt(mean(s^2, na.rm = TRUE))
}

#' Decision for one contrast (a day or the instrument) on one marker, with three states for rule 1
#' (Kruschke 2018, region of practical equivalence; Lakens 2017 for the frequentist analogue):
#'   detected      P(|effect| > margin) >= prob_threshold
#'   negligible    P(|effect| < margin) >= prob_threshold
#'   inconclusive  neither: the data cannot tell (most single-patient days)
#' `draws`: posterior draws of the effect. `rep_sd`: replicate SD on the same scale, kept for the
#' record. `margin`: the relevance margin (defaults to rep_sd, the earlier rule). `flips`: number of
#' positivity calls that change when the effect is removed (0 when not applicable, NA when unknown), or
#' the list returned by flip_summary(), whose noise baseline then decides rule 2:
#'   detected, correct                   expected flips above the 95th percentile of counting-noise flips
#'   detected, within counting noise     flips > 0, but no more than counting noise alone produces
#'   detected, not clinically relevant   no flips
drift_decision <- function(draws, rep_sd, flips, prob_threshold = 0.95, margin = rep_sd) {
  p_exceeds <- mean(abs(draws) > margin)
  p_within <- mean(abs(draws) < margin)
  detected <- is.finite(p_exceeds) && p_exceeds >= prob_threshold
  fs <- if (is.list(flips)) flips else list(mean = flips, lo = NA_real_, hi = NA_real_, noise_mean = NA_real_, noise_q95 = NA_real_)
  verdict <- if (!detected) { if (is.finite(p_within) && p_within >= prob_threshold) "negligible" else "inconclusive" }
             else if (is.na(fs$mean)) "detected, relevance unknown"
             else if (fs$mean == 0) "detected, not clinically relevant"
             else if (is.na(fs$noise_q95) || fs$mean > fs$noise_q95) "detected, correct"
             else "detected, within counting noise"
  data.frame(effect_mean = mean(draws),
             ci_low = unname(quantile(draws, 0.025)),
             ci_high = unname(quantile(draws, 0.975)),
             replicate_sd = rep_sd, margin = margin, p_exceeds_noise = p_exceeds, p_within_margin = p_within,
             positivity_flips = fs$mean, flips_lo = fs$lo, flips_hi = fs$hi,
             noise_flips_mean = fs$noise_mean, noise_flips_q95 = fs$noise_q95, verdict = verdict,
             stringsAsFactors = FALSE)
}

#' Relevance margin of rule 1 for a response of `layer` ("A", "BG", "B") on the count (logit) scale or
#' not. Counts: log(1 + cv) with cv = cfg$relevance$count_cv (default 0.35, the inter-laboratory CV of
#' ICS responses above 0.2 percent in proficiency panels, Jaimes et al. 2011; Maecker et al. 2005), since
#' a logit shift of a rare frequency is a log fold change. Other responses: the replicate SD, as before,
#' because no external standard exists for arcsinh medians or scatter. Returns list(margin, source).
relevance_margin <- function(cfg, counts, rep_sd) {
  cv <- if (is.null(cfg$relevance$count_cv)) 0.35 else as.numeric(cfg$relevance$count_cv)
  if (counts && is.finite(cv)) list(margin = log1p(cv), source = sprintf("count CV %.0f%%", 100 * cv))
  else list(margin = rep_sd, source = "replicate SD")
}

#' Rule 2 with uncertainty and a baseline. `flip_fun(effect)` returns the number of calls that change
#' when `effect` is removed; it is evaluated on `n_draws` posterior draws (expected flips and a 95 percent
#' interval). `noise_fun()` returns the flips produced by one binomial resample of the affected counts
#' at their observed proportions, with no effect at all; `n_noise` resamples give its mean and 95th
#' percentile. Returns list(mean, lo, hi, noise_mean, noise_q95).
flip_summary <- function(draws, flip_fun, noise_fun = NULL, n_draws = 200L, n_noise = 200L, seed = 1L) {
  idx <- unique(round(seq(1, length(draws), length.out = min(n_draws, length(draws)))))
  f <- vapply(draws[idx], flip_fun, numeric(1))
  nz <- if (is.null(noise_fun)) NA_real_ else { set.seed(seed); vapply(seq_len(n_noise), function(i) noise_fun(), numeric(1)) }
  list(mean = mean(f), lo = unname(quantile(f, 0.025)), hi = unname(quantile(f, 0.975)),
       noise_mean = mean(nz), noise_q95 = if (all(is.na(nz))) NA_real_ else unname(quantile(nz, 0.95)))
}

#' Calls among rows `which` of `fb` (one response) that change under one binomial resample of the
#' stimulated counts (`stim = TRUE`) and of the matched neg counts (`neg = TRUE`) at their observed
#' proportions: the counting-noise baseline for rule 2.
counting_noise_flips <- function(fb, which, min_subtracted_pct, min_fold_over_neg, stim = TRUE, neg = TRUE) {
  d <- fb[which, , drop = FALSE]
  if (!nrow(d)) return(0)
  before <- positivity_call(d$pct, d$neg_pct, min_subtracted_pct, min_fold_over_neg)$positive
  k <- if (stim) stats::rbinom(nrow(d), d$parent_count, d$count / d$parent_count) else d$count
  nn <- round(d$neg_parent_count)
  k0 <- if (neg) stats::rbinom(nrow(d), nn, pmin(1, d$neg_count / d$neg_parent_count)) else d$neg_count
  after <- positivity_call(100 * k / d$parent_count, 100 * k0 / d$neg_parent_count, min_subtracted_pct, min_fold_over_neg)$positive
  sum(before != after, na.rm = TRUE)
}

#' Per stimulated call, the probability that it changes under counting noise alone (binomial resampling
#' of the stimulated and matched neg counts at their observed proportions, `n` times), with the parent
#' count. Used for the call-stability table and the minimum-parent-count sensitivity.
call_flip_probability <- function(st, min_subtracted_pct, min_fold_over_neg, n = 500L, seed = 1L) {
  set.seed(seed)
  obs <- positivity_call(st$pct, st$neg_pct, min_subtracted_pct, min_fold_over_neg)$positive
  nn <- round(st$neg_parent_count)
  fl <- replicate(n, {
    k <- stats::rbinom(nrow(st), st$parent_count, st$count / st$parent_count)
    k0 <- stats::rbinom(nrow(st), nn, pmin(1, st$neg_count / st$neg_parent_count))
    positivity_call(100 * k / st$parent_count, 100 * k0 / st$neg_parent_count, min_subtracted_pct, min_fold_over_neg)$positive != obs
  })
  rowMeans(matrix(fl, nrow = nrow(st)), na.rm = TRUE)
}

#' Control (neg) matching as a policy. For every stimulated row of `fq` (the
#' step 5 table before subtraction) the neg of the same parent and gate is
#' found by exact match on `controls$match_keys`, then by the fallbacks in
#' `controls$fallbacks` in order: same_patient_day, same_patient_reference_instrument,
#' pooled_reference_instrument. A neg with an empty parent gate counts as
#' missing. Several candidates are pooled by summing count and parent_count.
#' Adds neg_count, neg_parent_count, neg_pct, neg_match ("exact", the fallback
#' name or "none"), neg_filename (comma-joined when pooled), neg_acq_day and
#' neg_instrument (NA when pooled across days / instruments). Neg rows get NA.
match_controls <- function(fq, controls, reference_instrument) {
  keys <- unlist(controls$match_keys); fallbacks <- unlist(controls$fallbacks)
  known <- c("same_patient_day", "same_patient_reference_instrument", "pooled_reference_instrument")
  if (length(bad <- setdiff(fallbacks, known))) stop("unknown control fallback: ", paste(bad, collapse = ", "))
  one <- function(x) if (length(unique(x)) == 1L) x[1] else NA
  neg <- fq[fq$is_neg & !is.na(fq$parent_count) & fq$parent_count > 0, ]
  res <- data.frame(neg_count = rep(NA_real_, nrow(fq)), neg_parent_count = NA_real_, neg_pct = NA_real_,
                    neg_match = NA_character_, neg_filename = NA_character_, neg_acq_day = NA_character_,
                    neg_instrument = NA_character_, stringsAsFactors = FALSE)
  key <- function(d) do.call(paste, c(d[keys], sep = "|"))
  for (i in which(!fq$is_neg)) {
    r <- fq[i, ]
    cand <- neg[neg$parent == r$parent & neg$gate == r$gate, ]
    hit <- cand[key(cand) == key(r), ]; label <- "exact"
    if (!nrow(hit)) for (fb in fallbacks) {
      hit <- switch(fb,
        same_patient_day = cand[cand$patient == r$patient & cand$acq_day == r$acq_day, ],
        same_patient_reference_instrument = cand[cand$patient == r$patient & cand$instrument == reference_instrument, ],
        pooled_reference_instrument = cand[cand$instrument == reference_instrument, ])
      if (nrow(hit)) { label <- fb; break }
    }
    if (!nrow(hit)) { res$neg_match[i] <- "none"; next }
    res$neg_count[i] <- sum(hit$count); res$neg_parent_count[i] <- sum(hit$parent_count)
    res$neg_pct[i] <- 100 * res$neg_count[i] / res$neg_parent_count[i]
    res$neg_match[i] <- label; res$neg_filename[i] <- paste(hit$filename, collapse = ",")
    res$neg_acq_day[i] <- one(hit$acq_day); res$neg_instrument[i] <- one(hit$instrument)
  }
  cbind(fq, res)
}

#' Level of technical factor `term` of the neg matched to each row of `fb`, looked up
#' through neg_filename in the tube table `tt`; NA when the pooled negs disagree.
#' "instrument" and "day" use the neg_instrument and neg_acq_day columns.
matched_neg_level <- function(fb, term, tt) {
  if (term == "instrument") return(fb$neg_instrument)
  if (term == "day") return(fb$neg_acq_day)
  lv <- setNames(as.character(tt[[term]]), tt$filename)
  vapply(strsplit(as.character(fb$neg_filename), ","), function(f) {
    x <- unique(lv[f]); if (length(x) == 1L) unname(x) else NA_character_ }, character(1))
}

#' Positivity calls that change when a background (neg) effect on the logit
#' scale is removed from the matched negs. `fb` is a layer B frame for ONE
#' response; `which_neg` marks the rows whose matched neg carries the effect
#' (its neg_acq_day or neg_instrument is the affected level). Only neg_pct is
#' shifted, the stimulated frequency stays as measured.
count_background_flips <- function(fb, effect, which_neg, min_subtracted_pct = 0.05, min_fold_over_neg = 3) {
  before <- positivity_call(fb$pct, fb$neg_pct, min_subtracted_pct, min_fold_over_neg)$positive
  pn <- pmin(pmax(fb$neg_pct / 100, 1e-6), 1 - 1e-6)
  adjn <- pn; adjn[which_neg] <- plogis(qlogis(pn[which_neg]) - effect)
  after <- positivity_call(fb$pct, 100 * adjn, min_subtracted_pct, min_fold_over_neg)$positive
  sum(before != after, na.rm = TRUE)
}

#' Responder calls of every stimulated row of `st` (gated_frequencies rows with a matched neg and
#' a `response` column) under the lab rule and a one-sided Fisher's exact test (stimulated vs matched
#' neg counts, Benjamini-Hochberg within each response, q < `alpha`). A `mimosa` column is added as NA.
calling_rules <- function(st, cfg, alpha = 0.05) {
  pos <- cfg$positivity
  lab <- positivity_call(st$pct, st$neg_pct, pos$min_subtracted_pct, pos$min_fold_over_neg)$positive
  p <- mapply(function(k, n, k0, n0) stats::fisher.test(matrix(c(k, n - k, k0, n0 - k0), 2), alternative = "greater")$p.value,
              st$count, st$parent_count, st$neg_count, st$neg_parent_count)
  q <- ave(p, st$response, FUN = function(x) stats::p.adjust(x, "BH"))
  data.frame(filename = st$filename, response = st$response, patient = st$patient, week = st$week,
             pct = st$pct, neg_pct = st$neg_pct, lab = lab, fisher_q = q, fisher = q < alpha,
             mimosa = NA, stringsAsFactors = FALSE)
}

#' Per response: number of positive calls under each rule and Cohen's kappa of every pair of rules.
calls_agreement <- function(calls) {
  rules <- intersect(c("lab", "fisher", "mimosa"), names(calls))
  rules <- rules[vapply(rules, function(r) any(!is.na(calls[[r]])), logical(1))]
  kappa <- function(x, y) { ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; po <- mean(x == y)
    pe <- mean(x) * mean(y) + (1 - mean(x)) * (1 - mean(y)); if (pe == 1) NA else (po - pe) / (1 - pe) }
  do.call(rbind, lapply(split(calls, calls$response), function(d) {
    out <- data.frame(response = d$response[1], n = nrow(d), stringsAsFactors = FALSE)
    for (r in rules) out[[paste0("positive_", r)]] <- sum(d[[r]], na.rm = TRUE)
    if (length(rules) > 1) for (pr in utils::combn(rules, 2, simplify = FALSE))
      out[[paste0("kappa_", pr[1], "_", pr[2])]] <- kappa(d[[pr[1]]], d[[pr[2]]])
    out
  }))
}

mimosa_available <- function() requireNamespace("MIMOSA", quietly = TRUE)

#' MIMOSA calls (Finak et al. 2014) for the stimulated rows of `st`: one pair per stimulated tube
#' (its counts and its matched neg's counts), all pools of a response in one EM fit, since pools are
#' patient-specific and MIMOSA shares information across the pairs of one stimulation. Returns
#' data.frame(prob, fdr, positive) in the row order of `st`; positive = fdr < `alpha`.
mimosa_calls <- function(st, alpha = 0.05) {
  out <- data.frame(prob = rep(NA_real_, nrow(st)), fdr = NA_real_, positive = NA)
  st$UID <- sprintf("u%05d", seq_len(nrow(st)))
  for (resp in unique(st$response)) {
    i <- which(st$response == resp); s <- st[i, ]
    long <- rbind(data.frame(UID = s$UID, ANTIGEN = "POOL", CYTNUM = s$count, NSUB = s$parent_count - s$count),
                  data.frame(UID = s$UID, ANTIGEN = "negctrl", CYTNUM = round(s$neg_count), NSUB = round(s$neg_parent_count - s$neg_count)))
    long$CYTOKINE <- "c"; long$TCELLSUBSET <- "t"
    res <- tryCatch(suppressMessages(suppressWarnings({
      # MIMOSA evaluates `reference` and `.variables` unevaluated (plyr style), so they are written out here
      E <- MIMOSA::ConstructMIMOSAExpressionSet(long, reference = ANTIGEN %in% "negctrl", measure.columns = c("CYTNUM", "NSUB"),
             other.annotations = c("CYTOKINE", "TCELLSUBSET", "ANTIGEN", "UID"),
             default.cast.formula = component ~ UID + ANTIGEN + CYTOKINE + TCELLSUBSET,
             .variables = plyr::.(TCELLSUBSET, CYTOKINE, UID), featureCols = 1, ref.append.replace = "_REF")
      fit <- MIMOSA::MIMOSA(NSUB + CYTNUM ~ UID + TCELLSUBSET + CYTOKINE | ANTIGEN, data = E, method = "EM",
                            subset = RefTreat %in% "Treatment" & ANTIGEN %in% "POOL", ref = ANTIGEN %in% "POOL" & RefTreat %in% "Reference")
      z <- MIMOSA::getZ(fit); q <- MIMOSA::fdr(fit)   # one row per pair, in the order of pData(fit[[1]])
      uid <- as.character(Biobase::pData(fit[[1]])$UID)
      list(prob = setNames(z[, "Pr.response"], uid), fdr = setNames(as.numeric(q), uid))
    })), error = function(e) { cat("  MIMOSA failed for", resp, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(res)) next
    out$prob[i] <- res$prob[s$UID]; out$fdr[i] <- res$fdr[s$UID]
  }
  out$positive <- out$fdr < alpha
  out
}
