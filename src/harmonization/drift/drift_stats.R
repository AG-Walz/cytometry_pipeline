# src/harmonization/drift/drift_stats.R — pure statistics helpers, base R
# Distribution summaries for neg tubes (layer A), 1-D Earth Mover's Distance,
# bimodal valley finder, background subtraction and the positivity rule
# (layer B), and the decision logic that combines model posteriors with the
# rule. No cytometry packages, so every function here is unit-tested locally.

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

#' Decision for one contrast (a day or the instrument) on one marker.
#' `draws`: posterior draws of the effect. `rep_sd`: replicate SD on the same
#' scale. `flips`: number of positivity calls that change when the effect is
#' removed (0 when not applicable, NA when unknown).
drift_decision <- function(draws, rep_sd, flips, prob_threshold = 0.95) {
  p_exceeds <- mean(abs(draws) > rep_sd)
  detected  <- is.finite(p_exceeds) && p_exceeds >= prob_threshold
  verdict <- if (!detected) "not detected"
             else if (is.na(flips)) "detected, relevance unknown"
             else if (flips > 0) "detected, correct"
             else "detected, not clinically relevant"
  data.frame(effect_mean = mean(draws),
             ci_low = unname(quantile(draws, 0.025)),
             ci_high = unname(quantile(draws, 0.975)),
             replicate_sd = rep_sd, p_exceeds_noise = p_exceeds,
             positivity_flips = flips, verdict = verdict,
             stringsAsFactors = FALSE)
}
