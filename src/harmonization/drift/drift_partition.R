# src/harmonization/drift/drift_partition.R — where do the differences between tubes come from?
# Base R (lme4 and brms only inside the fit helpers). docs/plans/drift_three_questions_plan.md, part D.
# Every readout (an MFI of a marker in a population, or a cytokine gate in a parent) is fitted once per
# data representation with the same factors, and the variance of the latent predictor is split into
# one share per factor:
#   unwanted: the technical factors (instrument, operator, ...) and the acquisition day
#   wanted:   the unit (patient), time (linear + spline), protected covariates, stimulated or not, pool
#   residual: tube-level variance left over
# Fixed factors contribute the variance of their part of the linear predictor across tubes, random
# factors their variance (Nakagawa & Schielzeth 2013, Nakagawa, Johnson & Schielzeth 2017).
# Count representations are binomial with one normal random effect per tube on the logit scale (in
# lme4 and brms alike); that tube-level variance is the residual, so binomial counting noise is not
# part of the total (as in Finak et al. 2016, which models logit proportions). A beta-binomial was
# tried first: its exact logit-scale tube variance explodes for rare events (unstimulated tubes near
# 0.01 percent) and is not comparable with lme4. For calls the residual is the latent-logistic pi^2/3.

#' The data representations. `kind` picks the readout table, `tubes` whether the unstimulated
#' tubes are rows of the model (then `stimulated` is a factor), `family` the likelihood.
#' The grid: MFI {linear, arcsinh, FlowJo display scale (each sample's own transformation from the
#' workspace)} x {plain, stim minus neg}; frequency {%, arcsinh %, logit counts} x
#' {plain, stim minus neg}; the responder call. `short` is the axis label of the comparison figures.
PARTITION_REPS <- list(
  mfi_raw        = list(kind = "mfi",  tubes = "all",  family = "gaussian",  label = "MFI, linear",                       short = "MFI linear"),
  mfi_raw_sub    = list(kind = "mfi",  tubes = "stim", family = "gaussian",  label = "MFI, linear, stim minus neg",        short = "MFI linear\nstim-neg"),
  mfi_asinh      = list(kind = "mfi",  tubes = "all",  family = "gaussian",  label = "MFI, arcsinh",                      short = "MFI arcsinh"),
  mfi_sub        = list(kind = "mfi",  tubes = "stim", family = "gaussian",  label = "MFI, arcsinh, stim minus neg",       short = "MFI arcsinh\nstim-neg"),
  mfi_flowjo     = list(kind = "mfi",  tubes = "all",  family = "gaussian",  label = "MFI, FlowJo display scale",         short = "MFI FlowJo"),
  mfi_flowjo_sub = list(kind = "mfi",  tubes = "stim", family = "gaussian",  label = "MFI, FlowJo scale, stim minus neg",  short = "MFI FlowJo\nstim-neg"),
  freq_pct       = list(kind = "freq", tubes = "all",  family = "gaussian",  label = "frequency, %",                      short = "freq %"),
  freq_sub       = list(kind = "freq", tubes = "stim", family = "gaussian",  label = "frequency, stim % minus neg %",      short = "freq %\nstim-neg"),
  freq_asinh     = list(kind = "freq", tubes = "all",  family = "gaussian",  label = "frequency, arcsinh %",              short = "freq arcsinh"),
  freq_asinh_sub = list(kind = "freq", tubes = "stim", family = "gaussian",  label = "frequency, arcsinh %, stim minus neg", short = "freq arcsinh\nstim-neg"),
  freq_raw       = list(kind = "freq", tubes = "all",  family = "binomial",  label = "frequency, logit (gated counts)",   short = "freq logit"),
  freq_ratio     = list(kind = "freq", tubes = "stim", family = "binomial",  label = "frequency, logit, stim over neg",    short = "freq logit\nstim/neg"),
  call           = list(kind = "freq", tubes = "stim", family = "bernoulli", label = "responder call (lab rule)",          short = "call",
                        compare = FALSE)
)

#' Settings that enter the comparison figures and tables. The responder call is a decision made from the
#' data, not a way of scaling it, so it is fitted (for the record) but not compared with the settings.
compared_reps <- function() names(Filter(function(x) !identical(x$compare, FALSE), PARTITION_REPS))

#' Cofactor (in percent) of the arcsinh frequency representations: linear below it, logarithmic
#' above. Default 0.05 percent, the minimum subtracted frequency of the lab positivity rule.
freq_asinh_cofactor <- function(cfg) if (is.null(cfg$freq_asinh_cofactor)) 0.05 else cfg$freq_asinh_cofactor

elogit <- function(k, n) qlogis((k + 0.5) / (n + 1))

#' Matched neg value per stimulated tube for MFI: the mean over the negs of the same
#' `match_keys` (default patient + week). NA when there is none.
mfi_matched_neg <- function(d, tt, match_keys = c("patient", "week")) {
  key <- function(x) do.call(paste, c(x[match_keys], sep = "|"))
  negs <- d[d$filename %in% tt$filename[tt$is_neg %in% TRUE], ]
  m <- tapply(negs$value, key(negs), mean)
  unname(m[key(d)])
}

#' Long model frame of one representation, all readouts stacked (column `readout`).
#' `mfi`: mfi.csv; `fq`: gated_frequencies.csv; `tt`: tubes.csv. Rows lacking the value
#' are dropped. Adds `stimulated` ("no"/"yes") and `pool` (cocktail, "neg" for negs).
partition_frame <- function(rep, mfi, fq, tt, cfg, ref = NULL, quality = NULL) {
  spec <- PARTITION_REPS[[rep]]
  if (is.null(spec)) stop("unknown representation: ", rep, call. = FALSE)
  keys <- unlist(cfg$controls$match_keys); if (!length(keys)) keys <- c("patient", "week")
  if (spec$kind == "mfi") {
    if (is.null(mfi)) return(NULL)
    d <- merge(mfi, tt[, intersect(unique(c("filename", "patient", "week", "is_neg", "cocktail", "acq_day", "patient_day",
                                            design_columns(cfg))), names(tt))], by = "filename")
    d$readout <- paste(d$population, d$marker, sep = ":")
    if (rep %in% c("mfi_flowjo", "mfi_flowjo_sub") && !"median_flowjo" %in% names(d)) return(NULL)
    d$value <- if (rep %in% c("mfi_raw", "mfi_raw_sub")) d$median_linear
               else if (rep %in% c("mfi_flowjo", "mfi_flowjo_sub")) d$median_flowjo else d$median_asinh
    if (rep %in% c("mfi_sub", "mfi_raw_sub", "mfi_flowjo_sub")) {
      d <- do.call(rbind, lapply(split(d, d$readout), function(x) { x$neg_value <- mfi_matched_neg(x, tt, keys); x }))
      d <- d[!d$is_neg & !is.na(d$neg_value), ]
      d$value <- d$value - d$neg_value
    }
    d$y <- d$value
  } else {
    d <- fq[!is.na(fq$parent_count) & fq$parent_count > 0, ]
    d <- with_tube_columns(d, NULL, cfg, tt)
    d$readout <- paste(d$parent, d$gate, sep = ":")
    if (spec$tubes == "stim") d <- d[!d$is_neg & !is.na(d$neg_pct), ]
    cf <- freq_asinh_cofactor(cfg)
    if (rep == "freq_pct") d$y <- d$pct
    if (rep == "freq_sub") d$y <- d$pct - d$neg_pct
    if (rep == "freq_asinh") d$y <- asinh(d$pct / cf)
    if (rep == "freq_asinh_sub") d$y <- asinh(d$pct / cf) - asinh(d$neg_pct / cf)
    if (rep == "freq_ratio") { d <- d[!is.na(d$neg_parent_count) & d$neg_parent_count > 0, ]; d$neg_offset <- elogit(d$neg_count, d$neg_parent_count) }
    if (rep == "call") d$y <- as.integer(d$positive %in% TRUE)
  }
  if (spec$tubes == "all") d <- d[!is.na(d$is_neg), ] else d <- d[d$is_neg %in% FALSE, ]
  d$stimulated <- factor(ifelse(d$is_neg, "no", "yes"), c("no", "yes"))
  d$pool <- ifelse(d$is_neg, "neg", as.character(d$cocktail))
  if (spec$family == "gaussian") d <- d[is.finite(d$y), ]
  # visit = one unit at one time point (patient x week): the sample all tubes of that visit come from
  ds <- cfg$design
  if (!is.null(ds$time) && ds$time %in% names(d)) d$visit <- paste(d[[ds$unit]], d[[ds$time]], sep = "_")
  d <- with_visit_quality(d, quality, tt, keys)
  d <- design_frame(d, cfg, ref)
  d$pool <- factor(d$pool)
  if ("visit" %in% names(d)) d$visit <- factor(d$visit)
  d$representation <- rep
  d
}

#' Sample quality per visit, measured on the visit's unstimulated tube (the peptide does not change it, so it
#' describes the thawed sample): `quality` holds filename, viability, live_cells (sample_quality.csv of the
#' acquisition audit). Adds qual_viability (logit, standardised) and qual_cells (log, standardised); a visit
#' without an unstimulated tube takes the mean of its tubes. Unchanged when `quality` is NULL.
with_visit_quality <- function(d, quality, tt, keys = c("patient", "week")) {
  if (is.null(quality) || !nrow(d)) return(d)
  q <- merge(tt[, intersect(c("filename", "is_neg", keys), names(tt))], quality[, c("filename", "viability", "live_cells")], by = "filename")
  key <- function(x) do.call(paste, c(x[keys], sep = "|"))
  q$k <- key(q)
  negq <- q[q$is_neg %in% TRUE, ]; negq <- negq[!duplicated(negq$k), ]
  allq <- aggregate(cbind(viability, live_cells) ~ k, q, mean)
  pick <- function(col) { v <- negq[[col]][match(key(d), negq$k)]; ifelse(is.na(v), allq[[col]][match(key(d), allq$k)], v) }
  via <- qlogis(pmin(pmax(pick("viability"), 0.01), 0.99)); cells <- log(pmax(pick("live_cells"), 1))
  z <- function(x) { s <- stats::sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s }   # constant -> 0
  d$qual_viability <- z(via); d$qual_cells <- z(cells)
  d[is.finite(d$qual_viability) & is.finite(d$qual_cells), ]
}

#' Factors of the partition model for one representation and frame: list(unwanted, wanted)
#' of term names, and the formula pieces. A factor with fewer than two levels in the frame
#' is left out (and so is reported as absent, not as zero).
partition_terms <- function(cfg, d, rep) {
  ds <- cfg$design; spec <- PARTITION_REPS[[rep]]
  lv <- function(v) v %in% names(d) && nlevels(droplevels(factor(d[[v]]))) >= 2L
  tech_fixed <- Filter(lv, design_technical_factors(cfg))
  day <- design_day_term(cfg); day <- if (!is.na(day) && lv(day)) day else character(0)
  prot <- Filter(lv, unlist(ds$protected_fixed))
  time <- if (!is.null(ds$time) && "week_c" %in% names(d)) c("week_c", time_nl_terms(cfg)) else character(0)
  stim <- if (spec$tubes == "all" && lv("stimulated")) "stimulated" else character(0)
  pool <- if (lv("pool")) "pool" else character(0)
  qual <- if (all(c("qual_viability", "qual_cells") %in% names(d))) c("qual_viability", "qual_cells") else character(0)
  # visit: only with several tubes per visit, otherwise it is the residual
  visit <- if (lv("visit") && nlevels(droplevels(factor(d$visit))) < nrow(d) * 0.8) "visit" else character(0)
  list(fixed = c(tech_fixed, qual, prot, time, stim), random = c(day, ds$unit, pool, visit),
       groups = c(setNames(rep("unwanted", length(tech_fixed) + length(day)), c(tech_fixed, day)),
                  if (length(qual)) c(quality = "unwanted"),
                  setNames(rep("wanted", 1L + length(prot) + length(stim) + length(pool)), c(ds$unit, prot, stim, pool)),
                  if (length(time)) c(week = "wanted"),
                  if (length(visit)) c(visit = "visit")))
}

#' Model formula of the partition fit. engine "lme4" or "brms"; counts get the observation-level
#' (1 | filename) term in both engines, which is then the residual (fit them as binomial).
partition_formula <- function(terms, rep, engine = c("lme4", "brms")) {
  engine <- match.arg(engine); fam <- PARTITION_REPS[[rep]]$family
  lhs <- switch(fam, gaussian = "y", bernoulli = "y",
                binomial = if (engine == "brms") "count | trials(parent_count)" else "cbind(count, parent_count - count)")
  rhs <- c("1", terms$fixed, sprintf("(1 | %s)", terms$random),
           if (fam == "binomial") "(1 | filename)",
           if (rep == "freq_ratio") "offset(neg_offset)")
  as.formula(paste(lhs, "~", paste(rhs, collapse = " + ")), env = globalenv())
}

#' Which fixed term each column of a model matrix belongs to: the time columns map to
#' "week", every other column to the factor whose name it starts with.
fixed_term_of <- function(colnames_x, fixed) {
  vapply(colnames_x, function(cn) {
    if (cn == "(Intercept)") return(NA_character_)
    if (cn == "week_c" || grepl("^week_nl[0-9]+$", cn)) return("week")
    if (grepl("^qual_", cn)) return("quality")
    hit <- fixed[startsWith(cn, fixed)]
    if (length(hit)) hit[which.max(nchar(hit))] else NA_character_
  }, character(1))
}

#' Variance shares from one set of components. `contrib`: named list of per-tube vectors (the
#' fixed part of each factor); `rand`: named random-effect variances; `resid`: residual variance;
#' `groups`: factor -> "unwanted"/"wanted". Returns data.frame(factor, group, variance, share)
#' including the residual row; shares sum to 1.
variance_shares <- function(contrib, rand, resid, groups) {
  v <- c(vapply(contrib, function(x) stats::var(x), numeric(1)), rand)
  v <- c(v[names(v) %in% names(groups)], residual = resid)
  tot <- sum(v)
  data.frame(factor = names(v), group = ifelse(names(v) == "residual", "residual", groups[names(v)]),
             variance = unname(v), share = unname(v / tot), stringsAsFactors = FALSE)
}

#' Per-factor fixed contributions from a model matrix X and coefficients b (named like X).
fixed_contrib <- function(X, b, terms) {
  ft <- fixed_term_of(colnames(X), terms$fixed)
  out <- list()
  for (f in unique(ft[!is.na(ft)])) {
    cols <- colnames(X)[ft %in% f & colnames(X) %in% names(b)]
    out[[f]] <- as.numeric(X[, cols, drop = FALSE] %*% b[cols])
  }
  out
}

#' Shares from an lme4 fit (lmer or glmer). Returns the variance_shares() table.
lme4_shares <- function(fit, terms, rep) {
  X <- lme4::getME(fit, "X"); b <- lme4::fixef(fit)
  vc <- as.data.frame(lme4::VarCorr(fit)); vc <- vc[is.na(vc$var2), ]
  rand <- setNames(vc$vcov, vc$grp)
  fam <- PARTITION_REPS[[rep]]$family
  resid <- switch(fam, gaussian = rand[["Residual"]], bernoulli = pi^2 / 3,
                  binomial = if ("filename" %in% names(rand)) rand[["filename"]] else 0)
  rand <- rand[!names(rand) %in% c("Residual", "filename")]
  variance_shares(fixed_contrib(X, b, terms), rand, resid, terms$groups)
}

#' Shares per posterior draw from a brms fit; returns one variance_shares() table per draw,
#' stacked with a `draw` column. Count models (binomial with a per-tube random effect) use that
#' random effect's variance as the residual.
brms_shares <- function(fit, terms, rep, d, n_draws = 1000L) {
  dr <- as.data.frame(fit)
  keep <- unique(round(seq(1, nrow(dr), length.out = min(n_draws, nrow(dr)))))
  dr <- dr[keep, , drop = FALSE]
  X <- stats::model.matrix(as.formula(paste("~", paste(c("1", terms$fixed), collapse = " + "))), d)
  clean <- function(x) gsub("[^[:alnum:]_]", "", x)
  bcols <- setNames(paste0("b_", clean(colnames(X))), colnames(X))
  fam <- PARTITION_REPS[[rep]]$family
  rows <- lapply(seq_len(nrow(dr)), function(i) {
    b <- setNames(vapply(bcols, function(cc) if (cc %in% names(dr)) dr[i, cc] else 0, numeric(1)), colnames(X))
    rand <- vapply(terms$random, function(g) { cc <- sprintf("sd_%s__Intercept", g); if (cc %in% names(dr)) dr[i, cc]^2 else 0 }, numeric(1))
    resid <- switch(fam, gaussian = dr[i, "sigma"]^2, bernoulli = pi^2 / 3,
                    binomial = dr[i, "sd_filename__Intercept"]^2)
    cbind(draw = i, variance_shares(fixed_contrib(X, b, terms), rand, resid, terms$groups))
  })
  do.call(rbind, rows)
}

#' Summary rows (mean and 95 percent interval of each share, mean variance) of a stacked draws
#' table, plus the unwanted and wanted totals per draw. `scale_var` multiplies the variances back to
#' the response scale when the model was fitted on a standardised response.
summarise_shares <- function(sh, scale_var = 1) {
  if (!"draw" %in% names(sh)) sh$draw <- 1L
  tot <- aggregate(cbind(share, variance) ~ draw + group, sh[sh$group != "residual", ], sum)
  tot$factor <- paste0("total_", tot$group)
  all <- rbind(sh[, c("draw", "factor", "group", "share", "variance")], tot[, c("draw", "factor", "group", "share", "variance")])
  s <- split(all$share, all$factor); v <- split(all$variance, all$factor)
  grp <- tapply(all$group, all$factor, `[`, 1)
  data.frame(factor = names(s), group = unname(grp[names(s)]), share = vapply(s, mean, numeric(1)),
             variance = vapply(v, mean, numeric(1)) * scale_var,
             lo = vapply(s, function(x) unname(quantile(x, 0.025)), numeric(1)),
             hi = vapply(s, function(x) unname(quantile(x, 0.975)), numeric(1)), stringsAsFactors = FALSE, row.names = NULL)
}

#' Technical CV implied by a variance on a log-like scale (arcsinh of bright signals, logit of
#' small proportions): sqrt(exp(v) - 1), the log-normal CV. Used against the CLSI H62 / O'Hara targets.
technical_cv <- function(v) sqrt(exp(v) - 1)

#' The partition table of a results folder: the brms one when it exists, else lme4. Full variant only
#' unless `variant` is NULL.
read_partition <- function(out, variant = "full") {
  for (e in c("brms", "lme4")) {
    p <- file.path(out, sprintf("variance_partition_%s.csv", e))
    if (file.exists(p)) { v <- read.csv(p, stringsAsFactors = FALSE); if (!is.null(variant)) v <- v[v$variant == variant, ]; return(v) }
  }
  NULL
}

#' Question 2: does an alternative representation of the same readout remove unwanted variance
#' while keeping the wanted variance? One row per readout x alternative with the unwanted and wanted
#' totals of both and a verdict. "keeps wanted" means the alternative's upper interval bound of the
#' wanted share reaches the base's estimate (lme4 tables have lo = hi = share).
subtraction_table <- function(vp, base = "freq_raw", alts = c("freq_pct", "freq_sub", "freq_asinh", "freq_asinh_sub", "freq_ratio")) {
  tot <- vp[vp$factor %in% c("total_unwanted", "total_wanted"), ]
  get <- function(r, ro, f, col = "share") { x <- tot[tot$representation == r & tot$readout == ro & tot$factor == f, col]; if (length(x)) x[1] else NA }
  rows <- list()
  for (ro in unique(tot$readout[tot$representation == base])) for (alt in intersect(alts, unique(tot$representation))) {
    bu <- get(base, ro, "total_unwanted"); bw <- get(base, ro, "total_wanted")
    au <- get(alt, ro, "total_unwanted"); aw <- get(alt, ro, "total_wanted"); aw_hi <- get(alt, ro, "total_wanted", "hi")
    if (anyNA(c(bu, bw, au, aw))) next
    less_unwanted <- au < bu; keeps <- isTRUE(aw_hi >= bw)
    verdict <- if (less_unwanted && keeps) "removes unwanted, keeps wanted" else if (less_unwanted) "removes unwanted, loses wanted" else "no reduction of unwanted"
    rows[[length(rows) + 1L]] <- data.frame(readout = ro, base = base, alternative = alt, unwanted_base = bu, unwanted_alt = au,
                                            wanted_base = bw, wanted_alt = aw, verdict = verdict, stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

#' Question 1: per technical term and layer, how many contrasts exceed replicate noise (rule 1) and
#' how many of those flip responder calls (rule 2). `dec` is decisions.csv.
correction_summary <- function(dec) {
  dec$detected <- grepl("^detected", dec$verdict)
  dec$flipping <- dec$detected & dec$positivity_flips %in% seq_len(1e6)
  agg <- aggregate(cbind(contrasts = 1, detected = detected, flipping = flipping,
                         calls_flipped = ifelse(flipping, positivity_flips, 0)) ~ term + layer, dec, sum)
  agg$action <- ifelse(agg$flipping > 0, "correct (see the contrasts that flip calls)",
                       ifelse(agg$detected > 0, "keep as covariate, no correction", "no correction"))
  agg[order(agg$term, agg$layer), ]
}

#' Counting (binomial) noise of a frequency representation: the mean over tubes of the variance a frequency
#' has from counting k positive cells out of n by chance alone, on the representation's own scale (percent,
#' or arcsinh of percent by the delta method; stimulated minus unstimulated adds both tubes). 0 for the logit
#' count models (their tube term already excludes it) and for MFI (medians of thousands of cells).
counting_noise_var <- function(d, rep, cfg) {
  vb <- function(k, n) 1e4 * (k / n) * (1 - k / n) / n                     # variance of a percentage
  cf <- freq_asinh_cofactor(cfg); dv <- function(pct) 1 / (cf * sqrt(1 + (pct / cf)^2))   # d asinh(pct / cf) / d pct
  v <- switch(rep,
    freq_pct = vb(d$count, d$parent_count),
    freq_sub = vb(d$count, d$parent_count) + vb(d$neg_count, d$neg_parent_count),
    freq_asinh = dv(d$pct)^2 * vb(d$count, d$parent_count),
    freq_asinh_sub = dv(d$pct)^2 * vb(d$count, d$parent_count) + dv(d$neg_pct)^2 * vb(d$neg_count, d$neg_parent_count),
    0)
  mean(v, na.rm = TRUE)
}

#' Split the residual of a summarise_shares() table into counting noise (`cn_var`, absolute, on the response
#' scale) and the rest. The variance column must be on the response scale (scale_var applied).
split_counting_noise <- function(s, cn_var) {
  if (!is.finite(cn_var) || cn_var <= 0) return(s)
  tot <- sum(s$variance[!grepl("^total_", s$factor)])
  cn <- min(cn_var / tot, s$share[s$factor == "residual"])
  r <- s$factor == "residual"
  s$share[r] <- s$share[r] - cn; s$lo[r] <- pmax(s$lo[r] - cn, 0); s$hi[r] <- pmax(s$hi[r] - cn, 0); s$variance[r] <- s$variance[r] - cn * tot
  rbind(s, data.frame(factor = "counting_noise", group = "residual", share = cn, variance = cn * tot, lo = cn, hi = cn, stringsAsFactors = FALSE))
}
