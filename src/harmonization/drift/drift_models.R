# src/harmonization/drift/drift_models.R — model frames and posterior handling
# Base R. Builds the model frames (layer A: neg summaries; layer BG: cytokine
# background in negs; layer B: gated counts) from the step 4/5 tables, assembles
# the model formulas from the `design` block of the config, and turns posterior
# draws of the day and instrument effects into the decision table. brms, Stan,
# lme4 and BayesFlow all consume/produce the same shapes, so they are comparable
# in step 9.

#' Design block used when a config predates the `design` field: the MTB layout
#' without protected covariates or a time slope. `time_df` adds that many non-linear time
#' columns (natural spline on log(1 + week)) next to the linear slope; 1 or less keeps the line
#' only. Default 5: on a synthetic cohort with the MTB design, a response peak around week 12
#' and no day drift, the straight line reported day variance of 0.13 to 0.94 logit^2 and an
#' inflated operator effect, the log-time spline 0.002 to 0.018 (2026-10-12, 9 fits each).
default_design <- function() {
  list(unit = "patient", technical = c("acq_day", "instrument"),
       protected_fixed = character(0), time = NULL, stimulation = "cocktail", time_df = 5)
}

#' Names of the non-linear time columns design_frame() adds (character(0) when
#' the design has no time column or time_df <= 1).
time_nl_terms <- function(cfg) {
  ds <- cfg$design; k <- if (is.null(ds$time_df)) 1L else as.integer(ds$time_df)
  if (is.null(ds$time) || k <= 1L) character(0) else paste0("week_nl", seq_len(k))
}

#' The acquisition-day term: the first technical term that is not `instrument`.
#' It is the random intercept whose levels are the day effects. NA when there is none.
design_day_term <- function(cfg) setdiff(unlist(cfg$design$technical), "instrument")[1]

#' The other technical terms (instrument and any tech_* column of the sample sheet,
#' e.g. operator): fixed factors, each contrasted against its reference level and
#' judged by the decision rule like the instrument.
design_technical_factors <- function(cfg) setdiff(unlist(cfg$design$technical), design_day_term(cfg))

#' Design columns of `cfg` (unlisted, in a fixed order) that a frame must carry.
design_columns <- function(cfg) {
  ds <- cfg$design
  unique(c(ds$unit, unlist(ds$technical), unlist(ds$protected_fixed), ds$time, ds$stimulation))
}

#' Turn the design columns of `df` into model-ready factors: unit, the technical
#' terms (instrument relevelled to `ref` when given, any other technical factor to
#' `design$reference_levels[[name]]` when given), every protected covariate,
#' and the stimulation term when present. Adds `week_c`, the time column centred
#' on the frame mean and in years, when the design names a time column.
design_frame <- function(df, cfg, ref = NULL) {
  ds <- cfg$design
  as_factor <- c(ds$unit, unlist(ds$technical), unlist(ds$protected_fixed))
  if (!is.null(ds$stimulation) && ds$stimulation %in% names(df)) as_factor <- c(as_factor, ds$stimulation)
  for (v in as_factor) {
    if (!v %in% names(df)) stop("design column `", v, "` is missing from the frame", call. = FALSE)
    df[[v]] <- factor(df[[v]])
  }
  if (!is.null(ref) && "instrument" %in% names(df) && ref %in% levels(df$instrument))
    df$instrument <- relevel(df$instrument, ref = ref)
  for (v in intersect(names(ds$reference_levels), setdiff(names(df), "instrument")))
    if (ds$reference_levels[[v]] %in% levels(df[[v]])) df[[v]] <- relevel(df[[v]], ref = ds$reference_levels[[v]])
  if (!is.null(ds$time)) df$week_c <- (df[[ds$time]] - mean(df[[ds$time]])) / 52   # keeps slopes O(1)
  nl <- time_nl_terms(cfg)
  if (length(nl)) {
    # natural spline of log(1 + week), made orthogonal to the intercept and the line, standardised:
    # the curvature (rise after a boost, waning) that a straight line would push into the day term.
    # Log time because visits crowd the first weeks and thin out over years.
    w <- df[[ds$time]]; ok <- is.finite(w)
    for (j in seq_along(nl)) df[[nl[j]]] <- 0
    if (length(unique(w[ok])) > length(nl) + 1L) {
      r <- qr.resid(qr(cbind(1, df$week_c[ok])), splines::ns(log1p(pmax(w[ok], 0)), df = length(nl)))
      for (j in seq_along(nl)) if (stats::var(r[, j]) > 1e-10) { df[[nl[j]]][ok] <- as.numeric(scale(r[, j])) }
      df[[nl[1]]][!ok] <- NA_real_
    }
  }
  df
}

#' Model formula for one response class, assembled from cfg$design.
#'   technical            (layer A):  instrument fixed + (1 | unit) + (1 | day terms); no biology
#'   background_counts    (layer BG): + protected_fixed, (1 + week_c | unit); beta-binomial counts
#'   background_intensity (layer BG): + protected_fixed + week_c fixed; gaussian
#'   endpoint             (layer B):  as background_counts + (1 | stimulation)
#' The day term (design_day_term()) is a random intercept; every other technical
#' term (instrument, operator, ...) a fixed factor with K-1 contrasts against the
#' reference level of the frame. lme4 count models get the cbind() response and an
#' observation-level (1 | filename) term for overdispersion.
design_formula <- function(cfg, class = c("technical", "background_counts", "background_intensity", "endpoint"),
                           engine = c("brms", "lme4"), response_lhs = NULL) {
  class <- match.arg(class); engine <- match.arg(engine)
  ds <- cfg$design
  counts <- class %in% c("background_counts", "endpoint")
  if (is.null(response_lhs))
    response_lhs <- if (!counts) "y" else if (engine == "brms") "count | trials(parent_count)" else "cbind(count, parent_count - count)"
  fixed <- c("1", design_technical_factors(cfg))
  unit_re <- sprintf("(1 | %s)", ds$unit)
  if (class != "technical") {
    fixed <- c(fixed, unlist(ds$protected_fixed))
    if (!is.null(ds$time)) {
      # population week slope as a fixed term (plus its spline curvature); the count classes also
      # get a per-unit slope around it
      fixed <- c(fixed, "week_c", time_nl_terms(cfg))
      if (counts) unit_re <- sprintf("(1 + week_c | %s)", ds$unit)
    }
  }
  day <- design_day_term(cfg)
  terms <- c(fixed, unit_re, if (!is.na(day)) sprintf("(1 | %s)", day))
  if (class == "endpoint" && !is.null(ds$stimulation)) terms <- c(terms, sprintf("(1 | %s)", ds$stimulation))
  if (counts && engine == "lme4") terms <- c(terms, "(1 | filename)")
  as.formula(paste(response_lhs, "~", paste(terms, collapse = " + ")), env = globalenv())
}

#' Population-level coefficient names of the fixed design terms `vars` in
#' `frame`, as brms and model.matrix() spell them: data.frame(coef, term, level).
#' Factors contribute one row per non-reference level, numeric terms one row.
fixed_terms_table <- function(frame, vars) {
  rows <- lapply(unlist(vars), function(v) {
    x <- frame[[v]]
    if (is.factor(x)) {
      lv <- levels(droplevels(x))[-1]
      data.frame(coef = paste0(v, lv), term = v, level = lv, stringsAsFactors = FALSE)
    } else data.frame(coef = v, term = v, level = NA_character_, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame(coef = character(0), term = character(0), level = character(0)) else out
}

#' Fixed design terms of one class, as they appear in the model formula:
#' the technical factors, the protected covariates and week_c. Order matches
#' design_formula().
design_fixed_vars <- function(cfg, class) {
  ds <- cfg$design
  vars <- design_technical_factors(cfg)
  if (class != "technical") {
    vars <- c(vars, unlist(ds$protected_fixed))
    if (!is.null(ds$time)) vars <- c(vars, "week_c", time_nl_terms(cfg))
  }
  vars
}

#' Layer A frame: one row per neg tube x marker x statistic.
#' Statistics used: median, spread (q90-q10), valley and frac_pos (bimodal
#' markers only), emd. `summaries`, `emd`, `tubes` are the CSVs of steps 4 and 1.
#' `cfg` is only needed for the design block; without it the MTB default is used.
layer_a_frame <- function(summaries, emd, tubes, technical_markers, cfg = NULL, ref = NULL) {
  s <- summaries[summaries$marker %in% technical_markers, ]
  long <- rbind(
    data.frame(filename = s$filename, marker = s$marker, statistic = "median", value = s$median),
    data.frame(filename = s$filename, marker = s$marker, statistic = "spread", value = s$spread))
  if ("valley" %in% names(s)) {
    v <- s[!is.na(s$valley), ]
    long <- rbind(long,
      data.frame(filename = v$filename, marker = v$marker, statistic = "valley", value = v$valley),
      data.frame(filename = v$filename, marker = v$marker, statistic = "frac_pos", value = v$frac_pos))
  }
  e <- emd[emd$marker %in% technical_markers, ]
  long <- rbind(long, data.frame(filename = e$filename, marker = e$marker, statistic = "emd", value = e$emd))
  fr <- merge_tubes(long, tubes, cfg)
  fr <- design_frame(fr, frame_cfg(cfg, tubes), ref)
  fr[order(fr$marker, fr$statistic, fr$patient, fr$acq_day), ]
}

#' Config used by the frame builders when none is passed: the default design
#' plus prior_vacc and week when tubes.csv carries them (today's MTB frames).
frame_cfg <- function(cfg, tubes) {
  if (!is.null(cfg$design)) return(cfg)
  ds <- default_design()
  if ("prior_vacc" %in% names(tubes)) ds$protected_fixed <- "prior_vacc"
  if ("week" %in% names(tubes)) ds$time <- "week"
  list(design = ds)
}

#' Join the tube metadata (identity, biology and design columns) onto `long`.
merge_tubes <- function(long, tubes, cfg = NULL) {
  cols <- unique(c("filename", "patient", "week", "is_neg", "acq_day", "instrument",
                   "patient_day", "cocktail", if (!is.null(cfg$design)) design_columns(cfg)))
  merge(long, tubes[, intersect(cols, names(tubes))], by = "filename")
}

#' Add the design columns of `cfg` that `df` lacks from tubes.csv in `out` (or from the
#' tube table `tt` when given), by filename (a step 5 table written before a technical
#' factor was added lacks it).
with_tube_columns <- function(df, out, cfg, tt = NULL) {
  miss <- setdiff(design_columns(cfg), names(df))
  if (!length(miss)) return(df)
  if (is.null(tt)) {
    tp <- file.path(out, "tubes.csv")
    if (!file.exists(tp)) return(df)
    tt <- read.csv(tp, stringsAsFactors = FALSE)
  }
  merge(df, tt[, c("filename", intersect(miss, names(tt))), drop = FALSE], by = "filename")
}

#' Layer B frame: cocktail tubes only, one row per tube x parent x gate.
layer_b_frame <- function(freqs, cfg = NULL, ref = NULL) {
  fb <- freqs[!freqs$is_neg & !is.na(freqs$neg_pct), ]
  fb <- design_frame(fb, frame_cfg(cfg, freqs), ref)
  fb$response <- paste(fb$parent, fb$gate, sep = ":")
  fb
}

#' The frames from the step 4/5 tables in `out`. With `ref` the instrument
#' factor is relevelled so that the fixed effect is the contrast against it.
read_layer_a <- function(out, cfg, tubes, ref = NULL) {
  layer_a_frame(read.csv(file.path(out, "neg_summaries.csv")),
                read.csv(file.path(out, "neg_emd.csv")), tubes, cfg$technical_markers, cfg, ref)
}
read_layer_b <- function(out, ref = NULL, cfg = NULL) {
  layer_b_frame(with_tube_columns(read.csv(file.path(out, "gated_frequencies.csv")), out, cfg), cfg, ref)
}

#' Layer BG counts: cytokine-positive counts in the neg tubes, one row per
#' neg tube x parent x gate (parent gate non-empty). Response key parent:gate.
read_layer_bg_counts <- function(out, cfg, ref = NULL) {
  fq <- with_tube_columns(read.csv(file.path(out, "gated_frequencies.csv"), stringsAsFactors = FALSE), out, cfg)
  bg <- fq[fq$is_neg & !is.na(fq$parent_count) & fq$parent_count > 0, ]
  cols <- unique(c("filename", "patient", "week", "acq_day", "instrument", "patient_day",
                   "parent", "gate", "count", "parent_count", "pct", design_columns(cfg)))
  bg <- bg[, intersect(cols, names(bg))]
  bg$response <- paste(bg$parent, bg$gate, sep = ":")
  design_frame(bg, cfg, ref)
}

#' Layer BG intensity: median of every cytokine marker in the neg tubes, in the
#' layer A column layout. Response key marker:median.
read_layer_bg_intensity <- function(out, cfg, tubes, ref = NULL) {
  s <- read.csv(file.path(out, "neg_summaries.csv"), stringsAsFactors = FALSE)
  s <- s[s$marker %in% cfg$cytokine_markers, ]
  long <- data.frame(filename = s$filename, marker = s$marker, statistic = "median", value = s$median,
                     stringsAsFactors = FALSE)
  fr <- merge_tubes(long, tubes, cfg)
  fr <- fr[fr$is_neg, ]
  fr$response <- paste(fr$marker, fr$statistic, sep = ":")
  design_frame(fr, cfg, ref)
}

#' Replicate SD per marker x statistic from ALL tubes of the same sample (patient,
#' day and visit; technical markers should not depend on the peptide). Needs step 4 run with
#' --all_tubes; with negs only, returns NA and step 9 falls back to the
#' posterior residual SD.
layer_a_replicate_sd <- function(frame_all_tubes) {
  sp <- split(frame_all_tubes, list(frame_all_tubes$marker, frame_all_tubes$statistic), drop = TRUE)
  do.call(rbind, lapply(sp, function(d) data.frame(
    marker = d$marker[1], statistic = d$statistic[1],
    replicate_sd = replicate_sd(d$value, same_sample(d)), stringsAsFactors = FALSE)))
}

#' Same-sample key: tubes of one patient, acquisition day and visit (week). On a day that
#' carries several visits of a patient, patient_day alone would mix samples.
same_sample <- function(d) if ("week" %in% names(d)) paste(d$patient_day, d$week) else d$patient_day

#' Count positivity calls that change when an effect on the logit scale is
#' removed from the tubes it applies to. `fb` is a layer B frame for ONE
#' response; `effect` is a scalar (posterior mean) applied to rows in `which`.
count_positivity_flips <- function(fb, effect, which, min_subtracted_pct, min_fold_over_neg) {
  before <- positivity_call(fb$pct, fb$neg_pct, min_subtracted_pct, min_fold_over_neg)$positive
  p <- fb$pct / 100; p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  adj <- p
  adj[which] <- plogis(qlogis(p[which]) - effect)
  # the neg tube of the same day/instrument is shifted the same way
  pn <- fb$neg_pct / 100; pn <- pmin(pmax(pn, 1e-6), 1 - 1e-6)
  adjn <- pn; adjn[which] <- plogis(qlogis(pn[which]) - effect)
  after <- positivity_call(100 * adj, 100 * adjn, min_subtracted_pct, min_fold_over_neg)$positive
  sum(before != after, na.rm = TRUE)
}

#' Standard shape for posterior draws of the effects, shared by every fitter:
#' data.frame(model, layer, response, term, level, draw, value)
#' term is "instrument" (level = non-reference instrument), "day" (level = acq_day) or
#' the name of another technical factor (level = its non-reference level).
#' `instrument_draws` is a matrix or named list with one column per non-reference
#' instrument (column names are the levels), or a single vector whose level is
#' `instrument_level`. `technical_draws` is a named list (term -> matrix with one
#' column per non-reference level) for the other technical factors.
effects_long <- function(model, layer, response, instrument_draws, day_draws, instrument_level = NA,
                         technical_draws = list()) {
  if (is.list(instrument_draws) && !is.data.frame(instrument_draws)) instrument_draws <- do.call(cbind, instrument_draws)
  if (is.null(dim(instrument_draws))) {
    instrument_draws <- matrix(instrument_draws, ncol = 1L, dimnames = list(NULL, instrument_level))
  }
  one <- function(term, level, v) data.frame(model = model, layer = layer, response = response, term = term,
                                            level = level, draw = seq_along(v), value = as.numeric(v),
                                            stringsAsFactors = FALSE)
  out <- lapply(seq_len(ncol(instrument_draws)), function(j)
    one("instrument", if (is.null(colnames(instrument_draws))) NA else colnames(instrument_draws)[j], instrument_draws[, j]))
  for (tm in names(technical_draws)) { m <- technical_draws[[tm]]
    out <- c(out, lapply(seq_len(ncol(m)), function(j) one(tm, colnames(m)[j], m[, j]))) }
  out <- c(out, lapply(colnames(day_draws), function(d) one("day", d, day_draws[, d])))
  do.call(rbind, out)
}

#' Posterior summaries of the protected covariates and time slope, one row per
#' coefficient: data.frame(model, layer, response, covariate, level, mean, lo, hi).
#' `draws_by_covariate` is a named list; each element a vector (level NA) or a
#' matrix with one column per factor level. These are protected biology: they
#' are reported, never treated as drift.
covariates_long <- function(model, layer, response, draws_by_covariate) {
  rows <- list()
  for (cv in names(draws_by_covariate)) {
    m <- draws_by_covariate[[cv]]
    if (is.null(dim(m))) m <- matrix(m, ncol = 1L, dimnames = list(NULL, NA_character_))
    for (j in seq_len(ncol(m))) {
      v <- m[, j]
      rows[[length(rows) + 1L]] <- data.frame(model = model, layer = layer, response = response, covariate = cv,
        level = colnames(m)[j], mean = mean(v), lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(model = character(0), layer = character(0), response = character(0),
    covariate = character(0), level = character(0), mean = numeric(0), lo = numeric(0), hi = numeric(0)))
  do.call(rbind, rows)
}

#' Split a draws matrix whose columns are named by `ft$coef` (a fixed_terms_table)
#' into the instrument matrix (columns = levels), the other `technical` factors
#' (named list of such matrices, for effects_long()) and the covariate list that
#' covariates_long() takes. Columns absent from `draws` are skipped.
split_fixed_draws <- function(draws, ft, technical = "instrument") {
  draws <- as.matrix(draws)
  ft <- ft[ft$coef %in% colnames(draws), ]
  by_term <- function(tm) { r <- ft[ft$term == tm, ]; m <- draws[, r$coef, drop = FALSE]; colnames(m) <- r$level; m }
  inst <- by_term("instrument")
  tech <- lapply(setNames(nm = setdiff(intersect(technical, ft$term), "instrument")), by_term)
  covs <- list()
  for (tm in unique(ft$term[!ft$term %in% c("instrument", technical)])) {
    r <- ft[ft$term == tm, ]
    m <- draws[, r$coef, drop = FALSE]; colnames(m) <- r$level
    covs[[tm]] <- if (all(is.na(r$level))) as.numeric(m) else m
  }
  list(instrument = inst, technical = tech, covariates = covs)
}

#' All effects_*.csv written by the step 6 fitters, or NULL when there are none.
read_effects <- function(out) {
  files <- list.files(out, "^effects_.*\\.csv$", full.names = TRUE)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

#' Variance components and the instrument fixed effect of one lme4 fit as a
#' summary row (one per non-reference instrument level). `resid_grp` is the
#' VarCorr group holding the residual variance ("Residual" for lmer, the
#' observation-level term for glmer). `unit` and `day` name the grouping terms.
lme_variance_row <- function(m, layer, response, n, resid_grp, unit = "patient", day = "acq_day") {
  vc <- as.data.frame(lme4::VarCorr(m))
  vc <- vc[is.na(vc$var2) & (is.na(vc$var1) | vc$var1 == "(Intercept)"), ]   # intercept variances only
  v <- setNames(vc$vcov, vc$grp)
  fe <- lme4::fixef(m); se <- sqrt(diag(as.matrix(vcov(m))))
  ic <- grep("^instrument", names(fe))
  data.frame(layer = layer, response = response, n = n,
    var_patient = v[unit], var_day = v[day], var_resid = v[resid_grp],
    icc_day = v[day] / sum(v), drift_to_biology = sqrt(v[day] / v[unit]),
    instrument_effect = if (length(ic)) fe[ic] else NA, instrument_se = if (length(ic)) se[ic] else NA,
    instrument_level = if (length(ic)) sub("^instrument", "", names(fe)[ic]) else NA,
    stringsAsFactors = FALSE, row.names = NULL)
}

#' Number of divergent transitions of a brms fit.
n_divergent <- function(fit) {
  np <- brms::nuts_params(fit)
  sum(np$Value[np$Parameter == "divergent__"])
}

#' Summary row of one brms fit: group SDs and sigma (multiplied by `scale` to
#' undo a standardised response; sigma is NA for families without one) plus
#' convergence diagnostics. `unit` and `day` name the grouping terms.
brms_fit_row <- function(fit, layer, response, n, scale = 1, unit = "patient", day = "acq_day") {
  vc <- brms::VarCorr(fit)
  s <- as.data.frame(vc[[day]]$sd); sp <- as.data.frame(vc[[unit]]$sd)
  spec <- summary(fit)$spec_pars
  data.frame(layer = layer, response = response, n = n,
    tau_day = s$Estimate * scale, tau_patient = sp$Estimate[1] * scale,
    sigma = if ("sigma" %in% rownames(spec)) spec["sigma", "Estimate"] * scale else NA,
    max_rhat = max(brms::rhat(fit), na.rm = TRUE), n_divergent = n_divergent(fit),
    stringsAsFactors = FALSE)
}

#' Agreement of two fitters on the same effects: does the posterior mean of B
#' fall inside the 95 percent interval of A, per response x term x level.
compare_effects <- function(eff, model_a, model_b) {
  key <- c("layer", "response", "term", "level")
  agg <- function(m) {
    d <- eff[eff$model == m, ]
    sp <- split(d, d[key], drop = TRUE)
    do.call(rbind, lapply(sp, function(x) cbind(x[1, key], mean = mean(x$value),
      lo = unname(quantile(x$value, 0.025)), hi = unname(quantile(x$value, 0.975)))))
  }
  a <- agg(model_a); b <- agg(model_b)
  names(b)[5:7] <- c("mean_b", "lo_b", "hi_b")
  m <- merge(a, b, by = key)
  m$model_a <- model_a; m$model_b <- model_b
  m$agree <- m$mean_b >= m$lo & m$mean_b <= m$hi
  m
}

#' Agreement counts "agreed/total" per model x layer x term.
agreement_table <- function(cmp) {
  aggregate(agree ~ model_b + layer + term, cmp, function(x) sprintf("%d/%d", sum(x), length(x)))
}
