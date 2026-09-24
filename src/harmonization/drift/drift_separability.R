# src/harmonization/drift/drift_separability.R — are the two clocks separable?
# Pre-flight check for a cohort design (docs/plans/drift_usability_plan.md,
# section 3.4): the biological clock (`week`, weeks since first vaccination) and
# the technical clock (`acq_day`) are confounded when every subject is acquired
# on one day per visit. Three tools, all working on tubes.csv alone:
#   confounding_index()   base R, instant
#   simulate_endpoint()   base R, one synthetic endpoint on the real design
#   design_power()        lme4 (requireNamespace), tens of glmer fits
#   separability_verdict() combines the two into one label

#' Within-unit share of the variance of `week` explained by `acq_day`.
#' Per patient, SS_p = sum (week - patient mean)^2 and SS_pd = sum over that
#' patient's days of (week - patient-day mean)^2; the index is
#' 1 - sum(SS_pd) / sum(SS_p), pooled over patients with more than one tube.
#' At 1 the day term is a relabelling of the week term for every subject and is
#' identifiable only through the unstimulated layer; the index drops below 1
#' only when one acquisition day carries tubes of more than one week.
#' Returns list(index, n_multi_week_days, n_subjects_multi_day,
#' n_subjects_multi_instrument). `index` is NA when no patient has more than
#' one week.
confounding_index <- function(tt) {
  per_pat <- split(tt, tt$patient)
  ss_p <- 0; ss_pd <- 0
  for (d in per_pat) {
    if (nrow(d) < 2L) next
    ss_p <- ss_p + sum((d$week - mean(d$week))^2)
    ss_pd <- ss_pd + sum(unlist(lapply(split(d$week, d$acq_day), function(w) (w - mean(w))^2)))
  }
  pd <- unique(tt[, c("patient", "acq_day", "week")])
  n_weeks_per_day <- table(paste(pd$patient, pd$acq_day))
  list(index = if (ss_p > 0) 1 - ss_pd / ss_p else NA_real_,
       n_multi_week_days = sum(n_weeks_per_day > 1L),
       n_subjects_multi_day = sum(vapply(per_pat, function(d) length(unique(d$acq_day)) > 1L, logical(1))),
       n_subjects_multi_instrument = sum(vapply(per_pat, function(d) length(unique(d$instrument)) > 1L, logical(1))))
}

#' Counts of ONE synthetic endpoint on the real design: the stimulated tubes of
#' `tt`, one row each. On the logit scale
#'   -6 + a_patient + week_slope * week_c + c_day[acq_day] + instr_shift * (instrument != ref) + N(0, 0.2)
#' with a_patient ~ N(0, 0.5), c_day ~ N(0, day_sd), week_c the centred week in
#' years (design_frame), parent_count ~ U(20000, 60000). `week_slope` is in
#' logit units per year. Returns the design frame (factors, week_c) plus
#' count and parent_count, ready for the lme4 endpoint formula.
simulate_endpoint <- function(tt, cfg, week_slope, day_sd, instr_shift, seed = 1, ref = cfg$reference_instrument) {
  set.seed(seed)
  s <- tt[!(tt$is_neg %in% TRUE), ]
  s <- design_frame(s, cfg, ref)
  pats <- levels(s[[cfg$design$unit]]); days <- sort(unique(s$acq_day))
  a_p <- setNames(rnorm(length(pats), 0, 0.5), pats)
  c_d <- setNames(rnorm(length(days), 0, day_sd), days)
  eta <- -6 + a_p[as.character(s[[cfg$design$unit]])] + week_slope * s$week_c + c_d[as.character(s$acq_day)] +
         instr_shift * (as.character(s$instrument) != ref) + rnorm(nrow(s), 0, 0.2)
  s$parent_count <- round(runif(nrow(s), 20000, 60000))
  s$count <- rbinom(nrow(s), s$parent_count, plogis(eta))
  rownames(s) <- NULL
  s
}

#' Fit the lme4 endpoint model to one simulated frame and read off the week
#' slope, the day SD and the instrument shift. The endpoint formula carries the
#' week only as a per-unit random slope, so the population slope `week_c` is
#' added as a fixed term. Returns list(week, day_sd, instr, converged), NA on
#' failure.
fit_simulated_endpoint <- function(d, cfg) {
  f <- update(design_formula(cfg, "endpoint", "lme4"), . ~ . + week_c)
  warned <- FALSE
  m <- tryCatch(withCallingHandlers(
         suppressMessages(lme4::glmer(f, data = d, family = binomial, control = lme4::glmerControl(optimizer = "bobyqa"))),
         warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") }),
       error = function(e) NULL)
  if (is.null(m)) return(list(week = NA_real_, day_sd = NA_real_, instr = NA_real_, converged = FALSE))
  vc <- as.data.frame(lme4::VarCorr(m))
  day_term <- setdiff(unlist(cfg$design$technical), "instrument")[1]
  dv <- vc$vcov[vc$grp == day_term & is.na(vc$var2) & vc$var1 %in% "(Intercept)"]
  fe <- lme4::fixef(m); ic <- grep("^instrument", names(fe))
  list(week = unname(fe["week_c"]), day_sd = sqrt(dv[1]),
       instr = if (length(ic)) unname(fe[ic[1]]) else NA_real_,
       converged = !warned && length(m@optinfo$conv$lme4$messages) == 0L)
}

#' Design power by simulation, `n_rep` repetitions of three scenarios:
#'   a  both effects (week_slope, day_sd, instr_shift): leakage = |bias| / truth
#'      for the week slope and the day SD
#'   b  week effect only (day_sd = 0): a false day detection is an estimated
#'      day SD above 0.1 logit (about a 10 percent multiplicative shift of the
#'      response; well below any drift worth correcting)
#'   c  day shift only (week_slope = 0): bias of the week slope in logit per year
#' `week_slope` is in logit per year (default 0.02 per week). Returns one row
#' per scenario x rep: scenario, rep, week_slope_true, week_slope_hat,
#' day_sd_true, day_sd_hat, instr_true, instr_hat, converged. Needs lme4.
design_power <- function(tt, cfg, n_rep = 10, week_slope = 0.02 * 52, day_sd = 0.5, instr_shift = 0.4, seed = 1) {
  if (!requireNamespace("lme4", quietly = TRUE)) stop("design_power() needs the lme4 package", call. = FALSE)
  scen <- data.frame(scenario = c("a_both", "b_week_only", "c_day_only"),
                     week = c(week_slope, week_slope, 0), day = c(day_sd, 0, day_sd), stringsAsFactors = FALSE)
  rows <- list()
  for (i in seq_len(nrow(scen))) for (r in seq_len(n_rep)) {
    d <- simulate_endpoint(tt, cfg, scen$week[i], scen$day[i], instr_shift, seed = seed + 1000L * i + r)
    fit <- fit_simulated_endpoint(d, cfg)
    rows[[length(rows) + 1L]] <- data.frame(scenario = scen$scenario[i], rep = r,
      week_slope_true = scen$week[i], week_slope_hat = fit$week, day_sd_true = scen$day[i], day_sd_hat = fit$day_sd,
      instr_true = instr_shift, instr_hat = fit$instr, converged = fit$converged, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Scenario summaries of a design_power() table (NA fits dropped):
#' leakage_week, leakage_day (scenario a), false_day_rate (b, day SD > 0.1),
#' week_bias (c, logit per year), n_fit (successful fits over all scenarios).
power_summary <- function(power) {
  a <- power[power$scenario == "a_both", ]; b <- power[power$scenario == "b_week_only", ]
  c <- power[power$scenario == "c_day_only", ]
  rel <- function(hat, truth) abs(mean(hat, na.rm = TRUE) - truth[1]) / truth[1]
  list(leakage_week = rel(a$week_slope_hat, a$week_slope_true),
       leakage_day = rel(a$day_sd_hat, a$day_sd_true),
       false_day_rate = mean(b$day_sd_hat > 0.1, na.rm = TRUE),
       week_bias = mean(c$week_slope_hat, na.rm = TRUE),
       n_fit = sum(!is.na(power$week_slope_hat)))
}

#' One label plus a one-paragraph explanation from the confounding index `ci`
#' and the design_power() table `power` (NULL when lme4 was not available):
#'   "not separable: add bridging acquisitions"          no subject on more than one day
#'   "separable"                                          leakage below 25 percent (a) and
#'                                                       false day detection below 20 percent (b)
#'   "separable only through the unstimulated layer"     otherwise, or index above 0.95
#' Returns list(label, explanation).
separability_verdict <- function(ci, power = NULL) {
  idx <- sprintf("confounding index %s, %d multi-week patient-day(s), %d subject(s) on more than one day, %d on more than one instrument",
                 format(round(ci$index, 3)), ci$n_multi_week_days, ci$n_subjects_multi_day, ci$n_subjects_multi_instrument)
  if (ci$n_subjects_multi_day == 0L)
    return(list(label = "not separable: add bridging acquisitions",
      explanation = paste0("No subject was acquired on more than one day, so a day effect cannot be told from a subject effect ",
                           "in any layer. Re-acquire a few subjects on other days (bridging acquisitions) before fitting. (", idx, ")")))
  if (is.null(power)) {
    lab <- if (!is.na(ci$index) && ci$index > 0.95) "separable only through the unstimulated layer" else "separable"
    return(list(label = lab, explanation = paste0("Judged from the confounding index alone (lme4 not available, no simulation). ", idx, ".")))
  }
  ps <- power_summary(power)
  sim <- sprintf("Simulation (%d fits): week-slope leakage %.0f percent, day-SD leakage %.0f percent, false day detection %.0f percent, week bias under a pure day shift %.2f logit/year.",
                 ps$n_fit, 100 * ps$leakage_week, 100 * ps$leakage_day, 100 * ps$false_day_rate, ps$week_bias)
  ok <- isTRUE(ps$leakage_week < 0.25) && isTRUE(ps$leakage_day < 0.25) && isTRUE(ps$false_day_rate < 0.20)
  # The simulation is the direct test on this design and outranks the index: a high index says the
  # clocks are nearly collinear within subjects, but a few multi-week days can still identify both.
  if (ok)
    return(list(label = "separable",
      explanation = paste0("The stimulated tubes alone recover a week slope and a day SD injected together, and a pure week effect ",
                           "is not mistaken for day drift. ", sim, " (", idx,
                           if (isTRUE(ci$index > 0.95)) "; the high index means this identification rests on the few multi-week acquisition days)" else ")")))
  why <- "The simulation did not recover both clocks from the stimulated tubes"
  list(label = "separable only through the unstimulated layer",
       explanation = paste0(why, ", so the day term is identified mainly by the unstimulated tubes (layers A and BG), ",
                            "where no biology is expected; endpoint estimates of drift lean on that assumption. ", sim, " (", idx, ")"))
}
