# src/harmonization/drift/drift_biology.R — the biological week slope and its sensitivity
# Base R. Formulas for the endpoint model with and without the technical terms, and the
# assembly of the week-slope tables written by scripts/drift/06_biology.R and read by the report.
# design_formula() keeps week_c in the per-patient random slope only; the population slope is
# added here as a fixed effect so the biological effect can be reported.

#' `cfg` with every technical term removed: no instrument contrast, no day random intercept.
biology_cfg <- function(cfg) { cfg$design$technical <- character(0); cfg }

#' lme4 endpoint formula with the fixed week slope. `technical = FALSE` drops the technical terms.
biology_formula <- function(cfg, technical = TRUE) {
  if (is.null(cfg$design$time)) stop("the design has no `time` field: no week slope to estimate", call. = FALSE)
  design_formula(if (technical) cfg else biology_cfg(cfg), "endpoint", "lme4")
}

#' Fixed week slope of an lme4 fit as c(est, se); NA for a failed fit (NULL).
slope_of <- function(fit) {
  if (is.null(fit)) return(c(NA_real_, NA_real_))
  fe <- lme4::fixef(fit); se <- sqrt(diag(as.matrix(vcov(fit))))
  c(unname(fe["week_c"]), unname(se[names(fe) == "week_c"]))
}

#' week_slope_sensitivity.csv from records list(response, model, slope = c(est, se), n).
sensitivity_table <- function(recs) {
  rows <- lapply(recs, function(r) data.frame(response = r$response, model = r$model, week_slope = r$slope[1],
                                              se = r$slope[2], n = r$n, stringsAsFactors = FALSE))
  if (length(rows)) do.call(rbind, rows)
  else data.frame(response = character(0), model = character(0), week_slope = numeric(0), se = numeric(0), n = integer(0))
}

#' week_slope_lodo.csv from records list(response, dropped_day, slope = c(est, se), n_tubes_dropped).
lodo_table <- function(recs) {
  rows <- lapply(recs, function(r) data.frame(response = r$response, dropped_day = r$dropped_day, week_slope = r$slope[1],
                                              se = r$slope[2], n_tubes_dropped = r$n_tubes_dropped, stringsAsFactors = FALSE))
  if (length(rows)) do.call(rbind, rows)
  else data.frame(response = character(0), dropped_day = character(0), week_slope = numeric(0), se = numeric(0), n_tubes_dropped = integer(0))
}

#' Per-patient slopes (fixed + random) from coef(fit)[[unit]]: a data frame with the units as row names
#' and a `week_c` column.
subject_slopes <- function(unit_coef, response) {
  data.frame(response = response, patient = rownames(unit_coef), slope = unit_coef$week_c, stringsAsFactors = FALSE, row.names = NULL)
}

#' Full and no_technical slopes side by side, with the change in units of the full SE and relative
#' to the full slope.
sensitivity_pairs <- function(sens) {
  full <- sens[sens$model == "full", ]; nt <- sens[sens$model == "no_technical", ]
  p <- merge(full[, c("response", "week_slope", "se")], nt[, c("response", "week_slope", "se")], by = "response", suffixes = c("", "_nt"))
  p$shift_se <- (p$week_slope_nt - p$week_slope) / p$se
  p$rel_change <- (p$week_slope_nt - p$week_slope) / abs(p$week_slope)
  p
}

#' One sentence: the endpoint whose slope moves most (in units of the full-model SE) without the
#' technical terms. The relative change is quoted only when the full slope exceeds its SE.
sensitivity_sentence <- function(sens) {
  p <- sensitivity_pairs(sens); p <- p[is.finite(p$shift_se), ]
  if (!nrow(p)) return("No endpoint has both a full and a no_technical slope.")
  i <- which.max(abs(p$shift_se))
  rel <- if (abs(p$week_slope[i]) > p$se[i]) sprintf("%.0f percent of the full slope, %.1f SE", 100 * abs(p$rel_change[i]), abs(p$shift_se[i]))
         else sprintf("%.1f SE; the full slope is within one SE of zero, so the relative change is not meaningful", abs(p$shift_se[i]))
  sprintf("Removing the technical terms changes the week slope most for %s: from %.2f to %.2f logit per year (%s); every other endpoint moves less.",
          p$response[i], p$week_slope[i], p$week_slope_nt[i], rel)
}

#' Leave-one-day-out shifts in units of the full-model SE; `flag` marks days beyond one SE.
lodo_shifts <- function(lodo, sens) {
  full <- sens[sens$model == "full", c("response", "week_slope", "se")]
  m <- merge(lodo, full, by = "response", suffixes = c("", "_full"))
  m$shift_se <- (m$week_slope - m$week_slope_full) / m$se_full
  m$flag <- is.finite(m$shift_se) & abs(m$shift_se) > 1
  m
}

#' One sentence: the day whose removal moves a slope most.
lodo_sentence <- function(lodo, sens) {
  m <- lodo_shifts(lodo, sens); m <- m[is.finite(m$shift_se), ]
  if (!nrow(m)) return("No leave-one-day-out refit converged.")
  i <- which.max(abs(m$shift_se))
  sprintf("The day whose removal moves a slope most is %s (%d tubes) for %s: from %.2f to %.2f logit per year (%.1f SE); %d of %d refits move a slope by more than one SE.",
          m$dropped_day[i], m$n_tubes_dropped[i], m$response[i], m$week_slope_full[i], m$week_slope[i], abs(m$shift_se[i]), sum(m$flag), nrow(m))
}
