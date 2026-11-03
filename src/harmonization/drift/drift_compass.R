# src/harmonization/drift/drift_compass.R — Boolean cytokine subsets, COMPASS inputs and the subset model
# Base R. Step 16 (scripts/drift/16_compass.R, docs/plans/compass_plan_2026-10-14.md):
#   - per-cell Boolean coding of the cytokine gates into the 2^M subsets, named as COMPASS names them
#     ("IFNg&!IL2&TNFa&!CD107a"), the all-negative subset last;
#   - the COMPASS category filter and the rebuilt all-negative column;
#   - functionality and polyfunctionality scores per posterior draw (Lin et al. 2015);
#   - the long frame of the two-sided subset model and the hypothesis summaries H1 to H3;
#   - a simulator of Boolean counts on the real design, for the synthetic checks.
# flowWorkspace, COMPASS and brms are only used by the script, so everything here is testable without them.

#' Defaults of the `compass:` config block.
default_compass_cfg <- function(cfg) {
  list(parents = names(cfg$parent_populations), markers = unlist(cfg$cytokine_markers),
       category_filter = list(min_count = 5, min_tubes = 2), low_info_parent = 1000,
       iterations = 10000, replications = 8, seeds = c(100, 200, 300),
       visit_column = "visit", dose_regex = "^V(\\d+)", reference_dose = 1, late_doses = c(3, 4),
       weeks_df = 3, chains = 4, iter = 2000, warmup = 1000, threads = 2, adapt_delta = 0.9, seed = 20261014)
}

compass_cfg <- function(cfg) {
  d <- default_compass_cfg(cfg)
  if (is.null(cfg$compass)) d else modifyList(d, cfg$compass)
}

#' FlowJo gate leaf names -> canonical cytokine names of the config. The same mapping as step 5
#' (05_gated_frequencies.R, canon_gate): the first non-channel marker pattern that matches wins;
#' unmatched names are returned unchanged.
canon_cytokine_gate <- function(g, cfg) {
  out <- g
  for (cy in cfg$cytokine_markers) for (pat in cfg$markers[[cy]]) {
    if (grepl("^\\^", pat)) next   # channel-name patterns do not apply to gate names
    out[grepl(pat, g, ignore.case = TRUE, perl = TRUE) & out == g] <- cy
  }
  out
}

#' All 2^M Boolean subsets of `markers`, COMPASS-named, ordered by degree (number of positive
#' markers) from high to low, ties in bit order; the all-negative subset is last.
compass_subset_names <- function(markers) {
  M <- length(markers)
  bits <- as.matrix(expand.grid(rep(list(0:1), M)))   # row r = code r - 1, column j = bit of marker j
  nm <- apply(bits, 1, function(b) paste(ifelse(b == 1, markers, paste0("!", markers)), collapse = "&"))
  deg <- rowSums(bits)
  nm[order(-deg, seq_along(nm))]
}

#' Number of positive markers in each COMPASS subset name.
subset_degree <- function(subsets) {
  vapply(strsplit(subsets, "&", fixed = TRUE), function(p) sum(!startsWith(p, "!")), integer(1))
}

#' Name of the all-negative subset of `markers`.
null_subset <- function(markers) paste(paste0("!", markers), collapse = "&")

#' Boolean subset counts of one parent population. `parent` is a logical vector over all events
#' (membership of the parent gate), `hits` a named list of logical vectors of the same length
#' (membership of each marker gate; names are the markers). Events outside the parent are ignored,
#' so a marker gate need not be nested. Returns integer counts named by compass_subset_names().
boolean_counts <- function(parent, hits) {
  markers <- names(hits)
  if (is.null(markers) || any(!nzchar(markers))) stop("`hits` must be a named list", call. = FALSE)
  if (any(lengths(hits) != length(parent))) stop("every marker vector must have the length of `parent`", call. = FALSE)
  code <- integer(sum(parent))
  for (j in seq_along(hits)) code <- code + as.integer(hits[[j]][parent]) * 2L^(j - 1L)
  tab <- tabulate(code + 1L, nbins = 2L^length(markers))
  bits <- as.matrix(expand.grid(rep(list(0:1), length(markers))))
  names(tab) <- apply(bits, 1, function(b) paste(ifelse(b == 1, markers, paste0("!", markers)), collapse = "&"))
  tab[compass_subset_names(markers)]
}

#' Marginal count of each marker from Boolean counts (named by compass_subset_names(markers)).
boolean_marginals <- function(counts, markers) {
  vapply(markers, function(m) {
    pos <- !grepl(paste0("(^|&)!", m, "(&|$)"), names(counts))
    sum(counts[pos])
  }, numeric(1))
}

#' COMPASS category filter (COMPASS(): category_filter = colSums(x > 5) > 2) on a stimulated count
#' matrix (tubes x subsets). Returns the names of the non-null subsets kept.
compass_filter <- function(n_s, null_name, min_count = 5, min_tubes = 2) {
  cand <- setdiff(colnames(n_s), null_name)
  cand[colSums(n_s[, cand, drop = FALSE] > min_count) > min_tubes]
}

#' Count matrix restricted to the subsets `keep`, with the all-negative column rebuilt as
#' parent_count minus the kept subsets (cells of dropped subsets go to the null), placed last as
#' COMPASS needs it. `m` has one column per subset; `parent_count` one value per row.
collapse_to_kept <- function(m, keep, null_name, parent_count) {
  out <- m[, keep, drop = FALSE]
  null <- parent_count - rowSums(out)
  if (any(null < 0)) stop("kept subsets exceed the parent count", call. = FALSE)
  out <- cbind(out, null); colnames(out)[ncol(out)] <- null_name
  out
}

#' Functionality (FS) and polyfunctionality (PFS) score of every individual and draw, from COMPASS
#' gamma draws (individuals x categories x draws, the null category last), as COMPASS's
#' FunctionalityScore() and PolyfunctionalityScore() compute them from mean_gamma, with `n` the
#' number of markers: FS = sum(gamma[non-null]) / (2^n - 1); PFS = sum(gamma * d / choose(n, d)) / n * 2 / (n + 1).
compass_score_draws <- function(gamma, degree, n) {
  K <- dim(gamma)[2]
  w_pfs <- degree / choose(n, degree) / n * 2 / (n + 1)
  fs <- apply(gamma[, -K, , drop = FALSE], c(1, 3), sum) / (2^n - 1)
  pfs <- apply(gamma, c(1, 3), function(g) sum(g * w_pfs))
  list(FS = fs, PFS = pfs)   # individuals x draws
}

#' Mean of `mean_gamma` over several COMPASS fits of the same data with different seeds, the
#' "mean model" of Lin et al. 2015. COMPASS::COMPASSMCMCDiagnosis() also advertises one, but its
#' mean_model returns the FIRST fit unchanged: it calls abind() through Map(), so each fit is bound
#' into its own 3-d array and `[[1]]` keeps only seed 1 (verified against COMPASS 1.48.0). This
#' averages across the seed dimension as intended.
mean_gamma_across_seeds <- function(fits) {
  mg <- lapply(fits, function(f) f$fit$mean_gamma)
  d <- dim(mg[[1]])
  if (any(vapply(mg, function(m) !identical(dim(m), d), logical(1)))) stop("fits differ in shape", call. = FALSE)
  out <- Reduce(`+`, mg) / length(mg)
  dimnames(out) <- dimnames(mg[[1]])
  out
}

#' Gamma draws (individuals x categories x draws) of several seeds concatenated along the draw
#' dimension, so interval summaries use every chain rather than one.
pooled_gamma_draws <- function(fits) {
  g <- lapply(fits, function(f) f$fit$gamma)
  d <- lapply(g, dim)
  if (any(vapply(d, function(x) !identical(x[1:2], d[[1]][1:2]), logical(1)))) stop("fits differ in shape", call. = FALSE)
  array(unlist(g), dim = c(d[[1]][1:2], sum(vapply(d, function(x) x[3], numeric(1)))),
        dimnames = dimnames(g[[1]])[1:2])
}

#' How many stimulated tubes share each tube's unstimulated control. COMPASS pairs every individual
#' with its own control and assumes individuals are independent; where one neg serves several
#' stimulated tubes (several peptide pools of one visit) that assumption does not hold, so the count
#' travels with the results.
neg_sharing <- function(neg_filename) {
  as.integer(table(neg_filename)[neg_filename])
}

#' Row summaries (mean, 2.5 % and 97.5 %) of a matrix of draws.
draw_summary <- function(m, prefix) {
  q <- t(apply(m, 1, stats::quantile, c(0.025, 0.975), names = FALSE))
  out <- data.frame(rowMeans(m), q[, 1], q[, 2])
  names(out) <- paste0(prefix, c("", "_lo", "_hi"))
  out
}

#' Vaccine dose number from the visit label (e.g. "V3W315" -> 3).
visit_dose <- function(x, regex = "^V(\\d+)") {
  m <- regmatches(x, regexec(regex, x))
  as.integer(vapply(m, function(v) if (length(v) >= 2) v[2] else NA_character_, character(1)))
}

#' Long frame of the subset model: one row per tube x parent x kept non-null subset, stimulated and
#' neg tubes alike. `bc` is boolean_counts.csv (filename, parent, subset, count, parent_count),
#' `tubes` the tube table with the design columns, `keep` the kept subsets. Adds stim (0/1), ps
#' (parent:subset), degree, dose and the visit factor (reference `ref_dose`), and obs (row id).
subset_model_frame <- function(bc, tubes, keep, ccfg) {
  d <- bc[bc$subset %in% keep, c("filename", "parent", "subset", "count", "parent_count")]
  cols <- intersect(c("filename", "patient", "week", "cocktail", "is_neg", "acq_day", "instrument", "operator",
                      "prior_vacc", ccfg$visit_column), names(tubes))
  d <- merge(d, tubes[, cols], by = "filename")
  d$stim <- as.integer(!as.logical(d$is_neg))
  d$ps <- paste(d$parent, d$subset, sep = ":")
  d$degree <- subset_degree(d$subset)
  d$dose <- visit_dose(as.character(d[[ccfg$visit_column]]), ccfg$dose_regex)
  lv <- sort(unique(d$dose)); lv <- c(ccfg$reference_dose, setdiff(lv, ccfg$reference_dose))
  d$visit <- factor(paste0("V", d$dose), levels = paste0("V", lv))
  d$visit_id <- paste(d$patient, d[[ccfg$visit_column]], sep = "_")   # one sample: patient x visit
  d <- d[order(d$parent, d$subset, d$filename), ]
  d$obs <- seq_len(nrow(d))
  rownames(d) <- NULL
  d
}

#' Natural-spline basis of log(1 + week) with `df` columns, named w1..wdf, plus the basis object
#' for prediction at new weeks (attribute "basis").
week_basis <- function(week, df) {
  b <- splines::ns(log1p(pmax(week, 0)), df = df)
  out <- as.data.frame(unclass(b)[, seq_len(df), drop = FALSE]); names(out) <- paste0("w", seq_len(df))
  attr(out, "basis") <- b
  out
}

week_basis_at <- function(basis, week) {
  loadNamespace("splines")   # registers predict.ns; a basis read back from an .rds does not load it
  out <- as.data.frame(stats::predict(basis, log1p(pmax(week, 0))))
  names(out) <- paste0("w", seq_len(ncol(out)))
  out
}

#' Hypothesis quantities of one parent at one time, per draw. `p_stim`, `p_neg`: draws x subsets
#' matrices of the probability of a cell falling in each kept non-null subset (stimulated and neg),
#' `degree` the degree of each column. Returns a data.frame with one row per draw:
#'   response  sum over subsets of p_stim - p_neg (net induced cells, fraction of the parent; equals
#'             the drop of the all-negative subset)
#'   induced   sum of the positive parts max(p_stim - p_neg, 0)
#'   share_poly share of the positive part carried by subsets with >= 2 markers (NA when nothing is induced)
#'   n_up      number of subsets with p_stim > p_neg
#'   n_down    number of subsets with p_stim < p_neg
hypothesis_draws <- function(p_stim, p_neg, degree) {
  dlt <- p_stim - p_neg
  pos <- pmax(dlt, 0)
  induced <- rowSums(pos)
  data.frame(response = rowSums(dlt), induced = induced,
             share_poly = ifelse(induced > 0, rowSums(pos[, degree >= 2, drop = FALSE]) / induced, NA_real_),
             n_up = rowSums(dlt > 0), n_down = rowSums(dlt < 0))
}

#' Posterior summary of a vector of draws: mean, 95 % interval, P(> 0).
post_summary <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(data.frame(mean = NA_real_, lo = NA_real_, hi = NA_real_, p_gt0 = NA_real_))
  data.frame(mean = mean(x), lo = unname(stats::quantile(x, 0.025)), hi = unname(stats::quantile(x, 0.975)), p_gt0 = mean(x > 0))
}

#' Simulate Boolean counts on a real design. `d` is a subset_model_frame() (stim, ps, visit, acq_day,
#' visit_id, cocktail, parent_count); `base` a named vector of baseline logits per ps; `effect` a
#' function(ps, visit) giving the stimulation log-odds ratio of each row. Random effects: visit x ps
#' (sd_visit), acquisition day (sd_day, shared by all subsets), pool x ps on stimulated rows (sd_pool),
#' observation level (sd_obs). Returns `d` with `count` replaced.
simulate_subset_counts <- function(d, base, effect, sd_visit = 0.4, sd_day = 0, sd_pool = 0.3, sd_obs = 0.3, seed = 1) {
  set.seed(seed)
  re <- function(g, sd) { lv <- unique(g); stats::setNames(stats::rnorm(length(lv), 0, sd), lv)[g] }
  eta <- base[d$ps] + re(paste(d$ps, d$visit_id), sd_visit) + re(d$acq_day, sd_day) +
    d$stim * (effect(d$ps, as.character(d$visit)) + re(paste(d$ps, d$cocktail), sd_pool)) + stats::rnorm(nrow(d), 0, sd_obs)
  d$count <- stats::rbinom(nrow(d), d$parent_count, stats::plogis(unname(eta)))
  d
}
