# scripts/drift/partition_plots.R — figures comparing the data representations (ggplot2)
# Sourced by 14_partition_figures.R and 10_report.Rmd. Input: a variance partition table
# (read_partition() in src/harmonization/drift/drift_partition.R). "Biological" is the wanted
# group (patient, time, protected covariates, cocktail used or not and which pool); "batch" the
# unwanted group (acquisition day and the technical factors such as instrument and operator).

# Factor colours, validated as one adjacent sequence in stack order (dataviz validate_palette.js, light
# mode: all checks pass, worst adjacent CVD dE 15.2, normal-vision dE 15.6). Batch sources warm and violet,
# biology blue/green/magenta, the visit (handling or biology, cannot be split) dark grey, counting noise mid
# grey, the rest light grey. Some colours sit below 3:1 contrast on white, so segments carry % labels.
PARTITION_COLOURS <- c(acq_day = "#eb6834", instrument = "#4a3aa7", operator = "#e34948", quality = "#eda100",
                       patient = "#2a78d6", week = "#1baf7a", prior_vacc = "#008300", cocktail = "#e87ba4",
                       visit = "#52514e", counting_noise = "#9a9a94", residual = "#d9d9d4")
PARTITION_LABELS <- c(acq_day = "acquisition day", instrument = "instrument", operator = "operator", quality = "sample quality (viability, cells)",
                      patient = "patient", week = "weeks since vaccination", prior_vacc = "prior vaccination",
                      cocktail = "cocktail (used or not, pool)", visit = "visit (handling or biology, cannot be split)",
                      counting_noise = "counting noise", residual = "unexplained rest")

# Lineage, viability and scatter markers: a negative control (no peptide response expected), shown per
# readout but left out of the summaries. Set from cfg$technical_markers by the callers.
PARTITION_CONTROL_MARKERS <- c("FSC_A", "SSC_A", "LiveDead", "CD4", "CD8")
is_control_readout <- function(readout) sub("^.*:", "", readout) %in% PARTITION_CONTROL_MARKERS

#' Gating level of a readout, ordered from no gate to the cytokine gate. MFI readouts are named
#' "<population>:<marker>" (all = all events, live = live cells, else a parent such as CD4); frequency
#' readouts "<parent>:<cytokine>" (the manual FlowJo cytokine gate).
GATING_LEVELS <- c("no gate (all events)", "live-cell gate", "CD4 / CD8 gate", "manual cytokine gate")
gating_level <- function(readout, kind) {
  pop <- sub(":.*$", "", readout)
  lv <- ifelse(kind == "mfi", ifelse(pop == "all", GATING_LEVELS[1], ifelse(pop == "live", GATING_LEVELS[2], GATING_LEVELS[3])),
               GATING_LEVELS[4])
  factor(lv, GATING_LEVELS)
}

#' Readouts that enter the summaries: every readout except the negative controls.
summary_readouts <- function(sh) unique(sh$readout[!is_control_readout(sh$readout)])

#' Long table of shares per readout, representation and factor, the stimulated and pool terms merged
#' into "cocktail", any further factor kept under its own name (coloured grey if unknown).
partition_long <- function(vp) {
  vp <- vp[vp$representation %in% compared_reps(), ]
  sh <- vp[!grepl("^total_", vp$factor), c("readout", "representation", "factor", "group", "share")]
  sh$factor[sh$factor %in% c("stimulated", "pool")] <- "cocktail"
  sh <- aggregate(share ~ readout + representation + factor + group, sh, sum)
  sh$kind <- vapply(sh$representation, function(r) PARTITION_REPS[[r]]$kind, "")
  sh$setting <- factor(vapply(sh$representation, function(r) PARTITION_REPS[[r]]$short, ""),
                       vapply(PARTITION_REPS, `[[`, "", "short"))
  sh$gating <- gating_level(sh$readout, sh$kind)
  sh$setting_name <- sh$representation
  sh$representation <- paste(sh$representation, sh$gating, sep = " | ")   # one combination = setting x gating level
  sh
}

#' The comparison table: per readout and setting the biological and batch percentages, their
#' difference and ratio. With `summary = TRUE` the medians over readouts per setting (MFI settings
#' over the readouts present in all of them, so the four MFI settings compare the same readouts).
partition_grid <- function(vp, summary = FALSE) {
  sh <- partition_long(vp)
  g <- aggregate(share ~ readout + representation + setting + gating + kind + group, sh, sum)
  w <- reshape(g, idvar = c("readout", "representation", "setting", "gating", "kind"), timevar = "group", direction = "wide")
  names(w) <- sub("^share\\.", "", names(w))
  for (k in c("wanted", "unwanted", "visit", "residual")) if (!k %in% names(w)) w[[k]] <- 0
  w$biological_pct <- 100 * w$wanted; w$batch_pct <- 100 * w$unwanted; w$visit_pct <- 100 * w$visit; w$unexplained_pct <- 100 * w$residual
  w$bio_minus_batch <- w$biological_pct - w$batch_pct
  w$bio_to_batch <- w$biological_pct / pmax(w$batch_pct, 0.1)
  w <- w[, c("readout", "kind", "gating", "representation", "setting", "biological_pct", "batch_pct", "visit_pct", "unexplained_pct", "bio_minus_batch", "bio_to_batch")]
  if (!summary) return(w[order(w$kind, w$gating, w$readout, w$setting), ])
  keep <- w$readout %in% summary_readouts(w)
  s <- aggregate(cbind(biological_pct, batch_pct, visit_pct, unexplained_pct, bio_minus_batch, bio_to_batch) ~ kind + gating + representation + setting,
                 w[keep, ], median)
  s$n_readouts <- as.integer(table(w$representation[keep])[s$representation])
  # difference and ratio of the medians (not medians of the per-readout values), so they match the two medians shown
  s$bio_minus_batch <- s$biological_pct - s$batch_pct; s$bio_to_batch <- s$biological_pct / pmax(s$batch_pct, 0.1)
  s <- decide_combinations(s)
  s[order(s$kind, s$gating, s$setting), ]
}

#' Which combinations of gating level and setting to keep. A combination is dropped when it has more
#' batch than biological variance, or when another combination has at least as much biology and at most
#' as much batch (it is dominated). The kept ones form the trade-off front. Adds `verdict` and `better`.
decide_combinations <- function(s) {
  s$verdict <- "keep"; s$better <- NA_character_
  lab <- paste(gsub("\n", " ", as.character(s$setting)), "/", s$gating)
  for (i in seq_len(nrow(s))) {
    dom <- which(s$biological_pct >= s$biological_pct[i] & s$batch_pct <= s$batch_pct[i] &
                 (s$biological_pct > s$biological_pct[i] | s$batch_pct < s$batch_pct[i]))
    if (length(dom)) { j <- dom[which.max(s$bio_minus_batch[dom])]; s$verdict[i] <- "drop: another combination is better on both"; s$better[i] <- lab[j] }
    if (s$batch_pct[i] > s$biological_pct[i]) s$verdict[i] <- if (length(dom)) "drop: more batch than biology, and dominated" else "drop: more batch than biology"
  }
  s
}

partition_theme <- function() ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "#f0f0ee", colour = NA), legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(size = 8, lineheight = .9))

#' Median shares per combination (gating level x setting) over the summary readouts, each factor's median
#' scaled so that a group's segments add up to the median of that group's total (medians of parts do not add
#' up to the median of the sum); the numbers then match partition_grid_summary.csv.
summary_shares <- function(sh, groups = c("wanted", "unwanted", "visit")) {
  sh <- sh[sh$readout %in% summary_readouts(sh) & sh$group %in% groups, ]
  tot_med <- aggregate(share ~ representation + group, aggregate(share ~ readout + representation + group, sh, sum), median)
  sh <- aggregate(share ~ representation + setting + gating + kind + factor + group, sh, median)
  part <- aggregate(share ~ representation + group, sh, sum)
  k <- paste(sh$representation, sh$group)
  sh$share <- sh$share * tot_med$share[match(k, paste(tot_med$representation, tot_med$group))] /
    part$share[match(k, paste(part$representation, part$group))]
  sh
}

#' Paired stacked bars: in every combination a biological bar and a batch bar, each stacked by its factors,
#' in percent of the variance. `kind` "mfi" or "freq". Summary: one row per gating level (medians over
#' readouts); `by_readout`: one row per readout.
plot_partition_paired <- function(vp, kind, by_readout = FALSE) {
  sh <- partition_long(vp); sh <- sh[sh$group %in% c("wanted", "unwanted") & sh$kind == kind, ]   # the visit share is in the full bars
  if (!by_readout) { sh <- summary_shares(sh); sh$row <- sh$gating } else sh$row <- sh$readout
  sh$bar <- factor(c(wanted = "bio", unwanted = "batch", visit = "visit")[sh$group], c("bio", "batch", "visit"))
  sh$setting <- droplevels(sh$setting)
  sh$factor <- factor(sh$factor, rev(intersect(names(PARTITION_COLOURS), unique(sh$factor))))
  sh$pct <- 100 * sh$share
  tot <- aggregate(pct ~ row + setting + bar, sh, sum)
  ggplot2::ggplot(sh, ggplot2::aes(bar, pct, fill = factor)) +
    ggplot2::geom_col(width = 0.85, colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(pct >= 5, sprintf("%.0f", pct), "")), position = ggplot2::position_stack(vjust = .5),
                       size = 2.6, colour = "white") +
    ggplot2::geom_text(data = tot, ggplot2::aes(bar, pct, label = sprintf("%.0f%%", pct)), inherit.aes = FALSE, vjust = -0.4, size = 2.7, colour = "#0b0b0b") +
    ggplot2::facet_grid(row ~ setting, switch = "x", labeller = ggplot2::labeller(row = ggplot2::label_wrap_gen(16))) +
    ggplot2::scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), labels = function(v) paste0(v, "%")) +
    ggplot2::scale_fill_manual(values = PARTITION_COLOURS, labels = PARTITION_LABELS, breaks = names(PARTITION_COLOURS), drop = TRUE) +
    partition_theme() + ggplot2::theme(strip.placement = "outside", panel.spacing.x = grid::unit(2, "mm"),
                                       strip.text = ggplot2::element_text(size = 8)) +
    ggplot2::labs(x = NULL, y = "share of variance", fill = NULL,
      title = paste(if (kind == "mfi") "MFI:" else "Frequencies:", "biological against batch variance per gating level and setting"),
      subtitle = paste(if (by_readout) "per readout;" else "rows = gating level; bar height = median over readouts of the group total;",
                       "higher bio and lower batch is better; visit = handling or biology, cannot be split; the rest of each 100% is unexplained"))
}

#' Full 100% stacked bar per combination: every factor in its own colour plus the unexplained rest.
#' `kind` "mfi" or "freq"; rows = gating level, medians over readouts renormalised to 100%.
plot_partition_full <- function(vp, kind) {
  sh <- partition_long(vp); sh <- sh[sh$kind == kind, ]
  sh <- summary_shares(sh, c("wanted", "unwanted", "visit", "residual"))
  s <- aggregate(share ~ representation, sh, sum); sh$share <- sh$share / s$share[match(sh$representation, s$representation)]
  sh$factor <- factor(sh$factor, rev(intersect(names(PARTITION_COLOURS), unique(sh$factor))))
  sh$setting <- droplevels(sh$setting)
  sh$pct <- 100 * sh$share
  ggplot2::ggplot(sh, ggplot2::aes(setting, pct, fill = factor)) +
    ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(pct >= 5, sprintf("%.0f", pct), ""), colour = ifelse(factor == "residual", "#0b0b0b", "white")),
                       position = ggplot2::position_stack(vjust = .5), size = 2.6) + ggplot2::scale_colour_identity() +
    ggplot2::facet_wrap(~gating, nrow = 1) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 25), labels = function(v) paste0(v, "%")) +
    ggplot2::scale_fill_manual(values = PARTITION_COLOURS, labels = PARTITION_LABELS, breaks = names(PARTITION_COLOURS), drop = TRUE) +
    partition_theme() + ggplot2::labs(x = NULL, y = "share of variance", fill = NULL,
      title = paste(if (kind == "mfi") "MFI:" else "Frequencies:", "every source of variance per gating level and setting"),
      subtitle = "median over readouts, renormalised to 100%")
}

#' Which combinations of gating level and setting are worth keeping: biology minus batch per cell, with
#' the verdict of decide_combinations() (crossed cells are dropped). One panel per kind.
plot_combinations <- function(vp) {
  s <- partition_grid(vp, summary = TRUE)
  s$setting_lab <- droplevels(s$setting)   # two-line labels, so neighbouring settings do not collide
  s$panel <- factor(ifelse(s$kind == "mfi", "MFI", "gated frequency"), c("MFI", "gated frequency"))
  s$lab <- sprintf("bio %.0f%%\nbatch %.0f%%", s$biological_pct, s$batch_pct)
  s$keep <- s$verdict == "keep"
  lim <- max(abs(s$bio_minus_batch))
  ggplot2::ggplot(s, ggplot2::aes(setting_lab, gating)) +
    ggplot2::geom_tile(ggplot2::aes(fill = bio_minus_batch), colour = "white", linewidth = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = lab, fontface = ifelse(keep, "bold", "plain")), size = 3, lineheight = .9,
                       colour = ifelse(s$keep, "#0b0b0b", "#52514e")) +
    ggplot2::geom_point(data = s[!s$keep, ], shape = 4, size = 14, colour = "#52514e", alpha = .35) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::scale_fill_gradient2(low = "#eb6834", mid = "#f0f0ee", high = "#2a78d6", midpoint = 0, limits = c(-lim, lim),
                                  name = "biology minus batch\n(percentage points)") +
    partition_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9), panel.grid = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::labs(x = NULL, y = NULL, title = "Which combinations of gating level and data setting are worth keeping?",
      subtitle = "blue = more biology than batch; crossed = dropped (more batch than biology, or another combination is better on both); bold = kept")
}

#' brms against lme4: the same model fitted by the Bayesian and the classical engine. One point per
#' readout x setting x factor (stimulated and pool merged into cocktail); on the diagonal they agree.
#' Returns list(plot, table) where the table holds, per setting, the correlation and the median absolute
#' difference of the shares in percentage points.
plot_engine_agreement <- function(out) {
  rd <- function(e) { p <- file.path(out, sprintf("variance_partition_%s.csv", e)); if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) }
  b <- rd("brms"); l <- rd("lme4")
  if (is.null(b) || is.null(l)) return(NULL)
  prep <- function(v) { v <- v[v$variant == "full" & !grepl("^total_", v$factor) & v$representation %in% compared_reps(), ]; v$factor[v$factor %in% c("stimulated", "pool")] <- "cocktail"
    aggregate(share ~ readout + representation + factor, v, sum) }
  m <- merge(prep(b), prep(l), by = c("readout", "representation", "factor"), suffixes = c("_brms", "_lme4"))
  m$setting <- factor(vapply(m$representation, function(r) gsub("\n", " ", PARTITION_REPS[[r]]$short), ""),
                      gsub("\n", " ", vapply(PARTITION_REPS, `[[`, "", "short")))
  m$factor <- factor(m$factor, intersect(names(PARTITION_COLOURS), unique(m$factor)))
  tab <- do.call(rbind, lapply(split(m, m$setting, drop = TRUE), function(x) data.frame(setting = x$setting[1], n = nrow(x),
    correlation = cor(x$share_brms, x$share_lme4), median_abs_diff_pp = 100 * median(abs(x$share_brms - x$share_lme4)),
    max_abs_diff_pp = 100 * max(abs(x$share_brms - x$share_lme4)))))
  p <- ggplot2::ggplot(m, ggplot2::aes(100 * share_lme4, 100 * share_brms, colour = factor)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#9a9a94") +
    ggplot2::geom_point(size = 2, alpha = .85) + ggplot2::facet_wrap(~setting, ncol = 4) + ggplot2::coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
    ggplot2::scale_colour_manual(values = PARTITION_COLOURS, labels = PARTITION_LABELS, breaks = names(PARTITION_COLOURS), drop = TRUE) +
    partition_theme() + ggplot2::labs(x = "lme4 (classical), % of variance", y = "brms (Bayesian), % of variance", colour = NULL,
      title = "The same models fitted by two engines: points on the diagonal agree",
      subtitle = "one point per readout and source of variance")
  list(plot = p, table = tab[order(tab$setting), ])
}

# ---- explanatory figures ---------------------------------------------------------------------------

#' What was run, as a diagram: inputs -> readouts -> settings -> one model per readout and setting ->
#' shares -> comparison. No data needed.
plot_explain_workflow <- function() {
  box <- data.frame(
    x = c(1, 2.6, 4.2, 5.8, 7.4),
    title = c("1  Tubes", "2  Readouts", "3  Settings", "4  One model each", "5  Compare"),
    body = c("328 FCS files\n20 patients, 23 days\n2 instruments, 3 operators\nFlowJo gates\nsample sheet",
             "Frequency: % of CD4 or CD8\ncells in a cytokine gate\n(8 gates)\n\nMFI: brightness of a marker\n(13 marker x population)",
             "MFI: linear, arcsinh,\nFlowJo display scale\nFrequency: %, arcsinh, logit\n\neach plain and\nstimulated minus\nunstimulated tube",
             "Splits the differences\nbetween tubes into:\nbatch: day, instrument,\noperator, sample quality\nbiology: patient, weeks,\nprior vaccination, cocktail\nvisit (cannot be split)\ncounting noise + rest",
             "Per setting:\nbiological %  vs  batch %\n\nbest setting: most\nbiology, least batch\n\nbrms checked with lme4"),
    stringsAsFactors = FALSE)
  ggplot2::ggplot(box) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x - .7, xmax = x + .7, ymin = 0, ymax = 1), fill = "#f0f0ee", colour = "#c9c9c4") +
    ggplot2::geom_text(ggplot2::aes(x, .9, label = title), fontface = "bold", size = 4.6, colour = "#0b0b0b") +
    ggplot2::geom_text(ggplot2::aes(x, .43, label = body), size = 3.7, colour = "#52514e", lineheight = .95) +
    ggplot2::geom_segment(data = data.frame(x = box$x[-5] + .72, xend = box$x[-1] - .72), ggplot2::aes(x = x, xend = xend, y = .5, yend = .5),
                          arrow = grid::arrow(length = grid::unit(2.5, "mm")), colour = "#52514e") +
    ggplot2::coord_cartesian(xlim = c(.25, 8.15), ylim = c(-.02, 1.02), expand = FALSE) + ggplot2::theme_void() +
    ggplot2::labs(title = "What was run: from tubes to a comparison of data settings")
}

#' One readout tube by tube against acquisition date: the unstimulated tubes, the stimulated tubes and
#' stimulated minus unstimulated. Colour = operator, shape = instrument, black bar = median per day, on a
#' compressed axis so low and high responders are both visible. Shows what "batch" means (days whose tubes
#' sit together above or below the others) and what subtraction removes (the day's background).
plot_explain_example <- function(fq, tt, response = "CD4:IFNg") {
  op_col <- c("#2a78d6", "#eb6834", "#1baf7a")   # validated all-pairs (dataviz validator, light mode)
  d <- fq[paste(fq$parent, fq$gate, sep = ":") == response, c("filename", "is_neg", "pct", "neg_pct")]
  d <- merge(d, tt[, intersect(c("filename", "acq_date", "instrument", "operator"), names(tt))], by = "filename")
  if (!"operator" %in% names(d)) d$operator <- "not recorded"
  d$acq_date <- as.Date(d$acq_date)
  lv <- c("1  unstimulated tubes (background)", "2  stimulated tubes", "3  stimulated minus its unstimulated tube")
  long <- rbind(data.frame(d[d$is_neg, ], panel = lv[1], value = d$pct[d$is_neg]),
                data.frame(d[!d$is_neg, ], panel = lv[2], value = d$pct[!d$is_neg]),
                data.frame(d[!d$is_neg & !is.na(d$neg_pct), ], panel = lv[3], value = (d$pct - d$neg_pct)[!d$is_neg & !is.na(d$neg_pct)]))
  long$panel <- factor(long$panel, lv)
  dm <- aggregate(value ~ panel + acq_date, long, median)
  ggplot2::ggplot(long, ggplot2::aes(acq_date, value)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#c9c9c4") +
    ggplot2::geom_point(ggplot2::aes(colour = operator, shape = instrument), size = 1.9, alpha = .75,
                        position = ggplot2::position_jitter(width = 4, height = 0, seed = 1)) +
    ggplot2::geom_point(data = dm, shape = 95, size = 9, colour = "#0b0b0b") +
    ggplot2::facet_wrap(~panel, ncol = 1) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.05, base = 10), breaks = c(-1, -0.1, 0, 0.1, 1, 10)) +
    ggplot2::scale_colour_manual(values = setNames(op_col[seq_along(sort(unique(d$operator)))], sort(unique(d$operator)))) +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 3)[seq_along(unique(d$instrument))]) +
    partition_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9)) +
    ggplot2::labs(x = "acquisition date", y = paste(response, "(% of parent, compressed axis)"), colour = "operator", shape = "instrument",
      title = paste("What batch looks like:", response, "tube by tube"),
      subtitle = paste0("black bar = median per acquisition day.\n1: the background differs from day to day; 2: the stimulated tubes carry that background plus the response;\n",
                        "3: subtracting each tube's own background removes the day's share, at the price of adding the background tube's noise"))
}

#' Scorecard: per combination (setting / gating level) one line from the batch % (orange) to the biological %
#' (blue), sorted by the difference; dropped combinations are greyed. The number is biology minus batch.
plot_explain_scorecard <- function(vp) {
  s <- partition_grid(vp, summary = TRUE)
  s$lab <- paste0(gsub("\n", " ", as.character(s$setting)), "  /  ", s$gating)
  s$lab <- factor(s$lab, s$lab[order(s$bio_minus_batch)])
  s$panel <- factor(ifelse(s$kind == "mfi", "MFI", "gated frequency"), c("gated frequency", "MFI"))
  s$keep <- ifelse(s$verdict == "keep", "kept", "dropped")
  long <- rbind(data.frame(s[, c("panel", "lab", "keep")], what = "batch (day, instrument, operator)", pct = s$batch_pct),
                data.frame(s[, c("panel", "lab", "keep")], what = "biology (patient, weeks, prior vaccination, cocktail)", pct = s$biological_pct))
  ggplot2::ggplot(s, ggplot2::aes(y = lab)) +
    ggplot2::geom_segment(ggplot2::aes(x = batch_pct, xend = biological_pct, yend = lab), colour = "#c9c9c4", linewidth = 1.2) +
    ggplot2::geom_point(data = long, ggplot2::aes(x = pct, colour = what, alpha = keep), size = 3.2) +
    ggplot2::geom_text(ggplot2::aes(x = pmax(batch_pct, biological_pct) + 4, label = sprintf("%+.0f", round(bio_minus_batch) + 0), alpha = keep),
                       hjust = 0, size = 3, colour = "#0b0b0b") +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_manual(values = c("batch (day, instrument, operator)" = "#eb6834", "biology (patient, weeks, prior vaccination, cocktail)" = "#2a78d6")) +
    ggplot2::scale_alpha_manual(values = c(kept = 1, dropped = .3), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, 105), labels = function(v) paste0(v, "%")) +
    partition_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "#f0f0ee"), axis.text.x = ggplot2::element_text(size = 9),
                                       axis.text.y = ggplot2::element_text(size = 8)) +
    ggplot2::labs(x = "share of variance (median over readouts)", y = NULL, colour = NULL,
      title = "Scorecard: every combination of data setting and gating level, best at the top",
      subtitle = "number = biology minus batch, percentage points; faded = dropped (more batch than biology, or another combination is better on both)")
}

#' brms against the hand-written Stan model for the correction question: every acquisition-day,
#' instrument and operator effect of the decision step (fitter_agreement.csv), the Stan posterior mean
#' against the brms one, both in units of the brms posterior SD. Agreement = Stan mean inside the brms
#' 95 percent interval. Returns list(plot, table) or NULL.
plot_engine_stan <- function(out) {
  p <- file.path(out, "fitter_agreement.csv"); if (!file.exists(p)) return(NULL)
  fa <- read.csv(p, stringsAsFactors = FALSE); fa <- fa[fa$model_b == "stan", ]
  if (!nrow(fa)) return(NULL)
  sd_a <- (fa$hi - fa$lo) / 3.92; fa$za <- fa$mean / sd_a; fa$zb <- fa$mean_b / sd_a
  fa$term <- factor(ifelse(fa$term == "day", "acquisition day", fa$term), c("acquisition day", "instrument", "operator"))
  fa$layer <- factor(c(A = "unstimulated tubes:\nlineage, viability, scatter", BG = "unstimulated tubes:\ncytokine background", B = "stimulated tubes:\ncytokine response")[fa$layer])
  fa$agree_f <- ifelse(fa$agree %in% TRUE, "Stan inside the brms 95% interval", "outside")
  tab <- aggregate(agree ~ term, fa, function(x) sprintf("%d of %d", sum(x), length(x)))
  pl <- ggplot2::ggplot(fa, ggplot2::aes(za, zb, colour = agree_f)) + ggplot2::geom_abline(colour = "#9a9a94", linetype = 2) +
    ggplot2::geom_point(alpha = .6, size = 1.6) + ggplot2::facet_grid(layer ~ term) +
    ggplot2::scale_colour_manual(values = c("Stan inside the brms 95% interval" = "#2a78d6", "outside" = "#e34948"), drop = TRUE) +
    partition_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9)) +
    ggplot2::labs(x = "brms effect (in units of its posterior SD)", y = "Stan effect (same units)", colour = NULL,
      title = sprintf("Correction question: brms and the hand-written Stan model agree on %d of %d effects", sum(fa$agree), nrow(fa)),
      subtitle = "one point per acquisition day, instrument contrast or operator contrast and per response; on the diagonal = same estimate")
  list(plot = pl, table = tab)
}

#' BayesFlow against brms from an earlier run kept in `out`/bayesflow_<date>/ (its fitter_agreement.csv and
#' bf_v*_acceptance.csv): every day and instrument effect both estimated, in units of the brms posterior SD.
#' The earlier model had no operator term and a straight-line week term, so the panel is dated and not
#' comparable with the current fits. Returns list(plot, acceptance, date) or NULL.
plot_bayesflow_history <- function(out) {
  dirs <- sort(list.files(out, "^bayesflow_[0-9-]+$", full.names = TRUE), decreasing = TRUE)
  if (!length(dirs)) return(NULL)
  d <- dirs[1]; date <- sub("^bayesflow_", "", basename(d))
  fa <- read.csv(file.path(d, "fitter_agreement.csv"), stringsAsFactors = FALSE); fa <- fa[grepl("^bayesflow", fa$model_b), ]
  if (!nrow(fa)) return(NULL)
  sd_a <- (fa$hi - fa$lo) / 3.92; fa$za <- fa$mean / sd_a; fa$zb <- fa$mean_b / sd_a
  fa$variant <- factor(c(bayesflow_v1 = "v1: tube summaries", bayesflow_v2 = "v2: single cells")[fa$model_b])
  fa$term <- ifelse(fa$term == "day", "acquisition day", fa$term)
  fa$agree_f <- ifelse(fa$agree %in% TRUE, "inside the brms 95% interval", "outside")
  acc <- do.call(rbind, lapply(list.files(d, "^bf_v[0-9]_acceptance\\.csv$", full.names = TRUE), read.csv))
  lab <- aggregate(agree ~ variant + term, fa, function(x) sprintf("%d of %d agree", sum(x), length(x)))
  pl <- ggplot2::ggplot(fa, ggplot2::aes(za, zb, colour = agree_f)) + ggplot2::geom_abline(colour = "#9a9a94", linetype = 2) +
    ggplot2::geom_hline(yintercept = 0, colour = "#e8e8e5") + ggplot2::geom_point(alpha = .6, size = 1.6) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = -Inf, y = Inf, label = agree), inherit.aes = FALSE, hjust = -0.1, vjust = 1.5, size = 3.2) +
    ggplot2::facet_grid(variant ~ term, scales = "free") +
    ggplot2::scale_colour_manual(values = c("inside the brms 95% interval" = "#2a78d6", "outside" = "#e34948")) +
    partition_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9)) +
    ggplot2::labs(x = "brms effect (in units of its posterior SD)", y = "BayesFlow effect (same units)", colour = NULL,
      title = paste("BayesFlow against brms, run of", date, "(older model: no operator term, straight-line week)"),
      subtitle = "points on the diagonal agree; points along the horizontal line are effects BayesFlow shrank to zero")
  list(plot = pl, acceptance = acc, date = date)
}

# ---- acquisition audit figures ----------------------------------------------------------------------

#' How much each factor determines the others (factor_overlap.csv of step 11b): 1 = the row factor fully
#' determines the column factor, so the two cannot be separate model terms.
plot_factor_overlap <- function(M) {
  d <- expand.grid(a = rownames(M), b = colnames(M), stringsAsFactors = FALSE); d$v <- as.vector(M)
  d$a <- factor(d$a, rev(rownames(M))); d$b <- factor(d$b, colnames(M))
  ggplot2::ggplot(d, ggplot2::aes(b, a, fill = v)) + ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(v), "", sprintf("%.2f", v)), colour = ifelse(!is.na(v) & v > .6, "white", "#0b0b0b")), size = 3) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient(low = "#f0f0ee", high = "#2a78d6", limits = c(0, 1), na.value = "white", name = "determines\n(0 to 1)") +
    partition_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "... determines this factor", y = "Knowing this factor ...", title = "Which sources overlap? (stimulated tubes)",
      subtitle = "near 1: the two cannot be separate model terms; the row factor already contains the column factor")
}

#' Detector voltage per marker against acquisition date, coloured by operator, shaped by instrument.
plot_voltages <- function(acq, tt) {
  vc <- grep("^V_", names(acq), value = TRUE); vc <- vc[!grepl("FSC|SSC", vc)]
  d <- merge(acq[, c("filename", vc)], tt[, intersect(c("filename", "acq_date", "instrument", "operator"), names(tt))], by = "filename")
  long <- do.call(rbind, lapply(vc, function(v) data.frame(d[, intersect(c("acq_date", "instrument", "operator"), names(d))], marker = sub("^V_", "", v), voltage = d[[v]])))
  long <- unique(long[is.finite(long$voltage), ]); long$acq_date <- as.Date(long$acq_date)
  if (!"operator" %in% names(long)) long$operator <- "not recorded"
  op_col <- c("#2a78d6", "#eb6834", "#1baf7a")
  ggplot2::ggplot(long, ggplot2::aes(acq_date, voltage, colour = operator, shape = instrument)) + ggplot2::geom_point(size = 2.2, alpha = .85) +
    ggplot2::facet_wrap(~marker, scales = "free_y", ncol = 4) +
    ggplot2::scale_colour_manual(values = setNames(op_col[seq_along(sort(unique(long$operator)))], sort(unique(long$operator)))) +
    partition_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8)) +
    ggplot2::labs(x = "acquisition date", y = "detector voltage (V)", colour = "operator", shape = "instrument",
      title = "Detector voltages per marker over time", subtitle = "a step in voltage changes the brightness of every cell measured on that detector")
}

#' Raw-unit position of every cytokine gate per sample against acquisition date (gate_positions.csv), log scale.
plot_gate_positions <- function(gp, tt, cytokine_regex = "IFN|INF|IL|TNF|CD107") {
  g <- gp[grepl(cytokine_regex, sub(".*/", "", gp$gate), ignore.case = TRUE) & !grepl("FSC|SSC", gp$dimension), ]
  g$position <- ifelse(is.na(g$lower_raw), g$upper_raw, g$lower_raw)
  g$parent <- sub("^.*/(CD4|CD8)/.*$", "\\1", g$gate); g <- g[g$parent %in% c("CD4", "CD8") & is.finite(g$position) & g$position > 0, ]
  leaf <- sub(".*/", "", g$gate)   # one name per cytokine whatever the FlowJo spelling (IFNy / INFy, TNF / TNFa ...)
  g$gate_name <- paste(g$parent, ifelse(grepl("IFN|INF", leaf, ignore.case = TRUE), "IFNg", ifelse(grepl("IL", leaf), "IL2",
                       ifelse(grepl("TNF", leaf, ignore.case = TRUE), "TNFa", ifelse(grepl("CD107", leaf), "CD107a", leaf)))))
  g <- merge(g, tt[, intersect(c("filename", "acq_date", "instrument", "operator"), names(tt))], by = "filename"); g$acq_date <- as.Date(g$acq_date)
  if (!"operator" %in% names(g)) g$operator <- "not recorded"
  op_col <- c("#2a78d6", "#eb6834", "#1baf7a")
  ggplot2::ggplot(unique(g[, c("gate_name", "acq_date", "position", "operator", "instrument")]), ggplot2::aes(acq_date, position, colour = operator, shape = instrument)) +
    ggplot2::geom_point(size = 1.9, alpha = .8) + ggplot2::scale_y_log10() + ggplot2::facet_wrap(~gate_name, ncol = 4) +
    ggplot2::scale_colour_manual(values = setNames(op_col[seq_along(sort(unique(g$operator)))], sort(unique(g$operator)))) +
    partition_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8)) +
    ggplot2::labs(x = "acquisition date", y = "gate threshold (raw fluorescence, log scale)", colour = "operator", shape = "instrument",
      title = "Where the cytokine gates were placed, per sample", subtitle = "one point per sample; gates follow the detector voltages of the day")
}
