#!/usr/bin/env Rscript
# test_data/test_drift_partition.R — unit tests for src/harmonization/drift/drift_partition.R (variance
# partition per representation), the spline time term of drift_models.R and the calling rules of
# drift_stats.R. Base R; the share-recovery checks run only when lme4 is installed.
#   Rscript test_data/test_drift_partition.R
cat("=== drift partition tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_partition.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
for (m in list.files(file.path(base_dir, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)

n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") }
}
fstr <- function(f) paste(deparse(f, width.cutoff = 500), collapse = "")

cfg <- list(design = list(unit = "patient", technical = c("acq_day", "instrument", "operator"), protected_fixed = "prior_vacc",
                          time = "week", stimulation = "cocktail", time_df = 3, reference_levels = list(operator = "P")),
            controls = list(match_keys = c("patient", "week")), reference_instrument = "Canto",
            cytokine_markers = c("IFNg"), positivity = list(min_subtracted_pct = 0.05, min_fold_over_neg = 3), arcsinh_cofactor = 150)

# --- simulated cohort --------------------------------------------------------------------
set.seed(3)
np <- 12; weeks <- c(0, 6, 12, 30, 60)
tt <- expand.grid(patient = seq_len(np), week = weeks, pool = c("neg", "p1", "p2"), stringsAsFactors = FALSE)
tt$is_neg <- tt$pool == "neg"
tt$cocktail <- ifelse(tt$is_neg, "neg", paste0(tt$pool, "_", tt$patient))
tt$acq_day <- paste0("d", (tt$patient + match(tt$week, weeks)) %% 9)
tt$instrument <- ifelse(tt$acq_day %in% c("d1", "d2"), "Symph", "Canto")
tt$operator <- c("A", "P", "Y")[(as.integer(sub("d", "", tt$acq_day)) %% 3) + 1]
tt$prior_vacc <- 1 + (tt$patient %% 2)
tt$filename <- sprintf("t%03d.fcs", seq_len(nrow(tt)))
tt$patient_day <- paste(tt$patient, tt$acq_day, sep = "_")

# --- spline time term ------------------------------------------------------------------------
cat("spline time term\n")
fr <- design_frame(tt, cfg, ref = "Canto")
check("time_df non-linear columns (log-time spline), orthogonal to the line and standardised",
      all(c("week_nl1", "week_nl2", "week_nl3") %in% names(fr)) && abs(cor(fr$week_nl1, fr$week_c)) < 1e-8 &&
        abs(mean(fr$week_nl1)) < 1e-8 && abs(sd(fr$week_nl2) - 1) < 1e-8)
check("time_df 1 keeps the line only", identical(time_nl_terms(list(design = list(time = "week", time_df = 1))), character(0)))
check("endpoint formula carries the curvature", grepl("week_c \\+ week_nl1 \\+ week_nl2 \\+ week_nl3", fstr(design_formula(cfg, "endpoint", "lme4"))))
check("fixed vars carry the curvature", all(c("week_nl1", "week_nl2") %in% design_fixed_vars(cfg, "background_counts")))
check("a frame with too few distinct weeks gets zero curvature", all(design_frame(transform(tt, week = 5), cfg)$week_nl1 == 0))

# --- frames --------------------------------------------------------------------------------
cat("partition frames\n")
mfi <- data.frame(filename = tt$filename, population = "CD4", marker = "IFNg", n = 1000,
                  median_linear = 100 + 50 * !tt$is_neg + rnorm(nrow(tt), 0, 5), stringsAsFactors = FALSE)
mfi$median_asinh <- asinh(mfi$median_linear / 150)
k0 <- rbinom(nrow(tt), 20000, 0.0005); k1 <- rbinom(nrow(tt), 20000, 0.004)
fq <- data.frame(filename = tt$filename, parent = "CD4", gate = "IFNg", count = ifelse(tt$is_neg, k0, k1), parent_count = 20000,
                 is_neg = tt$is_neg, patient = tt$patient, week = tt$week, cocktail = tt$cocktail, acq_day = tt$acq_day,
                 instrument = tt$instrument, patient_day = tt$patient_day, stringsAsFactors = FALSE)
fq$pct <- 100 * fq$count / fq$parent_count
negk <- fq[fq$is_neg, ]; key <- function(d) paste(d$patient, d$week)
fq$neg_count <- negk$count[match(key(fq), key(negk))]; fq$neg_parent_count <- 20000; fq$neg_pct <- 100 * fq$neg_count / 20000
fq$neg_count[fq$is_neg] <- NA; fq$neg_pct[fq$is_neg] <- NA
fq$positive <- positivity_call(fq$pct, fq$neg_pct)$positive; fq$positive[fq$is_neg] <- NA
d_raw <- partition_frame("freq_raw", mfi, fq, tt, cfg, "Canto")
check("freq_raw keeps negs, stimulated factor, pool 'neg' for negs, operator back-filled from tubes",
      nrow(d_raw) == nrow(tt) && identical(levels(d_raw$stimulated), c("no", "yes")) && all(d_raw$pool[d_raw$is_neg] == "neg") &&
        "operator" %in% names(d_raw) && levels(d_raw$operator)[1] == "P")
d_sub <- partition_frame("freq_sub", mfi, fq, tt, cfg, "Canto")
check("freq_sub: stimulated only, y = pct - neg_pct", all(!d_sub$is_neg) && isTRUE(all.equal(d_sub$y, d_sub$pct - d_sub$neg_pct)))
d_msub <- partition_frame("mfi_sub", mfi, fq, tt, cfg, "Canto")
i <- match(d_msub$filename[1], mfi$filename); j <- which(tt$is_neg & tt$patient == tt$patient[i] & tt$week == tt$week[i])
check("mfi_sub: arcsinh stim minus the matched neg", isTRUE(all.equal(d_msub$y[1], mfi$median_asinh[i] - mfi$median_asinh[j])))
d_rsub <- partition_frame("mfi_raw_sub", mfi, fq, tt, cfg, "Canto")
check("mfi_raw_sub: linear stim minus the matched neg", isTRUE(all.equal(d_rsub$y[1], mfi$median_linear[match(d_rsub$filename[1], mfi$filename)] -
      mfi$median_linear[which(tt$is_neg & tt$patient == tt$patient[match(d_rsub$filename[1], tt$filename)] & tt$week == tt$week[match(d_rsub$filename[1], tt$filename)])])))
d_as <- partition_frame("freq_asinh", mfi, fq, tt, cfg, "Canto"); d_assub <- partition_frame("freq_asinh_sub", mfi, fq, tt, cfg, "Canto")
check("freq_asinh on all tubes with cofactor 0.05 percent; freq_asinh_sub = difference of arcsinh",
      nrow(d_as) == nrow(tt) && isTRUE(all.equal(d_as$y, asinh(d_as$pct / 0.05))) &&
        isTRUE(all.equal(d_assub$y, asinh(d_assub$pct / 0.05) - asinh(d_assub$neg_pct / 0.05))) && all(!d_assub$is_neg))
check("freq_pct is the plain percentage, all tubes", isTRUE(all.equal(partition_frame("freq_pct", mfi, fq, tt, cfg, "Canto")$y, partition_frame("freq_pct", mfi, fq, tt, cfg, "Canto")$pct)))
check("thirteen representations, each with a short label", length(PARTITION_REPS) == 13 && all(vapply(PARTITION_REPS, function(x) nzchar(x$short), logical(1))))
d_call <- partition_frame("call", mfi, fq, tt, cfg, "Canto")
check("call equals positivity_call()", identical(d_call$y, as.integer(d_call$positive %in% TRUE)) && all(d_call$y %in% 0:1))
d_ratio <- partition_frame("freq_ratio", mfi, fq, tt, cfg, "Canto")
check("freq_ratio carries the empirical-logit neg offset", isTRUE(all.equal(d_ratio$neg_offset, elogit(d_ratio$neg_count, d_ratio$neg_parent_count))))

cat("terms and formulas\n")
tm <- partition_terms(cfg, d_raw, "freq_raw")
check("unwanted: instrument, operator, day; wanted: patient, prior_vacc, stimulated, pool, week",
      setequal(names(tm$groups)[tm$groups == "unwanted"], c("instrument", "operator", "acq_day")) &&
        setequal(names(tm$groups)[tm$groups == "wanted"], c("patient", "prior_vacc", "stimulated", "pool", "week")))
check("lme4 count formula", fstr(partition_formula(tm, "freq_raw", "lme4")) ==
      "cbind(count, parent_count - count) ~ 1 + instrument + operator + prior_vacc + week_c + week_nl1 + week_nl2 + week_nl3 + stimulated + (1 | acq_day) + (1 | patient) + (1 | pool) + (1 | visit) + (1 | filename)")
check("brms ratio formula has trials and the offset", grepl("^count \\| trials\\(parent_count\\) ~ .*offset\\(neg_offset\\)$", fstr(partition_formula(partition_terms(cfg, d_ratio, "freq_ratio"), "freq_ratio", "brms"))))
check("stimulated absent when only stimulated tubes", !"stimulated" %in% partition_terms(cfg, d_sub, "freq_sub")$fixed)
check("single-level factor left out", !"instrument" %in% partition_terms(cfg, transform(d_sub, instrument = factor("Canto")), "freq_sub")$fixed)
check("fixed_term_of maps columns to factors", identical(unname(fixed_term_of(c("(Intercept)", "instrumentSymph", "operatorA", "week_c", "week_nl2", "prior_vacc2"),
      c("instrument", "operator", "prior_vacc"))), c(NA, "instrument", "operator", "week", "week", "prior_vacc")))

cat("shares\n")
vs <- variance_shares(list(instrument = c(0, 0, 1, 1)), c(patient = 2, acq_day = 0.5), 1, c(instrument = "unwanted", acq_day = "unwanted", patient = "wanted"))
check("variance_shares: fixed part as variance, shares sum to 1",
      isTRUE(all.equal(sum(vs$share), 1)) && isTRUE(all.equal(vs$variance[vs$factor == "instrument"], var(c(0, 0, 1, 1)))) && vs$group[vs$factor == "residual"] == "residual")
ss <- summarise_shares(rbind(cbind(draw = 1, vs), cbind(draw = 2, vs)), scale_var = 4)
check("summarise_shares adds totals and rescales variances",
      isTRUE(all.equal(ss$share[ss$factor == "total_unwanted"], sum(vs$share[vs$group == "unwanted"]))) &&
        isTRUE(all.equal(ss$variance[ss$factor == "patient"], 8)))
vp <- rbind(data.frame(readout = "r", representation = "freq_raw", factor = c("total_unwanted", "total_wanted"), share = c(.4, .4), lo = c(.4, .4), hi = c(.4, .4)),
            data.frame(readout = "r", representation = "freq_sub", factor = c("total_unwanted", "total_wanted"), share = c(.1, .35), lo = c(.1, .3), hi = c(.1, .45)),
            data.frame(readout = "r", representation = "freq_ratio", factor = c("total_unwanted", "total_wanted"), share = c(.5, .2), lo = c(.5, .1), hi = c(.5, .3)))
st <- subtraction_table(vp)
check("subtraction_table verdicts", st$verdict[st$alternative == "freq_sub"] == "removes unwanted, keeps wanted" &&
        st$verdict[st$alternative == "freq_ratio"] == "no reduction of unwanted" && !"call" %in% compared_reps())
check("technical_cv is the log-normal CV", isTRUE(all.equal(technical_cv(log(0.25^2 + 1)), 0.25)))

cat("calling rules\n")
st_rows <- fq[!fq$is_neg, ]; st_rows$response <- "CD4:IFNg"
cr <- calling_rules(st_rows, cfg)
check("lab column equals positivity_call, Fisher q within [0,1]", identical(cr$lab, st_rows$positive) && all(cr$fisher_q >= 0 & cr$fisher_q <= 1))
ag <- calls_agreement(cr)
check("agreement table: counts and kappa", ag$positive_lab == sum(cr$lab) && "kappa_lab_fisher" %in% names(ag))
dec <- data.frame(term = c("day", "day", "instrument"), layer = "B", verdict = c("detected, correct", "not detected", "detected, not clinically relevant"),
                  positivity_flips = c(2, 0, 0), stringsAsFactors = FALSE)
cs <- correction_summary(dec)
check("correction_summary", cs$flipping[cs$term == "day"] == 1 && cs$calls_flipped[cs$term == "day"] == 2 &&
        cs$action[cs$term == "instrument"] == "keep as covariate, no correction")

cat("visit, sample quality and counting noise\n")
qual <- data.frame(filename = tt$filename, viability = ifelse(tt$is_neg, 0.8, 0.5), live_cells = 1e4)
qual$viability[tt$is_neg & tt$patient == 1] <- 0.4
dq <- partition_frame("freq_pct", mfi, fq, tt, cfg, "Canto", qual)
check("visit-level quality comes from the unstimulated tube of the visit", all(c("qual_viability", "qual_cells", "visit") %in% names(dq)) &&
      length(unique(round(dq$qual_viability[dq$patient == 1], 6))) == 1 && mean(dq$qual_viability[dq$patient == 1]) < mean(dq$qual_viability[dq$patient != 1]))
tq <- partition_terms(cfg, dq, "freq_pct")
check("terms: sample quality is technical, visit its own group",
      all(c("qual_viability", "qual_cells") %in% tq$fixed) && "visit" %in% tq$random && tq$groups[["quality"]] == "unwanted" && tq$groups[["visit"]] == "visit")
check("fixed_term_of maps the quality columns to one source", identical(unname(fixed_term_of(c("qual_viability", "qual_cells"), tq$fixed)), c("quality", "quality")))
dc1 <- data.frame(count = c(10, 20), parent_count = c(1000, 1000), neg_count = c(5, 5), neg_parent_count = c(1000, 1000), pct = c(1, 2), neg_pct = c(.5, .5))
check("counting noise of a percentage is p(1-p)/n in percent squared; subtraction adds both tubes; logit models 0",
      isTRUE(all.equal(counting_noise_var(dc1, "freq_pct", cfg), mean(1e4 * c(.01 * .99, .02 * .98) / 1000))) &&
        counting_noise_var(dc1, "freq_sub", cfg) > counting_noise_var(dc1, "freq_pct", cfg) && counting_noise_var(dc1, "freq_raw", cfg) == 0)
s0 <- data.frame(factor = c("patient", "residual", "total_wanted"), group = c("wanted", "residual", "wanted"), share = c(.6, .4, .6), variance = c(6, 4, 6), lo = c(.6, .4, .6), hi = c(.6, .4, .6))
s1 <- split_counting_noise(s0, 1)
check("split_counting_noise moves part of the residual into counting noise, shares still sum to 1",
      isTRUE(all.equal(s1$share[s1$factor == "counting_noise"], .1)) && isTRUE(all.equal(s1$share[s1$factor == "residual"], .3)) &&
        isTRUE(all.equal(sum(s1$share[!grepl("^total_", s1$factor)]), 1)))
d_fj <- partition_frame("mfi_flowjo", transform(mfi, median_flowjo = median_asinh * 2), fq, tt, cfg, "Canto")
check("FlowJo-scale MFI setting uses median_flowjo; skipped when the column is missing",
      isTRUE(all.equal(d_fj$y, d_fj$median_flowjo)) && is.null(partition_frame("mfi_flowjo", mfi, fq, tt, cfg, "Canto")))

cat("gating levels and keep/drop rule\n")
source(file.path(base_dir, "scripts", "drift", "partition_plots.R"))
check("gating level from readout name and kind",
      identical(as.character(gating_level(c("all:IFNg", "live:IFNg", "CD4:IFNg", "CD4:IFNg"), c("mfi", "mfi", "mfi", "freq"))),
                c("no gate (all events)", "live-cell gate", "CD4 / CD8 gate", "manual cytokine gate")))
check("controls are lineage, viability and scatter readouts", identical(is_control_readout(c("live:CD4", "live:IFNg", "CD4:IFNg", "live:FSC_A")), c(TRUE, FALSE, FALSE, TRUE)))
dc <- decide_combinations(data.frame(setting = factor(c("a", "b", "c", "d")), gating = factor(c("g", "g", "g", "g")),
                                     biological_pct = c(40, 30, 45, 20), batch_pct = c(20, 25, 30, 50), bio_minus_batch = c(20, 5, 15, -30)))
check("keep the trade-off front, drop dominated and batch-heavy combinations",
      identical(dc$verdict, c("keep", "drop: another combination is better on both", "keep", "drop: more batch than biology, and dominated")) &&
        dc$better[2] == "a / g")

# --- share recovery with lme4 ------------------------------------------------------------------
if (requireNamespace("lme4", quietly = TRUE)) {
  cat("share recovery (lme4)\n")
  set.seed(11)
  big <- expand.grid(patient = 1:40, day = 1:30, rep = 1:3); big$filename <- seq_len(nrow(big))
  # effects scaled to exactly the intended variances, so the test checks the arithmetic, not sampling luck
  big$y <- (sqrt(2) * as.numeric(scale(rnorm(40))))[big$patient] + as.numeric(scale(rnorm(30)))[big$day] + rnorm(nrow(big), 0, 1)
  big$patient <- factor(big$patient); big$acq_day <- factor(big$day)
  terms <- list(fixed = character(0), random = c("acq_day", "patient"), groups = c(acq_day = "unwanted", patient = "wanted"))
  PARTITION_REPS$test <- list(kind = "mfi", tubes = "stim", family = "gaussian", label = "test")
  fit <- lme4::lmer(y ~ 1 + (1 | acq_day) + (1 | patient), data = big)
  sh <- lme4_shares(fit, terms, "test")
  check("known shares 0.5 / 0.25 / 0.25 recovered within 10 percent (absolute)",
        abs(sh$share[sh$factor == "patient"] - .5) < .1 && abs(sh$share[sh$factor == "acq_day"] - .25) < .1 && abs(sh$share[sh$factor == "residual"] - .25) < .1)
} else cat("share recovery: lme4 not installed, skipped\n")

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
