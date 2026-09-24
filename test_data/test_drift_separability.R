#!/usr/bin/env Rscript
# test_data/test_drift_separability.R — tests for drift_separability.R
# Base R; the lme4-dependent design_power() part runs only when lme4 is installed.
#   Rscript test_data/test_drift_separability.R
cat("=== drift separability tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_separability.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
for (m in list.files(file.path(base_dir, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)

n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") }
}
cfg <- list(reference_instrument = "Canto", seed = 3,
            design = list(unit = "patient", technical = c("acq_day", "instrument"), protected_fixed = "prior_vacc",
                          time = "week", stimulation = "cocktail"))

# a design: `n_pat` patients, each with weeks 0/26/52; by default one day per week
make_design <- function(n_pat = 4, weeks = c(0, 26, 52), cocktails = c("neg", "A", "B"), instr = "Canto") {
  rows <- expand.grid(patient = seq_len(n_pat), week = weeks, cocktail = cocktails, stringsAsFactors = FALSE)
  rows$acq_day <- sprintf("2025-01-%02d", rows$patient * 10 + match(rows$week, weeks))
  rows$instrument <- instr; rows$prior_vacc <- 1 + rows$patient %% 2
  rows$is_neg <- rows$cocktail == "neg"
  rows$filename <- sprintf("P%d_W%d_%s.fcs", rows$patient, rows$week, rows$cocktail)
  rows$patient_day <- paste(rows$patient, rows$acq_day, sep = "_")
  rows
}

cat("confounding_index\n")
tt <- make_design()
ci <- confounding_index(tt)
check("week a function of day within patient -> index 1", abs(ci$index - 1) < 1e-12)
check("no multi-week days, every patient on 3 days", ci$n_multi_week_days == 0 && ci$n_subjects_multi_day == 4 && ci$n_subjects_multi_instrument == 0)
tt2 <- tt; tt2$acq_day[tt2$patient == 1 & tt2$week == 26] <- tt2$acq_day[tt2$patient == 1 & tt2$week == 0][1]
ci2 <- confounding_index(tt2)
check("two weeks on one day -> index below 1", ci2$index < 1 && ci2$index > 0.8)
check("one multi-week day counted", ci2$n_multi_week_days == 1)
tt3 <- tt; tt3$acq_day <- sprintf("2025-02-%02d", tt3$patient); tt3$instrument[tt3$patient == 2] <- rep_len(c("Canto", "Symphony"), sum(tt3$patient == 2))
ci3 <- confounding_index(tt3)
check("one day per patient -> index 0, no multi-day subject, one multi-instrument", abs(ci3$index) < 1e-12 && ci3$n_subjects_multi_day == 0 && ci3$n_subjects_multi_instrument == 1)
check("single week per patient -> NA index", is.na(confounding_index(make_design(weeks = 10))$index))

cat("simulate_endpoint\n")
d <- simulate_endpoint(tt, cfg, week_slope = 1, day_sd = 0.5, instr_shift = 0.4, seed = 7)
check("one row per stimulated tube", nrow(d) == sum(!tt$is_neg) && !any(d$cocktail == "neg"))
check("count within parent_count in [20000, 60000]", all(d$count <= d$parent_count) && all(d$parent_count >= 20000 & d$parent_count <= 60000))
check("design factors and centred week", is.factor(d$patient) && is.factor(d$acq_day) && is.factor(d$instrument) && abs(mean(d$week_c)) < 1e-9 && abs(diff(range(d$week_c)) - 1) < 1e-9)
d2 <- simulate_endpoint(tt, cfg, week_slope = 1, day_sd = 0.5, instr_shift = 0.4, seed = 7)
d3 <- simulate_endpoint(tt, cfg, week_slope = 1, day_sd = 0.5, instr_shift = 0.4, seed = 8)
check("same seed reproduces, other seed differs", identical(d$count, d2$count) && !identical(d$count, d3$count))
big <- simulate_endpoint(make_design(n_pat = 30), cfg, week_slope = 3, day_sd = 0, instr_shift = 0, seed = 1)
r <- tapply(big$count / big$parent_count, big$week, mean)
check("positive week slope raises the frequency over the weeks", all(diff(r) > 0))
ti <- make_design(instr = "Canto"); ti$instrument[ti$patient > 2] <- "Symphony"
di <- simulate_endpoint(ti, cfg, week_slope = 0, day_sd = 0, instr_shift = 3, seed = 2)
check("instrument shift applies to the non-reference instrument", mean(di$count[di$instrument == "Symphony"] / di$parent_count[di$instrument == "Symphony"]) >
        mean(di$count[di$instrument == "Canto"] / di$parent_count[di$instrument == "Canto"]))

cat("separability_verdict\n")
pw <- function(wa, da, db, wc) rbind(
  data.frame(scenario = "a_both", rep = 1:2, week_slope_true = 1, week_slope_hat = wa, day_sd_true = 0.5, day_sd_hat = da, instr_true = 0.4, instr_hat = 0.4, converged = TRUE),
  data.frame(scenario = "b_week_only", rep = 1:2, week_slope_true = 1, week_slope_hat = 1, day_sd_true = 0, day_sd_hat = db, instr_true = 0.4, instr_hat = 0.4, converged = TRUE),
  data.frame(scenario = "c_day_only", rep = 1:2, week_slope_true = 0, week_slope_hat = wc, day_sd_true = 0.5, day_sd_hat = 0.5, instr_true = 0.4, instr_hat = 0.4, converged = TRUE))
good <- pw(wa = c(0.95, 1.05), da = c(0.45, 0.55), db = c(0.02, 0.05), wc = c(0.01, -0.01))
v <- separability_verdict(ci2, good)
check("index < 0.95, clean recovery -> separable", v$label == "separable" && grepl("week-slope leakage 0 percent", v$explanation))
v1 <- separability_verdict(ci, good)
check("index 1 but clean recovery -> separable, with the multi-week caveat", v1$label == "separable" && grepl("multi-week acquisition days", v1$explanation))
v3 <- separability_verdict(ci3, good)
check("no multi-day subject -> not separable", v3$label == "not separable: add bridging acquisitions" && grepl("bridging", v3$explanation))
leaky <- pw(wa = c(0.5, 0.6), da = c(0.5, 0.5), db = c(0.02, 0.05), wc = 0)
check("week leakage 45 percent -> unstimulated layer", separability_verdict(ci2, leaky)$label == "separable only through the unstimulated layer")
false_day <- pw(wa = 1, da = 0.5, db = c(0.3, 0.02), wc = 0)
check("false day detection 50 percent -> unstimulated layer", separability_verdict(ci2, false_day)$label == "separable only through the unstimulated layer")
ps <- power_summary(false_day)
check("power_summary fields", abs(ps$false_day_rate - 0.5) < 1e-12 && ps$leakage_week == 0 && ps$n_fit == 6)
check("without a power table the index decides", separability_verdict(ci2, NULL)$label == "separable" && separability_verdict(ci, NULL)$label == "separable only through the unstimulated layer")

if (requireNamespace("lme4", quietly = TRUE)) {
  cat("design_power (lme4)\n")
  tb <- make_design(n_pat = 8, weeks = c(0, 26, 52, 78), cocktails = c("neg", "A", "B", "C"))
  # bridge: patients 1 to 4 have their week-26 tubes acquired on the week-0 day of the next patient
  for (p in 1:4) tb$acq_day[tb$patient == p & tb$week == 26] <- tb$acq_day[tb$patient == p + 1 & tb$week == 0][1]
  tb$instrument[tb$patient > 5] <- "Symphony"
  p <- design_power(tb, cfg, n_rep = 2, seed = 5)
  check("tidy table: 3 scenarios x 2 reps, expected columns", nrow(p) == 6 && identical(sort(unique(p$scenario)), c("a_both", "b_week_only", "c_day_only")) &&
          all(c("week_slope_true", "week_slope_hat", "day_sd_true", "day_sd_hat", "instr_true", "instr_hat", "converged") %in% names(p)))
  check("truth columns follow the scenarios", all(p$day_sd_true[p$scenario == "b_week_only"] == 0) && all(p$week_slope_true[p$scenario == "c_day_only"] == 0))
  check("at least one fit returned estimates", sum(!is.na(p$week_slope_hat)) >= 1)
  a <- p[p$scenario == "a_both" & !is.na(p$week_slope_hat), ]
  check("week slope of scenario a within 0.5 logit/year of the truth", nrow(a) > 0 && all(abs(a$week_slope_hat - a$week_slope_true) < 0.5))
} else cat("design_power: lme4 not installed, skipped\n")

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
