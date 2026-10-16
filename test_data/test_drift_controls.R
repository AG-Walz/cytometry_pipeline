#!/usr/bin/env Rscript
# test_data/test_drift_controls.R: unit tests for control matching and
# background flips in src/harmonization/drift/drift_stats.R. Base R only.
#   Rscript test_data/test_drift_controls.R
cat("=== MTB drift control-matching tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_controls.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
for (m in list.files(file.path(base_dir, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)

n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") }
}

ref <- "Canto"; oth <- "Symphony"
controls <- list(match_keys = c("patient", "week"),
                 fallbacks = c("same_patient_day", "same_patient_reference_instrument", "pooled_reference_instrument"))
tube <- function(fn, patient, week, cocktail, day, instr, count, parent_count, parent = "CD4", gate = "IFNg")
  data.frame(filename = fn, parent = parent, gate = gate, count = count, parent_count = parent_count,
             patient = patient, week = week, cocktail = cocktail, is_neg = cocktail == "neg", acq_day = day,
             instrument = instr, patient_day = paste(patient, day, sep = "_"), pct = 100 * count / parent_count,
             stringsAsFactors = FALSE)
fq <- rbind(
  # patient 1: exact neg on week 1; week 2 has no neg but shares the day with the week 1 neg
  tube("p1w1_neg", 1, 1, "neg", "d1", ref, 10, 10000), tube("p1w1_a", 1, 1, "a", "d1", ref, 100, 10000),
  tube("p1w2_a", 1, 2, "a", "d1", ref, 80, 10000),
  # patient 1 week 3, another day, on the other instrument: falls to the patient's reference negs
  tube("p1w3_a", 1, 3, "a", "d2", oth, 90, 10000),
  # patient 2: the exact neg has an empty parent gate, no other neg of the patient -> pooled reference negs
  tube("p2w1_neg", 2, 1, "neg", "d3", ref, 0, 0), tube("p2w1_a", 2, 1, "a", "d3", ref, 50, 10000),
  # patient 3: two negs on the same day (two visits) and a stimulated tube without a week match
  tube("p3w1_neg", 3, 1, "neg", "d4", ref, 20, 10000), tube("p3w2_neg", 3, 2, "neg", "d4", ref, 40, 20000),
  tube("p3w5_a", 3, 5, "a", "d4", ref, 300, 10000),
  # patient 4 on the other instrument only, with no neg anywhere near
  tube("p4w1_a", 4, 1, "a", "d5", oth, 60, 10000))

cat("match_controls\n")
m <- match_controls(fq, controls, ref)
r <- function(fn) m[m$filename == fn, ]
check("adds the neg columns", all(c("neg_count", "neg_parent_count", "neg_pct", "neg_match", "neg_filename", "neg_acq_day", "neg_instrument") %in% names(m)))
check("row count unchanged", nrow(m) == nrow(fq))
check("exact match on patient and week", r("p1w1_a")$neg_match == "exact" && r("p1w1_a")$neg_filename == "p1w1_neg" && r("p1w1_a")$neg_pct == 0.1)
check("exact: neg day and instrument recorded", r("p1w1_a")$neg_acq_day == "d1" && r("p1w1_a")$neg_instrument == ref)
check("fallback same_patient_day", r("p1w2_a")$neg_match == "same_patient_day" && r("p1w2_a")$neg_filename == "p1w1_neg")
check("fallback same_patient_reference_instrument", r("p1w3_a")$neg_match == "same_patient_reference_instrument" && r("p1w3_a")$neg_instrument == ref && r("p1w3_a")$neg_acq_day == "d1")
check("empty-gate neg counts as missing -> pooled_reference_instrument", r("p2w1_a")$neg_match == "pooled_reference_instrument")
p2 <- r("p2w1_a")
check("pooling sums counts over all reference negs", p2$neg_count == 70 && p2$neg_parent_count == 40000 && abs(p2$neg_pct - 100 * 70 / 40000) < 1e-12)
check("pooled across days -> neg_acq_day NA, instrument kept", is.na(p2$neg_acq_day) && p2$neg_instrument == ref)
check("pooled neg_filename is comma-joined", identical(sort(strsplit(p2$neg_filename, ",")[[1]]), c("p1w1_neg", "p3w1_neg", "p3w2_neg")))
p3 <- r("p3w5_a")
check("same_patient_day pools the two negs of that day", p3$neg_match == "same_patient_day" && p3$neg_count == 60 && p3$neg_parent_count == 30000 && p3$neg_acq_day == "d4")
p4 <- r("p4w1_a")
check("other instrument, no patient neg -> pooled reference negs", p4$neg_match == "pooled_reference_instrument" && p4$neg_count == 70)
check("neg rows get NA", all(is.na(m$neg_match[m$is_neg])) && all(is.na(m$neg_pct[m$is_neg])))
m0 <- match_controls(fq, list(match_keys = c("patient", "week"), fallbacks = character(0)), ref)
check("without fallbacks: none", m0$neg_match[m0$filename == "p1w2_a"] == "none" && is.na(m0$neg_pct[m0$filename == "p1w2_a"]) && sum(m0$neg_match == "exact", na.rm = TRUE) == 1)
check("unknown fallback is refused", inherits(try(match_controls(fq, list(match_keys = "patient", fallbacks = "nearest_moon"), ref), silent = TRUE), "try-error"))
check("keys other than patient/week work", match_controls(fq, list(match_keys = "patient_day", fallbacks = character(0)), ref)$neg_match[fq$filename == "p1w2_a"] == "exact")

cat("count_background_flips\n")
# one response, three stimulated rows: two with a neg on the shifted day, one on another day
fb <- data.frame(pct = c(0.30, 0.30, 0.30), neg_pct = c(0.20, 0.05, 0.20), neg_acq_day = c("d1", "d1", "d2"),
                 neg_instrument = c(oth, oth, ref), stringsAsFactors = FALSE)
before <- positivity_call(fb$pct, fb$neg_pct, 0.05, 3)$positive
check("baseline calls: negative, positive, negative", identical(before, c(FALSE, TRUE, FALSE)))
# removing a background effect of +2 logits from the d1 negs lowers their neg_pct ~7-fold: row 1 becomes positive
fl <- count_background_flips(fb, effect = 2, which_neg = fb$neg_acq_day == "d1", 0.05, 3)
check("shifted background flips the row whose neg was inflated", fl == 1)
check("zero effect flips nothing", count_background_flips(fb, 0, fb$neg_acq_day == "d1", 0.05, 3) == 0)
fl2 <- count_background_flips(fb, effect = -5, which_neg = fb$neg_acq_day == "d1", 0.05, 3)
check("negative effect (neg too low) turns the positive call negative, other day untouched", fl2 == 1)
check("rows whose neg is on another day are never touched", count_background_flips(fb, 5, fb$neg_acq_day == "d9", 0.05, 3) == 0)
check("which_neg by instrument works the same way", count_background_flips(fb, 2, fb$neg_instrument == oth, 0.05, 3) == 1)
check("default thresholds are the plan's 0.05 percent and 3-fold", count_background_flips(fb, 2, fb$neg_acq_day == "d1") == 1)

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
