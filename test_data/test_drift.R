#!/usr/bin/env Rscript
# test_data/test_drift.R — unit tests for the base-R drift modules
# Covers src/harmonization/drift/{drift_metadata,drift_panel,drift_stats,
# drift_models}.R. Needs only base R plus test_data/fcs_fixtures.R.
#   Rscript test_data/test_drift.R
cat("=== MTB drift module tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
for (m in list.files(file.path(base_dir, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)
source(file.path(base_dir, "test_data", "fcs_fixtures.R"))

n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") }
}
tmp <- tempfile("drift_"); dir.create(tmp)

# --- file names ---------------------------------------------------------------
cfg <- list(
  tubes = list(filename_regex = "^MTB(?<patient>\\d+)_(?<prior_vacc>\\d+)_V(?<vaccine>\\d+)W(?<week>\\d+)_(?<cocktail>.+)$",
               integer_fields = c("patient", "prior_vacc", "vaccine", "week"), neg_regex = "neg$"),
  instrument_labels = list(Symphony = "FACSymphonyA1", Canto = "FACSCantoII"))
cat("parse_filename\n")
p <- parse_filename("MTB012_2_V4W70_215008.fcs", cfg$tubes)
check("patient/prior_vacc/vaccine/week/cocktail", p$patient == 12 && p$prior_vacc == 2 && p$vaccine == 4 && p$week == 70 && p$cocktail == "215008" && !p$is_neg)
p <- parse_filename("MTB012_2_V4W110_1_neg.fcs", cfg$tubes)
check("neg with suffix", p$is_neg && p$cocktail == "neg")
check("unparsable -> NA", is.na(parse_filename("foo.fcs", cfg$tubes)$patient))
check("unnamed groups are refused", inherits(try(parse_filename("x.fcs", list(filename_regex = "^(\\d+)$", neg_regex = "neg")), silent = TRUE), "try-error"))
p <- parse_filename("S7_stimA.fcs", list(filename_regex = "^(?<patient>S\\d+)_(?<cocktail>.+)$", integer_fields = list(), neg_regex = "^unstim$"))
check("other cohort: character patient, no neg", p$patient == "S7" && p$cocktail == "stimA" && !p$is_neg)

cat("parse_fcs_date / instrument_label\n")
check("Diva date", parse_fcs_date("04-JUN-2025") == as.Date("2025-06-04"))
check("bad date -> NA", is.na(parse_fcs_date("2025/06/04")))
check("instrument labels", identical(instrument_label(c("FACSCantoII", "FACSymphonyA1", NA, "CytoFLEX S"), cfg$instrument_labels), c("FACSCantoII", "FACSymphonyA1", NA, "CytoFLEXS")))

# --- keywords from a byte-level fixture ------------------------------------------
cat("fcs_technical_keywords / build_tube_table\n")
f1 <- file.path(tmp, "MTB011_1_V3W315_234216.fcs")
fcs_fixture(f1, par = 3L, tot = 4L, keywords = list("$DATE" = "04-JUN-2025", "$CYT" = "FACSCantoII", "CREATOR" = "BD FACSDiva Software Version 6.1.2",
  "$P1N" = "FSC-A", "$P2N" = "PE-A", "$P2S" = "INFy", "$P3N" = "APC-Cy7-A", "SPILL" = "1,PE-A,1"))
f2 <- file.path(tmp, "MTB011_1_V3W315_neg.fcs")
fcs_fixture(f2, par = 3L, tot = 4L, keywords = list("$DATE" = "07-MAR-2024", "$CYT" = "FACSymphonyA1", "CREATOR" = "BD FACSDiva Software Version 9.0.2",
  "$P1N" = "FSC-A", "$P2N" = "PE-A", "$P3N" = "APC-Cy7-A"))
k <- fcs_technical_keywords(f1, cfg$instrument_labels)
check("date/instrument/software/par", k$acq_date == "2025-06-04" && k$instrument == "FACSCantoII" && k$software == "6.1.2" && k$n_par == 3 && k$has_spill)
tt <- build_tube_table(tmp, cfg)
check("two rows, neg flagged, patient_day built", nrow(tt) == 2 && sum(tt$is_neg) == 1 && all(grepl("^11_", tt$patient_day)))
sheet <- data.frame(filename = c("MTB011_1_V3W315_234216.fcs", "MTB011_1_V3W315_neg.fcs"), patient = "P1", cocktail = c("cock", "neg"), week = c(3, 3))
ts <- build_tube_table(tmp, cfg, sheet)
check("sample sheet overrides file names", all(ts$patient == "P1") && sum(ts$is_neg) == 1 && ts$cocktail[ts$is_neg] == "neg" && all(grepl("^P1_", ts$patient_day)))
ts2 <- build_tube_table(tmp, cfg, data.frame(filename = "MTB011_1_V3W315_neg.fcs", patient = 1, cocktail = "x", is_neg = TRUE))
check("file missing from sheet -> NA row, explicit is_neg respected", sum(is.na(ts2$patient)) == 1 && ts2$cocktail[ts2$is_neg %in% TRUE] == "neg")
check("sheet without required fields is refused", inherits(try(build_tube_table(tmp, cfg, data.frame(filename = "a.fcs", cocktail = "neg")), silent = TRUE), "try-error"))
br <- bridge_report(tt)
check("bridge report: 1 patient, 2 days, both instruments", nrow(br) == 1 && br$n_days == 2 && br$instr_bridge && br$day_bridge)

# --- panel ---------------------------------------------------------------------
cat("panel reconciliation\n")
pats <- list(FSC_A = "^FSC-A$", IFNg = c("IFN", "INF", "^PE-A$"), CD4 = c("^CD4$", "^APC-Cy7-A$"), IL2 = c("IL-?2", "^APC-A$"))
pm <- build_panel_map(tt, pats)
m1 <- pm$map[pm$map$filename == basename(f1), ]; m2 <- pm$map[pm$map$filename == basename(f2), ]
check("marker name INFy -> IFNg", m1$canonical[m1$channel == "PE-A"] == "IFNg")
check("blank $PnS falls back to channel", m2$canonical[m2$channel == "APC-Cy7-A"] == "CD4" && m2$canonical[m2$channel == "PE-A"] == "IFNg")
g <- panel_gaps(pm$map, c("FSC_A", "IFNg", "CD4", "IL2"))
check("IL2 missing in both", nrow(g) == 2 && all(g$missing == "IL2"))
check("one signature", length(unique(pm$signature$panel_signature)) == 1)

# --- statistics -----------------------------------------------------------------
cat("statistics\n")
set.seed(1)
check("EMD of identical samples ~ 0", emd_1d(rnorm(5000), rnorm(5000)) < 0.05)
check("EMD of shift 1 ~ 1", abs(emd_1d(rnorm(20000), rnorm(20000) + 1) - 1) < 0.05)
x <- c(rnorm(3000, 0, 0.3), rnorm(3000, 3, 0.3))
v <- bimodal_valley(x); check("valley between the modes", !is.na(v) && v > 1 && v < 2)
check("unimodal -> NA valley", is.na(bimodal_valley(rnorm(3000))))
s <- marker_summary(x, valley = TRUE); check("frac_pos ~ 0.5", abs(s$frac_pos - 0.5) < 0.03 && s$n == 6000)
pc <- positivity_call(c(0.20, 0.10, 0.06, 0.02), c(0.05, 0.05, 0.01, 0.00))
check("positivity rule", identical(pc$positive, c(TRUE, FALSE, TRUE, FALSE)))
check("replicate_sd pools within group", abs(replicate_sd(c(1, 2, 3, 10, 12), c("a", "a", "a", "b", "b")) - sqrt(mean(c(1, 2)))) < 1e-9)
dd <- drift_decision(rnorm(4000, 1, 0.1), rep_sd = 0.5, flips = 3)
check("decision: detected, correct", dd$verdict == "detected, correct" && dd$p_exceeds_noise > 0.99)
dd <- drift_decision(rnorm(4000, 0, 0.1), rep_sd = 0.5, flips = 0)
check("decision: not detected", dd$verdict == "not detected")
dd <- drift_decision(rnorm(4000, 1, 0.1), rep_sd = 0.5, flips = 0)
check("decision: detected, not clinically relevant", dd$verdict == "detected, not clinically relevant")

# --- model frames and flips -----------------------------------------------------------
cat("model frames\n")
fq <- data.frame(filename = c("a", "b", "n"), parent = "CD4", gate = "IFNg+", count = c(60, 30, 5), parent_count = 10000,
  patient = 1, prior_vacc = 1, week = 6, cocktail = c("x", "y", "neg"), is_neg = c(FALSE, FALSE, TRUE), acq_day = "2024-01-01",
  instrument = c("FACSymphonyA1", "FACSCantoII", "FACSCantoII"), patient_day = "1_2024-01-01", pct = c(0.6, 0.3, 0.05),
  neg_pct = 0.05, stringsAsFactors = FALSE)
fb <- layer_b_frame(fq)
check("layer_b_frame drops negs, adds response", nrow(fb) == 2 && all(fb$response == "CD4:IFNg+"))
fl <- count_positivity_flips(fb, effect = 3, which = fb$instrument == "FACSymphonyA1", 0.05, 3)
check("large removed effect flips the Symphony call", fl == 1)
fl0 <- count_positivity_flips(fb, effect = 0, which = fb$instrument == "FACSymphonyA1", 0.05, 3)
check("zero effect flips nothing", fl0 == 0)
e1 <- effects_long("m1", "A", "CD4:median", rnorm(100, 0.5, 0.1), matrix(rnorm(200, 0, 0.1), 100, dimnames = list(NULL, c("d1", "d2"))), "FACSymphonyA1")
e2 <- effects_long("m2", "A", "CD4:median", rnorm(100, 0.5, 0.1), matrix(rnorm(200, 0, 0.1), 100, dimnames = list(NULL, c("d1", "d2"))), "FACSymphonyA1")
cmp <- compare_effects(rbind(e1, e2), "m1", "m2")
check("compare_effects: 3 rows, all agree", nrow(cmp) == 3 && all(cmp$agree))
e3 <- e2; e3$model <- "m3"; e3$value[e3$term == "instrument"] <- e3$value[e3$term == "instrument"] + 5
cmp3 <- compare_effects(rbind(e1, e3), "m1", "m3")
check("compare_effects flags a disagreeing instrument effect", !cmp3$agree[cmp3$term == "instrument"])

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
