#!/usr/bin/env Rscript
# test_data/test_drift_biology.R — unit tests for src/harmonization/drift/drift_biology.R
# Base R only: the glmer fits are not run, the formulas and the table assembly are tested.
#   Rscript test_data/test_drift_biology.R
cat("=== MTB drift biology tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_biology.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
for (m in list.files(file.path(base_dir, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)

n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") }
}
cfg <- list(design = list(unit = "patient", technical = c("acq_day", "instrument"), protected_fixed = "prior_vacc",
                          time = "week", stimulation = "cocktail"))

# --- a small layer B frame -------------------------------------------------------------------
cat("layer B frame\n")
freqs <- data.frame(filename = paste0("f", 1:8), parent = "CD4", gate = "IFNg", count = c(5, 9, 3, 12, 7, 6, 4, 10),
                    parent_count = 1000, patient = rep(1:4, each = 2), prior_vacc = rep(1:2, 4), week = rep(c(0, 52), 4),
                    cocktail = "A", is_neg = FALSE, acq_day = rep(c("d1", "d2"), each = 4),
                    instrument = rep(c("Canto", "Symph"), 4), neg_pct = 0.1, stringsAsFactors = FALSE)
fb <- layer_b_frame(freqs, cfg, "Canto")
check("frame has week_c and factor technical terms", "week_c" %in% names(fb) && is.factor(fb$acq_day) && levels(fb$instrument)[1] == "Canto")

# --- formulas ---------------------------------------------------------------------------------
cat("biology_formula\n")
rhs <- function(f) attr(terms(f, keep.order = TRUE), "term.labels")
f_full <- biology_formula(cfg)
check("full keeps the technical terms and adds the fixed week slope",
      all(c("instrument", "prior_vacc", "week_c", "1 + week_c | patient", "1 | acq_day", "1 | cocktail", "1 | filename") %in% rhs(f_full)))
check("full response is the lme4 count form", deparse1(f_full[[2]]) == "cbind(count, parent_count - count)")
f_nt <- biology_formula(cfg, technical = FALSE)
check("no_technical lacks instrument and acq_day", !any(c("instrument", "acq_day") %in% all.vars(f_nt)))
check("no_technical keeps week_c, prior_vacc, patient, cocktail, filename",
      all(c("week_c", "prior_vacc", "patient", "cocktail", "filename") %in% all.vars(f_nt)) && "week_c" %in% rhs(f_nt))
check("biology_cfg does not touch the other design fields",
      length(biology_cfg(cfg)$design$technical) == 0 && identical(biology_cfg(cfg)$design$protected_fixed, "prior_vacc"))
check("no time field is refused", inherits(try(biology_formula(list(design = list(unit = "patient", technical = "acq_day"))), silent = TRUE), "try-error"))
check("formula terms exist in the frame", all(setdiff(all.vars(f_full), c("count", "parent_count")) %in% names(fb)))

# --- table assembly from fake fit records ---------------------------------------------------------
cat("tables\n")
check("slope_of(NULL) is NA", all(is.na(slope_of(NULL))) && length(slope_of(NULL)) == 2)
recs <- list(list(response = "CD4:IFNg", model = "full", slope = c(0.4, 0.1), n = 8L),
             list(response = "CD4:IFNg", model = "no_technical", slope = c(0.5, 0.08), n = 8L),
             list(response = "CD8:IL2", model = "full", slope = c(-0.2, 0.2), n = 8L),
             list(response = "CD8:IL2", model = "no_technical", slope = c(NA, NA), n = 8L))
s <- sensitivity_table(recs)
check("sensitivity columns and rows", identical(names(s), c("response", "model", "week_slope", "se", "n")) && nrow(s) == 4 && s$week_slope[2] == 0.5)
check("empty sensitivity table keeps the columns", identical(names(sensitivity_table(list())), names(s)) && nrow(sensitivity_table(list())) == 0)
l <- lodo_table(list(list(response = "CD4:IFNg", dropped_day = "d1", slope = c(0.55, 0.12), n_tubes_dropped = 4L),
                     list(response = "CD4:IFNg", dropped_day = "d2", slope = c(0.42, 0.11), n_tubes_dropped = 4L)))
check("lodo columns", identical(names(l), c("response", "dropped_day", "week_slope", "se", "n_tubes_dropped")) && nrow(l) == 2)
uc <- data.frame(`(Intercept)` = c(-5, -4.5), week_c = c(0.3, 0.6), row.names = c("1", "2"), check.names = FALSE)
sj <- subject_slopes(uc, "CD4:IFNg")
check("subject slopes: one row per patient", identical(names(sj), c("response", "patient", "slope")) && identical(sj$patient, c("1", "2")) && sj$slope[2] == 0.6)

cat("sentences\n")
p <- sensitivity_pairs(s)
check("pairs: shift in SE and relative change", nrow(p) == 2 && isTRUE(all.equal(p$shift_se[p$response == "CD4:IFNg"], 1)) &&
      isTRUE(all.equal(p$rel_change[p$response == "CD4:IFNg"], 0.25)) && is.na(p$shift_se[p$response == "CD8:IL2"]))
check("sensitivity sentence names the largest shift with its relative change", grepl("CD4:IFNg.*0.40 to 0.50.*25 percent of the full slope, 1.0 SE", sensitivity_sentence(s)))
s0 <- s; s0$week_slope[1] <- 0.01
check("relative change withheld when the full slope is within one SE of zero", grepl("not meaningful", sensitivity_sentence(s0)) && !grepl("percent", sensitivity_sentence(s0)))
sh <- lodo_shifts(l, s)
check("lodo shifts: d1 flagged (1.5 SE), d2 not", sh$flag[sh$dropped_day == "d1"] && !sh$flag[sh$dropped_day == "d2"])
check("lodo sentence names d1", grepl("^The day whose removal moves a slope most is d1 \\(4 tubes\\) for CD4:IFNg", lodo_sentence(l, s)) && grepl("1 of 2 refits", lodo_sentence(l, s)))
check("sentences degrade without data", grepl("^No endpoint", sensitivity_sentence(sensitivity_table(list()))) && grepl("^No leave", lodo_sentence(lodo_table(list()), s)))

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
