#!/usr/bin/env Rscript
# test_data/test_drift_design.R — unit tests for the design-driven model builders in
# src/harmonization/drift/drift_models.R: design_frame, design_formula, the multi-instrument
# effects_long, covariates_long and the layer BG frame readers. Base R only.
#   Rscript test_data/test_drift_design.R
cat("=== MTB drift design tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_design.R"
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
tmp <- tempfile("design_"); dir.create(tmp)

cfg <- list(
  design = list(unit = "patient", technical = c("acq_day", "instrument"), protected_fixed = "prior_vacc",
                time = "week", stimulation = "cocktail"),
  technical_markers = c("FSC_A", "CD4"), cytokine_markers = c("IFNg", "IL2"), reference_instrument = "Canto")

# --- design_frame -------------------------------------------------------------------
cat("design_frame\n")
df <- data.frame(patient = c(1, 1, 2, 2), acq_day = c("d1", "d2", "d1", "d3"), instrument = c("Symph", "Canto", "Canto", "Aurora"),
                 prior_vacc = c(1, 1, 2, 2), week = c(0, 52, 0, 52), cocktail = c("a", "b", "a", "neg"), stringsAsFactors = FALSE)
fr <- design_frame(df, cfg, ref = "Canto")
check("unit, day, instrument, protected and stimulation are factors",
      all(vapply(fr[c("patient", "acq_day", "instrument", "prior_vacc", "cocktail")], is.factor, logical(1))))
check("instrument relevelled to the reference", levels(fr$instrument)[1] == "Canto" && nlevels(fr$instrument) == 3)
check("week_c centred and in years", isTRUE(all.equal(fr$week_c, c(-0.5, 0.5, -0.5, 0.5))))
fr0 <- design_frame(df[, setdiff(names(df), "cocktail")], cfg)
check("stimulation column optional, no ref keeps alphabetical levels", !"cocktail" %in% names(fr0) && levels(fr0$instrument)[1] == "Aurora")
cfg_min <- list(design = default_design())
fr1 <- design_frame(df, cfg_min)
check("default design: no week_c, prior_vacc left numeric", !"week_c" %in% names(fr1) && is.numeric(fr1$prior_vacc))
check("missing design column is an error", inherits(try(design_frame(df[, -1], cfg), silent = TRUE), "try-error"))
check("design_columns lists every named column", identical(design_columns(cfg), c("patient", "acq_day", "instrument", "prior_vacc", "week", "cocktail")))

# --- design_formula -------------------------------------------------------------------
cat("design_formula\n")
check("technical / brms", fstr(design_formula(cfg, "technical", "brms")) == "y ~ 1 + instrument + (1 | patient) + (1 | acq_day)")
check("technical / lme4 with lhs", fstr(design_formula(cfg, "technical", "lme4", "value")) == "value ~ 1 + instrument + (1 | patient) + (1 | acq_day)")
check("background_counts / brms", fstr(design_formula(cfg, "background_counts", "brms")) ==
      "count | trials(parent_count) ~ 1 + instrument + prior_vacc + week_c + (1 + week_c | patient) + (1 | acq_day)")
check("background_counts / lme4", fstr(design_formula(cfg, "background_counts", "lme4")) ==
      "cbind(count, parent_count - count) ~ 1 + instrument + prior_vacc + week_c + (1 + week_c | patient) + (1 | acq_day) + (1 | filename)")
check("background_intensity / brms", fstr(design_formula(cfg, "background_intensity", "brms")) ==
      "y ~ 1 + instrument + prior_vacc + week_c + (1 | patient) + (1 | acq_day)")
check("background_intensity / lme4", fstr(design_formula(cfg, "background_intensity", "lme4", "value")) ==
      "value ~ 1 + instrument + prior_vacc + week_c + (1 | patient) + (1 | acq_day)")
check("endpoint / brms", fstr(design_formula(cfg, "endpoint", "brms")) ==
      "count | trials(parent_count) ~ 1 + instrument + prior_vacc + week_c + (1 + week_c | patient) + (1 | acq_day) + (1 | cocktail)")
check("endpoint / lme4", fstr(design_formula(cfg, "endpoint", "lme4")) ==
      "cbind(count, parent_count - count) ~ 1 + instrument + prior_vacc + week_c + (1 + week_c | patient) + (1 | acq_day) + (1 | cocktail) + (1 | filename)")
check("default design: endpoint without biology", fstr(design_formula(cfg_min, "endpoint", "brms")) ==
      "count | trials(parent_count) ~ 1 + instrument + (1 | patient) + (1 | acq_day) + (1 | cocktail)")
cfg2 <- cfg; cfg2$design$protected_fixed <- c("prior_vacc", "sex"); cfg2$design$technical <- c("acq_day", "run", "instrument")
check("several protected and technical terms", fstr(design_formula(cfg2, "background_intensity", "brms")) ==
      "y ~ 1 + instrument + prior_vacc + sex + week_c + (1 | patient) + (1 | acq_day) + (1 | run)")
check("unknown class refused", inherits(try(design_formula(cfg, "nonsense"), silent = TRUE), "try-error"))
check("design_fixed_vars per class", identical(design_fixed_vars(cfg, "technical"), "instrument") &&
      identical(design_fixed_vars(cfg, "background_intensity"), c("instrument", "prior_vacc", "week_c")) &&
      identical(design_fixed_vars(cfg, "endpoint"), c("instrument", "prior_vacc", "week_c")))
ft <- fixed_terms_table(fr, c("instrument", "prior_vacc", "week_c"))
check("fixed_terms_table: K-1 instrument rows, factor levels, numeric term",
      identical(ft$coef, c("instrumentAurora", "instrumentSymph", "prior_vacc2", "week_c")) &&
      identical(ft$term, c("instrument", "instrument", "prior_vacc", "week_c")) && is.na(ft$level[4]))

# --- effects_long with three instruments, covariates_long ---------------------------------
cat("effects_long / covariates_long\n")
set.seed(2)
inst <- cbind(Aurora = rnorm(50, 1), Symph = rnorm(50, -1))
days <- matrix(rnorm(100), 50, dimnames = list(NULL, c("d1", "d2")))
e <- effects_long("m", "A", "FSC_A:median", inst, days)
check("matrix: one instrument block per column", identical(unique(e$level[e$term == "instrument"]), c("Aurora", "Symph")) &&
      sum(e$term == "instrument") == 100 && sum(e$term == "day") == 100)
check("column values kept", isTRUE(all.equal(e$value[e$level == "Symph"], inst[, "Symph"])))
el <- effects_long("m", "A", "FSC_A:median", list(Aurora = inst[, 1], Symph = inst[, 2]), days)
check("named list equals matrix", isTRUE(all.equal(el, e)))
ev <- effects_long("m", "A", "x", inst[, 1], days, "Aurora")
check("single vector + instrument_level (old call) still works", sum(ev$term == "instrument") == 50 && all(ev$level[ev$term == "instrument"] == "Aurora"))
e0 <- effects_long("m", "A", "x", inst[, 0, drop = FALSE], days)
check("zero instrument columns: day rows only", all(e0$term == "day") && nrow(e0) == 100)
check("columns in the shared order", identical(names(e), c("model", "layer", "response", "term", "level", "draw", "value")))
cv <- covariates_long("m", "BG", "CD4:IFNg", list(prior_vacc = matrix(rnorm(50, 2), 50, dimnames = list(NULL, "2")), week_c = rnorm(50, -3)))
check("covariates_long: one row per coefficient", nrow(cv) == 2 && identical(cv$covariate, c("prior_vacc", "week_c")) &&
      cv$level[1] == "2" && is.na(cv$level[2]) && abs(cv$mean[2] + 3) < 0.5 && cv$lo[1] < cv$mean[1] && cv$mean[1] < cv$hi[1])
check("covariates_long: empty list gives empty table with the columns",
      nrow(covariates_long("m", "A", "x", list())) == 0 && identical(names(covariates_long("m", "A", "x", list())), c("model", "layer", "response", "covariate", "level", "mean", "lo", "hi")))
dr <- cbind(b_instrumentAurora = inst[, 1], b_instrumentSymph = inst[, 2], b_prior_vacc2 = rnorm(50), b_week_c = rnorm(50), other = 1)
ft$coef <- paste0("b_", ft$coef)
sp <- split_fixed_draws(dr, ft)
check("split_fixed_draws: instrument matrix named by level, covariates by term",
      identical(colnames(sp$instrument), c("Aurora", "Symph")) && identical(names(sp$covariates), c("prior_vacc", "week_c")) &&
      is.matrix(sp$covariates$prior_vacc) && is.numeric(sp$covariates$week_c) && is.null(dim(sp$covariates$week_c)))

# --- layer BG readers on a tiny folder -------------------------------------------------------
cat("read_layer_bg_counts / read_layer_bg_intensity\n")
tubes <- data.frame(filename = c("p1_neg.fcs", "p1_stim.fcs", "p2_neg.fcs", "p3_neg.fcs"), patient = c(1, 1, 2, 3), prior_vacc = c(1, 1, 2, 2),
                    week = c(0, 0, 52, 104), cocktail = c("neg", "stim", "neg", "neg"), is_neg = c(TRUE, FALSE, TRUE, TRUE),
                    acq_day = c("d1", "d1", "d2", "d3"), instrument = c("Canto", "Canto", "Symph", "Aurora"),
                    patient_day = c("1_d1", "1_d1", "2_d2", "3_d3"), stringsAsFactors = FALSE)
gf <- merge(expand.grid(filename = tubes$filename, parent = c("CD4", "CD8"), gate = "IFNg", stringsAsFactors = FALSE), tubes, by = "filename")
gf$count <- 5; gf$parent_count <- 1000; gf$parent_count[gf$filename == "p3_neg.fcs" & gf$parent == "CD8"] <- 0
gf$pct <- 100 * gf$count / pmax(gf$parent_count, 1); gf$neg_pct <- 0.5
write.csv(gf, file.path(tmp, "gated_frequencies.csv"), row.names = FALSE)
bc <- read_layer_bg_counts(tmp, cfg, ref = "Canto")
check("negs only, empty parent dropped", nrow(bc) == 5 && all(!grepl("stim", bc$filename)) && !any(bc$parent_count == 0))
check("response key and design factors", all(bc$response %in% c("CD4:IFNg", "CD8:IFNg")) && is.factor(bc$prior_vacc) &&
      levels(bc$instrument)[1] == "Canto" && "week_c" %in% names(bc))
check("expected columns present", all(c("filename", "patient", "prior_vacc", "week", "week_c", "acq_day", "instrument", "patient_day",
                                         "parent", "gate", "count", "parent_count", "pct", "response") %in% names(bc)))
ns <- merge(expand.grid(filename = tubes$filename, marker = c("FSC_A", "CD4", "IFNg", "IL2"), stringsAsFactors = FALSE), data.frame(n = 100), by = NULL)
ns$median <- seq_len(nrow(ns)); ns$mad <- 1; ns$q10 <- 0; ns$q90 <- 2; ns$spread <- 2; ns$valley <- NA; ns$frac_pos <- NA
write.csv(ns, file.path(tmp, "neg_summaries.csv"), row.names = FALSE)
bi <- read_layer_bg_intensity(tmp, cfg, tubes, ref = "Canto")
check("cytokine medians of the negs only", nrow(bi) == 6 && all(bi$marker %in% c("IFNg", "IL2")) && all(bi$statistic == "median") && all(bi$is_neg))
check("layer A layout plus response key", all(c("filename", "marker", "statistic", "value", "patient", "acq_day", "instrument", "week_c", "response") %in% names(bi)) &&
      all(bi$response %in% c("IFNg:median", "IL2:median")))
check("values are the medians", isTRUE(all.equal(sort(bi$value), sort(ns$median[ns$marker %in% c("IFNg", "IL2") & ns$filename != "p1_stim.fcs"]))))
# layer B frame still behaves as before and takes the config
fb <- layer_b_frame(gf, cfg, ref = "Canto")
check("layer_b_frame: stimulated rows, response, relevelled instrument", nrow(fb) == 2 && all(fb$response %in% c("CD4:IFNg", "CD8:IFNg")) && levels(fb$instrument)[1] == "Canto")
fb_old <- layer_b_frame(gf)
check("layer_b_frame without cfg keeps prior_vacc factor and week_c", is.factor(fb_old$prior_vacc) && "week_c" %in% names(fb_old))

unlink(tmp, recursive = TRUE)
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
