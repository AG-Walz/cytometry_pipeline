#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_cytokine_combinations.R — pure parts of the cytokine
# combination bar plots
# =============================================================================
# Covers the name parsing, the Boolean combination counting, the background
# subtraction and the gate-path resolution behind
# src/analysis/cytokine_combinations.R. Needs neither CytoML nor ggplot2 nor a
# real workspace, so it runs on base R; what it cannot cover is the GatingSet
# side and the rendering, which need an analysis host.
#
# Run from the project root:
#   Rscript test_data/test_cytokine_combinations.R
# =============================================================================

cat("=== Cytokine combination tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_cytokine_combinations.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(base_dir, "src", "loading", "flowjo_workspace.R"))
source(file.path(base_dir, "src", "loading", "sample_annotations.R"))
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
source(file.path(base_dir, "src", "harmonization", "panel_comparison.R"))
source(file.path(base_dir, "src", "visualization", "umap_from_flowjo.R"))
source(file.path(base_dir, "src", "analysis", "cytokine_combinations.R"))

.pass <- 0L
.fail <- 0L
ok <- function(label, cond) {
  if (isTRUE(cond)) { .pass <<- .pass + 1L; cat("  PASS  ", label, "\n", sep = "") }
  else              { .fail <<- .fail + 1L; cat("  FAIL  ", label, "\n", sep = "") }
}
eq <- function(label, actual, expected) {
  same <- isTRUE(all.equal(actual, expected))
  if (!same) {
    cat("  FAIL  ", label, "\n        expected: ",
        paste(format(expected), collapse = " | "), "\n        actual:   ",
        paste(format(actual), collapse = " | "), "\n", sep = "")
    .fail <<- .fail + 1L
  } else { .pass <<- .pass + 1L; cat("  PASS  ", label, "\n", sep = "") }
}

# --- 1. sample names ---------------------------------------------------------
cat("\n[1] parse_mtb_sample\n")
nm <- c("MTB59_1_V0W0_224065.fcs_0", "MTB059_1_V1W5_neg.fcs",
        "MTB012_2_V4W110_1_neg.fcs", "MTB016_2_V2W15_194335_018.fcs",
        "AML-VAC 01-07 2WpostV2_AML Pool_001.fcs")
m <- suppressWarnings(parse_mtb_sample(nm))
eq("patient is zero-padded", m$patient, c("MTB059", "MTB059", "MTB012", "MTB016", NA))
eq("week and vaccine are integers", m$week, c(0L, 5L, 110L, 15L, NA))
eq("timepoint label", m$timepoint, c("V0 W0", "V1 W5", "V4 W110", "V2 W15", NA))
eq("cocktail: neg, 1_neg -> neg, _NNN stripped",
   m$cocktail, c("224065", "neg", "neg", "194335", NA))
eq("sample name loses the .fcs and cytoset suffix", m$sample[1], "MTB59_1_V0W0_224065")
ok("unparsable name warns",
   inherits(tryCatch(parse_mtb_sample(nm[5]), warning = function(w) w), "warning"))

# --- 2. combinations ---------------------------------------------------------
cat("\n[2] combination labels and counts\n")
mk <- c("CD107a", "IFNy", "IL-2", "TNF")
lv <- combination_levels(mk)
eq("16 levels", length(lv), 16L)
eq("none first, then singles in marker order", lv[1:5],
   c("none", "CD107a+", "IFNy+", "IL-2+", "TNF+"))
eq("all-positive last", lv[16], "CD107a+ IFNy+ IL-2+ TNF+")
mem <- rbind(c(FALSE, FALSE, FALSE, FALSE),
             c(FALSE, TRUE,  TRUE,  FALSE),
             c(FALSE, TRUE,  TRUE,  FALSE),
             c(TRUE,  FALSE, FALSE, TRUE),
             c(TRUE,  TRUE,  TRUE,  TRUE))
lab <- combination_labels(mem, mk)
eq("per-cell labels", as.character(lab),
   c("none", "IFNy+ IL-2+", "IFNy+ IL-2+", "CD107a+ TNF+", "CD107a+ IFNy+ IL-2+ TNF+"))
cnt <- combination_counts(mem, mk)
eq("every level present", nrow(cnt), 16L)
eq("counts sum to the number of cells", sum(cnt$count), nrow(mem))
eq("IFNy+ IL-2+ counted twice", cnt$count[cnt$combination == "IFNy+ IL-2+"], 2L)
eq("empty membership gives zero counts", sum(combination_counts(matrix(FALSE, 0, 4), mk)$count), 0L)
eq("NA membership counts as negative",
   as.character(combination_labels(matrix(c(NA, TRUE, NA, NA), 1), mk)), "IFNy+")

# --- 3. percent and background -----------------------------------------------
cat("\n[3] percent and background subtraction\n")
tbl <- expand.grid(combination = c("none", "IFNy+", "TNF+"),
                   cocktail = c("neg", "224065", "999999"), stringsAsFactors = FALSE)
tbl$workspace <- "MTB059_1"; tbl$patient <- "MTB059"; tbl$timepoint <- "V1 W5"
tbl$parent <- "CD4"; tbl$sample <- paste0("MTB59_1_V1W5_", tbl$cocktail)
tbl$count <- c(90L, 5L, 5L,   80L, 15L, 5L,   97L, 1L, 2L)
tbl <- combination_percent(tbl)
eq("percent closes at 100 per sample",
   as.numeric(tapply(tbl$percent, tbl$sample, sum)), c(100, 100, 100))
bg <- suppressWarnings(subtract_background(tbl))
ok("neg rows are gone", !any(bg$cocktail == "neg"))
eq("IFNy+ minus background", bg$percent[bg$cocktail == "224065" & bg$combination == "IFNy+"], 10)
eq("floored at zero", bg$percent[bg$cocktail == "999999" & bg$combination == "IFNy+"], 0)
eq("none closes the stack",
   as.numeric(tapply(bg$percent, bg$sample, sum)), c(100, 100))
tbl2 <- tbl; tbl2$timepoint[tbl2$cocktail == "999999"] <- "V2 W11"
w <- tryCatch(subtract_background(tbl2), warning = function(w) w)
ok("a cocktail without neg at its timepoint warns", inherits(w, "warning"))
eq("...and is dropped", unique(suppressWarnings(subtract_background(tbl2))$cocktail), "224065")

# --- 4. gate paths -----------------------------------------------------------
cat("\n[4] resolve_gate_paths\n")
pops <- c("/Lymphocytes", "/Lymphocytes/Single Cells", "/Lymphocytes/Single Cells/Live cells",
          "/Lymphocytes/Single Cells/Live cells/CD4",
          "/Lymphocytes/Single Cells/Live cells/CD4/CD107a",
          "/Lymphocytes/Single Cells/Live cells/CD4/INFy",
          "/Lymphocytes/Single Cells/Live cells/CD4/IL-2",
          "/Lymphocytes/Single Cells/Live cells/CD4/TNF",
          "/Lymphocytes/Single Cells/Live cells/CD8",
          "/Lymphocytes/Single Cells/Live cells/CD8/TNF")
r <- resolve_gate_paths(pops, "CD4", mk)
eq("parent path", r$parent, "/Lymphocytes/Single Cells/Live cells/CD4")
eq("INFy is accepted as IFNy", unname(r$markers["IFNy"]),
   "/Lymphocytes/Single Cells/Live cells/CD4/INFy")
eq("children only under the parent", unname(r$markers["TNF"]),
   "/Lymphocytes/Single Cells/Live cells/CD4/TNF")
r8 <- suppressWarnings(resolve_gate_paths(pops, "CD8", mk))
ok("missing child is NA", is.na(r8$markers["IFNy"]) && !is.na(r8$markers["TNF"]))
ok("missing child warns",
   inherits(tryCatch(resolve_gate_paths(pops, "CD8", mk), warning = function(w) w), "warning"))
ok("ambiguous parent is an error",
   inherits(tryCatch(resolve_gate_paths(pops, "TNF", mk), error = function(e) e), "error"))
eq("gate_alias leaves other names alone", gate_alias(c("INFy", "TNF")), c("IFNy", "TNF"))

# --- 5. palette and level recovery ------------------------------------------
cat("\n[5] palette\n")
pal <- combination_palette(lv)
eq("one colour per level", length(pal), 16L)
eq("none is grey", unname(pal["none"]), "grey85")
eq("distinct colours", length(unique(pal)), 16L)
eq("levels recovered from labels", combination_levels_from(c("IFNy+ IL-2+", "none", "TNF+")),
   c("none", "TNF+", "IFNy+ IL-2+"))

cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")
if (.fail > 0) quit(status = 1)
