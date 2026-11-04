#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_umap_polyfunctionality.R: polyfunctionality maps of a
# gate-combination run (src/visualization/umap_polyfunctionality.R)
# =============================================================================
# The positive-gate count read from combination labels, the MTB colours and
# sizes, the per-sample table, and, when ggplot2 is installed, the maps and
# plot index written for a small hand-made run of 4 cytokines under CD4/CD8.
#
# Run from the project root:
#   Rscript test_data/test_umap_polyfunctionality.R
# =============================================================================

cat("=== UMAP polyfunctionality tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_umap_polyfunctionality.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

source(file.path(base_dir, "src", "visualization", "umap_polyfunctionality.R"))

# --- tiny test harness -------------------------------------------------------
.pass <- 0L
.fail <- 0L
ok <- function(label, cond) {
  if (isTRUE(cond)) {
    .pass <<- .pass + 1L
    cat("  PASS  ", label, "\n", sep = "")
  } else {
    .fail <<- .fail + 1L
    cat("  FAIL  ", label, "\n", sep = "")
  }
}
eq <- function(label, actual, expected) {
  same <- isTRUE(all.equal(actual, expected))
  if (!same) {
    cat("  FAIL  ", label, "\n",
        "        expected: ", paste(format(expected), collapse = " | "), "\n",
        "        actual:   ", paste(format(actual), collapse = " | "), "\n",
        sep = "")
    .fail <<- .fail + 1L
  } else {
    .pass <<- .pass + 1L
    cat("  PASS  ", label, "\n", sep = "")
  }
}

# --- 1. Positive gates per label ---------------------------------------------
cat("polyfunctional_degree()\n")
lab <- c("CD4+ | none", "CD4+ | IFNy+", "CD8+ | IFNy+ TNFa+",
         "CD8+ | CD107a+ IFNy+ IL-2+ TNFa+", "CD4- CD8-", "CD4+ CD8+",
         "other (rare)", NA)
pf <- polyfunctional_degree(lab)
eq("two tiers: the gates below the lineage are counted, the rest is NA",
   pf$degree, c(0L, 1L, 2L, 4L, NA, NA, NA, NA))
eq("...the gates in the order the most polyfunctional label spells them",
   pf$markers, c("CD107a", "IFNy", "IL-2", "TNFa"))
eq("one tier (no ' | '): the whole label, its negatives as 0",
   polyfunctional_degree(c("IFNy- TNFa-", "IFNy+", "IFNy+ TNFa+", "(none)"))$degree,
   c(0L, 1L, 2L, NA))
eq("three tiers: the last tier, and NA where a label stops early",
   polyfunctional_degree(c("CD3+ | CD4+ | IFNy+", "CD3+ | CD4+ | none",
                           "CD3+ | CD4+"))$degree, c(1L, 0L, NA))
pf <- polyfunctional_degree(c("CD4+ | TNF a+ IFN y+", "CD8+ | CD107a+ TNF a+ IL2+ IFN y+",
                             "CD4+ | none", "CD8+ | HLA-DR+"))
eq("marker names with spaces or a hyphen are one gate (J445: 'TNF a', 'IFN y')",
   pf$markers, c("CD107a", "TNF a", "IL2", "IFN y", "HLA-DR"))
eq("...and counted once each", pf$degree, c(2L, 4L, 0L, 1L))
ok("no positive gate anywhere: NULL (nothing to draw)",
   is.null(polyfunctional_degree(c("CD4+ | none", "CD4- CD8-"))))

# --- 2. Colours, sizes, table ------------------------------------------------
cat("\npolyfunctionality_colours(), _sizes(), _table()\n")
eq("four gates: the MTB cohort's five blues",
   polyfunctionality_colours(4L),
   c("#DCE9F6", "#9CC3E5", "#4A8FC9", "#0B4F8A", "#062F52"))
ok("other counts: k + 1 colours from the same light to dark ramp",
   length(polyfunctionality_colours(2L)) == 3L &&
     polyfunctionality_colours(2L)[1] == "#DCE9F6" &&
     polyfunctionality_colours(2L)[3] == "#062F52")
eq("sizes 1x to 3x the base, as the MTB maps (0.4 -> 0.4 ... 1.2)",
   polyfunctionality_sizes(4L, 0.4), c(0.4, 0.6, 0.8, 1.0, 1.2))
cells_t <- data.frame(sample = c("a", "a", "a", "b", "b"),
                      weight = c(1, 1, 2, 1, 1),
                      annot_timepoint = c("W1", "W1", "W1", "W5", "W5"))
tab <- polyfunctionality_table(cells_t, c(0L, 1L, 1L, NA, 2L), 2L, "timepoint")
eq("table: weighted cells per sample and count", tab$cells[tab$sample == "a"], c(1, 3, 0))
eq("...percent of the sample's counted cells", tab$percent[tab$sample == "a"],
   c(25, 75, 0))
eq("...a sample's uncounted cells are not in its denominator",
   tab$percent[tab$sample == "b"], c(0, 0, 100))
eq("...with the annotation columns", unique(tab$timepoint), c("W1", "W5"))

# --- 3. Maps and index -------------------------------------------------------
cat("\nplot_polyfunctionality()\n")
can_draw <- requireNamespace("ggplot2", quietly = TRUE) &&
  isTRUE(tryCatch({
    source(file.path(base_dir, "src", "R_shared", "plot_utils.R"))
    source(file.path(base_dir, "src", "visualization", "umap_style.R"))
    TRUE
  }, error = function(e) FALSE))
if (!can_draw) {
  cat("  SKIP  ggplot2 or the plot helpers are not available\n")
} else {
  set.seed(1)
  n <- 600
  combos <- c("CD4+ | none", "CD4+ | IFNy+", "CD4+ | IFNy+ TNFa+",
              "CD8+ | CD107a+ IFNy+ IL-2+ TNFa+", "CD4- CD8-")
  cells <- data.frame(UMAP1 = stats::rnorm(n), UMAP2 = stats::rnorm(n),
                      sample = rep(c("s1", "s2", "s3"), length.out = n),
                      combination = sample(combos, n, replace = TRUE),
                      annot_timepoint = rep(c("W5", "W12", "W1"), length.out = n),
                      stringsAsFactors = FALSE)
  out <- file.path(tempdir(), "pf_test")
  unlink(out, recursive = TRUE)
  idx <- NULL
  invisible(capture.output(idx <- plot_polyfunctionality(
    cells, population = "Lymphocytes/T cells", output_dir = out,
    split_cols = c("timepoint", "absent"), output_formats = list(png = TRUE))))
  eq("one map of the run, one per level of the split column (unknown column skipped)",
     nrow(idx), 4L)
  eq("...indexed by kind", idx$kind,
     c("polyfunctionality", rep("polyfunctionality_split", 3)))
  eq("...levels in natural order without ordered_levels(): sorted",
     idx$level[-1], sort(c("W5", "W12", "W1")))
  ok("...every PNG written", all(file.exists(file.path(out, paste0(idx$file, ".png")))))
  ok("...and polyfunctionality_per_sample.csv",
     file.exists(file.path(out, "polyfunctionality_per_sample.csv")))
  idx2 <- NULL
  invisible(capture.output(idx2 <- plot_polyfunctionality(
    cells, population = "T cells", output_dir = out,
    gates = c("IFNy", "IL-2", "TNFa", "CD107a"))))
  eq("polyfunctional_gates(): the last tier's markers",
     polyfunctional_gates(c(CD4 = 1L, CD8 = 1L, IFNy = 2L, `IL-2` = 2L)), c("IFNy", "IL-2"))
  ok("gates given: the map is drawn with them", nrow(idx2) == 1L)
  r <- "unset"
  invisible(capture.output(r <- plot_polyfunctionality(
    transform(cells, combination = "CD4+ | none"), "T cells", out)))
  ok("no positive gate in any label: nothing drawn, NULL", is.null(r))
  unlink(out, recursive = TRUE)
}

cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")
if (.fail > 0) quit(status = 1)
