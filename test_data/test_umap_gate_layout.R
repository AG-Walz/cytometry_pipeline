#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_umap_gate_layout.R: gate-anchored UMAP layout tests
# =============================================================================
# Covers the pure functions behind "--layout gates": reading a decision
# (marker, polarity, cut) off each gate, the per-cell calls built from those
# decisions, the combination labels, stratified and balanced sampling, the gap
# features, the metadata filter and level order, and the density grid behind
# the difference maps. When src/visualization/umap_style.R is present, the
# palette and the combination colours are checked as well.
#
# Every input is built by hand so that the right answer is known exactly. The
# gating example mirrors the real MTB tree: a CD4 x CD8 quadrant of four
# rectangles (two of them named " " and "  "), and cytokine x FSC-A rectangles
# under CD4 and CD8 with different cuts, the IFNg gate spelled IFNy in one
# sample and INFy in another. Nothing here needs CytoML or an FCS file, so the
# file runs on any machine with base R; the kNN-based separation metrics run
# only when RANN or FNN is installed.
#
# Run from the project root:
#   Rscript test_data/test_umap_gate_layout.R
# =============================================================================

cat("=== Gate-anchored UMAP layout tests ===\n\n")

.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_umap_gate_layout.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()

source(file.path(base_dir, "src", "visualization", "umap_gate_layout.R"))

# The style file needs ggplot2 only inside its layer constructors, but source
# it defensively: a host without it should still run every other test here.
style_file <- file.path(base_dir, "src", "visualization", "umap_style.R")
has_style <- file.exists(style_file) &&
  isTRUE(tryCatch({ source(style_file); TRUE }, error = function(e) {
    cat("NOTE: could not source umap_style.R (", conditionMessage(e),
        "); its tests are skipped.\n", sep = "")
    FALSE
  }))

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

throws <- function(label, expr, pattern = NULL) {
  msg <- tryCatch({ force(expr); NULL },
                  error = function(e) conditionMessage(e))
  if (is.null(msg)) {
    ok(label, FALSE)
  } else if (!is.null(pattern) && !grepl(pattern, msg, fixed = FALSE)) {
    cat("  FAIL  ", label, "\n        message did not match '", pattern,
        "':\n        ", msg, "\n", sep = "")
    .fail <<- .fail + 1L
  } else {
    ok(label, TRUE)
  }
}

# Like throws(), for a warning; returns the value so the result can be checked.
warns <- function(label, expr, pattern = NULL) {
  msg <- NULL
  value <- withCallingHandlers(expr, warning = function(w) {
    msg <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
  if (is.null(msg)) {
    ok(label, FALSE)
  } else if (!is.null(pattern) && !grepl(pattern, msg, fixed = FALSE)) {
    cat("  FAIL  ", label, "\n        warning did not match '", pattern,
        "':\n        ", msg, "\n", sep = "")
    .fail <<- .fail + 1L
  } else {
    ok(label, TRUE)
  }
  invisible(value)
}

# A rectangle in the form gate_geometry() returns.
rect <- function(min, max) {
  list(type = "rectangle", dims = names(min), min = min, max = max)
}

# min <= x < max on every dimension, as flowCore applies a rectangle.
in_rect <- function(x, g) {
  Reduce(`&`, lapply(g$dims, function(d) {
    x[, d] >= g$min[[d]] & x[, d] < g$max[[d]]
  }))
}

# --- an MTB-like sample ------------------------------------------------------
# Live cells split by a CD4 x CD8 quadrant (thresholds t4, t8), with IFNy and
# TNFa gates under CD4 and under CD8 drawn as cytokine x FSC-A rectangles. The
# two lineages have DIFFERENT cytokine cuts, as in the real workspaces.
#
# Each row is a group of identical cells and the calls the gater made for it:
#   b  sits at IFNy 80: above the CD4 cut (70), so IFNy+;
#   f  sits at IFNy 80 too, but under the CD8 cut (85), so IFNy-;
#   d  is above both CD4 cytokine cuts in IFNy but outside the gate on FSC-A,
#      so the gater called it IFNy-;
#   e  is at TNFa 87, just under the CD4 TNFa cut (88);
#   dp, dn have no cytokine gate at all, so their cytokines are undecided.
POP <- "/Lymphocytes/Single Cells/Live cells"
CH  <- c(CD4 = "Comp-APC-Cy7-A", CD8 = "Comp-PE-Cy7-A", IFNy = "Comp-PE-A",
         TNFa = "Comp-Pacific Blue-A")
# Channel -> marker, as the channel dictionary would give it; FSC-A is none.
CM  <- c(stats::setNames(names(CH), CH), "FSC-A" = NA_character_)

MTB_CELLS <- data.frame(
  group = c("a", "b", "c", "d", "e", "f", "g", "h", "dp", "dn"),
  n     = c( 40,  10,  10,   5,   5,  30,  10,   5,   10,   15),
  cd4   = c(150, 150, 150, 150, 150,  60,  60,  60,  150,   60),
  cd8   = c( 60,  60,  60,  60,  60, 160, 160, 160,  160,   60),
  ifn   = c( 60,  80, 100, 100,  60,  80,  90,  60,  100,  100),
  tnf   = c( 60,  60, 100,  60,  87,  60,  82,  85,  100,  100),
  fsc   = c(3e4, 3e4, 3e4, 6e4, 3e4, 3e4, 3e4, 3e4,  3e4,  3e4),
  CD4   = c(  1,   1,   1,   1,   1,   0,   0,   0,    1,    0),
  CD8   = c(  0,   0,   0,   0,   0,   1,   1,   1,    1,    0),
  IFNy  = c(  0,   1,   1,   0,   0,   0,   1,   0,   NA,   NA),
  TNFa  = c(  0,   0,   1,   0,   0,   0,   1,   1,   NA,   NA),
  label = c("CD4+ | none", "CD4+ | IFNy+", "CD4+ | IFNy+ TNFa+",
            "CD4+ | none", "CD4+ | none", "CD8+ | none",
            "CD8+ | IFNy+ TNFa+", "CD8+ | TNFa+", "CD4+ CD8+", "CD4- CD8-"),
  stringsAsFactors = FALSE)

mtb_sample <- function(t4 = 90, t8 = 120,
                       ifn_cut = c(CD4 = 70, CD8 = 85),
                       tnf_cut = c(CD4 = 88, CD8 = 80),
                       ifn_gate = c(CD4 = "IFNy", CD8 = "INFy")) {
  cells <- MTB_CELLS[rep(seq_len(nrow(MTB_CELLS)), MTB_CELLS$n), ]
  x <- cbind(cells$cd4, cells$cd8, cells$ifn, cells$tnf, cells$fsc)
  colnames(x) <- c(CH[["CD4"]], CH[["CD8"]], CH[["IFNy"]], CH[["TNFa"]],
                   "FSC-A")
  quad <- function(lo4, hi4, lo8, hi8) {
    rect(stats::setNames(c(lo4, lo8), CH[c("CD4", "CD8")]),
         stats::setNames(c(hi4, hi8), CH[c("CD4", "CD8")]))
  }
  cyt <- function(ch, cut) {
    rect(stats::setNames(c(cut, 1e4), c(ch, "FSC-A")),
         stats::setNames(c(250, 5e4), c(ch, "FSC-A")))
  }
  paths <- c(POP, paste0(POP, "/", c(" ", "  ", "CD4", "CD8")),
             paste0(POP, "/CD4/", c(ifn_gate[["CD4"]], "TNF")),
             paste0(POP, "/CD8/", c(ifn_gate[["CD8"]], "TNF")))
  gates <- gate_descendants(c("root", "/Lymphocytes",
                              "/Lymphocytes/Single Cells", paths), POP)
  gates$geom <- list(
    quad(t4, Inf, t8, Inf),      # " "  CD4+CD8+
    quad(-Inf, t4, -Inf, t8),    # "  " CD4-CD8-
    quad(t4, Inf, -Inf, t8),     # CD4
    quad(-Inf, t4, t8, Inf),     # CD8
    cyt(CH[["IFNy"]], ifn_cut[["CD4"]]),
    cyt(CH[["TNFa"]], tnf_cut[["CD4"]]),
    cyt(CH[["IFNy"]], ifn_cut[["CD8"]]),
    cyt(CH[["TNFa"]], tnf_cut[["CD8"]]))
  mem <- matrix(FALSE, nrow(x), nrow(gates),
                dimnames = list(NULL, gates$path))
  for (i in seq_len(nrow(gates))) {
    par <- if (gates$parent[i] == POP) TRUE else mem[, gates$parent[i]]
    mem[, i] <- in_rect(x, gates$geom[[i]]) & par
  }
  list(x = x, mem = mem, gates = gates, group = cells$group,
       expect = cells)
}

# --- gate_descendants --------------------------------------------------------
cat("\ngate_descendants()\n")

tree <- c("root", "/Lymphocytes", "/Lymphocytes/Single Cells", POP,
          paste0(POP, "/ "), paste0(POP, "/  "),
          paste0(POP, "/CD4"), paste0(POP, "/CD8"),
          paste0(POP, "/CD4/IFNy"), paste0(POP, "/CD4/TNF"),
          paste0(POP, "/CD8/INFy"),
          "/Lymphocytes/Single Cells/Live cells 2",
          "/Lymphocytes/Single Cells/Live cells 2/CD4")
gd <- gate_descendants(tree, POP)
eq("only gates strictly below the population, in input order",
   gd$path, tree[5:11])
eq("depth: quadrant gates 1, cytokine gates 2",
   gd$depth, c(1L, 1L, 1L, 1L, 2L, 2L, 2L))
eq("a quadrant gate's parent is the population", gd$parent[1:4], rep(POP, 4))
eq("a cytokine gate's parent is its lineage gate", gd$parent[5:7],
   paste0(POP, c("/CD4", "/CD4", "/CD8")))
ok("gate names of one and of two spaces stay two separate gates",
   identical(gd$path[1:2], paste0(POP, c("/ ", "/  "))) &&
     identical(gd$parent[1:2], c(POP, POP)))
ok("a sibling whose name only starts like the population is not below it",
   !any(grepl("Live cells 2", gd$path, fixed = TRUE)))

gd_slash <- gate_descendants(tree, paste0(POP, "/"))
eq("a trailing slash on the population gives the same paths",
   gd_slash$path, gd$path)
eq("...the same depths", gd_slash$depth, gd$depth)
eq("...and the same parents", gd_slash$parent, gd$parent)

# A gate whose NAME contains "/" ("CD4/CD8 DP" directly under the population)
# looks like a child of CD4 to the string rule; the real tree says otherwise.
slashy <- c(POP, paste0(POP, "/CD4"), paste0(POP, "/CD4/CD8 DP"),
            paste0(POP, "/CD4/CD8 DP/IFNy"))
gd_str <- gate_descendants(slashy, POP)
eq("without the tree, the string rule nests a slashed name",
   gd_str$depth, c(1L, 2L, 3L))
gd_par <- gate_descendants(slashy, POP,
                           parents = c(stats::setNames(POP, slashy[3])))
eq("parents override: the slashed gate sits directly under the population",
   gd_par$parent[2], POP)
eq("...depths follow the overridden tree, also for its child",
   gd_par$depth, c(1L, 1L, 2L))
eq("...paths missing from `parents` keep the string rule",
   gd_par$parent[c(1, 3)], c(POP, paste0(POP, "/CD4/CD8 DP")))

# "root" is a valid population (gs_get_pop_paths lists it and
# resolve_population() accepts it), and every other path lies below it even
# though none starts with "root/".
gd_root <- gate_descendants(tree, "root")
eq("population 'root': every gate is below it",
   gd_root$path, tree[-1])
eq("...with depth counted from the root",
   gd_root$depth[gd_root$path %in% c("/Lymphocytes", POP, paste0(POP, "/CD4/IFNy"))],
   c(1L, 3L, 5L))

gd_none <- gate_descendants(tree, "/Nope")
ok("no descendants: an empty frame with path, parent, depth",
   nrow(gd_none) == 0L && identical(names(gd_none), c("path", "parent", "depth")))

# --- gate_decisions_one: rectangles ------------------------------------------
cat("\ngate_decisions_one(): rectangles\n")

# 1000 parent cells at A = 1..1000; FSC-A spread so that its bounds DO exclude.
xr <- cbind(A = as.numeric(1:1000), B = as.numeric(1:1000),
            "FSC-A" = rep(c(1e4, 5e4), 500))
cm_r <- c(A = "M", "FSC-A" = NA_character_)
all_in <- rep(TRUE, 1000)
dec1 <- function(geom, in_parent = all_in, ...) {
  gate_decisions_one(xr, rep(FALSE, 1000), in_parent, geom, cm_r, ...)
}

d <- dec1(rect(c(A = 100.5), c(A = Inf)))
eq("a binding lower bound is a positive decision at that bound",
   d[, c("channel", "marker", "polarity", "cut", "method")],
   data.frame(channel = "A", marker = "M", polarity = 1L, cut = 100.5,
              method = "rectangle", stringsAsFactors = FALSE))
eq("...and records the share of the parent it excludes", d$excluded, 0.1)

d <- dec1(rect(c(A = -Inf), c(A = 900.5)))
eq("a binding upper bound is a negative decision at that bound",
   d[, c("polarity", "cut", "excluded")],
   data.frame(polarity = -1L, cut = 900.5, excluded = 0.1))

d_band <- dec1(rect(c(A = 100.5), c(A = 900.5)))
eq("a band (both bounds binding, 10 % each side) makes no decision",
   sum(d_band$type == "decision"), 0L)
eq("...but is reported as a band axis", d_band$type, "band")
d_clip <- dec1(rect(c(A = 100.5), c(A = 990.5)))
eq("a far edge clipping only 1 % leaves the other edge as the decision",
   d_clip[, c("polarity", "cut", "type")],
   data.frame(polarity = 1L, cut = 100.5, type = "decision",
              stringsAsFactors = FALSE))
eq("bounds that exclude nothing make no decision",
   nrow(dec1(rect(c(A = 0.5), c(A = 1000.5)))), 0L)

# FSC-A [2e4, Inf) excludes half the parent, but FSC-A is not a marker.
d <- dec1(rect(c(A = 100.5, "FSC-A" = 2e4), c(A = Inf, "FSC-A" = Inf)))
eq("an FSC-A dimension is ignored even when its bound is binding",
   d$channel, "A")
d <- dec1(rect(c("FSC-A" = 0, A = 100.5), c("FSC-A" = 1e6, A = Inf)))
eq("...and a non-binding FSC-A dimension is ignored too", d$channel, "A")

eq("5 of 1000 below the bound is exactly min_excluded: a decision",
   nrow(dec1(rect(c(A = 5.5), c(A = Inf)))), 1L)
eq("4 of 1000 below the bound is under min_excluded: no decision",
   nrow(dec1(rect(c(A = 4.5), c(A = Inf)))), 0L)
eq("min_excluded is honoured when given",
   nrow(dec1(rect(c(A = 100.5), c(A = Inf)), min_excluded = 0.2)), 0L)

eq("a channel with no marker is ignored",
   nrow(dec1(rect(c(B = 100.5), c(B = Inf)))), 0L)
eq("a channel whose marker is an empty string is ignored",
   nrow(gate_decisions_one(xr, all_in, all_in, rect(c(A = 100.5), c(A = Inf)),
                           c(A = ""))), 0L)
eq("a gate channel missing from the data is ignored",
   nrow(dec1(rect(c(Z = 100.5), c(Z = Inf)))), 0L)
eq("a Boolean gate makes no decision",
   nrow(dec1(list(type = "boolean", dims = character(0)))), 0L)
eq("...even when it lists dimensions",
   nrow(dec1(list(type = "boolean", dims = "A"))), 0L)
eq("no geometry at all makes no decision", nrow(dec1(NULL)), 0L)
eq("an empty parent makes no decision",
   nrow(dec1(rect(c(A = 100.5), c(A = Inf)), in_parent = rep(FALSE, 1000))), 0L)

# Exclusion is measured among the PARENT's cells only.
d <- dec1(rect(c(A = 600.5), c(A = Inf)), in_parent = xr[, "A"] > 500)
eq("exclusion is the share of the parent's cells, not of all cells",
   d$excluded, 0.2)

# --- gate_decisions_one: polygons (membership) -------------------------------
cat("\ngate_decisions_one(): polygons and other shapes\n")

xp <- cbind(A = as.numeric(1:20))
poly <- list(type = "polygon", dims = c("A", "FSC-A"))
decp <- function(inside, x = xp, in_parent = rep(TRUE, nrow(x)), ...) {
  gate_decisions_one(x, inside, in_parent, poly, c(A = "M"), ...)
}

d <- decp(xp[, "A"] >= 11)
eq("members above non-members (AUC 1): a positive membership decision",
   d[, c("channel", "marker", "polarity", "method")],
   data.frame(channel = "A", marker = "M", polarity = 1L,
              method = "membership", stringsAsFactors = FALSE))
eq("...with the Youden cut midway across the gap", d$cut, 10.5)
eq("...and the non-member share as `excluded`", d$excluded, 0.5)

d <- decp(xp[, "A"] <= 10)
eq("members below non-members (AUC 0): a negative decision, same cut",
   d[, c("polarity", "cut")], data.frame(polarity = -1L, cut = 10.5))

eq("interleaved members (AUC 0.55) make no decision",
   nrow(decp(xp[, "A"] %% 2 == 0)), 0L)
eq("9 members (fewer than min_class) make no decision",
   nrow(decp(xp[, "A"] >= 12)), 0L)
eq("9 non-members make no decision either", nrow(decp(xp[, "A"] >= 10)), 0L)

# Youden: the cut is the midpoint of the gap, wherever the gap is.
xg <- cbind(A = c(1:10, 30:39))
d <- decp(xg[, "A"] >= 30, x = xg)
eq("Youden cut on separable data is the midpoint of the gap", d$cut, 20)

# Cells outside the parent do not move the cut.
xo <- cbind(A = c(1:10, 30:39, rep(15, 50)))
d <- decp(xo[, "A"] >= 30, x = xo, in_parent = c(rep(TRUE, 20), rep(FALSE, 50)))
eq("cells outside the parent do not enter the AUC or the cut", d$cut, 20)

# The AUC thresholds are symmetric: "at least 0.8 or at most 0.2". Ten
# non-members at 1..10; eight members beyond them on one side, two on the
# other, so the AUC is exactly 0.8 or exactly 0.2.
non <- as.numeric(1:10)
x08 <- cbind(A = c(10 + (1:8) / 10, 0.1, 0.2, non))
x02 <- cbind(A = c((1:8) / 10, 10.1, 10.2, non))
mem10 <- c(rep(TRUE, 10), rep(FALSE, 10))
d08 <- decp(mem10, x = x08)
eq("AUC exactly 0.8 is a positive decision at the Youden cut",
   d08[, c("polarity", "cut")], data.frame(polarity = 1L, cut = 10.05))
d02 <- decp(mem10, x = x02)
eq("AUC exactly 0.2 is a negative decision at the Youden cut",
   d02[, c("polarity", "cut")], data.frame(polarity = -1L, cut = 0.9))

# --- gate_decisions_sample ---------------------------------------------------
cat("\ngate_decisions_sample()\n")

s1 <- mtb_sample()
d1 <- gate_decisions_sample(s1$x, s1$mem, s1$gates, CM, POP)
P <- function(...) paste0(POP, ...)
exp_dec <- data.frame(
  gate = c(P("/ "), P("/ "), P("/  "), P("/  "), P("/CD4"), P("/CD4"),
           P("/CD8"), P("/CD8"), P("/CD4/IFNy"), P("/CD4/TNF"),
           P("/CD8/INFy"), P("/CD8/TNF")),
  marker = c("CD4", "CD8", "CD4", "CD8", "CD4", "CD8", "CD4", "CD8",
             "IFNy", "TNFa", "IFNy", "TNFa"),
  polarity = c(1L, 1L, -1L, -1L, 1L, -1L, -1L, 1L, 1L, 1L, 1L, 1L),
  cut = c(90, 120, 90, 120, 90, 120, 90, 120, 70, 88, 85, 80),
  stringsAsFactors = FALSE)
eq("columns in the documented order", names(d1),
   c("gate", "parent", "depth", "channel", "marker", "polarity", "cut",
     "method", "excluded", "n_parent", "n_member", "n_axes"))
eq("quadrant and cytokine gates give the expected (gate, marker, polarity, cut)",
   d1[, c("gate", "marker", "polarity", "cut")], exp_dec)
ok("FSC-A never becomes a decision", !any(d1$channel == "FSC-A"))
eq("every decision here comes from a rectangle", unique(d1$method), "rectangle")
eq("depth is carried from the tree", d1$depth, rep(c(1L, 2L), c(8, 4)))
eq("excluded shares over each gate's parent",
   d1$excluded,
   c(c(60, 85, 80, 55, 60, 55, 80, 85) / 140,
     45 / 70, 60 / 70, 35 / 45, 30 / 45))
eq("n_parent: all population cells, or the lineage gate's cells",
   d1$n_parent, c(rep(140L, 8), 70L, 70L, 45L, 45L))
eq("n_member counts gate members within the parent (FSC-A-excluded cells out)",
   d1$n_member, c(10L, 10L, 15L, 15L, 70L, 70L, 45L, 45L, 20L, 10L, 10L, 15L))

g_missing <- s1$gates
g_missing$path[3] <- P("/not in mem")
d_missing <- gate_decisions_sample(s1$x, s1$mem, g_missing, CM, POP)
ok("a gate with no membership column is skipped",
   !any(d_missing$gate == P("/not in mem")) && nrow(d_missing) == 10L)

# --- gate_marker_calls -------------------------------------------------------
cat("\ngate_marker_calls()\n")

c1 <- gate_marker_calls(s1$x, s1$mem, d1, POP)
ex <- s1$expect
eq("one column per decided marker", colnames(c1$calls),
   c("CD4", "CD8", "IFNy", "TNFa"))
eq("CD4 calls follow the quadrant (CD4 and DP positive)",
   c1$calls[, "CD4"], as.integer(ex$CD4))
eq("CD8 calls follow the quadrant (CD8 and DP positive)",
   c1$calls[, "CD8"], as.integer(ex$CD8))
eq("IFNy calls: each lineage against its own cut, FSC-A respected",
   c1$calls[, "IFNy"], as.integer(ex$IFNy))
eq("TNFa calls: each lineage against its own cut",
   c1$calls[, "TNFa"], as.integer(ex$TNFa))
ok("IFNy 80 is positive under CD4 (cut 70) and negative under CD8 (cut 85)",
   all(c1$calls[s1$group == "b", "IFNy"] == 1L) &&
     all(c1$calls[s1$group == "f", "IFNy"] == 0L))
ok("a CD4 cell above the IFNy cut but outside the gate on FSC-A is negative",
   all(c1$calls[s1$group == "d", "IFNy"] == 0L) &&
     all(c1$values[s1$group == "d", "IFNy"] > c1$cuts[s1$group == "d", "IFNy"]))
ok("double negatives have no cytokine decision, even far above every cut",
   all(is.na(c1$calls[s1$group == "dn", c("IFNy", "TNFa")])))
ok("double positives have no cytokine decision either",
   all(is.na(c1$calls[s1$group == "dp", c("IFNy", "TNFa")])))
lineage <- ifelse(s1$group %in% c("a", "b", "c", "d", "e"), "CD4",
           ifelse(s1$group %in% c("f", "g", "h"), "CD8", "none"))
eq("the IFNy cut per cell is its lineage's cut (undecided: sample median)",
   c1$cuts[, "IFNy"],
   unname(c(CD4 = 70, CD8 = 85, none = 77.5)[lineage]))
eq("the TNFa cut per cell is its lineage's cut (undecided: sample median)",
   c1$cuts[, "TNFa"],
   unname(c(CD4 = 88, CD8 = 80, none = 84)[lineage]))
eq("the quadrant markers share one cut per sample",
   c(unique(c1$cuts[, "CD4"]), unique(c1$cuts[, "CD8"])), c(90, 120))
eq("values are the raw channel values of each marker",
   unname(c1$values[, "IFNy"]), unname(s1$x[, CH[["IFNy"]]]))

# Mixed polarity: one gate keeps X >= 100, another keeps X < 50, both directly
# under the population. A cell inside neither falls back to x >= cut, with the
# cut the median of the two (75).
xm <- cbind(X = c(40, 120, 60, 80, 75))
mem_m <- cbind(G1 = xm[, "X"] >= 100, G2 = xm[, "X"] < 50)
dec_m <- data.frame(gate = c("G1", "G2"), parent = c(POP, POP), depth = 1L,
                    channel = "X", marker = "M", polarity = c(1L, -1L),
                    cut = c(100, 50), stringsAsFactors = FALSE)
cm_m <- gate_marker_calls(xm, mem_m, dec_m, POP)
eq("inside the positive gate: 1; inside the negative gate: 0",
   cm_m$calls[1:2, "M"], c(0L, 1L))
eq("outside both (mixed polarity): x >= median cut, the edge included",
   cm_m$calls[3:5, "M"], c(0L, 1L, 1L))
eq("...with the median cut recorded", unique(cm_m$cuts[, "M"]), 75)
# Only the negative gate applies: a cell outside it is positive.
cm_n <- gate_marker_calls(xm, mem_m[, "G2", drop = FALSE], dec_m[2, ], POP)
eq("outside every applicable negative-polarity gate means positive",
   cm_n$calls[, "M"], c(0L, 1L, 1L, 1L, 1L))

# One marker gated twice: under the population (cut 100) and again under CD4
# (cut 60). A CD4 cell is decided by the deeper, CD4-level gate; other cells
# by the population-level one.
xd <- cbind(X = c(80, 80, 120, 50))
cd4 <- c(TRUE, FALSE, TRUE, TRUE)
mem_d <- cbind(CD4 = cd4, Mhi = xd[, "X"] >= 100, Mlo = cd4 & xd[, "X"] >= 60)
dec_d <- data.frame(gate = c("Mhi", "Mlo"), parent = c(POP, "CD4"),
                    depth = c(1L, 2L), channel = "X", marker = "M",
                    polarity = 1L, cut = c(100, 60), stringsAsFactors = FALSE)
cm_d <- gate_marker_calls(xd, mem_d, dec_d, POP)
eq("a marker gated at two depths: the deeper gate's cut for its cells",
   cm_d$cuts[, "M"], c(60, 100, 60, 60))
eq("...and its call (X 80 is positive in CD4, negative elsewhere)",
   cm_d$calls[, "M"], c(1L, 0L, 1L, 0L))

# --- labels are named after markers, not gates -------------------------------
cat("\ncombination labels from the MTB-like samples\n")

s2 <- mtb_sample(t4 = 100, t8 = 130, ifn_cut = c(CD4 = 75, CD8 = 88),
                 ifn_gate = c(CD4 = "INFy", CD8 = "IFNy"))
d2 <- gate_decisions_sample(s2$x, s2$mem, s2$gates, CM, POP)
c2 <- gate_marker_calls(s2$x, s2$mem, d2, POP)
ok("the second sample spells its CD4 IFNg gate INFy",
   P("/CD4/INFy") %in% d2$gate)
tiers <- marker_tiers(rbind(d1, d2), panel_order = c("CD4", "CD8", "IFNy", "TNFa"))
eq("tiers: CD4, CD8 first; IFNy, TNFa second",
   tiers, c(CD4 = 1L, CD8 = 1L, IFNy = 2L, TNFa = 2L))
lab1 <- combination_labels(c1$calls, tiers)
lab2 <- combination_labels(c2$calls, tiers)
eq("labels per cell match the gater's combinations", lab1, ex$label)
eq("IFNy- and INFy-named gates give the same labels", lab2, lab1)

# --- calls_key ---------------------------------------------------------------
cat("\ncalls_key()\n")

ck <- matrix(c(1L, 0L, NA, 0L, 0L, 1L, 1L, 0L, NA), 3, 3, byrow = TRUE,
             dimnames = list(NULL, c("CD4", "CD8", "IFNy")))
eq("one string per cell, NA spelled x", calls_key(ck), c("10x", "001", "10x"))
eq("no columns: an empty key per cell",
   calls_key(matrix(NA_integer_, 2, 0)), c("", ""))

# --- stratified_sample -------------------------------------------------------
cat("\nstratified_sample()\n")

set.seed(11)
key <- sample(rep(c("A", "B", "C"), c(50, 7, 13)))
set.seed(1); st <- stratified_sample(key, cap = 10)
kept <- table(key[st$idx])
eq("up to `cap` cells per key, rare keys kept whole",
   as.vector(kept[c("A", "B", "C")]), c(10L, 7L, 10L))
eq("weights sum to the original count of every key",
   as.vector(tapply(st$weight, key[st$idx], sum)[c("A", "B", "C")]),
   c(50, 7, 13))
eq("weight = cells of the key / cells kept",
   as.vector(tapply(st$weight, key[st$idx], unique)[c("A", "B", "C")]),
   c(5, 1, 1.3))
ok("indices are sorted and unique",
   !is.unsorted(st$idx, strictly = TRUE))
set.seed(1); st_again <- stratified_sample(key, cap = 10)
ok("the same seed gives the same sample", identical(st, st_again))
eq("cap 0 keeps every cell with weight 1",
   stratified_sample(key, 0), list(idx = seq_along(key), weight = rep(1, 70)))
eq("a negative cap does the same",
   stratified_sample(key, -5), list(idx = seq_along(key), weight = rep(1, 70)))
eq("no cells: nothing kept",
   stratified_sample(character(0), 10), list(idx = integer(0), weight = numeric(0)))

# --- balanced_fit_index ------------------------------------------------------
cat("\nbalanced_fit_index()\n")

labs <- rep(c("a", "b", "c"), c(5, 100, 100))
set.seed(5); bf <- balanced_fit_index(labs, 60)
eq("each label capped at budget / labels (20), the small one whole",
   as.vector(table(labs[bf])[c("a", "b", "c")]), c(5L, 20L, 20L))
ok("...unused share is not handed to the large labels", length(bf) <= 60)
ok("...sorted and unique", !is.unsorted(bf, strictly = TRUE))
set.seed(5); bf_again <- balanced_fit_index(labs, 60)
ok("the same seed gives the same fit set", identical(bf, bf_again))
eq("budget 0 means every cell", balanced_fit_index(labs, 0), seq_along(labs))
# A run smaller than the budget is balanced too: that is what keeps the
# abundant combinations from taking all the room in a small cohort.
labs2 <- rep(c("rare", "big1", "big2"), c(5, 400, 400))
set.seed(6); bf2 <- balanced_fit_index(labs2, 1000)
eq("a budget above the cell count still caps every label (1000 / 3 = 333)",
   as.vector(table(labs2[bf2])[c("rare", "big1", "big2")]), c(5L, 333L, 333L))
eq("a budget below the number of labels keeps one cell per label",
   length(balanced_fit_index(c("x", "y", "z"), 2)), 3L)

# --- marker_tiers ------------------------------------------------------------
cat("\nmarker_tiers()\n")

dt <- data.frame(marker = c("CD8", "CD4", "TNFa", "IFNy", "IL2", "CD4"),
                 depth = c(1L, 1L, 2L, 2L, 2L, 1L), stringsAsFactors = FALSE)
eq("without a panel order: tier, then first appearance",
   marker_tiers(dt), c(CD8 = 1L, CD4 = 1L, TNFa = 2L, IFNy = 2L, IL2 = 2L))
eq("panel order sorts within a tier, never across tiers",
   marker_tiers(dt, panel_order = c("IL2", "IFNy", "TNFa", "CD4", "CD8")),
   c(CD4 = 1L, CD8 = 1L, IL2 = 2L, IFNy = 2L, TNFa = 2L))
eq("markers missing from a partial panel order follow in first appearance",
   marker_tiers(dt, panel_order = c("IFNy", "CD4")),
   c(CD4 = 1L, CD8 = 1L, IFNy = 2L, TNFa = 2L, IL2 = 2L))
eq("a marker's tier is its median depth",
   marker_tiers(data.frame(marker = c("X", "X", "X", "Y"),
                           depth = c(1L, 3L, 3L, 1L))),
   c(Y = 1L, X = 3L))
eq("no decisions: no tiers", marker_tiers(d1[0, ]), integer(0))

# --- combination_labels ------------------------------------------------------
cat("\ncombination_labels()\n")

tl <- c(CD4 = 1L, CD8 = 1L, IFNy = 2L, TNFa = 2L)
cl <- matrix(c(1L, 0L, 1L, 1L,
               1L, 0L, 0L, 0L,
               0L, 0L, NA, NA,
               1L, 1L, NA, NA,
               0L, 1L, 0L, 1L,
               NA, NA, NA, NA,
               0L, NA, NA, NA),
             ncol = 4, byrow = TRUE, dimnames = list(NULL, names(tl)))
eq("labels: positives per tier, 'none', tier 1 spelled negative, NA tiers omitted",
   combination_labels(cl, tl),
   c("CD4+ | IFNy+ TNFa+", "CD4+ | none", "CD4- CD8-", "CD4+ CD8+",
     "CD8+ | TNFa+", "(no decision)", "CD4-"))
shuffled <- cl[, c("TNFa", "CD8", "IFNy", "CD4")]
eq("label order follows the tiers, not the column order of the calls",
   combination_labels(shuffled, tl), combination_labels(cl, tl))
eq("columns not in the tiers are ignored",
   combination_labels(cbind(cl, ALD = 1L), tl), combination_labels(cl, tl))
eq("a single decision marker still labels",
   combination_labels(matrix(c(1L, 0L), 2, 1, dimnames = list(NULL, "CD4")),
                      c(CD4 = 1L)),
   c("CD4+", "CD4-"))

# --- pool_rare_labels --------------------------------------------------------
cat("\npool_rare_labels()\n")

pl <- rep(c("a", "b", "c"), c(5, 2, 3))
eq("labels under min_cells become 'other (rare)'",
   pool_rare_labels(pl, min_cells = 3),
   rep(c("a", "other (rare)", "c"), c(5, 2, 3)))
eq("weights count, not cells", pool_rare_labels(pl, rep(c(1, 10, 1), c(5, 2, 3)),
                                                min_cells = 4),
   rep(c("a", "b", "other (rare)"), c(5, 2, 3)))
eq("min_cells 0 pools nothing", pool_rare_labels(pl, min_cells = 0), pl)
eq("the constant is the documented label", .GL_RARE, "other (rare)")

# --- fill_missing_decisions --------------------------------------------------
cat("\nfill_missing_decisions()\n")

tf <- c(CD4 = 1L, CD8 = 1L, IFNy = 2L)
cf <- cbind(CD4  = c(1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L),
            CD8  = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L),
            IFNy = c(1L, 0L, 0L, 0L, 1L, NA, NA, NA, 1L, NA))
sf <- c("s1", "s1", "s1", "s2", "s2", "s3", "s1", "s2", "s1", "s3")
ff <- fill_missing_decisions(cf, tf, sf)
eq("a sample missing a decision most of its group has: negative there",
   unname(ff$calls[6, "IFNy"]), 0L)
eq("a group where nobody has the decision is left undecided",
   ff$calls[7:8, "IFNy"], c(NA_integer_, NA_integer_))
eq("a group with exactly half decided is not filled (needs a majority)",
   unname(ff$calls[10, "IFNy"]), NA_integer_)
eq("decided calls are untouched", ff$calls[-c(6, 10), ], cf[-c(6, 10), ])
eq("the fill is reported per sample, group and marker",
   ff$filled, data.frame(sample = "s3", group = "10", marker = "IFNy",
                         cells = 1L, stringsAsFactors = FALSE))
nf <- fill_missing_decisions(cf[, 1:2], c(CD4 = 1L, CD8 = 1L), sf)
ok("a single tier: nothing to fill",
   identical(nf$calls, cf[, 1:2]) && is.null(nf$filled))

# --- gate_scales -------------------------------------------------------------
cat("\ngate_scales()\n")

# Negatives at cut + (-5..5); positives and undecided cells far away.
vals <- cbind(M = c(1:11, rep(200, 5), 1000))
cuts <- cbind(M = rep(6, 17))
calls <- cbind(M = c(rep(0L, 11), rep(1L, 5), NA))
eq("MAD of (value - cut) over the negatives only, times 1.4826",
   gate_scales(vals, cuts, calls), c(M = 1.4826 * 3))

flat <- cbind(M = rep(50, 12))
eq("a zero spread is floored at 1 % of the channel range",
   gate_scales(flat, flat, cbind(M = rep(0L, 12))), c(M = 2.56))
eq("...of the range given",
   gate_scales(flat, flat, cbind(M = rep(0L, 12)), channel_range = 1000),
   c(M = 10))

# Weighted: d = 0..9, the upper five weighted 10x. Unweighted median 4 and
# MAD 2; weighted median 7 and MAD 1.
v10 <- cbind(M = as.numeric(0:9)); z10 <- cbind(M = rep(0, 10))
n10 <- cbind(M = rep(0L, 10))
eq("unweighted MAD", gate_scales(v10, z10, n10, floor_frac = 0.001),
   c(M = 1.4826 * 2))
eq("weighted MAD", gate_scales(v10, z10, n10, weights = rep(c(1, 10), each = 5),
                               floor_frac = 0.001),
   c(M = 1.4826 * 1))
set.seed(7)
vr <- cbind(M = round(rnorm(40, 50, 8), 1)); wr <- sample(1:4, 40, replace = TRUE)
cr <- cbind(M = rep(45, 40)); kr <- cbind(M = rep(0L, 40))
eq("integer weights equal repeating the cells",
   gate_scales(vr, cr, kr, weights = wr, floor_frac = 0),
   gate_scales(vr[rep(1:40, wr), , drop = FALSE], cr[rep(1:40, wr), , drop = FALSE],
               kr[rep(1:40, wr), , drop = FALSE], floor_frac = 0))

# Fewer than 10 negatives: every finite cell is used instead. d over all cells
# is (0 x5, 1..20): median 8, MAD 6. The negatives alone would give 0.
vs <- cbind(M = c(rep(0, 5), 1:20)); cs <- cbind(M = rep(0, 25))
ks <- cbind(M = c(rep(0L, 5), rep(1L, 20)))
eq("fewer than 10 negatives: the scale comes from every cell",
   gate_scales(vs, cs, ks), c(M = 1.4826 * 6))

# --- gate_features -----------------------------------------------------------
cat("\ngate_features()\n")

set.seed(9)
nf_ <- 400
fv <- cbind(A = c(rnorm(nf_ - 4, 60, 30), 1e6, -1e6, 1e6, -1e6),
            B = rnorm(nf_, 100, 50))
fc <- cbind(A = rep(70, nf_), B = runif(nf_, 80, 120))
fk <- cbind(A = c(rbinom(nf_ - 4, 1, 0.4), 0L, 1L, NA, 1L),
            B = ifelse(runif(nf_) < 0.2, NA, rbinom(nf_, 1, 0.5)))
storage.mode(fk) <- "integer"
sig <- c(A = 5, B = 3)
separated <- function(f, k, m) {
  pos <- !is.na(k[, m]) & k[, m] == 1L
  min(f[pos, m]) > max(f[!pos, m])
}
for (g in c(1.01, 1.5, 2, 3)) {
  fg <- gate_features(fv, fc, fk, sig, gap = g)
  ok(sprintf("gap %.2f: positives and negatives never overlap, on every marker", g),
     separated(fg, fk, "A") && separated(fg, fk, "B"))
}
f2 <- gate_features(fv, fc, fk, sig, gap = 2, softness = 2)
eq("the feature is tanh((x - cut) / (softness * sigma)) +/- gap",
   unname(f2[1:5, "A"]),
   tanh((fv[1:5, "A"] - 70) / (2 * 5)) + 2 * ifelse(fk[1:5, "A"] == 1L, 1, -1))
fk0 <- fk; fk0[is.na(fk0)] <- 0L
ok("NA calls sit on the negative side (exactly as a negative call would)",
   identical(f2, gate_features(fv, fc, fk0, sig, gap = 2, softness = 2)))
eq("a cell at the cut lands exactly on +gap or -gap",
   unname(gate_features(cbind(A = c(70, 70, 70)), cbind(A = rep(70, 3)),
                        cbind(A = c(1L, 0L, NA)), c(A = 5), gap = 2)[, "A"]),
   c(2, -2, -2))

ext <- cbind(X = rnorm(nf_, 10, 3), Y = rexp(nf_), K = rep(4, nf_))
fe <- gate_features(fv, fc, fk, sig, extra = ext)
eq("extra markers are appended after the decision markers",
   colnames(fe), c("A", "B", "X", "Y", "K"))
ok("extra markers stay inside (-1, 1)",
   all(abs(fe[, c("X", "Y")]) < 1))
eq("a constant extra marker becomes 0", unname(fe[, "K"]), rep(0, nf_))
ok("extra markers take no jump (both signs present)",
   any(fe[, "X"] > 0) && any(fe[, "X"] < 0))

# --- parse_sample_filter -----------------------------------------------------
cat("\nparse_sample_filter()\n")

eq("one '=' clause",
   parse_sample_filter("stimulation=stimulated"),
   list(list(column = "stimulation", negate = FALSE, values = "stimulated")))
eq("a '!=' clause with several values",
   parse_sample_filter("visit!=W0,W2"),
   list(list(column = "visit", negate = TRUE, values = c("W0", "W2"))))
eq("several clauses, whitespace everywhere, an empty trailing clause",
   parse_sample_filter("  time point = W0 , W12 ;patient!= P1 ; "),
   list(list(column = "time point", negate = FALSE, values = c("W0", "W12")),
        list(column = "patient", negate = TRUE, values = "P1")))
eq("no filter: no clauses", parse_sample_filter(NULL), list())
eq("a blank filter: no clauses", parse_sample_filter("   "), list())
throws("a clause without an operator is an error",
       parse_sample_filter("visit W0"), "not of the form")
throws("a clause without a column is an error",
       parse_sample_filter("=W0"), "not of the form")

# --- apply_sample_filter -----------------------------------------------------
cat("\napply_sample_filter()\n")

ann <- data.frame(sample = paste0("s", 1:5),
                  visit = c("W0", "W2", "W12", "Unassigned", "W2"),
                  stim = c("neg", "stim", "stim", "stim", "neg"),
                  stringsAsFactors = FALSE)
eq("'=' keeps the listed levels",
   apply_sample_filter(ann, parse_sample_filter("visit=W0,W2")),
   c(TRUE, TRUE, FALSE, FALSE, TRUE))
eq("an Unassigned sample fails an '=' clause",
   apply_sample_filter(ann, parse_sample_filter("visit=W0,W2,W12"))[4], FALSE)
eq("'!=' drops the listed levels",
   apply_sample_filter(ann, parse_sample_filter("visit!=W2")),
   c(TRUE, FALSE, TRUE, TRUE, FALSE))
eq("clauses combine with AND",
   apply_sample_filter(ann, parse_sample_filter("stim=stim; visit!=W12")),
   c(FALSE, TRUE, FALSE, TRUE, FALSE))
eq("no clauses keep everything", apply_sample_filter(ann, list()), rep(TRUE, 5))
throws("an unknown column is an error that lists the columns",
       apply_sample_filter(ann, parse_sample_filter("tube=A")),
       "tube.*Columns: sample, visit, stim")

# --- natural_order -----------------------------------------------------------
cat("\nnatural_order()\n")

nat <- function(x) x[natural_order(x)]
eq("digit runs compare as numbers", nat(c("W110", "W12", "W5")),
   c("W5", "W12", "W110"))
eq("mixed text and numbers, several runs",
   nat(c("V2W67", "V10W0", "V2W8", "V1W6")), c("V1W6", "V2W8", "V2W67", "V10W0"))
eq("text compares case-insensitively", nat(c("gamma", "Beta", "alpha")),
   c("alpha", "Beta", "gamma"))
eq("a prefix sorts before its extensions", nat(c("W5a", "W5", "W")),
   c("W", "W5", "W5a"))
eq("plain numbers", nat(c("10", "9", "100", "1")), c("1", "9", "10", "100"))
eq("nothing to order", natural_order(character(0)), integer(0))
eq("empty strings are still ordered (a permutation comes back)",
   sort(natural_order(c("", ""))), 1:2)

# --- parse_level_order and ordered_levels ------------------------------------
cat("\nparse_level_order() and ordered_levels()\n")

eq("columns and their levels, whitespace trimmed",
   parse_level_order(" visit = W0, W2 ,W12 ; stim=neg,stimulated;"),
   list(visit = c("W0", "W2", "W12"), stim = c("neg", "stimulated")))
eq("no order given", parse_level_order(""), list())
throws("an entry without '=' is an error", parse_level_order("visit W0,W2"),
       "not of the form")

vl <- c("W12", "W0", "W5", "Unassigned", "W5", "W110")
eq("default: natural order, Unassigned dropped", ordered_levels(vl),
   c("W0", "W5", "W12", "W110"))
eq("explicit levels first, the rest in natural order",
   ordered_levels(vl, explicit = c("W12", "W0")), c("W12", "W0", "W5", "W110"))
lv_w <- warns("an explicit level not in the data is a warning",
              ordered_levels(vl, explicit = c("W99", "W5")), "W99")
eq("...and is left out", lv_w, c("W5", "W0", "W12", "W110"))

# --- weighted_density_grid ---------------------------------------------------
cat("\nweighted_density_grid()\n")

set.seed(13)
dx <- rnorm(5000, 2, 0.3); dy <- rnorm(5000, 3.5, 0.3)
dg <- weighted_density_grid(dx, dy, xlim = c(0, 5), ylim = c(0, 5), n = 101L,
                            bandwidth = 0.03)
eq("the grid sums to 1", sum(dg$z), 1)
ok("grid axes span the limits",
   length(dg$x) == 101L && isTRUE(all.equal(range(dg$x), c(0, 5))) &&
     all(dim(dg$z) == c(101L, 101L)))
pk <- which(dg$z == max(dg$z), arr.ind = TRUE)[1, ]
ok("the peak is near the mode (rows are x, columns are y)",
   abs(dg$x[pk[1]] - 2) < 0.25 && abs(dg$y[pk[2]] - 3.5) < 0.25)
dz <- weighted_density_grid(dx, dy, w = rep(0, 5000), xlim = c(0, 5),
                            ylim = c(0, 5), n = 32L)
ok("zero weights give an all-zero grid", all(dz$z == 0))
d256 <- weighted_density_grid(dx, dy, xlim = c(0, 5), ylim = c(0, 5))
ok("the default grid is 256 x 256 and sums to 1",
   all(dim(d256$z) == c(256L, 256L)) && isTRUE(all.equal(sum(d256$z), 1)))
# Weights move the mass: two blobs, all weight on the second.
bx <- c(rep(1, 50), rep(4, 50)); by <- c(rep(1, 50), rep(4, 50))
dw <- weighted_density_grid(bx, by, w = rep(c(0, 1), each = 50),
                            xlim = c(0, 5), ylim = c(0, 5), n = 64L)
ok("weights decide where the mass goes",
   sum(dw$z[dw$x > 2.5, dw$y > 2.5]) > 0.999)

# --- separation_metrics ------------------------------------------------------
cat("\nseparation_metrics()\n")

if (is.null(.gl_knn_fun())) {
  cat("  SKIP  separation_metrics(): neither RANN nor FNN is installed\n")
} else {
  set.seed(17)
  blob <- function(n, cx, cy) cbind(rnorm(n, cx, 0.3), rnorm(n, cy, 0.3))
  co <- rbind(blob(200, 0, 0), blob(200, 10, 0), blob(100, 0, 10))
  lb <- rep(c("p", "q", "r"), c(200, 200, 100))
  sm <- separation_metrics(co, lb)
  eq("columns", names(sm), c("label", "cells_scored", "knn_purity", "isolated",
                             "silhouette", "group_mixing", "sample_mixing"))
  eq("first row is the macro average, then labels by size",
     sm$label, c("(macro average)", "p", "q", "r"))
  eq("perfectly separated blobs: purity 1 for every label and overall",
     sm$knn_purity, rep(1, 4))
  eq("...every cell isolated", sm$isolated, rep(1, 4))
  ok("...silhouette close to 1", all(sm$silhouette > 0.9))
  ok("...no mixing score without groups", all(is.na(sm$group_mixing)))

  # One blob: half of it is label p, half label q, interleaved at random.
  mixed <- separation_metrics(blob(400, 0, 0), sample(rep(c("p", "q"), 200)))
  ok("two labels in one cloud: purity near one half",
     abs(mixed$knn_purity[1] - 0.5) < 0.1)

  # Groups that mix at random score about 1; groups on separate sides of each
  # blob score about 0.
  grp_rand <- sample(rep(c("g1", "g2"), 250))
  grp_side <- ifelse(co[, 1] < ifelse(lb == "q", 10, 0), "g1", "g2")
  mr <- separation_metrics(co, lb, groups = grp_rand)
  ms <- separation_metrics(co, lb, groups = grp_side)
  ok("random groups mix about as chance would",
     abs(mr$group_mixing[1] - 1) < 0.2)
  ok("groups split by side barely mix", ms$group_mixing[1] < 0.3)

  st1 <- separation_metrics(co, lb, max_cells = 150L)
  eq("a cell budget balances the labels (50 each)", st1$cells_scored, c(150L, 50L, 50L, 50L))
}

# --- helpers added to the core after the plan --------------------------------
# Each is checked only when it exists, so a rename shows up as a SKIP line
# instead of stopping the whole file.
cat("\ntier_gaps(), is_viability_marker(), parse_split_reference()\n")

if (exists("tier_gaps")) {
  eq("tier gaps: gap * ratio^(deepest tier - tier)",
     tier_gaps(c(CD4 = 1L, CD8 = 1L, IFNy = 2L), gap = 2, ratio = 1.5),
     c(CD4 = 3, CD8 = 3, IFNy = 2))
  eq("ratio 1 gives one gap for all", tier_gaps(c(A = 1L, B = 3L), gap = 2),
     c(A = 2, B = 2))
  eq("no tiers, no gaps", tier_gaps(integer(0), 2), numeric(0))
  gp <- c(A = 3, B = 1.2)
  fgp <- gate_features(fv, fc, fk, sig, gap = gp)
  ok("per-marker gaps: each marker's sides still never overlap",
     separated(fgp, fk, "A") && separated(fgp, fk, "B"))
  eq("per-marker gaps: a cell at the cut lands on its marker's +/- gap",
     unname(gate_features(cbind(A = c(70, 70), B = c(10, 10)),
                          cbind(A = c(70, 70), B = c(10, 10)),
                          cbind(A = c(1L, 0L), B = c(1L, 0L)),
                          c(A = 5, B = 5), gap = gp)),
     cbind(c(3, -3), c(1.2, -1.2)))
} else {
  cat("  SKIP  tier_gaps() is not defined\n")
}

if (exists("is_viability_marker")) {
  eq("viability dyes are recognised by marker name, other markers are not",
     is_viability_marker(c("ALD", "Aqua LiveDead", "Zombie NIR", "7-AAD", "DAPI",
                           "PI", "CD4", "CD107a", "IL-2", "TNFa", NA)),
     c(rep(TRUE, 6), rep(FALSE, 5)))
  eq("an empty pattern recognises nothing",
     is_viability_marker(c("ALD", "Live"), pattern = ""), c(FALSE, FALSE))
} else {
  cat("  SKIP  is_viability_marker() is not defined\n")
}

if (exists("parse_split_reference")) {
  eq("one reference level per column, whitespace trimmed",
     parse_split_reference(" visit = W0 ; stim=neg;"),
     list(visit = "W0", stim = "neg"))
  eq("no reference given", parse_split_reference(NULL), list())
  throws("an entry without '=' is an error", parse_split_reference("visit"),
         "not of the form")
  throws("an entry without a level is an error",
         parse_split_reference("visit= "), "not of the form")
} else {
  cat("  SKIP  parse_split_reference() is not defined\n")
}

# --- review fixes: ancestors, two-axis gates, anchors, empty keys ----------
cat("\nreview fixes (2026-10-13)\n")

# A CD4 x CD8 quadrant under Live, and under CD4+ an IFNy gate drawn on
# IFNy x CD4 whose CD4 edge clips the CD4-dim tail of the CD4+ parent.
set.seed(11)
n_r <- 4000
x_r <- cbind(C4 = c(rnorm(2000, 120, 10), rnorm(2000, 50, 8)),
             C8 = c(rnorm(2000, 40, 8), rnorm(1000, 140, 10), rnorm(1000, 45, 8)),
             IF = rnorm(n_r, 60, 5))
x_r[sample(2000, 150), "IF"] <- rnorm(150, 110, 10)
cm_rv <- c(C4 = "CD4", C8 = "CD8", IF = "IFNy")
P_r <- "/L"
q4 <- 85; q8 <- 90
g_cd4 <- rect(c(C4 = q4, C8 = -Inf), c(C4 = Inf, C8 = q8))
g_cd8 <- rect(c(C4 = -Inf, C8 = q8), c(C4 = q4, C8 = Inf))
g_dp  <- rect(c(C4 = q4, C8 = q8), c(C4 = Inf, C8 = Inf))
g_dn  <- rect(c(C4 = -Inf, C8 = -Inf), c(C4 = q4, C8 = q8))
cd4_tail <- as.numeric(stats::quantile(x_r[in_rect(x_r, g_cd4), "C4"], 0.02))
g_if  <- rect(c(IF = 75, C4 = cd4_tail), c(IF = Inf, C4 = Inf))
mem_r <- cbind(`/L/CD4` = in_rect(x_r, g_cd4), `/L/CD8` = in_rect(x_r, g_cd8),
               `/L/ ` = in_rect(x_r, g_dp), `/L/  ` = in_rect(x_r, g_dn))
mem_r <- cbind(mem_r, `/L/CD4/IFNy` = mem_r[, "/L/CD4"] & in_rect(x_r, g_if))
gates_r <- gate_descendants(c("/L", colnames(mem_r)), P_r)
gates_r$geom <- list(g_cd4, g_cd8, g_dp, g_dn, g_if)[match(gates_r$path, colnames(mem_r))]
dec_r <- gate_decisions_sample(x_r, mem_r, gates_r, cm_rv, P_r)
ok("a gate never re-decides a marker its ancestor decides (CD4 edge of the IFNy box)",
   !any(dec_r$gate == "/L/CD4/IFNy" & dec_r$marker == "CD4"))
ok("...but still decides its own marker", any(dec_r$gate == "/L/CD4/IFNy" & dec_r$marker == "IFNy"))
eq("...and binds on one marker axis once the ancestor's axis is dropped",
   dec_r$n_axes[dec_r$gate == "/L/CD4/IFNy"], 1L)
cl_r <- gate_marker_calls(x_r, mem_r, dec_r, P_r)
in4 <- mem_r[, "/L/CD4"]
ok("every cell in the CD4+ gate stays CD4-positive", all(cl_r$calls[in4, "CD4"] == 1L))
ok("IFNy calls of CD4+ cells follow the gater's IFNy box exactly",
   all(cl_r$calls[in4, "IFNy"] == as.integer(mem_r[in4, "/L/CD4/IFNy"])))

# A lone CD4+CD8- rectangle (no CD8+, DP or DN gates): cells outside it are
# outside because of CD4 OR CD8, so their own values decide.
mem_l <- cbind(`/L/CD4` = mem_r[, "/L/CD4"])
gates_l <- gate_descendants(c("/L", "/L/CD4"), P_r); gates_l$geom <- list(g_cd4)
dec_l <- gate_decisions_sample(x_r, mem_l, gates_l, cm_rv, P_r)
eq("a CD4+CD8- box binds on two marker axes", unique(dec_l$n_axes), 2L)
cl_l <- gate_marker_calls(x_r, mem_l, dec_l, P_r)
dn_true <- x_r[, "C4"] < q4 & x_r[, "C8"] < q8
ok("double negatives outside the box are not called CD8-positive",
   all(cl_l$calls[dn_true, "CD8"] == 0L))
ok("...nor CD4-positive", all(cl_l$calls[dn_true, "CD4"] == 0L))
ok("CD8-high cells outside the box are CD8-positive",
   all(cl_l$calls[x_r[, "C8"] >= q8, "CD8"] == 1L))

# Anchors: per sample, the median of each marker over a combination's cells.
an_r <- combination_anchors(cl_r$values, cl_r$calls)
neg4 <- in4 & !mem_r[, "/L/CD4/IFNy"]
eq("anchor of IFNy-negative CD4+ cells is their own median",
   unique(an_r[neg4, "IFNy"]), stats::median(x_r[neg4, "IF"]))
pos4 <- in4 & mem_r[, "/L/CD4/IFNy"]
eq("anchor of IFNy-positive CD4+ cells is theirs",
   unique(an_r[pos4, "IFNy"]), stats::median(x_r[pos4, "IF"]))
cd8_only <- mem_r[, "/L/CD8"]; dn_only <- mem_r[, "/L/  "]
ok("CD8+ and double-negative cells are centred separately on CD4 (both are CD4-negative)",
   all(an_r[cd8_only, "CD4"] == stats::median(x_r[cd8_only, "C4"])) &&
   all(an_r[dn_only, "CD4"] == stats::median(x_r[dn_only, "C4"])))
# Moving the gate's cut without moving the cells leaves the features alone.
# No cell between the two cuts, so the calls are the same.
x_r[, "IF"] <- ifelse(x_r[, "IF"] < 80, pmin(x_r[, "IF"], 70), pmax(x_r[, "IF"], 90))
mem_r[, "/L/CD4/IFNy"] <- mem_r[, "/L/CD4"] & in_rect(x_r, g_if)
cl_r <- gate_marker_calls(x_r, mem_r, gate_decisions_sample(x_r, mem_r, gates_r, cm_rv, P_r), P_r)
an_r <- combination_anchors(cl_r$values, cl_r$calls)
g_if2 <- rect(c(IF = 85, C4 = cd4_tail), c(IF = Inf, C4 = Inf))
mem_r2 <- mem_r; mem_r2[, "/L/CD4/IFNy"] <- mem_r[, "/L/CD4"] & in_rect(x_r, g_if2)
gates_r2 <- gates_r; gates_r2$geom <- list(g_cd4, g_cd8, g_dp, g_dn, g_if2)[match(gates_r$path, colnames(mem_r))]
same_calls <- function(a, b) all(a$calls == b$calls, na.rm = TRUE)
cl_r2 <- gate_marker_calls(x_r, mem_r2, gate_decisions_sample(x_r, mem_r2, gates_r2, cm_rv, P_r), P_r)
if (same_calls(cl_r, cl_r2)) {
  ok("a cut moved within the gap changes no anchor (features follow the data)",
     isTRUE(all.equal(an_r, combination_anchors(cl_r2$values, cl_r2$calls))))
} else {
  cat("  SKIP  moved cut changed a call in this draw\n")
}

set.seed(9)
ok(".gl_auc() does not overflow with 50000 x 50000 cells",
   is.finite(.gl_auc(rnorm(5e4, 2), rnorm(5e4))))

# A sample whose population has no deciding gate keeps every cell.
st0 <- stratified_sample(calls_key(matrix(NA_integer_, 5, 0)), 2)
eq("an empty stratification key is a key like any other (cap 2 keeps 2)",
   length(st0$idx), 2L)
eq("...weighted back to the 5 cells", sum(st0$weight), 5)

# The missing-decision fill weighs cells, not kept rows.
calls_f <- cbind(CD4 = rep(1L, 20), IFNy = c(rep(1L, 12), rep(NA_integer_, 8)))
w_f <- c(rep(1, 12), rep(9, 8))      # the undecided rows stand for 9 cells each
fm <- fill_missing_decisions(calls_f, c(CD4 = 1L, IFNy = 2L), rep(c("a", "b"), each = 10),
                             weights = w_f)
ok("a marker decided for 14 % of the group's cells is not filled",
   all(is.na(fm$calls[13:20, "IFNy"])))
fm2 <- fill_missing_decisions(calls_f, c(CD4 = 1L, IFNy = 2L), rep(c("a", "b"), each = 10))
ok("...while unweighted rows (60 %) would have filled it", all(fm2$calls[13:20, "IFNy"] == 0L))

if (!is.null(.gl_knn_fun())) {
  set.seed(4)
  co_m <- rbind(cbind(rnorm(400, 0), rnorm(400, 0)), cbind(rnorm(400, 20), rnorm(400, 0)))
  lb_m <- rep(c("p", "q"), each = 400)
  smp_split <- c(ifelse(co_m[1:400, 1] < 0, "s1", "s2"), rep(c("s1", "s2"), 200))
  sm <- separation_metrics(co_m, lb_m, samples = smp_split)
  ok("sample_mixing is low for an island split by sample",
     sm$sample_mixing[sm$label == "p"] < 0.3)
  ok("...and near 1 for a mixed one", abs(sm$sample_mixing[sm$label == "q"] - 1) < 0.2)
} else {
  cat("  SKIP  sample_mixing: neither RANN nor FNN is installed\n")
}

cat("\npool_gate_layout()\n")

if (exists("pool_gate_layout")) {
  # The two MTB-like samples, the second without the 5 CD8+ TNFa+ cells (h),
  # sampled at 4 cells per call pattern as the extraction would.
  # `group` is extra: pool_gate_layout() ignores it, the checks below use it.
  mk_per_sample <- function(smp, cm, keep, seed) {
    calls <- cm$calls[keep, , drop = FALSE]
    key <- calls_key(calls)
    set.seed(seed); st <- stratified_sample(key, 4)
    sub <- function(m) m[keep, , drop = FALSE][st$idx, , drop = FALSE]
    list(calls = calls[st$idx, , drop = FALSE], cuts = sub(cm$cuts),
         anchors = sub(combination_anchors(cm$values, cm$calls)),
         values = sub(cm$values), mtx = sub(cm$values),
         weight = st$weight, key = key[st$idx],
         key_counts = c(table(key)), population_cells = sum(keep),
         group = smp$group[keep][st$idx])
  }
  keep2 <- s2$group != "h"
  per <- list(s1 = mk_per_sample(s1, c1, rep(TRUE, 140), 1),
              s2 = mk_per_sample(s2, c2, keep2, 2))
  decs <- rbind(cbind(sample = "s1", d1), cbind(sample = "s2", d2))
  label_of <- function(g) MTB_CELLS$label[match(g, MTB_CELLS$group)]
  # The true counts over ALL cells, with combinations under `min_cells` in
  # total pooled, as list(s1 = table, s2 = table, levels).
  true_counts <- function(min_cells) {
    l1 <- label_of(s1$group); l2 <- label_of(s2$group[keep2])
    tot <- table(c(l1, l2))
    rare <- names(tot)[tot < min_cells]
    l1[l1 %in% rare] <- .GL_RARE; l2[l2 %in% rare] <- .GL_RARE
    list(s1 = table(l1), s2 = table(l2), levels = union(l1, l2))
  }
  check_counts <- function(gl, tc) {
    all(vapply(c("s1", "s2"), function(s) {
      cs <- gl$counts[gl$counts$sample == s, ]
      truth <- ifelse(cs$combination %in% names(tc[[s]]),
                      as.numeric(tc[[s]][cs$combination]), 0)
      setequal(cs$combination, tc$levels) && all(cs$cells == truth)
    }, logical(1)))
  }
  smp_of_row <- rep(c("s1", "s2"), c(length(per$s1$weight), length(per$s2$weight)))

  gl0 <- pool_gate_layout(per, decs, min_combination_cells = 0)
  eq("tiers from the pooled decisions", gl0$tiers,
     c(CD4 = 1L, CD8 = 1L, IFNy = 2L, TNFa = 2L))
  eq("one row per kept cell, samples in order",
     nrow(gl0$calls), length(smp_of_row))
  eq("kept cells carry their gater's combination", gl0$combination,
     label_of(c(per$s1$group, per$s2$group)))
  ok("counts cover ALL cells of every combination, zero-filled per sample",
     check_counts(gl0, true_counts(0)))
  ok("a combination a sample lacks is a zero row, not a missing one",
     any(gl0$counts$sample == "s2" & gl0$counts$combination == "CD8+ | TNFa+" &
           gl0$counts$cells == 0))
  eq("per sample, the counts add up to the sample's cells",
     as.vector(tapply(gl0$counts$cells, gl0$counts$sample, sum)), c(140, 135))
  eq("population_cells is carried per sample",
     as.vector(tapply(gl0$counts$population_cells, gl0$counts$sample, unique)),
     c(140, 135))
  w_ok <- TRUE
  for (r in which(gl0$counts$cells > 0)) {
    sel <- smp_of_row == gl0$counts$sample[r] &
      gl0$combination == gl0$counts$combination[r]
    if (!isTRUE(all.equal(sum(gl0$weight[sel]), gl0$counts$cells[r]))) w_ok <- FALSE
  }
  ok("kept weights add up to the count of each sample and combination", w_ok)
  gl25 <- pool_gate_layout(per, decs, min_combination_cells = 25)
  ok("combinations under 25 cells in total are pooled, in the counts too",
     check_counts(gl25, true_counts(25)) &&
       "other (rare)" %in% gl25$counts$combination)
  ok("no decision was missing, so nothing was filled", is.null(gl0$filled))
} else {
  cat("  SKIP  pool_gate_layout() is not defined\n")
}

# --- umap_style.R: palette and combination colours ---------------------------
cat("\numap_style.R\n")

if (!has_style) {
  cat("  SKIP  src/visualization/umap_style.R is not present\n")
} else {
  is_colour <- function(x) {
    isTRUE(tryCatch({ grDevices::col2rgb(x); TRUE }, error = function(e) FALSE))
  }
  # Chroma proxy in 0..255: a grey has red = green = blue.
  spread <- function(x) {
    rgb <- grDevices::col2rgb(x)
    apply(rgb, 2, function(v) max(v) - min(v))
  }
  lightness <- function(x) {
    grDevices::convertColor(t(grDevices::col2rgb(x)) / 255, from = "sRGB",
                            to = "Luv")[, 1]
  }
  hue <- function(x) {
    luv <- grDevices::convertColor(t(grDevices::col2rgb(x)) / 255,
                                   from = "sRGB", to = "Luv")
    (atan2(luv[, 3], luv[, 2]) * 180 / pi) %% 360
  }
  hue_dist <- function(a, b) { d <- abs(a - b) %% 360; pmin(d, 360 - d) }

  eq("umap_muted_palette(0) is empty", umap_muted_palette(0), character(0))
  for (n in c(1, 5, 12, 13, 30)) {
    p <- umap_muted_palette(n)
    ok(sprintf("umap_muted_palette(%d): %d distinct valid colours", n, n),
       length(p) == n && is_colour(p) && !anyDuplicated(toupper(p)))
  }
  ok("a longer palette extends a shorter one",
     identical(umap_muted_palette(30)[1:12], umap_muted_palette(12)))

  labs_c <- c("CD4+ | none", "CD4+ | IFNy+", "CD4+ | IFNy+ TNFa+", "CD4+ | TNFa+",
              "CD8+ | none", "CD8+ | IFNy+", "CD8+ | IFNy+ TNFa+",
              "CD4+ CD8+", "CD4- CD8-", "other (rare)")
  cc <- combination_colours(labs_c)
  ok("one named colour per distinct label",
     setequal(names(cc), labs_c) && length(cc) == length(labs_c))
  ok("all valid colours, all different", is_colour(cc) &&
       !anyDuplicated(toupper(cc)))
  ok("repeated labels do not change the result",
     identical(combination_colours(c(labs_c, labs_c)), cc))
  cc_rev <- combination_colours(rev(labs_c))
  ok("the colours depend on the set of labels, not their order",
     identical(cc_rev[sort(names(cc_rev))], cc[sort(names(cc))]))
  ok("'other (rare)' is grey", spread(cc[["other (rare)"]]) <= 6)
  ok("an all-negative first tier ('CD4- CD8-') is a neutral grey",
     spread(cc[["CD4- CD8-"]]) <= 15)
  fam4 <- grep("^CD4\\+ \\|", labs_c, value = TRUE)
  fam8 <- grep("^CD8\\+ \\|", labs_c, value = TRUE)
  ok("a group's '| none' is the lightest colour of its family (CD4+)",
     all(lightness(cc[["CD4+ | none"]]) >
           lightness(cc[setdiff(fam4, "CD4+ | none")])))
  ok("...and of the CD8+ family",
     all(lightness(cc[["CD8+ | none"]]) >
           lightness(cc[setdiff(fam8, "CD8+ | none")])))
  h4 <- hue(cc[["CD4+ | none"]]); h8 <- hue(cc[["CD8+ | none"]])
  ok("each family shares a hue: CD4+ colours sit nearer CD4+ than CD8+",
     all(hue_dist(hue(cc[fam4]), h4) < hue_dist(hue(cc[fam4]), h8)))
  ok("...and CD8+ colours nearer CD8+ than CD4+",
     all(hue_dist(hue(cc[fam8]), h8) < hue_dist(hue(cc[fam8]), h4)))
  ok("no labels: an empty named vector",
     length(combination_colours(character(0))) == 0L)
}

# =============================================================================
# umap_style.R and umap_split_plots.R: colour distances, change maps, marker
# colour limits. Kept as one section at the end so that it can change without
# touching the sections above.
# =============================================================================
cat("\numap_style.R colour distances, umap_split_plots.R change tables\n")

split_file <- file.path(base_dir, "src", "visualization", "umap_split_plots.R")
has_split <- has_style && file.exists(split_file) &&
  isTRUE(tryCatch({ source(split_file); TRUE }, error = function(e) {
    cat("NOTE: could not source umap_split_plots.R (", conditionMessage(e),
        "); its tests are skipped.\n", sep = "")
    FALSE
  }))

if (!has_style) {
  cat("  SKIP  src/visualization/umap_style.R is not present\n")
} else {
  # Smallest CIE-Lab distance between any two of the colours.
  min_de <- function(x) {
    lab <- grDevices::convertColor(t(grDevices::col2rgb(x)) / 255,
                                   from = "sRGB", to = "Lab")
    d <- as.matrix(stats::dist(lab)); diag(d) <- Inf
    min(d)
  }
  base_de <- min_de(umap_muted_palette(12))
  ok("umap_muted_palette(13..40): no two colours closer than the base twelve",
     all(vapply(13:40, function(n) min_de(umap_muted_palette(n)), 0) >= base_de - 0.5))
  ok("umap_muted_palette(): a longer palette extends a shorter one (beyond 12)",
     identical(umap_muted_palette(40)[1:20], umap_muted_palette(20)))
  # The single-tier case of the review: TNF+ and IFN+ islands came out in the
  # same green once there were more than twelve labels.
  flat16 <- c("none", "CD107a+", "IFNy+", "IL-2+", "TNFa+",
              combn(c("CD107a+", "IFNy+", "IL-2+", "TNFa+"), 2, paste, collapse = " "),
              combn(c("CD107a+", "IFNy+", "IL-2+", "TNFa+"), 3, paste, collapse = " "),
              "CD107a+ IFNy+ IL-2+ TNFa+")
  ok("16 labels without ' | ' (a single-tier population) are at least dE 13 apart",
     min_de(combination_colours(flat16)) >= 13)

  cyt4 <- c("CD107a", "IFNy", "IL-2", "TNFa")
  family <- function(lin) {
    pos <- unlist(lapply(1:4, function(k) combn(cyt4, k, function(p)
      paste0(p, "+", collapse = " "))))
    c(paste(lin, "| none"), paste(lin, "|", pos))
  }
  for (n_fam in 1:3) {
    lins <- c("CD4+", "CD8+", "CD56+")[seq_len(n_fam)]
    cc16 <- combination_colours(c(unlist(lapply(lins, family)), "CD4+ CD8+",
                                  "CD4- CD8-", "other (rare)"))
    within <- vapply(lins, function(l)
      min_de(cc16[startsWith(names(cc16), paste(l, "|"))]), 0)
    ok(sprintf("%d famil%s of 16 labels: colours within a family at least dE 9.5 apart (min %.1f)",
               n_fam, if (n_fam == 1) "y" else "ies", min(within)),
       all(within >= 9.5))
    fam_l <- function(l) {
      f <- cc16[startsWith(names(cc16), paste(l, "|"))]
      grDevices::convertColor(t(grDevices::col2rgb(f)) / 255, from = "sRGB",
                              to = "Lab")[, 1]
    }
    ok(sprintf("...and each family's '| none' is still its lightest colour (%d)", n_fam),
       all(vapply(lins, function(l) which.max(fam_l(l)) == 1L, logical(1))))
  }
}

if (!has_split) {
  cat("  SKIP  src/visualization/umap_split_plots.R is not present\n")
} else {
  fq <- data.frame(column = "visit", level = rep(c("V2", "V0"), each = 4),
                   combination = rep(c("up", "gone", "never", "same"), 2),
                   mean_pct = c(10, 0, 0, 2, 2.5, 1, 0, 2), sd_pct = NA_real_,
                   n_samples = rep(c(5L, 3L), each = 4), stringsAsFactors = FALSE)
  ct <- combination_change_table(fq, "V2", "V0", pseudo_pct = 0.5)
  eq("combination_change_table(): columns",
     names(ct), c("column", "level", "reference", "combination", "mean_pct_level",
                  "mean_pct_reference", "log2_fold_change", "n_samples_level",
                  "n_samples_reference", "pseudo_pct"))
  eq("log2 fold change of the means plus the pseudo-count",
     ct$log2_fold_change[1], log2(10.5 / 3))
  ok("a combination gone at the level is finite and negative",
     is.finite(ct$log2_fold_change[2]) && ct$log2_fold_change[2] < 0)
  ok("a combination absent at both levels is NA", is.na(ct$log2_fold_change[3]))
  eq("an unchanged combination is 0", ct$log2_fold_change[4], 0)
  ok("sample counts per side",
     all(ct$n_samples_level == 5L) && all(ct$n_samples_reference == 3L))

  # The population each sample's percentages are taken of: the weighted
  # plotted cells without counts, population_cells with them.
  pc <- data.frame(sample = rep(c("a", "b"), c(3, 2)),
                   combination = c("x", "x", "y", "x", "y"),
                   weight = c(10, 10, 20, 5, 5), stringsAsFactors = FALSE)
  pm <- .sp_percent_matrix(pc, NULL, c("x", "y"))
  eq("percent matrix without counts: rows sum to 100", unname(rowSums(pm)), c(100, 100))
  eq("...and the population is the weighted cells", attr(pm, "population"),
     c(a = 40, b = 10))
  cnt <- data.frame(sample = c("a", "a", "b"), combination = c("x", "z", "y"),
                    cells = c(30, 5, 4), population_cells = c(200, 200, 50),
                    stringsAsFactors = FALSE)
  pm2 <- .sp_percent_matrix(pc, cnt, c("x", "y"))
  eq("percent matrix with counts: % of population_cells",
     unname(pm2[, "x"]), c(15, 0))
  eq("...and the population is population_cells", attr(pm2, "population"),
     c(a = 200, b = 50))

  # Marker colour limits: the full finite range, as color_limits() in
  # umap_from_flowjo.R gives the per-file marker plots.
  eq(".sp_colour_limits(): full finite range",
     .sp_colour_limits(c(3, NA, -1, Inf, 7)), c(-1, 7))
  ok(".sp_colour_limits(): NULL when constant or empty",
     is.null(.sp_colour_limits(c(2, 2))) && is.null(.sp_colour_limits(c(NA, NaN))))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    sc <- umap_marker_colour(c(0, 10), "M")
    ok("umap_marker_colour(): the viridis ends, squished beyond the limits",
       identical(toupper(sc$map(c(-5, 0, 10, 50))),
                 toupper(substr(viridisLite::viridis(256)[c(1, 1, 256, 256)], 1, 7))))
    dv <- umap_diverging_colour(3, "Fold change")
    ok("umap_diverging_colour(): near-white at 0, squished beyond the limit",
       toupper(dv$map(0)) == "#F6F4F1" && identical(dv$map(10), dv$map(3)))
    ok("...and absent at both levels (NA) is a grey unlike the centre",
       dv$map(NA_real_) != dv$map(0) && min_de(c(dv$map(NA_real_), dv$map(0))) > 15)
  } else {
    cat("  SKIP  ggplot2 is not installed; the colour scales are not checked\n")
  }
}

cat("\nfit subsets: run_umap(fit_matrix =) and supercell_fit_matrix()\n")
source(file.path(base_dir, "src", "visualization", "umap_supercells.R"))
source(file.path(base_dir, "src", "R_shared", "stats_utils.R"))
set.seed(7)
fx <- rbind(matrix(rnorm(600, 0), ncol = 3), matrix(rnorm(600, 6), ncol = 3))
colnames(fx) <- c("A", "B", "C")
fs <- rep(c("s1", "s2"), each = 200)
if (requireNamespace("uwot", quietly = TRUE)) {
  fm <- fx[seq(1, 400, by = 4), ]
  uc <- suppressMessages(run_umap(fx, n_neighbors = 10, fit_matrix = fm, seed = 1))
  ok("run_umap(fit_matrix): one finite coordinate pair per cell",
     identical(dim(uc), c(400L, 2L)) && all(is.finite(uc)))
  ok("...and projected cells keep the two groups apart",
     abs(mean(uc[1:200, 1]) - mean(uc[201:400, 1])) +
       abs(mean(uc[1:200, 2]) - mean(uc[201:400, 2])) > 2)
  seen <- NULL
  uc2 <- suppressMessages(run_umap(fx, n_neighbors = 10, fit_matrix = fm, seed = 1,
                                   on_fit = function(m) seen <<- m))
  ok("run_umap(on_fit): called with the fitted model (fit rows only)",
     inherits(seen, "list") && identical(nrow(seen$embedding), nrow(fm)))
  ok("...and the model rides along on the coordinates",
     identical(nrow(attr(uc2, "model")$embedding), nrow(fm)))
  # A PCA start on values spread like FlowJo's 0-256 scale: unscaled, the
  # points start ~100 units apart and never come together.
  wide <- rbind(matrix(rnorm(1500, 60, 8), ncol = 3), matrix(rnorm(1500, 180, 8), ncol = 3))
  colnames(wide) <- c("A", "B", "C")
  up <- suppressMessages(run_umap(wide, n_neighbors = 15, init = "pca", seed = 1))
  ok("run_umap(init = \"pca\"): a scaled start, so a normal-width layout",
     max(apply(up, 2, function(v) diff(range(v)))) < 60)
  ok(".warn_spread(): warns on a layout hundreds of units wide, not on a normal one",
     inherits(tryCatch(.warn_spread(cbind(c(-150, 150), c(0, 1))),
                       warning = function(w) w), "warning") &&
       !inherits(tryCatch(.warn_spread(cbind(c(-10, 12), c(0, 1))),
                          warning = function(w) w), "warning"))
  # Resume: a projection from the saved model reproduces the uninterrupted
  # run (the projection is seeded per chunk, whatever the fit consumed).
  ua <- suppressMessages(run_umap(fx, n_neighbors = 10, fit_matrix = fm, seed = 1))
  ub <- suppressMessages(run_umap(fx, n_neighbors = 10, fit_matrix = fm, seed = 1,
                                  model = attr(ua, "model")))
  ok("run_umap(model =): projecting from a saved model gives the same coordinates",
     identical(unclass(ua)[, 1:2], unclass(ub)[, 1:2]))
  um <- suppressMessages(run_umap(fx, n_neighbors = 10, model = attr(ua, "model"),
                                  seed = 1))
  ok("run_umap(model =) without a fit set: every row projected, nothing refitted",
     identical(dim(unclass(um)[, 1:2]), c(400L, 2L)) && all(is.finite(um)) &&
       identical(attr(um, "model")$embedding, attr(ua, "model")$embedding))
  ok("...and refuses a model fitted on another fit set",
     inherits(try(suppressMessages(run_umap(fx, fit_matrix = fm[1:50, ],
                                            model = attr(ua, "model"))),
                  silent = TRUE), "try-error"))
  ok("...refuses a fit matrix with other columns",
     inherits(try(run_umap(fx, fit_matrix = fm[, 3:1]), silent = TRUE), "try-error"))
  ok("...and refuses it together with fit_idx",
     inherits(try(run_umap(fx, fit_matrix = fm, fit_idx = 1:10), silent = TRUE),
              "try-error"))
} else {
  cat("  SKIP  uwot is not installed; run_umap(fit_matrix) is not checked\n")
}
if (requireNamespace("SuperCellCyto", quietly = TRUE)) {
  # A third sample of 10 cells is too small to group (< 2 * gamma) and is
  # passed through cell by cell.
  sx <- rbind(fx, matrix(rnorm(30, 3), ncol = 3, dimnames = list(NULL, colnames(fx))))
  sc <- suppressWarnings(supercell_fit_matrix(sx, c(fs, rep("s3", 10)), gamma = 20))
  ok("supercell_fit_matrix(): the marker columns, in order",
     identical(colnames(sc$matrix), colnames(fx)))
  ok("...about cells / gamma supercells from the grouped samples",
     sc$n_supercells >= 10 && sc$n_supercells <= 40)
  eq("...the small sample's cells are passed through", sc$passed_through, 10L)
  eq("...fit rows = supercells + passed-through cells", nrow(sc$matrix),
     sc$n_supercells + 10L)
  # s1 sits around 0 and s2 around 6: a supercell mixing them would average 3.
  rm_sc <- rowMeans(sc$matrix[seq_len(sc$n_supercells), , drop = FALSE])
  ok("...each supercell stays within one sample (none averages s1 and s2)",
     all(rm_sc < 2 | rm_sc > 4))
  ok("...every grouped cell has a supercell, the passed-through ones none",
     !anyNA(sc$cell_supercell[1:400]) && all(is.na(sc$cell_supercell[401:410])))
  eq("...one supercells row per fit row, sizes summing to the cells",
     c(nrow(sc$supercells), sum(sc$supercells$cells)), c(nrow(sc$matrix), 410L))
  g1 <- sc$supercells$sample == "s1" & !is.na(sc$supercells$supercell)
  k1 <- sc$supercells$supercell[g1][1]
  ss <- c(fs, rep("s3", 10))
  eq("...a supercell's row is the mean of the cells mapped to it",
     unname(sc$matrix[which(g1)[1], ]),
     unname(colMeans(sx[which(ss == "s1" & sc$cell_supercell %in% k1), ,
                        drop = FALSE])))
  if (.Platform$OS.type == "unix") {
    sc2 <- suppressWarnings(supercell_fit_matrix(sx, c(fs, rep("s3", 10)),
                                                 gamma = 20, n_workers = 2))
    ok("...the same fit matrix with 2 workers (rows in sample order)",
       identical(sc2$matrix, sc$matrix))
  }
  ok("supercell_fit_matrix(): gamma must be above 1",
     inherits(try(supercell_fit_matrix(fx, fs, gamma = 1), silent = TRUE), "try-error"))
} else {
  cat("  SKIP  SuperCellCyto is not installed; supercell_fit_matrix() is not checked\n")
  ok("require_supercellcyto(): says how to install it",
     grepl("install_github", tryCatch(require_supercellcyto(),
                                      error = conditionMessage)))
}

cat("\nbackground subtraction: parse_background(), background_pairs(), subtract_background()\n")
eq("parse_background(): column and level",
   parse_background(" stimulation = unstimulated "),
   list(column = "stimulation", level = "unstimulated"))
ok("...NULL for none, an error for a malformed one",
   is.null(parse_background(NULL)) &&
     inherits(try(parse_background("stimulation"), silent = TRUE), "try-error"))
ann_bg <- data.frame(patient = c("P1", "P1", "P1", "P2", "P2", "P3"),
                     visit = c("V0", "V0", "V0", "V1", "V1", "V1"),
                     stim = c("s", "s", "neg", "s", "neg", "s"),
                     stringsAsFactors = FALSE)
smp <- c("a", "b", "a_neg", "c", "c_neg", "d")
bp <- suppressWarnings(background_pairs(smp, ann_bg, "stim", "neg",
                                        c("patient", "visit")))
eq("background_pairs(): each stimulated tube to the neg of its patient and visit",
   bp, c(a = "a_neg", b = "a_neg", c = "c_neg"))
ok("...warns about a sample without a background (d)",
   inherits(tryCatch(background_pairs(smp, ann_bg, "stim", "neg",
                                      c("patient", "visit")),
                     warning = function(w) w), "warning"))
ok("...and stops when no background sample is left (filtered out)",
   inherits(try(background_pairs(smp[-c(3, 5)], ann_bg[-c(3, 5), ], "stim", "neg",
                                 c("patient", "visit")), silent = TRUE), "try-error"))
pct_bg <- rbind(a = c(x = 5, y = 1), b = c(x = 2, y = 0.5), a_neg = c(x = 3, y = 0.8),
                c = c(x = 1, y = 2), c_neg = c(x = 0.5, y = 2.5))
attr(pct_bg, "population") <- c(a = 100, b = 100, a_neg = 100, c = 100, c_neg = 100)
sb <- subtract_background(pct_bg, bp)
eq("subtract_background(): sample minus its background, negatives 0",
   as.vector(sb), c(2, 0, 0.5, 0.2, 0, 0))
eq("...only the stimulated samples remain", rownames(sb), c("a", "b", "c"))
sb2 <- subtract_background(pct_bg, bp, subtract = "x")
ok("...a combination not subtracted is NA (not a response)",
   all(is.na(sb2[, "y"])) && identical(as.vector(sb2[, "x"]), c(2, 0, 0.5)))
eq("response_combinations(): cytokine-positive combinations only",
   response_combinations(c("CD4+ | none", "CD4+ | IFNy+", "CD4- CD8-",
                           "CD8+ | TNFa+ CD107a+", "other (rare)")),
   c("CD4+ | IFNy+", "CD8+ | TNFa+ CD107a+"))
eq("...and every label when they are not gate combinations",
   response_combinations(c("CD4", "CD8")), c("CD4", "CD8"))

cat("\ncheckpoints: ckpt_key(), ckpt_save(), ckpt_load()\n")
ff_env <- new.env()
ff_ok <- isTRUE(tryCatch({
  sys.source(file.path(base_dir, "src", "visualization", "umap_from_flowjo.R"),
             envir = ff_env)
  TRUE
}, error = function(e) FALSE))
if (ff_ok) {
  k1 <- ff_env$ckpt_key(list("extract", a = 1, b = "x"))
  ok("ckpt_key(): the same settings give the same key",
     identical(k1, ff_env$ckpt_key(list("extract", a = 1, b = "x"))))
  ok("...and any change gives another", !identical(k1,
     ff_env$ckpt_key(list("extract", a = 2, b = "x"))))
  cd <- file.path(tempdir(), "ckpt_test")
  unlink(cd, recursive = TRUE)
  invisible(capture.output(ff_env$ckpt_save(cd, "extract", k1, list(m = 1:3))))
  got <- NULL
  invisible(capture.output(got <- ff_env$ckpt_load(cd, "extract", k1)))
  eq("ckpt_load(): the saved value back under its key", got, list(m = 1:3))
  got2 <- "unset"
  invisible(capture.output(got2 <- ff_env$ckpt_load(cd, "extract", "other")))
  ok("...NULL (recompute) under another key", is.null(got2))
  got3 <- "unset"
  invisible(capture.output(got3 <- ff_env$ckpt_load(cd, "supercells", k1)))
  ok("...NULL when the stage was never saved", is.null(got3))
  ok("ckpt_save(): no .part file left behind",
     !length(list.files(file.path(cd, "checkpoints"), pattern = "\\.part$")))
  unlink(cd, recursive = TRUE)
} else {
  cat("  SKIP  umap_from_flowjo.R could not be sourced here\n")
}

cat("\n=== ", .pass, " passed, ", .fail, " failed ===\n", sep = "")

# --- summary -----------------------------------------------------------------
if (.fail > 0) quit(status = 1)
