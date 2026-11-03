#!/usr/bin/env Rscript
# test_data/test_drift_compass.R — unit tests for step 16 (src/harmonization/drift/drift_compass.R):
# Boolean coding, COMPASS subset names, the category filter and null rebuild, FS/PFS per draw,
# the subset model frame and the hypothesis summaries. Base R only.
#   Rscript test_data/test_drift_compass.R
cat("=== drift COMPASS tests ===\n\n")
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_drift_compass.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) base_dir <- getwd()
source(file.path(base_dir, "src", "R_shared", "fcs_keywords.R"))
for (m in list.files(file.path(base_dir, "src", "harmonization", "drift"), "\\.R$", full.names = TRUE)) source(m)

n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  ok  ", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL", label, "\n") }
}
mk <- c("IFNg", "IL2", "TNFa", "CD107a")

# --- subset names -----------------------------------------------------------------------
cat("subset names\n")
nm <- compass_subset_names(mk)
check("16 subsets for 4 markers", length(nm) == 16 && !anyDuplicated(nm))
check("all-negative subset last", nm[16] == "!IFNg&!IL2&!TNFa&!CD107a" && nm[16] == null_subset(mk))
check("all-positive subset first", nm[1] == "IFNg&IL2&TNFa&CD107a")
check("degrees run from 4 down to 0", identical(subset_degree(nm), c(4L, rep(3L, 4), rep(2L, 6), rep(1L, 4), 0L)))
check("COMPASS name format", "IFNg&!IL2&TNFa&!CD107a" %in% nm)

# --- Boolean coding ----------------------------------------------------------------------
cat("boolean_counts\n")
# 8 events: 6 in the parent; event 7 is a marker-positive event outside the parent
parent <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
hits <- list(IFNg = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE),
             IL2 = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
             TNFa = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
             CD107a = rep(FALSE, 8))
bc <- boolean_counts(parent, hits)
check("names in COMPASS order", identical(names(bc), nm))
check("counts sum to the parent", sum(bc) == sum(parent))
check("triple producer counted once", bc[["IFNg&IL2&TNFa&!CD107a"]] == 1)
check("IFNg+TNFa double", bc[["IFNg&!IL2&TNFa&!CD107a"]] == 1)
check("IFNg single", bc[["IFNg&!IL2&!TNFa&!CD107a"]] == 1)
check("TNFa single", bc[["!IFNg&!IL2&TNFa&!CD107a"]] == 1)
check("null holds the two negative parent events", bc[[null_subset(mk)]] == 2)
check("events outside the parent ignored", bc[["IFNg&IL2&!TNFa&!CD107a"]] == 0)
check("marginals equal the marker gate counts within the parent",
      identical(unname(boolean_marginals(bc, mk)), c(3, 1, 3, 0)))
check("unnamed hits rejected", inherits(try(boolean_counts(parent, unname(hits)), silent = TRUE), "try-error"))
check("length mismatch rejected", inherits(try(boolean_counts(parent[-1], hits), silent = TRUE), "try-error"))

# --- gate name mapping -------------------------------------------------------------------
cat("canon_cytokine_gate\n")
cfg <- list(cytokine_markers = mk, markers = list(IFNg = c("IFN", "INF", "^PE-A$"), IL2 = c("IL-?2", "^APC-A$"),
                                                   TNFa = c("TNF", "^Pacific Blue-A$"), CD107a = c("CD107", "^FITC-A$")))
check("FlowJo spellings map onto the config names",
      identical(canon_cytokine_gate(c("IFNy", "INFy", "IL-2", "TNF", "CD107a", "Foo"), cfg),
                c("IFNg", "IFNg", "IL2", "TNFa", "CD107a", "Foo")))

# --- filter and null rebuild ---------------------------------------------------------------
cat("compass_filter / collapse_to_kept\n")
n_s <- matrix(0, 4, 16, dimnames = list(paste0("t", 1:4), nm))
n_s[, "IFNg&!IL2&!TNFa&!CD107a"] <- c(10, 20, 30, 1)   # > 5 in 3 tubes: kept
n_s[, "IFNg&!IL2&TNFa&!CD107a"] <- c(10, 20, 1, 1)     # > 5 in 2 tubes: dropped
n_s[, null_subset(mk)] <- 1000
keep <- compass_filter(n_s, null_subset(mk))
check("filter keeps subsets with > 5 cells in > 2 tubes", identical(keep, "IFNg&!IL2&!TNFa&!CD107a"))
check("filter never keeps the null", !null_subset(mk) %in% keep)
pc <- rowSums(n_s)
cm <- collapse_to_kept(n_s, keep, null_subset(mk), pc)
check("null column last", colnames(cm)[ncol(cm)] == null_subset(mk) && ncol(cm) == 2)
check("dropped cells move into the null", isTRUE(all.equal(unname(cm[, 2]), unname(pc - n_s[, keep]))) && all(rowSums(cm) == pc))
check("negative null rejected", inherits(try(collapse_to_kept(n_s, keep, null_subset(mk), pc - 2000), silent = TRUE), "try-error"))

# --- FS / PFS ------------------------------------------------------------------------------------
cat("compass_score_draws\n")
deg <- c(3L, 2L, 1L, 0L)   # 3 non-null categories and the null, n = 3 markers
g <- array(0, c(2, 4, 2)); g[1, , 1] <- c(1, 1, 1, 0); g[1, , 2] <- c(1, 0, 0, 0); g[2, , ] <- 0
sc <- compass_score_draws(g, deg, 3)
check("FS divides by 2^n - 1", isTRUE(all.equal(sc$FS[1, ], c(3 / 7, 1 / 7))) && all(sc$FS[2, ] == 0))
pfs1 <- (3 / choose(3, 3) + 2 / choose(3, 2) + 1 / choose(3, 1)) / 3 * 2 / 4
check("PFS weights degree / choose(n, degree)", isTRUE(all.equal(sc$PFS[1, 1], pfs1)) && isTRUE(all.equal(sc$PFS[1, 2], 3 / choose(3, 3) / 3 * 2 / 4)))
# matches COMPASS's own formula on mean_gamma for a single draw
mg <- g[, , 1]
check("same as COMPASS PolyfunctionalityScore on mean_gamma",
      isTRUE(all.equal(unname(sc$PFS[, 1]), apply(mg, 1, function(row) sum(row * deg / choose(3, deg)) / 3 * (2 / 4)))))

# --- mean model over seeds, pooled draws, shared controls ---------------------------------------
cat("mean model and control sharing\n")
f <- function(m, g) list(fit = list(mean_gamma = m, gamma = g))
m1 <- matrix(c(0.2, 0.8, 0.4, 0.6), 2, 2, dimnames = list(c("t1", "t2"), c("A", "null")))
m2 <- matrix(c(0.4, 0.6, 0.8, 0.2), 2, 2, dimnames = dimnames(m1))
fits <- list(f(m1, array(1, c(2, 2, 3))), f(m2, array(0, c(2, 2, 5))))
mm <- mean_gamma_across_seeds(fits)
check("mean model averages the seeds", isTRUE(all.equal(mm, (m1 + m2) / 2)))
check("mean model keeps the dimnames", identical(dimnames(mm), dimnames(m1)))
check("mean model refuses mismatched fits",
      inherits(try(mean_gamma_across_seeds(list(fits[[1]], f(m1[1, , drop = FALSE], array(0, c(1, 2, 2))))), silent = TRUE), "try-error"))
pg <- pooled_gamma_draws(fits)
check("pooled draws concatenate along the draw dimension", identical(dim(pg), c(2L, 2L, 8L)))
check("pooled draws keep both seeds' values", all(pg[, , 1:3] == 1) && all(pg[, , 4:8] == 0))
check("control sharing counts tubes per control",
      identical(neg_sharing(c("a", "b", "a", "c", "a")), c(3L, 1L, 3L, 1L, 3L)))
check("control sharing is 1 when every tube has its own", identical(neg_sharing(c("a", "b")), c(1L, 1L)))

# --- visit dose and model frame -----------------------------------------------------------------
cat("subset_model_frame\n")
check("visit dose parsed", identical(visit_dose(c("V3W315", "V0W0", "x")), c(3L, 0L, NA)))
tubes <- data.frame(filename = c("a", "b", "c"), patient = 1, week = c(6, 6, 150), cocktail = c("p1", "neg", "p1"),
                    is_neg = c(FALSE, TRUE, FALSE), acq_day = "d", instrument = "Canto", operator = "X",
                    prior_vacc = 1, visit = c("V1W6", "V1W6", "V3W150"), stringsAsFactors = FALSE)
bcl <- expand.grid(filename = c("a", "b", "c"), parent = c("CD4", "CD8"), subset = nm[c(1, 16)], stringsAsFactors = FALSE)
bcl$count <- 1; bcl$parent_count <- 100
ccfg <- default_compass_cfg(list(parent_populations = list(CD4 = "CD4", CD8 = "CD8"), cytokine_markers = mk))
fr <- subset_model_frame(bcl, tubes, nm[1], ccfg)
check("only kept non-null subsets", all(fr$subset == nm[1]) && nrow(fr) == 6)
check("stim indicator from is_neg", identical(sort(unique(fr$stim[fr$filename == "b"])), 0L) && all(fr$stim[fr$filename != "b"] == 1))
check("visit factor with the reference dose first", levels(fr$visit)[1] == "V1" && "V3" %in% levels(fr$visit))
check("ps and degree", all(fr$ps %in% c(paste0("CD4:", nm[1]), paste0("CD8:", nm[1]))) && all(fr$degree == 4))
check("visit_id groups the stimulated and neg tube of a sample", length(unique(fr$visit_id[fr$filename %in% c("a", "b")])) == 1)

# --- week basis ---------------------------------------------------------------------------------
cat("week_basis\n")
wb <- week_basis(c(0, 6, 15, 150, 300), 3)
check("three named columns", identical(names(wb), c("w1", "w2", "w3")))
check("prediction at the data reproduces the basis", isTRUE(all.equal(as.matrix(week_basis_at(attr(wb, "basis"), c(0, 6, 15, 150, 300))), as.matrix(wb), check.attributes = FALSE)))

# --- hypothesis summaries ------------------------------------------------------------------------
cat("hypothesis_draws\n")
ps_ <- rbind(c(0.010, 0.002, 0.001), c(0.010, 0.001, 0.001))
pn_ <- rbind(c(0.002, 0.001, 0.002), c(0.012, 0.001, 0.001))
h <- hypothesis_draws(ps_, pn_, c(1L, 2L, 3L))
check("net response is the sum of differences", isTRUE(all.equal(h$response, c(0.008, -0.002))))
check("induced counts only increases", isTRUE(all.equal(h$induced, c(0.009, 0))))
check("polyfunctional share of the induced cells", isTRUE(all.equal(h$share_poly[1], 0.001 / 0.009)) && is.na(h$share_poly[2]))
check("subsets up and down", identical(h$n_up, c(2, 0)) && identical(h$n_down, c(1, 1)))
check("posterior summary", isTRUE(all.equal(post_summary(c(-1, 1, 2, 3))$p_gt0, 0.75)))

# --- simulator ----------------------------------------------------------------------------------------
cat("simulate_subset_counts\n")
d <- fr; d$parent_count <- 1e5
base <- setNames(rep(qlogis(0.001), 2), unique(d$ps))
s0 <- simulate_subset_counts(d, base, function(ps, v) 0, sd_visit = 0, sd_pool = 0, sd_obs = 0, seed = 1)
s1 <- simulate_subset_counts(d, base, function(ps, v) ifelse(grepl("^CD4", ps), 2, 0), sd_visit = 0, sd_pool = 0, sd_obs = 0, seed = 1)
check("null simulation near the baseline", all(abs(s0$count / 1e5 - 0.001) < 0.0005))
check("stimulation effect raises CD4 stimulated rows only",
      all(s1$count[s1$stim == 1 & s1$parent == "CD4"] > 3 * s0$count[s0$stim == 1 & s0$parent == "CD4"]) &&
      all(abs(s1$count[s1$parent == "CD8"] / 1e5 - 0.001) < 0.0005))

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
