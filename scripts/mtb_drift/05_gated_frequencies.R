#!/usr/bin/env Rscript
# Step 5, layer B: counts of cytokine-positive gates within CD4 and CD8 from the FlowJo gates,
# background subtraction against the matched neg tube (same patient-week), positivity calls.
#   Rscript scripts/mtb_drift/05_gated_frequencies.R --wsp_dir DIR --fcs_dir DIR --out DIR
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "05_gated_frequencies.R --wsp_dir DIR --fcs_dir DIR --out DIR [--config FILE]"
a <- parse_args(list(wsp_dir = "", fcs_dir = "", out = "",
                     config = file.path(REPO, "config", "mtb_drift.yaml")), usage = usage)
cfg <- read_drift_config(a$config)
need_packages(c("CytoML", "flowWorkspace"))

tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
excl <- if (file.exists(file.path(a$out, "excluded_tubes.csv")))
  read.csv(file.path(a$out, "excluded_tubes.csv"))$filename else character(0)
wsps <- list.files(a$wsp_dir, "\\.wsp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
sigs <- read.csv(file.path(a$out, "panel_signatures.csv"), stringsAsFactors = FALSE)

counts <- list()
for (w in wsps) {
  cat("\nWorkspace", basename(w), "\n")
  gsl <- tryCatch(load_workspace_gatingsets(w, a$fcs_dir, signatures = sigs), error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n"); list() })
  for (gs in gsl) counts[[length(counts) + 1L]] <- gatingset_counts(gs)
}
ct <- do.call(rbind, counts)
ct <- ct[!ct$filename %in% excl, ]
write_table(ct, file.path(a$out, "all_gate_counts.csv"))

# cytokine gates directly under a CD4 or CD8 parent
leaf <- function(p) sub("^.*/", "", p)
ct$parent_leaf <- leaf(ct$Parent); ct$gate_leaf <- leaf(ct$Population)
rows <- list()
for (pp in names(cfg$parent_populations)) {
  sel <- ct[grepl(paste0("^", cfg$parent_populations[[pp]], "$"), ct$parent_leaf, ignore.case = TRUE) &
            grepl(cfg$cytokine_gate_regex, ct$gate_leaf, ignore.case = TRUE), ]
  if (nrow(sel)) rows[[pp]] <- data.frame(filename = sel$filename, parent = pp, gate = sel$gate_leaf,
                                          count = sel$Count, parent_count = sel$ParentCount,
                                          stringsAsFactors = FALSE)
}
fq <- do.call(rbind, rows)
# FlowJo gate names differ between workspaces (IFNy, INFy, TNF ...); map them onto the canonical
# cytokine names of the config so one response spans all patients
canon_gate <- function(g) {
  out <- g
  for (cy in cfg$cytokine_markers) for (pat in cfg$markers[[cy]]) {
    if (grepl("^\\^", pat)) next   # channel-name patterns do not apply to gate names
    out[grepl(pat, g, ignore.case = TRUE, perl = TRUE) & out == g] <- cy
  }
  out
}
if (!is.null(fq)) { fq$gate_flowjo <- fq$gate; fq$gate <- canon_gate(fq$gate) }
if (is.null(fq) || !nrow(fq)) stop("no cytokine gates found under ", paste(cfg$parent_populations, collapse = "/"),
                                   "; inspect all_gate_counts.csv and adjust config parent_populations / cytokine_gate_regex")
fq <- merge(fq, tt[, c("filename", "patient", "prior_vacc", "week", "cocktail", "is_neg", "acq_day",
                       "instrument", "patient_day")], by = "filename")
fq$pct <- 100 * fq$count / fq$parent_count

# matched neg per patient-week-parent-gate
neg <- fq[fq$is_neg, c("patient", "week", "parent", "gate", "count", "parent_count", "pct")]
names(neg)[5:7] <- c("neg_count", "neg_parent_count", "neg_pct")
fq <- merge(fq, neg, by = c("patient", "week", "parent", "gate"), all.x = TRUE)
pc <- positivity_call(fq$pct, fq$neg_pct, cfg$positivity$min_subtracted_pct, cfg$positivity$min_fold_over_neg)
fq <- cbind(fq, pc)
fq$positive[fq$is_neg] <- NA
write_table(fq, file.path(a$out, "gated_frequencies.csv"))
cat("\nGates found:\n"); print(table(fq$parent, fq$gate))
cat("\nMissing matched neg:", sum(is.na(fq$neg_pct) & !fq$is_neg), "cocktail rows\n")
