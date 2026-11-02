#!/usr/bin/env Rscript
# Step 13: responder calls under different rules, for question 2 (docs/plans/drift_three_questions_plan.md).
# Per stimulated tube and cytokine gate:
#   lab     background-subtracted frequency >= min_subtracted_pct and >= min_fold_over_neg times the neg
#           (config `positivity`, the rule already in gated_frequencies.csv)
#   fisher  one-sided Fisher's exact test, stimulated vs matched neg counts, Benjamini-Hochberg across
#           all tubes of the response, q < 0.05 (HVTN practice, Horton et al. 2007)
#   mimosa  MIMOSA posterior probability of response >= 0.5 with FDR < 0.05 (Finak et al. 2014), when the
#           MIMOSA package is installed. Pools are patient-specific, so all pools of a response share
#           one MIMOSA model (stimulation: peptide pool).
# Writes calls_compare.csv (one row per tube x response) and calls_agreement.csv (per response: calls
# per rule, pairwise agreement and Cohen's kappa).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "13_calling_rules.R --out DIR --config FILE"
a <- parse_args(list(out = "", config = ""), usage = usage)
cfg <- read_drift_config(a$config)
fq <- read.csv(file.path(a$out, "gated_frequencies.csv"), stringsAsFactors = FALSE)
st <- fq[!fq$is_neg & !is.na(fq$neg_count) & fq$parent_count > 0 & fq$neg_parent_count > 0, ]
st$response <- paste(st$parent, st$gate, sep = ":")
calls <- calling_rules(st, cfg)
if (mimosa_available()) {
  suppressPackageStartupMessages(library(MIMOSA))   # its formula interface needs the package attached
  mc <- mimosa_calls(st); calls$mimosa <- mc$positive; calls$mimosa_prob <- mc$prob; calls$mimosa_fdr <- mc$fdr
} else cat("MIMOSA not installed: mimosa column left NA\n")
write_table(calls, file.path(a$out, "calls_compare.csv"))
write_table(calls_agreement(calls), file.path(a$out, "calls_agreement.csv"))
print(calls_agreement(calls), digits = 2, row.names = FALSE)
