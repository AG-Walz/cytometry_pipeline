#!/usr/bin/env Rscript
# =============================================================================
# test_data/mtb_drift_synthetic.R — synthetic step-4/5 outputs with KNOWN effects
# =============================================================================
# Reads a real tubes.csv (step 1) and writes neg_summaries.csv, neg_emd.csv and
# gated_frequencies.csv into the same folder, generated from the layer A / B
# models with an injected instrument shift and day drift. Lets steps 6 to 9 be
# exercised, and their recovered effects compared with the truth in
# synthetic_truth.csv, without CytoML. Base R only.
#   Rscript test_data/mtb_drift_synthetic.R --out DIR [--instr_shift 0.4] [--tau_day 0.15]
# =============================================================================
args <- commandArgs(TRUE)
get <- function(k, d) { i <- match(k, args); if (is.na(i)) d else type.convert(args[i + 1], as.is = TRUE) }
out <- get("--out", NA); if (is.na(out)) stop("--out required")
instr_shift <- get("--instr_shift", 0.4); tau_day <- get("--tau_day", 0.15); set.seed(get("--seed", 1))

tt <- read.csv(file.path(out, "tubes.csv"), stringsAsFactors = FALSE)
days <- sort(unique(tt$acq_day)); pats <- sort(unique(tt$patient))
c_d <- setNames(rnorm(length(days), 0, tau_day), days)
markers <- c("FSC_A", "SSC_A", "LiveDead", "CD4", "CD8")
summ <- list(); emd <- list()
for (mk in markers) {
  a_p <- setNames(rnorm(length(pats), 0, 0.3), pats)
  base <- c(FSC_A = 100000, SSC_A = 60000, LiveDead = 1.5, CD4 = 3.5, CD8 = 3.2)[mk]
  sc <- if (mk %in% c("FSC_A", "SSC_A")) 8000 else 1
  for (i in seq_len(nrow(tt))) {
    shift <- a_p[as.character(tt$patient[i])] + c_d[tt$acq_day[i]] +
             instr_shift * (tt$instrument[i] != "FACSCantoII") + rnorm(1, 0, 0.08)
    med <- base + sc * shift
    summ[[length(summ) + 1L]] <- data.frame(filename = tt$filename[i], marker = mk, n = 20000,
      median = med, mad = sc * (0.6 + rnorm(1, 0, 0.02)), q10 = med - sc * 0.8, q90 = med + sc * 0.8,
      spread = sc * (1.6 + 0.2 * (tt$instrument[i] != "FACSCantoII") + rnorm(1, 0, 0.05)),
      valley = if (mk %in% c("CD4", "CD8")) med - sc * 1.2 else NA,
      frac_pos = if (mk %in% c("CD4", "CD8")) plogis(qlogis(0.4) + shift * 0.3) else NA)
    emd[[length(emd) + 1L]] <- data.frame(filename = tt$filename[i], marker = mk,
      ref_kind = "pooled_reference_instrument", emd = sc * abs(shift) + sc * abs(rnorm(1, 0, 0.05)))
  }
}
write.csv(do.call(rbind, summ), file.path(out, "neg_summaries.csv"), row.names = FALSE)
# per-neg-tube cells (arcsinh scale) for BayesFlow variant 2: mixture per marker with the same shifts
dir.create(file.path(out, "neg_cells"), showWarnings = FALSE)
sm <- do.call(rbind, summ)
for (fn in tt$filename[tt$is_neg]) {
  n <- 2000; m <- matrix(NA_real_, n, length(markers), dimnames = list(NULL, markers))
  for (mk in markers) {
    med <- sm$median[sm$filename == fn & sm$marker == mk]
    if (mk %in% c("FSC_A", "SSC_A")) m[, mk] <- rnorm(n, med, 8000)
    else if (mk %in% c("CD4", "CD8")) { pos <- runif(n) < 0.45; m[, mk] <- ifelse(pos, rnorm(n, med, 0.35), rnorm(n, med - 2.5, 0.4)) }
    else m[, mk] <- rnorm(n, med, 0.4)
  }
  con <- gzfile(file.path(out, "neg_cells", paste0(fn, ".csv.gz")), "w"); write.csv(round(m, 4), con, row.names = FALSE); close(con)
}
write.csv(do.call(rbind, emd), file.path(out, "neg_emd.csv"), row.names = FALSE)

# layer B: counts per tube x parent x gate on the logit scale, matched negs
gates <- c("IFNg+", "IL2+", "TNFa+", "CD107a+"); rows <- list()
b_instr_B <- instr_shift; a_pB <- setNames(rnorm(length(pats), 0, 0.5), pats)
for (pp in c("CD4", "CD8")) for (g in gates) {
  for (i in seq_len(nrow(tt))) {
    n_parent <- round(runif(1, 20000, 60000))
    eta <- -7 + a_pB[as.character(tt$patient[i])] + c_d[tt$acq_day[i]] +
           b_instr_B * (tt$instrument[i] != "FACSCantoII") + 0.3 * (tt$prior_vacc[i] == 2) +
           (if (tt$is_neg[i]) 0 else 1.2 + 0.02 * tt$week[i]) + rnorm(1, 0, 0.2)
    k <- rbinom(1, n_parent, plogis(eta))
    rows[[length(rows) + 1L]] <- data.frame(filename = tt$filename[i], parent = pp, gate = g, count = k,
      parent_count = n_parent, patient = tt$patient[i], prior_vacc = tt$prior_vacc[i], week = tt$week[i],
      cocktail = tt$cocktail[i], is_neg = tt$is_neg[i], acq_day = tt$acq_day[i], instrument = tt$instrument[i],
      patient_day = tt$patient_day[i], pct = 100 * k / n_parent)
  }
}
fq <- do.call(rbind, rows)
neg <- fq[fq$is_neg, c("patient", "week", "parent", "gate", "count", "parent_count", "pct")]
names(neg)[5:7] <- c("neg_count", "neg_parent_count", "neg_pct")
fq <- merge(fq, neg, by = c("patient", "week", "parent", "gate"), all.x = TRUE)
fq$subtracted_pct <- fq$pct - fq$neg_pct; fq$fold_over_neg <- ifelse(fq$neg_pct > 0, fq$pct / fq$neg_pct, Inf)
fq$positive <- fq$subtracted_pct >= 0.05 & fq$fold_over_neg >= 3; fq$positive[fq$is_neg] <- NA
write.csv(fq, file.path(out, "gated_frequencies.csv"), row.names = FALSE)
write.csv(data.frame(term = c("instrument", rep("day", length(days))), level = c("FACSymphonyA1", days),
                     truth = c(instr_shift, c_d)), file.path(out, "synthetic_truth.csv"), row.names = FALSE)
cat("synthetic layer A/B tables written to", out, "\n")
