#!/usr/bin/env Rscript
# test_data/drift_synthetic.R — synthetic step-4/5 outputs with KNOWN effects
# Reads a real tubes.csv (step 1) and writes neg_summaries.csv, neg_emd.csv and
# gated_frequencies.csv into the same folder, generated from the layer A / BG / B
# models with an injected instrument shift, day drift and background shifts.
# Lets steps 6 to 9 be exercised, and their recovered effects compared with the
# truth in synthetic_truth.csv, without CytoML. Base R only.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "scripts", "drift", "common.R"))
usage <- "drift_synthetic.R --out DIR --config FILE [--instr_shift 0.4] [--tau_day 0.15] [--bg_shift 0.8] [--seed 1]"
a <- parse_args(list(out = "", config = "", instr_shift = 0.4, tau_day = 0.15, bg_shift = 0.8, seed = 1), usage = usage)
cfg <- read_drift_config(a$config)
out <- a$out; ref <- cfg$reference_instrument
instr_shift <- a$instr_shift; tau_day <- a$tau_day; bg_shift <- a$bg_shift; set.seed(a$seed)

tt <- read.csv(file.path(out, "tubes.csv"), stringsAsFactors = FALSE)
other <- setdiff(unique(tt$instrument), ref)[1]   # the instrument the shift is injected on
days <- sort(unique(tt$acq_day)); pats <- sort(unique(tt$patient))
c_d <- setNames(rnorm(length(days), 0, tau_day), days)
markers <- c("FSC_A", "SSC_A", "LiveDead", "CD4", "CD8")
cyto <- c("IFNg", "IL2", "TNFa", "CD107a")   # layer BG intensity: day drift, TNFa instrument shift, prior_vacc 0.3
summ <- list(); emd <- list()
for (mk in c(markers, cyto)) {
  a_p <- setNames(rnorm(length(pats), 0, 0.3), pats)
  base <- c(FSC_A = 100000, SSC_A = 60000, LiveDead = 1.5, CD4 = 3.5, CD8 = 3.2, IFNg = 0.8, IL2 = 0.9, TNFa = 1.0, CD107a = 1.1)[mk]
  sc <- if (mk %in% c("FSC_A", "SSC_A")) 8000 else 1
  for (i in seq_len(nrow(tt))) {
    shift <- a_p[as.character(tt$patient[i])] + c_d[tt$acq_day[i]] + rnorm(1, 0, 0.08) +
             (if (mk %in% cyto) (mk == "TNFa") * instr_shift * (tt$instrument[i] != ref) + 0.3 * (tt$prior_vacc[i] == 2)
              else instr_shift * (tt$instrument[i] != ref))
    med <- base + sc * shift
    summ[[length(summ) + 1L]] <- data.frame(filename = tt$filename[i], marker = mk, n = 20000,
      median = med, mad = sc * (0.6 + rnorm(1, 0, 0.02)), q10 = med - sc * 0.8, q90 = med + sc * 0.8,
      spread = sc * (1.6 + 0.2 * (tt$instrument[i] != ref) + rnorm(1, 0, 0.05)),
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

# layer B: counts per tube x parent x gate on the logit scale, matched negs. Layer BG: the neg
# counts additionally carry a background shift on the other instrument (CD8:TNFa+) and on one day (CD4:IFNg+)
gates <- c("IFNg+", "IL2+", "TNFa+", "CD107a+"); rows <- list()
b_instr_B <- instr_shift; a_pB <- setNames(rnorm(length(pats), 0, 0.5), pats)
bg_day <- names(which.max(table(tt$acq_day[tt$is_neg]))); bg_day_shift <- 1.0
for (pp in c("CD4", "CD8")) for (g in gates) {
  for (i in seq_len(nrow(tt))) {
    n_parent <- round(runif(1, 20000, 60000))
    bg <- if (!tt$is_neg[i]) 0 else bg_shift * (pp == "CD8" && g == "TNFa+" && tt$instrument[i] != ref) +
                                    bg_day_shift * (pp == "CD4" && g == "IFNg+" && tt$acq_day[i] == bg_day)
    eta <- -7 + a_pB[as.character(tt$patient[i])] + c_d[tt$acq_day[i]] + bg +
           b_instr_B * (tt$instrument[i] != ref) + 0.3 * (tt$prior_vacc[i] == 2) +
           (if (tt$is_neg[i]) 0 else 1.2 + 0.02 * tt$week[i]) + rnorm(1, 0, 0.2)
    k <- rbinom(1, n_parent, plogis(eta))
    rows[[length(rows) + 1L]] <- data.frame(filename = tt$filename[i], parent = pp, gate = g, count = k,
      parent_count = n_parent, patient = tt$patient[i], prior_vacc = tt$prior_vacc[i], week = tt$week[i],
      cocktail = tt$cocktail[i], is_neg = tt$is_neg[i], acq_day = tt$acq_day[i], instrument = tt$instrument[i],
      patient_day = tt$patient_day[i], pct = 100 * k / n_parent)
  }
}
fq <- do.call(rbind, rows)
fq <- match_controls(fq, cfg$controls, ref)
fq <- cbind(fq, positivity_call(fq$pct, fq$neg_pct, cfg$positivity$min_subtracted_pct, cfg$positivity$min_fold_over_neg))
fq$positive[fq$is_neg] <- NA
write.csv(fq, file.path(out, "gated_frequencies.csv"), row.names = FALSE)
# truth per layer; response NA = shared by every response of the layer (the day drift c_d, the instrument shift)
tr <- function(layer, response, term, level, truth) data.frame(layer = layer, response = response, term = term, level = level, truth = truth)
truth <- rbind(
  tr("A", NA, "instrument", other, instr_shift), tr("A", NA, "day", days, c_d),
  tr("BG", NA, "day", days, c_d), tr("BG", "TNFa:median", "instrument", other, instr_shift),
  tr("BG", NA, "instrument", other, b_instr_B),   # neg counts carry the layer B instrument shift ...
  tr("BG", "CD8:TNFa+", "instrument", other, b_instr_B + bg_shift),   # ... plus the background shift
  tr("BG", "CD4:IFNg+", "day", bg_day, c_d[bg_day] + bg_day_shift),
  tr("BG", NA, "prior_vacc", "2", 0.3),
  tr("B", NA, "instrument", other, b_instr_B), tr("B", NA, "day", days, c_d), tr("B", NA, "prior_vacc", "2", 0.3))
write.csv(truth, file.path(out, "synthetic_truth.csv"), row.names = FALSE)
cat("synthetic layer A/BG/B tables written to", out, "\n")
