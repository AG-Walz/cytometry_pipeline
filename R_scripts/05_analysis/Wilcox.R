source("Load_exprs_csv.R")

# ================================================================
# ----------------------------------------------------------------
# 4. Statistical testing function
# ----------------------------------------------------------------
wilcox_safe <- function(x, g) {
  tryCatch(wilcox.test(x ~ g, exact = FALSE)$p.value,
           warning = function(e) NA, error = function(e) NA)
}

out_list <- pblapply(sel_markers, function(mk) {  
  # --- read only the current marker + all selected gates ---
  cols_needed <- c(mk, sel_gates)
  dt          <- rbindlist(
    lapply(files, fread, select = cols_needed, showProgress = FALSE),
    use.names = TRUE)
  
  expr_vec    <- dt[[mk]]
  res_marker  <- vector("list", length(sel_gates))
  
  for (i in seq_along(sel_gates)) {
    gate       <- sel_gates[i]
    g          <- dt[[gate]]
    
    n_in  <- sum(g == 1L);           if (n_in == 0 || n_in == length(g)) next
    n_out <- length(g) - n_in
    w     <- wilcox.test(expr_vec ~ g, exact = FALSE)
    
    U     <- as.numeric(w$statistic)
    AUC   <- U / (as.numeric(n_in) * as.numeric(n_out))   # avoid int overflow
    
    med_in  <- median(expr_vec[g == 1L], na.rm = TRUE)
    med_out <- median(expr_vec[g == 0L], na.rm = TRUE)
    
    res_marker[[i]] <- data.table(
      gate,
      marker        = mk,
      p_value_raw   = max(w$p.value, .Machine$double.xmin),   # clip >0
      AUC,
      delta_median  = med_in - med_out,
      log2FC        = log2((med_in + 1e-9) / (med_out + 1e-9))
    )
  }
  
  rbindlist(res_marker)
})

results <- rbindlist(out_list, use.names = TRUE, fill = TRUE)

# -------- Benjamini–Hochberg per gate --------
results[, p_adj_BH := p.adjust(p_value_raw, "BH"), by = gate]

# -------- write to CSV --------
readr::write_csv(results, file.path(folder, "gate_stats.csv"))

# ----------------------------------------------------------------
# 5. Save & report
# ----------------------------------------------------------------
outfile <- file.path(folder, "gate_stats.csv")
write_csv(results, outfile)
message("\n✓ Analysis complete.")
message("  Results written to: ", outfile)

# Quick console peek
print(
  results %>% arrange(gate, p_adj_BH) %>% 
    group_by(gate) %>% slice_head(n = 5) %>% 
    ungroup()
)
