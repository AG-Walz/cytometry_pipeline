#!/usr/bin/env Rscript
# =============================================================================
# analysis/differential_expression.R — Per-Marker Differential Expression
# =============================================================================
# Per-marker differential expression between conditions using mixed-effects
# models, limma pseudobulk, or Wilcoxon rank-sum tests.
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Run per-marker differential expression between conditions
#'
#' Input from: clusters_to_gating_csv() (clustering/leiden.R),
#'             load_csv_gating() (loading/load_csv_gating.R)
#' Output to:  plot_volcano() (analysis/differential_abundance.R),
#'             plot_heatmap() (visualization/heatmaps.R)
#'
#' @param input_dir Directory with per-cluster CSVs (expression columns)
#' @param conditions_csv Path to CSV mapping sample_id to condition (columns: sample_id, condition)
#' @param output_dir Output directory
#' @param cluster_regex Regex for cluster column detection
#' @param markers Marker columns to test (NULL = auto-detect)
#' @param test_type "lmer" (mixed-effects), "limma" (pseudobulk), "wilcox" (Wilcoxon)
#' @param p_adjust_method P-value adjustment method
#' @param seed Random seed
#' @return data.frame(cluster, marker, log2fc, pvalue, padj, test)
run_differential_expression <- function(input_dir,
                                         conditions_csv,
                                         output_dir,
                                         cluster_regex = "Clustering annotated",
                                         markers = NULL,
                                         test_type = "lmer",
                                         p_adjust_method = "BH",
                                         seed = 123) {
  require(data.table)
  set.seed(seed)
  make_outdir(output_dir)

  # Read conditions mapping
  conditions <- fread(conditions_csv)
  stopifnot(all(c("sample_id", "condition") %in% names(conditions)))
  cond_levels <- unique(conditions$condition)
  if (length(cond_levels) != 2) {
    stop("Exactly 2 conditions required, found: ", paste(cond_levels, collapse = ", "))
  }
  cat("Conditions:", paste(cond_levels, collapse = " vs "), "\n")

  # List and load CSV files
  csv_info <- list_csvs(input_dir)
  cat("Found", nrow(csv_info), "CSV files\n")

  all_data <- list()
  for (i in seq_len(nrow(csv_info))) {
    dt <- fread(csv_info$file_path[i])
    dt$sample_id <- csv_info$file_name[i]
    all_data[[i]] <- dt
    if (i %% 10 == 0) cat("  Loaded", i, "/", nrow(csv_info), "files\n")
  }
  combined <- rbindlist(all_data, fill = TRUE)
  cat("Combined data:", nrow(combined), "cells\n")

  # Merge conditions
  combined <- merge(combined, conditions, by = "sample_id", all.x = TRUE)
  combined <- combined[!is.na(condition)]
  cat("Cells with condition:", nrow(combined), "\n")

  # Detect cluster columns
  cluster_cols <- grep(cluster_regex, names(combined), value = TRUE)
  if (length(cluster_cols) == 0) {
    stop("No cluster columns matching '", cluster_regex, "' found")
  }
  cat("Cluster columns:", length(cluster_cols), "\n")

  # Auto-detect marker columns if not provided
  if (is.null(markers)) {
    col_classes <- classify_csv_columns(names(combined))
    markers <- col_classes$markers
    # Remove added columns
    markers <- setdiff(markers, c("sample_id", "condition"))
  }
  cat("Markers to test:", length(markers), "\n")

  # Assign each cell to its cluster (use first matching cluster column = 1)
  combined$cluster_label <- NA_character_
  for (cl in cluster_cols) {
    mask <- !is.na(combined[[cl]]) & combined[[cl]] == 1 & is.na(combined$cluster_label)
    if (any(mask)) {
      combined$cluster_label[mask] <- clean_cluster_name(cl)
    }
  }
  combined <- combined[!is.na(cluster_label)]
  clusters <- unique(combined$cluster_label)
  cat("Clusters to test:", length(clusters), "\n")

  # Run differential expression per cluster x marker
  results <- list()
  for (cl in clusters) {
    cl_data <- combined[cluster_label == cl]
    n_samples <- length(unique(cl_data$sample_id))
    cat("  Cluster:", cl, "(", nrow(cl_data), "cells,", n_samples, "samples) ")

    if (n_samples < 3) {
      cat("- skipped (too few samples)\n")
      next
    }

    for (marker in markers) {
      if (!marker %in% names(cl_data)) next
      vals <- cl_data[[marker]]
      if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) next

      result_row <- tryCatch({
        switch(test_type,
          "wilcox" = .de_wilcox(cl_data, marker, cond_levels),
          "lmer"   = .de_lmer(cl_data, marker, cond_levels),
          "limma"  = .de_limma(cl_data, marker, cond_levels),
          stop("Unknown test_type: ", test_type)
        )
      }, error = function(e) {
        data.frame(
          log2fc = NA, pvalue = NA, test = paste0(test_type, "_error"),
          stringsAsFactors = FALSE
        )
      })

      result_row$cluster <- cl
      result_row$marker <- marker
      results[[length(results) + 1]] <- result_row
    }
    cat("\n")
  }

  if (length(results) == 0) {
    cat("No results produced. Check data and parameters.\n")
    return(data.frame())
  }

  results_df <- do.call(rbind, results)
  results_df$padj <- p.adjust(results_df$pvalue, method = p_adjust_method)
  results_df <- results_df[, c("cluster", "marker", "log2fc", "pvalue", "padj", "test")]

  # Save results
  fwrite(results_df, file.path(output_dir, "differential_expression_results.csv"), sep = ";")
  cat("\nDifferential expression analysis complete.\n")
  cat("Total tests:", nrow(results_df), "\n")
  cat("Significant (padj < 0.05):", sum(results_df$padj < 0.05, na.rm = TRUE), "\n")

  results_df
}

# --- Internal: Wilcoxon pseudobulk test ---
.de_wilcox <- function(cl_data, marker, cond_levels) {
  require(data.table)

  # Pseudobulk: median per sample
  pseudobulk <- cl_data[, .(median_expr = median(get(marker), na.rm = TRUE)),
                          by = .(sample_id, condition)]
  grp1 <- pseudobulk[condition == cond_levels[1]]$median_expr
  grp2 <- pseudobulk[condition == cond_levels[2]]$median_expr

  if (length(grp1) < 2 || length(grp2) < 2) {
    return(data.frame(log2fc = NA, pvalue = NA, test = "wilcox_insufficient",
                      stringsAsFactors = FALSE))
  }

  wt <- wilcox.test(grp1, grp2, exact = FALSE)
  log2fc <- log2((mean(grp1, na.rm = TRUE) + 1e-10) /
                 (mean(grp2, na.rm = TRUE) + 1e-10))

  data.frame(log2fc = log2fc, pvalue = wt$p.value, test = "wilcox",
             stringsAsFactors = FALSE)
}

# --- Internal: Mixed-effects model ---
.de_lmer <- function(cl_data, marker, cond_levels) {
  require(lmerTest)

  model_data <- data.frame(
    value = cl_data[[marker]],
    condition = factor(cl_data$condition, levels = cond_levels),
    sample_id = factor(cl_data$sample_id)
  )

  model <- tryCatch({
    lmerTest::lmer(value ~ condition + (1 | sample_id), data = model_data)
  }, error = function(e) {
    # Fallback: simple linear model without random effects
    return(NULL)
  })

  if (is.null(model)) {
    # Fallback to Wilcoxon
    return(.de_wilcox(cl_data, marker, cond_levels))
  }

  coefs <- summary(model)$coefficients
  cond_row <- grep("condition", rownames(coefs))
  if (length(cond_row) == 0) {
    return(data.frame(log2fc = NA, pvalue = NA, test = "lmer_no_effect",
                      stringsAsFactors = FALSE))
  }

  est <- coefs[cond_row[1], "Estimate"]
  p_val <- coefs[cond_row[1], "Pr(>|t|)"]

  # Compute log2FC from group means
  grp_means <- tapply(model_data$value, model_data$condition, mean, na.rm = TRUE)
  log2fc <- log2((grp_means[cond_levels[1]] + 1e-10) /
                 (grp_means[cond_levels[2]] + 1e-10))

  data.frame(log2fc = unname(log2fc), pvalue = p_val, test = "lmer",
             stringsAsFactors = FALSE)
}

# --- Internal: limma pseudobulk ---
.de_limma <- function(cl_data, marker, cond_levels) {
  require(data.table)

  if (!requireNamespace("limma", quietly = TRUE)) {
    cat("limma not available, falling back to wilcox... ")
    return(.de_wilcox(cl_data, marker, cond_levels))
  }

  # Pseudobulk: median per sample
  pseudobulk <- cl_data[, .(median_expr = median(get(marker), na.rm = TRUE)),
                          by = .(sample_id, condition)]
  if (nrow(pseudobulk) < 4) {
    return(data.frame(log2fc = NA, pvalue = NA, test = "limma_insufficient",
                      stringsAsFactors = FALSE))
  }

  # Design matrix
  design <- model.matrix(~ factor(pseudobulk$condition, levels = cond_levels))
  colnames(design) <- c("Intercept", paste0("condition_", cond_levels[2]))

  # Fit limma
  expr_mat <- matrix(pseudobulk$median_expr, nrow = 1)
  colnames(expr_mat) <- pseudobulk$sample_id
  rownames(expr_mat) <- marker

  fit <- limma::lmFit(expr_mat, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")

  data.frame(log2fc = tt$logFC[1], pvalue = tt$P.Value[1], test = "limma",
             stringsAsFactors = FALSE)
}

# --- Convenience: detect marker columns from a cluster CSV ---
.detect_markers <- function(col_names) {
  col_classes <- classify_csv_columns(col_names)
  col_classes$markers
}

# --- CLI execution ---
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat("Usage: Rscript differential_expression.R <input_dir> <conditions_csv> <output_dir> [test_type]\n")
    quit(status = 1)
  }
  run_differential_expression(
    input_dir = args[1],
    conditions_csv = args[2],
    output_dir = args[3],
    test_type = if (length(args) >= 4) args[4] else "lmer"
  )
}
