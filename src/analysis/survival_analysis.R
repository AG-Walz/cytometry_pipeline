#!/usr/bin/env Rscript
# =============================================================================
# analysis/survival_analysis.R — Survival Analysis with Cluster Abundances
# =============================================================================
# Performs Cox proportional hazards regression and Kaplan-Meier survival
# analysis using cluster abundance proportions as predictors. Supports
# median and optimal (maxstat) cutpoints for dichotomizing abundance.
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Run survival analysis using cluster abundances as predictors
#'
#' Input from: build_tracking_matrix() (analysis/abundance_plotting.R),
#'             clinical CSV with sample_id, time, event columns
#' Output to:  Standalone Kaplan-Meier plots + Cox results CSV
#'
#' @param abundance_matrix data.frame with sample rows x cluster columns (proportions)
#' @param clinical_csv Path to CSV with sample_id, time, event columns
#' @param output_dir Output directory
#' @param method "cox" (Cox PH), "km" (Kaplan-Meier per cluster), or "both"
#' @param split_method "median" or "optimal" (maximally selected rank statistics)
#' @param p_adjust_method P-value adjustment
#' @return data.frame(cluster, hr, ci_lower, ci_upper, pvalue, padj)
run_survival_analysis <- function(abundance_matrix,
                                   clinical_csv,
                                   output_dir,
                                   method = "both",
                                   split_method = "median",
                                   p_adjust_method = "BH") {
  require(data.table)

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  require(survival)

  make_outdir(output_dir)

  # Read clinical data
  clinical <- fread(clinical_csv)
  stopifnot(all(c("sample_id", "time", "event") %in% names(clinical)))
  cat("Clinical data:", nrow(clinical), "patients\n")

  # Identify sample ID column in abundance matrix
  if (!is.data.frame(abundance_matrix)) {
    abundance_matrix <- as.data.frame(abundance_matrix)
  }
  id_col <- names(abundance_matrix)[1]
  cluster_cols <- names(abundance_matrix)[-1]

  # Ensure cluster columns are numeric
  for (col in cluster_cols) {
    abundance_matrix[[col]] <- as.numeric(abundance_matrix[[col]])
  }

  cat("Abundance matrix:", nrow(abundance_matrix), "samples x",
      length(cluster_cols), "clusters\n")

  # Merge by sample_id
  abundance_matrix$sample_id <- abundance_matrix[[id_col]]
  merged <- merge(abundance_matrix, clinical, by = "sample_id")
  merged$time <- as.numeric(merged$time)
  merged$event <- as.integer(merged$event)

  n_matched <- nrow(merged)
  cat("Matched samples:", n_matched, "/", nrow(abundance_matrix), "\n")

  if (n_matched < 10) {
    stop("Too few matched samples (", n_matched, "). Check sample_id matching.")
  }

  # Remove rows with missing survival data
  merged <- merged[!is.na(merged$time) & !is.na(merged$event), ]
  cat("Samples with complete survival data:", nrow(merged), "\n")

  # Run analysis per cluster
  results <- list()
  km_plots_dir <- file.path(output_dir, "km_plots")
  if (method %in% c("km", "both")) make_outdir(km_plots_dir)

  for (cl in cluster_cols) {
    cl_values <- merged[[cl]]

    # Skip clusters with zero variance
    if (sd(cl_values, na.rm = TRUE) == 0) {
      cat("  Skipping", cl, "- zero variance\n")
      next
    }

    safe_cl <- create_safe_filename(cl)
    cat("  Cluster:", cl, "...")

    # Determine cutpoint for high/low groups
    cutpoint <- tryCatch({
      .determine_cutpoint(merged, cl, split_method)
    }, error = function(e) {
      median(cl_values, na.rm = TRUE)
    })

    merged$group <- ifelse(cl_values >= cutpoint, "high", "low")
    merged$group <- factor(merged$group, levels = c("low", "high"))

    # Check both groups have events
    grp_events <- tapply(merged$event, merged$group, sum, na.rm = TRUE)
    if (any(is.na(grp_events)) || any(grp_events == 0)) {
      cat(" no events in both groups, skipping\n")
      next
    }

    result_row <- data.frame(
      cluster = cl, cutpoint = cutpoint,
      n_high = sum(merged$group == "high"),
      n_low = sum(merged$group == "low"),
      hr = NA, ci_lower = NA, ci_upper = NA,
      pvalue = NA, cox_concordance = NA,
      km_pvalue = NA,
      stringsAsFactors = FALSE
    )

    # Cox proportional hazards
    if (method %in% c("cox", "both")) {
      cox_result <- tryCatch({
        .run_cox(merged, cl)
      }, error = function(e) {
        cat(" cox failed:", e$message)
        NULL
      })

      if (!is.null(cox_result)) {
        result_row$hr <- cox_result$hr
        result_row$ci_lower <- cox_result$ci_lower
        result_row$ci_upper <- cox_result$ci_upper
        result_row$pvalue <- cox_result$pvalue
        result_row$cox_concordance <- cox_result$concordance
      }
    }

    # Kaplan-Meier
    if (method %in% c("km", "both")) {
      km_result <- tryCatch({
        .run_km(merged, cl, km_plots_dir, safe_cl)
      }, error = function(e) {
        cat(" km failed:", e$message)
        NULL
      })

      if (!is.null(km_result)) {
        result_row$km_pvalue <- km_result$pvalue
        # Use KM p-value if Cox was not run
        if (is.na(result_row$pvalue)) {
          result_row$pvalue <- km_result$pvalue
        }
      }
    }

    results[[cl]] <- result_row
    cat(" HR =", round(result_row$hr, 3),
        " p =", round(result_row$pvalue, 4), "\n")
  }

  if (length(results) == 0) {
    cat("No results produced. Check data and parameters.\n")
    return(data.frame())
  }

  results_df <- do.call(rbind, results)
  results_df$padj <- p.adjust(results_df$pvalue, method = p_adjust_method)

  # Reorder columns
  results_df <- results_df[, c("cluster", "hr", "ci_lower", "ci_upper",
                                "pvalue", "padj", "km_pvalue",
                                "cutpoint", "n_high", "n_low",
                                "cox_concordance")]

  # Save results
  fwrite(results_df, file.path(output_dir, "survival_results.csv"), sep = ";")
  cat("\nSurvival analysis complete.\n")
  cat("Total clusters tested:", nrow(results_df), "\n")
  cat("Significant (padj < 0.05):", sum(results_df$padj < 0.05, na.rm = TRUE), "\n")

  # Summary forest plot
  tryCatch({
    .plot_forest(results_df, output_dir)
  }, error = function(e) {
    cat("Forest plot failed:", e$message, "\n")
  })

  results_df
}

# --- Internal: Determine cutpoint for dichotomization ---
.determine_cutpoint <- function(merged, cluster_col, split_method) {
  if (split_method == "optimal") {
    # Maximally selected rank statistics via survminer
    if (requireNamespace("survminer", quietly = TRUE)) {
      sc <- survminer::surv_cutpoint(
        merged, time = "time", event = "event",
        variables = cluster_col
      )
      return(sc$cutpoint[[cluster_col]]$cutpoint)
    } else {
      cat(" survminer not available, using median... ")
    }
  }

  # Default: median split
  median(merged[[cluster_col]], na.rm = TRUE)
}

# --- Internal: Cox PH regression ---
.run_cox <- function(merged, cluster_col) {
  require(survival)

  formula <- as.formula(paste0("Surv(time, event) ~ `", cluster_col, "`"))
  cox_fit <- coxph(formula, data = merged)
  cox_summary <- summary(cox_fit)

  list(
    hr = cox_summary$conf.int[1, "exp(coef)"],
    ci_lower = cox_summary$conf.int[1, "lower .95"],
    ci_upper = cox_summary$conf.int[1, "upper .95"],
    pvalue = cox_summary$coefficients[1, "Pr(>|z|)"],
    concordance = cox_summary$concordance["C"]
  )
}

# --- Internal: Kaplan-Meier analysis with plot ---
.run_km <- function(merged, cluster_col, plots_dir, safe_cl) {
  require(survival)

  # Log-rank test
  formula_diff <- as.formula("Surv(time, event) ~ group")
  surv_diff <- survdiff(formula_diff, data = merged)
  km_pvalue <- 1 - pchisq(surv_diff$chisq, df = length(surv_diff$n) - 1)

  # KM fit for plotting
  km_fit <- survfit(formula_diff, data = merged)

  # Plot KM curves
  tryCatch({
    if (requireNamespace("survminer", quietly = TRUE)) {
      p <- survminer::ggsurvplot(
        km_fit, data = merged,
        pval = TRUE, conf.int = TRUE,
        risk.table = TRUE,
        xlab = "Time", ylab = "Survival Probability",
        title = paste0("Kaplan-Meier: ", cluster_col),
        legend.labs = c("Low", "High"),
        palette = c("#2166AC", "#B2182B"),
        ggtheme = ggplot2::theme_minimal()
      )
      # survminer returns a list; save the combined plot
      grDevices::pdf(file.path(plots_dir, paste0("km_", safe_cl, ".pdf")),
                     width = 10, height = 8)
      print(p)
      grDevices::dev.off()

      grDevices::png(file.path(plots_dir, paste0("km_", safe_cl, ".png")),
                     width = 10, height = 8, units = "in", res = 200)
      print(p)
      grDevices::dev.off()
    } else {
      # Fallback: base R survival plot
      save_base_plot_formats(
        file.path(plots_dir, paste0("km_", safe_cl)),
        plot_fn = function() {
          plot(km_fit, col = c("#2166AC", "#B2182B"), lwd = 2,
               xlab = "Time", ylab = "Survival Probability",
               main = paste0("Kaplan-Meier: ", cluster_col))
          legend("topright", legend = c("Low", "High"),
                 col = c("#2166AC", "#B2182B"), lwd = 2)
          text(max(merged$time) * 0.7, 0.1,
               paste0("p = ", formatC(km_pvalue, digits = 4, format = "f")))
        },
        width = 10, height = 8
      )
    }
  }, error = function(e) {
    cat(" plot failed:", e$message)
  })

  list(pvalue = km_pvalue, km_fit = km_fit)
}

# --- Internal: Forest plot of hazard ratios ---
.plot_forest <- function(results_df, output_dir) {
  require(ggplot2)

  # Filter to results with valid HRs
  plot_data <- results_df[!is.na(results_df$hr), ]
  if (nrow(plot_data) == 0) return(invisible(NULL))

  # Sort by HR
  plot_data <- plot_data[order(plot_data$hr), ]
  plot_data$cluster <- factor(plot_data$cluster, levels = plot_data$cluster)
  plot_data$significant <- !is.na(plot_data$padj) & plot_data$padj < 0.05

  p <- ggplot(plot_data, aes(x = hr, y = cluster, color = significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("grey60", "red"), name = "Significant\n(padj < 0.05)") +
    scale_x_log10() +
    theme_minimal() +
    labs(title = "Hazard Ratios by Cluster Abundance",
         x = "Hazard Ratio (log scale)", y = NULL) +
    theme(axis.text.y = element_text(size = 8))

  save_plot_formats(p, file.path(output_dir, "forest_plot_hr"),
                    width = 12, height = max(6, nrow(plot_data) * 0.4))
  cat("Forest plot saved to:", output_dir, "\n")
  invisible(p)
}

# --- CLI execution ---
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat("Usage: Rscript survival_analysis.R <abundance_csv> <clinical_csv> <output_dir> [method] [split_method]\n")
    quit(status = 1)
  }

  require(data.table)
  abundance <- as.data.frame(fread(args[1]))
  run_survival_analysis(
    abundance_matrix = abundance,
    clinical_csv = args[2],
    output_dir = args[3],
    method = if (length(args) >= 4) args[4] else "both",
    split_method = if (length(args) >= 5) args[5] else "median"
  )
}
