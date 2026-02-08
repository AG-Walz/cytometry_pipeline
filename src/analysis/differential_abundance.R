#!/usr/bin/env Rscript
# =============================================================================
# analysis/differential_abundance.R — GLMM-based Differential Abundance
# =============================================================================
# Refactored from: R_scripts/05_analysis/Differential_Abundance_GLMM.R (1,326 lines)
# Keeps: GLMM fitting, overdispersion checks, results assembly
# Removed: duplicated I/O, memory mgmt, GUI, visualization (-> R_shared)
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Run differential abundance analysis with negative binomial GLMMs
#'
#' @param data_dir Directory with per-cell CSV files
#' @param conditions_csv Path to CSV mapping file_name -> condition
#' @param output_dir Output directory for results
#' @param cluster_regex Regex to identify cluster columns
#' @param condition_levels Character vector of 2 condition labels
#' @param markers_of_interest Character vector (NULL = interactive)
#' @param min_cluster_size Minimum cells per cluster
#' @param p_adjust_method Multiple testing correction method
#' @param max_cells_per_file Max cells to load per file
#' @param seed Random seed
#' @return Data frame of GLMM results per cluster
run_differential_abundance <- function(data_dir,
                                        conditions_csv,
                                        output_dir,
                                        cluster_regex = "Clustering annotated",
                                        condition_levels = c("pCR", "npCR"),
                                        markers_of_interest = NULL,
                                        min_cluster_size = 10,
                                        p_adjust_method = "BH",
                                        max_cells_per_file = 100000,
                                        seed = 123) {
  require(data.table)
  require(lme4)
  set.seed(seed)
  make_outdir(output_dir)

  # Load conditions mapping
  conditions <- fread(conditions_csv)
  stopifnot(all(c("file_name", "condition") %in% names(conditions)))
  cat("Conditions:", paste(unique(conditions$condition), collapse = ", "), "\n")

  # List and load CSV files
  csv_info <- list_csvs(data_dir)
  cat("Found", nrow(csv_info), "CSV files\n")

  # Build abundance matrix: rows = samples, cols = clusters
  abundance_list <- list()

  for (i in seq_len(nrow(csv_info))) {
    file_path <- csv_info$file_path[i]
    file_name <- csv_info$file_name[i]

    dt <- fread(file_path, nrows = max_cells_per_file)
    cluster_cols <- grep(cluster_regex, names(dt), value = TRUE)

    if (length(cluster_cols) == 0) next

    # Count cells per cluster
    counts <- sapply(cluster_cols, function(col) sum(dt[[col]] == 1, na.rm = TRUE))
    names(counts) <- sapply(names(counts), clean_cluster_name)

    abundance_list[[file_name]] <- c(
      total_cells = nrow(dt),
      counts
    )

    if (i %% 10 == 0) check_memory()
  }

  # Convert to data frame
  all_clusters <- unique(unlist(lapply(abundance_list, function(x) names(x)[-1])))
  abundance_df <- data.frame(
    file_name = names(abundance_list),
    total_cells = sapply(abundance_list, function(x) x["total_cells"])
  )
  for (cl in all_clusters) {
    abundance_df[[cl]] <- sapply(abundance_list, function(x) {
      if (cl %in% names(x)) x[[cl]] else 0
    })
  }

  # Merge with conditions
  abundance_df <- merge(abundance_df, conditions, by = "file_name")
  abundance_df$condition <- factor(abundance_df$condition, levels = condition_levels)
  abundance_df$log_total <- log(abundance_df$total_cells)

  cat("Abundance matrix:", nrow(abundance_df), "samples x",
      length(all_clusters), "clusters\n")

  # Filter clusters with too few cells
  cluster_totals <- colSums(abundance_df[, all_clusters, drop = FALSE])
  valid_clusters <- names(cluster_totals[cluster_totals >= min_cluster_size])
  cat("Clusters passing filter:", length(valid_clusters), "/",
      length(all_clusters), "\n")

  # Run GLMM for each cluster
  results <- list()
  for (cl in valid_clusters) {
    cat("  Testing:", cl, "...")
    count_data <- data.frame(
      count = abundance_df[[cl]],
      condition = abundance_df$condition,
      sample_id = abundance_df$file_name,
      log_total = abundance_df$log_total
    )

    model <- fit_nb_glmm(
      count_data,
      formula = count ~ condition + offset(log_total) + (1 | sample_id)
    )

    if (!is.null(model)) {
      # Extract results
      coefs <- summary(model)$coefficients
      cond_row <- grep("condition", rownames(coefs))
      if (length(cond_row) > 0) {
        est <- coefs[cond_row, "Estimate"]
        se <- coefs[cond_row, "Std. Error"]
        z <- coefs[cond_row, "z value"]
        p <- coefs[cond_row, "Pr(>|z|)"]
        results[[cl]] <- data.frame(
          cluster = cl, estimate = est, std_error = se,
          z_value = z, p_value = p,
          log2fc = est / log(2),
          odds_ratio = exp(est),
          method = class(model)[1]
        )
        cat(" p =", round(p, 4), "\n")
      } else {
        cat(" no condition effect\n")
      }
    } else {
      cat(" model failed\n")
      results[[cl]] <- data.frame(
        cluster = cl, estimate = NA, std_error = NA,
        z_value = NA, p_value = NA, log2fc = NA,
        odds_ratio = NA, method = "failed"
      )
    }
  }

  # Combine results
  results_df <- do.call(rbind, results)
  results_df$p_adjusted <- p.adjust(results_df$p_value, method = p_adjust_method)
  results_df$significant <- results_df$p_adjusted < 0.05

  # Save results
  fwrite(results_df, file.path(output_dir, "glmm_results.csv"), sep = ";")

  # Create volcano plot
  tryCatch({
    p <- ggplot2::ggplot(results_df,
                         ggplot2::aes(x = log2fc, y = -log10(p_adjusted),
                                      color = significant)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      ggplot2::scale_color_manual(values = c("grey60", "red")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Differential Abundance (GLMM)",
                    x = "log2 Fold Change", y = "-log10(adj. p-value)")
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      sig_results <- results_df[results_df$significant == TRUE, ]
      if (nrow(sig_results) > 0) {
        p <- p + ggrepel::geom_text_repel(
          data = sig_results,
          ggplot2::aes(label = cluster), size = 3, max.overlaps = 20
        )
      }
    }
    save_plot_formats(p, file.path(output_dir, "volcano_plot"),
                      width = 14, height = 10)
    cat("Volcano plot saved.\n")
  }, error = function(e) warning("Volcano plot failed: ", e$message))

  cat("\nDifferential abundance analysis complete.\n")
  cat("Significant clusters:", sum(results_df$significant, na.rm = TRUE), "\n")

  results_df
}

# --- CLI execution ---
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    cat("Usage: Rscript differential_abundance.R <data_dir> <conditions_csv> <output_dir>\n")
    quit(status = 1)
  }
  run_differential_abundance(
    data_dir = args[1], conditions_csv = args[2], output_dir = args[3]
  )
}
