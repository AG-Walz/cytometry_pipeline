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

#' Build a sample-by-cluster abundance matrix from per-cell CSVs
#'
#' Reads CSV files with binary cluster indicator columns, counts cells per
#' cluster per sample, and returns a data.frame ready for statistical testing
#' or abundance plotting.
#'
#' Input from: clusters_to_gating_csv() (clustering/leiden.R),
#'             load_csv_gating() (loading/load_csv_gating.R)
#' Output to:  run_differential_abundance() (this file),
#'             run_cluster_enrichment() (analysis/cluster_enrichment.R),
#'             plot_abundance() (analysis/abundance_plotting.R)
#'
#' @param data_dir Directory with per-cell CSV files
#' @param conditions_csv Path to CSV mapping file_name -> condition (optional)
#' @param cluster_regex Regex to identify cluster columns
#' @param condition_levels Character vector of condition labels (optional)
#' @param min_cluster_size Minimum total cells to retain a cluster
#' @param max_cells_per_file Max cells to load per file
#' @return List with: abundance_df (data.frame), all_clusters (character),
#'         valid_clusters (character filtered by min_cluster_size)
build_abundance_matrix <- function(data_dir,
                                    conditions_csv = NULL,
                                    cluster_regex = "Clustering annotated",
                                    condition_levels = NULL,
                                    min_cluster_size = 10,
                                    max_cells_per_file = 100000) {
  require(data.table)

  csv_info <- list_csvs(data_dir)
  cat("Found", nrow(csv_info), "CSV files\n")

  abundance_list <- list()
  for (i in seq_len(nrow(csv_info))) {
    file_path <- csv_info$file_path[i]
    file_name <- csv_info$file_name[i]

    dt <- fread(file_path, nrows = max_cells_per_file)
    cluster_cols <- grep(cluster_regex, names(dt), value = TRUE)
    if (length(cluster_cols) == 0) next

    counts <- sapply(cluster_cols, function(col) sum(dt[[col]] == 1, na.rm = TRUE))
    names(counts) <- sapply(names(counts), clean_cluster_name)

    abundance_list[[file_name]] <- c(total_cells = nrow(dt), counts)
    if (i %% 10 == 0) check_memory()
  }

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

  # Remove samples with zero total cells
  if (any(abundance_df$total_cells == 0)) {
    warning("Removing ", sum(abundance_df$total_cells == 0), " samples with zero total cells")
    abundance_df <- abundance_df[abundance_df$total_cells > 0, ]
  }

  # Merge with conditions if provided
  if (!is.null(conditions_csv)) {
    conditions <- fread(conditions_csv)
    stopifnot(all(c("file_name", "condition") %in% names(conditions)))
    abundance_df <- merge(abundance_df, conditions, by = "file_name")
    if (!is.null(condition_levels)) {
      abundance_df$condition <- factor(abundance_df$condition, levels = condition_levels)
    }
    cat("Conditions:", paste(unique(abundance_df$condition), collapse = ", "), "\n")
  }
  abundance_df$log_total <- log(abundance_df$total_cells)

  cat("Abundance matrix:", nrow(abundance_df), "samples x",
      length(all_clusters), "clusters\n")

  # Filter clusters with too few cells
  cluster_totals <- colSums(abundance_df[, all_clusters, drop = FALSE])
  valid_clusters <- names(cluster_totals[cluster_totals >= min_cluster_size])
  cat("Clusters passing filter:", length(valid_clusters), "/",
      length(all_clusters), "\n")

  list(
    abundance_df = abundance_df,
    all_clusters = all_clusters,
    valid_clusters = valid_clusters
  )
}

#' Create a volcano plot from differential analysis results
#'
#' Generic volcano plot that can visualize results from any statistical test
#' producing log2 fold-change and adjusted p-values.
#'
#' Input from: run_differential_abundance() (this file),
#'             run_cluster_enrichment() (analysis/cluster_enrichment.R)
#' Output to:  standalone visualization output
#'
#' @param results_df Data frame with columns: cluster, log2fc (or log2_enrichment),
#'                   p_adjusted, significant
#' @param output_dir Output directory
#' @param fc_col Column name for fold change (default "log2fc")
#' @param p_col Column name for adjusted p-value (default "p_adjusted")
#' @param label_col Column name for point labels (default "cluster")
#' @param title Plot title
#' @param fc_threshold Fold-change threshold lines (NULL = none)
#' @param p_threshold P-value threshold (default 0.05)
#' @return ggplot object (invisibly)
plot_volcano <- function(results_df, output_dir,
                          fc_col = "log2fc", p_col = "p_adjusted",
                          label_col = "cluster",
                          title = "Volcano Plot",
                          fc_threshold = NULL, p_threshold = 0.05) {
  require(ggplot2)
  make_outdir(output_dir)

  results_df$`.fc` <- results_df[[fc_col]]
  results_df$`.p` <- results_df[[p_col]]
  results_df$`.sig` <- !is.na(results_df$`.p`) & results_df$`.p` < p_threshold

  p <- ggplot(results_df, aes(x = `.fc`, y = -log10(`.p`), color = `.sig`)) +
    geom_point(size = 2) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
    scale_color_manual(values = c("grey60", "red"), name = "Significant") +
    theme_minimal() +
    labs(title = title, x = "log2 Fold Change", y = "-log10(adj. p-value)")

  if (!is.null(fc_threshold)) {
    p <- p + geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                         linetype = "dotted", color = "grey40")
  }

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    sig_results <- results_df[results_df$`.sig` == TRUE, ]
    if (nrow(sig_results) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = sig_results,
        aes(label = .data[[label_col]]), size = 3, max.overlaps = 20
      )
    }
  }

  save_plot_formats(p, file.path(output_dir, "volcano_plot"),
                    width = 14, height = 10)
  cat("Volcano plot saved to:", output_dir, "\n")
  invisible(p)
}

#' Run differential abundance analysis with negative binomial GLMMs
#'
#' Orchestrates the full GLMM-based differential abundance workflow:
#' build_abundance_matrix() -> fit per-cluster GLMMs -> plot_volcano().
#'
#' Input from: clusters_to_gating_csv() (clustering/leiden.R),
#'             load_csv_gating() (loading/load_csv_gating.R)
#' Output to:  plot_volcano() (this file),
#'             plot_abundance() (analysis/abundance_plotting.R)
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

  # Build abundance matrix (reusable step)
  abund <- build_abundance_matrix(
    data_dir = data_dir,
    conditions_csv = conditions_csv,
    cluster_regex = cluster_regex,
    condition_levels = condition_levels,
    min_cluster_size = min_cluster_size,
    max_cells_per_file = max_cells_per_file
  )
  abundance_df <- abund$abundance_df
  valid_clusters <- abund$valid_clusters

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
      # Check for overdispersion if performance package is available
      if (requireNamespace("performance", quietly = TRUE)) {
        tryCatch({
          od <- performance::check_overdispersion(model)
          if (od$p_value < 0.05) {
            cat(" [overdispersed, ratio=", round(od$overdispersion_ratio, 2), "]")
          }
        }, error = function(e) NULL)
      }

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

  # Create volcano plot (reusable function)
  tryCatch({
    plot_volcano(results_df, output_dir,
                 title = "Differential Abundance (GLMM)")
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
