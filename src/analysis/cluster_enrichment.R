#!/usr/bin/env Rscript
# =============================================================================
# analysis/cluster_enrichment.R — Fisher/Mixed-Effects Cluster Enrichment
# =============================================================================
# Refactored from: R_scripts/05_analysis/Fine_clustering_cytolution.R (1,252 lines)
# Keeps: Fisher/mixed-effects analysis, condition handling
# Removed: duplicated code (-> R_shared)
# Resume functionality built in via --resume flag
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Run cluster enrichment analysis between conditions
#'
#' @param data_dir Directory with per-cell CSV files
#' @param conditions_csv Path to CSV with file_name, condition columns
#' @param output_dir Output directory
#' @param cluster_regex Regex for cluster columns
#' @param condition_levels Character vector of 2 condition labels
#' @param stat_method "pooled" (Fisher only) or "mixed_effects"
#' @param min_cluster_frac Minimum fraction to keep a cluster
#' @param p_adjust_method Correction method
#' @param resume Resume from previous run (loads intermediate results)
#' @param seed Random seed
#' @return Data frame with enrichment results per cluster
run_cluster_enrichment <- function(data_dir,
                                    conditions_csv,
                                    output_dir,
                                    cluster_regex = "Clustering annotated",
                                    condition_levels = c("pCR", "npCR"),
                                    stat_method = "mixed_effects",
                                    min_cluster_frac = 0.001,
                                    p_adjust_method = "BH",
                                    resume = FALSE,
                                    seed = 123) {
  require(data.table)
  set.seed(seed)
  make_outdir(output_dir)

  # Check for resume
  intermediate_file <- file.path(output_dir, "intermediate_results.rds")
  if (resume && file.exists(intermediate_file)) {
    cat("Resuming from intermediate results...\n")
    prev <- readRDS(intermediate_file)
    completed_clusters <- prev$completed_clusters
    results <- prev$results
  } else {
    completed_clusters <- character(0)
    results <- list()
  }

  # Load conditions
  conditions <- fread(conditions_csv)
  stopifnot(all(c("file_name", "condition") %in% names(conditions)))

  # Load and combine cell data
  csv_info <- list_csvs(data_dir)
  cat("Loading", nrow(csv_info), "files...\n")

  cells_list <- list()
  for (i in seq_len(nrow(csv_info))) {
    dt <- fread(csv_info$file_path[i])
    dt$file_name <- csv_info$file_name[i]
    cells_list[[i]] <- dt
    if (i %% 10 == 0) { check_memory(); cat(".") }
  }
  cat("\n")

  cells_data <- rbindlist(cells_list, fill = TRUE)
  rm(cells_list); gc()

  # Merge conditions
  cells_data <- merge(cells_data, conditions, by = "file_name", all.x = TRUE)
  cells_data <- cells_data[!is.na(condition)]
  cat("Total cells with conditions:", nrow(cells_data), "\n")

  # Identify cluster columns
  cluster_cols <- grep(cluster_regex, names(cells_data), value = TRUE)
  cat("Cluster columns:", length(cluster_cols), "\n")

  # Filter by minimum fraction
  total_cells <- nrow(cells_data)
  valid_clusters <- cluster_cols[sapply(cluster_cols, function(col) {
    sum(cells_data[[col]] == 1, na.rm = TRUE) / total_cells >= min_cluster_frac
  })]
  cat("Clusters passing filter:", length(valid_clusters), "\n")

  # Skip already-completed clusters (resume mode)
  clusters_to_run <- setdiff(valid_clusters, completed_clusters)

  # Run analysis
  for (cl in clusters_to_run) {
    cl_clean <- clean_cluster_name(cl)
    cat("  Analyzing:", cl_clean, "...")

    # Fisher's exact test
    fisher_result <- calculate_enrichment(cl, cells_data, condition_levels)

    # Mixed-effects (optional)
    if (stat_method == "mixed_effects") {
      me_result <- do_mixed_effects(cl, cells_data, condition_levels)
      combined <- cbind(
        data.frame(cluster = cl_clean),
        fisher_result,
        me_result %>% dplyr::rename_with(~ paste0("me_", .x))
      )
    } else {
      combined <- cbind(data.frame(cluster = cl_clean), fisher_result)
    }

    results[[cl_clean]] <- combined
    completed_clusters <- c(completed_clusters, cl)

    # Save intermediate results periodically
    if (length(results) %% 10 == 0) {
      saveRDS(list(completed_clusters = completed_clusters, results = results),
              intermediate_file)
    }

    cat(" p =", round(fisher_result$p_value, 4), "\n")
  }

  # Combine results
  results_df <- do.call(rbind, results)
  results_df$p_adjusted <- p.adjust(results_df$p_value, method = p_adjust_method)
  results_df$significant <- results_df$p_adjusted < 0.05

  # Log2 enrichment
  results_df$log2_enrichment <- log2(
    (results_df$in_A / (results_df$in_A + results_df$not_in_A)) /
    (results_df$in_B / (results_df$in_B + results_df$not_in_B))
  )

  # Save final results
  fwrite(results_df, file.path(output_dir, "enrichment_results.csv"), sep = ";")

  # Clean up intermediate
  if (file.exists(intermediate_file)) file.remove(intermediate_file)

  cat("\nEnrichment analysis complete.\n")
  cat("Significant clusters:", sum(results_df$significant, na.rm = TRUE), "\n")

  results_df
}
