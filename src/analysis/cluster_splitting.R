#!/usr/bin/env Rscript
# =============================================================================
# analysis/cluster_splitting.R — Split Cells by Cluster
# =============================================================================
# Refactored from: R_scripts/05_analysis/01_split_clusters.R (1,479 lines)
# Keeps: cluster discovery, file splitting, tracking matrix
# Removed: duplicated utils (-> R_shared), GUI (-> gui_utils), stats (-> stats_utils)
# =============================================================================

# ---- cytometry src bootstrap (Nextflow-safe); see src/R_shared/paths.R ----
if (!exists(".cyto_paths_loaded", inherits = TRUE)) {
  local({
    root <- Sys.getenv("CYTOMETRY_SRC", "")
    if (!nzchar(root) || !dir.exists(file.path(root, "R_shared"))) {
      of <- NULL
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        f <- frames[[i]]$ofile
        if (!is.null(f) && nzchar(f)) { of <- f; break }
      }
      d <- if (!is.null(of)) dirname(normalizePath(of)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))

#' Split cells by cluster assignment into separate CSV files
#'
#' Input from: clusters_to_gating_csv() (clustering/leiden.R),
#'             load_csv_gating() (loading/load_csv_gating.R)
#' Output to:  tracking_matrix_to_long() (this file),
#'             run_cluster_enrichment() (analysis/cluster_enrichment.R),
#'             run_differential_abundance() (analysis/differential_abundance.R),
#'             create_expression_matrices() (analysis/expression_matrices.R)
#'
#' @param input_dir Directory with per-cell CSV files
#' @param output_dir Directory for cluster-split output files
#' @param cluster_regex Regex to identify cluster columns
#' @param markers_of_interest Character vector (NULL = interactive selection)
#' @param min_cells_per_cluster Minimum cells to retain a cluster
#' @param chunk_size Processing chunk size
#' @param overwrite_existing Overwrite existing output files
#' @return List with tracking_matrix, cluster_summary, output_files
split_by_cluster <- function(input_dir,
                              output_dir,
                              cluster_regex = "Clustering annotated",
                              markers_of_interest = NULL,
                              min_cells_per_cluster = 10,
                              chunk_size = 10000,
                              overwrite_existing = TRUE) {
  require(data.table)
  require(progress)

  make_outdir(output_dir)
  csv_info <- list_csvs(input_dir)
  cat("Found", nrow(csv_info), "CSV files in", input_dir, "\n")

  # Read first file to detect columns
  first_file <- fread(csv_info$file_path[1], nrows = 0)
  all_cols <- names(first_file)
  cluster_cols <- grep(cluster_regex, all_cols, value = TRUE)
  expr_cols <- setdiff(all_cols, c(cluster_cols, "index", "Index"))

  cat("Detected", length(cluster_cols), "cluster columns\n")
  cat("Detected", length(expr_cols), "expression columns\n")

  # Interactive marker selection if needed
  if (is.null(markers_of_interest)) {
    markers_of_interest <- select_columns_gui(expr_cols, "Select Markers")
    if (is.null(markers_of_interest) || length(markers_of_interest) == 0) {
      markers_of_interest <- expr_cols
    }
  }

  # Initialize tracking matrix
  all_clusters <- character(0)
  tracking_list <- list()
  output_files <- character(0)

  pb <- progress_bar$new(
    format = "Processing [:bar] :percent :current/:total",
    total = nrow(csv_info), clear = FALSE
  )

  for (i in seq_len(nrow(csv_info))) {
    file_path <- csv_info$file_path[i]
    file_name <- csv_info$file_name[i]
    pb$tick()

    # Read CSV
    dt <- fread(file_path)
    file_clusters <- grep(cluster_regex, names(dt), value = TRUE)

    if (length(file_clusters) == 0) {
      cat("  No cluster columns in:", file_name, "\n")
      next
    }

    # Process each cluster column
    for (clust_col in file_clusters) {
      cells_in_cluster <- which(dt[[clust_col]] == 1)
      cluster_name <- clean_cluster_name(clust_col)

      if (length(cells_in_cluster) < min_cells_per_cluster) next

      # Extract expression data for cells in this cluster
      keep_cols <- intersect(markers_of_interest, names(dt))
      cluster_data <- dt[cells_in_cluster, ..keep_cols]
      cluster_data[, file_origin := file_name]
      cluster_data[, cluster := cluster_name]

      # Write output
      out_file <- file.path(output_dir, paste0(file_name, "_", cluster_name, ".csv"))
      if (!file.exists(out_file) || overwrite_existing) {
        fwrite(cluster_data, out_file, sep = ";")
        output_files <- c(output_files, out_file)
      }

      # Track
      if (!cluster_name %in% all_clusters) {
        all_clusters <- c(all_clusters, cluster_name)
      }
      tracking_list[[paste0(file_name, "_", cluster_name)]] <- list(
        file = file_name, cluster = cluster_name,
        n_cells = length(cells_in_cluster)
      )
    }

    check_memory()
  }

  # Build tracking matrix
  tracking_df <- do.call(rbind, lapply(tracking_list, as.data.frame))
  if (nrow(tracking_df) > 0) {
    tracking_matrix <- tidyr::pivot_wider(
      tracking_df,
      id_cols = file, names_from = cluster, values_from = n_cells,
      values_fill = 0
    )
    fwrite(tracking_matrix,
           file.path(output_dir, "cluster_tracking_matrix.csv"), sep = ";")
  } else {
    tracking_matrix <- data.frame()
  }

  # Summary
  cluster_summary <- data.frame(
    cluster = names(table(tracking_df$cluster)),
    total_cells = as.integer(tapply(tracking_df$n_cells, tracking_df$cluster, sum)),
    n_files = as.integer(tapply(tracking_df$n_cells, tracking_df$cluster, length))
  )

  cat("\nSplit complete:", length(output_files), "files written\n")
  cat("Clusters found:", length(all_clusters), "\n")

  list(
    tracking_matrix = tracking_matrix,
    cluster_summary = cluster_summary,
    output_files = output_files,
    markers_used = markers_of_interest
  )
}

#' Convert a wide tracking matrix to long format for plotting
#'
#' Pivots a sample-by-cluster tracking matrix from wide format (one column
#' per cluster) to long format (one row per sample-cluster combination).
#'
#' Input from: split_by_cluster() (this file),
#'             build_tracking_matrix() (analysis/abundance_plotting.R)
#' Output to:  plot_composition_barplot() (visualization/barplots.R)
#'
#' @param tracking_matrix Data frame with an ID column and one numeric column per cluster
#' @param id_col Name of the ID column (default "file")
#' @param names_to Name for the cluster column in output (default "cluster")
#' @param values_to Name for the count column in output (default "count")
#' @return Data frame in long format with columns: id_col, names_to, values_to
tracking_matrix_to_long <- function(tracking_matrix, id_col = "file",
                                     names_to = "cluster", values_to = "count") {
  require(tidyr)
  long_data <- tidyr::pivot_longer(
    tracking_matrix,
    cols = -dplyr::all_of(id_col),
    names_to = names_to,
    values_to = values_to
  )
  as.data.frame(long_data)
}
