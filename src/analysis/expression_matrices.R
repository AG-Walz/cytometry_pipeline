#!/usr/bin/env Rscript
# =============================================================================
# analysis/expression_matrices.R — Create Labelled Expression Matrices
# =============================================================================
# Refactored from: R_scripts/05_analysis/Create_labelled_exprmat.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))

#' Create labelled expression matrices grouped by cluster
#'
#' Input from: split_by_cluster() (analysis/cluster_splitting.R)
#' Output to:  plot_heatmap() (visualization/heatmaps.R),
#'             plot_cluster_heatmap() (visualization/heatmaps.R)
#'
#' @param input_dir Directory with cluster-split CSV files
#' @param output_dir Output directory for expression matrices
#' @param markers Character vector of marker columns to include
#' @param group_by_pattern How to group files (regex for cluster group extraction)
#' @return List of group-level expression matrices
create_expression_matrices <- function(input_dir,
                                        output_dir,
                                        markers = NULL,
                                        group_by_pattern = NULL) {
  require(data.table)
  make_outdir(output_dir)

  csv_info <- list_csvs(input_dir)
  cat("Found", nrow(csv_info), "CSV files\n")

  # Detect markers from first file
  if (is.null(markers)) {
    first_dt <- fread(csv_info$file_path[1], nrows = 1)
    markers <- setdiff(names(first_dt), c("file_origin", "cluster", "index", "Index"))
  }

  # Group files by cluster
  groups <- list()
  for (i in seq_len(nrow(csv_info))) {
    fname <- csv_info$file_name[i]
    # Extract cluster group from filename
    parts <- strsplit(fname, "_")[[1]]
    cluster_name <- parts[length(parts)]
    group <- extract_cluster_group(cluster_name)

    if (!group %in% names(groups)) groups[[group]] <- character(0)
    groups[[group]] <- c(groups[[group]], csv_info$file_path[i])
  }

  cat("Found", length(groups), "cluster groups\n")

  # Build expression matrix per group
  result_matrices <- list()
  for (grp in names(groups)) {
    cat("  Building matrix for:", grp, "...")
    files <- groups[[grp]]

    data_list <- lapply(files, function(f) {
      dt <- fread(f, select = intersect(markers, names(fread(f, nrows = 0))))
      dt
    })
    combined <- rbindlist(data_list, fill = TRUE)
    cat(" ", nrow(combined), "cells x", ncol(combined), "markers\n")

    # Save
    out_file <- file.path(output_dir, paste0("exprmat_", grp, ".csv"))
    fwrite(combined, out_file, sep = ";")
    result_matrices[[grp]] <- combined

    # Calculate summary statistics
    stats <- calculate_statistics(combined, cluster_col = NULL,
                                  marker_cols = names(combined))
  }

  # Save summary
  summary_df <- data.frame(
    group = names(groups),
    n_files = sapply(groups, length),
    n_cells = sapply(result_matrices, nrow)
  )
  fwrite(summary_df, file.path(output_dir, "group_matrices_summary.csv"), sep = ";")

  cat("\nExpression matrices created for", length(groups), "groups\n")
  result_matrices
}

#' Combine expression data.tables into a single matrix
#'
#' Takes a named list of data.table objects (one per sample/file), subsets to
#' the intersection of requested markers, row-binds them, and optionally
#' subsamples rows for memory-efficient downstream operations.
#'
#' Input from: load_csv_expression() (loading/load_csv_expression.R)
#' Output to:  plot_heatmap() (visualization/heatmaps.R),
#'             run_pca() (R_shared/stats_utils.R)
#'
#' @param expr_data_list Named list of data.table objects with expression data
#' @param markers Character vector of marker column names to include
#' @param max_cells Maximum cells to retain (random subsample). NULL = no limit.
#' @param seed Random seed for subsampling (default 42)
#' @return Numeric matrix (cells x markers), or NULL if no markers matched
combine_expression_matrices <- function(expr_data_list, markers,
                                         max_cells = NULL, seed = 42) {
  require(data.table)

  combined_list <- lapply(names(expr_data_list), function(nm) {
    dt <- expr_data_list[[nm]]
    marker_cols <- intersect(markers, names(dt))
    if (length(marker_cols) > 0) {
      as.matrix(dt[, ..marker_cols])
    } else {
      NULL
    }
  })
  combined_list <- Filter(Negate(is.null), combined_list)

  if (length(combined_list) == 0) {
    cat("No matching markers found in any file.\n")
    return(NULL)
  }

  combined_matrix <- do.call(rbind, combined_list)
  cat("Combined matrix:", nrow(combined_matrix), "cells x",
      ncol(combined_matrix), "markers from", length(combined_list), "files\n")

  if (!is.null(max_cells) && nrow(combined_matrix) > max_cells) {
    set.seed(seed)
    combined_matrix <- combined_matrix[sample(nrow(combined_matrix), max_cells), ]
    cat("Subsampled to", max_cells, "cells\n")
  }

  combined_matrix
}
