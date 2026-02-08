#!/usr/bin/env Rscript
# =============================================================================
# analysis/abundance_plotting.R — Cluster Abundance Plots
# =============================================================================
# Refactored from: pipelines/run_ICK_abundance.R + R_scripts/05_analysis/Plot_Abundance.R
#
# Provides:
#   build_tracking_matrix()       — count cells per cluster per file
#   plot_abundance_stacked()      — barplot sorted by total cell count
#   plot_abundance_proportional() — proportional barplot with dendrogram
#   plot_abundance_totals()       — total cells per cluster
#   plot_abundance_heatmap()      — proportion heatmap (biclustered)
#   plot_abundance()              — orchestrator calling all four
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Build a sample-by-cluster tracking matrix from CSVs
#'
#' @param input_dir Directory containing per-sample CSVs
#' @param cluster_regex Regex to identify cluster columns
#' @param exclude_patterns Character vector of regex patterns to exclude
#' @param output_dir Output directory (tracking matrix CSV saved here)
#' @return data.frame with columns: file, cluster1, cluster2, ...
build_tracking_matrix <- function(input_dir, cluster_regex, exclude_patterns,
                                  output_dir) {
  require(data.table)
  make_outdir(output_dir)

  csv_info <- list_csvs(input_dir)
  cat("Found", nrow(csv_info), "CSV files in", input_dir, "\n")

  # Discover columns from first file
  first_header <- names(fread(csv_info$file_path[1], nrows = 0))
  all_cluster_cols <- grep(cluster_regex, first_header, value = TRUE, perl = TRUE)

  # Exclude sub-cluster patterns to keep parent-level only
  parent_cols <- all_cluster_cols
  for (pat in exclude_patterns) {
    parent_cols <- parent_cols[!grepl(pat, parent_cols)]
  }
  cat("Total clustering columns:", length(all_cluster_cols), "\n")
  cat("Parent-level clusters:  ", length(parent_cols), "\n")
  cat("Clusters to track:\n")
  for (col in parent_cols) cat("  -", clean_cluster_name(col), "\n")

  # Count cells per cluster per file
  n_files <- nrow(csv_info)
  tracking_list <- vector("list", n_files)
  batch_time <- Sys.time()

  n_skipped <- 0
  for (i in seq_len(n_files)) {
    # Check that the file actually contains the expected columns
    file_header <- names(fread(csv_info$file_path[i], nrows = 0, showProgress = FALSE))
    cols_present <- intersect(parent_cols, file_header)
    if (length(cols_present) == 0) {
      n_skipped <- n_skipped + 1
      next
    }

    dt <- fread(csv_info$file_path[i], select = cols_present, showProgress = FALSE)
    counts <- colSums(dt, na.rm = TRUE)
    names(counts) <- sapply(names(counts), clean_cluster_name)
    tracking_list[[i]] <- c(file = csv_info$file_name[i], counts)

    if (i %% 10 == 0 || i == n_files) {
      elapsed_s <- as.numeric(difftime(Sys.time(), batch_time, units = "secs"))
      rate <- i / elapsed_s
      eta <- (n_files - i) / rate
      cat(sprintf("  [%d/%d] %s | %.1f files/sec | ETA: %.0fs\n",
                  i, n_files, csv_info$file_name[i], rate, eta))
      flush(stdout())
    }
  }
  if (n_skipped > 0) cat("  Skipped", n_skipped, "files without cluster columns\n")

  # Remove NULLs (skipped files)
  tracking_list <- Filter(Negate(is.null), tracking_list)
  tracking_matrix <- as.data.frame(do.call(rbind, tracking_list),
                                   stringsAsFactors = FALSE)
  for (col in names(tracking_matrix)[-1]) {
    tracking_matrix[[col]] <- as.integer(tracking_matrix[[col]])
  }

  tracking_path <- file.path(output_dir, "cluster_tracking_matrix.csv")
  fwrite(tracking_matrix, tracking_path, sep = ";")
  cat("\nTracking matrix saved:", tracking_path, "\n")
  cat("  ", nrow(tracking_matrix), "samples x",
      ncol(tracking_matrix) - 1, "clusters\n")

  tracking_matrix
}

#' Stacked barplot sorted by total cell count (descending)
#'
#' @param tracking_matrix data.frame from build_tracking_matrix()
#' @param output_dir Output directory
plot_abundance_stacked <- function(tracking_matrix, output_dir) {
  require(ggplot2)
  require(tidyr)
  make_outdir(output_dir)

  long_df <- tidyr::pivot_longer(tracking_matrix, cols = -file,
                                  names_to = "cluster", values_to = "count")
  long_df <- as.data.frame(long_df)

  sample_totals <- rowSums(tracking_matrix[, -1])
  names(sample_totals) <- tracking_matrix$file
  stacked_order <- names(sort(sample_totals, decreasing = TRUE))
  long_df$file <- factor(long_df$file, levels = stacked_order)

  p <- ggplot(long_df, aes(x = file, y = count, fill = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    theme_minimal() +
    labs(title = "Cluster Abundance per Sample (sorted by total cell count)",
         x = NULL, y = "Cell Count") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          panel.grid.major.x = element_blank()) +
    scale_fill_viridis_d(option = "turbo") +
    scale_y_continuous(labels = scales::comma)

  save_plot_formats(p, file.path(output_dir, "abundance_stacked"),
                    width = 28, height = 10, png_dpi = 150)
  invisible(p)
}

#' Proportional barplot with Ward.D2 dendrogram on top
#'
#' @param tracking_matrix data.frame from build_tracking_matrix()
#' @param output_dir Output directory
plot_abundance_proportional <- function(tracking_matrix, output_dir) {
  require(ggplot2)
  require(tidyr)
  require(ggdendro)
  require(patchwork)
  make_outdir(output_dir)

  # Remove samples with zero total cells
  row_totals <- rowSums(tracking_matrix[, -1])
  tracking_matrix <- tracking_matrix[row_totals > 0, , drop = FALSE]

  long_df <- tidyr::pivot_longer(tracking_matrix, cols = -file,
                                  names_to = "cluster", values_to = "count")
  long_df <- as.data.frame(long_df)

  # Compute dendrogram order from proportions
  mat <- as.matrix(tracking_matrix[, -1])
  rownames(mat) <- tracking_matrix$file
  prop_mat <- mat / rowSums(mat)
  prop_mat[is.nan(prop_mat)] <- 0
  hc_samples <- hclust(dist(prop_mat), method = "ward.D2")
  dendro_order <- hc_samples$labels[hc_samples$order]
  long_df$file <- factor(long_df$file, levels = dendro_order)

  # Build dendrogram plot
  dendro_data <- ggdendro::dendro_data(as.dendrogram(hc_samples), type = "rectangle")
  dendro_segments <- ggdendro::segment(dendro_data)

  p_dendro <- ggplot(dendro_segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.3) +
    scale_x_continuous(expand = c(0.005, 0.005)) +
    theme_void()

  # Build proportional barplot
  p_prop <- ggplot(long_df, aes(x = file, y = count, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    theme_minimal() +
    labs(x = NULL, y = "Proportion") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(0, 5, 5, 5)) +
    scale_fill_viridis_d(option = "turbo") +
    scale_y_continuous(expand = c(0, 0))

  # Combine dendrogram + barplot
  p_combined <- p_dendro / p_prop +
    patchwork::plot_layout(heights = c(1, 5)) +
    patchwork::plot_annotation(
      title = "Cluster Proportions per Sample (hierarchically clustered)"
    )

  save_plot_formats(p_combined, file.path(output_dir, "abundance_proportional"),
                    width = 28, height = 12, png_dpi = 150)
  invisible(p_combined)
}

#' Total cells per cluster across all samples
#'
#' @param tracking_matrix data.frame from build_tracking_matrix()
#' @param output_dir Output directory
plot_abundance_totals <- function(tracking_matrix, output_dir) {
  require(ggplot2)
  require(tidyr)
  make_outdir(output_dir)

  long_df <- tidyr::pivot_longer(tracking_matrix, cols = -file,
                                  names_to = "cluster", values_to = "count")
  long_df <- as.data.frame(long_df)

  cluster_totals <- aggregate(count ~ cluster, data = long_df, sum)
  cluster_totals <- cluster_totals[order(cluster_totals$count, decreasing = TRUE), ]

  p <- ggplot(cluster_totals,
               aes(x = reorder(cluster, -count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Total Cells per Cluster (All Samples)",
         x = NULL, y = "Cell Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    scale_y_continuous(labels = scales::comma)

  save_plot_formats(p, file.path(output_dir, "total_cells_per_cluster"),
                    width = 14, height = 8)
  invisible(p)
}

#' Proportion heatmap (biclustered via Ward.D2)
#'
#' @param tracking_matrix data.frame from build_tracking_matrix()
#' @param output_dir Output directory
plot_abundance_heatmap <- function(tracking_matrix, output_dir) {
  require(ggplot2)
  require(reshape2)
  make_outdir(output_dir)

  # Remove samples with zero total cells
  row_totals <- rowSums(tracking_matrix[, -1])
  tracking_matrix <- tracking_matrix[row_totals > 0, , drop = FALSE]

  mat <- as.matrix(tracking_matrix[, -1])
  rownames(mat) <- tracking_matrix$file
  prop_mat <- mat / rowSums(mat)
  prop_mat[is.nan(prop_mat)] <- 0

  hc_samples  <- hclust(dist(prop_mat), method = "ward.D2")
  hc_clusters <- hclust(dist(t(prop_mat)), method = "ward.D2")
  hm_mat <- prop_mat[hc_samples$order, hc_clusters$order]

  long_hm <- reshape2::melt(hm_mat)
  names(long_hm) <- c("sample", "cluster", "proportion")
  long_hm$sample  <- factor(long_hm$sample,  levels = rownames(hm_mat))
  long_hm$cluster <- factor(long_hm$cluster, levels = colnames(hm_mat))

  p <- ggplot(long_hm, aes(x = cluster, y = sample, fill = proportion)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno", name = "Proportion") +
    theme_minimal() +
    labs(title = "Cluster Proportion per Sample (hierarchically clustered)",
         x = "Cluster", y = "Sample") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 3)
    )

  save_plot_formats(p, file.path(output_dir, "cluster_proportion_heatmap"),
                    width = 16, height = 24, png_dpi = 150)
  invisible(p)
}

#' Generate all abundance plots
#'
#' @param tracking_matrix_path Path to cluster_tracking_matrix.csv (semicolon-separated)
#' @param output_dir Output directory
#' @param tracking_matrix Optional: pre-loaded data.frame (if NULL, reads from path)
#' @param metadata_file Optional: path to metadata CSV for grouped plots
#' @param file_col Column name in metadata containing the file identifiers
#' @param metadata_sep Separator for the metadata CSV
#' @param skip_prefixes Prefix(es) for metadata columns to skip (default "*")
plot_abundance <- function(tracking_matrix_path = NULL, output_dir,
                           tracking_matrix = NULL,
                           metadata_file = NULL,
                           file_col = "*file_name",
                           metadata_sep = ";",
                           skip_prefixes = "\\*") {
  require(data.table)

  if (is.null(tracking_matrix)) {
    stopifnot(!is.null(tracking_matrix_path), file.exists(tracking_matrix_path))
    tracking_matrix <- as.data.frame(fread(tracking_matrix_path))
  }

  cat("Generating abundance plots:", nrow(tracking_matrix), "samples x",
      ncol(tracking_matrix) - 1, "clusters\n")

  cat("  Plotting: stacked barplot\n")
  plot_abundance_stacked(tracking_matrix, output_dir)

  cat("  Plotting: proportional barplot (dendrogram-ordered)\n")
  plot_abundance_proportional(tracking_matrix, output_dir)

  cat("  Plotting: total cells per cluster\n")
  plot_abundance_totals(tracking_matrix, output_dir)

  cat("  Plotting: cluster proportion heatmap\n")
  plot_abundance_heatmap(tracking_matrix, output_dir)

  # Metadata-grouped plots
  if (!is.null(metadata_file) && file.exists(metadata_file)) {
    cat("  Plotting: metadata-grouped abundance\n")
    plot_abundance_by_metadata(tracking_matrix, metadata_file, output_dir,
                               file_col = file_col,
                               metadata_sep = metadata_sep,
                               skip_prefixes = skip_prefixes)
  }

  cat("Abundance plots saved to:", output_dir, "\n")
  invisible(NULL)
}

# =============================================================================
# Metadata-grouped abundance plots
# =============================================================================

#' Generate abundance plots grouped by each metadata column
#'
#' Merges the tracking matrix with a metadata CSV and, for each usable metadata
#' column, produces:
#'   1. Faceted stacked barplot (one panel per group value, samples as bars)
#'   2. Aggregated proportional barplot (one bar per group value)
#'
#' @param tracking_matrix data.frame with file + cluster count columns
#' @param metadata_file Path to metadata CSV (e.g. export_cytolution.csv)
#' @param output_dir Base output directory
#' @param file_col Column name in metadata with file identifiers (default "*file_name")
#' @param metadata_sep CSV separator for metadata file (default ";")
#' @param skip_prefixes Regex prefix(es) for columns to skip (default "\\*")
#' @param max_unique Max unique values for a column to be plotted (default 30)
plot_abundance_by_metadata <- function(tracking_matrix, metadata_file, output_dir,
                                       file_col = "*file_name",
                                       metadata_sep = ";",
                                       skip_prefixes = "\\*",
                                       max_unique = 30) {
  require(data.table)
  require(ggplot2)
  require(tidyr)

  metadata <- as.data.frame(fread(metadata_file, sep = metadata_sep))
  stopifnot(file_col %in% names(metadata))

  # Match file names: strip extensions from metadata, tracking matrix already has no ext
  metadata$`.file_match` <- tools::file_path_sans_ext(metadata[[file_col]])

  # Merge
  merged <- merge(tracking_matrix, metadata, by.x = "file", by.y = ".file_match")
  n_matched <- nrow(merged)
  cat("  Metadata merge:", n_matched, "/", nrow(tracking_matrix), "samples matched\n")
  if (n_matched == 0) {
    warning("No samples matched between tracking matrix and metadata. Check file names.")
    return(invisible(NULL))
  }

  # Identify usable metadata columns
  all_meta_cols <- names(metadata)
  meta_cols <- all_meta_cols[!grepl(paste0("^(", skip_prefixes, ")"), all_meta_cols)]
  meta_cols <- setdiff(meta_cols, c(file_col, ".file_match"))

  # Filter: skip single-value and high-cardinality columns
  meta_cols <- Filter(function(col) {
    vals <- merged[[col]]
    vals <- vals[!is.na(vals) & vals != ""]
    n_unique <- length(unique(vals))
    n_unique > 1 && n_unique <= max_unique
  }, meta_cols)

  if (length(meta_cols) == 0) {
    cat("  No metadata columns suitable for grouping.\n")
    return(invisible(NULL))
  }
  cat("  Metadata columns to plot:", paste(meta_cols, collapse = ", "), "\n")

  cluster_cols <- setdiff(names(tracking_matrix), "file")
  meta_out_dir <- file.path(output_dir, "by_metadata")
  make_outdir(meta_out_dir)

  for (meta_col in meta_cols) {
    safe_name <- create_safe_filename(meta_col)
    cat("    Group:", meta_col, "\n")

    # Clean up group values: replace empty/NA with "Unknown"
    group_vals <- as.character(merged[[meta_col]])
    group_vals[is.na(group_vals) | group_vals == ""] <- "Unknown"
    merged[[meta_col]] <- group_vals

    # --- 1. Faceted stacked barplot (samples within each group) ---
    long_df <- tidyr::pivot_longer(merged, cols = all_of(cluster_cols),
                                    names_to = "cluster", values_to = "count")
    long_df <- as.data.frame(long_df)

    # Sort samples by total count within each group
    sample_totals <- rowSums(merged[, cluster_cols, drop = FALSE])
    names(sample_totals) <- merged$file
    long_df$file <- factor(long_df$file,
                           levels = names(sort(sample_totals, decreasing = TRUE)))

    n_groups <- length(unique(group_vals))
    facet_cols <- min(n_groups, 5)
    facet_rows <- ceiling(n_groups / facet_cols)
    plot_h <- max(8, 6 * facet_rows)

    p_facet <- ggplot(long_df, aes(x = file, y = count, fill = cluster)) +
      geom_bar(stat = "identity", width = 1) +
      facet_wrap(as.formula(paste("~", paste0("`", meta_col, "`"))),
                 scales = "free_x", ncol = facet_cols) +
      theme_minimal() +
      labs(title = paste("Cluster Abundance by", meta_col),
           x = NULL, y = "Cell Count") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 2),
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.4, "cm"),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(size = 10, face = "bold")) +
      scale_fill_viridis_d(option = "turbo") +
      scale_y_continuous(labels = scales::comma)

    save_plot_formats(p_facet,
                      file.path(meta_out_dir, paste0("stacked_by_", safe_name)),
                      width = 28, height = plot_h, png_dpi = 150)

    # --- 2. Aggregated proportional barplot (one bar per group value) ---
    agg_df <- aggregate(
      long_df[, "count", drop = FALSE],
      by = list(group = long_df[[meta_col]], cluster = long_df$cluster),
      FUN = sum
    )

    # Order groups naturally (try numeric first, then alpha)
    group_levels <- unique(agg_df$group)
    numeric_check <- suppressWarnings(as.numeric(group_levels))
    if (!any(is.na(numeric_check))) {
      group_levels <- group_levels[order(numeric_check)]
    } else {
      group_levels <- sort(group_levels)
    }
    # Push "Unknown" to end if present
    if ("Unknown" %in% group_levels) {
      group_levels <- c(setdiff(group_levels, "Unknown"), "Unknown")
    }
    agg_df$group <- factor(agg_df$group, levels = group_levels)

    p_agg <- ggplot(agg_df, aes(x = group, y = count, fill = cluster)) +
      geom_bar(stat = "identity", position = "fill", width = 0.8) +
      theme_minimal() +
      labs(title = paste("Cluster Proportions by", meta_col),
           x = meta_col, y = "Proportion") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.text = element_text(size = 8)) +
      scale_fill_viridis_d(option = "turbo") +
      scale_y_continuous(expand = c(0, 0))

    save_plot_formats(p_agg,
                      file.path(meta_out_dir, paste0("proportions_by_", safe_name)),
                      width = max(10, n_groups * 1.5), height = 10)

    # --- 3. Aggregated stacked barplot (total counts per group) ---
    p_agg_stacked <- ggplot(agg_df, aes(x = group, y = count, fill = cluster)) +
      geom_bar(stat = "identity", width = 0.8) +
      theme_minimal() +
      labs(title = paste("Total Cluster Abundance by", meta_col),
           x = meta_col, y = "Cell Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.text = element_text(size = 8)) +
      scale_fill_viridis_d(option = "turbo") +
      scale_y_continuous(labels = scales::comma)

    save_plot_formats(p_agg_stacked,
                      file.path(meta_out_dir, paste0("totals_by_", safe_name)),
                      width = max(10, n_groups * 1.5), height = 10)
  }

  cat("  Metadata-grouped plots saved to:", meta_out_dir, "\n")
  invisible(NULL)
}
