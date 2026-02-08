#!/usr/bin/env Rscript
# =============================================================================
# metadata_qc/check_metadata.R — Metadata Barplots & Heatmaps
# =============================================================================
# Refactored from: R_scripts/00_utilities/Check_metadata.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Create compact metadata barplots for all clinical variables
#'
#' @param metadata_file Path to metadata CSV
#' @param output_dir Output directory
#' @param sep Separator for CSV (default ",")
#' @param skip_columns Character vector of columns to exclude
#' @param only_cols Character vector of columns to include (NULL = all non-skipped)
#' @param sort_column Column name to sort rows by before splitting (e.g. "Datum")
#' @param sort_desc Logical; sort descending?
#' @param n_splits Number of splits for batch-colored plot (default 3)
#' @param placeholders Character vector of strings to treat as "Unassigned"
#' @return Invisible NULL; saves plots and CSV to output_dir
check_metadata <- function(metadata_file, output_dir,
                           sep = ",",
                           skip_columns = NULL,
                           only_cols = NULL,
                           sort_column = NULL,
                           sort_desc = FALSE,
                           n_splits = 3,
                           placeholders = c("?", "na", "n/a", "null", "-")) {
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(forcats)
  require(stringr)
  require(scales)

  make_outdir(output_dir)
  stopifnot(file.exists(metadata_file))

  # --- Helper: convert blanks/placeholders to "Unassigned" ---
  to_unassigned <- function(x) {
    x_chr  <- as.character(x)
    x_trim <- stringr::str_trim(x_chr)
    blank  <- x_trim == "" | is.na(x_trim)
    low    <- tolower(x_trim)
    place  <- low %in% placeholders
    out <- x_trim
    out[blank | place] <- "Unassigned"
    out
  }

  md <- suppressWarnings(
    data.table::fread(metadata_file, sep = sep, header = TRUE, data.table = FALSE)
  )

  # --- Optional sorting ---
  if (!is.null(sort_column) && nzchar(sort_column) && sort_column %in% colnames(md)) {
    col_vals <- md[[sort_column]]
    sorted <- FALSE
    date_formats <- c("%d/%m/%Y", "%d-%m-%Y", "%d.%m.%Y", "%Y-%m-%d",
                      "%Y/%m/%d", "%m/%d/%Y", "%m-%d-%Y")
    for (fmt in date_formats) {
      parsed <- suppressWarnings(as.Date(col_vals, format = fmt))
      if (sum(!is.na(parsed)) > 0) {
        md <- md[order(parsed, decreasing = sort_desc, na.last = TRUE), , drop = FALSE]
        message("Sorted by '", sort_column, "' as Date (format: ", fmt, ")")
        sorted <- TRUE
        break
      }
    }
    if (!sorted) {
      md <- md[order(as.character(col_vals), decreasing = sort_desc, na.last = TRUE), , drop = FALSE]
      message("Sorted by '", sort_column, "' as character")
    }
  }

  # --- Column selection ---
  all_cols <- colnames(md)
  analyzable <- if (!is.null(skip_columns)) setdiff(all_cols, skip_columns) else all_cols
  analyze_cols <- if (!is.null(only_cols)) intersect(only_cols, analyzable) else analyzable
  if (!length(analyze_cols)) stop("No columns selected to analyze.")

  # --- Long format + counts ---
  long <- md |>
    dplyr::select(all_of(analyze_cols)) |>
    dplyr::mutate(across(everything(), ~ as.character(.))) |>
    tidyr::pivot_longer(everything(), names_to = "column", values_to = "value_raw") |>
    dplyr::mutate(value = to_unassigned(value_raw),
                  is_unassigned = value == "Unassigned")

  counts <- long |>
    dplyr::count(column, value, name = "Count") |>
    dplyr::group_by(column) |>
    dplyr::mutate(Total = sum(Count)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      is_unassigned = value == "Unassigned",
      value_tmp = dplyr::if_else(is_unassigned, "zzz__Unassigned__", value)
    ) |>
    dplyr::group_by(column) |>
    dplyr::mutate(value_tmp = forcats::fct_reorder(value_tmp, Count, .desc = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(value_tmp = forcats::fct_recode(value_tmp, Unassigned = "zzz__Unassigned__"))

  # Fill palette: grey for Unassigned
  counts <- counts |>
    dplyr::mutate(fill_key = dplyr::if_else(value_tmp == "Unassigned", "Unassigned", as.character(value_tmp)))
  fill_levels <- unique(counts$fill_key)
  other_levels <- setdiff(fill_levels, "Unassigned")
  pal <- c("Unassigned" = "grey70")
  if (length(other_levels) > 0) {
    pal <- c(pal, setNames(scales::hue_pal()(length(other_levels)), other_levels))
  }

  k_cols <- length(unique(counts$column))
  nrow_facets <- 2
  ncol_facets <- ceiling(k_cols / nrow_facets)

  p <- ggplot(counts, aes(x = value_tmp, y = Count, fill = fill_key)) +
    geom_col() +
    geom_text(aes(label = Count), vjust = -0.25, size = 3) +
    facet_wrap(~ column, scales = "free_x", nrow = nrow_facets) +
    scale_fill_manual(values = pal, name = "Category") +
    labs(x = NULL, y = "Count") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(size = 9, face = "bold")
    )

  plot_width  <- max(8, 2.1 * ncol_facets)
  plot_height <- max(6, 3.0 * nrow_facets)
  save_plot_formats(p, file.path(output_dir, "consensus_barplots_compact"),
                    width = plot_width, height = plot_height)

  # Save summary CSV
  data.table::fwrite(as.data.frame(counts[, c("column", "value", "Count")]),
                     file.path(output_dir, "metadata_summary_table.csv"))

  # --- Batch-split barplot ---
  if (n_splits > 1) {
    batch_base <- c("#D62728", "#2CA02C", "#1F77B4", "#FF7F0E", "#9467BD",
                    "#FFFF00", "#8C564B", "#17BECF", "#E377C2", "#90EE91")
    batch_levels <- paste0("Split ", seq_len(n_splits))
    batch_cols <- setNames(rep(batch_base, length.out = n_splits), batch_levels)

    md_split <- md |>
      dplyr::mutate(.row_id = dplyr::row_number(),
                    split = paste0("Split ", dplyr::ntile(.row_id, n_splits)))

    long_b <- md_split |>
      dplyr::select(all_of(analyze_cols), split) |>
      dplyr::mutate(across(all_of(analyze_cols), ~ as.character(.))) |>
      tidyr::pivot_longer(cols = all_of(analyze_cols),
                          names_to = "column", values_to = "value_raw") |>
      dplyr::mutate(value = to_unassigned(value_raw),
                    split = factor(split, levels = batch_levels))

    counts_b <- long_b |>
      dplyr::count(column, value, split, name = "Count")

    totals_b <- counts_b |>
      dplyr::group_by(column, value) |>
      dplyr::summarise(Total = sum(Count), .groups = "drop")

    counts_b <- counts_b |>
      dplyr::left_join(totals_b, by = c("column", "value")) |>
      dplyr::mutate(fill_key = as.character(split))

    p_b <- ggplot(counts_b, aes(x = value, y = Count, fill = fill_key)) +
      geom_col(position = position_stack(reverse = TRUE)) +
      facet_wrap(~ column, scales = "free_x", nrow = nrow_facets) +
      scale_fill_manual(values = batch_cols, name = "Batches") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
      labs(x = NULL, y = "Count") +
      theme_minimal(base_size = 9) +
      theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 9, face = "bold")
      )

    save_plot_formats(p_b, file.path(output_dir, sprintf("consensus_barplots_compact_bySplit_%d", n_splits)),
                      width = plot_width, height = plot_height)
  }

  cat("Metadata check plots saved to:", output_dir, "\n")
  invisible(NULL)
}

#' Create temporal heatmap of metadata features over time
#'
#' @param metadata_file Path to metadata CSV
#' @param output_dir Output directory
#' @param sep CSV separator
#' @param skip_columns Columns to exclude
#' @param date_column Name of the date column (default "Datum")
#' @param placeholders Strings to treat as "Unassigned"
plot_temporal_heatmap <- function(metadata_file, output_dir,
                                  sep = ",",
                                  skip_columns = NULL,
                                  date_column = "Datum",
                                  placeholders = c("?", "na", "n/a", "null", "-")) {
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(scales)

  make_outdir(output_dir)

  to_unassigned <- function(x) {
    x_chr <- stringr::str_trim(as.character(x))
    blank <- x_chr == "" | is.na(x_chr)
    place <- tolower(x_chr) %in% placeholders
    x_chr[blank | place] <- "Unassigned"
    x_chr
  }

  md <- suppressWarnings(
    data.table::fread(metadata_file, sep = sep, header = TRUE, data.table = FALSE)
  )

  if (!date_column %in% colnames(md)) {
    warning("Date column '", date_column, "' not found. Skipping temporal heatmap.")
    return(invisible(NULL))
  }

  all_cols <- colnames(md)
  analyze_cols <- if (!is.null(skip_columns)) setdiff(all_cols, skip_columns) else all_cols
  analyze_cols <- setdiff(analyze_cols, date_column)

  # Parse dates
  date_formats <- c("%d/%m/%Y", "%d-%m-%Y", "%d.%m.%Y", "%Y-%m-%d",
                    "%Y/%m/%d", "%m/%d/%Y", "%m-%d-%Y")
  parsed_dates <- NULL
  for (fmt in date_formats) {
    parsed_dates <- suppressWarnings(as.Date(md[[date_column]], format = fmt))
    if (sum(!is.na(parsed_dates)) > 0) break
  }

  md <- md[order(parsed_dates, na.last = TRUE), , drop = FALSE]
  timepoint_ids <- seq_len(nrow(md))
  timepoint_labels <- if (!is.null(parsed_dates) && !all(is.na(parsed_dates))) {
    format(sort(parsed_dates, na.last = TRUE), "%b %Y")
  } else {
    as.character(timepoint_ids)
  }

  heatmap_data <- md |>
    dplyr::mutate(timepoint_id = timepoint_ids, timepoint_label = timepoint_labels) |>
    dplyr::select(timepoint_id, timepoint_label, all_of(analyze_cols)) |>
    dplyr::mutate(across(all_of(analyze_cols), ~ as.character(.))) |>
    tidyr::pivot_longer(cols = all_of(analyze_cols),
                        names_to = "column", values_to = "value") |>
    dplyr::mutate(value = to_unassigned(value),
                  feature = paste0(column, ":", value))

  n_features <- length(unique(heatmap_data$feature))
  n_timepoints <- length(unique(heatmap_data$timepoint_id))

  hm <- heatmap_data |>
    dplyr::distinct(timepoint_id, feature) |>
    dplyr::mutate(present = 1)

  feature_order <- heatmap_data |>
    dplyr::count(column, feature, name = "freq") |>
    dplyr::arrange(column, desc(freq)) |>
    dplyr::pull(feature)

  hm <- hm |>
    dplyr::mutate(
      feature = factor(feature, levels = rev(feature_order)),
      column = sapply(strsplit(as.character(feature), ":"), `[`, 1),
      value = sapply(strsplit(as.character(feature), ":"), function(x) if (length(x) > 1) x[2] else ""),
      is_unassigned = value == "Unassigned"
    )

  unique_cols <- unique(hm$column)
  col_colors <- setNames(scales::hue_pal()(length(unique_cols)), unique_cols)
  hm$fill_color <- ifelse(hm$is_unassigned, "grey70", col_colors[hm$column])

  label_interval <- max(1, floor(n_timepoints / 20))
  x_breaks <- seq(1, n_timepoints, by = label_interval)

  plot_height_heat <- min(20, max(8, n_features * 0.25))

  p_heat <- ggplot(hm, aes(x = timepoint_id, y = feature, fill = fill_color)) +
    geom_tile(color = NA) +
    scale_fill_identity() +
    scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Timepoint", y = "Feature (Column:Value)",
         title = "Temporal Distribution of Metadata Features") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )

  save_plot_formats(p_heat, file.path(output_dir, "temporal_feature_heatmap"),
                    width = 14, height = plot_height_heat)

  cat("Temporal heatmap saved to:", output_dir, "\n")
  invisible(NULL)
}
