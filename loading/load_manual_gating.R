#!/usr/bin/env Rscript
# =============================================================================
# loading/load_manual_gating.R — Load Manual Gating from CSV
# =============================================================================
# Refactored from: R_scripts/01_loading/Load_manual_gating.R
# Generalized from CD3-specific to configurable gate columns
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "fcs_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

#' Load manual gating labels from CSV files
#'
#' @param csv_dir Directory with CSV files containing manual gating columns
#' @param fcs_dir Directory with matching FCS files
#' @param gating_pattern Regex pattern identifying gating columns
#' @param label_dict Optional named vector for label translation
#' @param gating_tree Optional named numeric vector of gating depths
#' @param output_dir Output directory for label files
#' @return Named list of label vectors (one per sample)
load_manual_gating <- function(csv_dir,
                               fcs_dir = csv_dir,
                               gating_pattern = "Manual gating",
                               label_dict = NULL,
                               gating_tree = NULL,
                               output_dir = NULL) {
  require(data.table)

  csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
  fcs_files <- list.files(fcs_dir, pattern = "\\.fcs$", full.names = TRUE)

  if (length(csv_files) == 0) stop("No CSV files found in: ", csv_dir)

  # Filter to files with gating columns
  valid_csvs <- filter_valid_csv_files(csv_files, pattern = gating_pattern)
  if (length(valid_csvs) == 0) stop("No CSV files with gating columns found")

  # Map CSV to FCS
  csv_to_fcs <- if (length(fcs_files) > 0) {
    create_csv_fcs_mapping(valid_csvs, fcs_files)
  } else NULL

  if (!is.null(output_dir)) make_outdir(output_dir)

  results <- list()
  for (csv_file in valid_csvs) {
    cat("Processing:", basename(csv_file), "...")
    csv_data <- read_csv_filtered(csv_file, pattern = gating_pattern)
    if (is.null(csv_data)) { cat(" skipped\n"); next }

    # Get event count from FCS if available
    num_events <- if (!is.null(csv_to_fcs)) {
      get_fcs_event_count(csv_to_fcs[[csv_file]])
    } else nrow(csv_data)

    # Extract labels
    label_lists <- labels_from_binary_dt(csv_data)

    # Flatten to single label per cell
    if (!is.null(gating_tree)) {
      flat_labels <- flatten_label_list(label_lists, gating_tree)
    } else {
      flat_labels <- sapply(label_lists, function(x) {
        if (length(x) == 0) "Unassigned" else x[1]
      })
    }

    # Translate labels if dictionary provided
    if (!is.null(label_dict)) {
      flat_labels <- translate_labels(flat_labels, label_dict)
    }

    sample_name <- tools::file_path_sans_ext(basename(csv_file))
    results[[sample_name]] <- flat_labels

    # Save if output dir specified
    if (!is.null(output_dir)) {
      out_file <- file.path(output_dir, paste0(sample_name, "_labels.csv"))
      fwrite(data.table(cell_id = seq_along(flat_labels), label = flat_labels),
             out_file, sep = ";")
    }

    cat(sprintf(" %d cells, %d populations\n",
                length(flat_labels), length(unique(flat_labels))))
  }

  results
}
