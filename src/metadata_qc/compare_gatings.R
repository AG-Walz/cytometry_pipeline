#!/usr/bin/env Rscript
# =============================================================================
# metadata_qc/compare_gatings.R — Compare Gating Strategies
# =============================================================================
# Refactored from: R_scripts/05_analysis/Compare_gatings.R
# =============================================================================

if (!exists("shared_dir")) shared_dir <- { .s <- Sys.getenv("CYTOMETRY_SRC", ""); if (nzchar(.s)) file.path(.s, "R_shared") else file.path(dirname(getwd()), "src", "R_shared") }
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

#' Compare gating strategies across FlowJo workspaces
#'
#' Loads label lists from multiple FlowJo workspaces, translates labels
#' using a user-provided dictionary, and generates Venn/Euler diagrams
#' for each population.
#'
#' @param pop_strings_list Named list of population string lists (from load_flowjo_wsp)
#' @param target_fcs FCS filename to compare across workspaces
#' @param label_dict Named character vector mapping gate paths to simplified labels
#' @param custom_name Character vector of custom names for each workspace
#' @param output_dir Output directory for plots and CSVs
#' @param idx_labels Optional named list of CSV-derived label vectors for comparison
#' @param scale_size Use Euler (TRUE) or Venn (FALSE) diagrams
#' @return Invisible data frame of combined label lists
compare_gatings <- function(pop_strings_list, target_fcs, label_dict,
                             custom_name = NULL, output_dir,
                             idx_labels = NULL, scale_size = TRUE) {
  require(data.table)

  # Source visualization functions using shared_dir pattern
  if (!exists("shared_dir")) {
    shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
  }
  viz_dir <- file.path(dirname(shared_dir), "visualization")
  source(file.path(viz_dir, "venn_diagrams.R"))

  make_outdir(output_dir)

  # Check target FCS presence in all workspaces
  missing <- sapply(pop_strings_list, function(x) {
    !any(grepl(target_fcs, names(x), fixed = TRUE))
  })
  if (any(missing)) {
    stop("Target FCS '", target_fcs, "' missing from: ",
         paste(names(pop_strings_list)[missing], collapse = ", "))
  }
  message("Target FCS '", target_fcs, "' found in all workspaces.")

  # Extract label lists for target FCS
  fj_label_lists <- sapply(pop_strings_list, function(x) {
    idx <- grep(target_fcs, names(x), fixed = TRUE)
    if (length(idx) > 1) {
      warning("Multiple matches for '", target_fcs, "', using first.")
      idx <- idx[1]
    }
    x[[idx]]
  })

  # Translate labels
  for (i in seq_along(fj_label_lists)) {
    fj_label_lists[[i]] <- translate_labels(fj_label_lists[[i]], label_dict)
  }

  # Validate that all vectors have the same length before converting to data frame
  lengths <- sapply(fj_label_lists, length)
  if (length(unique(lengths)) > 1) {
    warning("Workspaces have different cell counts: ", paste(lengths, collapse = ", "),
            ". Truncating to minimum length for comparison.")
    min_len <- min(lengths)
    fj_label_lists <- lapply(fj_label_lists, function(x) x[seq_len(min_len)])
  }

  # Save FlowJo label lists
  fj_df <- as.data.frame(fj_label_lists, stringsAsFactors = FALSE)
  data.table::fwrite(fj_df, file.path(output_dir, paste0("fj_label_lists_", target_fcs, ".csv")))

  label_lists <- fj_df

  # Optionally merge CSV-derived labels
  if (!is.null(idx_labels)) {
    target_base <- gsub("\\.fcs$", "", target_fcs)
    target_norm <- gsub("[-_]", "", target_base)
    csv_names_norm <- gsub("[-_]", "", names(idx_labels))
    match_idx <- which(grepl(target_norm, csv_names_norm, fixed = TRUE))

    if (length(match_idx) > 0) {
      cy_lists <- lapply(match_idx, function(idx) {
        tmp <- idx_labels[[idx]]
        tmp[tmp == ""] <- "Removed"
        tmp
      })
      label_lists <- cbind(label_lists, do.call(cbind, cy_lists))
    } else {
      warning("No CSV labels matching '", target_base, "'")
    }
  }

  # Save combined label lists
  data.table::fwrite(as.data.frame(label_lists),
                     file.path(output_dir, paste0("label_lists_", target_fcs, ".csv")))

  # Generate Venn/Euler for each population
  unique_pops <- unique(as.vector(unlist(label_lists)))
  for (population in unique_pops) {
    if (is.na(population) || population == "") next
    cat("Plotting:", population, "\n")
    plot_and_save_venn(
      lists = label_lists, population = population,
      filename = paste0("subsetted ", population, "_", target_fcs),
      output_dir = output_dir, scale_size = scale_size,
      custom_name = custom_name
    )
  }

  cat("Gating comparison complete. Results in:", output_dir, "\n")
  invisible(label_lists)
}
