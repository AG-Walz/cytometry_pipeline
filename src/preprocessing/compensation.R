#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/compensation.R — Fluorescence Compensation
# =============================================================================
# Applies spillover compensation to conventional flow cytometry data.
# Uses the compensation matrix from FCS $SPILL keyword or an external CSV.
# =============================================================================

#' Extract compensation matrix from FCS file
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R)
#' Output to:  apply_compensation() (this file)
#'
#' @param flowframe flowFrame object
#' @return Compensation matrix or NULL if not found
extract_compensation_matrix <- function(flowframe) {
  require(flowCore)

  # Try standard keyword names
  keywords <- c("SPILL", "SPILLOVER", "$SPILLOVER", "$SPILL",
                "COMP", "$COMP", "spillover")

  for (kw in keywords) {
    comp <- tryCatch(
      flowCore::keyword(flowframe, kw)[[1]],
      error = function(e) NULL
    )
    if (!is.null(comp) && is.matrix(comp)) {
      cat("Found compensation matrix from keyword '", kw, "': ",
          nrow(comp), "x", ncol(comp), "\n")
      return(comp)
    }
  }

  # Try spillover() convenience function
  comp <- tryCatch(flowCore::spillover(flowframe), error = function(e) NULL)
  if (!is.null(comp) && is.matrix(comp)) {
    cat("Found compensation matrix via spillover(): ",
        nrow(comp), "x", ncol(comp), "\n")
    return(comp)
  }

  warning("No compensation matrix found in flowFrame keywords.")
  NULL
}

#' Apply fluorescence compensation to a flowFrame
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R)
#' Output to:  quality_control() (preprocessing/quality_control.R),
#'             transform_flowframe() (preprocessing/transformation.R)
#'
#' @param flowframe flowFrame object
#' @param comp_matrix Compensation matrix (from FCS $SPILL or external CSV).
#'   If NULL and auto_extract = TRUE, extracts from flowFrame keywords.
#' @param auto_extract If TRUE, extract comp matrix from flowFrame's $SPILL keyword
#' @return Compensated flowFrame
apply_compensation <- function(flowframe, comp_matrix = NULL, auto_extract = TRUE) {
  require(flowCore)

  if (is.null(comp_matrix) && auto_extract) {
    comp_matrix <- extract_compensation_matrix(flowframe)
  }

  if (is.null(comp_matrix)) {
    stop("No compensation matrix provided and auto-extraction failed.")
  }

  # Validate matrix dimensions match available channels
  comp_channels <- colnames(comp_matrix)
  ff_channels <- colnames(flowframe)
  matched <- intersect(comp_channels, ff_channels)

  if (length(matched) == 0) {
    stop("No channels in compensation matrix match flowFrame channels.\n",
         "Matrix channels: ", paste(comp_channels[1:min(5, length(comp_channels))], collapse = ", "),
         "\nflowFrame channels: ", paste(ff_channels[1:min(5, length(ff_channels))], collapse = ", "))
  }

  if (length(matched) < length(comp_channels)) {
    warning("Only ", length(matched), " of ", length(comp_channels),
            " compensation channels found in flowFrame. Subsetting matrix.")
    idx <- which(comp_channels %in% matched)
    comp_matrix <- comp_matrix[idx, idx, drop = FALSE]
  }

  cat("Applying compensation (", nrow(comp_matrix), "x", ncol(comp_matrix),
      " matrix) to flowFrame with ", nrow(flowframe), " cells\n")

  compensated <- tryCatch(
    flowCore::compensate(flowframe, comp_matrix),
    error = function(e) {
      stop("Compensation failed: ", e$message)
    }
  )

  cat("Compensation complete.\n")
  compensated
}

#' Load compensation matrix from external CSV file
#'
#' @param csv_path Path to CSV file with compensation matrix
#' @param sep CSV separator (default ";")
#' @return Compensation matrix
load_compensation_matrix <- function(csv_path, sep = ";") {
  require(data.table)

  dt <- data.table::fread(csv_path, sep = sep, header = TRUE)
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[1]]
  colnames(mat) <- names(dt)[-1]

  cat("Loaded compensation matrix: ", nrow(mat), "x", ncol(mat), " from ", csv_path, "\n")
  mat
}
