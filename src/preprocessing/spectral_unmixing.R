#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/spectral_unmixing.R — Spectral Cytometry Unmixing
# =============================================================================
# Unmixes spectral cytometry data using reference spectra from single-stain
# controls. Uses non-negative least squares (NNLS) for robust unmixing.
# =============================================================================

#' Build unmixing matrix from single-stain controls
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R) — single-stain controls
#' Output to:  unmix_spectral() (this file)
#'
#' @param control_frames Named list of single-stain flowFrames (names = fluorochrome labels)
#' @param unstained_frame Unstained control flowFrame for autofluorescence baseline
#' @param channels Character vector of spectral channel names (NULL = auto-detect)
#' @return Unmixing matrix (channels x fluorochromes)
build_unmixing_matrix <- function(control_frames, unstained_frame, channels = NULL) {
  require(flowCore)

  if (is.null(channels)) {
    channels <- colnames(flowCore::exprs(unstained_frame))
    # Remove non-spectral channels
    exclude <- c("Time", "Event_length", "FSC-A", "FSC-H", "FSC-W",
                 "SSC-A", "SSC-H", "SSC-W")
    channels <- setdiff(channels, exclude)
  }

  cat("Building unmixing matrix from", length(control_frames),
      "single-stain controls across", length(channels), "channels\n")

  # Compute autofluorescence baseline (median of unstained)
  unstained_expr <- flowCore::exprs(unstained_frame)[, channels, drop = FALSE]
  autofluorescence <- apply(unstained_expr, 2, median, na.rm = TRUE)

  # Build reference spectra from each single-stain control
  n_fluors <- length(control_frames)
  unmix_mat <- matrix(0, nrow = length(channels), ncol = n_fluors)
  rownames(unmix_mat) <- channels
  colnames(unmix_mat) <- names(control_frames)

  for (i in seq_along(control_frames)) {
    fluor_name <- names(control_frames)[i]
    ff <- control_frames[[i]]
    expr <- flowCore::exprs(ff)[, channels, drop = FALSE]

    # Subtract autofluorescence and take median spectrum
    spectrum <- apply(expr, 2, median, na.rm = TRUE) - autofluorescence
    spectrum[spectrum < 0] <- 0

    # Normalize to unit sum
    spec_sum <- sum(spectrum)
    if (spec_sum > 0) {
      spectrum <- spectrum / spec_sum
    } else {
      warning("Zero spectrum for control '", fluor_name,
              "'. Check single-stain quality.")
    }

    unmix_mat[, i] <- spectrum
    cat("  ", fluor_name, ": peak channel =",
        channels[which.max(spectrum)], "\n")
  }

  cat("Unmixing matrix built:", nrow(unmix_mat), "channels x",
      ncol(unmix_mat), "fluorochromes\n")
  unmix_mat
}

#' Unmix spectral cytometry data
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R)
#' Output to:  quality_control() (preprocessing/quality_control.R),
#'             transform_flowframe() (preprocessing/transformation.R)
#'
#' @param flowframe flowFrame with spectral channels
#' @param unmixing_matrix Reference unmixing matrix (channels x fluorochromes)
#'   from build_unmixing_matrix()
#' @param autofluorescence Autofluorescence reference vector (optional, subtracted
#'   before unmixing). If NULL, no subtraction is performed.
#' @return Unmixed flowFrame with virtual marker channels
unmix_spectral <- function(flowframe, unmixing_matrix, autofluorescence = NULL) {
  require(flowCore)
  require(nnls)

  channels <- rownames(unmixing_matrix)
  fluors <- colnames(unmixing_matrix)

  # Validate channels exist in flowFrame
  ff_channels <- colnames(flowframe)
  missing <- setdiff(channels, ff_channels)
  if (length(missing) > 0) {
    stop("Missing channels in flowFrame: ", paste(missing, collapse = ", "))
  }

  expr <- flowCore::exprs(flowframe)[, channels, drop = FALSE]
  n_cells <- nrow(expr)
  n_fluors <- length(fluors)

  cat("Unmixing", n_cells, "cells across", length(channels),
      "channels into", n_fluors, "fluorochromes\n")

  # Subtract autofluorescence if provided
  if (!is.null(autofluorescence)) {
    expr <- sweep(expr, 2, autofluorescence[channels], "-")
    expr[expr < 0] <- 0
  }

  # NNLS unmixing: solve for each cell
  unmixed <- matrix(0, nrow = n_cells, ncol = n_fluors)
  colnames(unmixed) <- fluors

  # Process in chunks for memory efficiency
  chunk_size <- 10000
  n_chunks <- ceiling(n_cells / chunk_size)

  for (chunk in seq_len(n_chunks)) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, n_cells)
    chunk_data <- expr[start_idx:end_idx, , drop = FALSE]

    for (j in seq_len(nrow(chunk_data))) {
      result <- tryCatch(
        nnls::nnls(unmixing_matrix, chunk_data[j, ]),
        error = function(e) NULL
      )
      if (!is.null(result)) {
        unmixed[start_idx + j - 1, ] <- result$x
      }
    }

    if (n_chunks > 1) {
      cat("  Unmixed chunk", chunk, "/", n_chunks, "\n")
    }
  }

  # Preserve non-spectral channels (scatter, time, etc.)
  non_spectral <- setdiff(ff_channels, channels)
  if (length(non_spectral) > 0) {
    extra <- flowCore::exprs(flowframe)[, non_spectral, drop = FALSE]
    unmixed <- cbind(extra, unmixed)
  }

  # Build output flowFrame
  params <- Biobase::AnnotatedDataFrame(
    data.frame(
      name = colnames(unmixed),
      desc = colnames(unmixed),
      range = apply(unmixed, 2, function(x) diff(range(x, na.rm = TRUE))),
      minRange = apply(unmixed, 2, min, na.rm = TRUE),
      maxRange = apply(unmixed, 2, max, na.rm = TRUE),
      row.names = paste0("$P", seq_len(ncol(unmixed))),
      stringsAsFactors = FALSE
    )
  )

  result_ff <- flowCore::flowFrame(exprs = unmixed, parameters = params)
  cat("Unmixing complete:", ncol(unmixed), "output channels\n")
  result_ff
}
