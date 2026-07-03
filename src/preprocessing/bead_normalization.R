#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/bead_normalization.R — CyTOF Bead-Based Normalization
# =============================================================================
# Normalizes CyTOF/mass cytometry data using bead standards to correct
# for signal decay over acquisition time. Identifies bead events, computes
# per-channel normalization factors, and optionally removes bead events.
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

#' Normalize CyTOF data using bead standards
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R)
#' Output to:  quality_control() (preprocessing/quality_control.R),
#'             transform_flowframe() (preprocessing/transformation.R)
#'
#' @param flowframe flowFrame with bead and cell events
#' @param bead_channels Character vector of bead channel names
#'   (e.g., c("140Ce", "151Eu", "153Eu", "165Ho", "175Lu"))
#' @param bead_type "four_element" (standard EQ beads: Ce140, Eu151, Eu153,
#'   Ho165, Lu175) or "custom" (use bead_channels directly)
#' @param remove_beads Remove bead events after normalization (default TRUE)
#' @param baseline "median" (normalize to global median) or "start" (normalize
#'   to first time window)
#' @param window_size Number of events per smoothing window (default 500)
#' @param output_dir Optional directory for normalization diagnostic plots
#' @return list(normalized_ff, bead_indices, normalization_factors)
normalize_beads <- function(flowframe, bead_channels = NULL,
                            bead_type = "four_element",
                            remove_beads = TRUE,
                            baseline = "median",
                            window_size = 500,
                            output_dir = NULL) {
  require(flowCore)

  # Default bead channels for standard EQ Four Element beads
  if (is.null(bead_channels)) {
    bead_channels <- switch(bead_type,
      "four_element" = c("140Ce", "142Ce", "151Eu", "153Eu", "165Ho", "175Lu"),
      stop("bead_type must be 'four_element' or 'custom' (with bead_channels specified).")
    )
  }

  expr <- flowCore::exprs(flowframe)
  ff_channels <- colnames(expr)

  # Match bead channels (partial matching for metal names)
  matched_channels <- c()
  for (bc in bead_channels) {
    matches <- grep(bc, ff_channels, value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) {
      matched_channels <- c(matched_channels, matches[1])
    }
  }

  if (length(matched_channels) < 2) {
    stop("Found only ", length(matched_channels), " bead channels. Need at least 2.\n",
         "Searched for: ", paste(bead_channels, collapse = ", "), "\n",
         "Available: ", paste(ff_channels, collapse = ", "))
  }

  cat("Using", length(matched_channels), "bead channels:",
      paste(matched_channels, collapse = ", "), "\n")

  # Identify bead events using Mahalanobis distance on bead channels
  bead_expr <- expr[, matched_channels, drop = FALSE]
  bead_indices <- identify_beads(bead_expr)

  n_beads <- sum(bead_indices)
  cat("Identified", n_beads, "bead events (",
      round(100 * n_beads / nrow(expr), 1), "% of total)\n")

  if (n_beads < 100) {
    warning("Very few bead events (", n_beads, "). Normalization may be unreliable.")
  }

  # Compute normalization factors using sliding window on bead events
  bead_data <- expr[bead_indices, matched_channels, drop = FALSE]

  # Sort beads by acquisition order (row index = time proxy)
  bead_order <- which(bead_indices)

  # Compute running median per channel
  n_windows <- max(1, floor(n_beads / window_size))
  norm_factors <- matrix(1, nrow = nrow(expr), ncol = length(matched_channels))
  colnames(norm_factors) <- matched_channels

  for (ch_idx in seq_along(matched_channels)) {
    ch <- matched_channels[ch_idx]
    bead_vals <- expr[bead_order, ch]

    # Compute baseline
    baseline_val <- switch(baseline,
      "median" = median(bead_vals, na.rm = TRUE),
      "start"  = median(bead_vals[1:min(window_size, length(bead_vals))], na.rm = TRUE),
      median(bead_vals, na.rm = TRUE)
    )

    if (baseline_val <= 0) {
      warning("Non-positive baseline for channel '", ch, "'. Skipping.")
      next
    }

    # Sliding window medians
    window_medians <- numeric(n_beads)
    for (i in seq_len(n_beads)) {
      w_start <- max(1, i - window_size %/% 2)
      w_end <- min(n_beads, i + window_size %/% 2)
      window_medians[i] <- median(bead_vals[w_start:w_end], na.rm = TRUE)
    }

    # Interpolate normalization factors for all events
    bead_factors <- baseline_val / pmax(window_medians, 1e-10)
    all_factors <- approx(bead_order, bead_factors, xout = seq_len(nrow(expr)),
                          rule = 2)$y
    norm_factors[, ch_idx] <- all_factors
  }

  # Apply normalization to ALL channels (not just bead channels)
  # Use mean normalization factor across bead channels for non-bead channels
  mean_factor <- rowMeans(norm_factors)
  normalized_expr <- expr
  for (ch in ff_channels) {
    if (ch %in% matched_channels) {
      ch_idx <- which(matched_channels == ch)
      normalized_expr[, ch] <- expr[, ch] * norm_factors[, ch_idx]
    } else {
      normalized_expr[, ch] <- expr[, ch] * mean_factor
    }
  }

  # Remove beads if requested
  if (remove_beads) {
    normalized_expr <- normalized_expr[!bead_indices, , drop = FALSE]
    cat("Removed", n_beads, "bead events.",
        nrow(normalized_expr), "cells remaining.\n")
  }

  # Build output flowFrame
  params <- flowCore::parameters(flowframe)
  if (remove_beads) {
    result_ff <- flowCore::flowFrame(exprs = normalized_expr,
                                     parameters = params)
  } else {
    flowCore::exprs(flowframe) <- normalized_expr
    result_ff <- flowframe
  }

  # Save diagnostic plots if output dir specified
  if (!is.null(output_dir)) {
    make_outdir(output_dir)
    save_bead_diagnostics(expr, normalized_expr, bead_indices,
                          matched_channels, output_dir)
  }

  cat("Bead normalization complete.\n")
  list(
    normalized_ff = result_ff,
    bead_indices = bead_indices,
    normalization_factors = norm_factors
  )
}

#' Identify bead events using intensity thresholding
#'
#' @param bead_expr Matrix of bead channel expression values
#' @param method "mahalanobis" or "threshold"
#' @return Logical vector (TRUE = bead event)
identify_beads <- function(bead_expr, method = "mahalanobis") {
  n_cells <- nrow(bead_expr)

  if (method == "mahalanobis") {
    # Beads have high signal in ALL bead channels
    # Use Mahalanobis distance from the high-intensity cluster
    high_threshold <- apply(bead_expr, 2, quantile, probs = 0.9, na.rm = TRUE)
    candidate_beads <- apply(bead_expr, 1, function(x) all(x > high_threshold * 0.5))

    if (sum(candidate_beads) < 10) {
      # Fallback to simpler threshold
      warning("Mahalanobis bead identification found < 10 candidates. Using threshold method.")
      return(identify_beads(bead_expr, method = "threshold"))
    }

    # Refine using Mahalanobis distance from bead cluster center
    bead_center <- colMeans(bead_expr[candidate_beads, , drop = FALSE])
    bead_cov <- tryCatch(
      cov(bead_expr[candidate_beads, , drop = FALSE]),
      error = function(e) diag(ncol(bead_expr))
    )

    md <- tryCatch(
      mahalanobis(bead_expr, center = bead_center, cov = bead_cov),
      error = function(e) rep(Inf, n_cells)
    )

    # Chi-squared cutoff (p = 0.01)
    cutoff <- qchisq(0.99, df = ncol(bead_expr))
    return(md < cutoff)
  }

  # Threshold method: beads are bright in all bead channels
  thresholds <- apply(bead_expr, 2, quantile, probs = 0.95, na.rm = TRUE)
  apply(bead_expr, 1, function(x) all(x > thresholds * 0.3))
}

#' Save bead normalization diagnostic plots
#' @keywords internal
save_bead_diagnostics <- function(original_expr, normalized_expr,
                                   bead_indices, channels, output_dir) {
  tryCatch({
    require(ggplot2)

    for (ch in channels[1:min(3, length(channels))]) {
      if (!ch %in% colnames(original_expr)) next

      df <- data.frame(
        event = seq_len(nrow(original_expr)),
        original = original_expr[, ch],
        is_bead = bead_indices
      )

      p <- ggplot(df, aes(x = event, y = original, color = is_bead)) +
        geom_point(size = 0.1, alpha = 0.3) +
        scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
        theme_minimal() +
        labs(title = paste("Bead identification -", ch),
             x = "Event index", y = "Intensity")

      ggsave(file.path(output_dir, paste0("bead_diag_", ch, ".png")),
             plot = p, width = 10, height = 4, dpi = 150)
    }
    cat("Diagnostic plots saved to:", output_dir, "\n")
  }, error = function(e) {
    warning("Failed to save diagnostic plots: ", e$message)
  })
}
