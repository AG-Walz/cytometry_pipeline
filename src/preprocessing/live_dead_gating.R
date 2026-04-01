#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/live_dead_gating.R — Live/Dead Cell Gating
# =============================================================================
# Gates on live cells using a viability marker (e.g., Live/Dead Aqua,
# Zombie dyes, propidium iodide). Supports threshold-based and Gaussian
# mixture model-based methods.
# =============================================================================

#' Gate on live cells using viability marker
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             detect_doublets() (preprocessing/doublet_detection.R)
#' Output to:  quality_control() (preprocessing/quality_control.R),
#'             transform_flowframe() (preprocessing/transformation.R)
#'
#' @param flowframe flowFrame object
#' @param viability_channel Name of live/dead marker column
#' @param method "threshold" (fixed cutoff) or "mixture" (Gaussian mixture model)
#' @param threshold Manual threshold (for "threshold" method). If NULL with
#'   threshold method, auto-detects using valley between bimodal peaks.
#' @param dead_high If TRUE (default), dead cells are viability-HIGH
#'   (standard for amine-reactive dyes). Set FALSE for dyes where dead = low.
#' @return list(live_ff, dead_indices, viability_rate)
gate_live_cells <- function(flowframe, viability_channel,
                            method = "threshold",
                            threshold = NULL,
                            dead_high = TRUE) {
  require(flowCore)

  expr <- flowCore::exprs(flowframe)
  ff_channels <- colnames(expr)

  # Find viability channel
  if (!viability_channel %in% ff_channels) {
    via_match <- grep(viability_channel, ff_channels, value = TRUE,
                      ignore.case = TRUE)
    if (length(via_match) > 0) {
      viability_channel <- via_match[1]
      cat("Using", viability_channel, "as viability channel\n")
    } else {
      stop("Viability channel '", viability_channel, "' not found.\n",
           "Available channels: ", paste(ff_channels, collapse = ", "))
    }
  }

  via_data <- expr[, viability_channel]
  n_cells <- nrow(expr)

  cat("Live/dead gating (method = '", method, "') on ", n_cells,
      " events using '", viability_channel, "'\n", sep = "")

  cutoff <- switch(method,
    "threshold" = {
      if (is.null(threshold)) {
        # Auto-detect threshold from density valley
        threshold <- find_density_valley(via_data)
        cat("  Auto-detected threshold:", round(threshold, 2), "\n")
      }
      threshold
    },
    "mixture" = {
      # Gaussian mixture model with 2 components
      gmm_threshold <- fit_gmm_threshold(via_data)
      cat("  GMM-estimated threshold:", round(gmm_threshold, 2), "\n")
      gmm_threshold
    },
    stop("Invalid method: '", method, "'. Use 'threshold' or 'mixture'.")
  )

  # Classify cells
  if (dead_high) {
    is_dead <- via_data > cutoff
  } else {
    is_dead <- via_data < cutoff
  }

  n_dead <- sum(is_dead)
  viability_rate <- 1 - n_dead / n_cells

  cat("  Dead cells:", n_dead, "(", round(100 * n_dead / n_cells, 1), "%)\n")
  cat("  Viability rate:", round(100 * viability_rate, 1), "%\n")

  if (viability_rate < 0.5) {
    warning("Low viability rate (", round(100 * viability_rate, 1),
            "%). Check viability channel and threshold.")
  }

  # Filter to live cells
  live_expr <- expr[!is_dead, , drop = FALSE]
  params <- flowCore::parameters(flowframe)
  live_ff <- flowCore::flowFrame(exprs = live_expr, parameters = params)

  list(
    live_ff = live_ff,
    dead_indices = which(is_dead),
    viability_rate = viability_rate
  )
}

#' Find valley between two density peaks (for auto-threshold)
#' @param x Numeric vector
#' @return Threshold value at density valley
#' @keywords internal
find_density_valley <- function(x) {
  x <- x[is.finite(x)]
  d <- density(x, n = 512)

  # Find peaks
  peaks <- which(diff(sign(diff(d$y))) == -2) + 1

  if (length(peaks) < 2) {
    # No bimodal pattern: use median
    warning("No bimodal pattern detected. Using median as threshold.")
    return(median(x, na.rm = TRUE))
  }

  # Find valley between two highest peaks
  peak_heights <- d$y[peaks]
  top2 <- order(peak_heights, decreasing = TRUE)[1:2]
  peak_range <- sort(peaks[top2])

  valley_region <- d$y[peak_range[1]:peak_range[2]]
  valley_idx <- peak_range[1] + which.min(valley_region) - 1

  d$x[valley_idx]
}

#' Fit 2-component Gaussian mixture model for threshold estimation
#' @param x Numeric vector
#' @return Threshold at intersection of two Gaussians
#' @keywords internal
fit_gmm_threshold <- function(x) {
  x <- x[is.finite(x)]

  # Simple EM for 2-component mixture (avoid mclust dependency)
  # Initialize with k-means
  km <- tryCatch(
    kmeans(x, centers = 2, nstart = 10),
    error = function(e) NULL
  )

  if (is.null(km)) {
    warning("GMM fitting failed. Falling back to density valley.")
    return(find_density_valley(x))
  }

  mu <- sort(km$centers[, 1])
  sigma <- c(sd(x[km$cluster == which.min(km$centers)]),
             sd(x[km$cluster == which.max(km$centers)]))
  sigma[sigma < 1e-10] <- sd(x) / 10

  # EM iterations
  for (iter in 1:50) {
    # E-step
    d1 <- dnorm(x, mean = mu[1], sd = sigma[1])
    d2 <- dnorm(x, mean = mu[2], sd = sigma[2])
    gamma <- d2 / (d1 + d2 + 1e-300)

    # M-step
    mu[1] <- weighted.mean(x, 1 - gamma)
    mu[2] <- weighted.mean(x, gamma)
    sigma[1] <- sqrt(sum((1 - gamma) * (x - mu[1])^2) / sum(1 - gamma))
    sigma[2] <- sqrt(sum(gamma * (x - mu[2])^2) / sum(gamma))
    sigma[sigma < 1e-10] <- sd(x) / 10
  }

  # Threshold at intersection of the two Gaussians
  # Use midpoint weighted by standard deviations
  threshold <- (mu[1] * sigma[2] + mu[2] * sigma[1]) / (sigma[1] + sigma[2])
  threshold
}
