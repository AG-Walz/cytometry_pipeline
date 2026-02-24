#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/doublet_detection.R — Doublet/Multiplet Detection
# =============================================================================
# Detects and removes doublets using forward scatter area vs. height ratios
# or bivariate density gating on scatter parameters.
# =============================================================================

#' Detect and remove doublets using scatter parameters
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             apply_compensation() (preprocessing/compensation.R)
#' Output to:  quality_control() (preprocessing/quality_control.R),
#'             transform_flowframe() (preprocessing/transformation.R)
#'
#' @param flowframe flowFrame object
#' @param area_col FSC-A (forward scatter area) column name
#' @param height_col FSC-H (forward scatter height) column name
#' @param method "ratio" (FSC-A/FSC-H threshold) or "model" (bivariate density gating)
#' @param threshold Ratio threshold for "ratio" method (default 1.1).
#'   Events with FSC-A/FSC-H > threshold are classified as doublets.
#' @param width_col Optional FSC-W (forward scatter width) column for additional filtering
#' @return list(filtered_ff, doublet_indices, singlet_rate)
detect_doublets <- function(flowframe, area_col = "FSC-A", height_col = "FSC-H",
                            method = "ratio", threshold = 1.1,
                            width_col = NULL) {
  require(flowCore)

  expr <- flowCore::exprs(flowframe)
  ff_channels <- colnames(expr)

  # Validate required columns exist
  if (!area_col %in% ff_channels) {
    # Try partial matching
    area_match <- grep("FSC.*A", ff_channels, value = TRUE, ignore.case = TRUE)
    if (length(area_match) > 0) {
      area_col <- area_match[1]
      cat("Using", area_col, "as area channel\n")
    } else {
      stop("Area column '", area_col, "' not found. Available: ",
           paste(ff_channels, collapse = ", "))
    }
  }

  if (!height_col %in% ff_channels) {
    height_match <- grep("FSC.*H", ff_channels, value = TRUE, ignore.case = TRUE)
    if (length(height_match) > 0) {
      height_col <- height_match[1]
      cat("Using", height_col, "as height channel\n")
    } else {
      stop("Height column '", height_col, "' not found. Available: ",
           paste(ff_channels, collapse = ", "))
    }
  }

  area <- expr[, area_col]
  height <- expr[, height_col]
  n_cells <- nrow(expr)

  cat("Doublet detection (method = '", method, "') on ", n_cells, " events\n", sep = "")

  is_doublet <- switch(method,
    "ratio" = {
      # Ratio method: doublets have disproportionate area vs height
      ratio <- area / pmax(height, 1e-10)
      ratio > threshold
    },
    "model" = {
      # Bivariate density gating: fit linear model, flag outliers
      valid <- height > 0 & area > 0
      fit <- tryCatch({
        lm(area[valid] ~ height[valid])
      }, error = function(e) {
        warning("Linear model failed, falling back to ratio method.")
        return(NULL)
      })

      if (is.null(fit)) {
        ratio <- area / pmax(height, 1e-10)
        ratio > threshold
      } else {
        predicted <- predict(fit, newdata = data.frame(height = height))
        residuals <- area - predicted
        sd_resid <- sd(residuals[valid], na.rm = TRUE)
        abs(residuals) > 2 * sd_resid
      }
    },
    stop("Invalid method: '", method, "'. Use 'ratio' or 'model'.")
  )

  # Optional: additional width-based filtering
  if (!is.null(width_col) && width_col %in% ff_channels) {
    width <- expr[, width_col]
    width_q <- quantile(width, probs = c(0.01, 0.99), na.rm = TRUE)
    width_doublets <- width < width_q[1] | width > width_q[2]
    is_doublet <- is_doublet | width_doublets
    cat("  Width-based filtering removed", sum(width_doublets & !is_doublet),
        "additional events\n")
  }

  n_doublets <- sum(is_doublet)
  singlet_rate <- 1 - n_doublets / n_cells

  cat("  Detected", n_doublets, "doublets (",
      round(100 * (1 - singlet_rate), 1), "%)\n")
  cat("  Singlet rate:", round(100 * singlet_rate, 1), "%\n")

  # Filter to singlets
  filtered_expr <- expr[!is_doublet, , drop = FALSE]
  params <- flowCore::parameters(flowframe)
  filtered_ff <- flowCore::flowFrame(exprs = filtered_expr, parameters = params)

  list(
    filtered_ff = filtered_ff,
    doublet_indices = which(is_doublet),
    singlet_rate = singlet_rate
  )
}
