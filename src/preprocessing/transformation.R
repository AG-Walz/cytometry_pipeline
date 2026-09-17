#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/transformation.R — Data Transformation Functions
# =============================================================================
# Extracted from: R_scripts/02_preprocessing/Preprocessing_functions.R
# =============================================================================

#' Transform a flowFrame using specified method
#'
#' @param flowframe A flowFrame object
#' @param sample_name Sample identifier
#' @param markers_of_interest Character vector of marker names to transform
#' @param offset_list Numeric vector of offsets (one per marker). When NULL,
#'   defaults to 0 for all markers (cofactor-style arcsinh).
#' @param width_list Numeric vector of widths (one per marker). When NULL,
#'   defaults to a sensible per-transform value: 1/5 for arcsinh (cofactor=5),
#'   1 for log, ignored for logicle/biexponential.
#' @param transform_type "arcsinh", "logicle", "biexponential", "log", or "None"
#' @return Transformed flowFrame
transform_flowframe <- function(flowframe, sample_name,
                                markers_of_interest,
                                offset_list = NULL, width_list = NULL,
                                transform_type = "arcsinh") {
  valid_transforms <- c("arcsinh", "logicle", "biexponential", "log", "None")
  if (!transform_type %in% valid_transforms) {
    stop("Invalid transform_type: '", transform_type, "'. Must be one of: ",
         paste(valid_transforms, collapse = ", "))
  }
  if (is.null(offset_list)) {
    offset_list <- rep(0, length(markers_of_interest))
  }
  if (is.null(width_list)) {
    default_width <- switch(transform_type,
                            "arcsinh" = 1 / 5,    # cofactor = 5
                            "log"     = 1,
                            1)                     # ignored by logicle/biexponential
    width_list <- rep(default_width, length(markers_of_interest))
  }
  if (length(markers_of_interest) != length(offset_list)) {
    stop("Marker and offset lists must have the same length.")
  }
  if (length(markers_of_interest) != length(width_list)) {
    stop("Marker and width lists must have the same length.")
  }

  cat(transform_type, "transform sample", sample_name, "\n")
  tf_ff <- flowframe

  for (i in seq_along(markers_of_interest)) {
    marker <- markers_of_interest[i]
    offset <- offset_list[i]
    width <- width_list[i]

    if (!(marker %in% flowCore::colnames(flowframe))) {
      warning("Skipping marker '", marker, "' — not found in flowFrame.")
      next
    }

    log_offset <- offset  # capture for closure
    trans <- switch(transform_type,
      "arcsinh"        = flowCore::arcsinhTransform(a = -offset, b = 1 / width),
      "logicle"        = flowCore::logicleTransform(),
      "biexponential"  = flowCore::biexponentialTransform(),
      "log"            = function(x) log10(pmax(x + log_offset, 1e-10)),
      "None"           = NULL,
      stop("Invalid transform_type: '", transform_type, "'.")
    )

    if (!is.null(trans)) {
      tflist <- flowCore::transformList(from = marker, tfun = trans)
      tf_ff <- flowCore::transform(tf_ff, tflist, truncate_max_range = FALSE)
    }
  }

  tf_ff
}

#' Optimize arcsinh cofactor for bimodality
#'
#' @param data Numeric vector of expression values
#' @param target_bimodality Target bimodality coefficient
#' @return Optimal b cofactor value
find_optimal_b_cofactor <- function(data, target_bimodality_coefficient = 0.9) {
  if (!requireNamespace("mousetrap", quietly = TRUE)) {
    warning("mousetrap package not available; using default cofactor of 5")
    return(5)
  }
  calculate_bc <- function(b_cofactor, data, target) {
    transformed <- asinh(data * b_cofactor)
    bc <- mousetrap::bimodality_coefficient(transformed)
    abs(bc - target)
  }
  result <- optimize(
    f = calculate_bc,
    interval = c(0.001, 10000),
    data = data,
    target = target_bimodality_coefficient,
    maximum = FALSE
  )
  result$minimum
}
