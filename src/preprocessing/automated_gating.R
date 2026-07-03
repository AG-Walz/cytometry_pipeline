#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/automated_gating.R — Automated Gating Methods
# =============================================================================
# Data-driven automated gating using flowDensity (density-based cutpoints)
# and openCyto (template-based gating). Also provides programmatic gate
# construction for rectangle, polygon, threshold, and quadrant gates.
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

#' Automated gating using data-driven methods
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             apply_compensation() (preprocessing/compensation.R)
#' Output to:  quality_control() (preprocessing/quality_control.R),
#'             apply_gate_sequence() (preprocessing/gate_utils.R),
#'             downstream analysis
#'
#' @param flowframe flowFrame object
#' @param method "flowDensity" (density-based cutpoints) or "openCyto" (template-based)
#' @param gating_template Path to openCyto gating template CSV (for "openCyto" method)
#' @param channels Character vector of channels to gate on
#' @param population Target population name (e.g., "CD3+CD4+")
#' @return list(gated_ff, gate_coords, population_stats)
run_auto_gating <- function(flowframe, method = "flowDensity",
                            gating_template = NULL, channels,
                            population = NULL) {
  require(flowCore)

  expr <- flowCore::exprs(flowframe)
  ff_channels <- colnames(expr)
  n_cells <- nrow(expr)

  # --- Validate channels -----------------------------------------------------
  missing_ch <- setdiff(channels, ff_channels)
  if (length(missing_ch) > 0) {
    stop("Channels not found in flowFrame: ",
         paste(missing_ch, collapse = ", "),
         "\nAvailable: ", paste(ff_channels, collapse = ", "))
  }

  cat("Automated gating (method = '", method, "') on ", n_cells,
      " events using ", length(channels), " channel(s)\n", sep = "")

  result <- switch(method,
    "flowDensity" = .run_flowdensity_gating(flowframe, channels, population),
    "openCyto"    = .run_opencyto_gating(flowframe, gating_template, channels,
                                          population),
    stop("Invalid method: '", method,
         "'. Use 'flowDensity' or 'openCyto'.")
  )

  result
}

#' flowDensity-based sequential gating (internal)
#' @keywords internal
.run_flowdensity_gating <- function(flowframe, channels, population) {
  require(flowCore)
  require(flowDensity)

  n_cells <- nrow(flowframe)
  gate_coords <- list()
  current_ff <- flowframe
  included <- rep(TRUE, n_cells)

  for (ch in channels) {
    cat("  Channel '", ch, "': ", sep = "")

    fd_result <- tryCatch({
      flowDensity::flowDensity(
        obj    = current_ff,
        channels = ch,
        position = TRUE
      )
    }, error = function(e) {
      warning("flowDensity failed for channel '", ch, "': ", e$message,
              ". Falling back to density-valley threshold.")
      NULL
    })

    if (!is.null(fd_result)) {
      # Extract threshold from flowDensity result
      threshold <- fd_result@gates
      cat("threshold =", round(threshold, 2))

      # Get indices of positive cells relative to current_ff
      positive <- flowCore::exprs(current_ff)[, ch] > threshold
      gate_coords[[ch]] <- list(
        channel   = ch,
        threshold = threshold,
        direction = "positive"
      )
    } else {
      # Fallback: density valley
      values <- flowCore::exprs(current_ff)[, ch]
      threshold <- .find_density_cutpoint(values)
      cat("fallback threshold =", round(threshold, 2))

      positive <- values > threshold
      gate_coords[[ch]] <- list(
        channel   = ch,
        threshold = threshold,
        direction = "positive",
        method    = "density_valley_fallback"
      )
    }

    n_pos <- sum(positive)
    cat(" ->", n_pos, "/", nrow(current_ff), "positive (",
        round(100 * n_pos / nrow(current_ff), 1), "%)\n")

    # Update included mask on the original flowframe indices
    current_included <- which(included)
    included[current_included[!positive]] <- FALSE

    # Subset for next iteration
    current_ff <- flowframe[included, ]
  }

  n_gated <- sum(included)
  pct_parent <- 100 * n_gated / n_cells

  cat("  Final gated population: ", n_gated, " / ", n_cells,
      " (", round(pct_parent, 1), "%)\n", sep = "")

  gated_ff <- flowframe[included, ]

  list(
    gated_ff = gated_ff,
    gate_coords = gate_coords,
    population_stats = list(
      n_parent  = n_cells,
      n_gated   = n_gated,
      pct_parent = pct_parent,
      population = population,
      method     = "flowDensity"
    )
  )
}

#' openCyto template-based gating (internal)
#' @keywords internal
.run_opencyto_gating <- function(flowframe, gating_template, channels,
                                  population) {
  require(flowCore)
  require(flowWorkspace)
  require(openCyto)

  if (is.null(gating_template)) {
    stop("gating_template path is required for method 'openCyto'.")
  }
  if (!file.exists(gating_template)) {
    stop("Gating template file not found: ", gating_template)
  }

  cat("  Loading gating template:", gating_template, "\n")

  n_cells <- nrow(flowframe)

  result <- tryCatch({
    # Create GatingSet from flowFrame
    fs <- flowCore::flowSet(list(sample = flowframe))
    gs <- flowWorkspace::GatingSet(fs)

    # Load and apply template
    gt <- openCyto::gatingTemplate(gating_template)
    openCyto::gt_gating(gt, gs)

    # Extract population
    if (!is.null(population)) {
      node <- population
    } else {
      # Use the last node in the template
      nodes <- flowWorkspace::gs_get_pop_paths(gs)
      node <- nodes[length(nodes)]
      cat("  Using last gating node:", node, "\n")
    }

    gated_idx <- flowWorkspace::gh_pop_get_indices(gs[[1]], node)
    gated_ff <- flowframe[gated_idx, ]
    n_gated <- sum(gated_idx)

    # Collect gate coordinates from the GatingSet
    gate_obj <- flowWorkspace::gh_pop_get_gate(gs[[1]], node)
    gate_coords <- list()
    gate_coords[[node]] <- list(
      gate   = gate_obj,
      node   = node,
      method = "openCyto"
    )

    list(
      gated_ff = gated_ff,
      gate_coords = gate_coords,
      population_stats = list(
        n_parent   = n_cells,
        n_gated    = n_gated,
        pct_parent = 100 * n_gated / n_cells,
        population = node,
        method     = "openCyto"
      )
    )
  }, error = function(e) {
    stop("openCyto gating failed: ", e$message,
         "\nCheck gating template format and channel names.")
  })

  cat("  openCyto gated population: ", result$population_stats$n_gated,
      " / ", n_cells, " (", round(result$population_stats$pct_parent, 1),
      "%)\n", sep = "")

  result
}

#' Find density-based cutpoint between populations (internal fallback)
#' @param x Numeric vector
#' @return Numeric threshold
#' @keywords internal
.find_density_cutpoint <- function(x) {
  x <- x[is.finite(x)]
  d <- density(x, n = 512)

  # Find peaks
  peaks <- which(diff(sign(diff(d$y))) == -2) + 1

  if (length(peaks) < 2) {
    warning("No bimodal pattern detected. Using median as cutpoint.")
    return(median(x, na.rm = TRUE))
  }

  # Valley between the two highest peaks
  peak_heights <- d$y[peaks]
  top2 <- order(peak_heights, decreasing = TRUE)[1:2]
  peak_range <- sort(peaks[top2])

  valley_region <- d$y[peak_range[1]:peak_range[2]]
  valley_idx <- peak_range[1] + which.min(valley_region) - 1

  d$x[valley_idx]
}

#' Define manual gates programmatically
#'
#' Constructs flowCore gate objects from specified coordinates. Supports
#' rectangle, polygon, threshold, and quadrant gate types.
#'
#' Input from: user-defined gate coordinates
#' Output to:  apply_gate_sequence() (preprocessing/gate_utils.R),
#'             run_auto_gating() (this file)
#'
#' @param gate_type "rectangle", "polygon", "threshold", "quadrant"
#' @param channels Channel names (1 for threshold, 2 for rectangle/polygon/quadrant)
#' @param coords Gate coordinates. Format depends on gate_type:
#'   - "rectangle": list(x = c(min, max), y = c(min, max))
#'   - "polygon": matrix with columns x, y (one row per vertex)
#'   - "threshold": numeric scalar (cutpoint on the single channel)
#'   - "quadrant": list(x = numeric, y = numeric) — intersection point
#' @return flowCore gate object (or list of 4 gates for "quadrant")
define_gate <- function(gate_type, channels, coords) {
  require(flowCore)

  gate_type <- tolower(gate_type)

  result <- switch(gate_type,
    "rectangle" = {
      if (length(channels) != 2) {
        stop("Rectangle gate requires exactly 2 channels, got ", length(channels))
      }
      if (!all(c("x", "y") %in% names(coords))) {
        stop("Rectangle coords must be a list with x = c(min, max) and y = c(min, max)")
      }
      bounds <- list()
      bounds[[channels[1]]] <- coords$x
      bounds[[channels[2]]] <- coords$y
      gate <- flowCore::rectangleGate(.gate = bounds)
      cat("Rectangle gate on", channels[1], "[",
          coords$x[1], ",", coords$x[2], "] x",
          channels[2], "[", coords$y[1], ",", coords$y[2], "]\n")
      gate
    },

    "polygon" = {
      if (length(channels) != 2) {
        stop("Polygon gate requires exactly 2 channels, got ", length(channels))
      }
      if (!is.matrix(coords) || ncol(coords) != 2) {
        stop("Polygon coords must be a matrix with 2 columns (x, y vertices)")
      }
      colnames(coords) <- channels
      gate <- flowCore::polygonGate(.gate = coords)
      cat("Polygon gate on", channels[1], "/", channels[2],
          "with", nrow(coords), "vertices\n")
      gate
    },

    "threshold" = {
      if (length(channels) != 1) {
        stop("Threshold gate requires exactly 1 channel, got ", length(channels))
      }
      if (!is.numeric(coords) || length(coords) != 1) {
        stop("Threshold coords must be a single numeric value")
      }
      bounds <- list()
      bounds[[channels[1]]] <- c(coords, Inf)
      gate <- flowCore::rectangleGate(.gate = bounds)
      cat("Threshold gate on", channels[1], ">", coords, "\n")
      gate
    },

    "quadrant" = {
      if (length(channels) != 2) {
        stop("Quadrant gate requires exactly 2 channels, got ", length(channels))
      }
      if (!all(c("x", "y") %in% names(coords))) {
        stop("Quadrant coords must be a list with x and y (intersection point)")
      }
      x_cut <- coords$x
      y_cut <- coords$y

      # Q1: x-low, y-high (upper-left)
      q1_bounds <- list()
      q1_bounds[[channels[1]]] <- c(-Inf, x_cut)
      q1_bounds[[channels[2]]] <- c(y_cut, Inf)

      # Q2: x-high, y-high (upper-right)
      q2_bounds <- list()
      q2_bounds[[channels[1]]] <- c(x_cut, Inf)
      q2_bounds[[channels[2]]] <- c(y_cut, Inf)

      # Q3: x-low, y-low (lower-left)
      q3_bounds <- list()
      q3_bounds[[channels[1]]] <- c(-Inf, x_cut)
      q3_bounds[[channels[2]]] <- c(-Inf, y_cut)

      # Q4: x-high, y-low (lower-right)
      q4_bounds <- list()
      q4_bounds[[channels[1]]] <- c(x_cut, Inf)
      q4_bounds[[channels[2]]] <- c(-Inf, y_cut)

      gates <- list(
        Q1_xneg_ypos = flowCore::rectangleGate(.gate = q1_bounds),
        Q2_xpos_ypos = flowCore::rectangleGate(.gate = q2_bounds),
        Q3_xneg_yneg = flowCore::rectangleGate(.gate = q3_bounds),
        Q4_xpos_yneg = flowCore::rectangleGate(.gate = q4_bounds)
      )

      cat("Quadrant gate on", channels[1], "/", channels[2],
          "at intersection (", x_cut, ",", y_cut, ")\n")
      gates
    },

    stop("Invalid gate_type: '", gate_type,
         "'. Use 'rectangle', 'polygon', 'threshold', or 'quadrant'.")
  )

  result
}
