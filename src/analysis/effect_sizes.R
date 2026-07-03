#!/usr/bin/env Rscript
# =============================================================================
# analysis/effect_sizes.R — Effect Size Computation
# =============================================================================
# Computes Cohen's d, Cliff's delta, and Hedges' g for marker expression
# differences between two conditions. Returns effect sizes with confidence
# intervals and standard interpretations.
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

#' Compute effect sizes for marker expression between conditions
#'
#' Input from: Per-cluster CSVs from cluster_splitting or expression_matrices
#'             (analysis/cluster_splitting.R, analysis/expression_matrices.R)
#' Output to:  Standalone report, plot_heatmap() (visualization/heatmaps.R)
#'
#' @param expr_data data.table with marker columns + condition + sample_id columns
#' @param condition_col Column name for condition grouping
#' @param value_cols Marker column names to compute effect sizes for
#' @param methods Character vector: "cohens_d", "cliffs_delta", "hedges_g"
#' @return data.frame(marker, effect_size, ci_lower, ci_upper, method, interpretation)
compute_effect_sizes <- function(expr_data,
                                  condition_col,
                                  value_cols,
                                  methods = c("cohens_d", "cliffs_delta")) {
  require(data.table)

  # Validate input
  if (!condition_col %in% names(expr_data)) {
    stop("Condition column '", condition_col, "' not found in data")
  }
  conditions <- unique(expr_data[[condition_col]])
  conditions <- conditions[!is.na(conditions)]
  if (length(conditions) != 2) {
    stop("Exactly 2 conditions required, found: ", paste(conditions, collapse = ", "))
  }

  # Filter to valid value columns
  valid_cols <- intersect(value_cols, names(expr_data))
  if (length(valid_cols) == 0) {
    stop("No valid marker columns found")
  }
  cat("Computing effect sizes for", length(valid_cols), "markers\n")
  cat("Conditions:", conditions[1], "vs", conditions[2], "\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n")

  # Split data by condition
  grp1_data <- expr_data[expr_data[[condition_col]] == conditions[1], ]
  grp2_data <- expr_data[expr_data[[condition_col]] == conditions[2], ]
  cat("Group 1 (", conditions[1], "):", nrow(grp1_data), "cells\n")
  cat("Group 2 (", conditions[2], "):", nrow(grp2_data), "cells\n")

  results <- list()

  for (marker in valid_cols) {
    x <- grp1_data[[marker]]
    y <- grp2_data[[marker]]
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]

    if (length(x) < 2 || length(y) < 2) {
      cat("  Skipping", marker, "- insufficient data\n")
      next
    }

    for (method in methods) {
      result_row <- tryCatch({
        switch(method,
          "cohens_d"     = .compute_cohens_d(x, y),
          "cliffs_delta" = .compute_cliffs_delta(x, y),
          "hedges_g"     = .compute_hedges_g(x, y),
          stop("Unknown method: ", method)
        )
      }, error = function(e) {
        data.frame(
          effect_size = NA, ci_lower = NA, ci_upper = NA,
          method = method, interpretation = "error",
          stringsAsFactors = FALSE
        )
      })

      result_row$marker <- marker
      results[[length(results) + 1]] <- result_row
    }
  }

  if (length(results) == 0) {
    cat("No results produced.\n")
    return(data.frame())
  }

  results_df <- do.call(rbind, results)
  results_df <- results_df[, c("marker", "effect_size", "ci_lower", "ci_upper",
                                "method", "interpretation")]

  cat("\nEffect size computation complete:", nrow(results_df), "results\n")
  results_df
}

#' Compute effect sizes and save results to CSV
#'
#' Wrapper around compute_effect_sizes() that reads input from files and
#' writes output to a directory.
#'
#' Input from: combine_expression_matrices() (analysis/expression_matrices.R)
#' Output to:  Standalone CSV report
#'
#' @param input_csv Path to expression CSV with condition and sample columns
#' @param output_dir Output directory
#' @param condition_col Column name for condition
#' @param value_cols Marker column names (NULL = auto-detect)
#' @param methods Effect size methods to compute
#' @return data.frame with effect sizes
run_effect_sizes <- function(input_csv,
                              output_dir,
                              condition_col = "condition",
                              value_cols = NULL,
                              methods = c("cohens_d", "cliffs_delta")) {
  require(data.table)
  make_outdir(output_dir)

  expr_data <- fread(input_csv)
  cat("Loaded:", nrow(expr_data), "rows x", ncol(expr_data), "columns\n")

  # Auto-detect marker columns if not specified
  if (is.null(value_cols)) {
    col_classes <- classify_csv_columns(names(expr_data))
    value_cols <- col_classes$markers
    value_cols <- setdiff(value_cols, c(condition_col, "sample_id"))
  }

  results_df <- compute_effect_sizes(expr_data, condition_col, value_cols, methods)

  if (nrow(results_df) > 0) {
    fwrite(results_df, file.path(output_dir, "effect_sizes.csv"), sep = ";")
    cat("Results saved to:", file.path(output_dir, "effect_sizes.csv"), "\n")
  }

  results_df
}

# --- Internal: Cohen's d with noncentral t CI ---
.compute_cohens_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  mean_diff <- mean(x) - mean(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))

  if (pooled_sd == 0) {
    return(data.frame(
      effect_size = 0, ci_lower = 0, ci_upper = 0,
      method = "cohens_d", interpretation = "zero_variance",
      stringsAsFactors = FALSE
    ))
  }

  d <- mean_diff / pooled_sd

  # CI via noncentral t distribution
  se_d <- sqrt((nx + ny) / (nx * ny) + d^2 / (2 * (nx + ny)))
  ci <- tryCatch({
    ncp <- d * sqrt(nx * ny / (nx + ny))
    df <- nx + ny - 2
    ci_ncp <- c(
      suppressWarnings(qt(0.025, df, ncp)),
      suppressWarnings(qt(0.975, df, ncp))
    )
    ci_ncp / sqrt(nx * ny / (nx + ny))
  }, error = function(e) {
    c(d - 1.96 * se_d, d + 1.96 * se_d)
  })

  data.frame(
    effect_size = d,
    ci_lower = ci[1],
    ci_upper = ci[2],
    method = "cohens_d",
    interpretation = .interpret_cohens_d(abs(d)),
    stringsAsFactors = FALSE
  )
}

# --- Internal: Cliff's delta ---
.compute_cliffs_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)

  # Compute P(X > Y) - P(Y > X) via ranks
  # More efficient than all pairwise comparisons for large n
  if (nx * ny > 1e8) {
    # Subsample for very large datasets
    set.seed(42)
    x <- sample(x, min(nx, 10000))
    y <- sample(y, min(ny, 10000))
    nx <- length(x)
    ny <- length(y)
  }

  combined <- c(x, y)
  ranks <- rank(combined)
  rank_x <- ranks[1:nx]

  # Cliff's delta = (2 * mean(rank_x) - nx - ny - 1) / ny
  # Equivalently: (sum of dominance) / (nx * ny)
  delta <- (2 * sum(rank_x) - nx * (nx + 1)) / (nx * ny) - 1 +
           nx / ny  # correction
  # Simpler: delta = 2 * (mean(rank_x) - (nx + 1) / 2) / ny
  # But use the standard dominance matrix formulation:
  dominance <- outer(x, y, function(a, b) sign(a - b))
  delta <- mean(dominance)

  # Bootstrap CI (percentile method)
  n_boot <- 1000
  boot_deltas <- numeric(n_boot)
  set.seed(42)
  for (b in seq_len(n_boot)) {
    bx <- sample(x, nx, replace = TRUE)
    by <- sample(y, ny, replace = TRUE)
    # Fast dominance via rank-based approach for bootstrap
    bc <- c(bx, by)
    br <- rank(bc)
    brx <- br[1:nx]
    boot_deltas[b] <- 2 * mean(brx) / (nx + ny) - (nx + ny + 1) / (nx + ny) +
                       (ny - nx) / (nx + ny)
  }
  # Use simpler bootstrap with dominance for consistency
  boot_deltas <- replicate(n_boot, {
    bx <- sample(x, nx, replace = TRUE)
    by <- sample(y, ny, replace = TRUE)
    mean(sign(outer(bx, by, `-`)))
  })

  ci <- quantile(boot_deltas, c(0.025, 0.975))

  data.frame(
    effect_size = delta,
    ci_lower = ci[1],
    ci_upper = ci[2],
    method = "cliffs_delta",
    interpretation = .interpret_cliffs_delta(abs(delta)),
    stringsAsFactors = FALSE
  )
}

# --- Internal: Hedges' g (bias-corrected Cohen's d) ---
.compute_hedges_g <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  mean_diff <- mean(x) - mean(y)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))

  if (pooled_sd == 0) {
    return(data.frame(
      effect_size = 0, ci_lower = 0, ci_upper = 0,
      method = "hedges_g", interpretation = "zero_variance",
      stringsAsFactors = FALSE
    ))
  }

  d <- mean_diff / pooled_sd

  # Hedges' correction factor: J = 1 - 3 / (4 * df - 1)
  df <- nx + ny - 2
  j <- 1 - 3 / (4 * df - 1)
  g <- d * j

  # CI via corrected SE
  se_d <- sqrt((nx + ny) / (nx * ny) + d^2 / (2 * (nx + ny)))
  se_g <- se_d * j
  ci <- c(g - 1.96 * se_g, g + 1.96 * se_g)

  data.frame(
    effect_size = g,
    ci_lower = ci[1],
    ci_upper = ci[2],
    method = "hedges_g",
    interpretation = .interpret_cohens_d(abs(g)),
    stringsAsFactors = FALSE
  )
}

# --- Internal: Interpretation thresholds ---
.interpret_cohens_d <- function(abs_d) {
  if (is.na(abs_d)) return("NA")
  if (abs_d < 0.2) return("negligible")
  if (abs_d < 0.5) return("small")
  if (abs_d < 0.8) return("medium")
  return("large")
}

.interpret_cliffs_delta <- function(abs_delta) {
  if (is.na(abs_delta)) return("NA")
  if (abs_delta < 0.147) return("negligible")
  if (abs_delta < 0.33) return("small")
  if (abs_delta < 0.474) return("medium")
  return("large")
}
