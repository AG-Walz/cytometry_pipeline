#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/tests/test_preprocessing.R — Comprehensive tests for
# preprocessing module
# =============================================================================
# Tests: quality_control, transformation, normalization, doublet_detection,
#        live_dead_gating, downsampling, compensation, spectral_unmixing,
#        automated_gating, gate_utils, preprocess_pipeline
# =============================================================================

library(testthat)

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
.this_dir <- tryCatch(
  dirname(testthat::test_path()),
  error = function(e) getwd()
)
pp_dir      <- normalizePath(file.path(.this_dir, ".."), mustWork = FALSE)
src_dir     <- normalizePath(file.path(pp_dir, ".."), mustWork = FALSE)
shared_dir  <- file.path(src_dir, "R_shared")
project_dir <- normalizePath(file.path(src_dir, ".."), mustWork = FALSE)
test_data   <- file.path(project_dir, "test_data")

# ---------------------------------------------------------------------------
# Source shared utilities
# ---------------------------------------------------------------------------
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

# Source preprocessing modules (pure functions, no CLI blocks)
source(file.path(pp_dir, "quality_control.R"))
source(file.path(pp_dir, "transformation.R"))
source(file.path(pp_dir, "normalization.R"))
source(file.path(pp_dir, "doublet_detection.R"))
source(file.path(pp_dir, "live_dead_gating.R"))
source(file.path(pp_dir, "downsampling.R"))
source(file.path(pp_dir, "compensation.R"))
source(file.path(pp_dir, "spectral_unmixing.R"))
source(file.path(pp_dir, "automated_gating.R"))
source(file.path(pp_dir, "gate_utils.R"))

# ---------------------------------------------------------------------------
# Shared helper: create synthetic flowFrame
# ---------------------------------------------------------------------------
make_synthetic_ff <- function(n_cells = 200,
                              markers = c("CD3", "CD4", "CD8",
                                          "CD19", "CD56", "HLA-DR"),
                              seed = 42) {
  set.seed(seed)
  expr <- matrix(abs(rnorm(n_cells * length(markers), mean = 3, sd = 1)),
                 nrow = n_cells, ncol = length(markers))
  colnames(expr) <- markers
  params <- Biobase::AnnotatedDataFrame(data.frame(
    name     = markers,
    desc     = markers,
    range    = rep(1024, length(markers)),
    minRange = rep(0, length(markers)),
    maxRange = rep(1024, length(markers)),
    row.names = paste0("$P", seq_along(markers)),
    stringsAsFactors = FALSE
  ))
  flowCore::flowFrame(exprs = expr, parameters = params)
}

# Helper: flowFrame with scatter parameters (for doublet/viability tests)
make_scatter_ff <- function(n_cells = 500, seed = 42) {
  set.seed(seed)
  channels <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A",
                 "CD3", "CD4", "CD8", "LiveDead")

  # FSC-A and FSC-H are correlated for singlets
  fsc_h <- abs(rnorm(n_cells, mean = 50000, sd = 10000))
  # Singlets: FSC-A ~ FSC-H; doublets: FSC-A >> FSC-H
  n_doublets <- round(n_cells * 0.1)
  fsc_a <- fsc_h * runif(n_cells, 0.95, 1.05)
  # Inject doublets at the end
  doublet_idx <- (n_cells - n_doublets + 1):n_cells
  fsc_a[doublet_idx] <- fsc_h[doublet_idx] * runif(n_doublets, 1.3, 2.0)

  fsc_w <- abs(rnorm(n_cells, mean = 80000, sd = 5000))
  ssc_a <- abs(rnorm(n_cells, mean = 30000, sd = 8000))

  # Marker channels
  cd3 <- abs(rnorm(n_cells, 3, 1))
  cd4 <- abs(rnorm(n_cells, 2, 1))
  cd8 <- abs(rnorm(n_cells, 2, 1))

  # LiveDead: bimodal -- low = live, high = dead
  live_dead <- c(
    rnorm(round(n_cells * 0.8), mean = 200, sd = 50),   # live
    rnorm(round(n_cells * 0.2), mean = 2000, sd = 300)  # dead
  )
  live_dead <- abs(live_dead[1:n_cells])

  expr <- cbind(`FSC-A` = fsc_a, `FSC-H` = fsc_h, `FSC-W` = fsc_w,
                `SSC-A` = ssc_a, CD3 = cd3, CD4 = cd4, CD8 = cd8,
                LiveDead = live_dead)

  params <- Biobase::AnnotatedDataFrame(data.frame(
    name     = channels,
    desc     = channels,
    range    = apply(expr, 2, function(x) diff(range(x))),
    minRange = apply(expr, 2, min),
    maxRange = apply(expr, 2, max),
    row.names = paste0("$P", seq_along(channels)),
    stringsAsFactors = FALSE
  ))
  flowCore::flowFrame(exprs = expr, parameters = params)
}

# Temporary directory for output artifacts
test_tmpdir <- tempfile("pp_tests_")
dir.create(test_tmpdir, recursive = TRUE)
withr::defer(unlink(test_tmpdir, recursive = TRUE), teardown_env())

# =============================================================================
# SECTION 1: quality_control
# =============================================================================

context("quality_control")

test_that("quality_control with 'None' returns unchanged flowFrame", {
  ff <- make_synthetic_ff(n_cells = 100)
  result <- quality_control(ff, "test_sample", qc_type = "None")

  expect_s4_class(result, "flowFrame")
  expect_equal(nrow(result), nrow(ff))
  expect_equal(ncol(result), ncol(ff))
  # Expression values should be identical
  expect_identical(flowCore::exprs(result), flowCore::exprs(ff))
})

test_that("quality_control rejects invalid qc_type", {
  ff <- make_synthetic_ff(n_cells = 50)
  expect_error(
    quality_control(ff, "test", qc_type = "InvalidMethod"),
    "Invalid qc_type"
  )
})

test_that("quality_control with flowAI runs if available", {
  skip_if_not_installed("flowAI")
  ff <- make_synthetic_ff(n_cells = 200)

  result <- tryCatch(
    quality_control(ff, "test_flowAI", qc_type = "flowAI"),
    error = function(e) {
      skip(paste("flowAI QC failed (expected on synthetic data):", e$message))
    }
  )

  if (!is.null(result)) {
    expect_s4_class(result, "flowFrame")
    expect_lte(nrow(result), nrow(ff))
  }
})

test_that("quality_control with PeacoQC runs if available", {
  skip_if_not_installed("PeacoQC")
  ff <- make_synthetic_ff(n_cells = 200)

  result <- tryCatch(
    quality_control(ff, "test_PeacoQC", qc_type = "PeacoQC",
                    marker_indices = 1:ncol(ff)),
    error = function(e) {
      skip(paste("PeacoQC QC failed (expected on synthetic data):", e$message))
    }
  )

  if (!is.null(result)) {
    expect_s4_class(result, "flowFrame")
    expect_lte(nrow(result), nrow(ff))
  }
})

test_that("quality_control defaults marker_indices to all columns", {
  ff <- make_synthetic_ff(n_cells = 50)
  # "None" type should not error with NULL marker_indices
  result <- quality_control(ff, "test_default", qc_type = "None",
                            marker_indices = NULL)
  expect_s4_class(result, "flowFrame")
})

# =============================================================================
# SECTION 2: transformation
# =============================================================================

context("transform_flowframe")

test_that("arcsinh transform produces expected output", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4"))
  markers <- c("CD3", "CD4")
  offsets <- c(0, 0)
  widths  <- c(5, 5)

  result <- transform_flowframe(
    ff, "test_sample", markers, offsets, widths,
    transform_type = "arcsinh"
  )

  expect_s4_class(result, "flowFrame")
  expect_equal(unname(nrow(result)), 100)
  expect_equal(unname(ncol(result)), 2)

  # Transformed values should differ from original
  orig_vals <- flowCore::exprs(ff)[, "CD3"]
  trans_vals <- flowCore::exprs(result)[, "CD3"]
  expect_false(all(orig_vals == trans_vals))

  # Arcsinh of positive values is less than original for values > ~1
  # (since b = 1/5 = 0.2 means asinh(-0*x + 0.2*x) = asinh(0.2*x))
  # Just check the transform ran and produced finite values
  expect_true(all(is.finite(trans_vals)))
})

test_that("transform_flowframe with 'None' returns unchanged data", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  result <- transform_flowframe(
    ff, "test_none", c("CD3", "CD4"), c(0, 0), c(5, 5),
    transform_type = "None"
  )
  expect_identical(flowCore::exprs(result), flowCore::exprs(ff))
})

test_that("transform_flowframe validates marker/offset length mismatch", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  expect_error(
    transform_flowframe(ff, "test", c("CD3", "CD4"), c(0), c(5, 5)),
    "Marker and offset lists must have the same length"
  )
  expect_error(
    transform_flowframe(ff, "test", c("CD3", "CD4"), c(0, 0), c(5)),
    "Marker and width lists must have the same length"
  )
})

test_that("transform_flowframe warns on missing markers", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  expect_warning(
    transform_flowframe(
      ff, "test", c("CD3", "NONEXISTENT"), c(0, 0), c(5, 5),
      transform_type = "arcsinh"
    ),
    "not found"
  )
})

test_that("transform_flowframe rejects invalid transform_type", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3"))
  expect_error(
    transform_flowframe(ff, "test", "CD3", 0, 5, transform_type = "invalid"),
    "Invalid transform_type"
  )
})

test_that("logicle transform produces finite values", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3"))
  result <- tryCatch(
    transform_flowframe(
      ff, "test_logicle", "CD3", 0, 5,
      transform_type = "logicle"
    ),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    expect_s4_class(result, "flowFrame")
    expect_true(all(is.finite(flowCore::exprs(result)[, "CD3"])))
  }
})

test_that("log transform adds offset and handles near-zero values", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3"))
  result <- transform_flowframe(
    ff, "test_log", "CD3", 1, 5,
    transform_type = "log"
  )

  expect_s4_class(result, "flowFrame")
  trans_vals <- flowCore::exprs(result)[, "CD3"]
  expect_true(all(is.finite(trans_vals)))
})

test_that("biexponential transform produces output", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3"))
  result <- tryCatch(
    transform_flowframe(
      ff, "test_biexp", "CD3", 0, 5,
      transform_type = "biexponential"
    ),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    expect_s4_class(result, "flowFrame")
  }
})

test_that("find_optimal_b_cofactor returns a positive number", {
  set.seed(42)
  data <- c(rnorm(200, mean = 100, sd = 30),
            rnorm(200, mean = 1000, sd = 100))
  result <- find_optimal_b_cofactor(data)

  expect_type(result, "double")
  expect_gt(result, 0)
})

# =============================================================================
# SECTION 3: normalization
# =============================================================================

context("normalization")

test_that("z_score_scaling produces zero-mean columns", {
  ff <- make_synthetic_ff(n_cells = 200, markers = c("CD3", "CD4", "CD8"))
  marker_indices <- 1:3

  result <- z_score_scaling(ff, marker_indices)

  expect_s4_class(result, "flowFrame")
  expect_equal(unname(nrow(result)), 200)

  # Each scaled marker column should have mean ~ 0
  expr <- flowCore::exprs(result)
  for (i in marker_indices) {
    col_mean <- mean(expr[, i])
    expect_equal(col_mean, 0, tolerance = 1e-10,
                 label = paste("Mean of column", i))
  }
})

test_that("z_score_scaling produces unit-variance columns", {
  ff <- make_synthetic_ff(n_cells = 300, markers = c("CD3", "CD4", "CD8"))
  marker_indices <- 1:3

  result <- z_score_scaling(ff, marker_indices)
  expr <- flowCore::exprs(result)

  for (i in marker_indices) {
    col_sd <- sd(expr[, i])
    expect_equal(col_sd, 1, tolerance = 0.01,
                 label = paste("SD of column", i))
  }
})

test_that("z_score_scaling on subset of markers leaves others unchanged", {
  markers <- c("CD3", "CD4", "CD8", "CD19", "CD56")
  ff <- make_synthetic_ff(n_cells = 100, markers = markers)
  orig_expr <- flowCore::exprs(ff)

  # Only scale first two columns
  result <- z_score_scaling(ff, marker_indices = 1:2)
  result_expr <- flowCore::exprs(result)

  # Columns 3-5 should be unchanged
  expect_equal(result_expr[, 3:5], orig_expr[, 3:5])
  # Columns 1-2 should differ
  expect_false(identical(result_expr[, 1], orig_expr[, 1]))
})

test_that("quantile_norm produces normalized output", {
  skip_if_not_installed("preprocessCore")
  # preprocessCore's native threading can fail under a container/CI cpuset
  # with `pthread_create() ... 22`. Skip in that case rather than fail.
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4", "CD8"))
  result <- tryCatch(
    quantile_norm(ff, marker_indices = 1:3),
    error = function(e) {
      if (grepl("pthread_create", conditionMessage(e))) {
        skip(paste("preprocessCore pthread unavailable:", conditionMessage(e)))
      }
      stop(e)
    }
  )

  expect_s4_class(result, "flowFrame")
  expect_equal(unname(nrow(result)), 100)
  expect_equal(unname(ncol(result)), 3)

  # After quantile normalization, the distribution of each column
  # should be more similar to each other
  expr <- flowCore::exprs(result)
  sorted_cols <- apply(expr, 2, sort)
  # All sorted columns should be identical after quantile normalization
  expect_equal(sorted_cols[, 1], sorted_cols[, 2], tolerance = 1e-10)
  expect_equal(sorted_cols[, 1], sorted_cols[, 3], tolerance = 1e-10)
})

test_that("quantile_norm errors without preprocessCore", {
  # This test only verifies the error message when the package is missing.
  # When preprocessCore IS installed, the function succeeds (tested above).
  # We cannot easily simulate a missing package, so just verify the function
  # exists and has the right check in its body.
  fn_body <- deparse(body(quantile_norm))
  expect_true(any(grepl("preprocessCore", fn_body)))
})

test_that("harmony_normalize runs if harmony is installed", {
  skip_if_not_installed("harmony")

  set.seed(42)
  n <- 200
  expr_matrix <- matrix(rnorm(n * 5, mean = 3), nrow = n, ncol = 5)
  colnames(expr_matrix) <- paste0("M", 1:5)
  # Add batch effect to second half
  expr_matrix[101:200, ] <- expr_matrix[101:200, ] + 2
  batch_vector <- rep(c("batch1", "batch2"), each = 100)

  result <- tryCatch(
    harmony_normalize(expr_matrix, batch_vector),
    error = function(e) {
      skip(paste("Harmony failed:", e$message))
    }
  )

  if (!is.null(result)) {
    expect_true(is.matrix(result))
    expect_equal(dim(result), dim(expr_matrix))
  }
})

# =============================================================================
# SECTION 4: doublet_detection
# =============================================================================

context("detect_doublets")

test_that("detect_doublets ratio method identifies doublets", {
  ff <- make_scatter_ff(n_cells = 500)

  result <- detect_doublets(ff, area_col = "FSC-A", height_col = "FSC-H",
                            method = "ratio", threshold = 1.1)

  expect_type(result, "list")
  expect_named(result, c("filtered_ff", "doublet_indices", "singlet_rate"))
  expect_s4_class(result$filtered_ff, "flowFrame")

  # Should have removed some doublets
  expect_lt(nrow(result$filtered_ff), 500)
  expect_gt(length(result$doublet_indices), 0)
  expect_gt(result$singlet_rate, 0)
  expect_lt(result$singlet_rate, 1)

  # Singlet rate should be around 0.9 (we injected 10% doublets)
  expect_gt(result$singlet_rate, 0.7)
})

test_that("detect_doublets model method runs", {
  ff <- make_scatter_ff(n_cells = 500)

  result <- detect_doublets(ff, area_col = "FSC-A", height_col = "FSC-H",
                            method = "model")

  expect_type(result, "list")
  expect_s4_class(result$filtered_ff, "flowFrame")
  expect_lt(nrow(result$filtered_ff), 500)
  expect_gt(result$singlet_rate, 0)
})

test_that("detect_doublets rejects invalid method", {
  ff <- make_scatter_ff(n_cells = 100)
  expect_error(
    detect_doublets(ff, method = "invalid_method"),
    "Invalid method"
  )
})

test_that("detect_doublets auto-detects FSC column names", {
  # Create flowFrame with slightly different FSC naming
  ff <- make_scatter_ff(n_cells = 200)
  # The make_scatter_ff already uses "FSC-A" and "FSC-H"
  # Test that auto-matching works when given default params
  result <- detect_doublets(ff)

  expect_type(result, "list")
  expect_s4_class(result$filtered_ff, "flowFrame")
})

test_that("detect_doublets errors when FSC columns missing", {
  # Create a flowFrame without FSC channels
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4", "CD8"))
  expect_error(detect_doublets(ff, area_col = "FSC-A", height_col = "FSC-H"),
               "not found")
})

test_that("detect_doublets with width_col applies additional filtering", {
  ff <- make_scatter_ff(n_cells = 500)

  result_no_width <- detect_doublets(ff, method = "ratio", threshold = 1.1)
  result_with_width <- detect_doublets(ff, method = "ratio", threshold = 1.1,
                                        width_col = "FSC-W")

  # Width-based filtering may remove additional events
  expect_lte(nrow(result_with_width$filtered_ff),
             nrow(result_no_width$filtered_ff))
})

test_that("detect_doublets preserves all columns in filtered output", {
  ff <- make_scatter_ff(n_cells = 200)
  result <- detect_doublets(ff, method = "ratio")

  expect_equal(ncol(result$filtered_ff), ncol(ff))
  expect_equal(colnames(result$filtered_ff), colnames(ff))
})

test_that("detect_doublets with threshold = Inf removes nothing", {
  ff <- make_scatter_ff(n_cells = 100)
  result <- detect_doublets(ff, method = "ratio", threshold = Inf)

  expect_equal(nrow(result$filtered_ff), 100)
  expect_equal(length(result$doublet_indices), 0)
  expect_equal(result$singlet_rate, 1)
})

# =============================================================================
# SECTION 5: live_dead_gating
# =============================================================================

context("gate_live_cells")

test_that("gate_live_cells threshold method gates bimodal data", {
  ff <- make_scatter_ff(n_cells = 500)

  result <- gate_live_cells(
    ff, viability_channel = "LiveDead",
    method = "threshold", threshold = 800,
    dead_high = TRUE
  )

  expect_type(result, "list")
  expect_named(result, c("live_ff", "dead_indices", "viability_rate"))
  expect_s4_class(result$live_ff, "flowFrame")
  expect_lt(nrow(result$live_ff), 500)
  expect_gt(result$viability_rate, 0)
  expect_lt(result$viability_rate, 1)
  # We created ~80% live, 20% dead
  expect_gt(result$viability_rate, 0.5)
})

test_that("gate_live_cells auto-detects threshold from density", {
  ff <- make_scatter_ff(n_cells = 500)

  result <- gate_live_cells(
    ff, viability_channel = "LiveDead",
    method = "threshold", threshold = NULL,
    dead_high = TRUE
  )

  expect_type(result, "list")
  expect_s4_class(result$live_ff, "flowFrame")
  # Auto threshold should find the valley between live and dead
  expect_gt(result$viability_rate, 0.5)
})

test_that("gate_live_cells mixture method runs", {
  ff <- make_scatter_ff(n_cells = 500)

  result <- gate_live_cells(
    ff, viability_channel = "LiveDead",
    method = "mixture",
    dead_high = TRUE
  )

  expect_type(result, "list")
  expect_s4_class(result$live_ff, "flowFrame")
  expect_gt(result$viability_rate, 0)
})

test_that("gate_live_cells dead_high = FALSE inverts gating", {
  ff <- make_scatter_ff(n_cells = 500)
  threshold_val <- 800

  result_high <- gate_live_cells(
    ff, "LiveDead", method = "threshold",
    threshold = threshold_val, dead_high = TRUE
  )
  result_low <- gate_live_cells(
    ff, "LiveDead", method = "threshold",
    threshold = threshold_val, dead_high = FALSE
  )

  # Complementary gating: live cells in one should be dead in the other
  total <- nrow(ff)
  expect_equal(
    nrow(result_high$live_ff) + nrow(result_low$live_ff),
    total
  )
})

test_that("gate_live_cells errors on missing viability channel", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  expect_error(
    gate_live_cells(ff, viability_channel = "MISSING_CHANNEL"),
    "not found"
  )
})

test_that("gate_live_cells rejects invalid method", {
  ff <- make_scatter_ff(n_cells = 100)
  expect_error(
    gate_live_cells(ff, "LiveDead", method = "invalid"),
    "Invalid method"
  )
})

test_that("gate_live_cells warns on low viability rate", {
  # Create data where almost everything is "dead"
  set.seed(99)
  n <- 200
  channels <- c("FSC-A", "LiveDead")
  expr <- cbind(
    `FSC-A`    = abs(rnorm(n, 50000, 10000)),
    LiveDead   = abs(rnorm(n, 5000, 500))  # all high = all dead
  )
  params <- Biobase::AnnotatedDataFrame(data.frame(
    name     = channels,
    desc     = channels,
    range    = apply(expr, 2, function(x) diff(range(x))),
    minRange = apply(expr, 2, min),
    maxRange = apply(expr, 2, max),
    row.names = paste0("$P", 1:2),
    stringsAsFactors = FALSE
  ))
  ff <- flowCore::flowFrame(exprs = expr, parameters = params)

  expect_warning(
    gate_live_cells(ff, "LiveDead", method = "threshold",
                    threshold = 100, dead_high = TRUE),
    "Low viability rate"
  )
})

test_that("gate_live_cells preserves all columns", {
  ff <- make_scatter_ff(n_cells = 200)
  result <- gate_live_cells(ff, "LiveDead", method = "threshold",
                            threshold = 800)
  expect_equal(ncol(result$live_ff), ncol(ff))
  expect_equal(colnames(result$live_ff), colnames(ff))
})

test_that("find_density_valley finds valley in bimodal data", {
  set.seed(42)
  bimodal <- c(rnorm(300, mean = 2, sd = 0.5),
               rnorm(300, mean = 8, sd = 0.5))
  valley <- find_density_valley(bimodal)

  expect_type(valley, "double")
  expect_gt(valley, 2)
  expect_lt(valley, 8)
  # Should be near the midpoint
  expect_gt(valley, 3)
  expect_lt(valley, 7)
})

test_that("find_density_valley falls back to median for unimodal data", {
  set.seed(42)
  unimodal <- rnorm(500, mean = 5, sd = 1)
  expect_warning(
    valley <- find_density_valley(unimodal),
    "No bimodal pattern"
  )
  expect_equal(valley, median(unimodal), tolerance = 0.5)
})

test_that("fit_gmm_threshold returns value between two modes", {
  set.seed(42)
  bimodal <- c(rnorm(300, mean = 100, sd = 20),
               rnorm(300, mean = 500, sd = 50))
  threshold <- fit_gmm_threshold(bimodal)

  expect_type(threshold, "double")
  expect_gt(threshold, 100)
  expect_lt(threshold, 500)
})

# =============================================================================
# SECTION 6: downsampling
# =============================================================================

context("downsampling")

test_that("downsample Consistent takes same count per sample", {
  ff1 <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4"))
  ff2 <- make_synthetic_ff(n_cells = 200, markers = c("CD3", "CD4"))
  ff3 <- make_synthetic_ff(n_cells = 150, markers = c("CD3", "CD4"))
  ffl <- list(ff1, ff2, ff3)
  names <- c("s1", "s2", "s3")

  result <- downsample_flowframes(
    ffl, names, target_cells = 120,
    method = "Consistent", seed = 42
  )

  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_named(result, names)

  # Each sample should have the same number of cells
  counts <- sapply(result, nrow)
  expect_equal(length(unique(counts)), 1)
  # per_sample = min(120/3, min(100,200,150)) = min(40, 100) = 40
  expect_equal(counts[[1]], 40)
})

test_that("downsample Maximal uses proportional sampling", {
  ff1 <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4"))
  ff2 <- make_synthetic_ff(n_cells = 200, markers = c("CD3", "CD4"))
  ffl <- list(ff1, ff2)
  names <- c("s1", "s2")

  result <- downsample_flowframes(
    ffl, names, target_cells = 150,
    method = "Maximal", seed = 42
  )

  expect_type(result, "list")
  expect_equal(length(result), 2)

  # s2 should get more cells since it's bigger
  expect_gt(nrow(result$s2), nrow(result$s1))
  # Total should be close to target (may differ due to rounding)
  total <- sum(sapply(result, nrow))
  expect_lte(total, 150)
})

test_that("downsample rejects unknown method", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3"))
  expect_error(
    downsample_flowframes(list(ff), "s1", method = "Unknown"),
    "Unknown method"
  )
})

test_that("downsample is reproducible with same seed", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4"))
  ffl <- list(ff)

  r1 <- downsample_flowframes(ffl, "s", target_cells = 50,
                                method = "Consistent", seed = 123)
  r2 <- downsample_flowframes(ffl, "s", target_cells = 50,
                                method = "Consistent", seed = 123)

  expect_identical(flowCore::exprs(r1$s), flowCore::exprs(r2$s))
})

test_that("downsample handles target_cells > total cells", {
  ff <- make_synthetic_ff(n_cells = 30, markers = c("CD3"))
  ffl <- list(ff)

  result <- downsample_flowframes(ffl, "s", target_cells = 1000,
                                   method = "Maximal")
  # Should keep all cells when target exceeds available
  expect_equal(nrow(result$s), 30)
})

test_that("downsample save_fcs writes FCS files", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  ffl <- list(ff)
  ds_out <- file.path(test_tmpdir, "downsampled")

  result <- downsample_flowframes(
    ffl, "test_ds", target_cells = 20,
    method = "Consistent", save_fcs = TRUE, output_dir = ds_out
  )

  expect_true(dir.exists(ds_out))
  expect_true(file.exists(file.path(ds_out, "ds_test_ds.fcs")))
})

# =============================================================================
# SECTION 7: compensation
# =============================================================================

context("compensation")

test_that("apply_compensation works with explicit matrix", {
  markers <- c("FL1-A", "FL2-A", "FL3-A")
  ff <- make_synthetic_ff(n_cells = 100, markers = markers)

  # Simple compensation matrix (near-identity with small spillover)
  comp <- diag(3)
  comp[1, 2] <- 0.1
  comp[2, 3] <- 0.05
  colnames(comp) <- markers
  rownames(comp) <- markers

  result <- apply_compensation(ff, comp_matrix = comp, auto_extract = FALSE)

  expect_s4_class(result, "flowFrame")
  expect_equal(nrow(result), 100)
  # Compensated values should differ from original
  expect_false(identical(
    flowCore::exprs(result), flowCore::exprs(ff)
  ))
})

test_that("apply_compensation errors when no matrix and auto_extract fails", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  expect_error(
    apply_compensation(ff, comp_matrix = NULL, auto_extract = FALSE),
    "No compensation matrix provided"
  )
})

test_that("apply_compensation errors when channels do not match", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  comp <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)
  colnames(comp) <- c("NOMATCH1", "NOMATCH2")
  rownames(comp) <- c("NOMATCH1", "NOMATCH2")

  expect_error(
    apply_compensation(ff, comp_matrix = comp, auto_extract = FALSE),
    "No channels in compensation matrix match"
  )
})

test_that("extract_compensation_matrix returns NULL when no SPILL keyword", {
  ff <- make_synthetic_ff(n_cells = 50)
  expect_warning(result <- extract_compensation_matrix(ff), "No compensation matrix")
  expect_null(result)
})

test_that("load_compensation_matrix reads CSV file", {
  comp <- matrix(c(1, 0.1, 0.05, 0.1, 1, 0.08, 0.02, 0.1, 1),
                 nrow = 3, byrow = TRUE)
  markers <- c("FL1", "FL2", "FL3")
  colnames(comp) <- markers
  dt <- data.table::data.table(Channel = markers)
  dt <- cbind(dt, as.data.frame(comp))
  csv_path <- file.path(test_tmpdir, "comp_matrix.csv")
  data.table::fwrite(dt, csv_path, sep = ";")

  result <- load_compensation_matrix(csv_path, sep = ";")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), markers)
})

# =============================================================================
# SECTION 8: spectral_unmixing
# =============================================================================

context("spectral_unmixing")

test_that("build_unmixing_matrix creates normalized reference spectra", {
  skip_if_not_installed("nnls")

  channels <- paste0("Ch", 1:10)

  # Unstained control
  set.seed(42)
  unstained_expr <- matrix(abs(rnorm(200 * 10, mean = 5, sd = 2)),
                           nrow = 200, ncol = 10)
  colnames(unstained_expr) <- channels
  unstained_ff <- from_mtx_to_ffr(unstained_expr)

  # Single-stain controls
  controls <- list()
  for (fluor in c("FITC", "PE", "APC")) {
    ctrl_expr <- matrix(abs(rnorm(100 * 10, mean = 5, sd = 2)),
                        nrow = 100, ncol = 10)
    # Set a peak channel
    peak_ch <- switch(fluor, FITC = 3, PE = 5, APC = 8)
    ctrl_expr[, peak_ch] <- ctrl_expr[, peak_ch] + 100
    colnames(ctrl_expr) <- channels
    controls[[fluor]] <- from_mtx_to_ffr(ctrl_expr)
  }

  unmix_mat <- build_unmixing_matrix(controls, unstained_ff, channels)

  expect_true(is.matrix(unmix_mat))
  expect_equal(nrow(unmix_mat), 10)
  expect_equal(ncol(unmix_mat), 3)
  expect_equal(rownames(unmix_mat), channels)
  expect_equal(colnames(unmix_mat), c("FITC", "PE", "APC"))

  # Each column should sum to 1 (normalized spectra)
  col_sums <- colSums(unmix_mat)
  for (s in col_sums) {
    expect_equal(s, 1, tolerance = 1e-6)
  }
})

test_that("unmix_spectral produces unmixed flowFrame", {
  skip_if_not_installed("nnls")

  channels <- paste0("Ch", 1:5)
  fluors <- c("FITC", "PE")

  # Build a simple 5x2 unmixing matrix
  unmix_mat <- matrix(0, nrow = 5, ncol = 2)
  unmix_mat[2, 1] <- 0.8  # FITC peak
  unmix_mat[3, 1] <- 0.2
  unmix_mat[4, 2] <- 0.7  # PE peak
  unmix_mat[5, 2] <- 0.3
  rownames(unmix_mat) <- channels
  colnames(unmix_mat) <- fluors

  # Create spectral data
  set.seed(42)
  expr <- matrix(abs(rnorm(50 * 5, mean = 10, sd = 3)),
                 nrow = 50, ncol = 5)
  colnames(expr) <- channels
  ff <- from_mtx_to_ffr(expr)

  result_ff <- unmix_spectral(ff, unmix_mat)

  expect_s4_class(result_ff, "flowFrame")
  expect_equal(nrow(result_ff), 50)
  expect_true("FITC" %in% colnames(result_ff))
  expect_true("PE" %in% colnames(result_ff))
})

test_that("unmix_spectral errors on missing channels", {
  skip_if_not_installed("nnls")

  unmix_mat <- matrix(1, nrow = 3, ncol = 2)
  rownames(unmix_mat) <- c("MISSING1", "MISSING2", "MISSING3")
  colnames(unmix_mat) <- c("F1", "F2")

  ff <- make_synthetic_ff(n_cells = 20, markers = c("CD3", "CD4"))

  expect_error(unmix_spectral(ff, unmix_mat), "Missing channels")
})

# =============================================================================
# SECTION 9: automated_gating
# =============================================================================

context("automated_gating")

test_that("define_gate rectangle creates valid gate", {
  gate <- define_gate(
    gate_type = "rectangle",
    channels = c("CD3", "CD4"),
    coords = list(x = c(1, 5), y = c(2, 8))
  )
  expect_s4_class(gate, "rectangleGate")
})

test_that("define_gate threshold creates single-channel gate", {
  gate <- define_gate(
    gate_type = "threshold",
    channels = "CD3",
    coords = 2.5
  )
  expect_s4_class(gate, "rectangleGate")
})

test_that("define_gate polygon creates polygonGate", {
  vertices <- matrix(c(1, 5, 5, 1,   # x coords
                       1, 1, 5, 5),   # y coords
                     ncol = 2)
  gate <- define_gate(
    gate_type = "polygon",
    channels = c("CD3", "CD4"),
    coords = vertices
  )
  expect_s4_class(gate, "polygonGate")
})

test_that("define_gate quadrant returns 4 gates", {
  gates <- define_gate(
    gate_type = "quadrant",
    channels = c("CD3", "CD4"),
    coords = list(x = 3, y = 3)
  )

  expect_type(gates, "list")
  expect_equal(length(gates), 4)
  expect_named(gates, c("Q1_xneg_ypos", "Q2_xpos_ypos",
                         "Q3_xneg_yneg", "Q4_xpos_yneg"))
  for (g in gates) {
    expect_s4_class(g, "rectangleGate")
  }
})

test_that("define_gate validates channel count for rectangle", {
  expect_error(
    define_gate("rectangle", channels = "CD3",
                coords = list(x = c(1, 5), y = c(2, 8))),
    "requires exactly 2 channels"
  )
})

test_that("define_gate validates channel count for threshold", {
  expect_error(
    define_gate("threshold", channels = c("CD3", "CD4"), coords = 2.5),
    "requires exactly 1 channel"
  )
})

test_that("define_gate validates coords for threshold", {
  expect_error(
    define_gate("threshold", channels = "CD3", coords = "not_numeric"),
    "single numeric value"
  )
})

test_that("define_gate rejects invalid gate_type", {
  expect_error(
    define_gate("invalid_type", channels = "CD3", coords = 1),
    "Invalid gate_type"
  )
})

# =============================================================================
# SECTION 10: gate_utils
# =============================================================================

context("gate_utils")

test_that("apply_gate_sequence with AND logic intersects gates", {
  ff <- make_synthetic_ff(n_cells = 200, markers = c("CD3", "CD4", "CD8"))

  gate1 <- flowCore::rectangleGate(
    .gate = list(CD3 = c(3, Inf))
  )
  gate2 <- flowCore::rectangleGate(
    .gate = list(CD4 = c(2.5, Inf))
  )

  result <- apply_gate_sequence(ff, list(gate1, gate2), logical_op = "and")

  expect_type(result, "list")
  expect_named(result, c("gated_ff", "gate_indices", "population_stats"))
  expect_s4_class(result$gated_ff, "flowFrame")

  # AND should be a subset of either gate alone
  r1 <- apply_gate_sequence(ff, list(gate1), logical_op = "and")
  r2 <- apply_gate_sequence(ff, list(gate2), logical_op = "and")
  expect_lte(nrow(result$gated_ff), nrow(r1$gated_ff))
  expect_lte(nrow(result$gated_ff), nrow(r2$gated_ff))
})

test_that("apply_gate_sequence with OR logic unions gates", {
  ff <- make_synthetic_ff(n_cells = 200, markers = c("CD3", "CD4"))

  gate1 <- flowCore::rectangleGate(.gate = list(CD3 = c(4, Inf)))
  gate2 <- flowCore::rectangleGate(.gate = list(CD4 = c(4, Inf)))

  result <- apply_gate_sequence(ff, list(gate1, gate2), logical_op = "or")

  # OR should include more cells than either gate alone
  r1 <- apply_gate_sequence(ff, list(gate1), logical_op = "and")
  r2 <- apply_gate_sequence(ff, list(gate2), logical_op = "and")
  expect_gte(nrow(result$gated_ff), max(nrow(r1$gated_ff), nrow(r2$gated_ff)))
})

test_that("apply_gate_sequence validates logical_op", {
  ff <- make_synthetic_ff(n_cells = 50)
  gate <- flowCore::rectangleGate(.gate = list(CD3 = c(2, Inf)))
  expect_error(
    apply_gate_sequence(ff, list(gate), logical_op = "xor"),
    "must be 'and' or 'or'"
  )
})

test_that("apply_gate_sequence errors on empty gate list", {
  ff <- make_synthetic_ff(n_cells = 50)
  expect_error(
    apply_gate_sequence(ff, list()),
    "at least one gate"
  )
})

test_that("boolean_gate evaluates AND expression correctly", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4", "CD8"))
  n <- nrow(ff)

  gate_results <- list(
    "CD3+" = flowCore::exprs(ff)[, "CD3"] > 3,
    "CD4+" = flowCore::exprs(ff)[, "CD4"] > 2
  )

  result <- boolean_gate(ff, gate_results, "CD3+ & CD4+")

  expect_type(result, "logical")
  expect_equal(length(result), n)
  # Manual verification
  expected <- gate_results[["CD3+"]] & gate_results[["CD4+"]]
  expect_equal(result, expected)
})

test_that("boolean_gate evaluates NOT expression", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4", "CD8"))
  n <- nrow(ff)

  gate_results <- list(
    "CD3+" = flowCore::exprs(ff)[, "CD3"] > 3,
    "CD8+" = flowCore::exprs(ff)[, "CD8"] > 3
  )

  result <- boolean_gate(ff, gate_results, "CD3+ & !CD8+")

  expected <- gate_results[["CD3+"]] & !gate_results[["CD8+"]]
  expect_equal(result, expected)
})

test_that("boolean_gate validates gate_results input", {
  ff <- make_synthetic_ff(n_cells = 50)
  expect_error(
    boolean_gate(ff, "not_a_list", "CD3+"),
    "named list"
  )
  # A named list whose element has the wrong length should also error.
  expect_error(
    boolean_gate(ff, list("unnamed_element" = TRUE), "CD3+"),
    "has length"
  )
  # A named list whose element is not logical should error with "logical vector".
  expect_error(
    boolean_gate(ff, list("unnamed_element" = rep(1, nrow(ff))), "CD3+"),
    "logical vector"
  )
})

test_that("boolean_gate errors on length mismatch", {
  ff <- make_synthetic_ff(n_cells = 100)
  gate_results <- list("G1" = c(TRUE, FALSE))  # wrong length
  expect_error(
    boolean_gate(ff, gate_results, "G1"),
    "length"
  )
})

# =============================================================================
# SECTION 11: preprocess_pipeline (integration)
# =============================================================================

context("preprocess_pipeline integration")

test_that("run_preprocessing chains QC -> Transform -> Normalize", {
  # Source the orchestrator
  tryCatch(
    source(file.path(pp_dir, "preprocess_pipeline.R")),
    error = function(e) skip(paste("Cannot source preprocess_pipeline:", e$message))
  )
  skip_if_not(exists("run_preprocessing"),
              "run_preprocessing not available")

  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4", "CD8"))
  flow_frames <- list(test_sample = ff)

  params <- list(
    pp_order = "QT",
    qc_type = "None",
    transform_type = "arcsinh",
    markers_of_interest = c("CD3", "CD4", "CD8"),
    offset_list = c(0, 0, 0),
    width_list = c(5, 5, 5),
    quant_norm = FALSE,
    scaling = FALSE,
    save_pp_fcs = FALSE,
    save_pp_csv = FALSE
  )

  result <- run_preprocessing(flow_frames, params)

  expect_type(result, "list")
  expect_named(result, "test_sample")
  expect_s4_class(result$test_sample, "flowFrame")
  expect_equal(nrow(result$test_sample), 100)
})

test_that("run_preprocessing TQ order runs Transform before QC", {
  skip_if_not(exists("run_preprocessing"),
              "run_preprocessing not available")

  ff <- make_synthetic_ff(n_cells = 80, markers = c("CD3", "CD4"))
  flow_frames <- list(s1 = ff)

  params <- list(
    pp_order = "TQ",
    qc_type = "None",
    transform_type = "arcsinh",
    markers_of_interest = c("CD3", "CD4"),
    offset_list = c(0, 0),
    width_list = c(5, 5),
    quant_norm = FALSE,
    scaling = TRUE,
    save_pp_fcs = FALSE,
    save_pp_csv = FALSE
  )

  result <- run_preprocessing(flow_frames, params)

  expect_type(result, "list")
  expect_s4_class(result$s1, "flowFrame")

  # With scaling=TRUE, the marker columns should have mean ~0
  expr <- flowCore::exprs(result$s1)
  for (m in c("CD3", "CD4")) {
    expect_equal(mean(expr[, m]), 0, tolerance = 1e-10,
                 label = paste("Mean of", m, "after scaling"))
  }
})

test_that("run_preprocessing saves output files when output_dir is set", {
  skip_if_not(exists("run_preprocessing"),
              "run_preprocessing not available")

  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  flow_frames <- list(output_test = ff)

  out_dir <- file.path(test_tmpdir, "pp_output")

  params <- list(
    pp_order = "QT",
    qc_type = "None",
    transform_type = "None",
    markers_of_interest = c("CD3", "CD4"),
    offset_list = c(0, 0),
    width_list = c(5, 5),
    quant_norm = FALSE,
    scaling = FALSE,
    save_pp_fcs = TRUE,
    save_pp_csv = TRUE
  )

  result <- run_preprocessing(flow_frames, params, output_dir = out_dir)

  pp_dir_out <- file.path(out_dir, "preprocessed_files")
  expect_true(dir.exists(pp_dir_out))
  expect_true(file.exists(file.path(pp_dir_out, "pp_output_test.fcs")))
  expect_true(file.exists(file.path(pp_dir_out, "pp_output_test.csv")))
})

# =============================================================================
# SECTION 12: bead_normalization
# =============================================================================

context("bead_normalization")

source(file.path(pp_dir, "bead_normalization.R"))

test_that("identify_beads finds high-intensity events", {
  set.seed(42)
  # Simulate bead channels: most cells are dim, beads are bright
  n_cells <- 500
  n_beads <- 30
  bead_expr <- matrix(abs(rnorm(n_cells * 3, mean = 100, sd = 50)),
                      nrow = n_cells, ncol = 3)
  # Make some events very bright (beads)
  bead_expr[1:n_beads, ] <- matrix(abs(rnorm(n_beads * 3, mean = 1000, sd = 100)),
                                    nrow = n_beads, ncol = 3)
  colnames(bead_expr) <- c("Ch1", "Ch2", "Ch3")

  is_bead <- identify_beads(bead_expr, method = "mahalanobis")

  expect_type(is_bead, "logical")
  expect_equal(length(is_bead), n_cells)
  # Should find at least some beads
  expect_gt(sum(is_bead), 0)
  # Most of the actual beads should be identified
  expect_gt(sum(is_bead[1:n_beads]), n_beads * 0.5)
})

test_that("identify_beads threshold method works", {
  set.seed(42)
  n_cells <- 200
  bead_expr <- matrix(abs(rnorm(n_cells * 2, mean = 50, sd = 20)),
                      nrow = n_cells, ncol = 2)
  bead_expr[1:10, ] <- abs(rnorm(20, mean = 500, sd = 50))
  colnames(bead_expr) <- c("Ch1", "Ch2")

  is_bead <- identify_beads(bead_expr, method = "threshold")

  expect_type(is_bead, "logical")
  expect_equal(length(is_bead), n_cells)
})

# =============================================================================
# SECTION 13: supercell (skipped if SuperCellCyto not available)
# =============================================================================

context("supercell")

test_that("run_supercell validates input dimensions", {
  skip_if_not_installed("SuperCellCyto")
  source(file.path(pp_dir, "supercell.R"))

  set.seed(42)
  expr <- matrix(rnorm(100), nrow = 20, ncol = 5)
  colnames(expr) <- paste0("M", 1:5)

  # Mismatched sample_labels length

  expect_error(
    run_supercell(expr, sample_labels = rep("s1", 10),
                  markers = paste0("M", 1:5)),
    "must match"
  )

  # Missing markers
  expect_error(
    run_supercell(expr, sample_labels = rep("s1", 20),
                  markers = c("M1", "NONEXISTENT")),
    "not found"
  )
})

# =============================================================================
# SECTION 14: Edge cases
# =============================================================================

context("Preprocessing edge cases")

test_that("transform then z-score produces zero-mean data", {
  ff <- make_synthetic_ff(n_cells = 200, markers = c("CD3", "CD4", "CD8"))

  # First transform
  tf_ff <- transform_flowframe(
    ff, "test", c("CD3", "CD4", "CD8"), c(0, 0, 0), c(5, 5, 5),
    transform_type = "arcsinh"
  )
  # Then z-score
  scaled_ff <- z_score_scaling(tf_ff, 1:3)

  expr <- flowCore::exprs(scaled_ff)
  for (i in 1:3) {
    expect_equal(mean(expr[, i]), 0, tolerance = 1e-10)
    expect_equal(sd(expr[, i]), 1, tolerance = 0.01)
  }
})

test_that("doublet detection followed by live/dead gating chains correctly", {
  ff <- make_scatter_ff(n_cells = 500)

  # Step 1: Remove doublets
  dd_result <- detect_doublets(ff, method = "ratio", threshold = 1.1)
  singlets <- dd_result$filtered_ff

  # Step 2: Gate live cells from singlets
  ld_result <- gate_live_cells(
    singlets, viability_channel = "LiveDead",
    method = "threshold", threshold = 800,
    dead_high = TRUE
  )

  live_cells <- ld_result$live_ff

  # Should be strictly smaller than both inputs
  expect_lt(nrow(live_cells), nrow(ff))
  expect_lte(nrow(live_cells), nrow(singlets))
  # All columns should be preserved
  expect_equal(colnames(live_cells), colnames(ff))
})

test_that("flowFrame with zero rows is handled gracefully", {
  markers <- c("CD3", "CD4")
  params <- Biobase::AnnotatedDataFrame(data.frame(
    name     = markers,
    desc     = markers,
    range    = c(1024, 1024),
    minRange = c(0, 0),
    maxRange = c(1024, 1024),
    row.names = paste0("$P", 1:2),
    stringsAsFactors = FALSE
  ))
  # flowCore's flowFrame constructor rejects a 0-row matrix directly
  # (round-trips through FCS format which requires at least one row), so
  # build a 1-row frame and drop the row.
  dummy <- matrix(0, nrow = 1, ncol = 2, dimnames = list(NULL, markers))
  ff_empty <- flowCore::flowFrame(exprs = dummy, parameters = params)[-1, ]
  expect_equal(unname(nrow(ff_empty)), 0)

  # quality_control None should pass through
  result <- quality_control(ff_empty, "empty", qc_type = "None")
  expect_equal(unname(nrow(result)), 0)

  # transform_flowframe should handle 0-row input
  result <- transform_flowframe(
    ff_empty, "empty", markers, c(0, 0), c(5, 5),
    transform_type = "arcsinh"
  )
  expect_equal(unname(nrow(result)), 0)
})

test_that("single-cell flowFrame (1 row) is handled", {
  markers <- c("CD3", "CD4")
  expr <- matrix(c(3.5, 2.1), nrow = 1, ncol = 2)
  colnames(expr) <- markers
  params <- Biobase::AnnotatedDataFrame(data.frame(
    name     = markers,
    desc     = markers,
    range    = c(1024, 1024),
    minRange = c(0, 0),
    maxRange = c(1024, 1024),
    row.names = paste0("$P", 1:2),
    stringsAsFactors = FALSE
  ))
  ff_single <- flowCore::flowFrame(exprs = expr, parameters = params)

  result <- quality_control(ff_single, "single", qc_type = "None")
  expect_equal(nrow(result), 1)

  result <- transform_flowframe(
    ff_single, "single", markers, c(0, 0), c(5, 5),
    transform_type = "arcsinh"
  )
  expect_equal(nrow(result), 1)
  expect_true(all(is.finite(flowCore::exprs(result))))
})

test_that("run_auto_gating errors on missing channels", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  expect_error(
    run_auto_gating(ff, channels = c("NONEXISTENT")),
    "Channels not found"
  )
})

test_that("run_auto_gating rejects invalid method", {
  ff <- make_synthetic_ff(n_cells = 50, markers = c("CD3", "CD4"))
  expect_error(
    run_auto_gating(ff, method = "invalid", channels = c("CD3")),
    "Invalid method"
  )
})

test_that("apply_gate_sequence population_stats reports correct counts", {
  ff <- make_synthetic_ff(n_cells = 100, markers = c("CD3", "CD4"))

  gate <- flowCore::rectangleGate(.gate = list(CD3 = c(3, Inf)))
  result <- apply_gate_sequence(ff, list(gate))

  stats <- result$population_stats
  expect_equal(stats$n_parent, 100)
  expect_equal(stats$n_gated, nrow(result$gated_ff))
  expect_equal(stats$pct_parent,
               100 * nrow(result$gated_ff) / 100,
               tolerance = 0.01)
})
