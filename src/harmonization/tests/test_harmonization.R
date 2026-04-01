library(testthat)

# ---------------------------------------------------------------------------
# Helpers: resolve shared and harmonization source directories
# ---------------------------------------------------------------------------
shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
harm_dir   <- file.path(dirname(getwd()), "src", "harmonization")

# Source shared utilities
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

# ---------------------------------------------------------------------------
# Helper: create synthetic flowFrames for testing
# ---------------------------------------------------------------------------
make_test_flowframe <- function(n_cells, markers, marker_descs = markers,
                                 seed = 42) {
  skip_if_not_installed("flowCore")
  set.seed(seed)
  expr_data <- matrix(abs(rnorm(n_cells * length(markers), mean = 3, sd = 1)),
                      nrow = n_cells, ncol = length(markers))
  colnames(expr_data) <- markers
  params <- Biobase::AnnotatedDataFrame(
    data.frame(
      name = markers,
      desc = marker_descs,
      range = rep(1024, length(markers)),
      minRange = rep(0, length(markers)),
      maxRange = rep(1024, length(markers)),
      row.names = paste0("$P", seq_along(markers)),
      stringsAsFactors = FALSE
    )
  )
  flowCore::flowFrame(exprs = expr_data, parameters = params)
}

# =============================================================================
# Tests for panel_comparison.R — compare_panels_detailed()
# =============================================================================
test_that("compare_panels_detailed works with identical panels", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_comparison.R"))

  markers <- c("CD3", "CD4", "CD8", "CD19")
  ff1 <- make_test_flowframe(30, markers, seed = 1)
  ff2 <- make_test_flowframe(40, markers, seed = 2)

  result <- compare_panels_detailed(list(S1 = ff1, S2 = ff2))

  expect_type(result, "list")
  expect_true("common_markers" %in% names(result))
  expect_true("not_common" %in% names(result))
  expect_true("presence_matrix" %in% names(result))

  # All markers should be common
  expect_equal(sort(result$common_markers), sort(markers))
  expect_length(result$not_common, 0)
})

test_that("compare_panels_detailed detects different panels", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_comparison.R"))

  markers1 <- c("CD3", "CD4", "CD8", "CD19")
  markers2 <- c("CD3", "CD4", "CD56", "HLA-DR")
  ff1 <- make_test_flowframe(30, markers1, seed = 1)
  ff2 <- make_test_flowframe(30, markers2, seed = 2)

  result <- compare_panels_detailed(list(S1 = ff1, S2 = ff2))

  # Common = CD3, CD4 only
  expect_equal(sort(result$common_markers), c("CD3", "CD4"))

  # Not common should include CD8, CD19, CD56, HLA-DR
  expect_true("CD8" %in% result$not_common)
  expect_true("CD19" %in% result$not_common)
  expect_true("CD56" %in% result$not_common)
  expect_true("HLA-DR" %in% result$not_common)

  # Presence matrix should have correct dimensions
  expect_equal(ncol(result$presence_matrix), 2)
  expect_equal(nrow(result$presence_matrix), length(result$all_markers))
})

test_that("compare_panels_detailed presence matrix is boolean", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_comparison.R"))

  markers1 <- c("CD3", "CD4", "CD8")
  markers2 <- c("CD3", "CD4", "CD56")
  ff1 <- make_test_flowframe(20, markers1, seed = 1)
  ff2 <- make_test_flowframe(20, markers2, seed = 2)

  result <- compare_panels_detailed(list(A = ff1, B = ff2))

  pm <- result$presence_matrix
  expect_true(all(pm %in% c(TRUE, FALSE)))

  # CD8 present in A but not B
  expect_true(pm["CD8", "A"])
  expect_false(pm["CD8", "B"])

  # CD56 present in B but not A
  expect_false(pm["CD56", "A"])
  expect_true(pm["CD56", "B"])
})

# =============================================================================
# Tests for panel_homogenization.R — translate_vector()
# =============================================================================
test_that("translate_vector translates known keys", {
  source(file.path(harm_dir, "panel_homogenization.R"))

  vec <- c("A", "B", "C")
  dict <- list(A = "alpha", B = "beta")
  result <- translate_vector(vec, dict)
  expect_equal(unname(result), c("alpha", "beta", "Unknown"))
})

test_that("translate_vector returns Unknown for unmatched keys", {
  source(file.path(harm_dir, "panel_homogenization.R"))

  vec <- c("X", "Y")
  dict <- list(A = "alpha")
  result <- translate_vector(vec, dict)
  expect_true(all(result == "Unknown"))
})

# =============================================================================
# Tests for panel_homogenization.R — homogenize_flowframes()
# =============================================================================
test_that("homogenize_flowframes produces frames with identical columns", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers1 <- c("CD3", "CD4", "CD8", "CD19", "CD56")
  markers2 <- c("CD3", "CD4", "CD8", "CD19", "HLA-DR")
  ff1 <- make_test_flowframe(30, markers1, seed = 1)
  ff2 <- make_test_flowframe(40, markers2, seed = 2)

  target_markers <- c("CD3", "CD4", "CD8")  # common subset
  result <- homogenize_flowframes(
    flow_frames = list(ff1, ff2),
    sample_names = c("S1", "S2"),
    markers_of_interest = target_markers
  )

  expect_length(result, 2)
  expect_equal(names(result), c("S1", "S2"))

  # Both should have exactly the target markers as column names
  cols1 <- colnames(flowCore::exprs(result[["S1"]]))
  cols2 <- colnames(flowCore::exprs(result[["S2"]]))
  expect_equal(cols1, target_markers)
  expect_equal(cols2, target_markers)
})

test_that("homogenize_flowframes preserves cell counts", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers <- c("CD3", "CD4", "CD8")
  ff1 <- make_test_flowframe(25, markers, seed = 1)
  ff2 <- make_test_flowframe(35, markers, seed = 2)

  result <- homogenize_flowframes(
    flow_frames = list(ff1, ff2),
    sample_names = c("A", "B"),
    markers_of_interest = markers
  )

  expect_equal(nrow(flowCore::exprs(result[["A"]])), 25)
  expect_equal(nrow(flowCore::exprs(result[["B"]])), 35)
})

test_that("homogenize_flowframes fills NA for missing markers", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers1 <- c("CD3", "CD4")
  markers2 <- c("CD3", "CD8")
  ff1 <- make_test_flowframe(20, markers1, seed = 1)
  ff2 <- make_test_flowframe(20, markers2, seed = 2)

  # Request a marker that only one frame has
  target <- c("CD3", "CD4", "CD8")
  expect_warning(
    result <- homogenize_flowframes(
      flow_frames = list(ff1, ff2),
      sample_names = c("S1", "S2"),
      markers_of_interest = target
    ),
    "Markers filled with NA"
  )

  # S1 should have NA in CD8 column; S2 should have NA in CD4 column
  expr_s1 <- flowCore::exprs(result[["S1"]])
  expr_s2 <- flowCore::exprs(result[["S2"]])
  expect_true(all(is.na(expr_s1[, "CD8"])))
  expect_true(all(is.na(expr_s2[, "CD4"])))
  # Non-missing should be populated
  expect_false(any(is.na(expr_s1[, "CD3"])))
})

test_that("homogenize_flowframes detects markername collision", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers <- c("CD3", "CD4", "CD8")
  ff1 <- make_test_flowframe(20, markers, seed = 1)

  # Two old markers mapping to same new name
  bad_replace <- list("CD3" = "target", "CD4" = "target")
  expect_error(
    homogenize_flowframes(
      flow_frames = list(ff1),
      sample_names = c("S1"),
      markers_of_interest = markers,
      markername_replacing = bad_replace
    ),
    "markername_replacing collision"
  )
})

test_that("homogenize_flowframes saves CSV when requested", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers <- c("CD3", "CD4")
  ff1 <- make_test_flowframe(15, markers, seed = 1)
  out_dir <- file.path(tempdir(), "hg_csv_test")

  result <- homogenize_flowframes(
    flow_frames = list(ff1),
    sample_names = c("TestSample"),
    markers_of_interest = markers,
    save_csv = TRUE,
    output_dir = out_dir
  )

  expect_true(file.exists(file.path(out_dir, "hg_TestSample.csv")))

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

# =============================================================================
# Tests for batch_correction.R — run_batch_correction()
# =============================================================================
test_that("run_batch_correction rejects non-homogenized frames", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("harmony")
  source(file.path(harm_dir, "batch_correction.R"))

  # Create two frames with different column sets
  markers1 <- c("CD3", "CD4")
  markers2 <- c("CD3", "CD8")
  ff1 <- make_test_flowframe(20, markers1, seed = 1)
  ff2 <- make_test_flowframe(20, markers2, seed = 2)

  expect_error(
    run_batch_correction(
      flow_frames = list(S1 = ff1, S2 = ff2),
      sample_names = c("S1", "S2"),
      batch_labels = c("Batch1", "Batch2")
    ),
    "different markers"
  )
})

# =============================================================================
# Tests for combat.R — run_combat()
# =============================================================================
test_that("run_combat validates input types", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  # Non-matrix input
  expect_error(
    run_combat("not_a_matrix", c("A", "B")),
    "expr_matrix must be a numeric matrix"
  )
})

test_that("run_combat validates batch_vector length", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  mtx <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")

  # batch_vector length doesn't match nrow
  expect_error(
    run_combat(mtx, c("A", "B")),
    "batch_vector length"
  )
})

test_that("run_combat requires at least 2 batches", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  mtx <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")

  # Single batch
  expect_error(
    run_combat(mtx, rep("A", 20)),
    "At least 2 batches"
  )
})

test_that("run_combat rejects invalid method", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  mtx <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")
  batch <- rep(c("A", "B"), each = 10)

  expect_error(
    run_combat(mtx, batch, method = "unknown"),
    "method must be 'parametric' or 'nonparametric'"
  )
})

test_that("run_combat validates covariates dimensions", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  mtx <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")
  batch <- rep(c("A", "B"), each = 10)

  # Covariates with wrong number of rows
  bad_cov <- data.frame(condition = c("ctrl", "treat"))
  expect_error(
    run_combat(mtx, batch, covariates = bad_cov),
    "covariates nrow"
  )
})

test_that("run_combat rejects non-dataframe covariates", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  mtx <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")
  batch <- rep(c("A", "B"), each = 10)

  expect_error(
    run_combat(mtx, batch, covariates = "not_a_dataframe"),
    "covariates must be a data.frame"
  )
})

test_that("run_combat rejects covariates with NA values", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  mtx <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")
  batch <- rep(c("A", "B"), each = 10)
  cov <- data.frame(condition = c(rep("ctrl", 10), rep(NA, 10)))

  expect_error(
    run_combat(mtx, batch, covariates = cov),
    "Covariates contain NA"
  )
})

test_that("run_combat returns corrected matrix with same dimensions", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  set.seed(42)
  n <- 40
  mtx <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(mtx) <- c("CD3", "CD4", "CD8", "CD19")
  # Add batch effect
  batch <- rep(c("A", "B"), each = n / 2)
  mtx[batch == "B", ] <- mtx[batch == "B", ] + 2

  result <- run_combat(mtx, batch)

  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mtx))
  expect_equal(colnames(result), colnames(mtx))
})

test_that("run_combat accepts data.frame input", {
  skip_if_not_installed("sva")
  source(file.path(harm_dir, "combat.R"))

  set.seed(42)
  n <- 40
  df <- data.frame(
    CD3  = rnorm(n),
    CD4  = rnorm(n),
    CD8  = rnorm(n)
  )
  batch <- rep(c("A", "B"), each = n / 2)

  result <- run_combat(df, batch)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), 3)
})

# =============================================================================
# Tests for batch_correction_dispatch.R — run_batch_correction_dispatch()
# =============================================================================
test_that("batch_correction_dispatch with method='none' passes through matrix", {
  source(file.path(shared_dir, "batch_correction_dispatch.R"))

  set.seed(42)
  mtx <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")
  batch <- rep(c("A", "B"), each = 5)
  samples <- rep(c("S1", "S2"), each = 5)

  result <- run_batch_correction_dispatch(
    data = mtx,
    method = "none",
    batch_labels = batch,
    sample_names = samples
  )

  expect_type(result, "list")
  expect_true("corrected_matrix" %in% names(result))
  expect_true("sample_labels" %in% names(result))
  expect_true("batch_labels" %in% names(result))
  expect_true("markers" %in% names(result))

  # Passthrough: matrix should be identical
  expect_identical(result$corrected_matrix, mtx)
  expect_equal(result$sample_labels, samples)
  expect_equal(result$batch_labels, batch)
  expect_equal(result$markers, c("CD3", "CD4", "CD8"))
  expect_null(result$corrected_frames)
})

test_that("batch_correction_dispatch with method='none' passes through flowFrames", {
  skip_if_not_installed("flowCore")
  source(file.path(shared_dir, "batch_correction_dispatch.R"))

  markers <- c("CD3", "CD4", "CD8")
  ff1 <- make_test_flowframe(20, markers, seed = 1)
  ff2 <- make_test_flowframe(30, markers, seed = 2)
  ff_list <- list(S1 = ff1, S2 = ff2)

  result <- run_batch_correction_dispatch(
    data = ff_list,
    method = "none",
    batch_labels = c("BatchA", "BatchB"),
    markers_of_interest = markers,
    sample_names = c("S1", "S2")
  )

  expect_type(result, "list")
  expect_equal(nrow(result$corrected_matrix), 50)  # 20 + 30
  expect_equal(ncol(result$corrected_matrix), 3)
  expect_equal(length(result$sample_labels), 50)
  expect_equal(length(result$batch_labels), 50)
})

test_that("batch_correction_dispatch rejects invalid method", {
  source(file.path(shared_dir, "batch_correction_dispatch.R"))

  mtx <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")

  expect_error(
    run_batch_correction_dispatch(
      data = mtx,
      method = "invalid_method",
      batch_labels = rep("A", 10)
    ),
    "Unknown batch correction method"
  )
})

test_that("batch_correction_dispatch method is case-insensitive", {
  source(file.path(shared_dir, "batch_correction_dispatch.R"))

  mtx <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(mtx) <- c("CD3", "CD4", "CD8")
  batch <- rep(c("A", "B"), each = 5)
  samples <- rep(c("S1", "S2"), each = 5)

  # "NONE" in uppercase should work
  result <- run_batch_correction_dispatch(
    data = mtx,
    method = "NONE",
    batch_labels = batch,
    sample_names = samples
  )
  expect_identical(result$corrected_matrix, mtx)
})

# =============================================================================
# Tests for reorder_exprmtx()
# =============================================================================
test_that("reorder_exprmtx reorders columns to target order", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers <- c("CD3", "CD4", "CD8", "CD19")
  ff <- make_test_flowframe(15, markers, seed = 1)

  # Request reversed subset
  target <- c("CD19", "CD3")
  result <- reorder_exprmtx(ff, target)

  expect_equal(colnames(result), target)
  expect_equal(nrow(result), 15)
  expect_equal(ncol(result), 2)
})

test_that("reorder_exprmtx warns about missing markers and fills NA", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers <- c("CD3", "CD4")
  ff <- make_test_flowframe(10, markers, seed = 1)

  expect_warning(
    result <- reorder_exprmtx(ff, c("CD3", "CD4", "CD56")),
    "Markers filled with NA"
  )
  expect_true(all(is.na(result[, "CD56"])))
  expect_false(any(is.na(result[, "CD3"])))
})

# =============================================================================
# Tests for generate_cellID_list() (internal helper)
# =============================================================================
test_that("generate_cellID_list produces correct ranges", {
  skip_if_not_installed("flowCore")
  source(file.path(harm_dir, "panel_homogenization.R"))

  markers <- c("CD3", "CD4")
  ff1 <- make_test_flowframe(10, markers, seed = 1)
  ff2 <- make_test_flowframe(15, markers, seed = 2)

  ids <- generate_cellID_list(list(ff1, ff2))

  expect_length(ids, 2)
  expect_equal(ids[[1]], 1:10)
  expect_equal(ids[[2]], 11:25)
})
