library(testthat)

# ---------------------------------------------------------------------------
# Helpers: resolve shared and metadata_qc source directories
# ---------------------------------------------------------------------------
shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
mqc_dir    <- file.path(dirname(getwd()), "src", "metadata_qc")
test_data  <- file.path(dirname(getwd()), "test_data")

# Source shared utilities
source(file.path(shared_dir, "io_utils.R"))

# =============================================================================
# Tests for compute_fleiss_kappa() — overlap_statistics.R
# =============================================================================
test_that("compute_fleiss_kappa rejects non-dataframe input", {
  skip_if_not_installed("irr")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "overlap_statistics.R"))

  expect_error(
    compute_fleiss_kappa("not_a_dataframe", tempdir()),
    "label_lists must be a data frame"
  )
})

test_that("compute_fleiss_kappa rejects matrix input", {
  skip_if_not_installed("irr")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "overlap_statistics.R"))

  mtx <- matrix(1:10, nrow = 5, ncol = 2)
  expect_error(
    compute_fleiss_kappa(mtx, tempdir()),
    "label_lists must be a data frame"
  )
})

test_that("compute_fleiss_kappa requires at least 2 columns", {
  skip_if_not_installed("irr")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "overlap_statistics.R"))

  df <- data.frame(rater1 = c("A", "B", "C", "A", "B"))
  expect_error(
    compute_fleiss_kappa(df, tempdir()),
    "At least 2 raters"
  )
})

test_that("compute_fleiss_kappa works with perfect agreement", {
  skip_if_not_installed("irr")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "overlap_statistics.R"))

  # Perfect agreement: all raters assign same labels
  df <- data.frame(
    rater1 = c("A", "B", "C", "A", "B", "C", "A", "B"),
    rater2 = c("A", "B", "C", "A", "B", "C", "A", "B"),
    rater3 = c("A", "B", "C", "A", "B", "C", "A", "B")
  )
  out_dir <- file.path(tempdir(), "kappa_perfect_test")
  result <- compute_fleiss_kappa(df, out_dir)

  expect_type(result, "list")
  expect_length(result, 1)
  expect_equal(names(result), "Full data")

  # Perfect agreement should have kappa = 1
  expect_equal(result[["Full data"]]$value, 1)

  # Log file should be created
  expect_true(file.exists(file.path(out_dir, "fleiss_kappa_log.txt")))

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

test_that("compute_fleiss_kappa works with partial disagreement", {
  skip_if_not_installed("irr")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "overlap_statistics.R"))

  set.seed(42)
  n <- 30
  labels <- c("CD4+", "CD8+", "Other")
  df <- data.frame(
    rater1 = sample(labels, n, replace = TRUE),
    rater2 = sample(labels, n, replace = TRUE)
  )
  out_dir <- file.path(tempdir(), "kappa_partial_test")
  result <- compute_fleiss_kappa(df, out_dir)

  expect_type(result, "list")
  # Kappa should be a numeric value between -1 and 1
  kappa_val <- result[["Full data"]]$value
  expect_true(is.numeric(kappa_val))
  expect_true(kappa_val >= -1 && kappa_val <= 1)

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

test_that("compute_fleiss_kappa handles subset_specs", {
  skip_if_not_installed("irr")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "overlap_statistics.R"))

  df <- data.frame(
    r1 = c("A", "B", "A", "B", "A", "B", "A", "B"),
    r2 = c("A", "B", "A", "B", "A", "B", "A", "B"),
    r3 = c("A", "B", "B", "A", "A", "B", "A", "B")
  )
  out_dir <- file.path(tempdir(), "kappa_subset_test")

  result <- compute_fleiss_kappa(
    df, out_dir,
    subset_specs = list(c(1, 2), c(2, 3)),
    subset_labels = c("Raters 1-2", "Raters 2-3")
  )

  expect_length(result, 2)
  expect_equal(names(result), c("Raters 1-2", "Raters 2-3"))
  # Raters 1-2 have perfect agreement
  expect_equal(result[["Raters 1-2"]]$value, 1)

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

test_that("compute_fleiss_kappa respects row_range", {
  skip_if_not_installed("irr")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "overlap_statistics.R"))

  df <- data.frame(
    r1 = c("A", "B", "C", "A", "B", "C", "A", "B", "C", "A"),
    r2 = c("A", "B", "C", "A", "B", "C", "A", "B", "C", "A")
  )
  out_dir <- file.path(tempdir(), "kappa_rowrange_test")

  # Use only first 5 rows
  result <- compute_fleiss_kappa(df, out_dir, row_range = 1:5)

  expect_type(result, "list")
  # With perfect agreement on subset, kappa = 1
  expect_equal(result[["Full data"]]$value, 1)

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

# =============================================================================
# Tests for get_fcs_keyword() — check_acquisition.R
# =============================================================================
test_that("get_fcs_keyword extracts $DATE from test FCS files", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  fcs_path <- file.path(test_data, "fcs_samples", "sample_001.fcs")
  skip_if(!file.exists(fcs_path), "Test FCS file not found")

  result <- get_fcs_keyword(fcs_path, "$DATE")
  expect_type(result, "character")
  expect_length(result, 1)
  # Should return a date string or "Unknown"
  expect_true(nchar(result) > 0)
})

test_that("get_fcs_keyword returns Unknown for missing keyword", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  fcs_path <- file.path(test_data, "fcs_samples", "sample_001.fcs")
  skip_if(!file.exists(fcs_path), "Test FCS file not found")

  result <- get_fcs_keyword(fcs_path, "$NONEXISTENT_KEYWORD_XYZ")
  expect_equal(result, "Unknown")
})

test_that("get_fcs_keyword returns Error for nonexistent file", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  result <- get_fcs_keyword("/nonexistent/path/fake.fcs", "$DATE")
  expect_equal(result, "Error")
})

test_that("get_fcs_keyword extracts $TOT from test FCS files", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  fcs_path <- file.path(test_data, "fcs_samples", "sample_001.fcs")
  skip_if(!file.exists(fcs_path), "Test FCS file not found")

  result <- get_fcs_keyword(fcs_path, "$TOT")
  expect_type(result, "character")
  # $TOT should be a numeric string representing total events
  expect_true(result != "Unknown" && result != "Error")
})

# =============================================================================
# Tests for get_fcs_metadata() — check_acquisition.R
# =============================================================================
test_that("get_fcs_metadata returns named vector with expected keys", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  fcs_path <- file.path(test_data, "fcs_samples", "sample_001.fcs")
  skip_if(!file.exists(fcs_path), "Test FCS file not found")

  result <- get_fcs_metadata(fcs_path)

  expect_type(result, "character")
  expect_true("acquisition_date" %in% names(result))
  expect_true("file_size" %in% names(result))
  # file_size should be a positive number or "Unknown"
  if (result["file_size"] != "Unknown") {
    expect_true(as.numeric(result["file_size"]) > 0)
  }
})

test_that("get_fcs_metadata returns Error for invalid file", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  result <- get_fcs_metadata("/nonexistent/path/fake.fcs")
  expect_equal(unname(result["acquisition_date"]), "Error")
  expect_equal(unname(result["file_size"]), "Error")
})

# =============================================================================
# Tests for explore_fcs_keywords() — check_acquisition.R
# =============================================================================
test_that("explore_fcs_keywords returns keyword list from valid FCS", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  fcs_path <- file.path(test_data, "fcs_samples", "sample_001.fcs")
  skip_if(!file.exists(fcs_path), "Test FCS file not found")

  result <- explore_fcs_keywords(fcs_path)
  expect_type(result, "list")
  # Should contain standard FCS keywords
  expect_true("$TOT" %in% names(result) || "$DATATYPE" %in% names(result))
})

test_that("explore_fcs_keywords returns NULL for invalid file", {
  skip_if_not_installed("flowCore")
  source(file.path(mqc_dir, "check_acquisition.R"))

  result <- explore_fcs_keywords("/nonexistent/path/fake.fcs")
  expect_null(result)
})

# =============================================================================
# Tests for check_acquisition_dates() — check_acquisition.R
# =============================================================================
test_that("check_acquisition_dates processes test_data FCS directory", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "check_acquisition.R"))

  fcs_dir <- file.path(test_data, "fcs_samples")
  skip_if(!dir.exists(fcs_dir), "Test FCS directory not found")

  out_dir <- file.path(tempdir(), "acq_dates_test")
  result <- check_acquisition_dates(fcs_dir, out_dir)

  expect_true(is.data.frame(result))
  expect_true("file_relative_path" %in% colnames(result))
  expect_true("acquisition_date" %in% colnames(result))
  expect_true("file_size" %in% colnames(result))
  expect_true("unique_file_id" %in% colnames(result))

  # Should have one row per FCS file (3 test files)
  expect_equal(nrow(result), 3)

  # Output CSV should be created
  expect_true(file.exists(file.path(out_dir, "acquisition_dates.csv")))

  # unique_file_id should be integer
  expect_true(is.integer(result$unique_file_id))

  # Cleanup
  unlink(out_dir, recursive = TRUE)
})

test_that("check_acquisition_dates errors on nonexistent folder", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "check_acquisition.R"))

  expect_error(
    check_acquisition_dates("/nonexistent/folder", tempdir()),
    "Folder does not exist"
  )
})

test_that("check_acquisition_dates errors on folder with no FCS files", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("data.table")
  source(file.path(mqc_dir, "check_acquisition.R"))

  # Create empty temp directory
  empty_dir <- file.path(tempdir(), "empty_fcs_dir")
  dir.create(empty_dir, showWarnings = FALSE)

  expect_error(
    check_acquisition_dates(empty_dir, tempdir()),
    "No FCS files found"
  )

  # Cleanup
  unlink(empty_dir, recursive = TRUE)
})

# =============================================================================
# Tests for check_metadata.R — check_metadata()
# =============================================================================
test_that("check_metadata creates plots from synthetic metadata CSV", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("forcats")
  skip_if_not_installed("stringr")
  skip_if_not_installed("scales")

  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(mqc_dir, "check_metadata.R"))

  # Create synthetic metadata CSV
  meta_csv <- file.path(tempdir(), "test_metadata.csv")
  md <- data.frame(
    SampleID = paste0("S", 1:10),
    Condition = rep(c("Control", "Treatment"), each = 5),
    Sex = c("M", "F", "M", "F", "M", "F", "M", "F", "M", "F"),
    Age = c("20-30", "30-40", "20-30", "40-50", "30-40",
            "20-30", "30-40", "40-50", "20-30", "?"),
    stringsAsFactors = FALSE
  )
  write.csv(md, meta_csv, row.names = FALSE)

  out_dir <- file.path(tempdir(), "check_metadata_test")

  expect_no_error(
    check_metadata(meta_csv, out_dir, sep = ",")
  )

  # Barplot output files should be created
  # save_plot_formats creates pdf/png/svg files
  barplot_files <- list.files(out_dir, pattern = "consensus_barplots")
  expect_true(length(barplot_files) > 0)

  # Summary CSV should exist
  expect_true(file.exists(file.path(out_dir, "metadata_summary_table.csv")))

  # Cleanup
  unlink(out_dir, recursive = TRUE)
  unlink(meta_csv)
})

test_that("check_metadata errors when no columns selected", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("forcats")
  skip_if_not_installed("stringr")
  skip_if_not_installed("scales")

  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(mqc_dir, "check_metadata.R"))

  meta_csv <- file.path(tempdir(), "test_metadata_empty_cols.csv")
  md <- data.frame(A = 1:5, B = 6:10)
  write.csv(md, meta_csv, row.names = FALSE)

  out_dir <- file.path(tempdir(), "check_metadata_empty_test")

  # Skip all columns so nothing is left to analyze
  expect_error(
    check_metadata(meta_csv, out_dir, sep = ",",
                   skip_columns = c("A", "B")),
    "No columns selected"
  )

  # Cleanup
  unlink(meta_csv)
  unlink(out_dir, recursive = TRUE)
})

test_that("check_metadata warns about missing only_cols", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("forcats")
  skip_if_not_installed("stringr")
  skip_if_not_installed("scales")

  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(mqc_dir, "check_metadata.R"))

  meta_csv <- file.path(tempdir(), "test_metadata_only_cols.csv")
  md <- data.frame(A = 1:5, B = 6:10)
  write.csv(md, meta_csv, row.names = FALSE)

  out_dir <- file.path(tempdir(), "check_metadata_onlycols_test")

  # Request columns that don't exist (+ one that does)
  expect_warning(
    check_metadata(meta_csv, out_dir, sep = ",",
                   only_cols = c("A", "NONEXISTENT")),
    "Columns not found"
  )

  # Cleanup
  unlink(meta_csv)
  unlink(out_dir, recursive = TRUE)
})

# =============================================================================
# Tests for compare_gatings.R — compare_gatings()
# =============================================================================
test_that("compare_gatings module defines expected functions", {
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(mqc_dir, "compare_gatings.R"))

  expect_true(is.function(compare_gatings))
})

test_that("compare_gatings errors when target FCS is missing from workspace", {
  skip_if_not_installed("data.table")
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(mqc_dir, "compare_gatings.R"))

  pop_strings <- list(
    WS1 = list("file_A.fcs" = c("CD4+", "CD8+")),
    WS2 = list("file_A.fcs" = c("CD4+", "CD8+"))
  )

  expect_error(
    compare_gatings(
      pop_strings, target_fcs = "file_B.fcs",
      label_dict = c("CD4+" = "CD4+"),
      output_dir = tempdir()
    ),
    "Target FCS 'file_B.fcs' missing"
  )
})
