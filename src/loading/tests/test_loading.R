#!/usr/bin/env Rscript
# =============================================================================
# loading/tests/test_loading.R — Comprehensive tests for loading module
# =============================================================================
# Tests: load_fcs_raw, load_csv_expression, create_flowframes,
#        load_manual_gating, load_csv_gating, load_flowjo_wsp helpers
# =============================================================================

library(testthat)

# ---------------------------------------------------------------------------
# Resolve paths relative to this test file so the tests work regardless of
# whether they are run via testthat::test_dir("tests/") from within
# src/loading/ or from the project root.
# ---------------------------------------------------------------------------
.this_dir <- tryCatch(
  dirname(testthat::test_path()),
  error = function(e) getwd()
)
loading_dir <- normalizePath(file.path(.this_dir, ".."), mustWork = FALSE)
src_dir     <- normalizePath(file.path(loading_dir, ".."), mustWork = FALSE)
shared_dir  <- file.path(src_dir, "R_shared")
project_dir <- normalizePath(file.path(src_dir, ".."), mustWork = FALSE)
test_data   <- file.path(project_dir, "test_data")

# ---------------------------------------------------------------------------
# Source shared utilities (needed by all loading scripts)
# ---------------------------------------------------------------------------
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "fcs_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

# Source the loading scripts themselves (functions only, not CLI blocks)
source(file.path(loading_dir, "load_fcs_raw.R"))
source(file.path(loading_dir, "load_csv_expression.R"))
source(file.path(loading_dir, "create_flowframes.R"))

# Conditionally source scripts that need gui_utils
tryCatch({
  source(file.path(shared_dir, "gui_utils.R"))
  source(file.path(loading_dir, "load_csv_gating.R"))
  source(file.path(loading_dir, "load_manual_gating.R"))
}, error = function(e) {
  message("Skipping gui_utils-dependent scripts: ", e$message)
})

# ---------------------------------------------------------------------------
# Helpers for creating synthetic data
# ---------------------------------------------------------------------------
make_synthetic_flowframe <- function(n_cells = 100,
                                     markers = c("CD3", "CD4", "CD8",
                                                 "CD19", "CD56")) {
  set.seed(42)
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

# Create a temporary directory that is cleaned up after all tests
test_tmpdir <- tempfile("loading_tests_")
dir.create(test_tmpdir, recursive = TRUE)
withr::defer(unlink(test_tmpdir, recursive = TRUE), teardown_env())

# Write a few synthetic FCS files into the temp dir for FCS-based tests
.markers <- c("CD3", "CD4", "CD8", "CD19", "CD56")
.sample_names <- c("sample_001", "sample_002", "sample_003")
.fcs_dir <- file.path(test_tmpdir, "fcs")
dir.create(.fcs_dir)
for (i in seq_along(.sample_names)) {
  ff <- make_synthetic_flowframe(n_cells = 50 + i * 10, markers = .markers)
  flowCore::write.FCS(ff, file.path(.fcs_dir, paste0(.sample_names[i], ".fcs")))
}

# Write synthetic CSV expression files
.csv_dir <- file.path(test_tmpdir, "csv_expr")
dir.create(.csv_dir)
for (sname in .sample_names) {
  set.seed(which(.sample_names == sname))
  dt <- data.table::data.table(
    CD3  = abs(rnorm(40, 3, 1)),
    CD4  = abs(rnorm(40, 2, 1)),
    CD8  = abs(rnorm(40, 2, 1)),
    CD19 = abs(rnorm(40, 1, 0.5)),
    CD56 = abs(rnorm(40, 1.5, 0.8))
  )
  data.table::fwrite(dt, file.path(.csv_dir, paste0(sname, ".csv")), sep = ";")
}

# =============================================================================
# SECTION 1: load_fcs_raw
# =============================================================================

context("load_fcs_raw")

test_that("load_fcs_raw validates input directory", {
  expect_error(
    load_fcs_raw("/nonexistent/dir/that/does/not/exist"),
    "dir.exists"
  )
})

test_that("load_fcs_raw errors when directory has no FCS files", {
  empty_dir <- file.path(test_tmpdir, "empty_fcs")
  dir.create(empty_dir, showWarnings = FALSE)
  expect_error(load_fcs_raw(empty_dir), "No FCS files found")
})

test_that("load_fcs_raw loads synthetic FCS files from temp dir", {
  result <- load_fcs_raw(.fcs_dir)

  expect_type(result, "list")
  expect_named(result, c("flow_frames", "sample_names", "input_files"))
  expect_equal(length(result$flow_frames), 3)
  expect_equal(length(result$sample_names), 3)
  expect_equal(length(result$input_files), 3)

  # Sample names should be derived from filenames (without .fcs extension)
  expect_true(all(.sample_names %in% result$sample_names))

  # Each element should be a valid flowFrame with rows > 0

  for (ff in result$flow_frames) {
    expect_s4_class(ff, "flowFrame")
    expect_gt(nrow(ff), 0)
    expect_equal(ncol(ff), length(.markers))
  }
})

test_that("load_fcs_raw loads from test_data if it exists", {
  fcs_test_dir <- file.path(test_data, "fcs_samples")
  skip_if_not(dir.exists(fcs_test_dir), "test_data/fcs_samples not found")

  result <- load_fcs_raw(fcs_test_dir)
  expect_type(result, "list")
  expect_true("flow_frames" %in% names(result))
  expect_true("sample_names" %in% names(result))
  expect_equal(length(result$flow_frames), length(result$sample_names))

  for (ff in result$flow_frames) {
    expect_s4_class(ff, "flowFrame")
    expect_gt(nrow(ff), 0)
  }
})

test_that("load_fcs_raw respects custom sample_names", {
  custom_names <- c("alpha", "beta", "gamma")
  result <- load_fcs_raw(.fcs_dir, sample_names = custom_names)

  expect_equal(result$sample_names, custom_names)
  expect_equal(names(result$flow_frames), custom_names)
})

test_that("load_fcs_raw errors when sample_names length mismatches", {
  expect_error(
    load_fcs_raw(.fcs_dir, sample_names = c("only_one")),
    "Number of sample names"
  )
})

test_that("load_fcs_raw extract_csv saves CSV files", {
  csv_out <- file.path(test_tmpdir, "fcs_extracted_csvs")
  result <- load_fcs_raw(.fcs_dir, extract_csv = TRUE, output_dir = csv_out)

  expect_true(dir.exists(csv_out))
  csv_files <- list.files(csv_out, pattern = "\\.csv$")
  expect_equal(length(csv_files), 3)

  # Check one CSV has actual content
  one_csv <- data.table::fread(file.path(csv_out, csv_files[1]))
  expect_gt(nrow(one_csv), 0)
  expect_gt(ncol(one_csv), 0)
})

test_that("load_fcs_raw truncation parameter is forwarded", {
  # Should not error with truncation options
  result <- load_fcs_raw(.fcs_dir, truncation = TRUE)
  expect_type(result, "list")
  expect_equal(length(result$flow_frames), 3)
})

# =============================================================================
# SECTION 2: load_csv_expression
# =============================================================================

context("load_csv_expression")

test_that("load_csv_expression validates input directory", {
  expect_error(
    load_csv_expression("/nonexistent/csv/dir"),
    "dir.exists"
  )
})

test_that("load_csv_expression errors when no CSVs found", {
  empty_dir <- file.path(test_tmpdir, "empty_csv")
  dir.create(empty_dir, showWarnings = FALSE)
  expect_error(load_csv_expression(empty_dir), "No .csv files found")
})

test_that("load_csv_expression loads synthetic CSVs", {
  result <- load_csv_expression(.csv_dir)

  expect_type(result, "list")
  expect_equal(length(result), 3)

  for (nm in names(result)) {
    dt <- result[[nm]]
    expect_s3_class(dt, "data.table")
    expect_equal(nrow(dt), 40)
    expect_equal(ncol(dt), 5)
    expect_true(all(c("CD3", "CD4", "CD8", "CD19", "CD56") %in% names(dt)))
  }
})

test_that("load_csv_expression auto-detects separator", {
  # Write a comma-separated CSV
  comma_dir <- file.path(test_tmpdir, "comma_csv")
  dir.create(comma_dir, showWarnings = FALSE)
  dt <- data.table::data.table(A = 1:5, B = 6:10, C = 11:15)
  data.table::fwrite(dt, file.path(comma_dir, "test.csv"), sep = ",")

  result <- load_csv_expression(comma_dir)
  expect_equal(length(result), 1)
  expect_equal(nrow(result[[1]]), 5)
  expect_equal(ncol(result[[1]]), 3)
})

test_that("load_csv_expression loads test_data expression CSVs", {
  expr_dir <- file.path(test_data, "csv_expression")
  skip_if_not(dir.exists(expr_dir), "test_data/csv_expression not found")

  result <- load_csv_expression(expr_dir)
  expect_type(result, "list")
  expect_gt(length(result), 0)

  for (nm in names(result)) {
    expect_s3_class(result[[nm]], "data.table")
    expect_gt(nrow(result[[nm]]), 0)
    expect_gt(ncol(result[[nm]]), 0)
  }
})

test_that("load_csv_expression explicit sep overrides auto-detection", {
  result <- load_csv_expression(.csv_dir, sep = ";")
  expect_type(result, "list")
  expect_equal(length(result), 3)
})

# =============================================================================
# SECTION 3: create_flowframes
# =============================================================================

context("create_flowframes_from_matrices")

test_that("create_flowframes_from_matrices creates flowFrames from lists", {
  expr_list <- list(
    sampleA = matrix(rnorm(100), nrow = 20, ncol = 5,
                     dimnames = list(NULL, c("CD3", "CD4", "CD8", "CD19", "CD56"))),
    sampleB = matrix(rnorm(75), nrow = 15, ncol = 5,
                     dimnames = list(NULL, c("CD3", "CD4", "CD8", "CD19", "CD56")))
  )

  result <- create_flowframes_from_matrices(expr_list)

  expect_type(result, "list")
  expect_equal(length(result), 2)
  expect_named(result, c("sampleA", "sampleB"))

  expect_s4_class(result$sampleA, "flowFrame")
  expect_s4_class(result$sampleB, "flowFrame")
  expect_equal(nrow(result$sampleA), 20)
  expect_equal(nrow(result$sampleB), 15)
  expect_equal(ncol(result$sampleA), 5)
})

test_that("create_flowframes_from_matrices filters to markers_of_interest", {
  expr_list <- list(
    sample1 = matrix(rnorm(100), nrow = 20, ncol = 5,
                     dimnames = list(NULL, c("CD3", "CD4", "CD8", "CD19", "CD56")))
  )
  result <- create_flowframes_from_matrices(
    expr_list, markers_of_interest = c("CD3", "CD4")
  )

  expect_s4_class(result$sample1, "flowFrame")
  expect_equal(ncol(result$sample1), 2)
})

test_that("create_flowframes_from_matrices handles data.frames", {
  expr_list <- list(
    sample1 = data.frame(
      CD3  = rnorm(30),
      CD4  = rnorm(30),
      CD8  = rnorm(30)
    )
  )
  result <- create_flowframes_from_matrices(expr_list)

  expect_s4_class(result$sample1, "flowFrame")
  expect_equal(nrow(result$sample1), 30)
  expect_equal(ncol(result$sample1), 3)
})

test_that("create_flowframes_from_matrices save_fcs writes files", {
  expr_list <- list(
    testSample = matrix(rnorm(50), nrow = 10, ncol = 5,
                        dimnames = list(NULL, c("A", "B", "C", "D", "E")))
  )
  fcs_out <- file.path(test_tmpdir, "created_fcs")
  result <- create_flowframes_from_matrices(
    expr_list, save_fcs = TRUE, output_dir = fcs_out
  )

  expect_true(dir.exists(fcs_out))
  expect_true(file.exists(file.path(fcs_out, "testSample.fcs")))
})

test_that("create_flowframes_from_matrices ignores unmatched markers_of_interest", {
  expr_list <- list(
    s1 = matrix(rnorm(30), nrow = 10, ncol = 3,
                dimnames = list(NULL, c("CD3", "CD4", "CD8")))
  )
  # Request markers that partially overlap
  result <- create_flowframes_from_matrices(
    expr_list, markers_of_interest = c("CD3", "NONEXISTENT")
  )
  expect_s4_class(result$s1, "flowFrame")
  expect_equal(ncol(result$s1), 1)  # Only CD3 matched
})

# =============================================================================
# SECTION 4: create_merged_flowset
# =============================================================================

context("create_merged_flowset")

test_that("create_merged_flowset merges compatible flowFrames", {
  ff1 <- make_synthetic_flowframe(n_cells = 30, markers = c("A", "B", "C"))
  ff2 <- make_synthetic_flowframe(n_cells = 20, markers = c("A", "B", "C"))
  ffl <- list(s1 = ff1, s2 = ff2)

  fs <- create_merged_flowset(ffl)
  expect_s4_class(fs, "flowSet")
  expect_equal(length(fs), 2)
})

test_that("create_merged_flowset rejects incompatible column counts", {
  ff1 <- make_synthetic_flowframe(n_cells = 30, markers = c("A", "B", "C"))
  ff2 <- make_synthetic_flowframe(n_cells = 20, markers = c("A", "B"))
  ffl <- list(s1 = ff1, s2 = ff2)

  expect_error(create_merged_flowset(ffl), "different numbers of columns")
})

# =============================================================================
# SECTION 5: R_shared/fcs_utils helper functions
# =============================================================================

context("fcs_utils helpers")

test_that("create_csv_fcs_mapping creates correct CSV-to-FCS mapping", {
  csv_files <- file.path(.csv_dir, paste0(.sample_names, ".csv"))
  fcs_files <- list.files(.fcs_dir, pattern = "\\.fcs$", full.names = TRUE)

  # Write matching CSVs into the fcs dir for mapping test
  mapping_csv_dir <- file.path(test_tmpdir, "mapping_csv")
  dir.create(mapping_csv_dir, showWarnings = FALSE)
  for (sname in .sample_names) {
    file.copy(
      file.path(.csv_dir, paste0(sname, ".csv")),
      file.path(mapping_csv_dir, paste0(sname, ".csv"))
    )
  }
  mapping_csv_files <- list.files(mapping_csv_dir, pattern = "\\.csv$",
                                  full.names = TRUE)

  result <- create_csv_fcs_mapping(mapping_csv_files, fcs_files)

  expect_type(result, "character")
  expect_equal(length(result), 3)
  # Each CSV should map to a FCS with the same base name
  for (csv_path in mapping_csv_files) {
    fcs_path <- result[csv_path]
    csv_base <- gsub("[-_]", "", gsub("\\.csv$", "", basename(csv_path)))
    fcs_base <- gsub("[-_]", "", gsub("\\.fcs$", "", basename(fcs_path)))
    expect_true(grepl(fcs_base, csv_base, fixed = TRUE))
  }
})

test_that("get_fcs_event_count returns correct count", {
  fcs_file <- file.path(.fcs_dir, paste0(.sample_names[1], ".fcs"))
  count <- get_fcs_event_count(fcs_file)

  expect_type(count, "double")
  expect_equal(count, 60)  # 50 + 1*10
})

test_that("detect_csv_sep detects semicolon separator", {
  csv_path <- file.path(.csv_dir, paste0(.sample_names[1], ".csv"))
  sep <- detect_csv_sep(csv_path)
  expect_equal(sep, ";")
})

test_that("detect_csv_sep detects comma separator", {
  comma_file <- file.path(test_tmpdir, "comma_test.csv")
  writeLines("A,B,C\n1,2,3", comma_file)
  sep <- detect_csv_sep(comma_file)
  expect_equal(sep, ",")
})

test_that("list_csvs returns expected tibble structure", {
  result <- list_csvs(.csv_dir)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("file_path", "file_name", "file_size_mb") %in% names(result)))
  expect_equal(nrow(result), 3)
  expect_true(all(result$file_size_mb > 0))
})

test_that("check_tot validates $TOT on synthetic flowFrame", {
  fcs_path <- file.path(.fcs_dir, paste0(.sample_names[1], ".fcs"))
  ff <- flowCore::read.FCS(fcs_path, transformation = FALSE)
  # Should not error -- just produces a message or warning
  expect_message(check_tot(ff, index = 1), regexp = "sample 1")
})

# =============================================================================
# SECTION 6: R_shared/flowframe_utils helper functions
# =============================================================================

context("flowframe_utils helpers")

test_that("explore_data returns markers and column indices", {
  ff <- make_synthetic_flowframe(n_cells = 50, markers = c("CD3", "CD4", "CD8"))
  result <- explore_data(ff)

  expect_type(result, "list")
  expect_named(result, c("dyes", "markers", "mcols"))
  expect_equal(length(result$dyes), 3)
  expect_equal(length(result$markers), 3)
  expect_true(all(result$mcols %in% 1:3))
})

test_that("from_mtx_to_ffr creates valid flowFrame from matrix", {
  mtx <- matrix(rnorm(100), nrow = 20, ncol = 5,
                dimnames = list(NULL, c("A", "B", "C", "D", "E")))
  ff <- from_mtx_to_ffr(mtx)

  expect_s4_class(ff, "flowFrame")
  expect_equal(nrow(ff), 20)
  expect_equal(ncol(ff), 5)
  expect_equal(colnames(ff), c("A", "B", "C", "D", "E"))
})

test_that("from_mtx_to_ffr with custom marker names", {
  mtx <- matrix(rnorm(60), nrow = 20, ncol = 3)
  ff <- from_mtx_to_ffr(mtx, marker_names = c("X1", "X2", "X3"))

  expect_equal(colnames(ff), c("X1", "X2", "X3"))
})

test_that("build_combined_expr_matrix merges correctly", {
  ff1 <- make_synthetic_flowframe(n_cells = 30, markers = c("CD3", "CD4", "CD8"))
  ff2 <- make_synthetic_flowframe(n_cells = 20, markers = c("CD3", "CD4", "CD8"))
  flow_frames <- list(ff1, ff2)

  result <- build_combined_expr_matrix(
    flow_frames,
    sample_names = c("s1", "s2"),
    markers_of_interest = c("CD3", "CD4", "CD8")
  )

  expect_type(result, "list")
  expect_named(result, c("combined_matrix", "sample_labels"))
  expect_equal(nrow(result$combined_matrix), 50)  # 30 + 20
  expect_equal(ncol(result$combined_matrix), 3)
  expect_equal(length(result$sample_labels), 50)
  expect_equal(sum(result$sample_labels == "s1"), 30)
  expect_equal(sum(result$sample_labels == "s2"), 20)
})

test_that("generate_cellID creates correct IDs", {
  ff <- make_synthetic_flowframe(n_cells = 10)
  ids <- generate_cellID(ff, "test")

  expect_type(ids, "character")
  expect_equal(length(ids), 10)
  expect_true(all(startsWith(ids, "test_cell_")))
  expect_equal(ids[1], "test_cell_1")
  expect_equal(ids[10], "test_cell_10")
})

# =============================================================================
# SECTION 7: load_manual_gating (if available)
# =============================================================================

context("load_manual_gating")

test_that("load_manual_gating works with test_data CSV", {
  skip_if_not(exists("load_manual_gating"), "load_manual_gating not sourced")

  # Create synthetic manual gating CSV
  mg_dir <- file.path(test_tmpdir, "manual_gating")
  dir.create(mg_dir, showWarnings = FALSE)

  dt <- data.table::data.table(
    index = 1:50,
    `CD4+ T cells (Manual gating)` = c(rep(1L, 20), rep(0L, 30)),
    `CD8+ T cells (Manual gating)` = c(rep(0L, 20), rep(1L, 15), rep(0L, 15)),
    `B cells (Manual gating)` = c(rep(0L, 35), rep(1L, 15))
  )
  data.table::fwrite(dt, file.path(mg_dir, "sample_001.csv"), sep = ";")

  # Copy matching FCS so the mapping works
  file.copy(
    file.path(.fcs_dir, "sample_001.fcs"),
    file.path(mg_dir, "sample_001.fcs")
  )

  result <- load_manual_gating(
    csv_dir = mg_dir,
    gating_pattern = "Manual gating"
  )

  expect_type(result, "list")
  expect_gt(length(result), 0)

  labels <- result[[1]]
  expect_type(labels, "character")
  expect_equal(length(labels), 50)
  expect_true(all(unique(labels) %in%
    c("CD4+ T cells (Manual gating)",
      "CD8+ T cells (Manual gating)",
      "B cells (Manual gating)",
      "Unassigned")))
})

# =============================================================================
# SECTION 8: load_csv_gating (non-interactive mode)
# =============================================================================

context("load_csv_gating")

test_that("load_csv_gating works in manual mode with test_data", {
  skip_if_not(exists("load_csv_gating"), "load_csv_gating not sourced")
  gating_dir <- file.path(test_data, "csv_gating")
  skip_if_not(dir.exists(gating_dir), "test_data/csv_gating not found")

  out_dir <- file.path(test_tmpdir, "gating_output")

  result <- load_csv_gating(
    input_fcs_path = gating_dir,
    output_path = out_dir,
    processing_mode = "manual",
    output_format = "wide"
  )

  expect_type(result, "list")
  expect_true("n_processed" %in% names(result))
  expect_true("output_dir" %in% names(result))
  expect_gt(result$n_processed, 0)
})

test_that("load_csv_gating errors with no CSV files", {
  skip_if_not(exists("load_csv_gating"), "load_csv_gating not sourced")
  empty_dir <- file.path(test_tmpdir, "empty_gating")
  dir.create(empty_dir, showWarnings = FALSE)
  # Write a dummy FCS so the FCS check passes
  ff <- make_synthetic_flowframe(n_cells = 10)
  flowCore::write.FCS(ff, file.path(empty_dir, "dummy.fcs"))

  expect_error(load_csv_gating(empty_dir), "No CSV files found")
})

# =============================================================================
# SECTION 9: classify_csv_columns (from io_utils)
# =============================================================================

context("classify_csv_columns")

test_that("classify_csv_columns correctly classifies columns", {
  cols <- c(
    "Index",
    "FSC.A | FSC.A",
    "SSC.A | SSC.A",
    "TIME | TIME.A",
    "CD3 | BV421.A",
    "CD4 | BV510.A",
    "CD8 | BV605.A",
    "CD4+ (Clustering annotated)",
    "CD8+ (Clustering annotated)",
    "Clustering raw",
    "cyto_2d_x",
    "cyto_2d_y",
    "AF.P1 | AF.P1.A"
  )

  result <- classify_csv_columns(cols)

  expect_type(result, "list")
  expect_true("Index" %in% result$index)
  expect_true(any(grepl("FSC", result$scatter)))
  expect_true(any(grepl("SSC", result$scatter)))
  expect_true("CD4+ (Clustering annotated)" %in% result$gating)
  expect_true("Clustering raw" %in% result$gating)
  expect_true("cyto_2d_x" %in% result$dimreduction)
  expect_true("AF.P1 | AF.P1.A" %in% result$autofluorescence)
  # Markers should include CD3, CD4, CD8
  expect_true("CD3 | BV421.A" %in% result$markers)
  expect_true("CD4 | BV510.A" %in% result$markers)
})

# =============================================================================
# SECTION 10: cluster_utils helpers used by loading
# =============================================================================

context("cluster_utils helpers")

test_that("clean_cluster_name removes prefixes and suffixes", {
  expect_equal(clean_cluster_name("AI_Cluster_CD4"), "CD4")
  expect_equal(clean_cluster_name("CD8+ (Clustering annotated)"), "CD8+")
  expect_equal(clean_cluster_name(""), "Unknown")
  expect_equal(clean_cluster_name(NA), "Unknown")
})

test_that("labels_from_binary_df extracts labels correctly", {
  df <- data.frame(
    `Gate_A` = c(1, 0, 1, 0),
    `Gate_B` = c(0, 1, 1, 0),
    check.names = FALSE
  )
  result <- labels_from_binary_df(df)

  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_equal(result[[1]], "Gate_A")
  expect_equal(result[[2]], "Gate_B")
  expect_equal(sort(result[[3]]), c("Gate_A", "Gate_B"))
  expect_equal(result[[4]], character(0))
})

test_that("translate_labels translates using dictionary", {
  labels <- c("CD4+ T cells", "CD8+ T cells", "B cells", "Unknown")
  dict <- c(
    "CD4+ T cells" = "CD4_T",
    "CD8+ T cells" = "CD8_T",
    "B cells"      = "B"
  )
  result <- translate_labels(labels, dict)

  expect_equal(result, c("CD4_T", "CD8_T", "B", "Unknown"))
})

test_that("flatten_label_list respects gating tree depth", {
  label_lists <- list(
    c("CD3+", "CD4+"),
    c("CD3+", "CD8+"),
    c("B cells"),
    character(0)
  )
  gating_tree <- c("CD3+" = 1, "CD4+" = 2, "CD8+" = 2, "B cells" = 1)

  result <- flatten_label_list(label_lists, gating_tree)
  expect_equal(result, c("CD4+", "CD8+", "B cells", "root"))
})

# =============================================================================
# SECTION 11: load_flowjo_wsp (skipped unless CytoML available)
# =============================================================================

context("load_flowjo_wsp")

test_that("load_flowjo_wsp functions exist after sourcing", {
  skip_if_not_installed("CytoML")
  skip_if_not_installed("flowWorkspace")

  tryCatch(
    source(file.path(loading_dir, "load_flowjo_wsp.R")),
    error = function(e) skip(paste("Cannot source load_flowjo_wsp.R:", e$message))
  )

  expect_true(is.function(load_flowjo_wsp))
  expect_true(is.function(get_population_stats))
  expect_true(is.function(get_population_strings_list))
  expect_true(is.function(load_multiple_workspaces))
})

test_that("load_flowjo_wsp validates wsp_path existence", {
  skip_if_not_installed("CytoML")
  skip_if_not_installed("flowWorkspace")
  skip_if_not(exists("load_flowjo_wsp"), "load_flowjo_wsp not sourced")

  expect_error(load_flowjo_wsp("/nonexistent/file.wsp"), "file.exists")
})

test_that("load_multiple_workspaces validates matching lengths", {
  skip_if_not_installed("CytoML")
  skip_if_not_installed("flowWorkspace")
  skip_if_not(exists("load_multiple_workspaces"),
              "load_multiple_workspaces not sourced")

  expect_error(
    load_multiple_workspaces(
      wsp_paths = c("a.wsp", "b.wsp"),
      names = c("only_one")
    ),
    "length"
  )
})

# =============================================================================
# SECTION 12: Edge cases and integration
# =============================================================================

context("Loading edge cases")

test_that("load_fcs_raw with single-file directory", {
  single_dir <- file.path(test_tmpdir, "single_fcs")
  dir.create(single_dir, showWarnings = FALSE)
  ff <- make_synthetic_flowframe(n_cells = 25)
  flowCore::write.FCS(ff, file.path(single_dir, "only_one.fcs"))

  result <- load_fcs_raw(single_dir)
  expect_equal(length(result$flow_frames), 1)
  expect_equal(result$sample_names, "only_one")
  expect_equal(nrow(result$flow_frames[[1]]), 25)
})

test_that("create_flowframes roundtrip: CSV -> matrix -> flowFrame -> matrix", {
  # Load CSV expression data
  result_csv <- load_csv_expression(.csv_dir)
  # Convert to flowFrames
  flow_frames <- create_flowframes_from_matrices(result_csv)
  # Extract back to matrices
  for (nm in names(flow_frames)) {
    ff <- flow_frames[[nm]]
    mtx <- flowCore::exprs(ff)
    original <- as.matrix(result_csv[[nm]])
    # Dimensions should match
    expect_equal(nrow(mtx), nrow(original))
    expect_equal(ncol(mtx), ncol(original))
    # Values should be the same (within floating point tolerance)
    expect_equal(mtx, original, tolerance = 1e-6, ignore_attr = TRUE)
  }
})

test_that("create_safe_filename sanitizes special characters", {
  expect_equal(
    create_safe_filename("CD4+ T cells (gated)"),
    "CD4_T_cells_gated_"
  )
  expect_equal(create_safe_filename("simple"), "simple")
})
