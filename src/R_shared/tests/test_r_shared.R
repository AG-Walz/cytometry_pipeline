# =============================================================================
# Tests for R_shared/ — Shared Utility Library
# =============================================================================
# Run with: testthat::test_file("test_r_shared.R")
# =============================================================================

library(testthat)

# ---------------------------------------------------------------------------
# Helper: resolve shared_dir from test file location
# ---------------------------------------------------------------------------
# When run via testthat::test_dir("tests/") from R_shared/, test_path("..")
# points to R_shared. When run directly, we try sys.frame or fall back to
# the known absolute path.
shared_dir <- tryCatch(
  normalizePath(testthat::test_path(".."), mustWork = TRUE),
  error = function(e) {
    # Fallback 1: derive from this script's file path (works with Rscript/source)
    script_dir <- tryCatch(
      normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE),
      error = function(e2) NULL
    )
    if (!is.null(script_dir) && file.exists(file.path(script_dir, "io_utils.R"))) {
      return(script_dir)
    }
    # Fallback 2: assume cwd is R_shared/tests or R_shared
    for (candidate in c(file.path(getwd(), ".."), getwd())) {
      cand <- normalizePath(candidate, mustWork = FALSE)
      if (file.exists(file.path(cand, "io_utils.R"))) return(cand)
    }
    # Fallback 3: hardcoded absolute path
    "/mnt/volume/nacho/cytometry/src/R_shared"
  }
)
stopifnot(file.exists(file.path(shared_dir, "io_utils.R")))

# %||% polyfill (needed before sourcing anything)
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# Source all modules under test
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))
source(file.path(shared_dir, "clustering_dispatch.R"))
source(file.path(shared_dir, "dimred_dispatch.R"))
source(file.path(shared_dir, "batch_correction_dispatch.R"))

# Only source flowframe_utils if flowCore is available
has_flowcore <- requireNamespace("flowCore", quietly = TRUE) &&
  requireNamespace("Biobase", quietly = TRUE)
if (has_flowcore) {
  source(file.path(shared_dir, "flowframe_utils.R"))
}


# =============================================================================
# 1. io_utils.R tests
# =============================================================================
context("io_utils.R")

test_that("detect_csv_sep detects semicolon separator", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("a;b;c\n1;2;3", tmp)
  expect_equal(detect_csv_sep(tmp), ";")
})

test_that("detect_csv_sep detects comma separator", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("a,b,c\n1,2,3", tmp)
  expect_equal(detect_csv_sep(tmp), ",")
})

test_that("detect_csv_sep defaults to comma for tab-separated content", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("a\tb\tc\n1\t2\t3", tmp)
  # No semicolons found -> falls back to ","
  expect_equal(detect_csv_sep(tmp), ",")
})

test_that("make_outdir creates directory and is idempotent", {
  tmp_dir <- file.path(tempdir(), paste0("test_make_outdir_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  expect_false(dir.exists(tmp_dir))
  make_outdir(tmp_dir)
  expect_true(dir.exists(tmp_dir))

  # Calling again should not error (idempotent)
  expect_silent(make_outdir(tmp_dir))
  expect_true(dir.exists(tmp_dir))
})

test_that("make_outdir creates nested directories", {
  tmp_dir <- file.path(tempdir(), paste0("test_nested_", Sys.getpid()),
                       "level1", "level2")
  on.exit(unlink(file.path(tempdir(), paste0("test_nested_", Sys.getpid())),
                 recursive = TRUE), add = TRUE)

  make_outdir(tmp_dir)
  expect_true(dir.exists(tmp_dir))
})

test_that("classify_csv_columns classifies column types correctly", {
  col_names <- c(
    "Index",
    "SSC.A | SSC.A",
    "FSC.H | FSC.H",
    "TIME.1 | TIME.1",
    "CD3 | BV421.A",
    "CD4 | FITC.A",
    "CD8 | PE.A",
    "T-Cells (Manual gating)",
    "B-Cells (Clustering annotated)",
    "Clustering raw",
    "cyto_2d_x",
    "cyto_2d_y",
    "AF.P1 | AF.P1.A",
    "AF-P2 | AF-P2.A",
    "BUV496.A | BUV496.A"
  )

  result <- classify_csv_columns(col_names)

  expect_type(result, "list")
  expect_true("Index" %in% result$index)
  expect_true("SSC.A | SSC.A" %in% result$scatter)
  expect_true("FSC.H | FSC.H" %in% result$scatter)
  expect_true("TIME.1 | TIME.1" %in% result$scatter)
  expect_true("T-Cells (Manual gating)" %in% result$gating)
  expect_true("B-Cells (Clustering annotated)" %in% result$gating)
  expect_true("Clustering raw" %in% result$gating)
  expect_true("cyto_2d_x" %in% result$dimreduction)
  expect_true("cyto_2d_y" %in% result$dimreduction)
  expect_true("AF.P1 | AF.P1.A" %in% result$autofluorescence)
  expect_true("AF-P2 | AF-P2.A" %in% result$autofluorescence)
  expect_true("CD3 | BV421.A" %in% result$markers)
  expect_true("CD4 | FITC.A" %in% result$markers)
  # BUV496.A | BUV496.A has marker == channel_root, so should be autofluorescence
  expect_true("BUV496.A | BUV496.A" %in% result$autofluorescence)
})

test_that("classify_csv_columns returns all categories even when empty", {
  result <- classify_csv_columns(c("CD3 | BV421.A", "CD4 | PE.A"))
  expect_true(all(c("index", "scatter", "markers", "gating",
                     "dimreduction", "autofluorescence") %in% names(result)))
  expect_length(result$index, 0)
  expect_length(result$scatter, 0)
  expect_length(result$markers, 2)
})

test_that("list_csvs lists CSV files with metadata", {
  tmp_dir <- file.path(tempdir(), paste0("test_list_csvs_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  dir.create(tmp_dir, recursive = TRUE)

  writeLines("a;b\n1;2", file.path(tmp_dir, "sample_A.csv"))
  writeLines("x;y\n3;4", file.path(tmp_dir, "sample_B.csv"))
  writeLines("not a csv", file.path(tmp_dir, "readme.txt"))

  result <- list_csvs(tmp_dir)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true(all(c("file_path", "file_name", "file_size_mb") %in% names(result)))
  expect_true("sample_A" %in% result$file_name)
  expect_true("sample_B" %in% result$file_name)
})

test_that("list_csvs errors on empty directory", {
  tmp_dir <- file.path(tempdir(), paste0("test_empty_csvs_", Sys.getpid()))
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  dir.create(tmp_dir, recursive = TRUE)

  expect_error(list_csvs(tmp_dir), "No .csv files found")
})

test_that("die stops execution", {
  # die() calls quit(), which we cannot test directly in testthat.
  # Instead we test that it would call message() and quit().
  # We can mock quit or just check the function exists and is callable.
  expect_true(is.function(die))
})

test_that("is_directory works for existing directories and files", {
  expect_true(is_directory(tempdir()))
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines("test", tmp)
  expect_false(is_directory(tmp))
  expect_false(is_directory("/nonexistent/path/abc"))
})

test_that("time_formatter formats seconds correctly", {
  expect_equal(time_formatter(30), "30s")
  expect_equal(time_formatter(120), "2min")
  expect_equal(time_formatter(7200), "2.0h")
  # Vectorized
  result <- time_formatter(c(5, 300, 3600))
  expect_equal(result, c("5s", "5min", "1.0h"))
})

test_that("get_file_size_mb returns numeric", {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(rep("some data here", 100), tmp)
  size <- get_file_size_mb(tmp)
  expect_type(size, "double")
  expect_true(size > 0)
})

test_that("dir_create creates multiple directories", {
  base <- file.path(tempdir(), paste0("test_dir_create_", Sys.getpid()))
  dirs <- file.path(base, c("dir_a", "dir_b", "dir_c"))
  on.exit(unlink(base, recursive = TRUE), add = TRUE)

  dir_create(dirs)
  for (d in dirs) expect_true(dir.exists(d))
  # Idempotent
  expect_silent(dir_create(dirs))
})


# =============================================================================
# 2. cluster_utils.R tests
# =============================================================================
context("cluster_utils.R")

test_that("clean_cluster_name removes AI_Cluster_ prefix", {
  expect_equal(clean_cluster_name("AI_Cluster_CD4TCs_1"), "CD4TCs_1")
})

test_that("clean_cluster_name removes (Clustering annotated) suffix", {
  expect_equal(clean_cluster_name("BCs (Clustering annotated)"), "BCs")
})

test_that("clean_cluster_name handles both prefix and suffix", {
  expect_equal(clean_cluster_name("AI_Cluster_NKCs (Clustering annotated)"), "NKCs")
})

test_that("clean_cluster_name returns Unknown for NA or empty", {
  expect_equal(clean_cluster_name(NA), "Unknown")
  expect_equal(clean_cluster_name(""), "Unknown")
})

test_that("clean_cluster_name returns unchanged string if no pattern matches", {
  expect_equal(clean_cluster_name("MyCluster"), "MyCluster")
})

test_that("extract_cluster_group extracts B cell group", {
  expect_equal(extract_cluster_group("BCs_1"), "BCs_1")
  expect_equal(extract_cluster_group("BCs"), "BCs")
  expect_equal(extract_cluster_group("AI_Cluster_BCs_3"), "BCs_3")
})

test_that("extract_cluster_group extracts CD4 T cell group", {
  expect_equal(extract_cluster_group("CD4TCs_5"), "CD4TCs_5")
  expect_equal(extract_cluster_group("CD4TCs"), "CD4TCs")
})

test_that("extract_cluster_group extracts CD8 T cell group", {
  expect_equal(extract_cluster_group("CD8TCs_2"), "CD8TCs_2")
})

test_that("extract_cluster_group extracts NK cell group", {
  expect_equal(extract_cluster_group("NKCs"), "NKCs")
  expect_equal(extract_cluster_group("NKCs_1"), "NKCs_1")
})

test_that("extract_cluster_group returns cleaned name for unrecognized patterns", {
  expect_equal(extract_cluster_group("SomeOtherCluster"), "SomeOtherCluster")
})

test_that("extract_cluster_group handles gdTCs", {
  expect_equal(extract_cluster_group("gdTCs_4"), "gdTCs_4")
  expect_equal(extract_cluster_group("gdTCs"), "gdTCs")
})

test_that("translate_labels translates using dictionary", {
  dict <- c("BCs" = "B Cells", "NKCs" = "NK Cells", "CD4TCs" = "CD4 T Cells")
  input <- c("BCs", "NKCs", "Unknown", "CD4TCs")
  result <- translate_labels(input, dict)
  expect_equal(result, c("B Cells", "NK Cells", "Unknown", "CD4 T Cells"))
})

test_that("translate_labels preserves unmatched labels", {
  dict <- c("A" = "Alpha")
  result <- translate_labels(c("A", "B", "C"), dict)
  expect_equal(result, c("Alpha", "B", "C"))
})

test_that("translate_labels handles empty input", {
  dict <- c("A" = "Alpha")
  result <- translate_labels(character(0), dict)
  expect_equal(result, character(0))
})

test_that("labels_from_binary_df extracts labels from binary data frame", {
  df <- data.frame(
    CD4TCs = c(1, 0, 1, 0),
    CD8TCs = c(0, 1, 0, 0),
    BCs    = c(0, 0, 1, 0)
  )
  result <- labels_from_binary_df(df)
  expect_type(result, "list")
  expect_equal(result[[1]], "CD4TCs")
  expect_equal(result[[2]], "CD8TCs")
  expect_equal(sort(result[[3]]), c("BCs", "CD4TCs"))
  expect_equal(result[[4]], character(0))
})

test_that("labels_from_binary_df handles all-zero rows", {
  df <- data.frame(X = c(0, 0), Y = c(0, 0))
  result <- labels_from_binary_df(df)
  expect_equal(result[[1]], character(0))
  expect_equal(result[[2]], character(0))
})

test_that("convert_binary_to_single_cluster assigns clusters correctly", {
  df <- data.frame(
    value  = c(10, 20, 30, 40),
    BCs    = c(1, 0, 0, 0),
    CD4TCs = c(0, 1, 0, 0),
    NKCs   = c(0, 0, 1, 0)
  )
  cluster_cols <- c("BCs", "CD4TCs", "NKCs")
  result <- convert_binary_to_single_cluster(df, cluster_cols)

  expect_true("cluster_assignment" %in% names(result))
  expect_equal(result$cluster_assignment[1], "BCs")
  expect_equal(result$cluster_assignment[2], "CD4TCs")
  expect_equal(result$cluster_assignment[3], "NKCs")
  expect_equal(result$cluster_assignment[4], "Unassigned")
})

test_that("create_cluster_id truncates long names", {
  long_name <- paste(rep("A", 50), collapse = "")
  result <- create_cluster_id(long_name, max_length = 20)
  expect_true(nchar(result) <= 20)
})

test_that("create_cluster_id sanitizes special characters", {
  result <- create_cluster_id("CD4+ T Cells (activated)")
  expect_false(grepl("[^A-Za-z0-9_-]", result))
})

test_that("flatten_label_list selects deepest label", {
  labels_list <- list(
    c("root", "TCells", "CD4TCs"),
    c("root", "BCs"),
    character(0)
  )
  gating_tree <- c("root" = 0, "TCells" = 1, "CD4TCs" = 2, "BCs" = 1)
  result <- flatten_label_list(labels_list, gating_tree)
  expect_equal(result, c("CD4TCs", "BCs", "root"))
})

test_that("labels_from_binary_dt works with data.table", {
  skip_if_not_installed("data.table")
  dt <- data.table::data.table(
    index  = 1:3,
    CD4TCs = c(1, 0, 0),
    BCs    = c(0, 1, 0),
    NKCs   = c(0, 0, 1)
  )
  result <- labels_from_binary_dt(dt)
  expect_type(result, "list")
  expect_equal(result[[1]], "CD4TCs")
  expect_equal(result[[2]], "BCs")
  expect_equal(result[[3]], "NKCs")
})


# =============================================================================
# 3. stats_utils.R tests
# =============================================================================
context("stats_utils.R")

test_that("calculate_enrichment returns correct Fisher's test results", {
  cells_data <- data.frame(
    cluster_X = c(rep(1, 30), rep(0, 70), rep(1, 10), rep(0, 90)),
    condition = c(rep("A", 100), rep("B", 100))
  )
  result <- calculate_enrichment("cluster_X", cells_data, c("A", "B"))

  expect_s3_class(result, "tbl_df")
  expect_true("p_value" %in% names(result))
  expect_true("odds_ratio" %in% names(result))
  expect_equal(result$in_A, 30)
  expect_equal(result$not_in_A, 70)
  expect_equal(result$in_B, 10)
  expect_equal(result$not_in_B, 90)
  expect_equal(result$method, "fisher")
  # With 30% in A vs 10% in B, should be significant
  expect_true(result$p_value < 0.05)
  expect_true(result$odds_ratio > 1)
})

test_that("calculate_enrichment handles balanced data", {
  cells_data <- data.frame(
    cluster_X = c(rep(1, 50), rep(0, 50), rep(1, 50), rep(0, 50)),
    condition = c(rep("A", 100), rep("B", 100))
  )
  result <- calculate_enrichment("cluster_X", cells_data, c("A", "B"))
  expect_equal(result$method, "fisher")
  # Balanced -> p should be 1.0 (or very close)
  expect_true(result$p_value > 0.9)
})

test_that("calculate_enrichment handles degenerate case (empty condition)", {
  cells_data <- data.frame(
    cluster_X = c(1, 0, 1),
    condition = c("A", "A", "A")
  )
  result <- calculate_enrichment("cluster_X", cells_data, c("A", "B"))
  # condB has 0 cells, so cc + d == 0 -> fisher_failed
  expect_equal(result$method, "fisher_failed")
  expect_true(is.na(result$p_value))
})

test_that("run_pca returns correct structure and dimensions", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 50, ncol = 4,
                dimnames = list(NULL, paste0("M", 1:4)))

  result <- run_pca(mat, n_components = 2)

  expect_type(result, "list")
  expect_true(all(c("coordinates", "variance_explained", "loadings",
                     "pca_object") %in% names(result)))
  expect_equal(nrow(result$coordinates), 50)
  expect_equal(ncol(result$coordinates), 2)
  expect_length(result$variance_explained, 2)
  expect_true(all(result$variance_explained >= 0))
  expect_true(all(result$variance_explained <= 1))
  expect_equal(nrow(result$loadings), 4)
  expect_equal(ncol(result$loadings), 2)
})

test_that("run_pca with 3 components", {
  set.seed(123)
  mat <- matrix(rnorm(500), nrow = 100, ncol = 5,
                dimnames = list(NULL, paste0("V", 1:5)))
  result <- run_pca(mat, n_components = 3, scale = FALSE)
  expect_equal(ncol(result$coordinates), 3)
  expect_length(result$variance_explained, 3)
})

test_that("run_umap returns correct dimensions", {
  skip_if_not_installed("uwot")
  set.seed(42)
  # Need enough data points for UMAP (n > n_neighbors)
  mat <- matrix(rnorm(400), nrow = 40, ncol = 10)

  result <- run_umap(mat, n_neighbors = 10, n_components = 2, seed = 42)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 40)
  expect_equal(ncol(result), 2)
})

test_that("run_umap with 3 components", {
  skip_if_not_installed("uwot")
  set.seed(42)
  mat <- matrix(rnorm(400), nrow = 40, ncol = 10)
  result <- run_umap(mat, n_neighbors = 10, n_components = 3, seed = 42)
  expect_equal(ncol(result), 3)
})

test_that("get_top_features selects by variance", {
  mat <- matrix(c(
    1, 1, 1, 1,   # low variance
    1, 2, 3, 4,   # medium variance
    1, 10, 100, 1000  # high variance
  ), nrow = 4, ncol = 3, dimnames = list(NULL, c("low", "med", "high")))

  result <- get_top_features(mat, n_top = 2, method = "variance")
  expect_length(result, 2)
  expect_equal(result[1], "high")
})

test_that("get_top_features selects by mean", {
  mat <- matrix(c(
    1, 1, 1, 1,      # low mean
    10, 10, 10, 10,   # medium mean
    100, 100, 100, 100  # high mean
  ), nrow = 4, ncol = 3, dimnames = list(NULL, c("low", "med", "high")))

  result <- get_top_features(mat, n_top = 1, method = "mean")
  expect_equal(result, "high")
})

test_that("get_top_features handles n_top > ncol", {
  mat <- matrix(rnorm(12), nrow = 4, ncol = 3,
                dimnames = list(NULL, c("a", "b", "c")))
  result <- get_top_features(mat, n_top = 10)
  expect_length(result, 3)
})

test_that("get_top_features defaults to variance for unknown method", {
  mat <- matrix(c(1, 1, 1, 10, 1, 100), nrow = 3, ncol = 2,
                dimnames = list(NULL, c("lo", "hi")))
  result <- get_top_features(mat, n_top = 1, method = "invalid_method")
  # Should default to variance, which is higher for "hi"
  expect_equal(result, "hi")
})

test_that("check_overdispersion returns expected structure", {
  # Create a simple Poisson GLM to test
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- rpois(n, lambda = exp(1 + 0.5 * x))
  model <- glm(y ~ x, family = poisson)

  result <- check_overdispersion(model)

  expect_type(result, "list")
  expect_true(all(c("ratio", "p_value", "is_overdispersed") %in% names(result)))
  expect_type(result$ratio, "double")
})

test_that("check_overdispersion handles errors gracefully", {
  # Pass something that is not a model
  result <- check_overdispersion("not_a_model")
  expect_true(is.na(result$ratio))
  expect_true(is.na(result$p_value))
  expect_true(is.na(result$is_overdispersed))
})

test_that("calculate_statistics computes per-cluster per-marker stats", {
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("rlang")

  expr_data <- data.frame(
    cluster = rep(c("A", "B"), each = 20),
    marker1 = c(rnorm(20, mean = 5), rnorm(20, mean = 10)),
    marker2 = c(rnorm(20, mean = 1), rnorm(20, mean = 2))
  )

  result <- calculate_statistics(expr_data, "cluster", c("marker1", "marker2"))
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("mean", "median", "sd", "q25", "q75", "n") %in% names(result)))
  # Two clusters x two markers = 4 rows

  expect_equal(nrow(result), 4)
  # Cluster A marker1 should have mean around 5
  a_m1 <- result[result$cluster == "A" & result$marker == "marker1", ]
  expect_true(abs(a_m1$mean - 5) < 2)
})


# =============================================================================
# 4. flowframe_utils.R tests
# =============================================================================
context("flowframe_utils.R")

test_that("from_mtx_to_ffr converts matrix to flowFrame", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat <- matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
                nrow = 3, ncol = 2)

  ff <- from_mtx_to_ffr(mat, marker_names = c("CD3", "CD4"))

  expect_s4_class(ff, "flowFrame")
  expect_equal(nrow(flowCore::exprs(ff)), 3)
  expect_equal(ncol(flowCore::exprs(ff)), 2)
  expect_equal(unname(colnames(flowCore::exprs(ff))), c("CD3", "CD4"))
})

test_that("from_mtx_to_ffr preserves expression values", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat <- matrix(c(10.5, 20.3, 30.1, 40.7, 50.2, 60.9),
                nrow = 3, ncol = 2)
  colnames(mat) <- c("M1", "M2")

  ff <- from_mtx_to_ffr(mat)
  expr <- flowCore::exprs(ff)
  # flowFrame stores values as single-precision floats, so use a loose tolerance.
  expect_equal(unname(expr[1, 1]), 10.5, tolerance = 1e-5)
  expect_equal(unname(expr[3, 2]), 60.9, tolerance = 1e-5)
})

test_that("build_combined_expr_matrix combines multiple flowFrames", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat1 <- matrix(rnorm(20), nrow = 10, ncol = 2)
  mat2 <- matrix(rnorm(30), nrow = 15, ncol = 2)

  ff1 <- from_mtx_to_ffr(mat1, marker_names = c("CD3", "CD4"))
  ff2 <- from_mtx_to_ffr(mat2, marker_names = c("CD3", "CD4"))

  result <- build_combined_expr_matrix(
    flow_frames_list = list(ff1, ff2),
    sample_names = c("Sample_A", "Sample_B"),
    markers_of_interest = c("CD3", "CD4")
  )

  expect_type(result, "list")
  expect_true(all(c("combined_matrix", "sample_labels") %in% names(result)))
  expect_equal(nrow(result$combined_matrix), 25)  # 10 + 15
  expect_equal(ncol(result$combined_matrix), 2)
  expect_length(result$sample_labels, 25)
  expect_equal(sum(result$sample_labels == "Sample_A"), 10)
  expect_equal(sum(result$sample_labels == "Sample_B"), 15)
})

test_that("explore_data extracts marker information", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  ff <- from_mtx_to_ffr(mat, marker_names = c("CD3", "CD4", "CD8"))

  result <- explore_data(ff)
  expect_type(result, "list")
  expect_true(all(c("dyes", "markers", "mcols") %in% names(result)))
  expect_length(result$dyes, 3)
  expect_length(result$markers, 3)
})

test_that("generate_cellID creates correct IDs", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat <- matrix(rnorm(10), nrow = 5, ncol = 2)
  ff <- from_mtx_to_ffr(mat, marker_names = c("M1", "M2"))

  ids <- generate_cellID(ff, "Sample1")
  expect_length(ids, 5)
  expect_equal(ids[1], "Sample1_cell_1")
  expect_equal(ids[5], "Sample1_cell_5")
})

test_that("from_ffl_to_fSet converts list to flowSet", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat1 <- matrix(rnorm(10), nrow = 5, ncol = 2)
  mat2 <- matrix(rnorm(10), nrow = 5, ncol = 2)
  ff1 <- from_mtx_to_ffr(mat1, marker_names = c("CD3", "CD4"))
  ff2 <- from_mtx_to_ffr(mat2, marker_names = c("CD3", "CD4"))

  fs <- from_ffl_to_fSet(list(ff1, ff2), sample_names = c("S1", "S2"))
  expect_s4_class(fs, "flowSet")
  expect_equal(length(fs), 2)
})


# =============================================================================
# 5. plot_utils.R tests
# =============================================================================
context("plot_utils.R")

test_that("create_safe_filename removes special characters", {
  expect_equal(create_safe_filename("hello world!"), "hello_world")
  expect_equal(create_safe_filename("file/with\\slashes"), "file_with_slashes")
  expect_equal(create_safe_filename("dots.are.ok"), "dots.are.ok")
  expect_equal(create_safe_filename("under_scores-and-dashes"), "under_scores-and-dashes")
})

test_that("create_safe_filename collapses multiple underscores", {
  expect_equal(create_safe_filename("a   b   c"), "a_b_c")
  expect_equal(create_safe_filename("a!!!b"), "a_b")
})

test_that("create_safe_filename respects max_length", {
  long_name <- paste(rep("a", 200), collapse = "")
  result <- create_safe_filename(long_name, max_length = 50)
  expect_true(nchar(result) <= 50)
})

test_that("create_safe_filename strips leading/trailing underscores", {
  expect_equal(create_safe_filename("_leading"), "leading")
  expect_equal(create_safe_filename("trailing_"), "trailing")
  expect_equal(create_safe_filename("!special!"), "special")
})

test_that("create_safe_filename handles empty-like input", {
  result <- create_safe_filename("!!!")
  # All special chars removed, then underscores collapsed and stripped
  expect_equal(result, "")
})

test_that("create_compact_names truncates long names", {
  names_vec <- c("ShortName", "ThisIsAVeryLongNameThatShouldBeTruncated")
  result <- create_compact_names(names_vec, max_length = 20)
  expect_equal(result[["ShortName"]], "ShortName")
  expect_true(nchar(result[[2]]) == 20)
  expect_true(grepl("\\.\\.\\.$", result[[2]]))
})

test_that("create_compact_names preserves names within limit", {
  result <- create_compact_names(c("abc", "def"), max_length = 10)
  expect_equal(unname(result), c("abc", "def"))
})

test_that("create_compact_names handles exact boundary length", {
  result <- create_compact_names("12345", max_length = 5)
  expect_equal(unname(result), "12345")
})

test_that("create_compact_names adds ellipsis at exact boundary + 1", {
  result <- create_compact_names("123456", max_length = 5)
  expect_equal(nchar(unname(result)), 5)
  expect_true(grepl("\\.\\.\\.$", unname(result)))
  # "12..."
  expect_equal(unname(result), "12...")
})


# =============================================================================
# 6. memory_utils.R tests
# =============================================================================
context("memory_utils.R")

test_that("get_available_memory_gb returns numeric or NA", {
  result <- get_available_memory_gb()
  expect_true(is.numeric(result))
  # On Linux, should return a positive number; on unknown OS, NA
  if (!is.na(result)) {
    expect_true(result > 0)
  }
})

test_that("check_memory returns logical without error", {
  result <- check_memory()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("check_memory with custom threshold", {
  # Very low threshold should return FALSE (memory is "high")
  result_low <- check_memory(warn_threshold = 0.01)
  # On most systems, memory usage will be > 1%, so this should return FALSE
  # unless the system cannot determine memory at all (returns TRUE)
  expect_type(result_low, "logical")

  # Very high threshold should return TRUE
  result_high <- check_memory(warn_threshold = 0.999)
  expect_type(result_high, "logical")
})

test_that(".get_os returns a valid OS string", {
  os <- .get_os()
  expect_type(os, "character")
  expect_true(os %in% c("Linux", "Darwin", "Windows"))
})

test_that(".query_system_memory returns named vector or NULL", {
  mem <- .query_system_memory()
  if (!is.null(mem)) {
    expect_true(is.numeric(mem))
    expect_true(all(c("total_mb", "used_mb", "available_mb") %in% names(mem)))
    expect_true(mem["total_mb"] > 0)
  }
})

test_that(".count_lines counts file lines correctly", {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("line1", "line2", "line3", "line4", "line5"), tmp)
  count <- .count_lines(tmp)
  expect_equal(count, 5)
})

test_that("read_csv_chunked reads small files directly", {
  skip_if_not_installed("data.table")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  data.table::fwrite(data.table::data.table(a = 1:10, b = 11:20), tmp)

  result <- read_csv_chunked(tmp, max_memory_gb = 100)
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 10)
  expect_true("a" %in% names(result))
})


# =============================================================================
# 7. %||% polyfill (defined in clustering_dispatch.R)
# =============================================================================
test_that("%||% operator works correctly", {
  expect_equal(NULL %||% 42, 42)
  expect_equal(5 %||% 42, 5)
  expect_equal("a" %||% "b", "a")
  expect_equal(NULL %||% NULL, NULL)
})


# =============================================================================
# 8. Dispatch function tests
# =============================================================================
context("clustering_dispatch.R")

test_that("run_clustering rejects invalid method", {
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  expect_error(
    run_clustering(mat, method = "invalid_method"),
    "Unknown clustering method"
  )
})

test_that("run_clustering rejects non-matrix input", {
  expect_error(
    run_clustering(data.frame(a = 1:5, b = 6:10), method = "leiden"),
    "is.matrix"
  )
})

test_that("run_clustering rejects non-numeric matrix", {
  mat <- matrix(letters[1:12], nrow = 3, ncol = 4)
  expect_error(
    run_clustering(mat, method = "leiden"),
    "is.numeric"
  )
})

test_that("run_clustering valid methods list is complete", {
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  # Check that the error message includes all valid method names
  err <- tryCatch(run_clustering(mat, method = "xyz"), error = conditionMessage)
  expect_true(grepl("leiden", err))
  expect_true(grepl("flowsom", err))
  expect_true(grepl("phenograph", err))
  expect_true(grepl("dbscan", err))
})

context("dimred_dispatch.R")

test_that("run_dimred rejects invalid method", {
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  expect_error(
    run_dimred(mat, method = "invalid_method"),
    "Unknown dimred method"
  )
})

test_that("run_dimred rejects non-matrix input", {
  expect_error(
    run_dimred(data.frame(a = 1:5, b = 6:10), method = "umap"),
    "is.matrix"
  )
})

test_that("run_dimred rejects non-numeric matrix", {
  mat <- matrix(letters[1:12], nrow = 3, ncol = 4)
  expect_error(
    run_dimred(mat, method = "umap"),
    "is.numeric"
  )
})

test_that("run_dimred valid methods list is complete", {
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  err <- tryCatch(run_dimred(mat, method = "xyz"), error = conditionMessage)
  expect_true(grepl("umap", err))
  expect_true(grepl("tsne", err))
  expect_true(grepl("diffusion_map", err))
})

context("batch_correction_dispatch.R")

test_that("run_batch_correction_dispatch rejects invalid method", {
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  expect_error(
    run_batch_correction_dispatch(mat, method = "invalid_method",
                                   batch_labels = rep("A", 20)),
    "Unknown batch correction method"
  )
})

test_that("run_batch_correction_dispatch 'none' passes through matrix data", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5,
                dimnames = list(NULL, paste0("M", 1:5)))
  sample_labels <- rep(c("S1", "S2"), each = 10)
  batch_labels <- rep(c("B1", "B2"), each = 10)

  result <- run_batch_correction_dispatch(
    data = mat,
    method = "none",
    batch_labels = batch_labels,
    sample_names = sample_labels
  )

  expect_type(result, "list")
  expect_true(all(c("corrected_matrix", "sample_labels", "batch_labels",
                     "markers", "corrected_frames") %in% names(result)))
  # "none" should pass through the matrix unchanged
  expect_identical(result$corrected_matrix, mat)
  expect_equal(result$sample_labels, sample_labels)
  expect_equal(result$batch_labels, batch_labels)
  expect_equal(result$markers, paste0("M", 1:5))
  expect_null(result$corrected_frames)
})

test_that("run_batch_correction_dispatch 'none' passes through flowFrame data", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat1 <- matrix(rnorm(20), nrow = 10, ncol = 2)
  mat2 <- matrix(rnorm(20), nrow = 10, ncol = 2)
  ff1 <- from_mtx_to_ffr(mat1, marker_names = c("CD3", "CD4"))
  ff2 <- from_mtx_to_ffr(mat2, marker_names = c("CD3", "CD4"))

  result <- run_batch_correction_dispatch(
    data = list(ff1, ff2),
    method = "none",
    batch_labels = c("B1", "B2"),
    markers_of_interest = c("CD3", "CD4"),
    sample_names = c("S1", "S2")
  )

  expect_type(result, "list")
  expect_equal(nrow(result$corrected_matrix), 20)  # 10 + 10
  expect_equal(ncol(result$corrected_matrix), 2)
  expect_length(result$sample_labels, 20)
  expect_length(result$batch_labels, 20)
  expect_null(result$corrected_frames)
})

test_that("run_batch_correction_dispatch valid methods list is complete", {
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
  err <- tryCatch(
    run_batch_correction_dispatch(mat, method = "xyz",
                                   batch_labels = rep("A", 20)),
    error = conditionMessage
  )
  expect_true(grepl("harmony", err))
  expect_true(grepl("cytonorm", err))
  expect_true(grepl("combat", err))
  expect_true(grepl("none", err))
})

test_that("run_batch_correction_dispatch method is case-insensitive", {
  mat <- matrix(rnorm(100), nrow = 20, ncol = 5,
                dimnames = list(NULL, paste0("M", 1:5)))

  result <- run_batch_correction_dispatch(
    data = mat, method = "NONE",
    batch_labels = rep("B1", 20),
    sample_names = rep("S1", 20)
  )
  expect_identical(result$corrected_matrix, mat)
})


# =============================================================================
# 9. Additional edge-case and integration tests
# =============================================================================
context("Edge cases and integration")

test_that("save_and_cleanup writes semicolon-separated CSV", {
  skip_if_not_installed("data.table")

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  dt <- data.table::data.table(x = 1:3, y = c("a", "b", "c"))

  save_and_cleanup(dt, tmp)

  expect_true(file.exists(tmp))
  content <- readLines(tmp, n = 1)
  expect_true(grepl(";", content))
})

test_that("classify_csv_columns with only marker columns", {
  result <- classify_csv_columns(c("CD3 | BV421.A", "CD4 | FITC.A", "CD8 | PE.A"))
  expect_length(result$markers, 3)
  expect_length(result$index, 0)
  expect_length(result$scatter, 0)
  expect_length(result$gating, 0)
})

test_that("classify_csv_columns case-insensitive index detection", {
  result_lower <- classify_csv_columns(c("index", "marker1"))
  result_upper <- classify_csv_columns(c("Index", "marker1"))
  expect_true("index" %in% result_lower$index)
  expect_true("Index" %in% result_upper$index)
})

test_that("from_mtx_to_ffr handles single-column matrix", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  mat <- matrix(1:5, ncol = 1)
  ff <- from_mtx_to_ffr(mat, marker_names = c("CD3"))
  expect_equal(ncol(flowCore::exprs(ff)), 1)
  expect_equal(nrow(flowCore::exprs(ff)), 5)
})

test_that("run_pca handles constant columns gracefully", {
  mat <- matrix(c(rep(5, 20), rnorm(20)), nrow = 20, ncol = 2,
                dimnames = list(NULL, c("const", "varied")))
  # scale = TRUE with a constant column would produce NaN,
  # but prcomp handles it. Test with scale = FALSE.
  result <- run_pca(mat, n_components = 1, scale = FALSE)
  expect_equal(nrow(result$coordinates), 20)
  expect_equal(ncol(result$coordinates), 1)
})

test_that("extract_cluster_group handles NKTCells pattern", {
  expect_equal(extract_cluster_group("NKTCells"), "NKTCells")
  expect_equal(extract_cluster_group("NKTCells_1"), "NKTCells_1")
})

test_that("extract_cluster_group handles CD3- pattern", {
  expect_equal(extract_cluster_group("CD3-_1"), "CD3-_1")
  expect_equal(extract_cluster_group("CD3-"), "CD3-")
})

test_that("extract_cluster_group handles LymphoidCells", {
  expect_equal(extract_cluster_group("LymphoidCells"), "LymphoidCells")
  expect_equal(extract_cluster_group("LymphoidCells_2"), "LymphoidCells_2")
})

test_that("extract_cluster_group handles TCells", {
  expect_equal(extract_cluster_group("TCells_7"), "TCells_7")
  expect_equal(extract_cluster_group("TCells"), "TCells")
})

test_that("calculate_enrichment returns CI bounds", {
  cells_data <- data.frame(
    cluster_X = c(rep(1, 40), rep(0, 60), rep(1, 10), rep(0, 90)),
    condition = c(rep("A", 100), rep("B", 100))
  )
  result <- calculate_enrichment("cluster_X", cells_data, c("A", "B"))
  expect_true(!is.na(result$ci_lower))
  expect_true(!is.na(result$ci_upper))
  expect_true(result$ci_lower < result$ci_upper)
  expect_true(result$ci_lower < result$odds_ratio)
  expect_true(result$odds_ratio < result$ci_upper)
})


# =============================================================================
# Summary
# =============================================================================
cat("\n=== R_shared test suite complete ===\n")
