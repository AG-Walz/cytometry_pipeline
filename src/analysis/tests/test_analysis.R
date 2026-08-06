library(testthat)

# ---------------------------------------------------------------------------
# Helper: resolve paths relative to the project root so tests work both
# when run via testthat::test_dir() and from the repo root.
# ---------------------------------------------------------------------------
find_shared_dir <- function() {
  candidates <- c(
    file.path(dirname(getwd()), "src", "R_shared"),
    file.path(getwd(), "..", "..", "R_shared"),
    file.path(getwd(), "..", "R_shared")
  )
  for (d in candidates) {
    if (dir.exists(d)) return(normalizePath(d))
  }
  stop("Cannot locate src/R_shared from working directory: ", getwd())
}

shared_dir   <- find_shared_dir()
analysis_dir <- file.path(dirname(shared_dir), "analysis")

# Source base shared utilities
source(file.path(shared_dir, "io_utils.R"))

# ============================= calculate_enrichment() =======================
test_that("calculate_enrichment returns correct Fisher's exact test results", {
  source(file.path(shared_dir, "stats_utils.R"))

  # Build a data frame with known 2x2 table:
  #                condA   condB
  # in_cluster       20       5
  # not_in_cluster   30      45
  cells <- data.frame(
    cluster_col = c(rep(1, 20), rep(0, 30), rep(1, 5), rep(0, 45)),
    condition   = c(rep("A", 50), rep("B", 50)),
    stringsAsFactors = FALSE
  )

  result <- calculate_enrichment("cluster_col", cells, c("A", "B"))

  expect_s3_class(result, "tbl_df")
  expect_equal(result$in_A, 20)
  expect_equal(result$not_in_A, 30)
  expect_equal(result$in_B, 5)
  expect_equal(result$not_in_B, 45)
  expect_equal(result$method, "fisher")
  expect_true(!is.na(result$odds_ratio))
  expect_true(result$odds_ratio > 1)                  # cluster enriched in A
  expect_true(!is.na(result$p_value))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(!is.na(result$ci_lower))
  expect_true(!is.na(result$ci_upper))
  expect_true(result$ci_lower < result$ci_upper)
})

test_that("calculate_enrichment handles zero-count gracefully", {
  source(file.path(shared_dir, "stats_utils.R"))

  # All cells in cluster for condA, none in condB
  cells <- data.frame(
    cluster_col = c(rep(1, 10), rep(0, 0), rep(0, 10)),
    condition   = c(rep("A", 10), rep("B", 10)),
    stringsAsFactors = FALSE
  )

  result <- calculate_enrichment("cluster_col", cells, c("A", "B"))
  expect_s3_class(result, "tbl_df")
  # Should still return a tibble (either with fisher result or fisher_failed)
  expect_true(result$method %in% c("fisher", "fisher_failed"))
})

test_that("calculate_enrichment returns failed when one condition is empty", {
  source(file.path(shared_dir, "stats_utils.R"))

  # No cells in condition B at all
  cells <- data.frame(
    cluster_col = c(rep(1, 5), rep(0, 5)),
    condition   = rep("A", 10),
    stringsAsFactors = FALSE
  )

  result <- calculate_enrichment("cluster_col", cells, c("A", "B"))
  expect_equal(result$method, "fisher_failed")
  expect_true(is.na(result$odds_ratio))
})

# ============================= split_by_cluster() ===========================
test_that("split_by_cluster produces per-cluster CSV files", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("progress")

  source(file.path(shared_dir, "gui_utils.R"))
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "memory_utils.R"))
  source(file.path(analysis_dir, "cluster_splitting.R"))

  # Create synthetic per-cell CSVs with binary cluster columns
  tmp_in  <- tempfile("split_input_")
  tmp_out <- tempfile("split_output_")
  dir.create(tmp_in, recursive = TRUE)
  on.exit(unlink(c(tmp_in, tmp_out), recursive = TRUE), add = TRUE)

  set.seed(30)
  for (s in c("sample1", "sample2")) {
    n <- 50
    dt <- data.table::data.table(
      marker_A = rnorm(n),
      marker_B = rnorm(n),
      `Cluster_1 (Clustering annotated)` = c(rep(1L, 20), rep(0L, 30)),
      `Cluster_2 (Clustering annotated)` = c(rep(0L, 20), rep(1L, 15), rep(0L, 15))
    )
    data.table::fwrite(dt, file.path(tmp_in, paste0(s, ".csv")), sep = ";")
  }

  result <- split_by_cluster(
    input_dir = tmp_in,
    output_dir = tmp_out,
    cluster_regex = "Clustering annotated",
    markers_of_interest = c("marker_A", "marker_B"),
    min_cells_per_cluster = 5,
    overwrite_existing = TRUE
  )

  expect_type(result, "list")
  expect_true("tracking_matrix" %in% names(result))
  expect_true("cluster_summary" %in% names(result))
  expect_true("output_files" %in% names(result))

  # Output files should exist
  expect_true(length(result$output_files) > 0)
  expect_true(all(file.exists(result$output_files)))

  # Tracking matrix should have rows for each sample
  tm <- result$tracking_matrix
  expect_true(nrow(tm) >= 1)
})

# ======================== build_tracking_matrix (via split_by_cluster) ======
test_that("tracking_matrix from split_by_cluster is wide format with correct dims", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("progress")

  source(file.path(shared_dir, "gui_utils.R"))
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "memory_utils.R"))
  source(file.path(analysis_dir, "cluster_splitting.R"))

  tmp_in  <- tempfile("tm_input_")
  tmp_out <- tempfile("tm_output_")
  dir.create(tmp_in, recursive = TRUE)
  on.exit(unlink(c(tmp_in, tmp_out), recursive = TRUE), add = TRUE)

  # 3 samples, 2 clusters
  for (s in c("s1", "s2", "s3")) {
    dt <- data.table::data.table(
      m1 = rnorm(40),
      `Cluster_A (Clustering annotated)` = c(rep(1L, 15), rep(0L, 25)),
      `Cluster_B (Clustering annotated)` = c(rep(0L, 15), rep(1L, 10), rep(0L, 15))
    )
    data.table::fwrite(dt, file.path(tmp_in, paste0(s, ".csv")), sep = ";")
  }

  result <- split_by_cluster(
    input_dir = tmp_in, output_dir = tmp_out,
    cluster_regex = "Clustering annotated",
    markers_of_interest = "m1",
    min_cells_per_cluster = 5
  )

  tm <- result$tracking_matrix
  # Wide format: one row per sample, columns for file + cluster names
  expect_true("file" %in% names(tm) || ncol(tm) >= 2)
  # 3 samples -> 3 rows
  expect_equal(nrow(tm), 3)
  # At least 2 cluster columns + id column
  expect_true(ncol(tm) >= 3)
})

# ======================== tracking_matrix_to_long() =========================
test_that("tracking_matrix_to_long pivots correctly", {
  skip_if_not_installed("tidyr")

  source(file.path(shared_dir, "gui_utils.R"))
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "memory_utils.R"))
  source(file.path(analysis_dir, "cluster_splitting.R"))

  wide <- data.frame(
    file = c("s1", "s2", "s3"),
    Cluster_A = c(10, 20, 15),
    Cluster_B = c(5, 8, 12),
    stringsAsFactors = FALSE
  )

  long <- tracking_matrix_to_long(wide, id_col = "file",
                                    names_to = "cluster", values_to = "count")

  expect_s3_class(long, "data.frame")
  expect_equal(nrow(long), 6)                         # 3 samples x 2 clusters
  expect_true("file" %in% names(long))
  expect_true("cluster" %in% names(long))
  expect_true("count" %in% names(long))

  # Check a specific value
  val <- long$count[long$file == "s2" & long$cluster == "Cluster_A"]
  expect_equal(val, 20)
})

test_that("tracking_matrix_to_long preserves all rows", {
  skip_if_not_installed("tidyr")

  source(file.path(shared_dir, "gui_utils.R"))
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "memory_utils.R"))
  source(file.path(analysis_dir, "cluster_splitting.R"))

  wide <- data.frame(
    file = paste0("sample_", 1:5),
    c1 = 1:5,
    c2 = 6:10,
    c3 = 11:15
  )

  long <- tracking_matrix_to_long(wide)
  expect_equal(nrow(long), 15)                        # 5 x 3
  expect_equal(sum(long$count), sum(1:15))
})

# ======================== calculate_cell_ratios() ===========================
test_that("calculate_cell_ratios computes correct ratio", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "cell_ratio.R"))

  # Create a temporary abundance CSV
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv), add = TRUE)

  abundance <- data.table::data.table(
    sample   = c("s1", "s2", "s3"),
    CD4      = c(100, 200, 50),
    CD8      = c(50, 100, 150),
    NK       = c(30, 40, 20)
  )
  data.table::fwrite(abundance, tmp_csv, sep = ";")

  result <- calculate_cell_ratios(
    abundance_csv = tmp_csv,
    numerator_clusters = "CD4",
    denominator_clusters = "CD8"
  )

  expect_s3_class(result, "data.frame")
  expect_true("ratio" %in% names(result))
  expect_true("log2_ratio" %in% names(result))
  expect_true("sample" %in% names(result))
  expect_equal(nrow(result), 3)

  # s1: CD4=100, CD8=50 => ratio = 2.0
  expect_equal(result$ratio[result$sample == "s1"], 2.0)

  # s2: CD4=200, CD8=100 => ratio = 2.0
  expect_equal(result$ratio[result$sample == "s2"], 2.0)

  # s3: CD4=50, CD8=150 => ratio = 50/150 = 1/3
  expect_equal(result$ratio[result$sample == "s3"], 50 / 150, tolerance = 1e-6)
})

test_that("calculate_cell_ratios handles multi-column numerator/denominator", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "cell_ratio.R"))

  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv), add = TRUE)

  abundance <- data.table::data.table(
    sample = c("s1"),
    A = c(10),
    B = c(20),
    C = c(5),
    D = c(15)
  )
  data.table::fwrite(abundance, tmp_csv, sep = ";")

  result <- calculate_cell_ratios(
    abundance_csv = tmp_csv,
    numerator_clusters = c("A", "B"),
    denominator_clusters = c("C", "D")
  )

  # numerator = 10 + 20 = 30, denominator = 5 + 15 = 20, ratio = 1.5
  expect_equal(result$ratio, 1.5)
})

test_that("calculate_cell_ratios errors on missing columns", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "cell_ratio.R"))

  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv), add = TRUE)

  abundance <- data.table::data.table(sample = "s1", A = 10)
  data.table::fwrite(abundance, tmp_csv, sep = ";")

  expect_error(
    calculate_cell_ratios(tmp_csv, numerator_clusters = "MISSING",
                           denominator_clusters = "A"),
    "No numerator clusters found"
  )
})

# ======================== compute_effect_sizes() ============================
test_that("compute_effect_sizes returns Cohen's d for known groups", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "effect_sizes.R"))

  set.seed(40)
  n <- 50
  expr_data <- data.table::data.table(
    marker1   = c(rnorm(n, mean = 0, sd = 1), rnorm(n, mean = 2, sd = 1)),
    marker2   = c(rnorm(n, mean = 5, sd = 1), rnorm(n, mean = 5, sd = 1)),
    condition = rep(c("A", "B"), each = n)
  )

  result <- compute_effect_sizes(expr_data,
                                  condition_col = "condition",
                                  value_cols = c("marker1", "marker2"),
                                  methods = "cohens_d")

  expect_s3_class(result, "data.frame")
  expect_true("marker" %in% names(result))
  expect_true("effect_size" %in% names(result))
  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true("method" %in% names(result))
  expect_true("interpretation" %in% names(result))

  # marker1: means differ by 2 with sd=1 => Cohen's d ~ -2 (large)
  m1 <- result[result$marker == "marker1", ]
  expect_equal(m1$method, "cohens_d")
  expect_true(abs(m1$effect_size) > 1.0)               # large effect
  expect_equal(m1$interpretation, "large")

  # marker2: means are same => Cohen's d ~ 0 (negligible)
  m2 <- result[result$marker == "marker2", ]
  expect_true(abs(m2$effect_size) < 0.5)
})

test_that("compute_effect_sizes computes Cliff's delta", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "effect_sizes.R"))

  set.seed(41)
  expr_data <- data.table::data.table(
    m1        = c(rnorm(30, 0), rnorm(30, 3)),
    condition = rep(c("X", "Y"), each = 30)
  )

  result <- compute_effect_sizes(expr_data, "condition", "m1",
                                  methods = "cliffs_delta")

  expect_equal(nrow(result), 1)
  expect_equal(result$method, "cliffs_delta")
  expect_true(abs(result$effect_size) > 0.3)           # should be substantial
  expect_true(!is.na(result$ci_lower))
  expect_true(!is.na(result$ci_upper))
})

test_that("compute_effect_sizes computes Hedges' g", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "effect_sizes.R"))

  set.seed(42)
  expr_data <- data.table::data.table(
    m1        = c(rnorm(20, 0), rnorm(20, 1)),
    condition = rep(c("G1", "G2"), each = 20)
  )

  result <- compute_effect_sizes(expr_data, "condition", "m1",
                                  methods = "hedges_g")

  expect_equal(result$method, "hedges_g")
  # Hedges' g is slightly smaller than Cohen's d (bias correction)
  expect_true(is.numeric(result$effect_size))
  expect_true(!is.na(result$effect_size))
})

test_that("compute_effect_sizes errors with != 2 conditions", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "effect_sizes.R"))

  expr_data <- data.table::data.table(
    m1        = rnorm(30),
    condition = rep(c("A", "B", "C"), each = 10)
  )

  expect_error(
    compute_effect_sizes(expr_data, "condition", "m1"),
    "Exactly 2 conditions required"
  )
})

test_that("compute_effect_sizes handles zero-variance gracefully", {
  skip_if_not_installed("data.table")

  source(file.path(analysis_dir, "effect_sizes.R"))

  expr_data <- data.table::data.table(
    m1        = c(rep(5, 10), rep(5, 10)),
    condition = rep(c("A", "B"), each = 10)
  )

  result <- compute_effect_sizes(expr_data, "condition", "m1",
                                  methods = "cohens_d")

  expect_equal(result$effect_size, 0)
  expect_equal(result$interpretation, "zero_variance")
})

# ======================== run_correlation_analysis() ========================
test_that("run_correlation_analysis returns cor_matrix and significant_pairs", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(analysis_dir, "correlation_analysis.R"))

  set.seed(43)
  n <- 100
  x <- rnorm(n)
  mat <- matrix(c(x, x + rnorm(n, sd = 0.1), rnorm(n)), ncol = 3)
  colnames(mat) <- c("A", "B", "C")

  result <- run_correlation_analysis(mat, method = "spearman", output_dir = NULL)

  expect_type(result, "list")
  expect_true("cor_matrix" %in% names(result))
  expect_true("p_matrix" %in% names(result))
  expect_true("significant_pairs" %in% names(result))

  # cor_matrix dimensions
  expect_equal(dim(result$cor_matrix), c(3, 3))
  expect_equal(rownames(result$cor_matrix), c("A", "B", "C"))

  # A and B should be highly correlated
  expect_true(result$cor_matrix["A", "B"] > 0.9)

  # Diagonal should be 1 (unname strips marker names from the diag vector)
  expect_equal(unname(diag(result$cor_matrix)), c(1, 1, 1), tolerance = 1e-10)

  # Significant pairs should include A-B
  sig <- result$significant_pairs
  if (nrow(sig) > 0) {
    ab_pair <- sig[(sig$marker_1 == "A" & sig$marker_2 == "B") |
                    (sig$marker_1 == "B" & sig$marker_2 == "A"), ]
    expect_true(nrow(ab_pair) >= 1)
  }
})

test_that("run_correlation_analysis works with Pearson method", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(analysis_dir, "correlation_analysis.R"))

  set.seed(44)
  mat <- matrix(rnorm(200), ncol = 4)
  colnames(mat) <- paste0("V", 1:4)

  result <- run_correlation_analysis(mat, method = "pearson", output_dir = NULL)

  expect_equal(dim(result$cor_matrix), c(4, 4))
  expect_true(all(!is.na(diag(result$cor_matrix))))
})

test_that("run_correlation_analysis handles unnamed matrix", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(analysis_dir, "correlation_analysis.R"))

  mat <- matrix(rnorm(60), ncol = 3)
  # No colnames set

  result <- run_correlation_analysis(mat, method = "spearman", output_dir = NULL)

  # Should auto-name columns V1, V2, V3
  expect_equal(colnames(result$cor_matrix), c("V1", "V2", "V3"))
})

# ======================== build_abundance_matrix() ==========================
test_that("build_abundance_matrix filters zero-cell samples", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "gui_utils.R"))
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "memory_utils.R"))
  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(analysis_dir, "differential_abundance.R"))

  tmp_dir <- tempfile("abund_input_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Sample with cells
  dt1 <- data.table::data.table(
    marker = rnorm(30),
    `Cluster_1 (Clustering annotated)` = c(rep(1L, 10), rep(0L, 20)),
    `Cluster_2 (Clustering annotated)` = c(rep(0L, 10), rep(1L, 15), rep(0L, 5))
  )
  data.table::fwrite(dt1, file.path(tmp_dir, "sampleA.csv"), sep = ";")

  # Another sample with cells
  dt2 <- data.table::data.table(
    marker = rnorm(20),
    `Cluster_1 (Clustering annotated)` = c(rep(1L, 5), rep(0L, 15)),
    `Cluster_2 (Clustering annotated)` = c(rep(0L, 5), rep(1L, 8), rep(0L, 7))
  )
  data.table::fwrite(dt2, file.path(tmp_dir, "sampleB.csv"), sep = ";")

  result <- build_abundance_matrix(
    data_dir = tmp_dir,
    cluster_regex = "Clustering annotated",
    min_cluster_size = 5
  )

  expect_type(result, "list")
  expect_true("abundance_df" %in% names(result))
  expect_true("all_clusters" %in% names(result))
  expect_true("valid_clusters" %in% names(result))

  df <- result$abundance_df
  expect_true("file_name" %in% names(df))
  expect_true("total_cells" %in% names(df))
  expect_true("log_total" %in% names(df))

  # Both samples should be present (both have cells)
  expect_equal(nrow(df), 2)
  # All total_cells > 0
  expect_true(all(df$total_cells > 0))

  # Check cluster counts
  expect_equal(df[df$file_name == "sampleA", "Cluster_1"][[1]], 10)
  expect_equal(df[df$file_name == "sampleA", "Cluster_2"][[1]], 15)
})

test_that("build_abundance_matrix respects min_cluster_size", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "gui_utils.R"))
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "memory_utils.R"))
  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(shared_dir, "plot_utils.R"))
  source(file.path(analysis_dir, "differential_abundance.R"))

  tmp_dir <- tempfile("abund_min_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  dt <- data.table::data.table(
    m = rnorm(100),
    `Big (Clustering annotated)`   = c(rep(1L, 80), rep(0L, 20)),
    `Small (Clustering annotated)` = c(rep(0L, 97), rep(1L, 3))
  )
  data.table::fwrite(dt, file.path(tmp_dir, "s1.csv"), sep = ";")

  result <- build_abundance_matrix(tmp_dir, min_cluster_size = 10)

  # "Big" should pass, "Small" (only 3 cells total) should be filtered
  expect_true("Big" %in% result$valid_clusters)
  expect_false("Small" %in% result$valid_clusters)
  # But "Small" should still be in all_clusters
  expect_true("Small" %in% result$all_clusters)
})

# ======================== combine_expression_matrices() =====================
test_that("combine_expression_matrices merges and subsamples correctly", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(analysis_dir, "expression_matrices.R"))

  dt1 <- data.table::data.table(A = 1:10, B = 11:20, C = 21:30)
  dt2 <- data.table::data.table(A = 31:40, B = 41:50, C = 51:60)

  expr_list <- list(sample1 = dt1, sample2 = dt2)

  # Without subsampling
  mat <- combine_expression_matrices(expr_list, markers = c("A", "B"), max_cells = NULL)

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 20)                          # 10 + 10
  expect_equal(ncol(mat), 2)                           # A and B only
  expect_equal(colnames(mat), c("A", "B"))

  # With subsampling
  mat_sub <- combine_expression_matrices(expr_list, markers = c("A", "B"),
                                          max_cells = 5, seed = 1)
  expect_equal(nrow(mat_sub), 5)
  expect_equal(ncol(mat_sub), 2)
})

test_that("combine_expression_matrices returns NULL when no markers match", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(analysis_dir, "expression_matrices.R"))

  dt1 <- data.table::data.table(X = 1:5, Y = 6:10)
  expr_list <- list(s1 = dt1)

  result <- combine_expression_matrices(expr_list, markers = c("MISSING1", "MISSING2"))
  expect_null(result)
})

test_that("combine_expression_matrices skips samples missing markers", {
  skip_if_not_installed("data.table")

  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(analysis_dir, "expression_matrices.R"))

  dt1 <- data.table::data.table(A = 1:5, B = 6:10)
  dt2 <- data.table::data.table(A = 11:15, C = 16:20)  # no B column

  expr_list <- list(s1 = dt1, s2 = dt2)

  # Request A and B; s2 is missing B and should be skipped with a warning.
  expect_warning(
    mat <- combine_expression_matrices(expr_list, markers = c("A", "B")),
    "Sample 's2' is missing markers"
  )

  expect_true(is.matrix(mat))
  # Only s1 contributes, so 5 rows x 2 cols
  expect_equal(nrow(mat), 5)
  expect_equal(ncol(mat), 2)
  expect_equal(colnames(mat), c("A", "B"))
})

# ======================== run_pca() from stats_utils ========================
test_that("run_pca returns expected structure", {
  source(file.path(shared_dir, "stats_utils.R"))

  set.seed(45)
  mat <- matrix(rnorm(200), nrow = 40, ncol = 5)
  colnames(mat) <- paste0("m", 1:5)

  result <- run_pca(mat, n_components = 2, scale = TRUE)

  expect_type(result, "list")
  expect_true("coordinates" %in% names(result))
  expect_true("variance_explained" %in% names(result))
  expect_true("loadings" %in% names(result))

  expect_equal(dim(result$coordinates), c(40, 2))
  expect_length(result$variance_explained, 2)
  expect_true(all(result$variance_explained > 0))
  expect_true(sum(result$variance_explained) <= 1)
  expect_equal(dim(result$loadings), c(5, 2))
})

# ======================== get_top_features() ================================
test_that("get_top_features returns correct top markers by variance", {
  source(file.path(shared_dir, "stats_utils.R"))

  set.seed(46)
  mat <- matrix(0, nrow = 50, ncol = 4)
  colnames(mat) <- c("low_var", "med_var", "high_var", "zero_var")
  mat[, "low_var"]  <- rnorm(50, sd = 0.1)
  mat[, "med_var"]  <- rnorm(50, sd = 1.0)
  mat[, "high_var"] <- rnorm(50, sd = 10.0)
  mat[, "zero_var"] <- rep(5, 50)

  top2 <- get_top_features(mat, n_top = 2, method = "variance")
  expect_length(top2, 2)
  expect_equal(top2[1], "high_var")
  expect_equal(top2[2], "med_var")
})

test_that("get_top_features by mean returns correct order", {
  source(file.path(shared_dir, "stats_utils.R"))

  mat <- matrix(c(rep(1, 10), rep(5, 10), rep(10, 10)), ncol = 3)
  colnames(mat) <- c("low", "mid", "high")

  top <- get_top_features(mat, n_top = 2, method = "mean")
  expect_equal(top[1], "high")
  expect_equal(top[2], "mid")
})

# ======================== classify_csv_columns() ============================
test_that("classify_csv_columns separates markers from gating columns", {
  source(file.path(shared_dir, "io_utils.R"))

  cols <- c("Index", "SSC.A | SSC.A", "FSC.A | FSC.A",
            "CD3 | BV421.A", "CD4 | BUV496.A",
            "Cluster_1 (Clustering annotated)",
            "Gate_1 (Manual gating)")

  classified <- classify_csv_columns(cols)

  expect_true("Index" %in% classified$index)
  expect_true("SSC.A | SSC.A" %in% classified$scatter)
  expect_true("FSC.A | FSC.A" %in% classified$scatter)
  expect_true("CD3 | BV421.A" %in% classified$markers)
  expect_true("CD4 | BUV496.A" %in% classified$markers)
  expect_true("Cluster_1 (Clustering annotated)" %in% classified$gating)
  expect_true("Gate_1 (Manual gating)" %in% classified$gating)
})

# ======================== clean_cluster_name() ==============================
test_that("clean_cluster_name removes expected prefixes and suffixes", {
  source(file.path(shared_dir, "cluster_utils.R"))

  expect_equal(clean_cluster_name("Cluster_1 (Clustering annotated)"), "Cluster_1")
  expect_equal(clean_cluster_name("AI_Cluster_CD4TCs_1"), "CD4TCs_1")
  expect_equal(clean_cluster_name(""), "Unknown")
  expect_equal(clean_cluster_name(NA), "Unknown")
})

# ======================== Internal effect size helpers =======================
test_that(".interpret_cohens_d thresholds are correct", {
  source(file.path(analysis_dir, "effect_sizes.R"))

  expect_equal(.interpret_cohens_d(0.1), "negligible")
  expect_equal(.interpret_cohens_d(0.3), "small")
  expect_equal(.interpret_cohens_d(0.6), "medium")
  expect_equal(.interpret_cohens_d(1.0), "large")
  expect_equal(.interpret_cohens_d(NA), "NA")
})

test_that(".interpret_cliffs_delta thresholds are correct", {
  source(file.path(analysis_dir, "effect_sizes.R"))

  expect_equal(.interpret_cliffs_delta(0.1), "negligible")
  expect_equal(.interpret_cliffs_delta(0.2), "small")
  expect_equal(.interpret_cliffs_delta(0.4), "medium")
  expect_equal(.interpret_cliffs_delta(0.5), "large")
  expect_equal(.interpret_cliffs_delta(NA), "NA")
})
