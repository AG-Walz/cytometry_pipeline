library(testthat)

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

# ---------------------------------------------------------------------------
# Helper: resolve paths relative to the project root so tests work both
# when run via testthat::test_dir() and from the repo root.
# ---------------------------------------------------------------------------
find_shared_dir <- function() {
  cyto_shared_dir()
}

shared_dir  <- find_shared_dir()
cluster_dir <- file.path(dirname(shared_dir), "clustering")

# Source shared utilities (needed by every clustering module)
source(file.path(shared_dir, "io_utils.R"))

# ---------------------------------------------------------------------------
# Synthetic data factory — reused across many tests
# ---------------------------------------------------------------------------
make_test_matrix <- function(n_cells = 100, n_markers = 5, n_blobs = 3, seed = 42) {
  set.seed(seed)
  centers <- matrix(rnorm(n_blobs * n_markers, sd = 5), nrow = n_blobs)
  labels  <- sample(seq_len(n_blobs), n_cells, replace = TRUE)
  mat <- centers[labels, ] + matrix(rnorm(n_cells * n_markers, sd = 0.5),
                                     nrow = n_cells)
  colnames(mat) <- paste0("marker_", seq_len(n_markers))
  list(matrix = mat, labels = labels)
}

# ============================= run_leiden() =================================
test_that("run_leiden returns integer vector of correct length", {
  skip_if_not_installed("leiden")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("uwot")

  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(cluster_dir, "leiden.R"))

  syn <- make_test_matrix(n_cells = 80, n_markers = 5, seed = 1)
  # Use small n_neighbors since dataset is tiny
  clusters <- run_leiden(syn$matrix, n_neighbors = 10, resolution = 0.5, seed = 1)

  expect_type(clusters, "double")                     # leiden returns numeric
  ids <- as.integer(clusters)
  expect_length(ids, nrow(syn$matrix))
  expect_true(all(ids >= 1))
  expect_true(length(unique(ids)) >= 1)
})

# ============================= run_flowsom() ================================
test_that("run_flowsom returns integer vector of correct length", {
  skip_if_not_installed("FlowSOM")
  skip_if_not_installed("flowCore")
  skip_if_not_installed("Biobase")

  source(file.path(shared_dir, "flowframe_utils.R"))
  source(file.path(cluster_dir, "flowsom.R"))

  syn <- make_test_matrix(n_cells = 100, n_markers = 5, seed = 2)
  clusters <- run_flowsom(syn$matrix, n_clusters = 3,
                           grid_xdim = 5, grid_ydim = 5, seed = 2)

  expect_type(clusters, "integer")
  expect_length(clusters, nrow(syn$matrix))
  expect_true(all(clusters >= 1))
  expect_true(length(unique(clusters)) >= 1)
  expect_true(length(unique(clusters)) <= 25)         # grid_xdim * grid_ydim
})

# ============================= run_phenograph() =============================
test_that("run_phenograph returns integer vector of correct length", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("uwot")

  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(cluster_dir, "phenograph.R"))

  syn <- make_test_matrix(n_cells = 60, n_markers = 4, seed = 3)
  clusters <- run_phenograph(syn$matrix, k = 10, seed = 3)

  expect_type(clusters, "integer")
  expect_length(clusters, nrow(syn$matrix))
  expect_true(all(clusters >= 1))
  expect_true(length(unique(clusters)) >= 1)
})

test_that("run_phenograph auto-reduces k when k >= n_cells", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("uwot")

  source(file.path(shared_dir, "stats_utils.R"))
  source(file.path(cluster_dir, "phenograph.R"))

  syn <- make_test_matrix(n_cells = 20, n_markers = 3, seed = 4)
  # k=30 > 20 cells => should auto-reduce
  clusters <- run_phenograph(syn$matrix, k = 30, seed = 4)
  expect_length(clusters, 20)
})

# ============================= run_dbscan() =================================
test_that("run_dbscan returns list with expected keys", {
  skip_if_not_installed("dbscan")

  source(file.path(cluster_dir, "dbscan.R"))

  syn <- make_test_matrix(n_cells = 100, n_markers = 5, n_blobs = 3, seed = 5)
  result <- run_dbscan(syn$matrix, eps = 2.0, min_pts = 3, auto_eps = FALSE, seed = 5)

  expect_type(result, "list")
  expect_named(result,
               c("cluster_ids", "noise_indices", "n_clusters",
                 "eps_used", "noise_fraction"),
               ignore.order = TRUE)

  expect_type(result$cluster_ids, "integer")
  expect_length(result$cluster_ids, nrow(syn$matrix))
  expect_true(all(result$cluster_ids >= 0))            # 0 = noise
  expect_equal(result$eps_used, 2.0)
  expect_true(is.numeric(result$noise_fraction))
  expect_true(result$noise_fraction >= 0 && result$noise_fraction <= 1)
  expect_equal(result$n_clusters, max(result$cluster_ids))
})

test_that("run_dbscan auto_eps estimates epsilon when eps is NULL", {
  skip_if_not_installed("dbscan")

  source(file.path(cluster_dir, "dbscan.R"))

  syn <- make_test_matrix(n_cells = 80, n_markers = 4, seed = 6)
  result <- run_dbscan(syn$matrix, eps = NULL, auto_eps = TRUE, min_pts = 5, seed = 6)

  expect_type(result, "list")
  expect_true(!is.null(result$eps_used))
  expect_true(is.numeric(result$eps_used))
  expect_true(result$eps_used > 0)
})

test_that("run_dbscan NaN eps triggers fallback to 0.5", {
  skip_if_not_installed("dbscan")

  source(file.path(cluster_dir, "dbscan.R"))

  # Override estimate_eps to return NaN to test the fallback
  local_estimate_eps <- estimate_eps
  assign("estimate_eps", function(...) NaN, envir = globalenv())
  on.exit(assign("estimate_eps", local_estimate_eps, envir = globalenv()), add = TRUE)

  syn <- make_test_matrix(n_cells = 50, n_markers = 3, seed = 7)
  result <- run_dbscan(syn$matrix, eps = NULL, auto_eps = TRUE, min_pts = 3, seed = 7)

  expect_equal(result$eps_used, 0.5)
})

test_that("run_dbscan errors when eps=NULL and auto_eps=FALSE", {
  skip_if_not_installed("dbscan")

  source(file.path(cluster_dir, "dbscan.R"))

  syn <- make_test_matrix(n_cells = 30, n_markers = 3, seed = 8)
  expect_error(
    run_dbscan(syn$matrix, eps = NULL, auto_eps = FALSE),
    "eps must be provided"
  )
})

# ======================== estimate_eps() ====================================
test_that("estimate_eps returns a finite positive number", {
  skip_if_not_installed("dbscan")

  source(file.path(cluster_dir, "dbscan.R"))

  syn <- make_test_matrix(n_cells = 200, n_markers = 4, seed = 9)
  eps_val <- estimate_eps(syn$matrix, k = 5, max_sample = 100, seed = 9)

  expect_true(is.numeric(eps_val))
  expect_length(eps_val, 1)
  expect_true(is.finite(eps_val))
  expect_true(eps_val > 0)
})

# ====================== clusters_to_gating_csv() ============================
test_that("clusters_to_gating_csv writes correct per-sample CSV files", {
  source(file.path(cluster_dir, "leiden.R"))

  n_cells <- 30
  n_markers <- 3
  set.seed(10)
  mat <- matrix(rnorm(n_cells * n_markers), nrow = n_cells)
  colnames(mat) <- paste0("M", seq_len(n_markers))

  sample_labels <- rep(c("sampleA", "sampleB", "sampleC"), each = 10)
  cluster_ids   <- rep(c(1L, 2L, 3L), each = 10)

  tmp_dir <- tempfile("gating_csv_test_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  out_files <- clusters_to_gating_csv(mat, sample_labels, cluster_ids,
                                       output_dir = tmp_dir, sep = ";")

  # Three samples -> three files
  expect_length(out_files, 3)
  expect_true(all(file.exists(out_files)))

  # Read one and verify structure
  dt_a <- data.table::fread(out_files[1], sep = ";")
  expect_true("M1" %in% names(dt_a))
  expect_true("M2" %in% names(dt_a))
  expect_true("M3" %in% names(dt_a))

  # Binary cluster indicator columns present
  cluster_cols <- grep("Cluster_.*\\(Clustering annotated\\)", names(dt_a), value = TRUE)
  expect_true(length(cluster_cols) >= 1)

  # Each row has exactly one cluster column set to 1
  row_sums <- rowSums(dt_a[, ..cluster_cols])
  expect_true(all(row_sums == 1))

  # sample and cluster columns should NOT be in the output

  expect_false("sample" %in% names(dt_a))
  expect_false("cluster" %in% names(dt_a))

  # 10 rows per sample
  expect_equal(nrow(dt_a), 10)
})

# ======================== validate_clusters() ===============================
test_that("validate_clusters returns data.frame with expected metric names", {
  skip_if_not_installed("cluster")

  source(file.path(cluster_dir, "cluster_validation.R"))

  syn <- make_test_matrix(n_cells = 100, n_markers = 5, n_blobs = 3, seed = 11)
  cluster_ids <- syn$labels

  result <- validate_clusters(syn$matrix, cluster_ids,
                               metrics = c("silhouette", "calinski_harabasz"),
                               sample_size = 100, seed = 11)

  expect_s3_class(result, "data.frame")
  expect_true("metric" %in% names(result))
  expect_true("value" %in% names(result))
  expect_true("description" %in% names(result))

  metric_names <- result$metric
  expect_true("n_clusters" %in% metric_names)
  expect_true("n_cells" %in% metric_names)
  expect_true("silhouette" %in% metric_names)
  expect_true("calinski_harabasz" %in% metric_names)
})

test_that("validate_clusters handles single-cluster input gracefully", {
  source(file.path(cluster_dir, "cluster_validation.R"))

  mat <- matrix(rnorm(50), nrow = 10)
  colnames(mat) <- paste0("m", 1:5)
  single_cluster <- rep(1L, 10)

  expect_warning(
    result <- validate_clusters(mat, single_cluster, sample_size = 10),
    "Fewer than 2 clusters"
  )
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 1)
})

test_that("validate_clusters computes ARI when true_labels provided", {
  skip_if_not_installed("cluster")

  source(file.path(cluster_dir, "cluster_validation.R"))

  syn <- make_test_matrix(n_cells = 60, n_markers = 4, n_blobs = 3, seed = 12)
  result <- validate_clusters(syn$matrix, syn$labels,
                               true_labels = syn$labels,
                               metrics = c("silhouette"),
                               sample_size = 60, seed = 12)

  expect_true("adjusted_rand_index" %in% result$metric)
  ari_row <- result[result$metric == "adjusted_rand_index", ]
  # Perfect agreement since true_labels == cluster_ids
  expect_equal(ari_row$value, 1.0)
})

test_that("validate_clusters computes Davies-Bouldin index", {
  source(file.path(cluster_dir, "cluster_validation.R"))

  syn <- make_test_matrix(n_cells = 80, n_markers = 4, n_blobs = 3, seed = 13)
  result <- validate_clusters(syn$matrix, syn$labels,
                               metrics = c("davies_bouldin"),
                               sample_size = 80, seed = 13)

  expect_true("davies_bouldin" %in% result$metric)
  db_val <- result[result$metric == "davies_bouldin", "value"]
  expect_true(is.numeric(db_val))
  expect_true(db_val > 0)                             # DB index is always positive
})

# ======================== annotate_clusters() ===============================
test_that("annotate_clusters threshold method returns correct structure", {
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(cluster_dir, "cell_annotation.R"))

  set.seed(14)
  n <- 90
  mat <- matrix(0, nrow = n, ncol = 3)
  colnames(mat) <- c("CD3", "CD4", "CD8")

  # Cluster 1: CD3+CD4+CD8-
  mat[1:30, "CD3"] <- rnorm(30, mean = 5, sd = 0.5)
  mat[1:30, "CD4"] <- rnorm(30, mean = 4, sd = 0.5)
  mat[1:30, "CD8"] <- rnorm(30, mean = -1, sd = 0.5)

  # Cluster 2: CD3+CD4-CD8+
  mat[31:60, "CD3"] <- rnorm(30, mean = 5, sd = 0.5)
  mat[31:60, "CD4"] <- rnorm(30, mean = -1, sd = 0.5)
  mat[31:60, "CD8"] <- rnorm(30, mean = 4, sd = 0.5)

  # Cluster 3: CD3-CD4-CD8-
  mat[61:90, "CD3"] <- rnorm(30, mean = -2, sd = 0.5)
  mat[61:90, "CD4"] <- rnorm(30, mean = -2, sd = 0.5)
  mat[61:90, "CD8"] <- rnorm(30, mean = -2, sd = 0.5)

  cluster_ids <- rep(1:3, each = 30)

  result <- annotate_clusters(mat, cluster_ids,
                               markers = c("CD3", "CD4", "CD8"),
                               method = "threshold")

  expect_s3_class(result, "data.frame")
  expect_true("cluster" %in% names(result))
  expect_true("annotation" %in% names(result))
  expect_true("confidence" %in% names(result))
  expect_equal(nrow(result), 3)

  # Cluster 1 should have CD3+ and CD4+ in annotation
  cl1_annot <- result$annotation[result$cluster == 1]
  expect_true(grepl("CD3\\+", cl1_annot))
  expect_true(grepl("CD4\\+", cl1_annot))

  # Cluster 2 should have CD3+ and CD8+
  cl2_annot <- result$annotation[result$cluster == 2]
  expect_true(grepl("CD3\\+", cl2_annot))
  expect_true(grepl("CD8\\+", cl2_annot))

  # Median columns present
  expect_true("median_CD3" %in% names(result))
  expect_true("median_CD4" %in% names(result))
  expect_true("median_CD8" %in% names(result))

  # Confidence is between 0 and 1
  expect_true(all(result$confidence >= 0 & result$confidence <= 1))
})

test_that("annotate_clusters manual method returns (manual) annotations", {
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(cluster_dir, "cell_annotation.R"))

  mat <- matrix(rnorm(40), nrow = 10)
  colnames(mat) <- c("A", "B", "C", "D")
  ids <- rep(1:2, each = 5)

  result <- annotate_clusters(mat, ids, method = "manual")

  expect_equal(nrow(result), 2)
  expect_true(all(result$annotation == "(manual)"))
  expect_true(all(is.na(result$confidence)))
})

test_that("annotate_clusters correlation method needs reference matrix", {
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(cluster_dir, "cell_annotation.R"))

  mat <- matrix(rnorm(40), nrow = 10)
  colnames(mat) <- c("A", "B", "C", "D")
  ids <- rep(1:2, each = 5)

  # With reference="auto" and method="correlation", should fall back gracefully
  result <- annotate_clusters(mat, ids, method = "correlation", reference = "auto")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
})

test_that("annotate_clusters correlation method works with reference matrix", {
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(cluster_dir, "cell_annotation.R"))

  set.seed(50)
  mat <- matrix(0, nrow = 60, ncol = 3)
  colnames(mat) <- c("CD3", "CD19", "CD56")

  # Cluster 1: T cells (CD3 high)
  mat[1:30, "CD3"]  <- rnorm(30, 5, 0.3)
  mat[1:30, "CD19"] <- rnorm(30, -1, 0.3)
  mat[1:30, "CD56"] <- rnorm(30, -1, 0.3)

  # Cluster 2: B cells (CD19 high)
  mat[31:60, "CD3"]  <- rnorm(30, -1, 0.3)
  mat[31:60, "CD19"] <- rnorm(30, 5, 0.3)
  mat[31:60, "CD56"] <- rnorm(30, -1, 0.3)

  ids <- rep(1:2, each = 30)

  ref <- matrix(c(5, 0, 0,
                   0, 5, 0), nrow = 2, byrow = TRUE)
  colnames(ref) <- c("CD3", "CD19", "CD56")
  rownames(ref) <- c("T_cell", "B_cell")

  result <- annotate_clusters(mat, ids, method = "correlation", reference = ref)
  expect_equal(result$annotation[result$cluster == 1], "T_cell")
  expect_equal(result$annotation[result$cluster == 2], "B_cell")
})

# ======================== run_clustering() dispatcher =======================
test_that("run_clustering errors on invalid method", {
  source(file.path(shared_dir, "clustering_dispatch.R"))

  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  expect_error(
    run_clustering(mat, method = "kmeans_invalid"),
    "Unknown clustering method"
  )
})

test_that("run_clustering routes to dbscan correctly", {
  skip_if_not_installed("dbscan")

  source(file.path(shared_dir, "clustering_dispatch.R"))

  syn <- make_test_matrix(n_cells = 80, n_markers = 4, seed = 15)
  ids <- run_clustering(syn$matrix, method = "dbscan",
                         params = list(eps = 2.0, min_pts = 3))

  expect_type(ids, "integer")
  expect_length(ids, nrow(syn$matrix))

  # DBSCAN metadata attached as attribute
  meta <- attr(ids, "dbscan_metadata")
  expect_true(!is.null(meta))
  expect_true("noise_fraction" %in% names(meta))
  expect_true("eps_used" %in% names(meta))
})

test_that("run_clustering ensures integer output", {
  skip_if_not_installed("dbscan")

  source(file.path(shared_dir, "clustering_dispatch.R"))

  mat <- make_test_matrix(n_cells = 50, n_markers = 3, seed = 16)$matrix
  ids <- run_clustering(mat, method = "dbscan",
                         params = list(eps = 3.0, min_pts = 2))
  expect_true(is.integer(ids))
})

# ======================== cosine_similarity() ===============================
test_that("cosine_similarity returns correct values", {
  source(file.path(shared_dir, "cluster_utils.R"))
  source(file.path(cluster_dir, "cell_annotation.R"))

  # Identical vectors: similarity = 1
  expect_equal(cosine_similarity(c(1, 0, 0), c(1, 0, 0)), 1.0)

  # Orthogonal vectors: similarity = 0
  expect_equal(cosine_similarity(c(1, 0), c(0, 1)), 0.0)

  # Opposite vectors: similarity = -1
  expect_equal(cosine_similarity(c(1, 0), c(-1, 0)), -1.0)

  # Zero vector returns 0
  expect_equal(cosine_similarity(c(0, 0, 0), c(1, 2, 3)), 0.0)
})

# ======================== compute_calinski_harabasz_internal() ===============
test_that("Calinski-Harabasz index is higher for well-separated clusters", {
  source(file.path(cluster_dir, "cluster_validation.R"))

  set.seed(20)
  # Well-separated blobs
  good_mat <- rbind(
    matrix(rnorm(50, mean = 0), ncol = 5),
    matrix(rnorm(50, mean = 10), ncol = 5)
  )
  good_ids <- rep(1:2, each = 10)

  # Overlapping blobs
  bad_mat <- rbind(
    matrix(rnorm(50, mean = 0), ncol = 5),
    matrix(rnorm(50, mean = 0.1), ncol = 5)
  )
  bad_ids <- rep(1:2, each = 10)

  ch_good <- compute_calinski_harabasz_internal(good_mat, good_ids)
  ch_bad  <- compute_calinski_harabasz_internal(bad_mat, bad_ids)

  expect_true(is.numeric(ch_good))
  expect_true(is.numeric(ch_bad))
  expect_true(ch_good > ch_bad)
})

# ======================== compute_ari_manual() ==============================
test_that("compute_ari_manual returns 1 for identical labelings", {
  source(file.path(cluster_dir, "cluster_validation.R"))

  labels <- c(1, 1, 2, 2, 3, 3)
  ari <- compute_ari_manual(labels, labels)
  expect_equal(ari, 1.0)
})

test_that("compute_ari_manual returns 0 for independent labelings (approximately)", {
  source(file.path(cluster_dir, "cluster_validation.R"))

  set.seed(21)
  n <- 200
  labels1 <- sample(1:4, n, replace = TRUE)
  labels2 <- sample(1:4, n, replace = TRUE)
  ari <- compute_ari_manual(labels1, labels2)

  # Random labels should give ARI near 0
  expect_true(abs(ari) < 0.15)
})
