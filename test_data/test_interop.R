#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_interop.R — Method Interoperability Integration Tests
# =============================================================================
# Tests that all clustering, dimensionality reduction, and batch correction
# methods are interchangeable via the dispatcher layer. Covers ~29 test cases.
#
# Run from the project root:
#   Rscript test_data/test_interop.R
# =============================================================================

cat("=== Method Interoperability Integration Tests ===\n\n")

# --- Setup ---
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_interop.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) {
  base_dir <- getwd()
}

test_data_dir <- file.path(base_dir, "test_data")
output_dir    <- file.path(base_dir, "test_output_interop")
shared_dir    <- file.path(base_dir, "src", "R_shared")

# Point reticulate to the conda env Python (for leiden/leidenalg)
conda_python <- file.path("/mnt/volume/bin/miniconda/envs/cytometry_env/bin/python")
if (file.exists(conda_python)) {
  Sys.setenv(RETICULATE_PYTHON = conda_python)
}

# Set working directory so src/ scripts can find R_shared via dirname(getwd())
old_wd <- getwd()
setwd(file.path(base_dir, "pipelines"))

# Source shared utilities
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "fcs_utils.R"))

# Source dispatchers
source(file.path(shared_dir, "clustering_dispatch.R"))
source(file.path(shared_dir, "dimred_dispatch.R"))
source(file.path(shared_dir, "batch_correction_dispatch.R"))

# Clean output
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
make_outdir(output_dir)

# Check which optional packages are available
has_destiny <- requireNamespace("destiny", quietly = TRUE)
if (!has_destiny) cat("Note: 'destiny' package not available — diffusion_map tests will be skipped\n")

# Test counter
pass <- 0
fail <- 0
skip <- 0
test_results <- list()

run_test <- function(name, expr) {
  cat(sprintf("  [TEST] %s ... ", name))
  result <- tryCatch({
    force(expr)
    cat("PASS\n")
    test_results[[name]] <<- "PASS"
    pass <<- pass + 1
    TRUE
  }, error = function(e) {
    cat(sprintf("FAIL: %s\n", e$message))
    test_results[[name]] <<- paste("FAIL:", e$message)
    fail <<- fail + 1
    FALSE
  })
  result
}

skip_test <- function(name, reason) {
  cat(sprintf("  [TEST] %s ... SKIP: %s\n", name, reason))
  test_results[[name]] <<- paste("SKIP:", reason)
  skip <<- skip + 1
}

# --- Build test data ----------------------------------------------------------
markers <- c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR")

cat("Loading test data...\n")
source(file.path(base_dir, "src", "loading", "load_csv_expression.R"))
source(file.path(base_dir, "src", "analysis", "expression_matrices.R"))
source(file.path(base_dir, "src", "clustering", "leiden.R"))

expr_list <- load_csv_expression(file.path(test_data_dir, "csv_expression"))
combined_matrix <- combine_expression_matrices(expr_list, markers, max_cells = 200)
sample_labels <- rep(names(expr_list), each = ceiling(nrow(combined_matrix) / length(expr_list)))
sample_labels <- sample_labels[1:nrow(combined_matrix)]
cat("Test matrix:", nrow(combined_matrix), "cells x", ncol(combined_matrix), "markers\n\n")


# =============================================================================
# SECTION 1: Clustering Dispatcher — Each Method Individually
# =============================================================================
cat("\n--- Section 1: Clustering Dispatcher (individual methods) ---\n")

run_test("clustering_dispatch: leiden", {
  ids <- run_clustering(combined_matrix, method = "leiden",
                        params = list(n_neighbors = 5, resolution = 0.5, seed = 42))
  stopifnot(is.integer(ids))
  stopifnot(length(ids) == nrow(combined_matrix))
  stopifnot(max(ids) >= 1)
})

run_test("clustering_dispatch: flowsom", {
  ids <- run_clustering(combined_matrix, method = "flowsom",
                        params = list(n_clusters = 5, grid_xdim = 5, grid_ydim = 5, seed = 42))
  stopifnot(is.integer(ids))
  stopifnot(length(ids) == nrow(combined_matrix))
  stopifnot(max(ids) >= 1)
})

run_test("clustering_dispatch: phenograph", {
  ids <- run_clustering(combined_matrix, method = "phenograph",
                        params = list(k = 5, seed = 42))
  stopifnot(is.integer(ids))
  stopifnot(length(ids) == nrow(combined_matrix))
  stopifnot(max(ids) >= 1)
})

run_test("clustering_dispatch: dbscan", {
  ids <- run_clustering(combined_matrix, method = "dbscan",
                        params = list(min_pts = 3, seed = 42))
  stopifnot(is.integer(ids))
  stopifnot(length(ids) == nrow(combined_matrix))
  # DBSCAN can have cluster 0 (noise) — at least some assignments exist
  meta <- attr(ids, "dbscan_metadata")
  stopifnot(!is.null(meta))
  stopifnot("noise_fraction" %in% names(meta))
})


# =============================================================================
# SECTION 2: DimRed Dispatcher — Each Method Individually
# =============================================================================
cat("\n--- Section 2: DimRed Dispatcher (individual methods) ---\n")

run_test("dimred_dispatch: umap", {
  result <- run_dimred(combined_matrix, method = "umap",
                       sample_labels = sample_labels,
                       output_dir = file.path(output_dir, "dimred_umap"),
                       params = list(n_neighbors = 5, seed = 42))
  stopifnot(is.matrix(result$coords) || is.data.frame(result$coords))
  stopifnot(nrow(result$coords) == nrow(combined_matrix))
  stopifnot(result$method == "umap")
})

run_test("dimred_dispatch: tsne", {
  result <- run_dimred(combined_matrix, method = "tsne",
                       sample_labels = sample_labels,
                       output_dir = file.path(output_dir, "dimred_tsne"),
                       params = list(tsne_perplexity = 5, seed = 42))
  stopifnot(is.matrix(result$coords) || is.data.frame(result$coords))
  stopifnot(nrow(result$coords) == nrow(combined_matrix))
  stopifnot(result$method == "tsne")
})

if (has_destiny) {
  run_test("dimred_dispatch: diffusion_map", {
    result <- run_dimred(combined_matrix, method = "diffusion_map",
                         sample_labels = sample_labels,
                         output_dir = file.path(output_dir, "dimred_dm"),
                         params = list(dm_n_components = 3, seed = 42))
    stopifnot(is.matrix(result$coords) || is.data.frame(result$coords))
    stopifnot(nrow(result$coords) == nrow(combined_matrix))
    stopifnot(result$method == "diffusion_map")
    stopifnot(!is.null(result$extra$eigenvalues))
  })
} else {
  skip_test("dimred_dispatch: diffusion_map", "destiny package not installed")
}


# =============================================================================
# SECTION 3: Batch Correction Dispatcher
# =============================================================================
cat("\n--- Section 3: Batch Correction Dispatcher ---\n")

# Create batch labels for the matrix (2 batches)
batch_labels_cell <- rep(c("Batch_A", "Batch_B"), length.out = nrow(combined_matrix))

run_test("batch_dispatch: none (pass-through)", {
  result <- run_batch_correction_dispatch(
    data = combined_matrix, method = "none",
    batch_labels = batch_labels_cell,
    markers_of_interest = markers,
    sample_names = sample_labels
  )
  stopifnot(is.matrix(result$corrected_matrix))
  stopifnot(nrow(result$corrected_matrix) == nrow(combined_matrix))
  stopifnot(all(result$corrected_matrix == combined_matrix))
})

run_test("batch_dispatch: combat", {
  result <- run_batch_correction_dispatch(
    data = combined_matrix, method = "combat",
    batch_labels = batch_labels_cell,
    markers_of_interest = markers,
    sample_names = sample_labels,
    params = list(combat_method = "parametric")
  )
  stopifnot(is.matrix(result$corrected_matrix))
  stopifnot(nrow(result$corrected_matrix) == nrow(combined_matrix))
  stopifnot(ncol(result$corrected_matrix) == ncol(combined_matrix))
  stopifnot(!is.null(result$markers))
})

run_test("batch_dispatch: harmony (matrix input)", {
  # Harmony needs sample-level batch labels: assign half the samples per batch
  unique_samps <- unique(sample_labels)
  batch_map <- setNames(
    rep(c("Batch_A", "Batch_B"), length.out = length(unique_samps)),
    unique_samps
  )
  harmony_batch_cell <- batch_map[sample_labels]
  result <- run_batch_correction_dispatch(
    data = combined_matrix, method = "harmony",
    batch_labels = harmony_batch_cell,
    markers_of_interest = markers,
    sample_names = sample_labels,
    params = list(batch_do_pca = FALSE)
  )
  stopifnot(is.matrix(result$corrected_matrix))
  stopifnot(nrow(result$corrected_matrix) == nrow(combined_matrix))
  stopifnot(!is.null(result$sample_labels))
})


# =============================================================================
# SECTION 4: Clustering x DimRed Chains (4 x 3 = 12 tests)
# =============================================================================
cat("\n--- Section 4: Clustering x DimRed Chains ---\n")

cluster_methods <- c("leiden", "flowsom", "phenograph", "dbscan")
dimred_methods  <- c("umap", "tsne", "diffusion_map")

cluster_params <- list(
  leiden     = list(n_neighbors = 5, resolution = 0.5, seed = 42),
  flowsom    = list(n_clusters = 5, grid_xdim = 5, grid_ydim = 5, seed = 42),
  phenograph = list(k = 5, seed = 42),
  dbscan     = list(min_pts = 3, seed = 42)
)

dimred_params <- list(
  umap          = list(n_neighbors = 5, seed = 42),
  tsne          = list(tsne_perplexity = 5, seed = 42),
  diffusion_map = list(dm_n_components = 3, seed = 42)
)

for (cm in cluster_methods) {
  for (dm in dimred_methods) {
    test_name <- paste0("chain: ", cm, " + ", dm)
    out_sub <- file.path(output_dir, paste0("chain_", cm, "_", dm))

    # Skip diffusion_map tests if destiny not available
    if (dm == "diffusion_map" && !has_destiny) {
      skip_test(test_name, "destiny package not installed")
      next
    }

    run_test(test_name, {
      # Step 1: Cluster
      ids <- run_clustering(combined_matrix, method = cm, params = cluster_params[[cm]])
      stopifnot(length(ids) == nrow(combined_matrix))

      # Step 2: DimRed
      dr <- run_dimred(combined_matrix, method = dm,
                       sample_labels = sample_labels,
                       cluster_labels = ids,
                       output_dir = out_sub,
                       params = dimred_params[[dm]])
      stopifnot(nrow(dr$coords) == nrow(combined_matrix))

      # Step 3: Write gating CSV (verify bridge works with all cluster outputs)
      clusters_to_gating_csv(
        combined_matrix = combined_matrix,
        sample_labels   = sample_labels,
        cluster_ids     = ids,
        output_dir      = out_sub
      )
      csvs <- list.files(out_sub, pattern = "\\.csv$")
      stopifnot(length(csvs) > 0)
    })
  }
}


# =============================================================================
# SECTION 5: Batch Correction + Clustering Chains
# =============================================================================
cat("\n--- Section 5: Batch + Clustering Chains ---\n")

batch_methods_to_test <- c("combat", "none")
cluster_methods_to_test <- c("leiden", "flowsom")

for (bm in batch_methods_to_test) {
  for (cm in cluster_methods_to_test) {
    test_name <- paste0("batch_chain: ", bm, " + ", cm)

    run_test(test_name, {
      # Step 1: Batch correction
      bc <- run_batch_correction_dispatch(
        data = combined_matrix, method = bm,
        batch_labels = batch_labels_cell,
        markers_of_interest = markers,
        sample_names = sample_labels,
        params = list(combat_method = "parametric")
      )
      stopifnot(is.matrix(bc$corrected_matrix))

      # Step 2: Cluster on corrected data
      ids <- run_clustering(bc$corrected_matrix, method = cm,
                            params = cluster_params[[cm]])
      stopifnot(length(ids) == nrow(bc$corrected_matrix))
      stopifnot(max(ids) >= 1)
    })
  }
}


# =============================================================================
# SECTION 6: Interchangeability — Same Downstream Analysis
# =============================================================================
cat("\n--- Section 6: Interchangeability ---\n")

source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))

run_test("interchangeability: leiden vs flowsom → split + abundance", {
  # Cluster with Leiden
  ids_leiden <- run_clustering(combined_matrix, method = "leiden",
                               params = list(n_neighbors = 5, resolution = 0.5, seed = 42))
  gating_dir_l <- file.path(output_dir, "interop_leiden_gating")
  clusters_to_gating_csv(combined_matrix, sample_labels, ids_leiden, gating_dir_l)

  split_l <- split_by_cluster(
    input_dir = gating_dir_l,
    output_dir = file.path(output_dir, "interop_leiden_split"),
    cluster_regex = "Clustering annotated",
    markers_of_interest = markers
  )
  stopifnot(!is.null(split_l$tracking_matrix))

  abund_l <- build_abundance_matrix(
    data_dir = gating_dir_l,
    cluster_regex = "Clustering annotated",
    min_cluster_size = 1
  )
  stopifnot(nrow(abund_l$abundance_df) > 0)

  # Cluster with FlowSOM
  ids_fsom <- run_clustering(combined_matrix, method = "flowsom",
                              params = list(n_clusters = 5, grid_xdim = 5, grid_ydim = 5, seed = 42))
  gating_dir_f <- file.path(output_dir, "interop_flowsom_gating")
  clusters_to_gating_csv(combined_matrix, sample_labels, ids_fsom, gating_dir_f)

  split_f <- split_by_cluster(
    input_dir = gating_dir_f,
    output_dir = file.path(output_dir, "interop_flowsom_split"),
    cluster_regex = "Clustering annotated",
    markers_of_interest = markers
  )
  stopifnot(!is.null(split_f$tracking_matrix))

  abund_f <- build_abundance_matrix(
    data_dir = gating_dir_f,
    cluster_regex = "Clustering annotated",
    min_cluster_size = 1
  )
  stopifnot(nrow(abund_f$abundance_df) > 0)

  # Both should produce valid tracking matrices with the same structure
  stopifnot(identical(names(split_l$tracking_matrix)[1], names(split_f$tracking_matrix)[1]))
})

run_test("interchangeability: phenograph vs dbscan → gating CSV → abundance", {
  ids_pg <- run_clustering(combined_matrix, method = "phenograph",
                            params = list(k = 5, seed = 42))
  gating_dir_pg <- file.path(output_dir, "interop_pg_gating")
  clusters_to_gating_csv(combined_matrix, sample_labels, ids_pg, gating_dir_pg)
  abund_pg <- build_abundance_matrix(
    data_dir = gating_dir_pg, cluster_regex = "Clustering annotated", min_cluster_size = 1
  )
  stopifnot(nrow(abund_pg$abundance_df) > 0)

  ids_db <- run_clustering(combined_matrix, method = "dbscan",
                            params = list(min_pts = 3, seed = 42))
  gating_dir_db <- file.path(output_dir, "interop_db_gating")
  clusters_to_gating_csv(combined_matrix, sample_labels, ids_db, gating_dir_db)
  abund_db <- build_abundance_matrix(
    data_dir = gating_dir_db, cluster_regex = "Clustering annotated", min_cluster_size = 1
  )
  stopifnot(nrow(abund_db$abundance_df) > 0)
})


# =============================================================================
# SECTION 7: Config-Driven Pipeline Routing
# =============================================================================
cat("\n--- Section 7: Config-Driven Pipeline Routing ---\n")

run_test("config-driven: flowsom + tsne routing via params", {
  # Simulate config
  test_params <- list(
    clustering_method = "flowsom",
    dimred_method     = "tsne",
    n_clusters        = 5,
    grid_xdim         = 5,
    grid_ydim         = 5,
    tsne_perplexity   = 5,
    seed              = 42
  )

  # Use dispatcher with config-style params
  ids <- run_clustering(
    combined_matrix,
    method = test_params$clustering_method,
    params = test_params
  )
  stopifnot(is.integer(ids))
  stopifnot(length(ids) == nrow(combined_matrix))

  dr <- run_dimred(
    combined_matrix,
    method = test_params$dimred_method,
    sample_labels = sample_labels,
    cluster_labels = ids,
    output_dir = file.path(output_dir, "config_driven"),
    params = test_params
  )
  stopifnot(dr$method == "tsne")
  stopifnot(nrow(dr$coords) == nrow(combined_matrix))
})


# =============================================================================
# Summary
# =============================================================================
setwd(old_wd)

cat("\n=== Interop Test Summary ===\n")
cat(sprintf("  PASSED:  %d\n", pass))
cat(sprintf("  FAILED:  %d\n", fail))
cat(sprintf("  SKIPPED: %d\n", skip))
cat(sprintf("  TOTAL:   %d\n", pass + fail + skip))

if (fail > 0) {
  cat("\nFailed tests:\n")
  for (nm in names(test_results)) {
    if (grepl("^FAIL", test_results[[nm]])) {
      cat(sprintf("  - %s: %s\n", nm, test_results[[nm]]))
    }
  }
}

cat(sprintf("\nTest output: %s\n", output_dir))
cat("=== Done ===\n")

# Exit with error code if failures
if (fail > 0) quit(status = 1)
