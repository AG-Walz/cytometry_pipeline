#!/usr/bin/env Rscript
# =============================================================================
# test_data/test_cross_pipeline.R — Integration Tests for Cross-Pipeline Paths
# =============================================================================
# Tests all viable combinations of scripts from different pipelines using
# the synthetic test data. Run from the project root:
#   Rscript test_data/test_cross_pipeline.R
# =============================================================================

cat("=== Cross-Pipeline Integration Tests ===\n\n")

# --- Setup ---
.ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.ofile)) .ofile <- "test_data/test_cross_pipeline.R"
base_dir <- normalizePath(file.path(dirname(.ofile), ".."))
if (!dir.exists(file.path(base_dir, "src"))) {
  base_dir <- getwd()
}

test_data_dir <- file.path(base_dir, "test_data")
output_dir    <- file.path(base_dir, "test_output_cross_pipeline")
shared_dir    <- file.path(base_dir, "src", "R_shared")

# All src/ scripts compute shared_dir via: dirname(getwd()) + "src/R_shared"
# This works when getwd() is one level below the project root (e.g. pipelines/).
# Set getwd to pipelines/ so dirname(getwd()) = project root.
# Point reticulate to the conda env Python (for leiden/leidenalg)
conda_python <- file.path("/mnt/volume/bin/miniconda/envs/cytometry_env/bin/python")
if (file.exists(conda_python)) {
  Sys.setenv(RETICULATE_PYTHON = conda_python)
}

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

# Clean output
if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
make_outdir(output_dir)

# Test counter
pass <- 0
fail <- 0
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

markers <- c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR")

# =============================================================================
# PATH 1: CSV Expression → combine → matrix-based UMAP
# =============================================================================
cat("\n--- Path 1: CSV Expression → Matrix UMAP ---\n")

source(file.path(base_dir, "src", "loading", "load_csv_expression.R"))
source(file.path(base_dir, "src", "analysis", "expression_matrices.R"))
source(file.path(base_dir, "src", "visualization", "umap.R"))

expr_data <- NULL
run_test("load_csv_expression", {
  expr_data <<- load_csv_expression(file.path(test_data_dir, "csv_expression"))
  stopifnot(length(expr_data) == 6)
})

combined_matrix <- NULL
run_test("combine_expression_matrices", {
  combined_matrix <<- combine_expression_matrices(expr_data, markers, max_cells = 200)
  stopifnot(!is.null(combined_matrix))
  stopifnot(ncol(combined_matrix) == 6)
})

if (!is.null(combined_matrix)) {
  sample_labels_expr <- rep(names(expr_data), each = ceiling(nrow(combined_matrix) / length(expr_data)))
  sample_labels_expr <- sample_labels_expr[1:nrow(combined_matrix)]

  run_test("plot_umap_from_matrix (CSV expression)", {
    plot_umap_from_matrix(
      expr_matrix = combined_matrix,
      sample_labels = sample_labels_expr,
      markers_of_interest = markers,
      output_dir = file.path(output_dir, "path1_csv_umap"),
      n_neighbors = 5,
      seed = 42
    )
    stopifnot(file.exists(file.path(output_dir, "path1_csv_umap", "umap_sample.png")))
  })
}

# =============================================================================
# PATH 2: CSV Expression → combine → matrix-based Density
# =============================================================================
cat("\n--- Path 2: CSV Expression → Matrix Density ---\n")

source(file.path(base_dir, "src", "visualization", "density_plots.R"))

if (!is.null(combined_matrix)) {
  run_test("plot_density_from_matrix (CSV expression)", {
    plot_density_from_matrix(
      expr_matrix = combined_matrix,
      sample_labels = sample_labels_expr,
      markers = c("CD3", "CD4"),
      output_dir = file.path(output_dir, "path2_csv_density")
    )
    stopifnot(file.exists(file.path(output_dir, "path2_csv_density", "density_CD3.png")))
  })
}

# =============================================================================
# PATH 3: CSV Expression → combine → PCA
# =============================================================================
cat("\n--- Path 3: CSV Expression → PCA ---\n")

source(file.path(base_dir, "src", "visualization", "pca.R"))

if (!is.null(combined_matrix)) {
  run_test("plot_pca (CSV expression)", {
    plot_pca(
      expr_matrix = combined_matrix,
      sample_labels = sample_labels_expr,
      output_dir = file.path(output_dir, "path3_csv_pca")
    )
    stopifnot(file.exists(file.path(output_dir, "path3_csv_pca", "pca_plot.png")))
  })
}

# =============================================================================
# PATH 4: CSV Expression → combine → Heatmap
# =============================================================================
cat("\n--- Path 4: CSV Expression → Heatmap ---\n")

source(file.path(base_dir, "src", "visualization", "heatmaps.R"))

if (!is.null(combined_matrix)) {
  run_test("plot_heatmap (CSV expression)", {
    plot_heatmap(
      expr_matrix = t(combined_matrix),
      output_path = file.path(output_dir, "path4_csv_heatmap", "heatmap.png"),
      title = "Test: CSV Expression Heatmap"
    )
    stopifnot(file.exists(file.path(output_dir, "path4_csv_heatmap", "heatmap.png")))
  })
}

# =============================================================================
# PATH 5: CSV Gating → split → tracking matrix → long → barplot
# =============================================================================
cat("\n--- Path 5: CSV Gating → Split → Barplot ---\n")

source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
source(file.path(base_dir, "src", "visualization", "barplots.R"))

split_result <- NULL
run_test("split_by_cluster (CSV gating)", {
  split_result <<- split_by_cluster(
    input_dir = file.path(test_data_dir, "csv_gating"),
    output_dir = file.path(output_dir, "path5_split"),
    cluster_regex = "Clustering annotated",
    markers_of_interest = markers
  )
  stopifnot(!is.null(split_result$tracking_matrix))
  stopifnot(nrow(split_result$tracking_matrix) > 0)
})

if (!is.null(split_result) && nrow(split_result$tracking_matrix) > 0) {
  run_test("tracking_matrix_to_long", {
    long_data <- tracking_matrix_to_long(split_result$tracking_matrix)
    stopifnot(nrow(long_data) > 0)
    stopifnot(all(c("file", "cluster", "count") %in% names(long_data)))
  })

  run_test("plot_composition_barplot (from tracking matrix)", {
    long_data <- tracking_matrix_to_long(split_result$tracking_matrix)
    plot_composition_barplot(
      data = long_data,
      output_dir = file.path(output_dir, "path5_barplot"),
      x_col = "cluster",
      fill_col = "file",
      title = "Test: Composition Barplot"
    )
  })
}

# =============================================================================
# PATH 6: CSV Gating → build_abundance_matrix → abundance plot
# =============================================================================
cat("\n--- Path 6: CSV Gating → Abundance Matrix → Abundance Plot ---\n")

source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))
source(file.path(base_dir, "src", "analysis", "abundance_plotting.R"))

abund_result <- NULL
run_test("build_abundance_matrix (standalone)", {
  abund_result <<- build_abundance_matrix(
    data_dir = file.path(test_data_dir, "csv_gating"),
    cluster_regex = "Clustering annotated",
    min_cluster_size = 1
  )
  stopifnot(!is.null(abund_result$abundance_df))
  stopifnot(length(abund_result$valid_clusters) > 0)
})

if (!is.null(abund_result)) {
  run_test("build_tracking_matrix (from CSV gating)", {
    tracking_matrix <- build_tracking_matrix(
      input_dir = file.path(test_data_dir, "csv_gating"),
      cluster_regex = "Clustering annotated",
      exclude_patterns = character(0),
      output_dir = file.path(output_dir, "path6_abundance")
    )
    stopifnot(nrow(tracking_matrix) > 0)
  })

  run_test("plot_abundance (from tracking matrix)", {
    tryCatch({
      tracking_matrix <- build_tracking_matrix(
        input_dir = file.path(test_data_dir, "csv_gating"),
        cluster_regex = "Clustering annotated",
        exclude_patterns = character(0),
        output_dir = file.path(output_dir, "path6_abundance_plots")
      )
      plot_abundance(
        output_dir = file.path(output_dir, "path6_abundance_plots"),
        tracking_matrix = tracking_matrix
      )
    }, error = function(e) {
      # abundance plotting may need ggdendro; check if it's a package issue
      if (grepl("package", e$message, ignore.case = TRUE)) {
        cat("  (skipped: missing optional package)\n")
      } else {
        stop(e)
      }
    })
  })
}

# =============================================================================
# PATH 7: CSV Gating + conditions → GLMM → volcano
# =============================================================================
cat("\n--- Path 7: CSV Gating → GLMM → Volcano ---\n")

run_test("build_abundance_matrix with conditions", {
  abund_cond <- build_abundance_matrix(
    data_dir = file.path(test_data_dir, "csv_gating"),
    conditions_csv = file.path(test_data_dir, "csv_gating_conditions.csv"),
    cluster_regex = "Clustering annotated",
    condition_levels = c("Control", "Treatment"),
    min_cluster_size = 1
  )
  stopifnot("condition" %in% names(abund_cond$abundance_df))
})

run_test("plot_volcano (standalone)", {
  # Create mock results for volcano plot test
  mock_results <- data.frame(
    cluster = c("CD4+", "CD8+", "B_cells"),
    log2fc = c(1.5, -0.8, 0.2),
    p_adjusted = c(0.001, 0.03, 0.5),
    significant = c(TRUE, TRUE, FALSE)
  )
  plot_volcano(
    results_df = mock_results,
    output_dir = file.path(output_dir, "path7_volcano"),
    title = "Test: Volcano Plot"
  )
  stopifnot(file.exists(file.path(output_dir, "path7_volcano", "volcano_plot.png")))
})

# =============================================================================
# PATH 8: CSV Expression → flowFrame bridge → UMAP (flowFrame-based)
# =============================================================================
cat("\n--- Path 8: CSV Expression → flowFrame Bridge → UMAP ---\n")

source(file.path(base_dir, "src", "loading", "create_flowframes.R"))

if (!is.null(expr_data)) {
  ff_list <- NULL
  run_test("create_flowframes_from_matrices (bridge)", {
    ff_list <<- create_flowframes_from_matrices(
      expr_list = expr_data,
      markers_of_interest = markers
    )
    stopifnot(length(ff_list) == length(expr_data))
    stopifnot(inherits(ff_list[[1]], "flowFrame"))
  })

  if (!is.null(ff_list)) {
    run_test("plot_umap (flowFrame input from CSV bridge)", {
      plot_umap(
        flow_frames = ff_list,
        sample_names = names(ff_list),
        markers_of_interest = markers,
        output_dir = file.path(output_dir, "path8_ff_umap"),
        n_neighbors = 5,
        seed = 42
      )
      stopifnot(file.exists(file.path(output_dir, "path8_ff_umap", "umap_sample.png")))
    })
  }
}

# =============================================================================
# PATH 9: CSV Expression → clustering → gating CSV → split → expression matrices
# =============================================================================
cat("\n--- Path 9: CSV Expression → Cluster → Gating CSV → Split → Expr Matrices ---\n")

source(file.path(base_dir, "src", "clustering", "leiden.R"))

if (!is.null(combined_matrix)) {
  cluster_result <- NULL
  run_test("run_leiden (on CSV expression matrix)", {
    cluster_result <<- run_leiden(combined_matrix, n_neighbors = 5, resolution = 0.5, seed = 42)
    stopifnot(length(cluster_result) == nrow(combined_matrix))
  })

  if (!is.null(cluster_result)) {
    run_test("clusters_to_gating_csv (from expression + clusters)", {
      clusters_to_gating_csv(
        combined_matrix = combined_matrix,
        sample_labels = sample_labels_expr,
        cluster_ids = cluster_result,
        output_dir = file.path(output_dir, "path9_gating_csv")
      )
      csvs <- list.files(file.path(output_dir, "path9_gating_csv"), pattern = "\\.csv$")
      stopifnot(length(csvs) > 0)
    })

    run_test("split_by_cluster (from generated gating CSVs)", {
      split_result9 <- split_by_cluster(
        input_dir = file.path(output_dir, "path9_gating_csv"),
        output_dir = file.path(output_dir, "path9_split"),
        cluster_regex = "Clustering annotated",
        markers_of_interest = markers
      )
      stopifnot(!is.null(split_result9$tracking_matrix))
    })
  }
}

# =============================================================================
# PATH 10: CSV Expression → matrix UMAP + cluster overlay
# =============================================================================
cat("\n--- Path 10: CSV Expression → UMAP with Cluster Overlay ---\n")

if (!is.null(combined_matrix) && !is.null(cluster_result)) {
  run_test("plot_umap_from_matrix with cluster_labels", {
    plot_umap_from_matrix(
      expr_matrix = combined_matrix,
      sample_labels = sample_labels_expr,
      markers_of_interest = markers,
      cluster_labels = cluster_result,
      output_dir = file.path(output_dir, "path10_umap_clusters"),
      n_neighbors = 5,
      seed = 42
    )
    stopifnot(file.exists(file.path(output_dir, "path10_umap_clusters", "umap_cluster.png")))
  })
}

# =============================================================================
# PATH 11: CSV Expression → combine → Leiden + UMAP → gating CSV → abundance
# =============================================================================
cat("\n--- Path 11: Full CSV Expression → Cluster → Abundance Chain ---\n")

if (!is.null(combined_matrix) && !is.null(cluster_result)) {
  run_test("full chain: CSV → cluster → gating CSV → abundance matrix", {
    # Generate gating CSVs
    clusters_to_gating_csv(
      combined_matrix = combined_matrix,
      sample_labels = sample_labels_expr,
      cluster_ids = cluster_result,
      output_dir = file.path(output_dir, "path11_chain_gating")
    )

    # Build abundance matrix from generated CSVs
    abund11 <- build_abundance_matrix(
      data_dir = file.path(output_dir, "path11_chain_gating"),
      cluster_regex = "Clustering annotated",
      min_cluster_size = 1
    )
    stopifnot(nrow(abund11$abundance_df) > 0)
    stopifnot(length(abund11$valid_clusters) > 0)
  })
}

# =============================================================================
# PATH 12: CSV Expression → combine → Leiden → pre-computed UMAP reuse
# =============================================================================
cat("\n--- Path 12: Pre-computed UMAP Coordinates Reuse ---\n")

if (!is.null(combined_matrix)) {
  run_test("plot_umap_from_matrix with pre-computed coords", {
    # Compute UMAP once
    coords <- run_umap(combined_matrix, n_neighbors = 5, seed = 42)

    # Reuse in plot
    plot_umap_from_matrix(
      expr_matrix = combined_matrix,
      sample_labels = sample_labels_expr,
      markers_of_interest = c("CD3", "CD4"),
      umap_coords = coords,
      output_dir = file.path(output_dir, "path12_precomputed_umap"),
      seed = 42
    )
    stopifnot(file.exists(file.path(output_dir, "path12_precomputed_umap", "umap_sample.png")))
  })
}

# =============================================================================
# PATH 13: CSV Gating → cluster heatmap
# =============================================================================
cat("\n--- Path 13: CSV Gating → Cluster Heatmap ---\n")

run_test("plot_cluster_heatmap (from gating CSV)", {
  dt <- data.table::fread(file.path(test_data_dir, "csv_gating", "sample_001.csv"))
  # Add a cluster column from gating data
  dt$cluster <- ifelse(dt[["CD4+ (Clustering annotated)"]] == 1, "CD4+",
                        ifelse(dt[["CD8+ (Clustering annotated)"]] == 1, "CD8+", "Other"))
  plot_cluster_heatmap(
    expr_data = as.data.frame(dt),
    cluster_col = "cluster",
    marker_cols = markers,
    output_path = file.path(output_dir, "path13_cluster_heatmap", "heatmap.png"),
    title = "Test: Cluster Expression Heatmap"
  )
})

# =============================================================================
# Summary
# =============================================================================
cat("\n=== Test Summary ===\n")
cat(sprintf("  PASSED: %d\n", pass))
cat(sprintf("  FAILED: %d\n", fail))
cat(sprintf("  TOTAL:  %d\n", pass + fail))

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
