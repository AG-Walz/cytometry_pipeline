#!/usr/bin/env Rscript
# =============================================================================
# analysis/diffcyt_analysis.R — diffcyt Framework for Differential Discovery
# =============================================================================
# Wraps the diffcyt Bioconductor package for differential abundance (DA) and
# differential state (DS) analysis using FlowSOM clustering and edgeR/voom/GLMM
# or limma/LMM statistical frameworks.
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
      # Rscript src/foo/bar.R carries no source frame; recover from --file=.
      if (is.null(of)) {
        a <- commandArgs(FALSE); m <- grep("^--file=", a)
        if (length(m)) of <- sub("^--file=", "", a[m[[1]]])
      }
      d <- if (!is.null(of)) dirname(normalizePath(of, mustWork = FALSE)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

# Guarded, NOT unconditional: nextflow/bin/*_cli.R and scripts/*.R set shared_dir
# themselves (the --src_dir convention), and clobbering it here would send every
# source() below to whichever tree this file happens to live in rather than the
# one the caller chose.
if (!exists("shared_dir")) shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Run diffcyt framework for differential discovery
#'
#' Input from: flowFrame lists from load_fcs_raw() (loading/load_fcs_raw.R),
#'             build_combined_expr_matrix() (R_shared/flowframe_utils.R)
#' Output to:  Standalone results CSV, plot_volcano() (analysis/differential_abundance.R)
#'
#' @param flow_frames List of flowFrame objects
#' @param experiment_info data.frame with sample_id, group_id, patient_id columns
#' @param markers_of_interest Character vector of markers
#' @param analysis_type "DA" (differential abundance) or "DS" (differential state)
#' @param method_DA "diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM"
#' @param method_DS "diffcyt-DS-limma", "diffcyt-DS-LMM"
#' @param n_clusters Number of clusters (FlowSOM meta-clusters)
#' @param seed Random seed
#' @return list(results, top_table, cluster_assignments)
run_diffcyt <- function(flow_frames,
                         experiment_info,
                         markers_of_interest,
                         analysis_type = "DA",
                         method_DA = "diffcyt-DA-edgeR",
                         method_DS = "diffcyt-DS-limma",
                         n_clusters = 20,
                         seed = 42) {

  # Check diffcyt availability
  if (!requireNamespace("diffcyt", quietly = TRUE)) {
    stop("Package 'diffcyt' is required but not installed.\n",
         "Install with: BiocManager::install('diffcyt')")
  }
  require(diffcyt)
  require(flowCore)

  set.seed(seed)

  # Validate inputs
  stopifnot(is.list(flow_frames))
  stopifnot(all(sapply(flow_frames, is, "flowFrame")))
  stopifnot(is.data.frame(experiment_info))
  stopifnot(all(c("sample_id", "group_id") %in% names(experiment_info)))
  stopifnot(analysis_type %in% c("DA", "DS"))
  stopifnot(length(flow_frames) == nrow(experiment_info))

  cat("Running diffcyt", analysis_type, "analysis\n")
  cat("Samples:", length(flow_frames), "\n")
  cat("Groups:", paste(unique(experiment_info$group_id), collapse = ", "), "\n")
  cat("Markers:", length(markers_of_interest), "\n")
  cat("Target clusters:", n_clusters, "\n")

  # --- Step 1: Prepare data ---
  cat("  Step 1/5: Preparing data...\n")
  d_input <- tryCatch({
    diffcyt::prepareData(flow_frames, experiment_info)
  }, error = function(e) {
    stop("prepareData() failed: ", e$message)
  })

  # --- Step 2: Transform data ---
  cat("  Step 2/5: Transforming data...\n")
  d_transformed <- tryCatch({
    diffcyt::transformData(d_input)
  }, error = function(e) {
    stop("transformData() failed: ", e$message)
  })

  # --- Step 3: Generate clusters (FlowSOM) ---
  cat("  Step 3/5: Generating clusters (FlowSOM, k =", n_clusters, ")...\n")

  # Build marker_info for clustering
  all_markers <- colnames(flowCore::exprs(flow_frames[[1]]))
  marker_info <- data.frame(
    channel_name = all_markers,
    marker_name = all_markers,
    marker_class = ifelse(all_markers %in% markers_of_interest, "type", "state"),
    stringsAsFactors = FALSE
  )

  d_clustered <- tryCatch({
    diffcyt::generateClusters(d_transformed,
                               cols_clustering = markers_of_interest,
                               xdim = 10, ydim = 10,
                               meta_clustering = TRUE,
                               meta_k = n_clusters,
                               seed_clustering = seed)
  }, error = function(e) {
    stop("generateClusters() failed: ", e$message)
  })

  # Extract cluster assignments
  cluster_assignments <- tryCatch({
    as.data.frame(SummarizedExperiment::rowData(d_clustered))
  }, error = function(e) {
    cat("  Warning: could not extract cluster assignments\n")
    NULL
  })

  # --- Step 4: Calculate counts/medians ---
  cat("  Step 4/5: Calculating summary statistics...\n")

  if (analysis_type == "DA") {
    d_counts <- tryCatch({
      diffcyt::calculateCounts(d_clustered)
    }, error = function(e) {
      stop("calculateCounts() failed: ", e$message)
    })
  } else {
    d_medians <- tryCatch({
      diffcyt::calculateMedians(d_clustered)
    }, error = function(e) {
      stop("calculateMedians() failed: ", e$message)
    })
  }

  # --- Step 5: Statistical testing ---
  cat("  Step 5/5: Running statistical tests...\n")

  # Create design and contrast matrices
  design <- .create_design_matrix(experiment_info)
  contrast <- .create_contrast(design, experiment_info)

  if (analysis_type == "DA") {
    cat("  DA method:", method_DA, "\n")
    results <- tryCatch({
      switch(method_DA,
        "diffcyt-DA-edgeR" = diffcyt::testDA_edgeR(d_counts, design, contrast),
        "diffcyt-DA-voom"  = diffcyt::testDA_voom(d_counts, design, contrast),
        "diffcyt-DA-GLMM"  = diffcyt::testDA_GLMM(d_counts, experiment_info,
                                                     formula = ~ group_id + (1 | patient_id),
                                                     contrast = contrast),
        stop("Unknown DA method: ", method_DA)
      )
    }, error = function(e) {
      # Fallback: try edgeR if specified method fails
      if (method_DA != "diffcyt-DA-edgeR") {
        warning("DA method '", method_DA, "' failed: ", e$message,
                ". Falling back to diffcyt-DA-edgeR.")
        cat("  ", method_DA, "failed:", e$message, "\n")
        cat("  Falling back to diffcyt-DA-edgeR...\n")
        method_DA <<- "diffcyt-DA-edgeR"
        tryCatch({
          diffcyt::testDA_edgeR(d_counts, design, contrast)
        }, error = function(e2) {
          stop("Fallback DA-edgeR also failed: ", e2$message)
        })
      } else {
        stop("DA-edgeR failed: ", e$message)
      }
    })
  } else {
    cat("  DS method:", method_DS, "\n")
    results <- tryCatch({
      switch(method_DS,
        "diffcyt-DS-limma" = diffcyt::testDS_limma(d_clustered, d_medians,
                                                      design, contrast,
                                                      markers_test = markers_of_interest),
        "diffcyt-DS-LMM"   = diffcyt::testDS_LMM(d_clustered, d_medians,
                                                     experiment_info,
                                                     formula = ~ group_id + (1 | patient_id),
                                                     contrast = contrast,
                                                     markers_test = markers_of_interest),
        stop("Unknown DS method: ", method_DS)
      )
    }, error = function(e) {
      # Fallback: try limma if specified method fails
      if (method_DS != "diffcyt-DS-limma") {
        warning("DS method '", method_DS, "' failed: ", e$message,
                ". Falling back to diffcyt-DS-limma.")
        cat("  ", method_DS, "failed:", e$message, "\n")
        cat("  Falling back to diffcyt-DS-limma...\n")
        method_DS <<- "diffcyt-DS-limma"
        tryCatch({
          diffcyt::testDS_limma(d_clustered, d_medians,
                                 design, contrast,
                                 markers_test = markers_of_interest)
        }, error = function(e2) {
          stop("Fallback DS-limma also failed: ", e2$message)
        })
      } else {
        stop("DS-limma failed: ", e$message)
      }
    })
  }

  # Extract top table
  top_table <- tryCatch({
    tt <- diffcyt::topTable(results, all = TRUE)
    as.data.frame(tt)
  }, error = function(e) {
    cat("  Warning: topTable() extraction failed:", e$message, "\n")
    data.frame()
  })

  # Record the actual method used (may differ from requested if fallback occurred)
  actual_method <- if (analysis_type == "DA") method_DA else method_DS
  if (nrow(top_table) > 0) {
    top_table$method_used <- actual_method
  }

  n_sig <- sum(top_table$p_adj < 0.05, na.rm = TRUE)
  cat("\nDiffcyt analysis complete.\n")
  cat("Method used:", actual_method, "\n")
  cat("Total clusters tested:", nrow(top_table), "\n")
  cat("Significant (padj < 0.05):", n_sig, "\n")

  list(
    results = results,
    top_table = top_table,
    cluster_assignments = cluster_assignments,
    method_used = actual_method
  )
}

#' Run diffcyt and save results to files
#'
#' Convenience wrapper that runs run_diffcyt() and saves the top table
#' and cluster assignments to the output directory.
#'
#' @param flow_frames List of flowFrame objects
#' @param experiment_info data.frame with sample_id, group_id, patient_id columns
#' @param markers_of_interest Character vector of markers
#' @param output_dir Output directory
#' @param analysis_type "DA" or "DS"
#' @param method_DA DA method string
#' @param method_DS DS method string
#' @param n_clusters Number of FlowSOM meta-clusters
#' @param seed Random seed
#' @return list from run_diffcyt()
run_diffcyt_and_save <- function(flow_frames,
                                  experiment_info,
                                  markers_of_interest,
                                  output_dir,
                                  analysis_type = "DA",
                                  method_DA = "diffcyt-DA-edgeR",
                                  method_DS = "diffcyt-DS-limma",
                                  n_clusters = 20,
                                  seed = 42) {
  require(data.table)
  make_outdir(output_dir)

  result <- run_diffcyt(
    flow_frames = flow_frames,
    experiment_info = experiment_info,
    markers_of_interest = markers_of_interest,
    analysis_type = analysis_type,
    method_DA = method_DA,
    method_DS = method_DS,
    n_clusters = n_clusters,
    seed = seed
  )

  # Save top table
  if (nrow(result$top_table) > 0) {
    fwrite(result$top_table,
           file.path(output_dir, paste0("diffcyt_", tolower(analysis_type), "_results.csv")),
           sep = ";")
    cat("Results saved to:", output_dir, "\n")
  }

  # Save cluster assignments
  if (!is.null(result$cluster_assignments)) {
    fwrite(result$cluster_assignments,
           file.path(output_dir, "diffcyt_cluster_assignments.csv"),
           sep = ";")
  }

  result
}

# --- Internal: Create design matrix from experiment_info ---
.create_design_matrix <- function(experiment_info) {
  group_factor <- factor(experiment_info$group_id)
  design <- model.matrix(~ group_factor)
  colnames(design) <- c("Intercept", levels(group_factor)[-1])
  design
}

# --- Internal: Create contrast vector ---
.create_contrast <- function(design, experiment_info) {
  n_cols <- ncol(design)
  contrast <- rep(0, n_cols)
  # Test the second group vs reference (intercept)
  if (n_cols >= 2) {
    contrast[2] <- 1
  }
  contrast
}

# --- CLI execution ---
if (sys.nframe() == 0) {
  cat("diffcyt_analysis.R: This script is intended to be sourced and called\n")
  cat("programmatically with flowFrame objects. See run_diffcyt() documentation.\n")
  cat("\nExample usage:\n")
  cat("  source('diffcyt_analysis.R')\n")
  cat("  result <- run_diffcyt(flow_frames, experiment_info, markers)\n")
}
