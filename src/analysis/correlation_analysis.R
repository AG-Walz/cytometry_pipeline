#!/usr/bin/env Rscript
# =============================================================================
# analysis/correlation_analysis.R — Marker-Marker Correlation Analysis
# =============================================================================
# Computes pairwise marker correlations (Spearman or Pearson) with p-values,
# optionally stratified by cluster. Produces correlation heatmaps and
# identifies significant marker pairs.
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
source(file.path(shared_dir, "plot_utils.R"))

#' Compute marker-marker correlation matrices
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R)
#' Output to:  plot_heatmap() (visualization/heatmaps.R)
#'
#' @param expr_matrix Numeric matrix (cells x markers)
#' @param method "spearman" or "pearson"
#' @param output_dir Output directory (for plots + CSV). NULL = no saving.
#' @param cluster_labels Optional: compute per-cluster correlations
#' @param p_adjust_method P-value adjustment method
#' @return list(cor_matrix, p_matrix, significant_pairs)
run_correlation_analysis <- function(expr_matrix,
                                      method = "spearman",
                                      output_dir = NULL,
                                      cluster_labels = NULL,
                                      p_adjust_method = "BH") {
  require(data.table)

  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  n_markers <- ncol(expr_matrix)
  marker_names <- colnames(expr_matrix)
  if (is.null(marker_names)) {
    marker_names <- paste0("V", seq_len(n_markers))
    colnames(expr_matrix) <- marker_names
  }

  cat("Computing", method, "correlations for", n_markers, "markers\n")
  cat("Cells:", nrow(expr_matrix), "\n")

  if (!is.null(output_dir)) make_outdir(output_dir)

  # If no cluster labels, compute global correlations
  if (is.null(cluster_labels)) {
    result <- .compute_cor_with_pvalues(expr_matrix, method, p_adjust_method)

    if (!is.null(output_dir)) {
      .save_correlation_results(result, output_dir, prefix = "global", method = method)
    }

    cat("Significant pairs (padj < 0.05):", nrow(result$significant_pairs), "\n")
    return(result)
  }

  # Per-cluster correlations
  cluster_labels <- as.character(cluster_labels)
  stopifnot(length(cluster_labels) == nrow(expr_matrix))
  clusters <- sort(unique(cluster_labels))
  cat("Computing per-cluster correlations for", length(clusters), "clusters\n")

  # Global result first
  global_result <- .compute_cor_with_pvalues(expr_matrix, method, p_adjust_method)
  if (!is.null(output_dir)) {
    .save_correlation_results(global_result, output_dir, prefix = "global", method = method)
  }

  # Per-cluster results
  cluster_results <- list()
  for (cl in clusters) {
    cl_mask <- cluster_labels == cl
    cl_matrix <- expr_matrix[cl_mask, , drop = FALSE]

    if (nrow(cl_matrix) < 10) {
      cat("  Cluster", cl, "- skipped (", sum(cl_mask), "cells)\n")
      next
    }

    cat("  Cluster", cl, "(", sum(cl_mask), "cells)...")
    cl_result <- tryCatch({
      .compute_cor_with_pvalues(cl_matrix, method, p_adjust_method)
    }, error = function(e) {
      cat(" error:", e$message)
      NULL
    })

    if (!is.null(cl_result)) {
      cluster_results[[cl]] <- cl_result
      if (!is.null(output_dir)) {
        safe_cl <- create_safe_filename(cl)
        .save_correlation_results(cl_result, output_dir,
                                   prefix = paste0("cluster_", safe_cl),
                                   method = method)
      }
      cat(" ", nrow(cl_result$significant_pairs), "sig. pairs\n")
    } else {
      cat("\n")
    }
  }

  result <- list(
    cor_matrix = global_result$cor_matrix,
    p_matrix = global_result$p_matrix,
    significant_pairs = global_result$significant_pairs,
    cluster_results = cluster_results
  )

  cat("\nCorrelation analysis complete.\n")
  result
}

# --- Internal: Compute correlation matrix with pairwise p-values ---
.compute_cor_with_pvalues <- function(expr_matrix, method, p_adjust_method) {
  n <- ncol(expr_matrix)
  marker_names <- colnames(expr_matrix)

  # Correlation matrix
  cor_matrix <- cor(expr_matrix, method = method, use = "pairwise.complete.obs")

  # P-value matrix (from pairwise cor.test)
  p_matrix <- matrix(NA, nrow = n, ncol = n,
                      dimnames = list(marker_names, marker_names))

  pairs_list <- list()
  pair_idx <- 0

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      x <- expr_matrix[, i]
      y <- expr_matrix[, j]

      # Remove NAs pairwise
      valid <- !is.na(x) & !is.na(y)
      if (sum(valid) < 5) {
        p_matrix[i, j] <- NA
        p_matrix[j, i] <- NA
        next
      }

      ct <- tryCatch({
        cor.test(x[valid], y[valid], method = method)
      }, error = function(e) NULL)

      if (!is.null(ct)) {
        p_matrix[i, j] <- ct$p.value
        p_matrix[j, i] <- ct$p.value

        pair_idx <- pair_idx + 1
        pairs_list[[pair_idx]] <- data.frame(
          marker_1 = marker_names[i],
          marker_2 = marker_names[j],
          correlation = cor_matrix[i, j],
          pvalue = ct$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  diag(p_matrix) <- 0

  # Combine pairs and adjust p-values
  if (length(pairs_list) > 0) {
    all_pairs <- do.call(rbind, pairs_list)
    all_pairs$padj <- p.adjust(all_pairs$pvalue, method = p_adjust_method)
    all_pairs$abs_cor <- abs(all_pairs$correlation)

    # Significant pairs
    sig_pairs <- all_pairs[!is.na(all_pairs$padj) & all_pairs$padj < 0.05, ]
    sig_pairs <- sig_pairs[order(-sig_pairs$abs_cor), ]
  } else {
    all_pairs <- data.frame(
      marker_1 = character(0), marker_2 = character(0),
      correlation = numeric(0), pvalue = numeric(0),
      padj = numeric(0), abs_cor = numeric(0)
    )
    sig_pairs <- all_pairs
  }

  list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    all_pairs = all_pairs,
    significant_pairs = sig_pairs
  )
}

# --- Internal: Save correlation results and plot heatmap ---
.save_correlation_results <- function(result, output_dir, prefix, method) {
  require(data.table)
  require(ggplot2)

  # Save correlation matrix CSV
  cor_df <- as.data.frame(result$cor_matrix)
  cor_df <- cbind(marker = rownames(cor_df), cor_df)
  fwrite(cor_df, file.path(output_dir, paste0(prefix, "_correlation_matrix.csv")), sep = ";")

  # Save p-value matrix CSV
  p_df <- as.data.frame(result$p_matrix)
  p_df <- cbind(marker = rownames(p_df), p_df)
  fwrite(p_df, file.path(output_dir, paste0(prefix, "_pvalue_matrix.csv")), sep = ";")

  # Save significant pairs
  if (nrow(result$significant_pairs) > 0) {
    fwrite(result$significant_pairs,
           file.path(output_dir, paste0(prefix, "_significant_pairs.csv")), sep = ";")
  }

  # Plot correlation heatmap
  tryCatch({
    .plot_correlation_heatmap(result$cor_matrix, output_dir, prefix, method)
  }, error = function(e) {
    cat("  Heatmap plot failed:", e$message, "\n")
  })
}

# --- Internal: Correlation heatmap using ggplot2::geom_tile ---
.plot_correlation_heatmap <- function(cor_matrix, output_dir, prefix, method) {
  require(ggplot2)
  require(reshape2)

  # Hierarchical clustering for ordering
  hc <- tryCatch({
    hclust(as.dist(1 - abs(cor_matrix)), method = "ward.D2")
  }, error = function(e) NULL)

  if (!is.null(hc)) {
    ordered_names <- colnames(cor_matrix)[hc$order]
    cor_matrix <- cor_matrix[ordered_names, ordered_names]
  }

  long_cor <- reshape2::melt(cor_matrix)
  names(long_cor) <- c("Marker1", "Marker2", "Correlation")
  long_cor$Marker1 <- factor(long_cor$Marker1, levels = colnames(cor_matrix))
  long_cor$Marker2 <- factor(long_cor$Marker2, levels = rev(colnames(cor_matrix)))

  p <- ggplot(long_cor, aes(x = Marker1, y = Marker2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                          midpoint = 0, limits = c(-1, 1),
                          name = paste0(tools::toTitleCase(method), "\nCorrelation")) +
    theme_minimal() +
    labs(title = paste0("Marker Correlation Matrix (", method, ")"),
         x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    coord_fixed()

  save_plot_formats(p, file.path(output_dir, paste0(prefix, "_correlation_heatmap")),
                    width = max(10, ncol(cor_matrix) * 0.5),
                    height = max(10, ncol(cor_matrix) * 0.5))
}
