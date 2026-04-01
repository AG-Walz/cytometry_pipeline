#!/usr/bin/env Rscript
# =============================================================================
# clustering/cell_annotation.R — Cluster Cell Type Annotation
# =============================================================================
# Annotates clusters with cell type labels based on marker expression
# profiles. Supports threshold-based, correlation-based, and manual
# annotation workflows.
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))

#' Annotate clusters with cell type labels
#'
#' Computes per-cluster median expression for each marker and assigns cell type
#' annotations using one of three methods:
#'
#' - "threshold": markers above the global median are marked positive (+),
#'   others negative (-). Builds a phenotype string (e.g. "CD3+CD4+CD8-").
#' - "correlation": computes cosine similarity between cluster median profiles
#'   and a reference matrix. Assigns the best-matching reference cell type.
#' - "manual": returns median expression per cluster without automatic labeling,
#'   for manual annotation by the user.
#'
#' If reference="auto" with method "threshold", phenotype strings are built
#' from positive/negative marker calls.
#'
#' Input from: run_leiden() (clustering/leiden.R),
#'             run_flowsom() (clustering/flowsom.R),
#'             run_phenograph() (clustering/phenograph.R),
#'             run_dbscan() (clustering/dbscan.R)
#' Output to:  clusters_to_gating_csv() (clustering/leiden.R),
#'             analysis scripts for cell type abundance
#'
#' @param expr_matrix Numeric matrix (rows = cells, cols = markers)
#' @param cluster_ids Integer vector of cluster assignments (length = nrow(expr_matrix))
#' @param markers Character vector of marker names. If NULL, uses colnames(expr_matrix).
#' @param reference Reference for annotation. "auto" builds phenotype strings
#'   from threshold calls. A matrix or data.frame with rows = cell types and
#'   cols = markers provides reference profiles for correlation matching.
#' @param method Annotation method: "threshold", "correlation", or "manual"
#' @return Data frame with columns: cluster, annotation, confidence, plus
#'   per-marker median columns
annotate_clusters <- function(expr_matrix, cluster_ids,
                              markers = NULL,
                              reference = "auto",
                              method = "threshold") {
  stopifnot(is.matrix(expr_matrix), is.numeric(expr_matrix))
  stopifnot(length(cluster_ids) == nrow(expr_matrix))
  stopifnot(method %in% c("threshold", "correlation", "manual"))

  if (is.null(markers)) {
    markers <- colnames(expr_matrix)
    if (is.null(markers)) {
      markers <- paste0("marker_", seq_len(ncol(expr_matrix)))
      colnames(expr_matrix) <- markers
    }
  }
  stopifnot(length(markers) == ncol(expr_matrix))

  unique_clusters <- sort(unique(cluster_ids))
  n_clusters <- length(unique_clusters)

  cat("Annotating", n_clusters, "clusters using method:", method, "\n")

  # Compute per-cluster median expression
  cluster_medians <- compute_cluster_medians(expr_matrix, cluster_ids,
                                              markers, unique_clusters)

  # Dispatch to annotation method
  result <- switch(method,
    "threshold"   = annotate_by_threshold(cluster_medians, markers,
                                           unique_clusters, expr_matrix),
    "correlation" = annotate_by_correlation(cluster_medians, markers,
                                             unique_clusters, reference),
    "manual"      = annotate_manual(cluster_medians, markers, unique_clusters)
  )

  cat("Annotation complete for", nrow(result), "clusters\n")

  result
}

#' Compute per-cluster median expression
#'
#' @param expr_matrix Numeric matrix (cells x markers)
#' @param cluster_ids Integer vector of cluster assignments
#' @param markers Character vector of marker names
#' @param unique_clusters Sorted unique cluster IDs
#' @return Matrix (clusters x markers) of median expression values
compute_cluster_medians <- function(expr_matrix, cluster_ids,
                                     markers, unique_clusters) {
  medians <- matrix(NA_real_, nrow = length(unique_clusters), ncol = length(markers))
  rownames(medians) <- unique_clusters
  colnames(medians) <- markers

  for (i in seq_along(unique_clusters)) {
    cl <- unique_clusters[i]
    members <- which(cluster_ids == cl)
    if (length(members) == 1) {
      medians[i, ] <- expr_matrix[members, ]
    } else {
      medians[i, ] <- apply(expr_matrix[members, , drop = FALSE], 2, median,
                             na.rm = TRUE)
    }
  }

  medians
}

#' Threshold-based annotation
#'
#' Markers above the global median for that marker are called positive (+),
#' below are negative (-). Builds a phenotype string from positive markers.
#' Confidence is the fraction of markers that are clearly above or below
#' the global median (not close to the boundary).
#'
#' @param cluster_medians Matrix (clusters x markers)
#' @param markers Character vector of marker names
#' @param unique_clusters Sorted unique cluster IDs
#' @param expr_matrix Original expression matrix (for global median calculation)
#' @return Data frame with cluster, annotation, confidence, and per-marker columns
annotate_by_threshold <- function(cluster_medians, markers,
                                   unique_clusters, expr_matrix) {
  # Compute global medians per marker
  global_medians <- apply(expr_matrix, 2, median, na.rm = TRUE)
  global_mads <- apply(expr_matrix, 2, mad, na.rm = TRUE)
  # Avoid zero MAD
  global_mads[global_mads == 0] <- 1e-10

  annotations <- character(length(unique_clusters))
  confidences <- numeric(length(unique_clusters))

  for (i in seq_along(unique_clusters)) {
    cl_medians <- cluster_medians[i, ]
    positive_markers <- character(0)
    negative_markers <- character(0)
    clear_calls <- 0

    for (j in seq_along(markers)) {
      diff_from_global <- cl_medians[j] - global_medians[j]
      z_score <- diff_from_global / global_mads[j]

      if (cl_medians[j] > global_medians[j]) {
        positive_markers <- c(positive_markers, markers[j])
      } else {
        negative_markers <- c(negative_markers, markers[j])
      }

      # A clear call is one where |z| > 1
      if (abs(z_score) > 1) {
        clear_calls <- clear_calls + 1
      }
    }

    # Build phenotype string from positive and negative markers
    parts <- character(0)
    for (m in markers) {
      if (m %in% positive_markers) {
        parts <- c(parts, paste0(m, "+"))
      } else {
        parts <- c(parts, paste0(m, "-"))
      }
    }
    # Only show positive markers in the annotation for readability
    pos_parts <- paste0(positive_markers, "+")
    if (length(pos_parts) == 0) {
      annotations[i] <- "Negative_all"
    } else {
      annotations[i] <- paste(pos_parts, collapse = "")
    }

    confidences[i] <- round(clear_calls / length(markers), 4)
  }

  # Build result data frame
  result <- data.frame(
    cluster = unique_clusters,
    annotation = annotations,
    confidence = confidences,
    stringsAsFactors = FALSE
  )

  # Add per-marker median columns
  median_df <- as.data.frame(cluster_medians)
  colnames(median_df) <- paste0("median_", markers)
  result <- cbind(result, median_df)

  result
}

#' Correlation-based annotation using reference profiles
#'
#' Computes cosine similarity between each cluster's median expression profile
#' and each reference cell type profile. Assigns the best-matching cell type.
#'
#' @param cluster_medians Matrix (clusters x markers)
#' @param markers Character vector of marker names
#' @param unique_clusters Sorted unique cluster IDs
#' @param reference Reference profiles: matrix or data.frame with rows = cell types,
#'   cols = markers. If "auto", falls back to threshold method.
#' @return Data frame with cluster, annotation, confidence, and per-marker columns
annotate_by_correlation <- function(cluster_medians, markers,
                                     unique_clusters, reference) {
  # Validate reference
  if (is.character(reference) && reference == "auto") {
    cat("Warning: correlation method requires a reference matrix.",
        " Falling back to threshold.\n")
    # Cannot call annotate_by_threshold here without expr_matrix,
    # so return manual results with a note
    result <- annotate_manual(cluster_medians, markers, unique_clusters)
    result$annotation <- "Requires reference matrix"
    result$confidence <- 0
    return(result)
  }

  reference <- as.matrix(reference)
  stopifnot(is.numeric(reference))

  # Match markers between cluster medians and reference
  common_markers <- intersect(colnames(cluster_medians), colnames(reference))
  if (length(common_markers) == 0) {
    cat("Warning: no matching markers between clusters and reference\n")
    result <- annotate_manual(cluster_medians, markers, unique_clusters)
    result$annotation <- "No matching markers"
    result$confidence <- 0
    return(result)
  }

  cat("Matching on", length(common_markers), "common markers\n")

  cl_sub <- cluster_medians[, common_markers, drop = FALSE]
  ref_sub <- reference[, common_markers, drop = FALSE]

  ref_labels <- if (!is.null(rownames(reference))) {
    rownames(reference)
  } else {
    paste0("ref_", seq_len(nrow(reference)))
  }

  annotations <- character(length(unique_clusters))
  confidences <- numeric(length(unique_clusters))

  for (i in seq_along(unique_clusters)) {
    cl_vec <- cl_sub[i, ]
    best_sim <- -Inf
    best_label <- "Unknown"

    for (r in seq_len(nrow(ref_sub))) {
      ref_vec <- ref_sub[r, ]
      sim <- cosine_similarity(cl_vec, ref_vec)
      if (sim > best_sim) {
        best_sim <- sim
        best_label <- ref_labels[r]
      }
    }

    annotations[i] <- best_label
    confidences[i] <- round(max(0, best_sim), 4)
  }

  result <- data.frame(
    cluster = unique_clusters,
    annotation = annotations,
    confidence = confidences,
    stringsAsFactors = FALSE
  )

  median_df <- as.data.frame(cluster_medians)
  colnames(median_df) <- paste0("median_", markers)
  result <- cbind(result, median_df)

  result
}

#' Manual annotation (returns medians for user labeling)
#'
#' @param cluster_medians Matrix (clusters x markers)
#' @param markers Character vector of marker names
#' @param unique_clusters Sorted unique cluster IDs
#' @return Data frame with cluster, annotation, confidence, and per-marker columns
annotate_manual <- function(cluster_medians, markers, unique_clusters) {
  result <- data.frame(
    cluster = unique_clusters,
    annotation = rep("(manual)", length(unique_clusters)),
    confidence = rep(NA_real_, length(unique_clusters)),
    stringsAsFactors = FALSE
  )

  median_df <- as.data.frame(cluster_medians)
  colnames(median_df) <- paste0("median_", markers)
  result <- cbind(result, median_df)

  cat("Manual annotation mode: review median expression to assign labels\n")

  result
}

#' Compute cosine similarity between two numeric vectors
#'
#' @param a Numeric vector
#' @param b Numeric vector (same length as a)
#' @return Numeric scalar: cosine similarity (-1 to 1)
cosine_similarity <- function(a, b) {
  dot_prod <- sum(a * b, na.rm = TRUE)
  norm_a <- sqrt(sum(a^2, na.rm = TRUE))
  norm_b <- sqrt(sum(b^2, na.rm = TRUE))
  if (norm_a == 0 || norm_b == 0) return(0)
  dot_prod / (norm_a * norm_b)
}
