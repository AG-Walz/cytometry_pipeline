#!/usr/bin/env Rscript
# =============================================================================
# harmonization/cytonorm.R — CytoNorm Batch Correction
# =============================================================================
# Quantile-based, cytometry-specific batch correction using FlowSOM-stratified
# normalization. Falls back to simple quantile normalization if CytoNorm is
# unavailable or fails.
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))

#' Batch correction using CytoNorm
#'
#' Input from: load_fcs_raw() (loading/load_fcs_raw.R),
#'             apply_compensation() (preprocessing/compensation.R)
#' Output to:  transform_flowframe() (preprocessing/transformation.R),
#'             run_leiden() (clustering/leiden.R)
#'
#' @param train_frames List of flowFrames from reference batch(es)
#' @param validation_frames List of flowFrames to correct
#' @param batch_labels Character vector of batch assignments (one per frame, across both train+validation)
#' @param markers Character vector of markers to correct
#' @param n_clusters Number of FlowSOM clusters for stratified normalization (default 10)
#' @param goal "mean" (normalize to mean of reference batches) or specific batch name
#' @param seed Random seed
#' @return list(corrected_frames, model)
run_cytonorm <- function(train_frames, validation_frames, batch_labels,
                          markers, n_clusters = 10, goal = "mean", seed = 42) {
  require(flowCore)

  # --- Validate inputs -------------------------------------------------------
  n_total <- length(train_frames) + length(validation_frames)
  if (length(batch_labels) != n_total) {
    stop("batch_labels length (", length(batch_labels), ") must equal ",
         "length(train_frames) + length(validation_frames) (", n_total, ")")
  }

  # Check that markers exist in all frames
  all_frames <- c(train_frames, validation_frames)
  for (i in seq_along(all_frames)) {
    ff_data <- explore_data(all_frames[[i]])
    available <- ff_data$markers[ff_data$mcols]
    missing <- setdiff(markers, available)
    if (length(missing) > 0) {
      stop("Frame ", i, " is missing markers: ", paste(missing, collapse = ", "))
    }
  }

  batch_train <- batch_labels[seq_along(train_frames)]
  batch_val   <- batch_labels[(length(train_frames) + 1):n_total]
  n_batches <- length(unique(batch_labels))

  cat("CytoNorm batch correction\n")
  cat("  Train frames:", length(train_frames), "\n")
  cat("  Validation frames:", length(validation_frames), "\n")
  cat("  Batches:", n_batches, "(", paste(unique(batch_labels), collapse = ", "), ")\n")
  cat("  Markers to correct:", length(markers), "\n")
  cat("  FlowSOM clusters:", n_clusters, "\n")
  cat("  Goal:", goal, "\n")

  set.seed(seed)

  # --- Try CytoNorm ----------------------------------------------------------
  tryCatch({
    require(CytoNorm)
    require(FlowSOM)

    cat("Training CytoNorm model on reference frames...\n")

    # Convert train frames to a flowSet with batch info in pData
    train_fset <- flowCore::flowSet(train_frames)
    flowCore::pData(train_fset)$batch <- batch_train
    flowCore::pData(train_fset)$sample <- paste0("train_", seq_along(train_frames))

    # Identify channel names corresponding to the requested markers
    channels <- .markers_to_channels(train_frames[[1]], markers)

    # Train the CytoNorm model
    model <- CytoNorm::CytoNorm.train(
      files    = train_fset,
      labels   = batch_train,
      channels = channels,
      transformList = NULL,
      FlowSOM.params = list(nCells = NULL, xdim = 5, ydim = 5,
                            nClus = n_clusters, seed = seed),
      goal = goal,
      verbose = TRUE
    )
    cat("CytoNorm model trained.\n")

    # Normalize validation frames
    cat("Normalizing validation frames...\n")
    val_fset <- flowCore::flowSet(validation_frames)
    flowCore::pData(val_fset)$batch <- batch_val
    flowCore::pData(val_fset)$sample <- paste0("val_", seq_along(validation_frames))

    corrected_fset <- CytoNorm::CytoNorm.normalize(
      model    = model,
      files    = val_fset,
      labels   = batch_val,
      channels = channels,
      transformList = NULL,
      verbose  = TRUE
    )

    # Convert flowSet back to list
    corrected_frames <- lapply(seq_len(length(corrected_fset)), function(i) {
      corrected_fset[[i]]
    })
    names(corrected_frames) <- names(validation_frames)

    cat("CytoNorm correction complete for", length(corrected_frames), "frames.\n")

    list(corrected_frames = corrected_frames, model = model)

  }, error = function(e) {
    cat("CytoNorm failed:", conditionMessage(e), "\n")
    cat("Falling back to quantile normalization...\n")
    .fallback_quantile_norm(train_frames, validation_frames,
                            batch_train, batch_val, markers, goal)
  })
}

# =============================================================================
# Internal helpers
# =============================================================================

#' Map marker names to channel names in a flowFrame
#' @param ff A flowFrame object
#' @param markers Character vector of marker names
#' @return Character vector of channel names
#' @keywords internal
.markers_to_channels <- function(ff, markers) {
  ff_data <- explore_data(ff)
  channels <- character(length(markers))
  for (i in seq_along(markers)) {
    idx <- which(ff_data$markers == markers[i])
    if (length(idx) == 0) stop("Marker '", markers[i], "' not found in flowFrame")
    channels[i] <- ff_data$dyes[idx[1]]
  }
  channels
}

#' Fallback quantile normalization when CytoNorm is unavailable
#'
#' Computes per-marker quantile reference from train_frames, then maps each
#' validation frame's marker distribution onto the reference quantiles.
#'
#' @param train_frames List of flowFrames (reference)
#' @param validation_frames List of flowFrames (to correct)
#' @param batch_train Batch labels for train frames
#' @param batch_val Batch labels for validation frames
#' @param markers Character vector of marker names to correct
#' @param goal Not used in fallback (always normalizes to pooled reference)
#' @return list(corrected_frames, model = NULL)
#' @keywords internal
.fallback_quantile_norm <- function(train_frames, validation_frames,
                                     batch_train, batch_val, markers, goal) {
  require(flowCore)

  n_quantiles <- 99
  probs <- seq(0, 1, length.out = n_quantiles + 2)

  # Build reference quantiles from train frames
  cat("  Computing reference quantiles from", length(train_frames), "training frames...\n")
  ref_quantiles <- list()
  for (m in markers) {
    all_vals <- c()
    for (ff in train_frames) {
      ff_data <- explore_data(ff)
      col_idx <- which(ff_data$markers == m)
      if (length(col_idx) > 0) {
        all_vals <- c(all_vals, flowCore::exprs(ff)[, col_idx[1]])
      }
    }
    ref_quantiles[[m]] <- quantile(all_vals, probs = probs, na.rm = TRUE)
  }

  # Apply quantile normalization to each validation frame
  corrected_frames <- list()
  for (i in seq_along(validation_frames)) {
    ff <- validation_frames[[i]]
    expr <- flowCore::exprs(ff)
    ff_data <- explore_data(ff)

    for (m in markers) {
      col_idx <- which(ff_data$markers == m)
      if (length(col_idx) == 0) next

      vals <- expr[, col_idx[1]]
      sample_q <- quantile(vals, probs = probs, na.rm = TRUE)

      # Map through quantile function: find quantile rank in sample, apply reference
      ranked <- ecdf(vals)(vals)
      corrected_vals <- approx(
        x = probs, y = ref_quantiles[[m]],
        xout = ranked, rule = 2
      )$y

      expr[, col_idx[1]] <- corrected_vals
    }

    flowCore::exprs(ff) <- expr
    corrected_frames[[i]] <- ff
    cat("  Corrected frame", i, "of", length(validation_frames),
        "(", nrow(expr), "cells)\n")
  }

  names(corrected_frames) <- names(validation_frames)
  cat("Quantile normalization fallback complete.\n")

  list(corrected_frames = corrected_frames, model = NULL)
}
