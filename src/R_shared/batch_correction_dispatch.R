#!/usr/bin/env Rscript
# =============================================================================
# R_shared/batch_correction_dispatch.R — Batch Correction Dispatcher
# =============================================================================
# Routes to Harmony, CytoNorm, ComBat, or pass-through via a single config
# parameter. Bridges input types: accepts flowFrame list OR matrix, converts
# as needed for each method.
#
# Usage:
#   source(file.path(shared_dir, "batch_correction_dispatch.R"))
#   result <- run_batch_correction_dispatch(
#     data, method = "combat", batch_labels = c("A","A","B","B"),
#     markers_of_interest = markers, sample_names = samples)
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
      d <- if (!is.null(of)) dirname(normalizePath(of)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}

# %||% polyfill for R < 4.0
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' Dispatch batch correction to the specified method
#'
#' Accepts either a list of flowFrame objects or a numeric matrix. Converts
#' between formats as needed for each backend method.
#'
#' @param data Either a named list of flowFrame objects OR a numeric matrix
#'   (cells x markers). If matrix, sample_labels must also be provided.
#' @param method One of "harmony", "cytonorm", "combat", "none"
#' @param batch_labels Character vector of batch assignments. If data is
#'   flowFrames: one per frame. If data is matrix: one per row.
#' @param markers_of_interest Character vector of marker names
#' @param sample_names Character vector of sample identifiers. For flowFrames:
#'   one per frame. For matrix: one per row (cell-level sample labels).
#' @param params Named list of method-specific parameters
#' @return List with:
#'   \item{corrected_matrix}{Numeric matrix (cells x markers)}
#'   \item{sample_labels}{Character vector (length == nrow)}
#'   \item{batch_labels}{Character vector (length == nrow)}
#'   \item{markers}{Character vector of marker names}
#'   \item{corrected_frames}{List of corrected flowFrames (CytoNorm only; NULL otherwise)}
run_batch_correction_dispatch <- function(data, method = "harmony",
                                           batch_labels,
                                           markers_of_interest = NULL,
                                           sample_names = NULL,
                                           params = list()) {
  method <- tolower(method)
  valid_methods <- c("harmony", "cytonorm", "combat", "none")
  if (!method %in% valid_methods) {
    stop("Unknown batch correction method '", method, "'. Valid: ",
         paste(valid_methods, collapse = ", "))
  }

  cat("Dispatching batch correction: method =", method, "\n")

  # Detect input type
  is_flowframes <- is.list(data) && !is.matrix(data) &&
    length(data) > 0 && inherits(data[[1]], "flowFrame")

  # Resolve source directories
  if (exists("shared_dir")) {
    harm_dir <- file.path(dirname(shared_dir), "harmonization")
    shared <- shared_dir
  } else {
    harm_dir <- file.path(cyto_src_root(), "harmonization")
    shared <- cyto_shared_dir()
  }

  # --- "none" pass-through ---------------------------------------------------
  if (method == "none") {
    cat("Batch correction: none (pass-through)\n")
    if (is_flowframes) {
      source(file.path(shared, "flowframe_utils.R"))
      expr_data <- build_combined_expr_matrix(
        data, sample_names %||% names(data), markers_of_interest
      )
      cell_batch <- rep(batch_labels, sapply(data, function(ff) nrow(flowCore::exprs(ff))))
      return(list(
        corrected_matrix = expr_data$combined_matrix,
        sample_labels    = expr_data$sample_labels,
        batch_labels     = cell_batch,
        markers          = colnames(expr_data$combined_matrix),
        corrected_frames = NULL
      ))
    } else {
      return(list(
        corrected_matrix = data,
        sample_labels    = sample_names,
        batch_labels     = batch_labels,
        markers          = colnames(data),
        corrected_frames = NULL
      ))
    }
  }

  # --- Harmony (needs flowFrame list) ----------------------------------------
  if (method == "harmony") {
    tryCatch(source(file.path(harm_dir, "batch_correction.R")),
             error = function(e) stop("Failed to load Harmony module: ", e$message,
                                      "\nExpected at: ", file.path(harm_dir, "batch_correction.R")))

    if (is_flowframes) {
      # Direct: Harmony expects flowFrame list
      result <- run_batch_correction(
        flow_frames         = data,
        sample_names        = sample_names %||% names(data),
        batch_labels        = batch_labels,
        markers_of_interest = markers_of_interest,
        do_pca              = params$batch_do_pca %||% FALSE
      )
    } else {
      # Convert matrix to flowFrames first
      source(file.path(shared, "flowframe_utils.R"))
      unique_samples <- unique(sample_names)
      ff_list <- list()
      frame_batch_labels <- character(length(unique_samples))
      for (i in seq_along(unique_samples)) {
        s <- unique_samples[i]
        idx <- which(sample_names == s)
        mtx <- data[idx, , drop = FALSE]
        ff_list[[s]] <- from_mtx_to_ffr(mtx)
        # Find batch for this sample (first matching cell's batch)
        frame_batch_labels[i] <- batch_labels[idx[1]]
      }
      result <- run_batch_correction(
        flow_frames         = ff_list,
        sample_names        = unique_samples,
        batch_labels        = frame_batch_labels,
        markers_of_interest = markers_of_interest,
        do_pca              = params$batch_do_pca %||% FALSE
      )
    }
    return(list(
      corrected_matrix = result$corrected_matrix,
      sample_labels    = result$sample_labels,
      batch_labels     = result$batch_labels,
      markers          = result$markers,
      corrected_frames = NULL
    ))
  }

  # --- ComBat (needs matrix) -------------------------------------------------
  if (method == "combat") {
    tryCatch(source(file.path(harm_dir, "combat.R")),
             error = function(e) stop("Failed to load ComBat module: ", e$message,
                                      "\nExpected at: ", file.path(harm_dir, "combat.R")))

    if (is_flowframes) {
      # Convert flowFrames to matrix first
      source(file.path(shared, "flowframe_utils.R"))
      expr_data <- build_combined_expr_matrix(
        data, sample_names %||% names(data), markers_of_interest
      )
      mtx <- expr_data$combined_matrix
      cell_samples <- expr_data$sample_labels
      # Expand frame-level batch_labels to cell-level
      cell_batch <- rep(batch_labels, sapply(data, function(ff) nrow(flowCore::exprs(ff))))
    } else {
      mtx <- data
      cell_samples <- sample_names
      cell_batch <- batch_labels
    }

    corrected <- run_combat(
      expr_matrix  = mtx,
      batch_vector = cell_batch,
      covariates   = params$combat_covariates %||% NULL,
      method       = params$combat_method %||% "parametric"
    )
    return(list(
      corrected_matrix = corrected,
      sample_labels    = cell_samples,
      batch_labels     = cell_batch,
      markers          = colnames(corrected),
      corrected_frames = NULL
    ))
  }

  # --- CytoNorm (needs flowFrame lists split into train/validation) ----------
  if (method == "cytonorm") {
    tryCatch(source(file.path(harm_dir, "cytonorm.R")),
             error = function(e) stop("Failed to load CytoNorm module: ", e$message,
                                      "\nExpected at: ", file.path(harm_dir, "cytonorm.R")))

    if (!is_flowframes) {
      stop("CytoNorm requires flowFrame input. Provide a list of flowFrame objects.")
    }

    if (is.null(markers_of_interest)) {
      stop("markers_of_interest is required for CytoNorm method. Specify in config.")
    }

    train_indices <- params$cytonorm_train_indices %||% seq_along(data)
    val_indices <- setdiff(seq_along(data), train_indices)
    if (length(val_indices) == 0) {
      warning("No validation indices specified; using all frames for both training ",
              "and normalization. Consider specifying cytonorm_train_indices in config.")
      val_indices <- seq_along(data)
    }

    cyto_result <- run_cytonorm(
      train_frames      = data[train_indices],
      validation_frames = data[val_indices],
      batch_labels      = batch_labels,
      markers           = markers_of_interest,
      n_clusters        = params$cytonorm_n_clusters %||% 10,
      goal              = params$cytonorm_goal %||% "mean",
      seed              = params$seed %||% 42
    )

    # Build combined matrix from corrected frames
    source(file.path(shared, "flowframe_utils.R"))
    corrected_names <- (sample_names %||% names(data))[val_indices]
    expr_data <- build_combined_expr_matrix(
      cyto_result$corrected_frames, corrected_names, markers_of_interest
    )
    cell_batch <- rep(
      batch_labels[val_indices],
      sapply(cyto_result$corrected_frames, function(ff) nrow(flowCore::exprs(ff)))
    )

    return(list(
      corrected_matrix = expr_data$combined_matrix,
      sample_labels    = expr_data$sample_labels,
      batch_labels     = cell_batch,
      markers          = colnames(expr_data$combined_matrix),
      corrected_frames = cyto_result$corrected_frames
    ))
  }
}
