#!/usr/bin/env Rscript
# =============================================================================
# harmonization/combat.R — ComBat Batch Correction
# =============================================================================
# Parametric empirical Bayes batch correction via sva::ComBat.
# Operates on expression matrices (cells x markers) and preserves biological
# covariates when provided.
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

#' Batch correction using ComBat
#'
#' Input from: build_combined_expr_matrix() (R_shared/flowframe_utils.R),
#'             combine_expression_matrices() (analysis/expression_matrices.R)
#' Output to:  run_leiden() (clustering/leiden.R),
#'             plot_umap_from_matrix() (visualization/umap.R)
#'
#' @param expr_matrix Numeric matrix (cells x markers)
#' @param batch_vector Character/factor vector of batch assignments (one per row)
#' @param covariates Optional data.frame of biological covariates to preserve (e.g., condition)
#' @param method "parametric" or "nonparametric"
#' @return Corrected numeric matrix (same dimensions)
run_combat <- function(expr_matrix, batch_vector, covariates = NULL,
                        method = "parametric") {
  require(sva)

  # --- Validate inputs -------------------------------------------------------
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("expr_matrix must be a numeric matrix or data.frame")
  }
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
  }
  if (length(batch_vector) != nrow(expr_matrix)) {
    stop("batch_vector length (", length(batch_vector),
         ") must equal nrow(expr_matrix) (", nrow(expr_matrix), ")")
  }

  batch_vector <- as.factor(batch_vector)
  batch_levels <- levels(batch_vector)
  if (length(batch_levels) < 2) {
    stop("At least 2 batches are required for ComBat correction, found: ",
         length(batch_levels))
  }

  par_prior <- switch(method,
    "parametric"    = TRUE,
    "nonparametric" = FALSE,
    stop("method must be 'parametric' or 'nonparametric', got: ", method)
  )

  # --- Log batch summary ------------------------------------------------------
  batch_table <- table(batch_vector)
  cat("ComBat batch correction (", method, ")\n", sep = "")
  cat("  Input:", nrow(expr_matrix), "cells x", ncol(expr_matrix), "markers\n")
  cat("  Batches (", length(batch_levels), "):\n", sep = "")
  for (b in names(batch_table)) {
    cat("    ", b, ":", batch_table[b], "cells\n")
  }

  # --- Build model matrix from covariates ------------------------------------
  mod <- NULL
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) {
      stop("covariates must be a data.frame")
    }
    if (nrow(covariates) != nrow(expr_matrix)) {
      stop("covariates nrow (", nrow(covariates),
           ") must equal nrow(expr_matrix) (", nrow(expr_matrix), ")")
    }
    # Check for NA values in covariates
    na_cols <- names(covariates)[sapply(covariates, function(col) any(is.na(col)))]
    if (length(na_cols) > 0) {
      stop("Covariates contain NA values in columns: ",
           paste(na_cols, collapse = ", "),
           ". Remove or impute NAs before running ComBat.")
    }
    cat("  Preserving biological covariates:", paste(names(covariates), collapse = ", "), "\n")
    mod_formula <- as.formula(
      paste("~", paste(names(covariates), collapse = " + "))
    )
    mod <- model.matrix(mod_formula, data = covariates)
  }

  # --- Run ComBat ------------------------------------------------------------
  # ComBat expects features (markers) in rows and samples (cells) in columns
  cat("Transposing matrix (markers x cells)...\n")
  expr_t <- t(expr_matrix)

  cat("Running ComBat...\n")
  tryCatch({
    corrected_t <- sva::ComBat(
      dat    = expr_t,
      batch  = batch_vector,
      mod    = mod,
      par.prior = par_prior,
      prior.plots = FALSE
    )

    # Transpose back to cells x markers
    corrected <- t(corrected_t)
    rownames(corrected) <- rownames(expr_matrix)
    colnames(corrected) <- colnames(expr_matrix)

    cat("ComBat correction complete.\n")
    cat("  Output:", nrow(corrected), "cells x", ncol(corrected), "markers\n")

    corrected

  }, error = function(e) {
    stop("ComBat failed: ", conditionMessage(e),
         "\n  Hint: Ensure each batch has enough cells and that markers have ",
         "non-zero variance in every batch.")
  })
}
