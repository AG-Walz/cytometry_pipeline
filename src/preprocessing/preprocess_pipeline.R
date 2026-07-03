#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/preprocess_pipeline.R — Preprocessing Orchestrator
# =============================================================================
# Simplified orchestrator that chains QC -> Transform -> Normalize
# Refactored from: R_scripts/02_preprocessing/Preprocess.R
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

shared_dir <- cyto_shared_dir()
tryCatch(source(file.path(shared_dir, "io_utils.R")),
         error = function(e) stop("Failed to load io_utils: ", e$message))
tryCatch(source(file.path(shared_dir, "flowframe_utils.R")),
         error = function(e) stop("Failed to load flowframe_utils: ", e$message))
tryCatch(source(file.path(cyto_src_root(), "preprocessing", "quality_control.R")),
         error = function(e) stop("Failed to load quality_control module: ", e$message))
tryCatch(source(file.path(cyto_src_root(), "preprocessing", "transformation.R")),
         error = function(e) stop("Failed to load transformation module: ", e$message))
tryCatch(source(file.path(cyto_src_root(), "preprocessing", "normalization.R")),
         error = function(e) stop("Failed to load normalization module: ", e$message))

#' Run the full preprocessing pipeline on a list of flowFrames
#'
#' @param flow_frames Named list of flowFrame objects
#' @param params List of preprocessing parameters (from config/Parameters_template.R)
#' @param output_dir Output directory for preprocessed files
#' @return Named list of preprocessed flowFrame objects
run_preprocessing <- function(flow_frames, params, output_dir = NULL) {
  if (!is.null(output_dir)) {
    pp_dir <- file.path(output_dir, "preprocessed_files")
    qc_dir <- file.path(output_dir, "qc_results")
    make_outdir(pp_dir)
    make_outdir(qc_dir)
  } else {
    pp_dir <- qc_dir <- NULL
  }

  sample_names <- names(flow_frames)
  pp_frames <- list()

  for (i in seq_along(flow_frames)) {
    sname <- sample_names[i]
    cat("\n=== Preprocessing sample:", sname, "===\n")

    ff <- flow_frames[[i]]
    marker_indices <- explore_data(ff)$mcols

    # Step 1: QC or Transform (depending on pp_order)
    if (params$pp_order == "QT") {
      ff <- quality_control(ff, sname, qc_type = params$qc_type,
                           marker_indices = marker_indices, output_dir = qc_dir)
      ff <- transform_flowframe(ff, sname, params$markers_of_interest,
                                params$offset_list, params$width_list,
                                params$transform_type)
    } else {
      ff <- transform_flowframe(ff, sname, params$markers_of_interest,
                                params$offset_list, params$width_list,
                                params$transform_type)
      ff <- quality_control(ff, sname, qc_type = params$qc_type,
                           marker_indices = marker_indices, output_dir = qc_dir)
    }

    # Step 2: Quantile normalization (optional)
    if (isTRUE(params$quant_norm)) {
      cat("Quantile normalizing", sname, "\n")
      ff <- quantile_norm(ff, marker_indices)
    }

    # Step 3: Scaling (optional)
    if (isTRUE(params$scaling)) {
      cat("Scaling", sname, "\n")
      ff <- z_score_scaling(ff, marker_indices)
    }

    # Save outputs
    if (!is.null(pp_dir)) {
      if (isTRUE(params$save_pp_fcs)) {
        flowCore::write.FCS(ff, file.path(pp_dir, paste0("pp_", sname, ".fcs")))
      }
      if (isTRUE(params$save_pp_csv)) {
        write.csv(flowCore::exprs(ff), file.path(pp_dir, paste0("pp_", sname, ".csv")))
      }
    }

    pp_frames[[sname]] <- ff
  }

  cat("\nPreprocessing complete for", length(pp_frames), "samples.\n")
  pp_frames
}
