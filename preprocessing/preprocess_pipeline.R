#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/preprocess_pipeline.R — Preprocessing Orchestrator
# =============================================================================
# Simplified orchestrator that chains QC -> Transform -> Normalize
# Refactored from: R_scripts/02_preprocessing/Preprocess.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(dirname(getwd()), "preprocessing", "quality_control.R"))
source(file.path(dirname(getwd()), "preprocessing", "transformation.R"))
source(file.path(dirname(getwd()), "preprocessing", "normalization.R"))

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
