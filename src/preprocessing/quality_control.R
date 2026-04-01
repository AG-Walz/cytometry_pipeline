#!/usr/bin/env Rscript
# =============================================================================
# preprocessing/quality_control.R — Quality Control Functions
# =============================================================================
# Extracted from: R_scripts/02_preprocessing/Preprocessing_functions.R
# =============================================================================

#' Run quality control on a flowFrame
#'
#' @param flowframe A flowFrame object
#' @param sample_name Sample identifier string
#' @param qc_type Algorithm: "flowAI", "PeacoQC", or "None"
#' @param marker_indices Integer vector of marker column indices
#' @param output_dir Output directory for QC results
#' @return QC-filtered flowFrame
quality_control <- function(flowframe, sample_name,
                            qc_type = "PeacoQC",
                            marker_indices = NULL,
                            output_dir = NULL) {
  if (is.null(marker_indices)) {
    marker_indices <- seq_len(ncol(flowframe))
  }

  qc_result <- switch(qc_type,
    "flowAI" = {
      cat("Sample", sample_name, "QC with flowAI.\n")
      qc_ff <- flowAI::flow_auto_qc(
        flowframe,
        folder_results = if (!is.null(output_dir))
          file.path(output_dir, paste0("FlowAI_", sample_name))
        else tempdir()
      )
      if (!is.null(output_dir)) {
        flowCore::write.FCS(qc_ff,
          file.path(output_dir, paste0("FlowAI_", sample_name, ".fcs")))
      }
      cat(sprintf("After QC: %d/%d cells remained.\n", nrow(qc_ff), nrow(flowframe)))
      qc_ff
    },
    "PeacoQC" = {
      cat("Sample", sample_name, "QC with PeacoQC.\n")
      qc_result <- PeacoQC::PeacoQC(
        flowframe,
        channels = marker_indices,
        plot = TRUE,
        save_fcs = !is.null(output_dir),
        determine_good_cells = "all",
        output_directory = if (is.null(output_dir)) tempdir() else output_dir
      )
      cat(sprintf("After QC: %d/%d cells remained.\n",
                  nrow(qc_result$FinalFF), nrow(flowframe)))
      qc_result$FinalFF
    },
    "None" = {
      cat("Sample", sample_name, "— no QC applied.\n")
      flowframe
    },
    stop("Invalid qc_type: '", qc_type, "'. Use 'flowAI', 'PeacoQC', or 'None'.")
  )

  qc_result
}
