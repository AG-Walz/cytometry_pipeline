#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_rds", type = "character", help = "Input flowFrame RDS"),
  make_option("--sample_name", type = "character", help = "Sample name"),
  make_option("--output_rds", type = "character", default = "preprocessed.rds"),
  make_option("--transform_type", type = "character", default = "arcsinh"),
  make_option("--qc_type", type = "character", default = "None"),
  make_option("--pp_order", type = "character", default = "QT"),
  make_option("--markers_of_interest", type = "character", default = NULL),
  make_option("--quant_norm", type = "logical", default = FALSE),
  make_option("--scaling", type = "logical", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Input validation ---
if (!is.null(opt$input_rds) && !file.exists(opt$input_rds)) {
  stop("Input file not found: ", opt$input_rds)
}

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(src_dir, "preprocessing", "quality_control.R"))
source(file.path(src_dir, "preprocessing", "transformation.R"))
source(file.path(src_dir, "preprocessing", "normalization.R"))

tryCatch({
  ff <- readRDS(opt$input_rds)
  markers <- if (!is.null(opt$markers_of_interest)) {
    strsplit(opt$markers_of_interest, ",")[[1]]
  } else {
    NULL
  }

  explored <- explore_data(ff)
  marker_indices <- explored$mcols

  if (opt$pp_order == "QT") {
    ff <- quality_control(ff, opt$sample_name,
                          qc_type = opt$qc_type,
                          marker_indices = marker_indices,
                          output_dir = NULL)
    ff <- transform_flowframe(ff, opt$sample_name, markers,
                              NULL, NULL, opt$transform_type)
  } else {
    ff <- transform_flowframe(ff, opt$sample_name, markers,
                              NULL, NULL, opt$transform_type)
    ff <- quality_control(ff, opt$sample_name,
                          qc_type = opt$qc_type,
                          marker_indices = marker_indices,
                          output_dir = NULL)
  }

  if (isTRUE(opt$quant_norm)) ff <- quantile_norm(ff, marker_indices)
  if (isTRUE(opt$scaling))    ff <- z_score_scaling(ff, marker_indices)

  saveRDS(ff, opt$output_rds)
  cat("Preprocessed:", opt$sample_name, "->", opt$output_rds, "\n")
}, error = function(e) {
  cat(sprintf("\n[ERROR] Preprocessing failed: %s\n", conditionMessage(e)), file = stderr())
  cat("Stack trace:\n", file = stderr())
  traceback(e)
  quit(status = 1)
})
