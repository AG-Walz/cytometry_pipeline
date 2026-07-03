#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--input_rds_dir", type = "character",
              help = "Directory containing per-sample RDS files"),
  make_option("--method", type = "character", default = "harmony"),
  make_option("--batch_names", type = "character",
              help = "Comma-separated batch labels per sample"),
  make_option("--markers_of_interest", type = "character"),
  make_option("--batch_do_pca", type = "logical", default = FALSE),
  make_option("--combat_method", type = "character", default = "parametric"),
  make_option("--cytonorm_n_clusters", type = "integer", default = 10),
  make_option("--cytonorm_goal", type = "character", default = "mean"),
  make_option("--seed", type = "integer", default = 42),
  make_option("--output_dir", type = "character", default = "corrected")
)
opt <- parse_args(OptionParser(option_list = option_list))

suppressPackageStartupMessages(library(flowCore))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(shared_dir, "batch_correction_dispatch.R"))
source(file.path(src_dir, "loading", "create_flowframes.R"))

markers <- strsplit(opt$markers_of_interest, ",")[[1]]
batch_labels <- strsplit(opt$batch_names, ",")[[1]]

# Load all flowFrames
rds_files <- sort(list.files(opt$input_rds_dir, pattern = "\\.rds$",
                             full.names = TRUE))
flow_frames <- lapply(rds_files, readRDS)
sample_names <- tools::file_path_sans_ext(basename(rds_files))
names(flow_frames) <- sample_names

batch_result <- run_batch_correction_dispatch(
  data                = flow_frames,
  method              = opt$method,
  batch_labels        = batch_labels,
  markers_of_interest = markers,
  sample_names        = sample_names,
  params              = list(
    batch_do_pca        = opt$batch_do_pca,
    combat_method       = opt$combat_method,
    cytonorm_n_clusters = opt$cytonorm_n_clusters,
    cytonorm_goal       = opt$cytonorm_goal,
    seed                = opt$seed
  )
)

# Convert back to flowFrame list if needed
if (!is.null(batch_result$corrected_frames)) {
  corrected <- batch_result$corrected_frames
} else if (!is.null(batch_result$corrected_matrix)) {
  corrected <- create_flowframes_from_matrices(
    expr_list = split(as.data.frame(batch_result$corrected_matrix),
                      batch_result$sample_labels),
    markers_of_interest = batch_result$markers
  )
} else {
  corrected <- flow_frames
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
for (nm in names(corrected)) {
  saveRDS(corrected[[nm]], file.path(opt$output_dir, paste0(nm, ".rds")))
}
cat("Batch-corrected", length(corrected), "samples ->", opt$output_dir, "\n")
