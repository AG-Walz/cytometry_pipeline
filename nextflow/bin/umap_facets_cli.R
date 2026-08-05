#!/usr/bin/env Rscript
# =============================================================================
# umap_facets_cli.R — pipeline F, last step
# =============================================================================
# Colours ONE embedding many ways: by source file, by workspace, by the FlowJo
# subgate each cell falls in, by every marker, and by every sample annotation.
#
# Reads the coordinates a previous DIMRED step computed rather than embedding
# again. That is the whole reason this is its own step: the embedding is the
# expensive part and is paid for once, while re-colouring is cheap and is what
# anyone actually iterates on.
#
# Thin: everything below the argument parsing is
# plot_flowjo_umap_facets() in src/visualization/umap_from_flowjo.R.
# =============================================================================
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--expr_rds", type = "character", help = "expr_data.rds"),
  make_option("--labels_rds", type = "character", help = "flowjo_labels.rds"),
  make_option("--subpop_rds", type = "character", help = "subpop_ids.rds"),
  make_option("--coords_rds", type = "character",
              help = "Coordinates from dimred_cli.R --output_rds"),
  make_option("--output_dir", type = "character", default = "umap_facets"),
  make_option("--normalize", type = "character", default = "none",
              help = "Only used to label the colour scales honestly"),
  make_option("--point_size", type = "double", default = 0.4),
  make_option("--point_alpha", type = "double", default = 0.6),
  make_option("--plot_width", type = "double", default = 10),
  make_option("--plot_height", type = "double", default = 8),
  make_option("--output_formats", type = "character", default = "png")
)
opt <- parse_args(OptionParser(option_list = option_list))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))
source(file.path(shared_dir, "stats_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(shared_dir, "dimred_dispatch.R"))
source(file.path(src_dir, "loading", "flowjo_workspace.R"))
source(file.path(src_dir, "loading", "sample_annotations.R"))
source(file.path(src_dir, "harmonization", "panel_homogenization.R"))
source(file.path(src_dir, "harmonization", "panel_comparison.R"))
source(file.path(src_dir, "visualization", "umap_load_report.R"))
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))

expr   <- readRDS(opt$expr_rds)
labels <- readRDS(opt$labels_rds)
subpop <- readRDS(opt$subpop_rds)
coords <- readRDS(opt$coords_rds)$coords

parts <- trimws(strsplit(opt$output_formats, ",")[[1]])
output_formats <- list(pdf = "pdf" %in% parts, png = "png" %in% parts,
                       svg = "svg" %in% parts)

# Row alignment is the one thing that can silently ruin every figure here, so
# it is asserted rather than assumed — the coordinates and the labels came from
# two different processes and only this check ties them together.
stopifnot(nrow(coords) == nrow(expr$combined_matrix),
          length(subpop) == nrow(expr$combined_matrix),
          length(labels$file_labels) == nrow(expr$combined_matrix))

umap_df <- data.frame(
  UMAP1 = coords[, 1], UMAP2 = coords[, 2],
  file = labels$file_labels, workspace = labels$workspace_labels,
  subpop = subpop, stringsAsFactors = FALSE)

# Per-SAMPLE annotations expanded to per-CELL through the same order and counts
# the matrix was built with, exactly as umap_from_flowjo() does it.
if (length(labels$annot_cols)) {
  annot_of <- match(labels$samples, labels$sample_index$display)
  stopifnot(!anyNA(annot_of))
  for (cl in labels$annot_cols) {
    umap_df[[paste0("annot_", cl)]] <-
      rep(labels$sample_annotations[[cl]][annot_of], labels$n_cells_per)
  }
}

make_outdir(opt$output_dir)
plot_flowjo_umap_facets(
  umap_df = umap_df, combined_matrix = expr$combined_matrix,
  use_markers = labels$use_markers, subpop_labels = labels$subpop_levels,
  annot_cols = labels$annot_cols, population = labels$population,
  output_dir = opt$output_dir, n_workspaces = length(labels$workspaces),
  subpop_mode_tag = labels$subpop_mode_tag %||% "",
  normalize = opt$normalize, point_size = opt$point_size,
  point_alpha = opt$point_alpha, plot_width = opt$plot_width,
  plot_height = opt$plot_height, output_formats = output_formats)

cat("FlowJo UMAP facets ->", opt$output_dir, "\n")
