#!/usr/bin/env Rscript
# =============================================================================
# extract_flowjo_population_cli.R — pipeline F, step 1
# =============================================================================
# Turns a folder of FlowJo workspaces plus a folder of FCS files into the
# repo's ordinary expression-matrix RDS, so every step after this one is a step
# the pipeline already had.
#
# Thin, like every other wrapper here: it parses arguments and calls
# extract_flowjo_population() in src/visualization/umap_from_flowjo.R. Nothing
# about workspaces, gating, pairing or marker panels is decided in this file.
#
# Writes:
#   expr_data.rds     list(combined_matrix, sample_labels) — the same shape
#                     build_expr_matrix_cli.R emits, which is what lets
#                     DOWNSAMPLE and DIMRED consume it unchanged.
#   subpop_ids.rds    the FlowJo subgate of each cell. Handed to DIMRED as its
#                     `cluster_rds`: a gate IS a cluster as far as the plotting
#                     is concerned, and passing it there means the standard
#                     by-cluster figure comes out right with no new code.
#   flowjo_labels.rds workspace / annotation labels and the marker list, for
#                     the facet step.
#   load_report.json  the full load report (see umap_load_report.R).
# =============================================================================
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--wsp", type = "character",
              help = "A .wsp file, or a directory searched recursively"),
  make_option("--fcs_dir", type = "character",
              help = "Directory of FCS files, searched recursively"),
  make_option("--population", type = "character",
              default = "Lymphocytes/Single Cells/Viable cells"),
  make_option("--population_map", type = "character", default = NULL),
  make_option("--group", type = "character", default = "1"),
  make_option("--match_mode", type = "character", default = "content"),
  make_option("--allow_missing_samples", action = "store_true", default = FALSE),
  make_option("--marker_aliases", type = "character", default = NULL),
  make_option("--workspaces", type = "character", default = NULL),
  make_option("--no_recover_markers", action = "store_true", default = FALSE),
  make_option("--markers", type = "character", default = NULL),
  make_option("--subpop_mode", type = "character", default = "furthest"),
  make_option("--subpop_select", type = "character", default = NULL),
  make_option("--sample_metadata", type = "character", default = NULL),
  make_option("--sample_pattern", type = "character", default = NULL),
  make_option("--annotate", type = "character", default = NULL),
  make_option("--sample_id_column", type = "character", default = NULL),
  make_option("--normalize", type = "character", default = "none"),
  make_option("--downsample", type = "integer", default = 0),
  make_option("--seed", type = "integer", default = 42),
  make_option("--list_populations", action = "store_true", default = FALSE,
              help = "Report the load and exit without extracting anything"),
  make_option("--output_dir", type = "character", default = "flowjo_population"),
  make_option("--report_json", type = "character", default = NULL)
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

split_csv <- function(s) {
  if (is.null(s)) return(NULL)
  v <- trimws(strsplit(s, ",")[[1]])
  v <- v[nzchar(v)]
  if (length(v) == 0) NULL else v
}

make_outdir(opt$output_dir)
report_json <- opt$report_json %||% file.path(opt$output_dir, "load_report.json")

ext <- extract_flowjo_population(
  wsp_path = opt$wsp, fcs_dir = opt$fcs_dir, population = opt$population,
  group = opt$group, sample_metadata = opt$sample_metadata,
  sample_pattern = opt$sample_pattern, annotate = split_csv(opt$annotate),
  sample_id_column = opt$sample_id_column, markers = split_csv(opt$markers),
  subpop_mode = opt$subpop_mode, subpop_select = split_csv(opt$subpop_select),
  normalize = opt$normalize, downsample = opt$downsample, seed = opt$seed,
  list_populations = opt$list_populations,
  allow_missing_samples = opt$allow_missing_samples,
  population_map = opt$population_map, match_mode = opt$match_mode,
  marker_aliases = opt$marker_aliases, workspaces = opt$workspaces,
  recover_markers = !opt$no_recover_markers,
  report_json = report_json)

if (isTRUE(ext$inspected)) {
  cat("Inspect only — report written to ", report_json, "\n", sep = "")
  quit(status = 0)
}

# The pipeline contract, and nothing else in this file.
saveRDS(list(combined_matrix = ext$combined_matrix,
             sample_labels   = ext$sample_labels),
        file.path(opt$output_dir, "expr_data.rds"))
saveRDS(ext$subpop_labels_per_cell,
        file.path(opt$output_dir, "subpop_ids.rds"))
saveRDS(list(workspace_labels = ext$workspace_labels,
             subpop_levels    = ext$subpop_levels,
             subpop_mode_tag  = ext$subpop_mode_tag,
             file_labels      = ext$file_labels,
             n_cells_per      = ext$n_cells_per,
             samples          = ext$samples,
             sample_index     = ext$sample_index,
             sample_annotations = ext$sample_annotations,
             annot_cols       = ext$annot_cols,
             use_markers      = ext$use_markers,
             population       = ext$population,
             workspaces       = ext$workspaces),
        file.path(opt$output_dir, "flowjo_labels.rds"))

cat("Extracted ", nrow(ext$combined_matrix), " cells x ",
    ncol(ext$combined_matrix), " markers from ", length(ext$workspaces),
    " workspace(s) -> ", opt$output_dir, "\n", sep = "")
