#!/usr/bin/env Rscript
# =============================================================================
# scripts/cytokine_combinations.R — CLI wrapper for per-patient cytokine
# combination bar plots from FlowJo workspaces
# =============================================================================
# Thin wrapper: parses options, locates the repo 'src' directory, sources the
# shared utilities and the analysis module
# (src/analysis/cytokine_combinations.R) and calls cytokine_combinations().
#
# Example (one patient):
#   Rscript scripts/cytokine_combinations.R \
#     --wsp /mnt/volume/data/MTB_FACS/wsp/MTB059_1.wsp \
#     --fcs_dir /mnt/volume/data/MTB_FACS/fcs \
#     --output_dir /mnt/volume/data/MTB_FACS
#
# Every workspace in a folder, but plots for two patients only:
#   Rscript scripts/cytokine_combinations.R --wsp data/MTB_FACS/wsp \
#     --fcs_dir data/MTB_FACS/fcs --output_dir results/combos \
#     --patients MTB059,MTB074
#
# Re-plot from the CSV a previous run wrote (no CytoML needed):
#   Rscript scripts/cytokine_combinations.R --from_table results/combos/cytokine_combinations.csv \
#     --output_dir results/combos
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--wsp", type = "character", default = NULL,
              help = "A FlowJo .wsp file, or a folder searched recursively for them"),
  make_option("--fcs_dir", type = "character", default = NULL,
              help = "Folder searched recursively for the FCS files the workspaces name"),
  make_option("--output_dir", type = "character", default = "results/cytokine_combinations",
              help = "Where the PNGs and cytokine_combinations.csv go [default %default]"),
  make_option("--patients", type = "character", default = NULL,
              help = "Comma-separated patients to plot, e.g. MTB059,MTB074 [default: all]"),
  make_option("--parents", type = "character", default = "CD4,CD8",
              help = "Comma-separated parent gates (one facet row each) [default %default]"),
  make_option("--markers", type = "character", default = "CD107a,IFNy,IL-2,TNF",
              help = "Comma-separated cytokine gates under each parent [default %default]"),
  make_option("--group", type = "character", default = "All Samples",
              help = "FlowJo sample group to load [default %default]"),
  make_option("--output_formats", type = "character", default = "png",
              help = "Comma-separated: png, pdf, svg [default %default]"),
  make_option("--from_table", type = "character", default = NULL,
              help = "Re-plot from a cytokine_combinations.csv instead of loading workspaces"),
  make_option("--src_dir", type = "character", default = NULL,
              help = "Path to the repo 'src' directory [default: auto-detected]")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = paste0(
    "\nOne stacked bar per FCS file: the Boolean combinations of the cytokine ",
    "gates under CD4 and CD8 as % of parent, grouped by timepoint and faceted ",
    "by cocktail. One figure per patient, plus a neg-subtracted one.\n")
))

if (is.null(opt$from_table) && (is.null(opt$wsp) || is.null(opt$fcs_dir))) {
  stop("--wsp and --fcs_dir are required unless --from_table is given.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Locate the repo 'src' dir, set the shared_dir/src_dir globals (CLI convention),
# and source the shared utilities + analysis module in dependency order.
# -----------------------------------------------------------------------------
# find_script_path, resolve_src_dir, split_csv, parse_formats: shared by the CLI
# wrappers and kept next to this script, so located from Rscript's --file=.
script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file) == 1) {
  cli_dir <- dirname(sub("^--file=", "", script_file))
} else {
  cli_dir <- file.path(getwd(), "scripts")
}
source(file.path(cli_dir, "_cli_common.R"))

src_dir    <<- resolve_src_dir(opt$src_dir)
shared_dir <<- file.path(src_dir, "R_shared")
if (!dir.exists(shared_dir)) stop("R_shared directory not found at: ", shared_dir)

source(file.path(shared_dir, "io_utils.R"))        # make_outdir, emit_progress
source(file.path(shared_dir, "plot_utils.R"))      # save_plot_formats
source(file.path(shared_dir, "flowframe_utils.R")) # needed by umap_from_flowjo.R
# THE FlowJo workspace layer: discovery, pairing, GatingSet building
source(file.path(src_dir, "loading", "flowjo_workspace.R"))
# annotations_from_pattern (name -> patient/timepoint/cocktail)
source(file.path(src_dir, "loading", "sample_annotations.R"))
# clean_sample_label, get_subpop_membership (pure at source time)
source(file.path(src_dir, "harmonization", "panel_homogenization.R"))
source(file.path(src_dir, "harmonization", "panel_comparison.R"))
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))
source(file.path(src_dir, "analysis", "cytokine_combinations.R"))

cytokine_combinations(
  wsp            = opt$wsp,
  fcs_dir        = opt$fcs_dir,
  output_dir     = opt$output_dir,
  patients       = split_csv(opt$patients),
  parents        = split_csv(opt$parents),
  markers        = split_csv(opt$markers),
  group          = opt$group,
  output_formats = parse_formats(opt$output_formats),
  from_table     = opt$from_table
)
