#!/usr/bin/env Rscript
# =============================================================================
# scripts/umap_from_flowjo.R — CLI wrapper for the FlowJo-population UMAP
# =============================================================================
# Thin standalone wrapper (no Nextflow). Parses CLI options, locates the repo
# 'src' directory, sources the R_shared utilities and the analysis module
# (src/visualization/umap_from_flowjo.R), and calls umap_from_flowjo(). All of
# the analysis/plotting logic lives in that shared module — keep this file thin.
#
# Example:
#   Rscript scripts/umap_from_flowjo.R \
#     --wsp data/Killing_AML.wsp --fcs_dir data \
#     --population "Lymphocytes/Single Cells/Viable cells" \
#     --output_dir results/umap_viable
#
#   # Many workspaces in one shared embedding: --wsp is a FOLDER, searched
#   # recursively, and each workspace is paired with its FCS files by name.
#   Rscript scripts/umap_from_flowjo.R \
#     --wsp data/workspaces --fcs_dir data/fcs \
#     --group "All Samples" \
#     --population "Single Cells/Viable cells" \
#     --output_dir results/umap_cohort
#
#   # List the gating tree, then exit (so you can pick a population):
#   Rscript scripts/umap_from_flowjo.R --wsp data/Killing_AML.wsp \
#     --fcs_dir data --list_populations
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# -----------------------------------------------------------------------------
# CLI options
#
# Defaults for --wsp / --fcs_dir / --population point at the known on-disk data
# (data/Killing_AML.wsp, data/, and the target "Viable cells" population), so a
# bare `Rscript scripts/umap_from_flowjo.R` targets the intended analysis.
# -----------------------------------------------------------------------------
option_list <- list(
  make_option("--wsp", type = "character", default = "data/Killing_AML.wsp",
              help = paste0("A FlowJo .wsp file, OR a DIRECTORY searched ",
                            "recursively for .wsp files. With several ",
                            "workspaces all their cells go into ONE shared ",
                            "UMAP, and an extra 'by workspace' plot is ",
                            "written. [default %default]")),
  make_option("--fcs_dir", type = "character", default = "data",
              help = paste0("Directory holding the FCS files, searched ",
                            "RECURSIVELY. Each workspace is paired with its ",
                            "own files by file name, so subfolders need no ",
                            "particular layout. [default %default]")),
  make_option("--population", type = "character",
              default = "Lymphocytes/Single Cells/Viable cells",
              help = paste0("Population path to extract. A leading '/' is ",
                            "optional; a suffix (e.g. 'Single Cells/Viable cells') ",
                            "is matched against the gating tree. [default %default]")),
  make_option("--list_populations", action = "store_true", default = FALSE,
              help = "Print the population tree (gs_get_pop_paths) and exit."),
  make_option("--group", type = "character", default = "1",
              help = paste0("FlowJo sample group to load: an index (e.g. 1) ",
                            "or a group NAME (e.g. \"All Samples\"). A name ",
                            "is what lets one value cover several workspaces ",
                            "that order their groups differently. ",
                            "--list_populations prints each workspace's ",
                            "groups. [default %default]")),
  make_option("--expected_samples", type = "integer", default = 0,
              help = paste0("If > 0, warn when the run does not total exactly ",
                            "this many samples across all workspaces ",
                            "[default %default]")),
  make_option("--sample_metadata", type = "character", default = NULL,
              help = paste0("Optional CSV keyed on sample (a 'sample_id' ",
                            "column, or the first column). Every other column ",
                            "becomes a per-sample annotation — patient, ",
                            "timepoint, response, anything — and each one ",
                            "yields a plot coloured by it plus one plot per ",
                            "level, all on the SAME embedding.")),
  make_option("--sample_pattern", type = "character", default = NULL,
              help = paste0("Alternative to --sample_metadata: a regex with ",
                            "NAMED capture groups applied to the sample name, ",
                            "e.g. \"(?<patient>MTB[0-9]+)_(?<timepoint>TP",
                            "[0-9]+)\". Group names become the columns. ",
                            "Ignored when --sample_metadata is given.")),
  make_option("--annotate", type = "character", default = NULL,
              help = paste0("Comma-separated subset of annotation columns to ",
                            "plot. Default: all of them.")),
  make_option("--sample_id_column", type = "character", default = NULL,
              help = paste0("Name of the --sample_metadata column holding the ",
                            "sample identifier [default: 'sample_id', else ",
                            "the first column].")),
  make_option("--markers", type = "character", default = NULL,
              help = paste0("Optional comma-separated subset of marker names to ",
                            "use/plot. Default: all common non-scatter markers.")),
  make_option("--subpop_mode", type = "character", default = "furthest",
              help = paste0("Colour the shared embedding by the FlowJo subgates ",
                            "found UNDER --population. One of: 'furthest' (label ",
                            "each cell by the most specific / furthest gating-tree ",
                            "node it reaches, falling back to a parent gate when ",
                            "no deeper node is reached; lineage CD4+/CD8+ falls ",
                            "out of the node path), 'cd4cd8' (only the CD4+/CD8+ ",
                            "lineage split), 'direct' (every immediate child gate), ",
                            "'deepest' (deepest by path depth), 'none' (skip). ",
                            "Note: 'direct'/'deepest' collapse to '(multiple)' ",
                            "when sibling gates overlap. Cells in no selected ",
                            "subgate are labelled '(none)'. [default %default]")),
  make_option("--subpop_select", type = "character", default = NULL,
              help = paste0("Comma-separated list of gates (relative path under ",
                            "--population, or bare gate name) to use as the ",
                            "colour partition, overriding --subpop_mode. Use this ",
                            "for an explicit mutually-exclusive set, e.g. the 16 ",
                            "Boolean functional combos. A cell in >1 listed gate ",
                            "is '(multiple)'; in none is '(none)'.")),
  make_option("--output_dir", type = "character", default = "umap_flowjo",
              help = "Output directory for plots [default %default]"),
  make_option("--output_formats", type = "character", default = "svg",
              help = "Comma-separated formats: png,pdf,svg [default %default]"),
  make_option("--downsample", type = "integer", default = 0,
              help = paste0("Max cells per sample AFTER gating (0 = no ",
                            "downsampling) [default %default]")),
  make_option("--normalize", type = "character", default = "none",
              help = paste0("Per-marker normalization computed FROM the ",
                            "concatenated, gated cells themselves. One of: ",
                            "'none', 'zscore' (centre to mean 0, scale to unit ",
                            "SD), 'minmax' (rescale to [0,1]). Applied before ",
                            "UMAP and used for the marker colour scales. ",
                            "[default %default]")),
  make_option("--n_neighbors", type = "integer", default = 15,
              help = "UMAP n_neighbors [default %default]"),
  make_option("--min_dist", type = "double", default = 0.1,
              help = "UMAP min_dist [default %default]"),
  make_option("--metric", type = "character", default = "euclidean",
              help = "UMAP distance metric [default %default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed (UMAP + downsampling) [default %default]"),
  make_option("--point_size", type = "double", default = 0.4,
              help = "ggplot point size [default %default]"),
  make_option("--point_alpha", type = "double", default = 0.6,
              help = "ggplot point alpha [default %default]"),
  make_option("--plot_width", type = "double", default = 10,
              help = "Plot width (inches) [default %default]"),
  make_option("--plot_height", type = "double", default = 8,
              help = "Plot height (inches) [default %default]"),
  make_option("--src_dir", type = "character", default = NULL,
              help = paste0("Path to the repo 'src' directory. Default: ",
                            "auto-detected relative to this script."))
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = paste0(
    "\nCompute a single shared UMAP over one FlowJo-gated population across ",
    "samples and plot it by sample, by subpopulation, and by every marker.\n")
))

# -----------------------------------------------------------------------------
# Locate the repo 'src' dir, set the shared_dir/src_dir globals (CLI convention),
# and source the shared utilities + analysis module in dependency order.
# -----------------------------------------------------------------------------
find_script_path <- function() {
  # Works under Rscript (--file=) and source().
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE))
  }
  if (!is.null(sys.function(1)) && !is.null(attr(body(sys.function(1)), "srcfile"))) {
    sf <- attr(body(sys.function(1)), "srcfile")
    if (!is.null(sf$filename)) return(normalizePath(sf$filename, mustWork = FALSE))
  }
  NA_character_
}

resolve_src_dir <- function(user_src) {
  if (!is.null(user_src)) {
    sd <- normalizePath(user_src, mustWork = FALSE)
    if (!dir.exists(sd)) stop("--src_dir does not exist: ", sd)
    return(sd)
  }
  sp <- find_script_path()
  if (!is.na(sp)) {
    # script is repo/scripts/umap_from_flowjo.R -> repo/src
    cand <- file.path(dirname(dirname(sp)), "src")
    if (dir.exists(cand)) return(normalizePath(cand))
  }
  # Fallback: relative to current working directory (repo root).
  cand <- file.path(getwd(), "src")
  if (dir.exists(cand)) return(normalizePath(cand))
  stop("Could not locate the 'src' directory. Pass --src_dir explicitly.")
}

src_dir    <<- resolve_src_dir(opt$src_dir)
shared_dir <<- file.path(src_dir, "R_shared")
if (!dir.exists(shared_dir)) {
  stop("R_shared directory not found at: ", shared_dir)
}

source(file.path(shared_dir, "io_utils.R"))        # make_outdir
source(file.path(shared_dir, "plot_utils.R"))      # save_plot_formats, create_safe_filename
source(file.path(shared_dir, "stats_utils.R"))     # run_umap
source(file.path(shared_dir, "flowframe_utils.R")) # explore_data
# THE FlowJo workspace layer: wsp_open, wsp_to_gatingset, discovery, pairing
source(file.path(src_dir, "loading", "flowjo_workspace.R"))
# Per-sample annotations (patient, timepoint, ...) from a CSV or a name regex
source(file.path(src_dir, "loading", "sample_annotations.R"))
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))  # umap_from_flowjo

# -----------------------------------------------------------------------------
# Parse the comma-separated CLI strings into the structured values the analysis
# module expects, then hand off. All analysis/plotting logic lives in the module.
# -----------------------------------------------------------------------------
split_csv <- function(s) {
  if (is.null(s)) return(NULL)
  v <- trimws(strsplit(s, ",")[[1]])
  v <- v[nzchar(v)]
  if (length(v) == 0) NULL else v
}

parse_formats <- function(s) {
  parts <- trimws(strsplit(s, ",")[[1]])
  list(pdf = "pdf" %in% parts, png = "png" %in% parts, svg = "svg" %in% parts)
}

umap_from_flowjo(
  wsp_path         = opt$wsp,
  fcs_dir          = opt$fcs_dir,
  population       = opt$population,
  output_dir       = opt$output_dir,
  group            = opt$group,
  expected_samples = opt$expected_samples,
  sample_metadata  = opt$sample_metadata,
  sample_pattern   = opt$sample_pattern,
  annotate         = split_csv(opt$annotate),
  sample_id_column = opt$sample_id_column,
  markers          = split_csv(opt$markers),
  subpop_mode      = opt$subpop_mode,
  subpop_select    = split_csv(opt$subpop_select),
  normalize        = opt$normalize,
  downsample       = opt$downsample,
  n_neighbors      = opt$n_neighbors,
  min_dist         = opt$min_dist,
  metric           = opt$metric,
  seed             = opt$seed,
  point_size       = opt$point_size,
  point_alpha      = opt$point_alpha,
  plot_width       = opt$plot_width,
  plot_height      = opt$plot_height,
  output_formats   = parse_formats(opt$output_formats),
  list_populations = opt$list_populations
)
