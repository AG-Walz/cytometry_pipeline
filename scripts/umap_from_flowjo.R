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
  make_option("--population_map", type = "character", default = NULL,
              help = paste0("Choose the gate PER WORKSPACE, for studies whose ",
                            "workspaces name the same biology differently. ",
                            "Either a path to a two-column file (workspace, ",
                            "gate) or inline \"wsA=/gate;wsB=/other gate\". ",
                            "A workspace with no entry is EXCLUDED from the ",
                            "run. Overrides --population. Use ",
                            "--list_populations to see each workspace's own ",
                            "gates.")),
  make_option("--list_populations", action = "store_true", default = FALSE,
              help = "Print the population tree (gs_get_pop_paths) and exit."),
  make_option("--annotation_preview", type = "character", default = NULL,
              help = paste0("Resolve --sample_metadata / --sample_pattern ",
                            "against a list of sample names and write the ",
                            "result to this path as JSON, then exit. Reads no ",
                            "FCS file and opens no workspace, so it answers in ",
                            "milliseconds — it exists so a GUI can show which ",
                            "columns a metadata table offers and how many ",
                            "samples it matches WITHOUT re-running a load, and ",
                            "without a second implementation of the matching ",
                            "rules that would drift from this one. Sample ",
                            "names come from --preview_samples.")),
  make_option("--preview_samples", type = "character", default = NULL,
              help = paste0("Newline-separated file of sample display names, ",
                            "for --annotation_preview.")),
  make_option("--report_json", type = "character", default = NULL,
              help = paste0("With --list_populations, also write the whole ",
                            "load to this path as JSON: every workspace with ",
                            "its status, sample counts, how each sample was ",
                            "paired, its gating tree and its marker panel, ",
                            "plus what diverges between them. This is what the ",
                            "Dash app reads instead of parsing the text above.")),
  make_option("--group", type = "character", default = "1",
              help = paste0("FlowJo sample group to load: an index (e.g. 1) ",
                            "or a group NAME (e.g. \"All Samples\"). A name ",
                            "is what lets one value cover several workspaces ",
                            "that order their groups differently. ",
                            "--list_populations prints each workspace's ",
                            "groups. [default %default]")),
  make_option("--allow_missing_samples", action = "store_true", default = FALSE,
              help = paste0("Drop samples whose FCS file is not found under ",
                            "--fcs_dir instead of failing the run, and report ",
                            "exactly what was dropped. A workspace with no ",
                            "file left is skipped. An AMBIGUOUS file name is ",
                            "still fatal: the sample is not absent, it is ",
                            "unidentifiable, and guessing would put the wrong ",
                            "cells in the embedding. [default %default]")),
  make_option("--match_mode", type = "character", default = "content",
              help = paste0("How a workspace's samples are paired with FCS ",
                            "files. 'content' (default): try the file NAME ",
                            "first, then let CytoML pair the rest by the $FIL ",
                            "keyword stored inside each file plus its event ",
                            "count — which survives a rename on disk, the ",
                            "usual reason a real cohort will not load. ",
                            "'strict': a sample whose file name does not match ",
                            "is missing, full stop (the behaviour before this ",
                            "option existed). [default %default]")),
  make_option("--workspaces", type = "character", default = NULL,
              help = paste0("Comma-separated workspace names to USE; the rest ",
                            "are left out. Default: every workspace found. ",
                            "Applied before any GatingSet is built, so leaving ",
                            "one out makes the run shorter as well as smaller. ",
                            "Names are the ones --list_populations reports.")),
  make_option("--n_cores", type = "integer", default = 0,
              help = paste0("Worker processes for the per-file plots, which ",
                            "are the bulk of a run (7 markers x N samples). ",
                            "0 = decide automatically, leaving one core free; ",
                            "1 = no parallelism. Capped at the cgroup CPU ",
                            "quota, so asking for more than a container allows ",
                            "does nothing. [default %default]")),
  make_option("--umap_init", type = "character", default = "auto",
              help = paste0("How the embedding is initialised: 'auto' ",
                            "(spectral below 300k cells, PCA above), ",
                            "'spectral', 'pca', 'random'. Spectral is the best ",
                            "starting layout and is an eigendecomposition of a ",
                            "graph with one node per cell — SINGLE-THREADED, ",
                            "and on millions of cells it is the whole run. ",
                            "[default %default]")),
  make_option("--umap_subsample", type = "integer", default = 0,
              help = paste0("Fit the embedding on this many randomly chosen ",
                            "cells, then project the rest onto it. 0 = use ",
                            "every cell. Every stage of fitting is ",
                            "superlinear in the cell count while projection ",
                            "is not, so on a multi-million-cell run this is ",
                            "the difference between hours and minutes. The ",
                            "fitted cells get an exact embedding; the rest get ",
                            "an out-of-sample projection onto it. ",
                            "[default %default]")),
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
  make_option("--marker_aliases", type = "character", default = NULL,
              help = paste0("Treat differently-spelled marker names as the same ",
                            "marker. Either inline \"INFy=IFN;TNFa=TNF\" or a ",
                            "path to a two-column file (as written, to use). ",
                            "Needed when workspaces disagree about spelling, ",
                            "which makes the shared panel come back EMPTY. The ",
                            "run prints the spellings it saw and suggests ",
                            "pairings; nothing is ever merged on its own.")),
  make_option("--no_recover_markers", action = "store_true", default = FALSE,
              help = paste0("Do NOT fill in missing marker names from the ",
                            "channel they were measured on. By default a file ",
                            "with no $PnS borrows its marker names from the ",
                            "files that do have them, keyed on detector — ",
                            "without which a cohort where half the files were ",
                            "acquired from an unfilled template has no usable ",
                            "channels at all. Turn this off to see the panel ",
                            "exactly as the files carry it.")),
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
  make_option("--normalize_data", type = "character", default = "none",
              help = paste0("Per-marker normalisation of the values the ",
                            "EMBEDDING is computed from: 'none', 'zscore', ",
                            "'minmax'. Changing this moves the points — it is ",
                            "a different analysis. [default %default]")),
  make_option("--normalize_umap", type = "character", default = "none",
              help = paste0("Per-marker normalisation of the values the MARKER ",
                            "COLOURS are taken from, applied after the ",
                            "embedding exists. Changing this recolours the ",
                            "same coordinates and moves nothing, which is what ",
                            "you want when comparing scales on one embedding. ",
                            "[default %default]")),
  make_option("--subpop_levels_select", type = "character", default = NULL,
              help = paste0("Name the gates that make up each level, instead ",
                            "of letting the script guess. Levels separated by ",
                            "';', gates within a level by ',', e.g. ",
                            "\"CD4+,CD8+;CD4+/CD107a+,CD8+/CD107a+\". Use this ",
                            "when a workspace hangs lineage and functional ",
                            "gates off the same node — then no automatic ",
                            "grouping is a partition and the level plots come ",
                            "out one colour. The run warns when that happens.")),
  make_option("--no_subpop_levels", action = "store_true", default = FALSE,
              help = paste0("Do NOT write one subpopulation plot per gating ",
                            "level. By default, in addition to the single ",
                            "'by subpopulation' plot (which labels each cell by ",
                            "its deepest gate), one plot is written per level ",
                            "under --population: level 1 is the immediate ",
                            "children, level 2 theirs, and so on. Cells of the ",
                            "reference population in no gate at a level are ",
                            "drawn grey and left out of the legend.")),
  make_option("--normalize", type = "character", default = NULL,
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
if (!dir.exists(shared_dir)) {
  stop("R_shared directory not found at: ", shared_dir)
}

source(file.path(shared_dir, "io_utils.R"))        # make_outdir
source(file.path(shared_dir, "plot_utils.R"))      # save_plot_formats, create_safe_filename
source(file.path(shared_dir, "memory_utils.R"))    # resolve_cores, parallel_lapply
source(file.path(shared_dir, "stats_utils.R"))     # run_umap
source(file.path(shared_dir, "dimred_dispatch.R")) # run_dimred (umap/tsne/dm)
source(file.path(shared_dir, "flowframe_utils.R")) # explore_data
# THE FlowJo workspace layer: wsp_open, wsp_to_gatingset, discovery, pairing
source(file.path(src_dir, "loading", "flowjo_workspace.R"))
# Per-sample annotations (patient, timepoint, ...) from a CSV or a name regex
source(file.path(src_dir, "loading", "sample_annotations.R"))
# Marker-panel reconciliation: channel -> marker recovery, aliases, comparison.
# Sourced AFTER shared_dir is set, which is what their source guards look for.
source(file.path(src_dir, "harmonization", "panel_homogenization.R"))
source(file.path(src_dir, "harmonization", "panel_comparison.R"))
source(file.path(src_dir, "visualization", "umap_load_report.R"))  # build/write JSON
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))  # umap_from_flowjo

# -----------------------------------------------------------------------------
# Parse the comma-separated CLI strings into the structured values the analysis
# module expects, then hand off. All analysis/plotting logic lives in the module.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Annotation preview: answer and exit, before anything slow happens.
#
# This deliberately sits ahead of umap_from_flowjo() — it needs neither a
# workspace nor an FCS file, only the sample NAMES, and the whole reason it
# exists is to be fast enough to run while someone is typing a path into a form.
# The matching itself is resolve_sample_annotations(), the same function a run
# uses, so the preview and the run cannot disagree.
# -----------------------------------------------------------------------------
if (!is.null(opt$annotation_preview)) {
  preview_annotations(
    out_path  = opt$annotation_preview,
    samples   = opt$preview_samples,
    metadata  = opt$sample_metadata,
    pattern   = opt$sample_pattern,
    id_column = opt$sample_id_column)
  quit(status = 0)
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
  normalize_data   = opt$normalize_data,
  normalize_umap   = opt$normalize_umap,
  subpop_levels    = !opt$no_subpop_levels,
  subpop_levels_select = opt$subpop_levels_select,
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
  list_populations = opt$list_populations,
  allow_missing_samples = opt$allow_missing_samples,
  population_map   = opt$population_map,
  match_mode       = opt$match_mode,
  marker_aliases   = opt$marker_aliases,
  workspaces       = opt$workspaces,
  n_cores          = opt$n_cores,
  embed_cores      = ifelse(opt$n_cores == 0, 1, opt$n_cores),
  umap_init        = opt$umap_init,
  umap_subsample   = opt$umap_subsample,
  recover_markers  = !opt$no_recover_markers,
  report_json      = opt$report_json
)
