#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option("--wsp", type = "character", default = "data/Killing_AML.wsp",
              help = "Path to the FlowJo .wsp workspace file [default %default]"),
  make_option("--fcs_dir", type = "character", default = "data",
              help = "Directory containing the FCS files referenced by the wsp [default %default]"),
  make_option("--population", type = "character",
              default = "Lymphocytes/Single Cells/Viable cells",
              help = paste0("Population path to extract. A leading '/' is ",
                            "optional; a suffix (e.g. 'Single Cells/Viable cells') ",
                            "is matched against the gating tree. [default %default]")),
  make_option("--list_populations", action = "store_true", default = FALSE,
              help = "Print the population tree (gs_get_pop_paths) and exit."),
  make_option("--group", type = "integer", default = 1,
              help = paste0("FlowJo workspace group index to load ",
                            "(flowjo_to_gatingset name=) [default %default]")),
  make_option("--expected_samples", type = "integer", default = 0,
              help = paste0("If > 0, warn when the loaded group does not ",
                            "contain exactly this many samples [default %default]")),
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
              help = "Plot height (inches) [default %default]")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = paste0(
    "\nCompute a single shared UMAP over one FlowJo-gated population across ",
    "samples and plot it by sample, by subpopulation, and by every marker.\n")
))

src_dir    <- Sys.getenv("CYTOMETRY_SRC", "/opt/cytometry/src")
shared_dir <<- file.path(src_dir, "R_shared")

source(file.path(shared_dir, "io_utils.R"))        # make_outdir
source(file.path(shared_dir, "plot_utils.R"))      # save_plot_formats, create_safe_filename
source(file.path(shared_dir, "stats_utils.R"))     # run_umap
source(file.path(shared_dir, "flowframe_utils.R")) # explore_data
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))  # umap_from_flowjo

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
cat("FlowJo UMAP plots ->", opt$output_dir, "\n")
