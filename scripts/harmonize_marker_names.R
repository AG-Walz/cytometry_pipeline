#!/usr/bin/env Rscript
# =============================================================================
# scripts/harmonize_marker_names.R : CLI wrapper for the marker-name harmoniser
# =============================================================================
# Thin standalone wrapper (no Nextflow). Parses CLI options, sources the keyword
# layer, the pure alias functions and the tool module
# (src/metadata_qc/harmonize_marker_names.R), and calls
# harmonize_marker_names(). The logic lives in those modules; keep this thin.
#
# Base R apart from two optional packages: optparse for the option parser (a
# minimal --key value parser stands in when it is missing) and yaml, needed
# only for --canonical.
#
# WHAT IT IS FOR
# A cohort spells the same marker several ways in $PnS (INFy / IFN / IFNg,
# TNFa / TNF, ALD / Aqua LiveDead) and many files carry no $PnS at all. This
# writes one canonical spelling into every file, without decoding the events,
# and records the old values in a ledger so the run can be undone.
#
# Example:
#   # 1. Census: every channel with every spelling seen and how many files.
#   Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs --probe
#
#   # 2. Dry run against the drift config's `markers` block. Writes nothing.
#   Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs \
#     --canonical config/drift_mtb.yaml
#
#   # 3. Apply in place, confirming a low-confidence pair by hand.
#   Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs \
#     --canonical config/drift_mtb.yaml --aliases "INFy=IFNg" --apply --in_place
#
#   # Undo that run.
#   Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs --revert
# =============================================================================

option_spec <- list(
  fcs_dir   = list("data", "Folder holding the FCS files, searched recursively [default %default]"),
  probe     = list(FALSE, "Census only: channel x spelling x files. Writes nothing."),
  canonical = list(NULL, paste0("YAML file whose top level is a `markers` block (canonical -> ",
                                "regex patterns), or a drift config containing one.")),
  aliases   = list(NULL, paste0("Explicit renames, \"INFy=IFNg;TNF=TNFa\", or a two-column ",
                                "file (as written, to use). Take precedence over --canonical.")),
  apply     = list(FALSE, "Perform the writes. Without it the plan is printed and nothing changes."),
  in_place  = list(FALSE, "Patch the originals where they sit. Mutually exclusive with --out_dir."),
  out_dir   = list(NULL, "Write patched copies here, mirroring subfolders. Must not be inside --fcs_dir."),
  revert    = list(FALSE, "Undo the most recent run recorded in the ledger."),
  ledger    = list("marker_names_ledger.csv", "Undo ledger, written next to the data [default %default]"),
  no_recursive = list(FALSE, "Do not descend into subfolders."),
  src_dir   = list(NULL, "Path to the repo 'src' directory. Default: auto-detected.")
)
description <- paste0("\nWrite one canonical marker name into the $PnS keyword of every FCS ",
                      "file in a folder. Event data is never decoded. Start with --probe.\n")

if (requireNamespace("optparse", quietly = TRUE)) {
  opts <- lapply(names(option_spec), function(k) {
    d <- option_spec[[k]][[1]]; h <- option_spec[[k]][[2]]
    if (is.logical(d)) optparse::make_option(paste0("--", k), action = "store_true",
                                             default = FALSE, help = h)
    else optparse::make_option(paste0("--", k), type = "character", default = d, help = h)
  })
  opt <- optparse::parse_args(optparse::OptionParser(option_list = opts,
                                                     description = description))
} else {
  # Minimal --key value / --flag parser for a plain R install.
  opt <- lapply(option_spec, `[[`, 1L)
  usage <- paste0(description, "\nOptions:\n", paste(sprintf("  --%-13s %s", names(option_spec),
    vapply(option_spec, function(s) sub("%default", format(s[[1]]), s[[2]], fixed = TRUE),
           character(1))), collapse = "\n"))
  args <- commandArgs(trailingOnly = TRUE); i <- 1L
  while (i <= length(args)) {
    if (args[i] %in% c("-h", "--help")) { cat(usage, "\n"); quit(status = 0) }
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !key %in% names(option_spec)) {
      stop("unknown option: ", args[i], "\n", usage, call. = FALSE)
    }
    if (is.logical(opt[[key]])) { opt[[key]] <- TRUE; i <- i + 1L; next }
    if (i == length(args)) stop("--", key, " needs a value", call. = FALSE)
    opt[[key]] <- args[i + 1L]; i <- i + 2L
  }
}

# --- locate the repo and source the modules -----------------------------------
script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
cli_dir <- if (length(script_file) == 1) dirname(sub("^--file=", "", script_file)) else
  file.path(getwd(), "scripts")
source(file.path(cli_dir, "_cli_common.R"))
src_dir <- resolve_src_dir(opt$src_dir)
shared_dir <<- file.path(src_dir, "R_shared")
source(file.path(shared_dir, "fcs_keywords.R"))
source(file.path(src_dir, "harmonization", "panel_homogenization.R"))
source(file.path(src_dir, "harmonization", "drift", "drift_panel.R"))
source(file.path(src_dir, "metadata_qc", "harmonize_marker_names.R"))

# --- inputs --------------------------------------------------------------------
canonical <- NULL
if (!is.null(opt$canonical)) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("--canonical needs the 'yaml' package", call. = FALSE)
  }
  if (!file.exists(opt$canonical)) stop("--canonical file not found: ", opt$canonical, call. = FALSE)
  cfg <- yaml::read_yaml(opt$canonical)
  canonical <- if (!is.null(cfg$markers)) cfg$markers else cfg
  if (!is.list(canonical) || is.null(names(canonical))) {
    stop("--canonical must be a named `markers` block: canonical -> regex patterns", call. = FALSE)
  }
}
aliases <- parse_marker_aliases(opt$aliases)

if (opt$apply && !opt$in_place && is.null(opt$out_dir) && !opt$revert) {
  stop("refusing to guess where to write: pass --in_place to patch the ",
       "originals, or --out_dir <dir> to write copies.", call. = FALSE)
}
if (opt$in_place && !is.null(opt$out_dir)) {
  stop("--in_place and --out_dir are mutually exclusive.", call. = FALSE)
}

harmonize_marker_names(
  fcs_dir    = opt$fcs_dir,
  canonical  = canonical,
  aliases    = aliases,
  probe_only = opt$probe,
  apply      = opt$apply,
  in_place   = opt$in_place,
  out_dir    = opt$out_dir,
  revert     = opt$revert,
  ledger     = opt$ledger,
  recursive  = !opt$no_recursive
)
