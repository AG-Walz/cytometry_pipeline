#!/usr/bin/env Rscript
# =============================================================================
# scripts/sync_fcs_names.R — CLI wrapper for the FCS embedded-name sync
# =============================================================================
# Thin standalone wrapper (no Nextflow). Parses CLI options, locates the repo
# 'src' directory, sources the keyword layer and the tool module
# (src/metadata_qc/sync_fcs_names.R), and calls sync_fcs_names(). All of the
# logic lives in those two modules — keep this file thin.
#
# Unlike the other scripts here, the modules behind this one are BASE R ONLY:
# no conda, no Bioconductor, no flowCore. Only this wrapper needs optparse, so
# on a machine without it you can source the two modules and call
# sync_fcs_names() directly.
#
# WHAT IT IS FOR
# A run acquired under one file name (a plate loader, a Diva tube template)
# records that name inside every file, in $FIL. Rename the files on disk and
# $FIL still says the old thing — so FlowJo, which displays $FIL by default,
# shows a column of identical sample names. This writes each file's own name
# into its metadata WITHOUT decoding the event data.
#
# BEFORE YOU RUN IT ON A COHORT
# Open FlowJo -> Preferences -> Workspace -> "How should data files be named".
# If it says "Use Short File Name" you do not need this tool at all. If it says
# "Keyword Combination", find out which keywords, and pass them to --also.
#
# Example:
#   # 1. Census first. Read-only: layout, collisions, and which samples CytoML
#   #    is currently telling apart by luck.
#   Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --probe
#
#   # 2. Dry run: exactly what would change, per file. Still writes nothing.
#   Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed
#
#   # 3. Apply in place, recording an undo ledger next to the data.
#   Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --apply --in_place
#
#   # Undo that run.
#   Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --revert
#
#   # See what every file actually says about itself, without changing anything:
#   Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --probe \
#     --keyword_table keywords.csv
#
#   # A site whose FlowJo shows $SRC and Tube Name rather than $FIL:
#   Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --apply --in_place \
#     --also '$SRC,TUBE NAME'
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# -----------------------------------------------------------------------------
# CLI options
#
# Nothing is written unless --apply is given, and --apply still refuses to guess
# between editing the originals and copying them. The default mode is the
# read-only census, because "the conflict is $FIL" is an inference about
# FlowJo's behaviour and the census is how you check it before touching data.
# -----------------------------------------------------------------------------
option_list <- list(
  make_option("--fcs_dir", type = "character", default = "data",
              help = paste0("Folder holding the already-renamed FCS files, ",
                            "searched RECURSIVELY and case-insensitively ",
                            "(.fcs and .FCS) [default %default]")),
  make_option("--probe", action = "store_true", default = FALSE,
              help = paste0("Census only. Prints each file's layout, every ",
                            "identity keyword shared by more than one file, ",
                            "and any group of files sharing BOTH $FIL and ",
                            "$TOT — those are the ones CytoML currently tells ",
                            "apart by luck. Writes nothing.")),
  make_option("--list_keywords", action = "store_true", default = FALSE,
              help = "Print every keyword of the first few files, then exit."),
  make_option("--apply", action = "store_true", default = FALSE,
              help = paste0("Perform the writes. Without it the tool reports ",
                            "what it would do and stops. [default %default]")),
  make_option("--in_place", action = "store_true", default = FALSE,
              help = paste0("Patch the original files where they sit. Cheap ",
                            "even for multi-GB files: when the new name fits ",
                            "the existing space the event data is never read ",
                            "or moved. Mutually exclusive with --output_dir.")),
  make_option("--output_dir", type = "character", default = NULL,
              help = paste0("Write patched COPIES here, mirroring the ",
                            "subfolder layout, and leave the originals alone. ",
                            "Needs disk space for the whole dataset. Must not ",
                            "be inside --fcs_dir.")),
  make_option("--also", type = "character", default = NULL,
              help = paste0("Comma-separated extra keywords to give the same ",
                            "value, for sites whose FlowJo is set to \"Keyword ",
                            "Combination\", e.g. '$SRC,TUBE NAME'. Only $FIL ",
                            "is written by default.")),
  make_option("--set", type = "character", default = NULL,
              help = paste0("Semicolon-separated KEY=TEMPLATE pairs, e.g. ",
                            "\"$SRC=%stem%;TUBE NAME=%dir%\". Templates ",
                            "understand %name% %stem% %oldname% %relpath% ",
                            "%dir%.")),
  make_option("--unset", type = "character", default = NULL,
              help = "Comma-separated keywords to remove."),
  make_option("--guid", type = "character", default = "never",
              help = paste0("Whether to rewrite GUID: 'never', 'auto' (only ",
                            "when its current value is already an echo of the ",
                            "old file name) or 'always'. A real vendor GUID is ",
                            "left alone by 'auto', and a fresh UUID is never ",
                            "minted. [default %default]")),
  make_option("--strip_extension", action = "store_true", default = FALSE,
              help = paste0("Write $FIL without the .fcs suffix. The default ",
                            "keeps it, which is what $FIL means in the ",
                            "standard and what acquisition software writes.")),
  make_option("--keyword_table", type = "character", default = NULL,
              help = paste0("Write every keyword of every file to this CSV ",
                            "(long format: file, keyword, value, ",
                            "keyword_rank). Read-only, and done before any ",
                            "write, so it works with --probe too. The TEXT ",
                            "segment is a packed blob nothing else can read; ",
                            "this is how you see what the cohort says about ",
                            "itself.")),
  make_option("--no_original", action = "store_true", default = FALSE,
              help = paste0("Do NOT preserve the replaced $FIL in ORIGINALFIL. ",
                            "By default the name a file was acquired under is ",
                            "written to that keyword the first time $FIL is ",
                            "rewritten, and never overwritten afterwards.")),
  make_option("--revert", action = "store_true", default = FALSE,
              help = paste0("Undo a previous run by replaying its ledger. ",
                            "Restores every keyword this tool changed, and ",
                            "removes the ones it created.")),
  make_option("--allow", type = "character", default = NULL,
              help = paste0("Comma-separated refusal codes to downgrade, e.g. ",
                            "'crc_present'. Use only after reading why the ",
                            "file was refused.")),
  make_option("--ledger", type = "character", default = "fcs_name_sync_ledger.csv",
              help = paste0("Undo ledger, written next to the data (or in ",
                            "--output_dir). Semicolon-separated. ",
                            "[default %default]")),
  make_option("--no_recursive", action = "store_true", default = FALSE,
              help = "Do not descend into subfolders."),
  make_option("--src_dir", type = "character", default = NULL,
              help = paste0("Path to the repo 'src' directory. Default: ",
                            "auto-detected relative to this script."))
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = paste0("\nWrite each FCS file's own name into its $FIL ",
                       "keyword, so FlowJo stops showing a column of ",
                       "identically-named samples. Event data is never ",
                       "decoded. Start with --probe.\n")))

# -----------------------------------------------------------------------------
# Locate the repo and source the two modules
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
    # script is repo/scripts/sync_fcs_names.R -> repo/src
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

source(file.path(shared_dir, "io_utils.R"))      # make_outdir, die
# The FCS keyword layer: fcs_inspect, fcs_write_keywords, fcs_print_keywords
source(file.path(shared_dir, "fcs_keywords.R"))
source(file.path(src_dir, "metadata_qc", "sync_fcs_names.R"))  # sync_fcs_names

# -----------------------------------------------------------------------------
# Parse the comma/semicolon-separated CLI strings, then hand off.
# -----------------------------------------------------------------------------
split_csv <- function(s) {
  if (is.null(s)) return(character(0))
  v <- trimws(strsplit(s, ",")[[1]])
  v[nzchar(v)]
}

# KEY=VALUE pairs are semicolon-separated, because a keyword name may contain a
# space ("TUBE NAME") and a template may contain a comma.
parse_sets <- function(s) {
  if (is.null(s)) return(character(0))
  parts <- trimws(strsplit(s, ";")[[1]])
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(character(0))
  keys <- sub("=.*$", "", parts)
  vals <- sub("^[^=]*=", "", parts)
  bad <- !grepl("=", parts, fixed = TRUE)
  if (any(bad)) {
    stop("--set expects KEY=VALUE pairs separated by ';', got: ",
         paste(parts[bad], collapse = ", "), call. = FALSE)
  }
  structure(as.list(trimws(vals)), names = trimws(keys))
}

if (opt$apply && !opt$in_place && is.null(opt$output_dir) && !opt$revert) {
  die("refusing to guess where to write: pass --in_place to patch the ",
      "originals, or --output_dir <dir> to write copies.")
}
if (opt$in_place && !is.null(opt$output_dir)) {
  die("--in_place and --output_dir are mutually exclusive.")
}

sync_fcs_names(
  fcs_dir         = opt$fcs_dir,
  apply           = opt$apply,
  in_place        = opt$in_place,
  output_dir      = opt$output_dir,
  also            = split_csv(opt$also),
  set             = parse_sets(opt$set),
  unset           = split_csv(opt$unset),
  guid            = opt$guid,
  strip_extension = opt$strip_extension,
  recursive       = !opt$no_recursive,
  probe_only      = opt$probe,
  list_keywords   = opt$list_keywords,
  revert          = opt$revert,
  allow           = split_csv(opt$allow),
  ledger          = opt$ledger,
  keyword_table   = opt$keyword_table,
  write_original  = !opt$no_original
)
