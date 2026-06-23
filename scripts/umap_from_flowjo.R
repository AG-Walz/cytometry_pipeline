#!/usr/bin/env Rscript
# =============================================================================
# scripts/umap_from_flowjo.R — Shared UMAP from a FlowJo workspace population
# =============================================================================
# Standalone script (no Nextflow). Opens a FlowJo .wsp, builds a GatingSet over
# a directory of FCS files, extracts the cells of ONE selected population from
# each sample (data is already compensated + transformed in the FlowJo working
# space, so NO further transform is applied), concatenates the per-sample
# expression matrices, computes a SINGLE shared UMAP embedding, and plots that
# embedding coloured by sample and by every marker (conjoint + per-file).
#
# Reuses the repo's R_shared utilities (run_umap, save_plot_formats,
# create_safe_filename, make_outdir, explore_data).
#
# Index-alignment contract (CRITICAL):
#   combined_matrix   : rows = cells (sample 1 stacked on top of sample 2)
#   file_labels       : length == nrow(combined_matrix); file_labels[i] is the
#                       source file of combined_matrix row i
#   umap_coords       : same row order as combined_matrix (UMAP is row-stable)
#   per-marker colour : combined_matrix[, marker] -> same row order
# All subsetting for per-file plots is done by `which(file_labels == f)`, which
# preserves the shared embedding coordinates for those exact rows. We assert
# all(is.finite(combined_matrix)) BEFORE run_umap so that a channel with
# NA/NaN/Inf cannot make uwot drop rows and silently break this contract.
#
# Example:
#   Rscript scripts/umap_from_flowjo.R \
#     --wsp data/Killing_AML.wsp --fcs_dir data \
#     --population "Lymphocytes/Single Cells/Viable cells" \
#     --output_dir results/umap_viable
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
# bare `Rscript scripts/umap_from_flowjo.R` targets the intended analysis. The
# file.exists/dir.exists checks below still validate whatever is supplied.
# -----------------------------------------------------------------------------
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
  make_option("--output_dir", type = "character", default = "umap_flowjo",
              help = "Output directory for plots [default %default]"),
  make_option("--output_formats", type = "character", default = "png",
              help = "Comma-separated formats: png,pdf,svg [default %default]"),
  make_option("--downsample", type = "integer", default = 0,
              help = paste0("Max cells per sample AFTER gating (0 = no ",
                            "downsampling) [default %default]")),
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
    "samples and plot it by sample and by every marker.\n")
))

# -----------------------------------------------------------------------------
# Locate the repo 'src' dir, set the shared_dir/src_dir globals (CLI convention),
# and source the shared utilities in dependency order.
#
# The R_shared files we use (io_utils, plot_utils, stats_utils, flowframe_utils)
# have no problematic top-level source() calls, but we still set `shared_dir`
# and `src_dir` as globals exactly like nextflow/bin/*_cli.R do, so any helper
# that inspects them behaves consistently.
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

suppressPackageStartupMessages({
  library(CytoML)
  library(flowWorkspace)
  library(flowCore)
  library(ggplot2)
})
if (!requireNamespace("viridis", quietly = TRUE)) {
  stop("Package 'viridis' is required for marker colour scales.")
}

# -----------------------------------------------------------------------------
# Validate required args
# -----------------------------------------------------------------------------
if (!file.exists(opt$wsp)) stop("wsp file not found: ", opt$wsp)
if (!dir.exists(opt$fcs_dir)) stop("fcs_dir not found: ", opt$fcs_dir)

wsp_path <- normalizePath(opt$wsp, mustWork = TRUE)
fcs_dir  <- normalizePath(opt$fcs_dir, mustWork = TRUE)

parse_formats <- function(s) {
  parts <- trimws(strsplit(s, ",")[[1]])
  list(
    pdf = "pdf" %in% parts,
    png = "png" %in% parts,
    svg = "svg" %in% parts
  )
}
output_formats <- parse_formats(opt$output_formats)

# -----------------------------------------------------------------------------
# 1. Open workspace and build the GatingSet over the FCS dir.
# -----------------------------------------------------------------------------
cat("Opening FlowJo workspace:", wsp_path, "\n")
ws <- CytoML::open_flowjo_xml(wsp_path)
cat("Building GatingSet from FCS dir:", fcs_dir,
    "(group index =", opt$group, ")\n")
gs <- CytoML::flowjo_to_gatingset(ws, name = opt$group, path = fcs_dir)

sample_names <- flowWorkspace::sampleNames(gs)
cat("Samples in loaded group (", length(sample_names), "):\n", sep = "")
for (s in sample_names) cat("  -", s, "\n")
if (opt$expected_samples > 0 && length(sample_names) != opt$expected_samples) {
  warning(sprintf(
    "Loaded group %d has %d sample(s) but --expected_samples=%d. The wrong ",
    opt$group, length(sample_names), opt$expected_samples),
    "workspace group may be selected; re-check --group or --list_populations.")
}

all_pops <- flowWorkspace::gs_get_pop_paths(gs, path = "full")

# -----------------------------------------------------------------------------
# 2. --list_populations mode: print the tree and exit.
# -----------------------------------------------------------------------------
if (isTRUE(opt$list_populations)) {
  cat("\nPopulation tree (", length(all_pops), " nodes):\n", sep = "")
  for (p in all_pops) cat("  ", p, "\n", sep = "")
  cat("\nSelect one with --population (leading '/' optional), e.g.:\n")
  cat("  --population \"", sub("^/", "", all_pops[length(all_pops)]),
      "\"\n", sep = "")
  quit(status = 0)
}

if (is.null(opt$population) || !nzchar(opt$population)) {
  stop("--population is required (or use --list_populations to see options).")
}

# Resolve the requested population to a path that exists in the gating tree.
# Accept with/without leading '/', and also match by basename or suffix.
resolve_population <- function(requested, available) {
  req <- requested
  cand <- unique(c(req, sub("^/", "", req), paste0("/", sub("^/", "", req))))
  hit <- cand[cand %in% available]
  if (length(hit) >= 1) return(hit[1])
  # suffix match (e.g. user gave 'Lymphocytes/Single Cells/Viable cells')
  norm_req <- sub("^/", "", req)
  suffix_hits <- available[
    vapply(available, function(p) {
      endsWith(sub("^/", "", p), norm_req)
    }, logical(1))
  ]
  if (length(suffix_hits) == 1) return(suffix_hits[1])
  if (length(suffix_hits) > 1) {
    stop("Population '", requested, "' is ambiguous; matches:\n  ",
         paste(suffix_hits, collapse = "\n  "),
         "\nRe-run with the full path.")
  }
  stop("Population '", requested, "' not found. Run with ",
       "--list_populations to see the available paths.")
}

population <- resolve_population(opt$population, all_pops)
cat("\nTarget population resolved to:", population, "\n")

# Optional marker subset from the CLI.
markers_subset <- NULL
if (!is.null(opt$markers)) {
  markers_subset <- trimws(strsplit(opt$markers, ",")[[1]])
  markers_subset <- markers_subset[nzchar(markers_subset)]
}

# -----------------------------------------------------------------------------
# 3 & 4. Extract the population per sample and build per-sample marker matrices.
#
# We use explore_data(ff)$markers (per-column marker names, "" when none) for
# consistency with the rest of the repo. A channel is kept as a biological
# marker only when:
#   * it has a non-empty marker name, AND
#   * its CHANNEL name does not match the scatter/Time exclusion regex.
# The exclusion regex is applied to the CHANNEL name only (not the marker
# label), so a legitimate antibody whose marker string happens to contain
# "FSC"/"Time"/etc. as a substring is not silently dropped. Empty-marker
# scatter/Time channels are already removed by the marker != "" test.
# Kept columns are renamed to their MARKER name.
# -----------------------------------------------------------------------------
EXCLUDE_CHANNEL_PATTERN <- "Time|FSC|SSC|scatter"  # case-insensitive, CHANNEL only

extract_sample_marker_matrix <- function(gh) {
  # NOTE: data is already compensated + transformed in the FlowJo working
  # space. We take exprs() AS-IS, no further transform.
  ff   <- flowWorkspace::gh_pop_get_data(gh, y = population)
  expr <- flowCore::exprs(ff)
  info <- explore_data(ff)          # dyes (channels), markers, mcols

  channels <- colnames(expr)
  markers  <- info$markers          # length == ncol(expr)

  keep <- vapply(seq_along(channels), function(i) {
    m <- markers[i]
    if (is.na(m) || m == "") return(FALSE)                    # empty marker
    if (grepl(EXCLUDE_CHANNEL_PATTERN, channels[i],
              ignore.case = TRUE)) return(FALSE)              # scatter/Time channel
    TRUE
  }, logical(1))

  # Log excluded channels (channel + marker) so a wrongly-dropped biological
  # channel is visible to the user.
  if (any(!keep)) {
    dropped <- vapply(which(!keep), function(i) {
      sprintf("%s [%s]", channels[i],
              if (is.na(markers[i]) || markers[i] == "") "<no marker>" else markers[i])
    }, character(1))
    cat("    excluded channels:", paste(dropped, collapse = ", "), "\n")
  }

  if (!any(keep)) {
    stop("No biological marker columns found for a sample after excluding ",
         "scatter/time/empty channels.")
  }

  mtx <- expr[, keep, drop = FALSE]
  colnames(mtx) <- markers[keep]
  # Guard against duplicate marker names (keep first occurrence).
  if (anyDuplicated(colnames(mtx))) {
    dup <- colnames(mtx)[duplicated(colnames(mtx))]
    warning("Duplicate marker names dropped (kept first): ",
            paste(unique(dup), collapse = ", "))
    mtx <- mtx[, !duplicated(colnames(mtx)), drop = FALSE]
  }
  mtx
}

set.seed(opt$seed)  # governs downsampling; run_umap re-seeds internally for UMAP

per_sample_mtx <- list()
for (s in sample_names) {
  gh  <- gs[[s]]
  mtx <- extract_sample_marker_matrix(gh)

  # Optional per-sample downsampling AFTER gating.
  if (opt$downsample > 0 && nrow(mtx) > opt$downsample) {
    idx <- sort(sample(nrow(mtx), opt$downsample))
    mtx <- mtx[idx, , drop = FALSE]
  }

  cat(sprintf("  %s: %d cells x %d markers\n", s, nrow(mtx), ncol(mtx)))
  if (nrow(mtx) == 0) {
    warning("Sample '", s, "' has 0 cells in population '", population,
            "' after gating/downsampling; it will not contribute rows.")
  }
  per_sample_mtx[[s]] <- mtx
}

# Drop samples that contributed no cells (keeps alignment + plotting sane).
nonempty <- vapply(per_sample_mtx, nrow, integer(1)) > 0
if (!any(nonempty)) {
  stop("No sample has any cells in population '", population, "'.")
}
per_sample_mtx <- per_sample_mtx[nonempty]

# Intersect markers common to ALL (non-empty) samples.
common_markers <- Reduce(intersect, lapply(per_sample_mtx, colnames))
if (length(common_markers) == 0) {
  stop("No markers are common to all samples.")
}
# Preserve first-sample ordering for determinism.
common_markers <- colnames(per_sample_mtx[[1]])[
  colnames(per_sample_mtx[[1]]) %in% common_markers]

# Apply optional CLI marker subset (must be a subset of common markers).
if (!is.null(markers_subset)) {
  missing_m <- setdiff(markers_subset, common_markers)
  if (length(missing_m) > 0) {
    stop("These --markers are not common markers in the data: ",
         paste(missing_m, collapse = ", "),
         "\nAvailable common markers: ",
         paste(common_markers, collapse = ", "))
  }
  use_markers <- markers_subset
} else {
  use_markers <- common_markers
}

cat("\nMarkers used (", length(use_markers), "):\n  ",
    paste(use_markers, collapse = ", "), "\n", sep = "")

# Restrict every per-sample matrix to use_markers (same column order).
per_sample_mtx <- lapply(per_sample_mtx, function(m) m[, use_markers, drop = FALSE])

# -----------------------------------------------------------------------------
# 5. Concatenate. file_labels is index-aligned to combined_matrix rows.
# -----------------------------------------------------------------------------
combined_matrix <- do.call(rbind, per_sample_mtx)
file_labels     <- rep(names(per_sample_mtx),
                       vapply(per_sample_mtx, nrow, integer(1)))
stopifnot(length(file_labels) == nrow(combined_matrix))

cat(sprintf("\nConcatenated matrix: %d cells x %d markers (%d samples)\n",
            nrow(combined_matrix), ncol(combined_matrix),
            length(per_sample_mtx)))

# -----------------------------------------------------------------------------
# 6. Compute UMAP ONCE on the concatenated matrix.
#
# CRITICAL alignment guard: uwot may error or drop rows on non-finite input,
# which would silently break the row-for-row alignment between umap_coords,
# combined_matrix, and file_labels that the per-file subsetting relies on.
# Assert finiteness so any such problem fails loudly here.
# -----------------------------------------------------------------------------
if (!all(is.finite(combined_matrix))) {
  bad_cols <- colnames(combined_matrix)[
    apply(combined_matrix, 2, function(x) any(!is.finite(x)))]
  stop("combined_matrix contains non-finite (NA/NaN/Inf) values in marker(s): ",
       paste(bad_cols, collapse = ", "),
       ". Refusing to run UMAP because dropped rows would break index ",
       "alignment with the per-file plots.")
}

cat("Computing shared UMAP (n_neighbors=", opt$n_neighbors,
    ", min_dist=", opt$min_dist, ", metric=", opt$metric,
    ", seed=", opt$seed, ")...\n", sep = "")
umap_coords <- run_umap(
  combined_matrix,
  n_neighbors  = opt$n_neighbors,
  min_dist     = opt$min_dist,
  n_components = 2,
  metric       = opt$metric,
  seed         = opt$seed
)
stopifnot(nrow(umap_coords) == nrow(combined_matrix))

# Master data frame: one row per cell, aligned to combined_matrix.
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  file  = file_labels,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 7. Plotting.
# -----------------------------------------------------------------------------
make_outdir(opt$output_dir)
pt_size  <- opt$point_size
pt_alpha <- opt$point_alpha

pop_tag <- create_safe_filename(sub("^/", "", population))

# Shared axis limits so every plot (conjoint + per-file) uses the same frame.
xlim_shared <- range(umap_df$UMAP1, na.rm = TRUE)
ylim_shared <- range(umap_df$UMAP2, na.rm = TRUE)

save_one <- function(p, basename_no_ext) {
  save_plot_formats(
    p,
    file.path(opt$output_dir, basename_no_ext),
    width   = opt$plot_width,
    height  = opt$plot_height,
    formats = output_formats
  )
}

# 7a. Categorical plot coloured by file/sample (all cells).
p_file <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = file)) +
  geom_point(size = pt_size, alpha = pt_alpha) +
  coord_cartesian(xlim = xlim_shared, ylim = ylim_shared) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = paste0("UMAP - ", population, " - by file"),
       color = "File")
save_one(p_file, create_safe_filename(paste0("umap_", pop_tag, "_by_file")))
cat("Saved: by-file categorical plot\n")

files_unique <- unique(file_labels)

# 7b & 7c. For each marker: one conjoint plot (all cells) + one plot per file.
for (marker in use_markers) {
  # CONJOINT (all cells) — colour vector is combined_matrix[, marker], which is
  # index-aligned to umap_df rows.
  df_all <- umap_df
  df_all$expr <- combined_matrix[, marker]

  p_all <- ggplot(df_all, aes(x = UMAP1, y = UMAP2, color = expr)) +
    geom_point(size = pt_size, alpha = pt_alpha) +
    viridis::scale_color_viridis(name = marker) +
    coord_cartesian(xlim = xlim_shared, ylim = ylim_shared) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(title = paste0("UMAP - ", population, " - ", marker, " (all files)"))
  save_one(
    p_all,
    create_safe_filename(paste0("umap_", pop_tag, "_marker_", marker, "_ALL"))
  )

  # Consistent colour limits across the per-file plots = full marker range.
  # Guard against a degenerate (all-NA) column: range(..., na.rm=TRUE) would
  # return c(Inf, -Inf) and break scale_color_viridis(limits=...).
  expr_range <- range(combined_matrix[, marker], na.rm = TRUE)
  if (!all(is.finite(expr_range))) expr_range <- NULL

  # PER FILE — subset rows by file label; SAME coordinates, SAME colour scale.
  for (f in files_unique) {
    rows <- which(file_labels == f)               # index-aligned subset
    if (length(rows) == 0) next                    # nothing to plot for this file
    df_f <- umap_df[rows, , drop = FALSE]
    df_f$expr <- combined_matrix[rows, marker]    # aligned colour vector

    p_f <- ggplot(df_f, aes(x = UMAP1, y = UMAP2, color = expr)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      viridis::scale_color_viridis(name = marker, limits = expr_range) +
      coord_cartesian(xlim = xlim_shared, ylim = ylim_shared) +
      theme_minimal() +
      theme(panel.grid = element_blank()) +
      labs(title = paste0("UMAP - ", population, " - ", marker, " - ", f))
    save_one(
      p_f,
      create_safe_filename(
        paste0("umap_", pop_tag, "_marker_", marker, "_file_", f))
    )
  }
  cat("Saved marker plots for:", marker, "\n")
}

cat("\nDone. Plots written to:", normalizePath(opt$output_dir), "\n")