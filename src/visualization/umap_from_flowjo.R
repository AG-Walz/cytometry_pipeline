#!/usr/bin/env Rscript
# =============================================================================
# visualization/umap_from_flowjo.R — Shared UMAP from a FlowJo workspace population
# =============================================================================
# Core logic for the FlowJo-population UMAP analysis. Opens ONE OR MANY FlowJo
# .wsp workspaces, builds a GatingSet per workspace over the FCS files it
# references, extracts the cells of ONE selected population from each sample
# (data is already compensated + transformed in the FlowJo working space, so NO
# further transform is applied), concatenates the per-sample expression
# matrices, OPTIONALLY normalizes each marker using statistics computed from
# those concatenated cells, computes a SINGLE shared UMAP embedding, and plots
# that embedding coloured by sample, by workspace, by the FlowJo subgates found
# under the selected population (one colour per subgate; see `subpop_mode`), and
# by every marker (conjoint + per-file).
#
# Many workspaces go into ONE embedding, so their coordinates are comparable.
# Workspace discovery and the pairing of each workspace to its FCS files live in
# src/loading/flowjo_workspace.R; the two gating trees are reconciled
# through the gate paths RELATIVE to the selected population, which is the only
# key comparable across workspaces that name their ancestors differently.
# NOTHING here corrects for batch effects between workspaces (see
# src/harmonization/); the by-workspace plot is what makes them visible.
#
# Every plot uses an EQUAL-ASPECT coordinate system (coord_fixed, ratio = 1), so
# the embedding is never geometrically distorted and a wide legend (e.g. long
# file names) can no longer squish the panel — categorical figures are widened to
# give the legend its own space instead.
#
# This file holds only function definitions (no top-level side effects). The
# caller (scripts/umap_from_flowjo.R, or a CLI wrapper) is responsible for
# sourcing the R_shared utilities it relies on — run_umap (stats_utils.R),
# save_plot_formats + create_safe_filename (plot_utils.R), make_outdir
# (io_utils.R), explore_data (flowframe_utils.R) — before calling
# umap_from_flowjo().
#
# Index-alignment contract (CRITICAL):
#   combined_matrix   : rows = cells (sample 1 stacked on top of sample 2,
#                       workspace after workspace)
#   file_labels       : length == nrow(combined_matrix); file_labels[i] is the
#                       source file of combined_matrix row i
#   wsp_labels_per_cell : length == nrow(combined_matrix); owning workspace of
#                       row i. Built from the SAME per-sample order and row
#                       counts as file_labels, so the two cannot drift.
#   combined_subpop   : length == nrow(combined_matrix); subgate label of row i
#   umap_coords       : same row order as combined_matrix (UMAP is row-stable)
#   per-marker colour : combined_matrix[, marker] -> same row order
# All subsetting for per-file plots is done by `which(file_labels == f)`, which
# preserves the shared embedding coordinates for those exact rows. We assert
# all(is.finite(combined_matrix)) BEFORE run_umap so that a channel with
# NA/NaN/Inf cannot make uwot drop rows and silently break this contract.
# =============================================================================

# Catch-all subgate labels (a cell in no / more than one selected subgate).
.SUBPOP_NONE  <- "(none)"
.SUBPOP_MULTI <- "(multiple)"

# -----------------------------------------------------------------------------
# Resolve the requested population to a path that exists in the gating tree.
# Accept with/without leading '/', and also match by basename or suffix.
#
# Non-throwing: returns list(ok, path, reason). With several workspaces the
# same request has to be resolved against several trees, and every failure
# should be reported together rather than one stop() at a time — so the caller
# decides when to give up. resolve_population() below keeps the throwing
# behaviour its existing callers rely on.
# -----------------------------------------------------------------------------
try_resolve_population <- function(requested, available) {
  req <- requested
  cand <- unique(c(req, sub("^/", "", req), paste0("/", sub("^/", "", req))))
  hit <- cand[cand %in% available]
  if (length(hit) >= 1) return(list(ok = TRUE, path = hit[1], reason = ""))
  # suffix match (e.g. user gave 'Lymphocytes/Single Cells/Viable cells')
  norm_req <- sub("^/", "", req)
  suffix_hits <- available[
    vapply(available, function(p) {
      endsWith(sub("^/", "", p), norm_req)
    }, logical(1))
  ]
  if (length(suffix_hits) == 1) {
    return(list(ok = TRUE, path = suffix_hits[1], reason = ""))
  }
  if (length(suffix_hits) > 1) {
    return(list(ok = FALSE, path = NA_character_, reason = paste0(
      "population '", requested, "' is ambiguous; it matches:\n      ",
      paste(suffix_hits, collapse = "\n      "),
      "\n    Re-run with the full path.")))
  }
  list(ok = FALSE, path = NA_character_, reason = paste0(
    "population '", requested, "' not found."))
}

resolve_population <- function(requested, available) {
  res <- try_resolve_population(requested, available)
  if (!res$ok) {
    stop(sub("^population", "Population", res$reason),
         " Run with list_populations = TRUE to see the available paths.")
  }
  res$path
}

# -----------------------------------------------------------------------------
# Parse a per-workspace population mapping: which gate each workspace
# contributes to the shared embedding.
#
# Why this exists: one --population has to resolve in EVERY workspace, so a
# study whose analysts named the same biology differently ("Live/CD3" here,
# "Viable/CD3+" there) cannot be pooled at all — and with a generic "All
# Samples" group the gates common to every workspace can be three nodes deep.
# The mapping makes the pairing EXPLICIT: you say which gate in each workspace
# is the population, and workspaces you leave out are excluded from the run.
#
# That explicitness is the point. Pooling cells from gates that are not the same
# gate is a scientific error, and the only one who can say whether two
# differently-named gates are the same biology is the person who drew them.
#
# Two forms, because a FlowJo gate name may contain a comma but a file never
# has to be quoted:
#   * a path to a two-column file (TSV, or CSV/semicolon), workspace then gate,
#     an optional header row starting "workspace" being ignored
#   * inline "label=/gate;label2=/other gate", split on ';' and the FIRST '='
#
# Returns a named character vector (names = workspace labels), or NULL for "".
# -----------------------------------------------------------------------------
parse_population_map <- function(spec, labels = NULL) {
  if (is.null(spec) || !nzchar(spec)) return(NULL)

  if (file.exists(spec)) {
    txt <- readLines(spec, warn = FALSE)
    txt <- trimws(txt)
    txt <- txt[nzchar(txt) & !startsWith(txt, "#")]
    if (!length(txt)) stop("population_map file is empty: ", spec)
    sep <- if (any(grepl("\t", txt))) "\t" else if (any(grepl(";", txt))) ";" else ","
    parts <- strsplit(txt, sep, fixed = TRUE)
    bad <- vapply(parts, function(p) length(p) < 2L, logical(1))
    if (any(bad)) {
      stop("population_map file ", spec, " has ", sum(bad), " line(s) without ",
           "two '", sep, "'-separated fields, e.g.:\n  ", txt[which(bad)[1]])
    }
    ws  <- trimws(vapply(parts, function(p) p[1], character(1)))
    pop <- trimws(vapply(parts, function(p) paste(p[-1], collapse = sep),
                         character(1)))
    if (tolower(ws[1]) %in% c("workspace", "wsp", "label")) {
      ws <- ws[-1]; pop <- pop[-1]
    }
  } else {
    pairs <- trimws(strsplit(spec, ";", fixed = TRUE)[[1]])
    pairs <- pairs[nzchar(pairs)]
    if (!length(pairs)) return(NULL)
    has_eq <- grepl("=", pairs, fixed = TRUE)
    if (!all(has_eq)) {
      stop("population_map entries must be 'workspace=gate' separated by ';'. ",
           "This one has no '=':\n  ", pairs[!has_eq][1],
           "\n(If you meant a file, it does not exist: ", spec, ")")
    }
    ws  <- trimws(sub("=.*$", "", pairs))
    pop <- trimws(sub("^[^=]*=", "", pairs))
  }

  if (!length(ws)) return(NULL)
  if (any(!nzchar(pop))) {
    stop("population_map has no gate for workspace(s): ",
         paste(ws[!nzchar(pop)], collapse = ", "))
  }
  if (anyDuplicated(ws)) {
    dup <- unique(ws[duplicated(ws)])
    stop("population_map names the same workspace more than once: ",
         paste(dup, collapse = ", "),
         "\nOne gate per workspace — pool several gates by selecting their ",
         "common parent instead.")
  }
  if (!is.null(labels)) {
    unknown <- setdiff(ws, labels)
    if (length(unknown)) {
      stop("population_map names workspace(s) that were not loaded: ",
           paste(unknown, collapse = ", "),
           "\nLoaded workspaces are:\n  ", paste(labels, collapse = "\n  "))
    }
  }
  stats::setNames(pop, ws)
}

# -----------------------------------------------------------------------------
# Discover the candidate FlowJo subgates UNDER `population` to colour by.
#
# `pop_paths_union` should be the UNION of paths across all samples, because a
# FlowJo workspace can attach slightly different gating trees to different
# samples; the union guarantees the colour legend covers every subgate present
# in any sample. `all_pops` (the gs-level tree order) is used only to order the
# legend entries by tree position.
#
# Returns list(paths, labels, depths, tag), where labels are paths relative to
# `population` and depths are their depth below it. `mode`:
#   furthest - every descendant gate (assignment labels each cell by the most
#              specific / furthest node it reaches; see assign_subpop)
#   cd4cd8   - the CD4+/CD8+ lineage split (a clean partition)
#   direct   - every immediate child gate
#   deepest  - every descendant gate (assignment reduces to the deepest by path)
#   none     - no subgates
# `select` (non-NULL) overrides `mode`: an explicit list of gate identifiers
# matched against the relative path or the bare gate name.
# -----------------------------------------------------------------------------
flowjo_subgate_candidates <- function(pop_paths_union, all_pops, population,
                                      mode = "cd4cd8", select = NULL,
                                      warn = TRUE) {
  prefix    <- paste0(population, "/")
  desc      <- pop_paths_union[startsWith(pop_paths_union, prefix)]
  rel       <- substring(desc, nchar(prefix) + 1L)
  depth     <- lengths(strsplit(rel, "/", fixed = TRUE))
  base_name <- basename(rel)  # leaf gate name (last path component)

  if (!is.null(select)) {
    keep <- rel %in% select | base_name %in% select
    matched <- unique(c(rel[keep], base_name[keep]))
    missing_sel <- setdiff(select, matched)
    # warn = FALSE when several workspaces are being merged: a gate absent from
    # THIS workspace may well exist in another, so the warning is only correct
    # once the candidate sets have been merged (see merge_subgate_candidates).
    if (warn && length(missing_sel)) {
      warning("subpop_select entries not found under ", population, ": ",
              paste(missing_sel, collapse = ", "))
    }
    tag <- "select"
  } else {
    keep <- switch(mode,
                   furthest = rep(TRUE, length(desc)),
                   cd4cd8   = base_name %in% c("CD4+", "CD8+"),
                   direct   = depth == 1L,
                   deepest  = rep(TRUE, length(desc)),
                   none     = rep(FALSE, length(desc)))
    tag <- mode
  }

  paths  <- desc[keep]
  labels <- rel[keep]
  depths <- depth[keep]

  # Order legend entries by tree position (fall back to label for sample-specific
  # extras absent from the gs-level path list), so colours follow the gating tree.
  ord <- order(match(paths, all_pops, nomatch = length(all_pops) + 1L), labels)
  list(paths = paths[ord], labels = labels[ord], depths = depths[ord], tag = tag)
}

# -----------------------------------------------------------------------------
# Merge the per-workspace candidate sets into ONE colour partition.
#
# Across workspaces the ABSOLUTE gate paths cannot be compared: two workspaces
# may gate the same biology under differently-named ancestors, so --population
# legitimately resolves to '/Lymphocytes/Single Cells/Viable cells' in one and
# '/Lymph/Singlets/Viable cells' in another. What IS comparable — and what the
# colour actually means — is the path RELATIVE to the selected population. So
# the relative label is the merge key, and depth is derived from it.
#
# Returns list(labels, depths, paths, tag) where `paths` is a named list, one
# absolute-path vector per workspace, COLUMN-ALIGNED to `labels`. A gate that a
# workspace does not have is NA there; get_subpop_membership() leaves those
# columns all-FALSE, exactly as it already does for a gate missing from one
# sample's tree.
#
# Label order is first-appearance across workspaces, so each workspace's own
# gating-tree order is preserved as far as possible and the legend is stable.
# -----------------------------------------------------------------------------
merge_subgate_candidates <- function(cands) {
  cands <- cands[!vapply(cands, is.null, logical(1))]
  if (length(cands) == 0) {
    return(list(labels = character(0), depths = integer(0),
                paths = list(), tag = "none"))
  }

  labels <- unique(unlist(lapply(cands, function(c) c$labels), use.names = FALSE))
  if (is.null(labels)) labels <- character(0)
  depths <- lengths(strsplit(labels, "/", fixed = TRUE))

  paths <- lapply(cands, function(c) {
    out <- rep(NA_character_, length(labels))
    hit <- match(labels, c$labels)
    out[!is.na(hit)] <- c$paths[hit[!is.na(hit)]]
    out
  })
  names(paths) <- names(cands)

  list(labels = labels, depths = as.integer(depths), paths = paths,
       tag = cands[[1]]$tag)
}

# -----------------------------------------------------------------------------
# Per-cell subpopulation label from a (cells x candidate-gate) membership matrix.
# Fully vectorised. `assignment` selects the rule:
#
#   "specificity" (the 'furthest' mode) — label each cell with the MOST SPECIFIC
#       gate it reaches: the gate containing the fewest cells. A node deeper in
#       the gating tree is a SUBSET of its parents, so it has fewer cells; this
#       picks e.g. CD4+/CD107a-IFNg-IL2-TNF- over the overlapping single-marker
#       gates CD107a-/IFNg-/... that a cell sits in at the same path depth. A
#       cell that only reaches a shallower parent (and no deeper node) is
#       labelled by that parent; a cell in no gate is none_label. No "(multiple)"
#       — every reached cell gets the single furthest node. Path depth breaks
#       exact size ties toward the deeper node.
#
#   "depth" — reduce each cell to its deepest gate(s) BY PATH DEPTH, then a
#       single survivor -> that label, ties -> multiple_label, none -> none_label.
#       For 'direct'/'cd4cd8' candidates all depths are equal, so this collapses
#       to "exactly one gate -> that gate, else multiple".
# -----------------------------------------------------------------------------
assign_subpop <- function(mem, subpop_depths, subpop_labels,
                          assignment = "depth",
                          none_label = .SUBPOP_NONE,
                          multiple_label = .SUBPOP_MULTI) {
  n <- nrow(mem)
  if (n == 0L) return(character(0))
  if (ncol(mem) == 0L) return(rep(none_label, n))

  if (assignment == "specificity") {
    sizes <- colSums(mem)                         # per-sample gate sizes
    pri   <- sizes - 1e-6 * subpop_depths         # smaller gate wins; deeper breaks ties
    P <- matrix(pri, nrow = n, ncol = ncol(mem), byrow = TRUE)
    P[!mem] <- Inf                                # only member gates are eligible
    has <- rowSums(mem) > 0
    chosen <- max.col(-P, ties.method = "first")  # min priority among member gates
    out <- rep(none_label, n)
    out[has] <- subpop_labels[chosen[has]]
    return(out)
  }

  depth_mat <- matrix(subpop_depths, nrow = n, ncol = ncol(mem), byrow = TRUE)
  weighted  <- mem * depth_mat                 # gate depth where member, else 0
  max_depth <- apply(weighted, 1, max)         # 0 == cell in no candidate gate
  out <- rep(none_label, n)
  has <- max_depth > 0
  at_max <- mem & (weighted == matrix(max_depth, nrow = n, ncol = ncol(mem)))
  cnt <- rowSums(at_max)
  out[has & cnt > 1L] <- multiple_label
  singles <- which(has & cnt == 1L)
  if (length(singles)) {
    col_idx <- max.col(at_max[singles, , drop = FALSE], ties.method = "first")
    out[singles] <- subpop_labels[col_idx]
  }
  out
}

# -----------------------------------------------------------------------------
# (cells-in-`population` x candidate-gate) logical membership for one sample,
# row-aligned to gh_pop_get_data(gh, population). Columns are aligned to the
# `subpop_paths` order; gates absent from THIS sample's tree stay all-FALSE.
#
# An NA entry in `subpop_paths` means the gate belongs to a different workspace
# (see merge_subgate_candidates) and is likewise left all-FALSE.
# -----------------------------------------------------------------------------
get_subpop_membership <- function(gh, viable_idx, subpop_paths) {
  nv <- sum(viable_idx)
  M  <- matrix(FALSE, nrow = nv, ncol = length(subpop_paths))
  if (length(subpop_paths) == 0L || nv == 0L || all(is.na(subpop_paths))) {
    return(M)
  }
  gh_nodes <- flowWorkspace::gh_get_pop_paths(gh, path = "full")
  present  <- !is.na(subpop_paths) & subpop_paths %in% gh_nodes
  if (any(present)) {
    mat <- flowWorkspace::gh_pop_get_indices_mat(gh, y = subpop_paths[present])
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = sum(present))
    mat <- mat[viable_idx, , drop = FALSE]     # root order -> population order
    M[, present] <- mat
  }
  M
}

# -----------------------------------------------------------------------------
# Decide which columns of one flowFrame are biological markers.
#
# We use explore_data(ff)$markers (per-column marker names, "" when none) for
# consistency with the rest of the repo. A channel is kept as a biological
# marker only when it has a non-empty marker name AND its CHANNEL name does not
# match the scatter/Time exclusion regex (applied to the CHANNEL name only, so a
# legitimate antibody whose marker string contains "FSC"/"Time" is not dropped).
#
# Returns list(keep = logical per column, markers = per-column marker names,
# dropped = human-readable "channel [marker]" strings for the excluded columns).
# Shared by extract_population_sample() and gs_common_markers() so the marker set
# reported by list_populations mode can never drift from the one a run uses.
# -----------------------------------------------------------------------------
select_marker_columns <- function(channels, markers,
                                  exclude_pattern = "Time|FSC|SSC|scatter") {
  keep <- vapply(seq_along(channels), function(i) {
    m <- markers[i]
    if (is.na(m) || m == "") return(FALSE)                    # empty marker
    if (grepl(exclude_pattern, channels[i],
              ignore.case = TRUE)) return(FALSE)              # scatter/Time channel
    TRUE
  }, logical(1))

  dropped <- if (any(!keep)) {
    vapply(which(!keep), function(i) {
      sprintf("%s [%s]", channels[i],
              if (is.na(markers[i]) || markers[i] == "") "<no marker>" else markers[i])
    }, character(1))
  } else {
    character(0)
  }

  list(keep = keep, markers = markers, dropped = dropped)
}

# The flowFrame/cytoframe form. Reads parameters() rather than going through
# explore_data(), which materialises exprs() once per column — the same answer
# for a fraction of the work, and it keeps this and gs_common_markers() on one
# code path so the panel reported by list_populations cannot drift from the one
# a run uses.
marker_columns <- function(ff, exclude_pattern = "Time|FSC|SSC|scatter") {
  cm <- cf_channel_markers(ff)
  select_marker_columns(cm$channels, cm$markers, exclude_pattern)
}

# Channel and marker names of a cytoframe WITHOUT touching a single event.
#
# parameters() is metadata: $name is the channel, $desc the marker, already
# aligned and one row per column. The obvious-looking `cf[0, ]` is not a way to
# get here — on flowWorkspace 4.14 a zero-row cytoframe subset builds an index
# of UINT_MAX and dies with "The size of the new row index
# (4294967295,4294967295) is not within the original mat size". Reading the
# metadata directly is both correct and cheaper than materialising exprs().
cf_channel_markers <- function(cf) {
  pd <- flowCore::pData(flowCore::parameters(cf))
  list(channels = as.character(pd$name), markers = as.character(pd$desc))
}

# -----------------------------------------------------------------------------
# Marker names common to ALL samples of a GatingSet, evaluated at the ROOT node.
# Markers are a property of the channels, not of a gate, so this is the same set
# a run would end up with — it just does not require a population to be chosen
# yet. Order follows the first sample (deterministic), matching the run path.
# Used by list_populations mode so a caller (e.g. the Dash app) can offer a
# marker picker from a single workspace load.
# -----------------------------------------------------------------------------
gs_common_markers <- function(gs, sample_names,
                              exclude_pattern = "Time|FSC|SSC|scatter") {
  per_sample <- lapply(sample_names, function(s) {
    # Metadata only: no rows are read, and no column is subset. At the root node
    # a sample is the whole FCS file, so anything that materialises exprs() here
    # is both wasteful and — since flowWorkspace 4.14 — a crash. See
    # cf_channel_markers().
    cm  <- cf_channel_markers(flowWorkspace::gh_pop_get_data(gs[[s]], y = "root"))
    sel <- select_marker_columns(cm$channels, cm$markers, exclude_pattern)
    sel$markers[sel$keep]
  })
  if (length(per_sample) == 0) return(character(0))
  common <- Reduce(intersect, per_sample)
  # Preserve first-sample ordering for determinism (as the run path does).
  per_sample[[1]][per_sample[[1]] %in% common]
}

# -----------------------------------------------------------------------------
# Extract one sample: returns list(mtx = cells x markers, subpop = per-cell
# subgate label). Both are row-aligned to gh_pop_get_data(gh, population); the
# subpop vector is built from the SAME population indices so it stays aligned
# through downsampling and rbind. Kept columns are renamed to their MARKER name.
# -----------------------------------------------------------------------------
extract_population_sample <- function(gh, population, subpop_paths, subpop_depths,
                                      subpop_labels, assignment = "depth",
                                      exclude_pattern = "Time|FSC|SSC|scatter") {
  # NOTE: data is already compensated + transformed in the FlowJo working space.
  # We take exprs() AS-IS, no further transform.
  ff   <- flowWorkspace::gh_pop_get_data(gh, y = population)
  expr <- flowCore::exprs(ff)

  sel     <- marker_columns(ff, exclude_pattern)
  keep    <- sel$keep
  markers <- sel$markers

  # Log excluded channels (channel + marker) so a wrongly-dropped biological
  # channel is visible to the user.
  if (length(sel$dropped) > 0) {
    cat("    excluded channels:", paste(sel$dropped, collapse = ", "), "\n")
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

  # Per-cell subpopulation label, aligned to the SAME `population` cells as mtx.
  # viable_idx is root-relative; gh_pop_get_data() returns those events in root
  # order, so the membership matrix (also subset by viable_idx) is row-aligned.
  viable_idx <- flowWorkspace::gh_pop_get_indices(gh, population)
  stopifnot(sum(viable_idx) == nrow(mtx))   # length/order alignment guard
  subpop <- assign_subpop(
    get_subpop_membership(gh, viable_idx, subpop_paths),
    subpop_depths, subpop_labels, assignment = assignment)
  stopifnot(length(subpop) == nrow(mtx))

  list(mtx = mtx, subpop = subpop)
}

# -----------------------------------------------------------------------------
# Per-marker normalization computed FROM the supplied matrix itself (the pooled
# cells of the SELECTED population). Constant / zero-variance markers are handled
# explicitly so normalization can never introduce NaN/Inf.
# -----------------------------------------------------------------------------
normalize_columns <- function(m, method) {
  if (method == "none") return(m)
  for (j in seq_len(ncol(m))) {
    x <- m[, j]
    if (method == "zscore") {
      mu <- mean(x); s <- stats::sd(x)
      m[, j] <- if (is.finite(s) && s > 0) (x - mu) / s else x - mu
    } else if (method == "minmax") {
      lo <- min(x); hi <- max(x); span <- hi - lo
      m[, j] <- if (is.finite(span) && span > 0) (x - lo) / span else rep(0, length(x))
    }
  }
  m
}

# Strip the FCS extension and GatingSet/cytoset suffix so legends/filenames show
# a readable sample name (e.g. "...AML Pool_001.fcs_0" -> "...AML Pool_001").
clean_sample_label <- function(x) {
  x <- sub("\\.fcs_[0-9]+$", "", x, ignore.case = TRUE)  # cytoset suffix
  x <- sub("\\.fcs$",        "", x, ignore.case = TRUE)  # bare extension
  x
}

# Colour-scale limits for one marker = its full min-max range. Returns NULL
# (ggplot auto-scales) for a degenerate/constant vector so scale_color_viridis
# is never handed c(Inf,-Inf) or a zero-width range.
color_limits <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NULL)
  lims <- range(x)
  if (!all(is.finite(lims)) || lims[1] >= lims[2]) return(NULL)
  lims
}

# Extra figure width (inches) to reserve for a categorical legend so its keys and
# labels get dedicated space instead of compressing the (fixed-aspect) panel.
# Scales with the longest label — this is what stops long file names from
# squishing the plot.
legend_width_in <- function(strings) {
  longest <- max(nchar(as.character(strings)), 0L)
  0.9 + 0.085 * longest          # key glyph + padding + ~per-character width
}
# Number of legend columns so all keys fit within the FIXED plot height instead
# of overflowing a single tall column. Crucially we wrap sideways (more columns)
# rather than growing the figure height, so the data panel — which is height-
# constrained by coord_fixed — stays exactly the same size as every other plot;
# only the legend spreads into extra columns. ~0.3 in per legend row, reserving
# ~1 in for the legend title / margins.
legend_ncol <- function(n_keys, plot_height) {
  per_col <- max(1L, floor((plot_height - 1) / 0.3))
  max(1L, as.integer(ceiling(n_keys / per_col)))
}

# Maximally-distinct categorical colours for MANY (50+) labels, dependency-free.
# ggplot's default hue_pal() only rotates hue at fixed lightness/chroma, so with
# dozens of categories the hues are a few degrees apart and become impossible to
# tell apart. Instead we (1) spread hue by the GOLDEN ANGLE (137.508 deg) so
# consecutive labels land far apart on the colour wheel, and (2) cycle through
# several LIGHTNESS and CHROMA levels so even labels that eventually share a
# similar hue still differ in brightness/saturation. The result keeps adjacent
# legend entries visually separable far past the point where hue alone fails.
distinct_palette <- function(n) {
  if (n <= 0) return(character(0))
  i <- seq_len(n) - 1L
  h <- (i * 137.508) %% 360                       # golden-angle hue spread
  l_levels <- c(40, 58, 74)                       # dark / mid / light
  c_levels <- c(82, 55)                           # vivid / muted
  l  <- l_levels[(i %% length(l_levels)) + 1L]
  cc <- c_levels[((i %/% length(l_levels)) %% length(c_levels)) + 1L]
  grDevices::hcl(h = h, c = cc, l = l)
}

# =============================================================================
# Main entry point: compute and plot a shared UMAP over one FlowJo population,
# across ONE OR MANY FlowJo workspaces.
#
# With several workspaces every cell still goes into ONE embedding, so the
# coordinates are comparable across workspaces; workspace simply becomes another
# grouping level alongside sample, with its own categorical plot. Note that
# pooling workspaces from different acquisition batches will show batch effects
# — nothing here corrects for them (see src/harmonization/ for that); the
# by-workspace plot is what makes them visible.
#
# @param wsp_path A .wsp file, or a DIRECTORY searched recursively for .wsp
#   files. Each workspace is paired with its FCS files by file name.
# @param fcs_dir Directory holding the FCS files, searched recursively.
# @param population Population path to extract (leading '/' optional; a suffix
#   such as 'Single Cells/Viable cells' is matched against the gating tree).
#   Resolved independently in each workspace, so workspaces whose gating trees
#   differ above the population still line up.
# @param output_dir Output directory for plots.
# @param group FlowJo sample group to load: an integer index, or a group NAME
#   (which is what lets one value cover workspaces that order groups
#   differently).
# @param expected_samples If > 0, warn when the run does not total this many
#   samples across all workspaces.
# @param sample_metadata Optional CSV keyed on sample (a `sample_id` column, or
#   the first column) whose other columns become per-sample annotations —
#   patient, timepoint, response, anything. Each requested column yields one
#   plot coloured by it plus one plot per level, all on the SAME embedding.
# @param sample_pattern Optional regex with NAMED capture groups applied to the
#   sample name when no metadata file is given, e.g.
#   "(?<patient>MTB[0-9]+)_(?<timepoint>TP[0-9]+)".
# @param annotate Optional character vector naming which annotation columns to
#   plot (default: all of them).
# @param sample_id_column Optional name of the metadata column holding the
#   sample identifier (default: `sample_id`, else the first column).
# @param markers Optional character vector subset of marker names (NULL = all
#   common non-scatter markers).
# @param subpop_mode One of 'cd4cd8', 'direct', 'deepest', 'none'.
# @param subpop_select Optional character vector of gates (relative path or bare
#   name) used as the colour partition, overriding subpop_mode.
# @param normalize One of 'none', 'zscore', 'minmax'.
# @param downsample Max cells per sample AFTER gating (0 = no downsampling).
# @param n_neighbors,min_dist,metric,seed UMAP parameters / RNG seed.
# @param point_size,point_alpha ggplot point aesthetics.
# @param plot_width,plot_height Base plot dimensions (inches).
# @param output_formats Named logical list, e.g. list(pdf=,png=,svg=).
# @param list_populations If TRUE, print the gating tree and return it.
# @return Invisibly, a list with umap_coords, umap_df, output_dir (or, in
#   list_populations mode, the character vector of population paths).
# =============================================================================
umap_from_flowjo <- function(wsp_path, fcs_dir,
                             population = "Lymphocytes/Single Cells/Viable cells",
                             output_dir = "umap_flowjo",
                             group = "1", expected_samples = 0L,
                             sample_metadata = NULL, sample_pattern = NULL,
                             annotate = NULL, sample_id_column = NULL,
                             markers = NULL, subpop_mode = "furthest",
                             subpop_select = NULL, normalize = "none",
                             downsample = 0L, n_neighbors = 15L, min_dist = 0.1,
                             metric = "euclidean", seed = 42L,
                             point_size = 0.4, point_alpha = 0.6,
                             plot_width = 10, plot_height = 8,
                             output_formats = list(pdf = FALSE, png = FALSE, svg = TRUE),
                             list_populations = FALSE,
                             allow_missing_samples = FALSE,
                             population_map = NULL) {
  require(CytoML)
  require(flowWorkspace)
  require(flowCore)
  require(ggplot2)
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for marker colour scales.")
  }

  VALID_NORMALIZE <- c("none", "zscore", "minmax")
  if (!normalize %in% VALID_NORMALIZE) {
    stop("normalize must be one of: ", paste(VALID_NORMALIZE, collapse = ", "),
         " (got '", normalize, "').")
  }
  VALID_SUBPOP_MODE <- c("furthest", "cd4cd8", "direct", "deepest", "none")
  if (!subpop_mode %in% VALID_SUBPOP_MODE) {
    stop("subpop_mode must be one of: ",
         paste(VALID_SUBPOP_MODE, collapse = ", "),
         " (got '", subpop_mode, "').")
  }

  group <- as.character(group)
  if (!dir.exists(fcs_dir)) stop("fcs_dir not found: ", fcs_dir)
  fcs_dir <- normalizePath(fcs_dir, mustWork = TRUE)

  # --- 1. Discover the workspaces and pair each with its FCS files. ----------
  # Everything is resolved BEFORE any GatingSet is built, so a mislaid file or a
  # bad --group is reported for every workspace at once instead of one slow load
  # at a time. See src/loading/flowjo_workspace.R.
  wsp_paths <- find_workspace_files(wsp_path)
  wsp_names <- workspace_labels(wsp_paths)
  n_wsp     <- length(wsp_paths)

  # Names, not full paths: they share a root, and 33 repetitions of it push the
  # part that matters — what failed — off the top of the log.
  cat("Workspaces found (", n_wsp, ") under ", wsp_path, ":\n", sep = "")
  cat(strwrap(paste(wsp_names, collapse = ", "), width = 78,
              prefix = "  ", initial = "  "), sep = "\n")
  fcs_index <- index_fcs_dir(fcs_dir)
  cat("FCS files indexed: ", length(fcs_index), " distinct name(s) under ",
      fcs_dir, "\n", sep = "")

  cat("Resolving ", n_wsp, " workspace(s) (group = \"", group, "\") ...\n",
      sep = "")
  resolved <- resolve_workspace_specs(wsp_paths, wsp_names, fcs_index, group,
                                      allow_missing = allow_missing_samples)
  report   <- format_resolution_report(resolved$problems, fcs_root = fcs_dir)
  if (nzchar(report)) {
    # The report is written out BEFORE stop(), not passed to it. R truncates the
    # error message it PRINTS at getOption("warning.length") (1000 chars by
    # default), which silently swallowed most of a 33-workspace report — the
    # whole point of collecting every problem first. stop() then carries only a
    # short summary, which is also what a caller sees in conditionMessage().
    cat(report, "\n", sep = "", file = stderr())
    stop("could not build the workspace set: ", length(resolved$problems),
         " of ", n_wsp, " workspace(s) unusable. See the report above.",
         call. = FALSE)
  }
  specs <- resolved$specs

  # What --allow_missing_samples discarded. Said before the slow part, and again
  # in the run summary at the end, so it cannot be lost in the scrollback of a
  # long load — the embedding is short of cells and that has to stay visible.
  dropped_report <- format_dropped_report(resolved$dropped)
  if (nzchar(dropped_report)) {
    cat("\nWARNING: ", dropped_report, "\n", sep = "")
  }
  if (!length(specs)) {
    stop("every workspace was dropped: not one FCS file named by any of the ",
         n_wsp, " workspace(s) was found under ", fcs_dir,
         ". This is a wrong --fcs_dir, not a gap in the data.", call. = FALSE)
  }

  # What DID work, stated as plainly as what did not. Without this the log only
  # ever names failures, and a reader cannot tell a run that lost two controls
  # from one that lost half the study.
  n_samples <- sum(vapply(specs, function(s) nrow(s$files), integer(1)))
  cat(sprintf("\nGOING INTO THE UMAP: %d of %d workspace(s), %d sample(s).\n",
              length(specs), n_wsp, n_samples))

  # Overlapping workspaces double-count cells rather than failing, so this is
  # said loudly and early, before the slow part.
  shared_report <- format_shared_file_report(specs)
  if (nzchar(shared_report)) {
    cat("\nWARNING: ", shared_report, "\n", sep = "")
  }

  # --- 1b. Build one GatingSet per workspace. --------------------------------
  # A build can still fail after resolution succeeded — CytoML checks the event
  # count of each file it opens, so a name we resolved can turn out to be a
  # different acquisition. With --allow_missing_samples that costs the run one
  # workspace, not all of them: the flag says "load what is loadable", and a
  # workspace CytoML refuses is a workspace whose data is unusable here.
  gs_list <- list()
  build_failed <- list()
  for (label in names(specs)) {
    spec <- specs[[label]]
    cat(sprintf("\nBuilding GatingSet for '%s' — %s, %d sample(s)\n",
                label, describe_group(spec), nrow(spec$files)))
    gs <- tryCatch(build_workspace_gatingset(spec),
                   error = function(e) {
                     if (!allow_missing_samples) stop(e)
                     build_failed[[label]] <<- conditionMessage(e)
                     NULL
                   })
    if (!is.null(gs)) gs_list[[label]] <- gs
  }
  if (length(build_failed)) {
    cat("\nWARNING: ", length(build_failed),
        " workspace(s) resolved their files but could NOT be loaded, and are\n",
        "not in the analysis. Each reason is CytoML's own:\n", sep = "")
    for (label in names(build_failed)) {
      cat(sprintf("  %s\n    %s\n", label,
                  gsub("\n", "\n    ", build_failed[[label]])))
    }
  }
  if (!length(gs_list)) {
    stop("no workspace could be loaded; see the failures above.", call. = FALSE)
  }

  # A population map is also a workspace SELECTION: a workspace with no entry
  # contributes no gate and therefore no cells. Applied here, before anything
  # downstream is derived, so sample labels, annotations and the shared-file
  # accounting all describe the run that is actually happening.
  pop_map <- parse_population_map(population_map, names(gs_list))
  if (!is.null(pop_map) && !isTRUE(list_populations)) {
    dropped_by_map <- setdiff(names(gs_list), names(pop_map))
    gs_list <- gs_list[names(gs_list) %in% names(pop_map)]
    specs   <- specs[names(specs) %in% names(pop_map)]
    if (!length(gs_list)) {
      stop("the population map excludes every workspace that loaded.",
           call. = FALSE)
    }
    if (length(dropped_by_map)) {
      cat(sprintf(
        "\n%d workspace(s) are NOT in this UMAP because the population map gives\nthem no gate: %s\n",
        length(dropped_by_map), paste(dropped_by_map, collapse = ", ")))
    }
  }

  # Sample display labels across the WHOLE run. Names are kept exactly as they
  # were unless two workspaces genuinely collide, so a single-workspace run
  # produces the same legends and the same file names as it always did.
  sample_index <- do.call(rbind, lapply(names(gs_list), function(w) {
    raw <- flowWorkspace::sampleNames(gs_list[[w]])
    data.frame(workspace = rep(w, length(raw)), raw = raw,
               clean = clean_sample_label(raw), stringsAsFactors = FALSE)
  }))
  sample_index$display <- disambiguate_sample_labels(sample_index$clean,
                                                     sample_index$workspace)
  if (anyDuplicated(sample_index$display)) {
    warning("Sample labels still collide after cleaning; falling back to the ",
            "raw sample names for display.")
    sample_index$display <- paste0(sample_index$workspace, ": ",
                                   sample_index$raw)
  }

  # Per-sample annotations (patient, timepoint, ...). Resolved here, once, so a
  # bad column name or an unmatched metadata table is reported before the slow
  # part rather than after it.
  sample_annotations <- resolve_sample_annotations(
    sample_index, metadata_file = sample_metadata, pattern = sample_pattern,
    annotate = annotate, id_column = sample_id_column)
  annot_cols <- colnames(sample_annotations)
  if (length(annot_cols)) {
    cat("\nSample annotations (", length(annot_cols), "): ",
        paste(annot_cols, collapse = ", "), "\n", sep = "")
    for (cl in annot_cols) {
      lv <- unique(sample_annotations[[cl]])
      cat(sprintf("  %s: %d level(s) — %s\n", cl, length(lv),
                  paste(utils::head(lv, 8), collapse = ", ")))
    }
  }

  cat("\nSamples in loaded group (", nrow(sample_index), "):\n", sep = "")
  for (s in sample_index$display) cat("  -", s, "\n")
  if (expected_samples > 0 && nrow(sample_index) != expected_samples) {
    warning(sprintf(
      paste0("The %d loaded workspace(s) total %d sample(s) but ",
             "expected_samples=%d. The wrong sample group may be selected; ",
             "re-check group or list_populations."),
      n_wsp, nrow(sample_index), expected_samples))
  }

  # Gating tree of each workspace. They can differ, so they are kept apart and
  # only reconciled through the RELATIVE paths under the selected population.
  pops_by_wsp <- lapply(gs_list, function(gs)
    flowWorkspace::gs_get_pop_paths(gs, path = "full"))
  # Gates present in EVERY workspace, in the first workspace's tree order. This
  # is the set a run can actually use: --population must resolve everywhere.
  shared_pops <- Reduce(intersect, pops_by_wsp)
  shared_pops <- pops_by_wsp[[1]][pops_by_wsp[[1]] %in% shared_pops]

  # --- 2. list_populations mode: print the tree + markers, then return. ------
  # This is the "inspect" entry point: one workspace load reports everything a
  # caller needs to configure a run (samples above, populations and the common
  # marker panel here). The block headers below are a stable, parseable contract
  # relied on by app/ (the Dash front-end) — change them and update that too.
  # The two extra blocks at the end are indented FOUR spaces, so the existing
  # two-space item patterns cannot swallow them.
  if (isTRUE(list_populations)) {
    cat("\nPopulation tree (", length(shared_pops), " nodes):\n", sep = "")
    for (p in shared_pops) cat("  ", p, "\n", sep = "")
    if (length(shared_pops) == 0) {
      cat("  (no gate is shared by all workspaces — see the divergence list ",
          "below)\n", sep = "")
    }

    # Markers are a channel property, so they are evaluated at the root node and
    # do not depend on which population is eventually selected. Intersected
    # across every workspace as well as every sample.
    markers_by_wsp <- lapply(names(gs_list), function(w)
      gs_common_markers(gs_list[[w]], flowWorkspace::sampleNames(gs_list[[w]])))
    common_markers <- Reduce(intersect, markers_by_wsp)
    common_markers <- markers_by_wsp[[1]][markers_by_wsp[[1]] %in% common_markers]
    cat("\nCommon markers (", length(common_markers), "):\n", sep = "")
    for (m in common_markers) cat("  ", m, "\n", sep = "")
    if (length(common_markers) == 0) {
      cat("  (no marker common to all samples — check the panel / group)\n")
    }

    cat("\nSample annotations (", length(annot_cols), "):\n", sep = "")
    for (cl in annot_cols) {
      lv <- unique(sample_annotations[[cl]])
      cat(sprintf("    %s | %d level(s) | %s\n", cl, length(lv),
                  paste(lv, collapse = ", ")))
    }
    if (length(annot_cols) == 0) {
      cat("    (none — pass sample_metadata or sample_pattern to colour by ",
          "patient, timepoint, ...)\n", sep = "")
    }

    # Every workspace's OWN tree, which is what a per-workspace population
    # choice needs: the shared block above is an intersection and collapses to
    # almost nothing when the workspaces were gated independently. Six-space
    # indent for the gates so the two- and four-space item patterns in
    # app/services/rscript.py cannot swallow them; part of the same contract.
    cat("\nPopulations by workspace (", length(pops_by_wsp), "):\n", sep = "")
    for (w in names(pops_by_wsp)) {
      cat(sprintf("    %s | %d node(s)\n", w, length(pops_by_wsp[[w]])))
      for (p in pops_by_wsp[[w]]) cat("      ", p, "\n", sep = "")
    }

    cat("\nWorkspaces resolved (", length(specs), "):\n", sep = "")
    for (label in names(specs)) {
      spec <- specs[[label]]
      cat(sprintf("    %s | %s | %s | %d sample(s)\n",
                  label, spec$wsp, describe_group(spec), nrow(spec$files)))
    }

    divergent <- setdiff(unique(unlist(pops_by_wsp, use.names = FALSE)),
                         shared_pops)
    cat("\nPopulations not shared by all workspaces (", length(divergent),
        "):\n", sep = "")
    for (p in divergent) {
      absent <- names(pops_by_wsp)[!vapply(pops_by_wsp,
                                           function(x) p %in% x, logical(1))]
      cat(sprintf("    %s   [missing in: %s]\n", p,
                  paste(absent, collapse = ", ")))
    }
    if (length(divergent) == 0 && n_wsp > 1) {
      cat("    (none — every workspace has the same gating tree)\n")
    }

    cat("\nSelect one with population = (leading '/' optional), e.g.:\n")
    if (length(shared_pops)) {
      cat("  population = \"", sub("^/", "", shared_pops[length(shared_pops)]),
          "\"\n", sep = "")
    }
    return(invisible(list(populations = shared_pops,
                          populations_by_workspace = pops_by_wsp,
                          samples = sample_index$display,
                          markers = common_markers,
                          workspaces = names(specs),
                          annotations = sample_annotations)))
  }

  if (is.null(pop_map) && (is.null(population) || !nzchar(population))) {
    stop("population is required (or use list_populations = TRUE to see options).")
  }

  # --- 2a. Resolve the population in EVERY workspace. ------------------------
  # Resolved independently per workspace (the trees can differ above the
  # population). With ONE --population every failure is reported together and is
  # fatal: pooling cells from gates that are not the same gate would be a silent
  # scientific error, and a request that cannot be found somewhere may well have
  # landed on the wrong gate where it was found.
  #
  # With a population MAP the pairing is no longer inferred from a shared name —
  # it is stated per workspace — so differently-named gates pool by design. A
  # mapped gate that does not resolve is still fatal, because it is an explicit
  # instruction that could not be carried out.
  if (!is.null(pop_map)) {
    resolutions <- lapply(names(gs_list), function(w)
      try_resolve_population(pop_map[[w]], pops_by_wsp[[w]]))
    names(resolutions) <- names(gs_list)
  } else {
    resolutions <- lapply(pops_by_wsp, function(p)
      try_resolve_population(population, p))
  }
  failed <- names(resolutions)[!vapply(resolutions, function(r) r$ok, logical(1))]
  if (length(failed)) {
    detail <- vapply(failed, function(w) sprintf(
      "  %s (%s)\n    %s", w, specs[[w]]$wsp, resolutions[[w]]$reason),
      character(1))
    if (!is.null(pop_map)) {
      stop("the population map names a gate that does not exist in its ",
           "workspace.\n\n", paste(detail, collapse = "\n"),
           "\n\nRun with list_populations = TRUE to see each workspace's own ",
           "gates, and correct the entry (or remove that workspace).",
           call. = FALSE)
    }
    stop("--population could not be resolved in every workspace.\n\n",
         paste(detail, collapse = "\n"),
         "\n\nRun with list_populations = TRUE to see the gates shared by all ",
         "workspaces, or select a gate per workspace with --population_map.")
  }
  pop_by_wsp <- vapply(resolutions, function(r) r$path, character(1))
  population <- pop_by_wsp[[1]]   # display / file-name form
  if (!is.null(pop_map)) {
    # Named per workspace, so the pairing being pooled is on the record next to
    # the plots rather than only in whatever built the command.
    cat("\nPopulation per workspace (", length(pop_by_wsp),
        " workspace(s), paired by your choice):\n", sep = "")
    for (w in names(pop_by_wsp)) cat("     ", w, ": ", pop_by_wsp[[w]], "\n",
                                     sep = "")
    cat("  Plots and file names use '", population, "'.\n", sep = "")
  } else {
    cat("\nTarget population resolved to:", population, "\n")
    if (n_wsp > 1 && length(unique(pop_by_wsp)) > 1) {
      cat("  (it resolves to a different absolute path in some workspaces:\n")
      for (w in names(pop_by_wsp)) cat("     ", w, ": ", pop_by_wsp[[w]], "\n",
                                       sep = "")
      cat("   they are treated as the same population.)\n")
    }
  }

  # --- 2b. Discover the FlowJo subgates UNDER `population` to colour by. ------
  # Per workspace, then merged on the RELATIVE path (the only key comparable
  # across workspaces whose ancestors are named differently).
  cands <- lapply(names(gs_list), function(w) {
    gs <- gs_list[[w]]
    pop_paths_union <- unique(unlist(lapply(
      flowWorkspace::sampleNames(gs),
      function(s) flowWorkspace::gh_get_pop_paths(gs[[s]], path = "full"))))
    flowjo_subgate_candidates(pop_paths_union, pops_by_wsp[[w]],
                              pop_by_wsp[[w]], mode = subpop_mode,
                              select = subpop_select, warn = FALSE)
  })
  names(cands) <- names(gs_list)
  cand <- merge_subgate_candidates(cands)

  subpop_labels <- cand$labels
  subpop_depths <- cand$depths
  subpop_paths_by_wsp <- cand$paths
  subpop_mode_tag <- cand$tag

  # A --subpop_select entry is only genuinely missing when NO workspace has it.
  if (!is.null(subpop_select)) {
    matched <- unique(c(subpop_labels, basename(subpop_labels)))
    missing_sel <- setdiff(subpop_select, matched)
    if (length(missing_sel)) {
      warning("subpop_select entries not found under ", population,
              " in any workspace: ", paste(missing_sel, collapse = ", "))
    }
  }

  # 'furthest' labels each cell by the most specific (smallest) gate it reaches;
  # the other modes assign by path depth (ties -> "(multiple)").
  subpop_assignment <- if (subpop_mode == "furthest" && is.null(subpop_select)) {
    "specificity"
  } else {
    "depth"
  }

  do_subpop <- length(subpop_labels) > 0
  subpop_requested <- !is.null(subpop_select) || subpop_mode != "none"
  if (subpop_requested) {
    if (do_subpop) {
      cat(sprintf("\nSubpopulation colouring (%s): %d candidate subgate(s) under %s\n",
                  subpop_mode_tag, length(subpop_labels), population))
    } else {
      cat(sprintf(
        "\nSubpopulation colouring requested (%s) but no matching subgates exist under %s; skipping.\n",
        subpop_mode_tag, population))
    }
  }

  # --- 3 & 4. Extract the population per sample (matrix + subgate labels). ----
  set.seed(seed)  # governs downsampling; run_umap re-seeds internally for UMAP

  per_sample_mtx    <- list()
  per_sample_subpop <- list()
  per_sample_wsp    <- character(0)
  for (w in names(gs_list)) {
    gs <- gs_list[[w]]
    if (n_wsp > 1) cat(sprintf("\n[%s]\n", w))
    for (s in flowWorkspace::sampleNames(gs)) {
      # The display label is the key everywhere downstream (legends, file
      # names, per-file subsetting), so switch to it here and never look back.
      display <- sample_index$display[sample_index$workspace == w &
                                        sample_index$raw == s]
      res <- extract_population_sample(gs[[s]], pop_by_wsp[[w]],
                                       subpop_paths_by_wsp[[w]],
                                       subpop_depths, subpop_labels,
                                       assignment = subpop_assignment)
      mtx <- res$mtx
      sp  <- res$subpop

      # Optional per-sample downsampling AFTER gating. The SAME row index
      # subsets the subpopulation labels so they stay aligned to mtx.
      if (downsample > 0 && nrow(mtx) > downsample) {
        idx <- sort(sample(nrow(mtx), downsample))
        mtx <- mtx[idx, , drop = FALSE]
        sp  <- sp[idx]
      }

      cat(sprintf("  %s: %d cells x %d markers\n", display, nrow(mtx),
                  ncol(mtx)))
      if (nrow(mtx) == 0) {
        warning("Sample '", display, "' has 0 cells in population '",
                pop_by_wsp[[w]],
                "' after gating/downsampling; it will not contribute rows.")
      }
      per_sample_mtx[[display]]    <- mtx
      per_sample_subpop[[display]] <- sp
      per_sample_wsp[display]      <- w
    }
  }

  # Drop samples that contributed no cells (keeps alignment + plotting sane).
  nonempty <- vapply(per_sample_mtx, nrow, integer(1)) > 0
  if (!any(nonempty)) {
    stop("No sample has any cells in population '", population, "'.")
  }
  per_sample_mtx    <- per_sample_mtx[nonempty]
  per_sample_subpop <- per_sample_subpop[nonempty]
  per_sample_wsp    <- per_sample_wsp[names(per_sample_mtx)]

  # Intersect markers common to ALL (non-empty) samples.
  common_markers <- Reduce(intersect, lapply(per_sample_mtx, colnames))
  if (length(common_markers) == 0) {
    stop("No markers are common to all samples.")
  }
  # Preserve first-sample ordering for determinism.
  common_markers <- colnames(per_sample_mtx[[1]])[
    colnames(per_sample_mtx[[1]]) %in% common_markers]

  # Apply optional marker subset (must be a subset of common markers).
  if (!is.null(markers)) {
    missing_m <- setdiff(markers, common_markers)
    if (length(missing_m) > 0) {
      stop("These markers are not common markers in the data: ",
           paste(missing_m, collapse = ", "),
           "\nAvailable common markers: ",
           paste(common_markers, collapse = ", "))
    }
    use_markers <- markers
  } else {
    use_markers <- common_markers
  }

  cat("\nMarkers used (", length(use_markers), "):\n  ",
      paste(use_markers, collapse = ", "), "\n", sep = "")

  # Restrict every per-sample matrix to use_markers (same column order).
  per_sample_mtx <- lapply(per_sample_mtx, function(m) m[, use_markers, drop = FALSE])

  # --- 5. Concatenate. file_labels / combined_subpop align to matrix rows. ---
  # Names are already the run-wide display labels (see disambiguate_sample_
  # labels above), so no further cleaning happens here.
  sample_labels <- names(per_sample_mtx)
  n_cells_per   <- vapply(per_sample_mtx, nrow, integer(1))

  combined_matrix <- do.call(rbind, per_sample_mtx)
  file_labels     <- rep(sample_labels, n_cells_per)
  stopifnot(length(file_labels) == nrow(combined_matrix))

  # Workspace of each cell, built from the SAME per-sample order and counts, so
  # it is row-aligned to combined_matrix exactly as file_labels is.
  wsp_labels_per_cell <- rep(unname(per_sample_wsp[sample_labels]), n_cells_per)
  stopifnot(length(wsp_labels_per_cell) == nrow(combined_matrix))

  # Subpopulation labels, concatenated in the SAME sample order as
  # combined_matrix (the two lists are filtered/ordered identically).
  combined_subpop <- unlist(per_sample_subpop[sample_labels], use.names = FALSE)
  stopifnot(length(combined_subpop) == nrow(combined_matrix))

  cat(sprintf("\nConcatenated matrix: %d cells x %d markers (%d samples, %d workspace(s))\n",
              nrow(combined_matrix), ncol(combined_matrix),
              length(per_sample_mtx),
              length(unique(wsp_labels_per_cell))))

  # --- 5b. Optional per-marker normalization. --------------------------------
  if (normalize != "none") {
    cat(sprintf(
      "Normalizing %d markers per-column (method=%s) from the concatenated gated cells...\n",
      ncol(combined_matrix), normalize))
    combined_matrix <- normalize_columns(combined_matrix, normalize)
  }

  # --- 6. Compute UMAP ONCE on the concatenated matrix. ----------------------
  # CRITICAL alignment guard: uwot may error or drop rows on non-finite input,
  # which would silently break the row-for-row alignment. Fail loudly here.
  if (!all(is.finite(combined_matrix))) {
    bad_cols <- colnames(combined_matrix)[
      apply(combined_matrix, 2, function(x) any(!is.finite(x)))]
    stop("combined_matrix contains non-finite (NA/NaN/Inf) values in marker(s): ",
         paste(bad_cols, collapse = ", "),
         ". Refusing to run UMAP because dropped rows would break index ",
         "alignment with the per-file plots.")
  }

  cat("Computing shared UMAP (n_neighbors=", n_neighbors,
      ", min_dist=", min_dist, ", metric=", metric,
      ", seed=", seed, ")...\n", sep = "")
  umap_coords <- run_umap(
    combined_matrix,
    n_neighbors  = n_neighbors,
    min_dist     = min_dist,
    n_components = 2,
    metric       = metric,
    seed         = seed
  )
  stopifnot(nrow(umap_coords) == nrow(combined_matrix))

  # Master data frame: one row per cell, aligned to combined_matrix.
  umap_df <- data.frame(
    UMAP1     = umap_coords[, 1],
    UMAP2     = umap_coords[, 2],
    file      = file_labels,
    workspace = wsp_labels_per_cell,
    subpop    = combined_subpop,
    stringsAsFactors = FALSE
  )

  # Annotation columns, expanded from per-SAMPLE to per-CELL through the same
  # per-sample order and row counts as file_labels, so they cannot drift out of
  # alignment with the coordinates. Prefixed to keep an annotation called
  # "file" or "subpop" from colliding with the columns above.
  annot_of <- match(sample_labels, sample_index$display)
  stopifnot(!anyNA(annot_of))
  for (cl in annot_cols) {
    umap_df[[paste0("annot_", cl)]] <-
      rep(sample_annotations[[cl]][annot_of], n_cells_per)
  }
  stopifnot(nrow(umap_df) == nrow(combined_matrix))

  # --- 7. Plotting. ----------------------------------------------------------
  make_outdir(output_dir)
  pt_size  <- point_size
  pt_alpha <- point_alpha

  pop_tag <- create_safe_filename(sub("^/", "", population))

  # Suffix appended to marker colour-scale titles so a reader can tell whether
  # the colour encodes raw or normalized expression.
  norm_tag <- switch(normalize,
                     none   = "",
                     zscore = " (z-score)",
                     minmax = " (min-max)")

  # Shared axis limits so every plot (conjoint + per-file) uses the same frame.
  xlim_shared <- range(umap_df$UMAP1, na.rm = TRUE)
  ylim_shared <- range(umap_df$UMAP2, na.rm = TRUE)

  # Equal-aspect coordinate system shared by EVERY plot. UMAP1 and UMAP2 live in
  # the same embedding space, so one unit on each axis must render at the same
  # size (ratio = 1); otherwise the embedding is geometrically distorted. Because
  # the data aspect is fixed, the panel keeps its true shape no matter how wide
  # the figure or its legend is — a wide legend can no longer squish the panel,
  # it just gets letter-boxed padding.
  umap_coord <- function() {
    ggplot2::coord_fixed(ratio = 1, xlim = xlim_shared, ylim = ylim_shared)
  }
  base_umap_theme <- list(
    ggplot2::theme_minimal(),
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  )

  save_one <- function(p, basename_no_ext) {
    save_plot_formats(
      p,
      file.path(output_dir, basename_no_ext),
      width   = plot_width,
      height  = plot_height,
      formats = output_formats
    )
  }

  # Save a categorical (discrete-legend) plot. The figure HEIGHT stays at
  # plot_height (so the coord_fixed panel matches every other plot exactly) and
  # only the WIDTH grows — by ncol legend columns — so a long, many-label legend
  # spreads sideways instead of squishing or clipping the panel.
  save_categorical <- function(p, basename_no_ext, legend_strings, ncol = 1L) {
    save_plot_formats(
      p,
      file.path(output_dir, basename_no_ext),
      width   = plot_width + ncol * legend_width_in(legend_strings),
      height  = plot_height,
      formats = output_formats
    )
  }

  # 7a. Categorical plot coloured by file/sample (all cells).
  file_levels <- sort(unique(umap_df$file))
  file_ncol <- legend_ncol(length(file_levels), plot_height)
  p_file <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = file)) +
    geom_point(size = pt_size, alpha = pt_alpha) +
    umap_coord() +
    base_umap_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = file_ncol)) +
    labs(title = paste0("UMAP - ", population, " - by file"),
         color = "File")
  save_categorical(p_file, create_safe_filename(paste0("umap_", pop_tag, "_by_file")),
                   legend_strings = c(file_levels, "File"),
                   ncol = file_ncol)
  cat("Saved: by-file categorical plot\n")

  # 7a-ter. Categorical plot coloured by WORKSPACE (all cells, same embedding).
  # Only when more than one workspace contributed — a single-workspace run would
  # otherwise gain a one-colour plot and a file that did not exist before.
  #
  # This is also the plot to read for batch effects: nothing here corrects for
  # them, so workspaces separating into their own islands means the embedding is
  # partly driven by acquisition, not biology.
  wsp_levels <- unique(wsp_labels_per_cell)
  if (length(wsp_levels) > 1) {
    wsp_counts <- table(factor(umap_df$workspace, levels = wsp_levels))
    wsp_lab <- sprintf("%s (n=%d)", wsp_levels,
                       as.integer(wsp_counts[wsp_levels]))
    wsp_cols <- stats::setNames(distinct_palette(length(wsp_levels)), wsp_levels)

    df_ws <- umap_df
    df_ws$workspace <- factor(df_ws$workspace, levels = wsp_levels)

    wsp_ncol <- legend_ncol(length(wsp_levels), plot_height)
    p_ws <- ggplot(df_ws, aes(x = UMAP1, y = UMAP2, color = workspace)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      scale_color_manual(values = wsp_cols, breaks = wsp_levels,
                         labels = wsp_lab, drop = FALSE) +
      umap_coord() +
      base_umap_theme +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                  ncol = wsp_ncol)) +
      labs(title = paste0("UMAP - ", population, " - by workspace"),
           color = "Workspace")
    save_categorical(
      p_ws, create_safe_filename(paste0("umap_", pop_tag, "_by_workspace")),
      legend_strings = c(wsp_lab, "Workspace"),
      ncol = wsp_ncol)
    cat(sprintf("Saved: by-workspace categorical plot (%d workspaces)\n",
                length(wsp_levels)))
  }

  # 7a-bis. Categorical plot coloured by FlowJo subgate (all cells, same
  # embedding). Colours follow the gating tree; the catch-all
  # "(multiple)"/"(none)" levels get greys so real subgates keep the
  # distinguishable hues. Legend labels carry the per-level cell count.
  if (do_subpop) {
    present_real <- intersect(subpop_labels, unique(combined_subpop))
    catchall     <- intersect(c(.SUBPOP_MULTI, .SUBPOP_NONE), unique(combined_subpop))
    sp_levels    <- c(present_real, catchall)

    sp_counts <- table(factor(combined_subpop, levels = sp_levels))
    sp_lab    <- sprintf("%s (n=%d)", sp_levels, as.integer(sp_counts[sp_levels]))

    sp_cols <- character(0)
    if (length(present_real)) {
      sp_cols <- stats::setNames(
        distinct_palette(length(present_real)), present_real)
    }
    if (.SUBPOP_MULTI %in% sp_levels) sp_cols[.SUBPOP_MULTI] <- "grey55"
    if (.SUBPOP_NONE  %in% sp_levels) sp_cols[.SUBPOP_NONE]  <- "grey85"

    df_sp <- umap_df
    df_sp$subpop <- factor(df_sp$subpop, levels = sp_levels)

    sp_ncol <- legend_ncol(length(sp_levels), plot_height)
    p_sp <- ggplot(df_sp, aes(x = UMAP1, y = UMAP2, color = subpop)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      scale_color_manual(values = sp_cols, breaks = sp_levels, labels = sp_lab,
                         drop = FALSE) +
      umap_coord() +
      base_umap_theme +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                  ncol = sp_ncol)) +
      labs(title = paste0("UMAP - ", population, " - by subpopulation (",
                          subpop_mode_tag, ")"),
           color = "Subpopulation")
    save_categorical(
      p_sp, create_safe_filename(paste0("umap_", pop_tag, "_by_subpopulation")),
      legend_strings = c(sp_lab, "Subpopulation"),
      ncol = sp_ncol)
    cat(sprintf("Saved: by-subpopulation categorical plot (%d subgate(s) + %d catch-all level(s))\n",
                length(present_real), length(catchall)))
  }

  # 7a-quater. Sample annotations: for EACH column one plot coloured by it, and
  # then one plot per level with that level drawn over the whole embedding in
  # grey. One embedding, many views — the grey context is the point, because an
  # isolated subset gives you no idea where it sits relative to everything else,
  # which is exactly what a shared embedding is for.
  for (cl in annot_cols) {
    values <- umap_df[[paste0("annot_", cl)]]
    # Real levels first (sorted), the Unassigned catch-all last and in grey, so
    # it can never take a colour a real level needs.
    real   <- sort(setdiff(unique(values), .ANNOT_UNASSIGNED))
    levels_all <- c(real, intersect(.ANNOT_UNASSIGNED, unique(values)))

    counts <- table(factor(values, levels = levels_all))
    lab    <- sprintf("%s (n=%d)", levels_all, as.integer(counts[levels_all]))
    cols   <- stats::setNames(distinct_palette(length(real)), real)
    if (.ANNOT_UNASSIGNED %in% levels_all) {
      cols[.ANNOT_UNASSIGNED] <- "grey70"
    }

    df_a <- umap_df
    df_a$.annot <- factor(values, levels = levels_all)

    a_ncol <- legend_ncol(length(levels_all), plot_height)
    p_a <- ggplot(df_a, aes(x = UMAP1, y = UMAP2, color = .annot)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      scale_color_manual(values = cols, breaks = levels_all, labels = lab,
                         drop = FALSE) +
      umap_coord() +
      base_umap_theme +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                  ncol = a_ncol)) +
      labs(title = paste0("UMAP - ", population, " - by ", cl), color = cl)
    save_categorical(
      p_a,
      create_safe_filename(paste0("umap_", pop_tag, "_annot_", cl, "_ALL")),
      legend_strings = c(lab, cl), ncol = a_ncol)

    # One plot per level: the whole embedding in grey, this level on top.
    for (lv in levels_all) {
      rows <- which(values == lv)
      if (length(rows) == 0) next
      p_lv <- ggplot(mapping = aes(x = UMAP1, y = UMAP2)) +
        geom_point(data = umap_df, color = "grey88", size = pt_size,
                   alpha = pt_alpha) +
        geom_point(data = umap_df[rows, , drop = FALSE],
                   color = unname(cols[lv]), size = pt_size, alpha = pt_alpha) +
        umap_coord() +
        base_umap_theme +
        labs(title = paste0("UMAP - ", population, " - ", cl, ": ", lv),
             subtitle = sprintf("%d of %d cells (grey = all cells)",
                                length(rows), nrow(umap_df)))
      save_one(p_lv, create_safe_filename(
        paste0("umap_", pop_tag, "_annot_", cl, "_level_", lv)))
    }
    cat(sprintf("Saved: %s (1 coloured + %d level plot(s))\n", cl,
                length(levels_all)))
  }

  files_unique <- unique(file_labels)

  # 7b & 7c. For each marker: one conjoint plot (all cells) + one plot per file.
  for (marker in use_markers) {
    # ONE colour scale per marker, shared by the conjoint plot and every per-file
    # plot, so a given colour means the same expression value in every plot of
    # this marker. The limits span the full min-max range of the marker across
    # ALL pooled cells. Per-file subsets only ever cover part of that range.
    marker_limits <- color_limits(combined_matrix[, marker])

    # CONJOINT (all cells) — colour vector is combined_matrix[, marker], which is
    # index-aligned to umap_df rows.
    df_all <- umap_df
    df_all$expr <- combined_matrix[, marker]

    p_all <- ggplot(df_all, aes(x = UMAP1, y = UMAP2, color = expr)) +
      geom_point(size = pt_size, alpha = pt_alpha) +
      viridis::scale_color_viridis(name = paste0(marker, norm_tag),
                                   limits = marker_limits,
                                   oob = scales::squish) +
      umap_coord() +
      base_umap_theme +
      labs(title = paste0("UMAP - ", population, " - ", marker, " (all files)"))
    save_one(
      p_all,
      create_safe_filename(paste0("umap_", pop_tag, "_marker_", marker, "_ALL"))
    )

    # PER FILE — subset rows by file label; SAME coordinates and the SAME colour
    # scale as the conjoint plot, so colours are directly comparable.
    for (f in files_unique) {
      rows <- which(file_labels == f)               # index-aligned subset
      if (length(rows) == 0) next                    # nothing to plot for this file
      df_f <- umap_df[rows, , drop = FALSE]
      df_f$expr <- combined_matrix[rows, marker]    # aligned colour vector

      p_f <- ggplot(df_f, aes(x = UMAP1, y = UMAP2, color = expr)) +
        geom_point(size = pt_size, alpha = pt_alpha) +
        viridis::scale_color_viridis(name = paste0(marker, norm_tag),
                                     limits = marker_limits,
                                     oob = scales::squish) +
        umap_coord() +
        base_umap_theme +
        labs(title = paste0("UMAP - ", population, " - ", marker, " - ", f))
      save_one(
        p_f,
        create_safe_filename(
          paste0("umap_", pop_tag, "_marker_", marker, "_file_", f))
      )
    }
    cat("Saved marker plots for:", marker, "\n")
  }

  cat("\nDone. Plots written to:", normalizePath(output_dir), "\n")
  # Repeated here on purpose: the warning above scrolled past a long load, and
  # the plots must not be read as covering samples they do not contain.
  if (nzchar(dropped_report)) {
    cat("\nWARNING: ", dropped_report, "\n", sep = "")
  }
  invisible(list(umap_coords = umap_coords, umap_df = umap_df,
                 workspaces = names(gs_list),
                 dropped = resolved$dropped,
                 output_dir = normalizePath(output_dir)))
}
