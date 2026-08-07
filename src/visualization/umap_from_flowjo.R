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

  # LEAF: the rule of thumb — a cell is named by the deepest thing the gating
  # tree actually says about it. A cell inside CD4+/IFN+ is inside CD4+ too, and
  # the useful label is the former; the tree already encodes which is which, so
  # nothing has to be inferred from cell counts.
  #
  # This replaced size-based "specificity", which used the SMALLEST gate as a
  # proxy for the most specific one. That proxy breaks on exactly the gates
  # people care about: the all-negative boolean leaf (CD107a-IFNg-IL2-TNF-) is
  # usually the biggest population in the sample, so smallest-wins labelled
  # those cells by their parent CD4+ and the leaf never appeared. Size is not
  # specificity; position in the tree is.
  #
  # A leaf is a candidate with no other candidate beneath it. A cell is labelled
  # by the leaf it is in; a cell that reaches only interior gates (its branch
  # was not gated further) keeps the deepest of those, because "CD4+" is the
  # true answer for a CD4+ cell nobody sub-gated.
  if (assignment == "leaf") {
    is_leaf <- !vapply(subpop_labels, function(lab)
      any(startsWith(subpop_labels, paste0(lab, "/"))), logical(1))

    # A cell is in SEVERAL leaves at once, and that is normal rather than an
    # error: the fixture gates CD4+ into CD107a-, IFNg-, IL2-, TNF- AND the
    # four-way combination CD107a-IFNg-IL2-TNF-, all leaves at the same depth.
    # Calling that "(multiple)" throws away the answer — the four-way gate is
    # plainly the most informative of them.
    #
    # So: leaves beat interior gates (the rule), then depth, then the SMALLEST
    # gate. Size is a poor definition of specificity on its own — that is what
    # this replaced — but it is a good tie-break BETWEEN LEAVES, where the
    # smaller gate is the more restrictive statement about the same cell.
    sizes <- colSums(mem)
    pick <- function(cols) {
      if (!length(cols)) return(NA_integer_)
      cols <- cols[subpop_depths[cols] == max(subpop_depths[cols])]
      cols[which.min(sizes[cols])]
    }
    out <- rep(none_label, n)
    for (i in seq_len(n)) {
      inside <- which(mem[i, ])
      if (!length(inside)) next
      leaves <- inside[is_leaf[inside]]
      chosen <- if (length(leaves)) pick(leaves) else pick(inside)
      out[i] <- subpop_labels[chosen]
    }
    return(out)
  }

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
# -----------------------------------------------------------------------------
# One assignment PER GATING LEVEL, instead of one assignment overall.
#
# The default "furthest" mode labels each cell by the deepest gate it reaches,
# which means a CD4+ cell that also falls in CD4+/CD107a+ is reported as
# CD107a+ and disappears from the CD4+/CD8+ picture entirely. That is the right
# answer to "what is this cell, most specifically" and the wrong answer to
# "how do the lineages split" — and both questions get asked of the same
# embedding.
#
# So this returns one labelling per depth relative to the reference population:
# level 1 is the immediate children (CD4+ / CD8+), level 2 their children
# (CD107a+ / IL-2+ / …), and so on. A cell is labelled at a level only by the
# gates AT that level; a cell of the reference population in none of them is
# NA, not "(none)".
#
# NA rather than a label is the point: these cells are the parent population
# doing nothing at this level, they are not a category, and drawing them as one
# more coloured class in the legend makes the classes that matter unreadable.
# The plotting side renders NA grey and keeps it out of the legend.
#
# Returns a named list, "level 1", "level 2", …, each a factor with one entry
# per cell. Factors rather than character vectors because at seven million
# cells the difference is an integer vector versus millions of string pointers.
# -----------------------------------------------------------------------------
assign_subpop_levels <- function(mem, subpop_depths, subpop_labels,
                                 multiple_label = .SUBPOP_MULTI,
                                 spec = NULL) {
  n <- nrow(mem)
  out <- list()
  if (n == 0L || ncol(mem) == 0L) return(out)

  # AN EXPLICIT SPEC ALWAYS WINS. The automatic rule below is a guess about
  # which gates constitute a level, and on a tree that mixes lineage with
  # functional gates at the same depth there is no guess that works: the
  # reference fixture has CD4+, CD8+, bare CD107a+/IFNg-/TNF- and sixteen
  # boolean combinations all at depth 1, several of which also branch. A cell
  # is in a handful of them at once and every automatic grouping collapses to
  # "(multiple)". When that happens, say which gates you mean:
  #     --subpop_levels_select "CD4+,CD8+;CD4+/CD107a+,CD8+/CD107a+"
  # Levels separated by ';', gates within a level by ','. Names match the way
  # --subpop_select matches: relative path, or bare gate name.
  if (!is.null(spec) && length(spec)) {
    for (i in seq_along(spec)) {
      want <- spec[[i]]
      at <- which(subpop_labels %in% want | basename(subpop_labels) %in% want)
      if (!length(at)) {
        warning("subpop level ", i, " names no gate that exists: ",
                paste(want, collapse = ", "))
        next
      }
      out[[paste("level", i)]] <- .subpop_level_factor(
        mem[, at, drop = FALSE], basename(subpop_labels[at]), n, multiple_label)
    }
    return(out)
  }

  # WHICH GATES BELONG TO A LEVEL. Depth alone is not enough. A real workspace
  # hangs the lineage split AND a pile of flat functional gates off the same
  # node: under "Viable cells" this fixture has CD4+ and CD8+ next to CD107a+,
  # IFNg-, and sixteen boolean combinations — 27 gates at depth 1. Those are not
  # alternatives, they are several different questions, so every cell is inside
  # a handful of them at once and grouping by depth alone labels the entire
  # plot "(multiple)".
  #
  # A level is therefore built from the gates at that depth which BRANCH — the
  # ones with a deeper gate beneath them. Those are the structural nodes, the
  # ones the tree is actually organised by, and at depth 1 that is exactly
  # {CD4+, CD8+}. Where a depth has no branching gates (the bottom of the tree,
  # where the answers live) its leaves are used instead.
  has_child <- vapply(subpop_labels, function(lab) {
    any(startsWith(subpop_labels, paste0(lab, "/")))
  }, logical(1))

  for (d in sort(unique(subpop_depths))) {
    at <- which(subpop_depths == d)
    if (!length(at)) next
    branching <- at[has_child[at]]
    dropped <- length(at) - length(branching)
    if (length(branching)) {
      if (dropped > 0L) {
        message("level ", d, ": using the ", length(branching),
                " gate(s) that branch; ", dropped,
                " leaf gate(s) at this depth are shown at their own level.")
      }
      at <- branching
    }
    sub_mem <- mem[, at, drop = FALSE]
    # Label by the LEAF name at this level, not the relative path: at level 2
    # "CD4+/CD107a+" and "CD8+/CD107a+" are the same functional gate under
    # different parents, and showing the full path makes a legend nobody reads.
    # The parent split is exactly what level 1 already shows.
    lab <- basename(subpop_labels[at])

    out[[paste("level", d)]] <- .subpop_level_factor(sub_mem, lab, n,
                                                    multiple_label)
  }

  # A level where most cells are "(multiple)" is not a partition, and the plot
  # will be one colour. Said once, loudly, with the way out — this is the most
  # likely reason a level plot looks wrong.
  for (ln in names(out)) {
    v <- out[[ln]]
    frac <- sum(!is.na(v) & v == multiple_label) / max(1L, sum(!is.na(v)))
    if (frac > 0.5) {
      warning(ln, ": ", round(100 * frac), "% of gated cells are in more than ",
              "one gate at this level, so the plot will be mostly one colour. ",
              "These gates are not alternatives. Name the ones you mean with ",
              "--subpop_levels_select, e.g. \"CD4+,CD8+\".", call. = FALSE)
    }
  }
  out
}

# One level's labelling: the gate a cell is in, "(multiple)" when several, NA
# when none. Shared by the explicit and the automatic paths so they cannot
# disagree about what a level plot means.
.subpop_level_factor <- function(sub_mem, lab, n, multiple_label) {
  cnt <- rowSums(sub_mem)
  vals <- rep(NA_character_, n)
  multi <- cnt > 1L
  vals[multi] <- multiple_label
  singles <- which(cnt == 1L)
  if (length(singles)) {
    col_idx <- max.col(sub_mem[singles, , drop = FALSE], ties.method = "first")
    vals[singles] <- lab[col_idx]
  }
  factor(vals, levels = c(unique(lab), if (any(multi)) multiple_label))
}

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
# `dictionary` (channel -> marker, from channel_marker_dictionary()) fills in
# marker names the FILE does not carry. More than half of a real cohort can ship
# with no $PnS at all, and without this every one of those samples has zero
# biological columns — it contributes nothing to the panel intersection and
# takes the run down with it. `aliases` reconciles the spellings that remain.
# Both are applied HERE, at the one point both marker paths already go through,
# so the panel reported by list_populations mode still cannot drift from the one
# a run uses.
#
# Returns list(keep = logical per column, markers = per-column marker names,
# dropped = human-readable "channel [marker]" strings for the excluded columns,
# recovered = channels whose marker name came from the dictionary).
# Shared by extract_population_sample() and gs_common_markers().
# -----------------------------------------------------------------------------
select_marker_columns <- function(channels, markers,
                                  exclude_pattern = "Time|FSC|SSC|scatter",
                                  dictionary = NULL, aliases = NULL) {
  recovered <- character(0)
  if (!is.null(dictionary) && length(dictionary)) {
    rec       <- recover_marker_names(channels, markers, dictionary)
    markers   <- rec$markers
    recovered <- rec$recovered
  }
  markers <- apply_marker_aliases(markers, aliases)

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

  list(keep = keep, markers = markers, dropped = dropped, recovered = recovered)
}

# The flowFrame/cytoframe form. Reads parameters() rather than going through
# explore_data(), which materialises exprs() once per column — the same answer
# for a fraction of the work, and it keeps this and gs_common_markers() on one
# code path so the panel reported by list_populations cannot drift from the one
# a run uses.
marker_columns <- function(ff, exclude_pattern = "Time|FSC|SSC|scatter",
                           dictionary = NULL, aliases = NULL) {
  cm <- cf_channel_markers(ff)
  select_marker_columns(cm$channels, cm$markers, exclude_pattern,
                        dictionary = dictionary, aliases = aliases)
}

# -----------------------------------------------------------------------------
# Every loaded sample's channel/marker metadata, in one pass over the GatingSets.
#
# Metadata only — no events are read — so this is cheap even across a few
# hundred samples, and it is what channel_marker_dictionary() needs to learn
# what each detector carries before any sample is extracted.
# -----------------------------------------------------------------------------
gs_list_panels <- function(gs_list) {
  out <- list()
  for (w in names(gs_list)) {
    gs <- gs_list[[w]]
    for (s in flowWorkspace::sampleNames(gs)) {
      cm <- cf_channel_markers(flowWorkspace::gh_pop_get_data(gs[[s]], y = "root"))
      out[[length(out) + 1L]] <- list(workspace = w, sample = s,
                                      channels = cm$channels,
                                      markers = cm$markers)
    }
  }
  out
}

# How many samples carry no marker name at all, and in which workspaces. This is
# the number that explains an empty panel intersection, so it is reported before
# anything tries to use the markers.
format_unnamed_panel_report <- function(panels, max_show = 10L) {
  bare <- Filter(function(p) !any(!is.na(p$markers) & nzchar(p$markers)), panels)
  if (!length(bare)) return("")
  by_wsp <- table(vapply(bare, function(p) p$workspace, character(1)))
  lines <- c(sprintf(
    "%d of %d sample(s) carry NO marker names ($PnS) in the FCS file itself.",
    length(bare), length(panels)),
    "Their markers are being recovered from the channel each was measured on, using",
    "the samples that DO name theirs. Without that they would have no biological",
    "columns at all and could not join the shared embedding.",
    sprintf("    %-40s %d sample(s)", names(by_wsp), as.integer(by_wsp)))
  if (length(by_wsp) > max_show) {
    lines <- c(utils::head(lines, 4L + max_show),
               sprintf("    ... and %d more workspace(s)",
                       length(by_wsp) - max_show))
  }
  paste(lines, collapse = "\n")
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
                              exclude_pattern = "Time|FSC|SSC|scatter",
                              dictionary = NULL, aliases = NULL) {
  per_sample <- lapply(sample_names, function(s) {
    # Metadata only: no rows are read, and no column is subset. At the root node
    # a sample is the whole FCS file, so anything that materialises exprs() here
    # is both wasteful and — since flowWorkspace 4.14 — a crash. See
    # cf_channel_markers().
    cm  <- cf_channel_markers(flowWorkspace::gh_pop_get_data(gs[[s]], y = "root"))
    sel <- select_marker_columns(cm$channels, cm$markers, exclude_pattern,
                                 dictionary = dictionary, aliases = aliases)
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
                                      exclude_pattern = "Time|FSC|SSC|scatter",
                                      dictionary = NULL, aliases = NULL,
                                      levels = FALSE, level_spec = NULL) {
  # NOTE: data is already compensated + transformed in the FlowJo working space.
  # We take exprs() AS-IS, no further transform.
  ff   <- flowWorkspace::gh_pop_get_data(gh, y = population)
  expr <- flowCore::exprs(ff)

  sel     <- marker_columns(ff, exclude_pattern, dictionary = dictionary,
                            aliases = aliases)
  keep    <- sel$keep
  markers <- sel$markers

  # Log excluded channels (channel + marker) so a wrongly-dropped biological
  # channel is visible to the user.
  if (length(sel$dropped) > 0) {
    cat("    excluded channels:", paste(sel$dropped, collapse = ", "), "\n")
  }

  # NULL, not stop(): one sample whose channels carry no marker name — and whose
  # detectors the channel dictionary could not name either — is a sample to
  # leave out, not a reason to abandon the other three hundred. The caller
  # collects these and reports them together; if it turns out to be ALL of them,
  # that is where it becomes an error.
  if (!any(keep)) return(NULL)

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
  # Built ONCE and used twice: the overall assignment and the per-level ones
  # read the same membership matrix, so a cell can never be in a gate for one
  # picture and out of it for the other.
  mem <- get_subpop_membership(gh, viable_idx, subpop_paths)
  subpop <- assign_subpop(mem, subpop_depths, subpop_labels,
                          assignment = assignment)
  stopifnot(length(subpop) == nrow(mtx))

  # Character, not factor, per sample: the levels differ between samples (a gate
  # can be empty here and populated there) and unifying them once at
  # concatenation is cheaper and safer than reconciling factors sample by sample.
  levels_out <- NULL
  if (isTRUE(levels)) {
    levels_out <- lapply(assign_subpop_levels(mem, subpop_depths, subpop_labels,
                                              spec = level_spec),
                         as.character)
    for (v in levels_out) stopifnot(length(v) == nrow(mtx))
  }

  list(mtx = mtx, subpop = subpop, levels = levels_out)
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

# The run-wide display label for one sample of one GATINGSET.
#
# The whole point is the key translation. `w` names a GatingSet, and a workspace
# split by channel set has SEVERAL of those ("MTB016_1_2 [panel 2]"), while
# sample_index is keyed per WORKSPACE — everything a person reads is. `gs_source`
# is the bridge between the two, and crossing it in one named place is what stops
# a caller from matching a part name against the workspace column, getting
# character(0) back, and failing several lines later with an error about an index
# rather than about the split.
#
# A miss here is a bug in this file — sample_index is built from these same
# GatingSets — not anything the data can cause, so it is named as one.
sample_display_label <- function(sample_index, gs_source, w, s) {
  source_wsp <- gs_source[[w]]
  hit <- sample_index$display[sample_index$workspace == source_wsp &
                                sample_index$raw == s]
  if (length(hit) != 1L) {
    stop(sprintf(paste0("internal: sample '%s' of GatingSet '%s' matched %d ",
                        "row(s) of the sample index (workspace '%s')"),
                 s, w, length(hit), source_wsp), call. = FALSE)
  }
  hit
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
# =============================================================================
# extract_flowjo_population() — everything up to the embedding
# =============================================================================
# Discovers the workspaces, pairs each with its FCS files, builds a GatingSet
# per workspace, reconciles the marker panels, resolves the population in every
# workspace, extracts its cells sample by sample and concatenates them.
#
# It stops exactly where the analysis stops being FlowJo-specific. What comes
# back is the repo's ordinary expression-matrix contract —
# `list(combined_matrix, sample_labels)`, the same shape
# build_combined_expr_matrix() produces — plus the per-cell labels only a
# FlowJo gating tree can supply. Downstream that means the shared pipeline
# steps (homogenize, downsample, dimred) work on it unchanged, and the
# embedding is run_dimred()'s job rather than a private call to run_umap():
# t-SNE and diffusion maps come for free, and there is one fewer copy of the
# UMAP parameter list in the repo.
#
# `list_populations = TRUE` returns after the report instead of extracting
# anything, which is the inspect entry point the GUI and --list_populations use.
#
# @return list(combined_matrix, sample_labels, file_labels, workspace_labels,
#   subpop_labels_per_cell, subpop_levels, n_cells_per, sample_index,
#   sample_annotations, annot_cols, population, pop_by_wsp, use_markers,
#   workspaces, dropped, report). With list_populations, list(populations=...)
#   and `inspected = TRUE`.
# =============================================================================
extract_flowjo_population <- function(wsp_path, fcs_dir,
                             population = "Lymphocytes/Single Cells/Viable cells",
                             group = "1", expected_samples = 0L,
                             sample_metadata = NULL, sample_pattern = NULL,
                             annotate = NULL, sample_id_column = NULL,
                             markers = NULL, subpop_mode = "furthest",
                             subpop_select = NULL, normalize = NULL,
                             normalize_data = "none", normalize_umap = "none",
                             subpop_levels = TRUE,
                             subpop_levels_select = NULL,
                             downsample = 0L, seed = 42L,
                             list_populations = FALSE,
                             allow_missing_samples = FALSE,
                             population_map = NULL,
                             match_mode = "content",
                             marker_aliases = NULL,
                             recover_markers = TRUE,
                             workspaces = NULL,
                             report_json = NULL) {
  require(CytoML)
  require(flowWorkspace)
  require(flowCore)

  # NORMALISATION HAPPENS IN TWO INDEPENDENT PLACES, and conflating them was the
  # complaint: changing the scale to recolour a figure also moved every point,
  # because the same matrix fed the embedding and the colour scales.
  #
  #   normalize_data  the values the EMBEDDING is computed from. Changing it is
  #                   a different analysis and the points move — as they should.
  #   normalize_umap  the values the MARKER COLOURS are taken from, applied
  #                   after the embedding exists. Changing it recolours the same
  #                   coordinates and nothing moves.
  #
  # `normalize` is kept as a deprecated alias for normalize_data so existing
  # commands, params files and the app's saved form state keep working.
  VALID_NORMALIZE <- c("none", "zscore", "minmax")
  if (!is.null(normalize) && nzchar(normalize)) {
    if (!identical(normalize_data, "none")) {
      stop("Pass either --normalize (deprecated) or --normalize_data, not both.")
    }
    message("--normalize is deprecated; it now means --normalize_data. ",
            "Use --normalize_umap to change only the marker colour scales.")
    normalize_data <- normalize
  }
  for (nm in c("normalize_data", "normalize_umap")) {
    v <- get(nm)
    if (!v %in% VALID_NORMALIZE) {
      stop(nm, " must be one of: ", paste(VALID_NORMALIZE, collapse = ", "),
           " (got '", v, "').")
    }
  }
  VALID_SUBPOP_MODE <- c("furthest", "cd4cd8", "direct", "deepest", "none")
  if (!subpop_mode %in% VALID_SUBPOP_MODE) {
    stop("subpop_mode must be one of: ",
         paste(VALID_SUBPOP_MODE, collapse = ", "),
         " (got '", subpop_mode, "').")
  }
  VALID_MATCH_MODE <- c("content", "strict")
  if (!match_mode %in% VALID_MATCH_MODE) {
    stop("match_mode must be one of: ", paste(VALID_MATCH_MODE, collapse = ", "),
         " (got '", match_mode, "').")
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
  n_found   <- length(wsp_paths)

  # Names, not full paths: they share a root, and 33 repetitions of it push the
  # part that matters — what failed — off the top of the log.
  cat("Workspaces found (", n_found, ") under ", wsp_path, ":\n", sep = "")
  cat(strwrap(paste(wsp_names, collapse = ", "), width = 78,
              prefix = "  ", initial = "  "), sep = "\n")

  # Deselecting happens here, before the expensive part, so leaving a workspace
  # out makes the run shorter as well as smaller.
  keep_wsp <- select_workspaces(wsp_names, workspaces)
  if (!all(keep_wsp)) {
    cat("Excluded by --workspaces (", sum(!keep_wsp), "): ",
        paste(wsp_names[!keep_wsp], collapse = ", "), "\n", sep = "")
    wsp_paths <- wsp_paths[keep_wsp]
    wsp_names <- wsp_names[keep_wsp]
    if (!length(wsp_paths)) {
      stop("--workspaces excluded every workspace; nothing to load.",
           call. = FALSE)
    }
  }
  n_wsp <- length(wsp_paths)
  fcs_index <- index_fcs_dir(fcs_dir)
  cat("FCS files indexed: ", length(fcs_index), " distinct name(s) under ",
      fcs_dir, "\n", sep = "")

  cat("Resolving ", n_wsp, " workspace(s) (group = \"", group, "\", matching = \"",
      match_mode, "\") ...\n", sep = "")
  resolved <- resolve_workspace_specs(wsp_paths, wsp_names, fcs_index, group,
                                      allow_missing = allow_missing_samples,
                                      match_mode = match_mode,
                                      fcs_root = fcs_dir)
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
  # What "content" mode has NOT decided yet. Said before the build so a long
  # load is not silent about the samples it is still hoping to pair, and phrased
  # as pending rather than lost — the load records below are what settle it.
  deferred_report <- format_deferred_report(resolved$deferred)
  if (nzchar(deferred_report)) {
    cat("\n", deferred_report, "\n", sep = "")
  }
  if (!length(specs)) {
    stop("every workspace was dropped: not one FCS file named by any of the ",
         n_wsp, " workspace(s) was found under ", fcs_dir,
         ". This is a wrong --fcs_dir, not a gap in the data.", call. = FALSE)
  }

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
  #
  # A workspace normally yields ONE GatingSet named after itself. It yields
  # several when its samples were acquired with different channel sets, which
  # CytoML cannot hold in one cytoset — see build_workspace_gatingsets(). The
  # parts are separate entries in gs_list because everything downstream works
  # per GatingSet, and `gs_source` maps each back to the workspace it came from
  # so counts, specs and the population map stay per WORKSPACE.
  gs_list <- list()
  gs_source <- character(0)
  build_failed <- list()
  load_records <- list()
  for (label in names(specs)) {
    spec <- specs[[label]]
    cat(sprintf("\nBuilding GatingSet for '%s' — %s, %d listed sample(s)\n",
                label, describe_group(spec), length(spec$listed)))
    emit_progress("load", which(names(specs) == label) - 1L, length(specs), label)
    parts <- tryCatch(build_workspace_gatingsets(spec),
                      error = function(e) {
                        if (!allow_missing_samples) stop(e)
                        build_failed[[label]] <<- conditionMessage(e)
                        NULL
                      })
    if (!is.null(parts)) {
      for (part in names(parts)) {
        gs_list[[part]]   <- parts[[part]]
        gs_source[[part]] <- label
      }
      # A part that still refused after the split costs its own samples and
      # nothing else, but it is a real loss and is reported as one.
      failed_parts <- attr(parts, "failed_parts")
      for (part in names(failed_parts)) {
        build_failed[[part]] <- failed_parts[[part]]
      }
      # The only place that knows what ARRIVED rather than what was asked for.
      # Every count from here on is taken from these records — over ALL parts,
      # so a split workspace is still one workspace with one sample tally.
      load_records[[label]] <- workspace_load_record(
        spec, unlist(lapply(parts, flowWorkspace::sampleNames),
                     use.names = FALSE))
    }
  }
  if (length(build_failed)) {
    cat("\nWARNING: ", length(build_failed),
        " workspace(s) or channel-set part(s) resolved their files but could\n",
        "NOT be loaded, and are not in the analysis. Each reason is CytoML's own:\n",
        sep = "")
    for (label in names(build_failed)) {
      cat(sprintf("  %s\n    %s\n", label,
                  gsub("\n", "\n    ", build_failed[[label]])))
      if (is_channel_inconsistency(build_failed[[label]])) {
        cat("    (its samples do not all carry the same channels and could not be\n",
            "     split into consistent groups either — the FCS files whose names\n",
            "     the workspace records were not all found, so their channels are\n",
            "     unknown. Check --fcs_dir covers every acquisition.)\n", sep = "")
      }
    }
  }
  if (!length(gs_list)) {
    stop("no workspace could be loaded; see the failures above.", call. = FALSE)
  }

  # What DID work, stated as plainly as what did not, and counted from the
  # GatingSets themselves. Counting the RESOLVED files here instead — as this
  # used to — understates a run by every sample CytoML paired by content: on a
  # renamed cohort that was the difference between a reported 26 and an actual
  # 311, with the headline claiming the smaller number.
  n_samples <- sum(vapply(load_records, function(r) r$n_loaded, integer(1)))
  cat(sprintf("\nGOING INTO THE UMAP: %d of %d workspace(s), %d sample(s).\n",
              length(load_records), n_wsp, n_samples))
  if (length(gs_list) > length(load_records)) {
    cat(sprintf("  (%d of them were split by channel set into %d parts; see above.)\n",
                sum(table(gs_source) > 1L), length(gs_list)))
  }
  routes_line <- format_match_routes(load_records)
  if (nzchar(routes_line)) cat(routes_line, "\n", sep = "")
  records_report <- format_load_records(load_records)
  if (nzchar(records_report)) cat("\n", records_report, "\n", sep = "")

  # --- 1c. Reconcile the marker panels. --------------------------------------
  # Done once, over every loaded sample's channel metadata, BEFORE anything is
  # extracted — so the panel the run uses and the panel list_populations reports
  # are built from exactly the same dictionary.
  marker_alias_map <- parse_marker_aliases(marker_aliases)
  panels <- gs_list_panels(gs_list)
  unnamed_report <- format_unnamed_panel_report(panels)
  if (nzchar(unnamed_report)) cat("\n", unnamed_report, "\n", sep = "")

  channel_dict <- if (isTRUE(recover_markers)) {
    channel_marker_dictionary(panels, aliases = marker_alias_map)
  } else {
    list(dictionary = character(0),
         conflicts = data.frame(channel = character(0), names = character(0),
                                chosen = character(0), counts = character(0),
                                stringsAsFactors = FALSE))
  }
  dict_report <- format_channel_dictionary(channel_dict)
  if (nzchar(dict_report)) cat("\n", dict_report, "\n", sep = "")

  # What each workspace's panel looks like AFTER recovery and aliasing, and what
  # is still not shared. An empty intersection is the most common reason a
  # multi-workspace run cannot proceed, so the divergence is named rather than
  # left to be inferred from "Common markers (0)".
  # Keyed by the SOURCE workspace: a workspace split by channel set is one panel
  # to the person reading this, and its parts differing is exactly the kind of
  # thing the marker dictionary above exists to reconcile.
  panel_by_wsp <- lapply(split(panels, vapply(panels, function(p)
                                 gs_source[[p$workspace]], character(1))),
    function(ps) {
      unique(unlist(lapply(ps, function(p) {
        sel <- select_marker_columns(p$channels, p$markers,
                                     dictionary = channel_dict$dictionary,
                                     aliases = marker_alias_map)
        sel$markers[sel$keep]
      }), use.names = FALSE))
    })
  panel_cmp <- compare_marker_sets(panel_by_wsp)
  if (length(panel_cmp$not_common)) {
    cat(sprintf("\n%d marker(s) are NOT shared by every workspace and cannot enter a joint UMAP:\n",
                length(panel_cmp$not_common)))
    for (m in panel_cmp$not_common) {
      have <- colnames(panel_cmp$presence)[panel_cmp$presence[m, ]]
      cat(sprintf("    %-16s in %d of %d: %s\n", m, length(have),
                  ncol(panel_cmp$presence), paste(have, collapse = ", ")))
    }
    sugg <- suggest_marker_aliases(panel_by_wsp)
    if (nrow(sugg)) {
      cat("\n  Some of those look like the same marker spelled differently. Nothing is\n",
          "  merged automatically — pass the ones you agree with as marker_aliases:\n",
          sep = "")
      for (i in seq_len(nrow(sugg))) {
        cat(sprintf("    %-24s -> %-12s (%s)\n", sugg$names[i],
                    sugg$suggested[i], sugg$confidence[i]))
      }
      others <- lapply(strsplit(sugg$names, ", ", fixed = TRUE), identity)
      pairs <- unlist(lapply(seq_len(nrow(sugg)), function(i)
        sprintf("%s=%s", setdiff(others[[i]], sugg$suggested[i]),
                sugg$suggested[i])))
      cat("    e.g. --marker_aliases \"", paste(pairs, collapse = ";"), "\"\n",
          sep = "")
    }
  }

  # A population map is also a workspace SELECTION: a workspace with no entry
  # contributes no gate and therefore no cells. Applied here, before anything
  # downstream is derived, so sample labels, annotations and the shared-file
  # accounting all describe the run that is actually happening.
  # Keyed on WORKSPACE labels, not on gs_list names: a workspace split by
  # channel set must still be selected by its own name, not by "[panel 2]".
  pop_map <- parse_population_map(population_map, names(load_records))
  if (!is.null(pop_map) && !isTRUE(list_populations)) {
    dropped_by_map <- setdiff(names(load_records), names(pop_map))
    gs_list   <- gs_list[gs_source[names(gs_list)] %in% names(pop_map)]
    gs_source <- gs_source[names(gs_list)]
    load_records <- load_records[names(load_records) %in% names(pop_map)]
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
    data.frame(workspace = rep(gs_source[[w]], length(raw)), raw = raw,
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

  # Everything above is keyed per GATINGSET, because that is the unit a gate has
  # to resolve in. Everything a person reads is keyed per WORKSPACE. When a
  # workspace was split by channel set the two differ, and the bridge is an
  # INTERSECTION, never a union: a gate only one part carries is not a gate that
  # workspace can offer the run, and claiming otherwise would put it in the
  # health panel as available and then fail at extraction.
  by_source <- function(x) {
    if (!length(x)) return(x)
    grouped <- split(names(x), unname(gs_source[names(x)]))
    out <- lapply(grouped, function(ks) Reduce(intersect, x[ks]))
    out[intersect(names(load_records), names(out))]
  }

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
      gs_common_markers(gs_list[[w]], flowWorkspace::sampleNames(gs_list[[w]]),
                        dictionary = channel_dict$dictionary,
                        aliases = marker_alias_map))
    names(markers_by_wsp) <- names(gs_list)
    markers_by_source <- by_source(markers_by_wsp)
    common_markers <- Reduce(intersect, markers_by_wsp)
    common_markers <- markers_by_wsp[[1]][markers_by_wsp[[1]] %in% common_markers]
    cat("\nCommon markers (", length(common_markers), "):\n", sep = "")
    for (m in common_markers) cat("  ", m, "\n", sep = "")
    if (length(common_markers) == 0) {
      cat("  (no marker common to all samples — see the per-workspace panels ",
          "below)\n", sep = "")
    }

    # Each workspace's own panel. The intersection above collapses to nothing
    # the moment one workspace spells a marker differently, and then the only
    # way to see WHY is to compare the panels side by side. Four-space indent,
    # same contract as the other per-workspace blocks.
    cat("\nMarkers by workspace (", length(markers_by_source), "):\n", sep = "")
    for (w in names(markers_by_source)) {
      cat(sprintf("    %s | %d marker(s)\n", w, length(markers_by_source[[w]])))
      for (m in markers_by_source[[w]]) cat("      ", m, "\n", sep = "")
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
    pops_by_source <- by_source(pops_by_wsp)
    cat("\nPopulations by workspace (", length(pops_by_source), "):\n", sep = "")
    for (w in names(pops_by_source)) {
      cat(sprintf("    %s | %d node(s)\n", w, length(pops_by_source[[w]])))
      for (p in pops_by_source[[w]]) cat("      ", p, "\n", sep = "")
    }

    # Counted from the GatingSets, not from the resolved file list: a sample
    # CytoML paired by content is in the analysis and has to be in the count.
    cat("\nWorkspaces resolved (", length(load_records), "):\n", sep = "")
    for (label in names(load_records)) {
      spec  <- specs[[label]]
      parts <- names(gs_source)[gs_source == label]
      cat(sprintf("    %s | %s | %s | %d sample(s)%s\n",
                  label, spec$wsp, describe_group(spec),
                  load_records[[label]]$n_loaded,
                  if (length(parts) > 1L)
                    sprintf(" | %d channel-set parts", length(parts)) else ""))
    }

    divergent <- setdiff(unique(unlist(pops_by_source, use.names = FALSE)),
                         shared_pops)
    cat("\nPopulations not shared by all workspaces (", length(divergent),
        "):\n", sep = "")
    for (p in divergent) {
      absent <- names(pops_by_source)[!vapply(pops_by_source,
                                              function(x) p %in% x, logical(1))]
      cat(sprintf("    %s   [missing in: %s]\n", p,
                  paste(absent, collapse = ", ")))
    }
    if (length(divergent) == 0 && n_wsp > 1) {
      cat("    (none — every workspace has the same gating tree)\n")
    }

    # The same load, as data. Written after the prose so that a failure to write
    # it cannot cost anyone the report they were reading. See
    # visualization/umap_load_report.R for why this exists.
    write_load_report(build_load_report(
      inputs = list(wsp = wsp_path, fcs_dir = fcs_dir, group = group,
                    match_mode = match_mode,
                    allow_missing_samples = isTRUE(allow_missing_samples),
                    recover_markers = isTRUE(recover_markers)),
      workspaces_found = wsp_names,
      load_records = load_records, specs = specs,
      problems = resolved$problems, dropped = resolved$dropped,
      deferred = resolved$deferred, build_failed = build_failed,
      pops_by_wsp = pops_by_source, shared_pops = shared_pops,
      markers_by_wsp = markers_by_source, common_markers = common_markers,
      panel_cmp = panel_cmp,
      suggestions = suggest_marker_aliases(panel_by_wsp),
      channel_dict = channel_dict, samples = sample_index,
      annotations = stats::setNames(
        lapply(annot_cols, function(cl) unique(sample_annotations[[cl]])),
        annot_cols),
      unnamed_samples = sum(vapply(panels, function(p)
        !any(!is.na(p$markers) & nzchar(p$markers)), logical(1))),
      panel_parts = as.list(table(gs_source))
    ), report_json)

    cat("\nSelect one with population = (leading '/' optional), e.g.:\n")
    if (length(shared_pops)) {
      cat("  population = \"", sub("^/", "", shared_pops[length(shared_pops)]),
          "\"\n", sep = "")
    }
    return(invisible(list(inspected = TRUE, populations = shared_pops,
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
      try_resolve_population(pop_map[[gs_source[[w]]]], pops_by_wsp[[w]]))
    names(resolutions) <- names(gs_list)
  } else {
    resolutions <- lapply(pops_by_wsp, function(p)
      try_resolve_population(population, p))
  }
  failed <- names(resolutions)[!vapply(resolutions, function(r) r$ok, logical(1))]
  if (length(failed)) {
    detail <- vapply(failed, function(w) sprintf(
      "  %s (%s)\n    %s", w, specs[[gs_source[[w]]]]$wsp,
      resolutions[[w]]$reason),
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

  # 'furthest' names each cell by the deepest gate the TREE puts it in — its
  # leaf. The other modes assign by path depth (ties -> "(multiple)").
  subpop_assignment <- if (subpop_mode == "furthest" && is.null(subpop_select)) {
    "leaf"
  } else {
    "depth"
  }

  # "CD4+,CD8+;CD4+/CD107a+,CD8+/CD107a+" -> list(c("CD4+","CD8+"), c(...)).
  # ';' between levels and ',' within one, the same separators --population_map
  # and --marker_aliases already use.
  level_spec <- NULL
  if (!is.null(subpop_levels_select) && nzchar(trimws(subpop_levels_select))) {
    level_spec <- lapply(strsplit(trimws(subpop_levels_select), ";",
                                  fixed = TRUE)[[1]], function(g) {
      v <- trimws(strsplit(g, ",", fixed = TRUE)[[1]]); v[nzchar(v)]
    })
    level_spec <- Filter(length, level_spec)
    cat("Subpopulation levels given explicitly: ", length(level_spec), "\n",
        sep = "")
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
  no_marker_samples <- character(0)
  per_sample_levels <- list()
  # Progress is reported against the sample count the GatingSets actually hold,
  # which is the only total that is known before the loop and true after it.
  n_total_samples <- sum(vapply(gs_list, function(g)
    length(flowWorkspace::sampleNames(g)), integer(1)))
  n_done <- 0L
  for (w in names(gs_list)) {
    gs <- gs_list[[w]]
    if (n_wsp > 1) cat(sprintf("\n[%s]\n", w))
    for (s in flowWorkspace::sampleNames(gs)) {
      # The display label is the key everywhere downstream (legends, file
      # names, per-file subsetting), so switch to it here and never look back.
      # `w` is a GatingSet name and sample_index is keyed per workspace — see
      # sample_display_label(), which owns that translation.
      display <- sample_display_label(sample_index, gs_source, w, s)
      res <- extract_population_sample(gs[[s]], pop_by_wsp[[w]],
                                       subpop_paths_by_wsp[[w]],
                                       subpop_depths, subpop_labels,
                                       assignment = subpop_assignment,
                                       dictionary = channel_dict$dictionary,
                                       aliases = marker_alias_map,
                                       levels = isTRUE(subpop_levels),
                                       level_spec = level_spec)
      # No usable channel on this sample. Collected rather than fatal — see
      # extract_population_sample() — and reported together below.
      if (is.null(res)) {
        no_marker_samples <- c(no_marker_samples, display)
        cat(sprintf("  %s: SKIPPED — no named marker channel\n", display))
        next
      }
      mtx <- res$mtx
      sp  <- res$subpop

      # Optional per-sample downsampling AFTER gating. The SAME row index
      # subsets the subpopulation labels — overall and per level — so they stay
      # aligned to mtx. Every aligned vector is subset HERE; carrying the index
      # out of this block would let a stale one from a previous sample apply.
      lv <- res$levels
      if (downsample > 0 && nrow(mtx) > downsample) {
        idx <- sort(sample(nrow(mtx), downsample))
        mtx <- mtx[idx, , drop = FALSE]
        sp  <- sp[idx]
        if (!is.null(lv)) lv <- lapply(lv, function(v) v[idx])
      }

      cat(sprintf("  %s: %d cells x %d markers\n", display, nrow(mtx),
                  ncol(mtx)))
      n_done <- n_done + 1L
      emit_progress("extract", n_done, n_total_samples, display)
      if (nrow(mtx) == 0) {
        warning("Sample '", display, "' has 0 cells in population '",
                pop_by_wsp[[w]],
                "' after gating/downsampling; it will not contribute rows.")
      }
      per_sample_mtx[[display]]    <- mtx
      per_sample_subpop[[display]] <- sp
      if (!is.null(lv)) per_sample_levels[[display]] <- lv
      # The SOURCE workspace, not the GatingSet: a workspace split by channel
      # set is still one workspace to everyone reading a legend.
      per_sample_wsp[display]      <- gs_source[[w]]
    }
  }

  # Samples with no usable channel at all. Said once, in full, rather than as
  # one line lost among three hundred — these cells are not in the embedding.
  if (length(no_marker_samples)) {
    cat(sprintf(
      paste0("\nWARNING: %d sample(s) are NOT in the embedding: no channel of theirs ",
             "carries a\nmarker name, and the channel dictionary could not name one ",
             "either. Their\ndetectors are named differently from the rest of the ",
             "cohort:\n"),
      length(no_marker_samples)))
    for (s in no_marker_samples) cat("    ", s, "\n", sep = "")
  }
  if (!length(per_sample_mtx)) {
    stop("no sample had a usable marker channel: every one of them is missing ",
         "$PnS and no other sample's channels could supply the names. Check the ",
         "panel, or pass marker_aliases.", call. = FALSE)
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

  # Per-level labels, concatenated in the SAME sample order as the matrix and
  # turned into factors only now — the level set is the union across samples,
  # so a gate that was empty in one sample still gets a stable colour. A sample
  # that produced no vector for a level (its gates were absent) contributes NA
  # for its rows rather than shifting everything after it.
  combined_levels <- list()
  if (length(per_sample_levels)) {
    level_names <- unique(unlist(lapply(per_sample_levels, names),
                                 use.names = FALSE))
    for (ln in level_names) {
      parts <- lapply(sample_labels, function(s) {
        v <- per_sample_levels[[s]][[ln]]
        if (is.null(v)) rep(NA_character_, n_cells_per[[s]]) else v
      })
      vec <- unlist(parts, use.names = FALSE)
      stopifnot(length(vec) == nrow(combined_matrix))
      lv <- sort(setdiff(unique(vec[!is.na(vec)]), .SUBPOP_MULTI))
      combined_levels[[ln]] <- factor(vec,
                                      levels = c(lv, .SUBPOP_MULTI[
                                        .SUBPOP_MULTI %in% vec]))
    }
    cat(sprintf("Subpopulation levels: %d (%s)\n", length(combined_levels),
                paste(names(combined_levels), collapse = ", ")))
  }

  cat(sprintf("\nConcatenated matrix: %d cells x %d markers (%d samples, %d workspace(s))\n",
              nrow(combined_matrix), ncol(combined_matrix),
              length(per_sample_mtx),
              length(unique(wsp_labels_per_cell))))

  # --- 5b. Per-marker normalisation, twice and independently. -----------------
  # colour_matrix is what the marker plots are coloured by; combined_matrix is
  # what the embedding is computed from. They start identical and are scaled
  # separately, so --normalize_umap can be changed to recolour a figure without
  # moving a single point. When the two methods agree (including the usual case
  # where both are "none") they stay the same object, so nothing is copied.
  colour_matrix <- combined_matrix
  if (normalize_data != "none") {
    cat(sprintf(
      "Normalising %d markers for the EMBEDDING (method=%s) from the concatenated gated cells...\n",
      ncol(combined_matrix), normalize_data))
    combined_matrix <- normalize_columns(combined_matrix, normalize_data)
  }
  if (identical(normalize_umap, normalize_data)) {
    colour_matrix <- combined_matrix
  } else if (normalize_umap != "none") {
    cat(sprintf(
      "Normalising %d markers for the COLOUR SCALES (method=%s); the embedding is unaffected.\n",
      ncol(colour_matrix), normalize_umap))
    colour_matrix <- normalize_columns(colour_matrix, normalize_umap)
  }


  # The pipeline contract (combined_matrix / sample_labels), plus the per-cell
  # labels that only a FlowJo gating tree can supply. Everything downstream —
  # homogenize, downsample, dimred — reads the first two and ignores the rest.
  list(
    combined_matrix = combined_matrix,
    # Same shape and row order as combined_matrix; the marker plots read this
    # one. Identical object unless normalize_umap differs from normalize_data.
    colour_matrix   = colour_matrix,
    normalize_data  = normalize_data,
    normalize_umap  = normalize_umap,
    sample_labels   = file_labels,          # per CELL, as the pipeline expects
    samples         = sample_labels,        # per SAMPLE, in matrix order
    file_labels     = file_labels,
    workspace_labels = wsp_labels_per_cell,
    subpop_labels_per_cell = combined_subpop,
    subpop_by_level        = combined_levels,
    subpop_levels   = subpop_labels,
    subpop_mode_tag = subpop_mode_tag,
    n_cells_per     = n_cells_per,
    sample_index    = sample_index,
    sample_annotations = sample_annotations,
    annot_cols      = annot_cols,
    population      = population,
    pop_by_wsp      = pop_by_wsp,
    use_markers     = use_markers,
    workspaces      = names(load_records),
    dropped         = resolved$dropped,
    dropped_report  = dropped_report,
    inspected       = FALSE
  )
}


# =============================================================================
# umap_from_flowjo() — extraction, then embedding, then the FlowJo-specific
# plots. Kept as the one-call entry point the standalone script and the tests
# use; the pipeline calls the two halves as separate steps.
# =============================================================================
umap_from_flowjo <- function(wsp_path, fcs_dir,
                             population = "Lymphocytes/Single Cells/Viable cells",
                             output_dir = "umap_flowjo",
                             group = "1", expected_samples = 0L,
                             sample_metadata = NULL, sample_pattern = NULL,
                             annotate = NULL, sample_id_column = NULL,
                             markers = NULL, subpop_mode = "furthest",
                             subpop_select = NULL, normalize = NULL,
                             normalize_data = "none", normalize_umap = "none",
                             subpop_levels = TRUE,
                             subpop_levels_select = NULL,
                             downsample = 0L, n_neighbors = 15L, min_dist = 0.1,
                             metric = "euclidean", seed = 42L,
                             point_size = 0.4, point_alpha = 0.6,
                             plot_width = 10, plot_height = 8,
                             output_formats = list(pdf = FALSE, png = FALSE, svg = TRUE),
                             list_populations = FALSE,
                             allow_missing_samples = FALSE,
                             population_map = NULL,
                             match_mode = "content",
                             marker_aliases = NULL,
                             recover_markers = TRUE,
                             workspaces = NULL,
                             report_json = NULL,
                             method = "umap",
                             n_cores = 1L,
                             embed_cores = 0L,
                             embed_sgd_cores = 1L,
                             umap_init = "auto",
                             umap_subsample = 0L,
                             nn_method = NULL) {
  require(ggplot2)
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for marker colour scales.")
  }

  # Resolved ONCE, here, so every phase and every log line agrees on the number
  # — and so "0" (decide for me) becomes a concrete count the run can report.
  n_cores <- resolve_cores(n_cores)
  # TWO embedding knobs, because uwot has two thread pools and only one of them
  # can change the answer (see run_umap in R_shared/stats_utils.R):
  #
  #   embed_cores      the NEIGHBOUR SEARCH — Annoy index, kNN query, smooth-kNN
  #                    calibration. Deterministic at any thread count, so it
  #                    defaults to the whole machine like everything else.
  #   embed_sgd_cores  the OPTIMISER. Threading it means uwot's batch optimiser,
  #                    which is a DIFFERENT optimiser: same seed, different
  #                    coordinates. Defaults to serial, and stays there unless
  #                    someone asks for it by name.
  #
  # These used to be one number defaulting to 1, on the grounds that the
  # embedding was "worth ~1% of a run — the plots are the cost". That holds at a
  # few hundred thousand cells and is badly wrong past a few million: a 9.8M-cell
  # run spent 22 minutes building the Annoy index on ONE thread of a 28-core
  # machine, before a single-threaded kNN query over the same 9.8M points. The
  # reproducibility that bought was never at risk from the search.
  embed_cores     <- resolve_cores(embed_cores)
  embed_sgd_cores <- resolve_cores(embed_sgd_cores)
  cat("Cores: ", n_cores, " for plotting, ", embed_cores,
      " for the neighbour search, ", embed_sgd_cores,
      " for the optimiser.\n", sep = "")

  ext <- extract_flowjo_population(
    wsp_path = wsp_path, fcs_dir = fcs_dir, population = population,
    group = group, expected_samples = expected_samples,
    sample_metadata = sample_metadata, sample_pattern = sample_pattern,
    annotate = annotate, sample_id_column = sample_id_column,
    markers = markers, subpop_mode = subpop_mode,
    subpop_select = subpop_select, normalize = normalize,
    normalize_data = normalize_data, normalize_umap = normalize_umap,
    subpop_levels = subpop_levels,
    subpop_levels_select = subpop_levels_select,
    downsample = downsample, seed = seed,
    list_populations = list_populations,
    allow_missing_samples = allow_missing_samples,
    population_map = population_map, match_mode = match_mode,
    marker_aliases = marker_aliases, recover_markers = recover_markers,
    workspaces = workspaces, report_json = report_json)

  # Inspect mode already printed and returned everything it has.
  if (isTRUE(ext$inspected)) return(invisible(ext))

  combined_matrix     <- ext$combined_matrix
  sample_labels       <- ext$samples
  file_labels         <- ext$file_labels
  wsp_labels_per_cell <- ext$workspace_labels
  combined_subpop     <- ext$subpop_labels_per_cell
  subpop_labels       <- ext$subpop_levels
  subpop_mode_tag     <- ext$subpop_mode_tag
  n_cells_per         <- ext$n_cells_per
  sample_index        <- ext$sample_index
  sample_annotations  <- ext$sample_annotations
  annot_cols          <- ext$annot_cols
  population          <- ext$population
  use_markers         <- ext$use_markers
  gs_list             <- stats::setNames(vector("list", length(ext$workspaces)),
                                         ext$workspaces)
  resolved            <- list(dropped = ext$dropped)
  dropped_report      <- ext$dropped_report

  # --- 6. Compute the embedding ONCE on the concatenated matrix. -------------
  # CRITICAL alignment guard: uwot may error or drop rows on non-finite input,
  # which would silently break the row-for-row alignment. Fail loudly here.
  if (!all(is.finite(combined_matrix))) {
    bad_cols <- colnames(combined_matrix)[
      apply(combined_matrix, 2, function(x) any(!is.finite(x)))]
    stop("combined_matrix contains non-finite (NA/NaN/Inf) values in marker(s): ",
         paste(bad_cols, collapse = ", "),
         ". Refusing to embed because dropped rows would break index ",
         "alignment with the per-file plots.")
  }

  # Through run_dimred(), not run_umap(): the dispatcher is what the rest of the
  # repo embeds with, it re-checks that the coordinate count still matches the
  # matrix, and going through it means `method` can be tsne or diffusion_map
  # here without another line of code. draw_plots = FALSE because the figures
  # below are this path's own — the standard by-sample plot is not one of them.
  cat("Computing shared ", method, " over ", nrow(combined_matrix),
      " cells x ", ncol(combined_matrix), " markers (n_neighbors=", n_neighbors,
      ", min_dist=", min_dist, ", metric=", metric, ", seed=", seed,
      ", search threads=", embed_cores, ")...\n", sep = "")
  if (embed_sgd_cores > 1L) {
    cat("  The OPTIMISER is threaded (", embed_sgd_cores, "), which means ",
        "uwot's batch optimiser: reproducible from the seed, but NOT identical ",
        "to a serial-optimiser run. The neighbour search's threads do not ",
        "affect the coordinates.\n", sep = "")
  }
  # Said BEFORE the call, so the bar names the phase the moment it is entered
  # rather than when uwot first speaks. From here on the sub-phases report
  # themselves: emit_uwot_stage() (R_shared/stats_utils.R) translates uwot's own
  # stage announcements into embed_nn / embed_init, and the epoch callback
  # carries the optimisation. This used to be a single "0 / 100" that never
  # moved until the first epoch, which on millions of cells meant an empty bar
  # for most of the run.
  emit_progress("embed_nn", 0L, 3L, "starting the neighbour search")
  dimred <- run_dimred(
    expr_matrix = combined_matrix,
    method      = method,
    output_dir  = output_dir,
    params = list(n_neighbors = n_neighbors, umap_min_dist = min_dist,
                  metric = metric, seed = seed, draw_plots = FALSE,
                  n_threads = embed_cores, sgd_threads = embed_sgd_cores,
                  nn_method = nn_method,
                  progress = TRUE, verbose = TRUE,
                  umap_init = umap_init, umap_subsample = umap_subsample,
                  output_formats = output_formats,
                  plot_width = plot_width, plot_height = plot_height)
  )
  umap_coords <- dimred$coords
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
  plot_flowjo_umap_facets(
    umap_df = umap_df, combined_matrix = combined_matrix,
    use_markers = use_markers, subpop_labels = subpop_labels,
    annot_cols = annot_cols, population = population,
    output_dir = output_dir, n_workspaces = length(gs_list),
    subpop_mode_tag = subpop_mode_tag, n_cores = n_cores,
    colour_matrix = ext$colour_matrix,
    subpop_by_level = ext$subpop_by_level %||% list(),
    normalize = ext$normalize_data, normalize_umap = ext$normalize_umap,
    point_size = point_size,
    point_alpha = point_alpha, plot_width = plot_width,
    plot_height = plot_height, output_formats = output_formats)

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


# =============================================================================
# plot_flowjo_umap_facets() — the figures only a FlowJo gating tree can produce
# =============================================================================
# One embedding, coloured many ways: by source file, by workspace, by the FlowJo
# subgate each cell falls in, by every marker, and by every sample annotation.
#
# Split out from umap_from_flowjo() so the pipeline can run it as its own step
# against coordinates a previous step computed — the embedding is expensive and
# should be paid for once, while re-colouring is cheap and is what people
# actually iterate on.
#
# `umap_df` must be row-aligned to `combined_matrix`; every caller builds it the
# same way and the guards below assume it.
# =============================================================================
plot_flowjo_umap_facets <- function(umap_df, combined_matrix, use_markers,
                                    subpop_labels, annot_cols, population,
                                    output_dir, n_workspaces = 1L,
                                    subpop_mode_tag = "",
                                    normalize = "none", n_cores = 1L,
                                    colour_matrix = NULL,
                                    normalize_umap = NULL,
                                    subpop_by_level = list(),
                                    point_size = 0.4, point_alpha = 0.6,
                                    plot_width = 10, plot_height = 8,
                                    output_formats = list(pdf = FALSE,
                                                          png = FALSE,
                                                          svg = TRUE)) {
  require(ggplot2)
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for marker colour scales.")
  }
  stopifnot(nrow(umap_df) == nrow(combined_matrix))
  # Marker colours come from colour_matrix; the embedding came from
  # combined_matrix. They are the same object unless --normalize_umap differs.
  if (is.null(colour_matrix)) colour_matrix <- combined_matrix
  stopifnot(nrow(colour_matrix) == nrow(umap_df),
            identical(colnames(colour_matrix), colnames(combined_matrix)))
  if (is.null(normalize_umap)) normalize_umap <- normalize
  n_wsp <- n_workspaces

  # umap_df already carries these, one row per cell, and they are the ONLY
  # copy now that extraction happens in a different process: bound back to the
  # names the plotting below uses so there is one source of truth for each and
  # no chance of a second, unaligned vector being passed in beside the frame.
  file_labels         <- umap_df$file
  wsp_labels_per_cell <- umap_df$workspace
  combined_subpop     <- umap_df$subpop

  # Whether there is a subpopulation partition to colour by is a property of the
  # gate set itself, so it is derived here rather than passed — one fewer
  # argument that a caller could set inconsistently with the labels it supplies.
  do_subpop <- length(subpop_labels) > 0

  # --- 7. Plotting. ----------------------------------------------------------
  make_outdir(output_dir)
  pt_size  <- point_size
  pt_alpha <- point_alpha

  pop_tag <- create_safe_filename(sub("^/", "", population))

  # Suffix appended to marker colour-scale titles so a reader can tell whether
  # the colour encodes raw or normalized expression.
  norm_tag <- switch(normalize_umap,
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
  # ggplot2 refuses a canvas over 50 inches. A per-sample legend reaches that on
  # a real cohort — 311 samples with names 48 characters long want a 14-column
  # legend about 70 inches wide — and the run died there rather than producing
  # the figure. Past a few dozen keys the legend has stopped being readable
  # anyway, so it is dropped and the plot says so: a figure whose colours are
  # unlabelled is worth more than no figure, as long as it does not pretend
  # otherwise. The other colourings (subpopulation, workspace, annotation) have
  # few enough levels that this never triggers.
  MAX_CANVAS_IN <- 45
  save_categorical <- function(p, basename_no_ext, legend_strings, ncol = 1L) {
    width <- plot_width + ncol * legend_width_in(legend_strings)
    if (width > MAX_CANVAS_IN) {
      n_keys <- max(0L, length(legend_strings) - 1L)
      cat(sprintf(
        "    %s: %d legend entries would need a %.0f-inch canvas; drawing it without a legend.\n",
        basename_no_ext, n_keys, width))
      p <- p + theme(legend.position = "none") +
        labs(caption = sprintf(
          "%d levels - legend omitted, it would not fit on the canvas", n_keys))
      width <- plot_width
    }
    save_plot_formats(
      p,
      file.path(output_dir, basename_no_ext),
      width   = width,
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

    # 7b-ii. ONE PLOT PER GATING LEVEL.
    #
    # The plot above answers "what is this cell, most specifically", which is
    # why a CD4+ cell that also fires IL-2 appears as IL-2 and vanishes from the
    # lineage picture. These answer "how does the population divide HERE": level
    # 1 is the immediate children, level 2 theirs, and so on.
    #
    # Cells of the reference population that are in no gate at a level are NA
    # and drawn grey, with na.translate = FALSE keeping them out of the legend.
    # They are the parent doing nothing at this level, not a category — giving
    # them a swatch is what made the earlier "(none)" class dominate the key.
    for (ln in names(subpop_by_level)) {
      lv_vec <- subpop_by_level[[ln]]
      if (all(is.na(lv_vec))) next
      keys   <- levels(droplevels(lv_vec))
      lv_col <- stats::setNames(distinct_palette(length(keys)), keys)
      lv_ncol <- legend_ncol(length(keys), plot_height)
      df_lv <- umap_df
      df_lv$level <- lv_vec
      # NA first so the grey parent cells sit UNDER the coloured gates rather
      # than speckling over them.
      df_lv <- df_lv[order(!is.na(df_lv$level)), , drop = FALSE]

      p_lv <- ggplot(df_lv, aes(x = UMAP1, y = UMAP2, color = level)) +
        geom_point(size = pt_size, alpha = pt_alpha) +
        scale_color_manual(values = lv_col, breaks = keys,
                           na.value = "grey85", na.translate = FALSE,
                           drop = FALSE) +
        umap_coord() +
        base_umap_theme +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                    ncol = lv_ncol)) +
        labs(title = paste0("UMAP - ", population, " - ", ln,
                            " (", sum(!is.na(lv_vec)), " of ", length(lv_vec),
                            " cells gated)"),
             color = ln)
      save_categorical(
        p_lv,
        create_safe_filename(paste0("umap_", pop_tag, "_by_",
                                    gsub(" ", "", ln))),
        legend_strings = c(keys, ln), ncol = lv_ncol)
    }
    cat("Saved: ", length(subpop_by_level), " per-level subpopulation plot(s)\n",
        sep = "")
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
  n_marker_done <- 0L
  for (marker in use_markers) {
    n_marker_done <- n_marker_done + 1L
    emit_progress("plot", n_marker_done, length(use_markers), marker)
    # ONE colour scale per marker, shared by the conjoint plot and every per-file
    # plot, so a given colour means the same expression value in every plot of
    # this marker. The limits span the full min-max range of the marker across
    # ALL pooled cells. Per-file subsets only ever cover part of that range.
    marker_limits <- color_limits(colour_matrix[, marker])

    # CONJOINT (all cells) — colour vector is combined_matrix[, marker], which is
    # index-aligned to umap_df rows.
    df_all <- umap_df
    df_all$expr <- colour_matrix[, marker]

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
    #
    # This is the run: 7 markers x 311 files is 2,177 of the 2,187 plots, and on
    # the study this was written for it was 2h54m of a 3h15m run. Every
    # iteration is independent — it reads two shared read-only objects, builds a
    # ggplot and writes one file whose name is unique to (marker, file) — so it
    # is spread across forked workers. Fork rather than a socket cluster because
    # umap_df and combined_matrix are then shared copy-on-write instead of being
    # serialised to each worker.
    #
    # The OUTER loop over markers stays serial: marker_limits has to be computed
    # once per marker anyway, progress reads sensibly per marker, and 311 items
    # is already far more parallelism than any core count here will absorb.
    plot_one_file <- function(f) {
      rows <- which(file_labels == f)               # index-aligned subset
      if (length(rows) == 0) return(invisible(NULL)) # nothing for this file
      df_f <- umap_df[rows, , drop = FALSE]
      df_f$expr <- colour_matrix[rows, marker]      # aligned colour vector

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
      invisible(NULL)
    }
    # This one really does return NULL, so a NULL is not evidence of a dead
    # worker. Everything else parallel_lapply() runs is checked for it.
    attr(plot_one_file, "returns_null") <- TRUE
    invisible(parallel_lapply(files_unique, plot_one_file, n_cores = n_cores))
    cat("Saved marker plots for:", marker, "\n")
  }


  invisible(normalizePath(output_dir))
}

