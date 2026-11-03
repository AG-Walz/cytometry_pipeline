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
# Folder layout of a run. A cohort run writes two thousand files; flat, the
# figure someone wants is lost among the per-sample marker plots. Each kind of
# output has one subfolder, by metadata column where there is one:
#
#   run_info.json, plot_index.json    what was run, and every plot written
#   overview/                         by combination, file, workspace, subgate
#   markers/                          one map per marker, all cells
#   markers/per_sample/<marker>/      the per-file marker plots
#   metadata/<column>/                coloured by an annotation, one per level
#   splits/<column>/                  per-level plots and the panel
#   changes/<column>/                 fold-change maps, frequency and change tables
#   data/                             cells, supercells, counts, decisions, pairs
#   model/                            umap_model.uwot
#   checkpoints/                      resume (see ckpt_save)
#
# File names are unchanged, so a name still says what a plot is wherever it is.
# run_subdir() is the only place the folder names are spelt.
# -----------------------------------------------------------------------------
run_subdir <- function(kind, column = NULL, marker = NULL) {
  part <- function(x) create_safe_filename(as.character(x))
  switch(kind,
         overview   = "overview",
         marker     = "markers",
         marker_file = file.path("markers", "per_sample", part(marker)),
         annot      = file.path("metadata", part(column)),
         split      = file.path("splits", part(column)),
         split_panel = file.path("splits", part(column)),
         diff       = file.path("changes", part(column)),
         table      = file.path("changes", part(column)),
         data       = "data",
         model      = "model",
         "")
}

# Path of a run's output file `name` in `kind`'s folder, created if needed.
run_path <- function(output_dir, kind, name, column = NULL, marker = NULL) {
  d <- file.path(output_dir, run_subdir(kind, column, marker))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, name)
}

# Where an EXISTING run keeps `name`: its folder in the layout, else the run
# folder itself (runs before 2026-10-15 wrote everything flat), else NA.
run_find <- function(run_dir, kind, name, column = NULL) {
  for (p in c(file.path(run_dir, run_subdir(kind, column), name),
              file.path(run_dir, name))) {
    if (file.exists(p)) return(p)
  }
  NA_character_
}

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
# Samples whose OWN gating tree lacks the resolved population.
#
# gs_get_pop_paths(gs) reads one sample's tree, so a population can resolve for
# a workspace and still be absent from some of its samples: an "All Samples"
# group happily holds a sample that was never gated, or was gated with another
# template. Left alone, extraction stops at that sample with flowWorkspace's
# bare "<path> not found!", minutes into the run and without saying which
# sample. Checked up front instead, so every such sample is named at once.
#
# `paths_by_sample` is a named list (sample name -> its full pop paths).
# Returns the names of the samples missing `population`, in input order.
# -----------------------------------------------------------------------------
samples_missing_population <- function(paths_by_sample, population) {
  has <- vapply(paths_by_sample, function(p) population %in% p, logical(1))
  names(paths_by_sample)[!has]
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

# The marker colour scale a run uses. "same" (the default) follows the scale
# the embedding was computed on, so the colours show exactly the values the
# map was built from: FlowJo values for a plain run, z-scores for a z-scored
# one. Any explicit method overrides it and only recolours.
resolve_colour_scale <- function(normalize_umap, normalize_data) {
  if (is.null(normalize_umap) || identical(normalize_umap, "same")) {
    normalize_data
  } else {
    normalize_umap
  }
}

# Colour values for a re-plot from what a run saved.
#
# `values` are the value_<marker> columns of cells.csv.gz on scale
# `values_scale`: "none" (FlowJo values; runs since 2026-11-03 always save
# these) or, for older runs, the scale they were embedded or coloured on.
# z-score and min-max are both affine per marker, so either can be computed
# from the other exactly; only the FlowJo values cannot be recovered from a
# normalised copy, which is said and the saved scale kept.
#
# @return list(matrix, scale): the colour matrix and the scale it is on.
colour_values_for <- function(values, values_scale, target) {
  if (identical(target, values_scale)) return(list(matrix = values, scale = target))
  if (identical(target, "none")) {
    warning("this run saved its marker values already normalised (",
            values_scale, "), so they cannot be shown as FlowJo values; ",
            "colouring by ", values_scale, " instead.", call. = FALSE)
    return(list(matrix = values, scale = values_scale))
  }
  list(matrix = normalize_columns(values, target), scale = target)
}

# DRAW ORDER. ggplot draws rows in order, so whatever comes last is on top.
# Cells used to be drawn sample by sample, so each sample's cells covered the
# ones before -- and a population of a few hundred cells (CD4/IL2 in J445: 326
# of 2.86M) was painted over by the common cells of every later sample.
#
# Categories: largest first, rarest last (on top); catch-alls (`bottom`, e.g.
# "(none)") under everything; shuffled within a category so no sample sits on
# top of the others. Seeded, so a re-plot draws the same picture.
draw_order_by_category <- function(x, bottom = character(0), seed = 42L) {
  x <- as.character(x)
  freq <- table(x, useNA = "no")
  rank <- match(x, names(sort(freq, decreasing = TRUE)))
  rank[is.na(x) | x %in% bottom] <- 0L
  order(rank, .draw_shuffle(length(x), seed))
}

# Continuous values (marker colours): low first, high on top, so the few
# positive cells of a mostly-negative marker are not buried. NA first.
draw_order_by_value <- function(v) order(!is.na(v), v)

# A random order (by sample, by workspace): nobody's cells on top by rank.
draw_order_shuffled <- function(n, seed = 42L) order(.draw_shuffle(n, seed))

# GATE SEPARATION: a middle ground between the markers layout (one continuum,
# rare gated populations smeared along the big ones) and the gates layout
# (every gate combination its own island).
#
# Each cell's gate label becomes one indicator column, multiplied by `weight`
# times the typical spread of the embedded markers (the median of their MADs),
# and appended to what UMAP embeds. Two cells in different gates are then that
# much further apart; cells in the same gate are unaffected. 0 is the plain
# markers layout; tested on J445 (300k cells, 2026-11-03): 0.5-1 detaches the
# gated populations into compact clusters that stay next to their parents,
# and 2 or more is already the gate islands. A continuous knob, which uwot's
# label supervision is not: there, any weight from 0.02 up gave full islands.
#
# @return The indicator columns (cells x labels), named "gate:<label>", or
#   NULL when weight is 0 or there is only one label.
gate_separation_columns <- function(labels, embed_matrix, weight) {
  weight <- as.numeric(weight)
  if (!length(weight) || !is.finite(weight) || weight <= 0) return(NULL)
  labels <- as.character(labels)
  labels[is.na(labels)] <- "(none)"
  levs <- sort(unique(labels))
  if (length(levs) < 2L) return(NULL)
  spread <- stats::median(apply(embed_matrix, 2, stats::mad), na.rm = TRUE)
  if (!is.finite(spread) || spread <= 0) {
    spread <- stats::median(apply(embed_matrix, 2, stats::sd), na.rm = TRUE)
  }
  if (!is.finite(spread) || spread <= 0) spread <- 1
  step <- weight * spread
  out <- matrix(0, nrow = length(labels), ncol = length(levs),
                dimnames = list(NULL, paste0("gate:", levs)))
  out[cbind(seq_along(labels), match(labels, levs))] <- step
  out
}

.draw_shuffle <- function(n, seed) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  stats::runif(n)
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
#   sample identifier (default: the first of sample_id, sample, sampleID,
#   sample_name, filename, file; else the first column).
# @param markers Optional character vector subset of marker names (NULL = all
#   common non-scatter markers).
# @param subpop_mode One of 'furthest', 'cd4cd8', 'direct', 'deepest',
#   'combinations' (markers layout: colour by the combination of gate calls,
#   as the gates layout labels cells; the embedding is not affected), 'none'.
# @param subpop_select Optional character vector of gates (relative path or bare
#   name) used as the colour partition, overriding subpop_mode.
# @param normalize One of 'none', 'zscore', 'minmax'.
# @param downsample Max cells per sample AFTER gating (0 = no downsampling).
# @param n_neighbors,min_dist,metric,seed UMAP parameters / RNG seed.
# @param point_size,point_alpha ggplot point aesthetics. NULL: 0.4 and 0.6 in
#   the marker, file and annotation plots; scaled to the number of cells in
#   the combination, split, fold-change and conjoint marker plots.
# @param plot_width,plot_height Base plot dimensions (inches).
# @param output_formats Named logical list, e.g. list(pdf=,png=,svg=).
# @param list_populations If TRUE, print the gating tree and return it.
# @param layout "markers" (embed the FlowJo values) or "gates" (embed
#   gate-anchored features so every combination of gate calls gets its own
#   island; see docs/plans/umap_gate_layout_plan.md).
# @param gate_gap,gate_tier_ratio,gate_softness Jump at each gate cut of the
#   deepest tier, how much larger it is per tier up, and the tanh softness of
#   the gate-anchored features (gates layout).
# @param viability_pattern Regex of marker names never used as decisions
#   (NULL: VIABILITY_MARKER_PATTERN; "" turns it off).
# @param supervise uwot target_weight with the combination as label; 0 = off.
# @param cells_per_combination,min_combination_cells Per-sample cap of kept
#   cells per gate-call pattern; combinations with fewer cells in total are
#   pooled as "other (rare)" (gates layout).
# @param filter Samples to keep, "col=a,b;col2!=c" on the sample metadata
#   (either layout).
# @param split_by,level_order,split_reference Metadata columns to draw one
#   plot per level of, their level order ("col=a,b;col2=..."), and the
#   reference level of each for the fold-change maps ("col=level;...").
# @param fit_subset Markers layout: what UMAP is fitted on before every cell
#   is projected. "random" (umap_subsample random cells; 0 = fit on all; a
#   negative value = AUTO_FIT_CELLS when there are more cells, else all) or
#   "supercell" (SuperCellCyto supercell means, see umap_supercells.R).
# @param supercell_gamma,supercell_k_knn SuperCellCyto's cells per supercell
#   and kNN size (fit_subset = "supercell").
# @param background,pair_by Background subtraction: "column=level" of the
#   background samples (e.g. "stimulation=unstimulated") and the metadata
#   columns a sample shares with its background (e.g. c("patient", "visit")).
#   Frequencies become sample minus background, negatives 0.
# @param umap_model Path of a saved model (umap_model.uwot of an earlier run):
#   the fit is skipped and every cell is projected onto that map.
# @param resume Reuse the steps a previous run into output_dir saved
#   (checkpoints/: extraction, supercells, model, coordinates) whose settings
#   match this run's; FALSE recomputes everything.
# @param supercell_max Most supercells UMAP is fitted on: gamma is raised to
#   cells / supercell_max when the cohort would give more (fit_subset =
#   "supercell"). Keeps the fit small enough for the spectral start.
# @return Invisibly, a list with umap_coords, umap_df, workspaces, dropped,
#   output_dir (or, in list_populations mode, what extract_flowjo_population()
#   returned for the inspection).
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
                             normalize_data = "none", normalize_umap = "same",
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
                             report_json = NULL,
                             layout = "markers",
                             sample_filter = NULL,
                             cells_per_combination = 2000L,
                             min_combination_cells = 50L,
                             viability_pattern = VIABILITY_MARKER_PATTERN) {
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
    valid <- if (nm == "normalize_umap") c("same", VALID_NORMALIZE) else VALID_NORMALIZE
    if (!v %in% valid) {
      stop(nm, " must be one of: ", paste(valid, collapse = ", "),
           " (got '", v, "').")
    }
  }
  VALID_SUBPOP_MODE <- c("furthest", "cd4cd8", "direct", "deepest",
                         "combinations", "none")
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

  # LAYOUT. "markers" is the plain embedding of FlowJo values. "gates" reads
  # the gates below the population, labels every cell by its gate calls and
  # keeps the per-cell calls, cuts and values the gate-anchored features are
  # built from (src/visualization/umap_gate_layout.R and
  # docs/plans/umap_gate_layout_plan.md). Everything the markers layout does is
  # left exactly as it was; the gates layout branches off at extraction.
  VALID_LAYOUT <- c("markers", "gates")
  if (!layout %in% VALID_LAYOUT) {
    stop("layout must be one of: ", paste(VALID_LAYOUT, collapse = ", "),
         " (got '", layout, "').")
  }
  gates_layout <- identical(layout, "gates")
  # The markers layout coloured by gate combinations: the cells are labelled
  # exactly as the gates layout labels them (pool_gate_layout()), but the
  # embedding still sees only the marker values. --subpop_select overrides it,
  # as it overrides every other mode.
  combo_colour <- !gates_layout && identical(subpop_mode, "combinations") &&
    is.null(subpop_select)
  if (identical(subpop_mode, "combinations") && !combo_colour) {
    subpop_mode <- "none"
  }
  gate_labels <- gates_layout || combo_colour
  if (gate_labels && !exists("gate_layout_extract_sample", mode = "function")) {
    stop(if (gates_layout) "layout = \"gates\"" else "subpop_mode = \"combinations\"",
         " needs src/visualization/umap_gate_layout.R to be sourced first.",
         call. = FALSE)
  }

  # SAMPLE FILTER, as parsed clauses or the raw "col=a,b;col2!=c" string. It
  # selects samples by their annotations, so it is meaningless without a
  # metadata table or a name pattern, and that is said before the slow load.
  filter_clauses <- sample_filter
  if (is.character(filter_clauses)) {
    if (!exists("parse_sample_filter", mode = "function")) {
      stop("a sample filter needs src/visualization/umap_gate_layout.R to be ",
           "sourced first.", call. = FALSE)
    }
    filter_clauses <- parse_sample_filter(filter_clauses)
  }
  if (is.null(filter_clauses)) filter_clauses <- list()
  if (length(filter_clauses) && is.null(sample_metadata) &&
      is.null(sample_pattern)) {
    stop("--filter selects samples by their metadata, but neither ",
         "--sample_metadata nor --sample_pattern was given.", call. = FALSE)
  }

  if (combo_colour) {
    # Labelled by gate calls below; the subgate candidates are not needed.
    subpop_mode   <- "none"
    subpop_levels <- FALSE
  }
  if (gates_layout) {
    # Cells are labelled by their gate calls here, so the subgate colouring of
    # the markers layout is switched off rather than computed and ignored.
    subpop_mode   <- "none"
    subpop_select <- NULL
    subpop_levels <- FALSE
    # Uniform downsampling is what loses the rare cytokine-positive cells; the
    # gates layout samples per call pattern instead (cells_per_combination).
    if (downsample > 0) {
      warning("--downsample is ignored in the gates layout; cells are sampled ",
              "per gate-call pattern instead (--cells_per_combination = ",
              cells_per_combination, ").", call. = FALSE)
      downsample <- 0L
    }
    # The embedding is computed from gate-anchored features, not from these
    # values, so normalising them for the embedding would do nothing.
    if (!identical(normalize_data, "none")) {
      warning("--normalize_data is ignored in the gates layout: the embedding ",
              "is built from gate-anchored features. --normalize_umap still ",
              "applies to the marker colours.", call. = FALSE)
      normalize_data <- "none"
    }
  }
  # "same": colour by the scale the embedding was computed on. Resolved after
  # the gates layout has set its embedding scale, so it follows that too.
  normalize_umap <- resolve_colour_scale(normalize_umap, normalize_data)

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

  # Which samples the run keeps, by their annotations. Decided here, before a
  # single event is read, so a filtered-out sample costs no extraction time.
  # The filter sees every metadata column, not only the --annotate subset: a
  # run may well filter on a column it does not want plotted.
  sample_keep <- stats::setNames(rep(TRUE, nrow(sample_index)),
                                 sample_index$display)
  if (length(filter_clauses) && !isTRUE(list_populations)) {
    filter_ann <- sample_annotations
    if (!is.null(annotate) && length(annotate)) {
      filter_ann <- suppressWarnings(resolve_sample_annotations(
        sample_index, metadata_file = sample_metadata, pattern = sample_pattern,
        annotate = NULL, id_column = sample_id_column))
    }
    sample_keep[] <- apply_sample_filter(filter_ann, filter_clauses)
    filter_txt <- paste(vapply(filter_clauses, function(cl)
      paste0(cl$column, if (cl$negate) "!=" else "=",
             paste(cl$values, collapse = ",")), character(1)), collapse = "; ")
    cat(sprintf("\nSample filter (%s): %d of %d sample(s) kept, %d left out.\n",
                filter_txt, sum(sample_keep), length(sample_keep),
                sum(!sample_keep)))
    if (!any(sample_keep)) {
      stop("the sample filter (", filter_txt, ") leaves no sample. Check the ",
           "values against the metadata levels listed above.", call. = FALSE)
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

  # Every sample's own tree, per gating set: the population has to be in each
  # of them, and 2b needs their union.
  paths_by_sample <- lapply(gs_list, function(gs) {
    nms <- flowWorkspace::sampleNames(gs)
    stats::setNames(lapply(nms, function(s)
      flowWorkspace::gh_get_pop_paths(gs[[s]], path = "full")), nms)
  })
  # A sample the --filter leaves out is never read, so its tree does not matter.
  lacking <- lapply(names(gs_list), function(w) {
    miss <- samples_missing_population(paths_by_sample[[w]], pop_by_wsp[[w]])
    miss[vapply(miss, function(s)
      sample_keep[[sample_display_label(sample_index, gs_source, w, s)]],
      logical(1))]
  })
  names(lacking) <- names(gs_list)
  lacking <- lacking[lengths(lacking) > 0]
  # Left out, loudly, rather than stopping the run: one stray sample in an "All
  # Samples" group is the usual case, and the rest of the run is sound without
  # it. Its cells must not be pooled under some other gate, so it is excluded,
  # named here and recorded as samples_left_out in run_info.json.
  if (length(lacking)) {
    detail <- vapply(names(lacking), function(w) {
      rows <- vapply(lacking[[w]], function(s) sprintf(
        "      %s\n        its gates: %s",
        sub("(\\.fcs)_[0-9]+$", "\\1", s, ignore.case = TRUE),
        paste(utils::head(setdiff(paths_by_sample[[w]][[s]], "root"), 6),
              collapse = ", ")), character(1))
      sprintf("  %s (%s): %d sample(s) without %s\n%s", w,
              specs[[gs_source[[w]]]]$wsp,
              length(lacking[[w]]), pop_by_wsp[[w]],
              paste(rows, collapse = "\n"))
    }, character(1))
    cat("\nWARNING: ", sum(lengths(lacking)), " sample(s) are LEFT OUT: they ",
        "do not have the population in their own gating tree, although the ",
        "rest of their workspace does. Usually a sample that was never gated, ",
        "or gated with a different template, inside a generic group such as ",
        "'All Samples'.\n\n", paste(detail, collapse = "\n"),
        "\n\nTo include them, gate them in FlowJo. To silence this, remove ",
        "them from the workspace or load a sample group that leaves them ",
        "out.\n", sep = "")
    for (w in names(lacking)) for (s in lacking[[w]]) {
      sample_keep[[sample_display_label(sample_index, gs_source, w, s)]] <- FALSE
    }
    if (!any(sample_keep)) {
      stop("no sample is left: every one lacks the population in its own ",
           "gating tree.", call. = FALSE)
    }
  }

  # --- 2b. Discover the FlowJo subgates UNDER `population` to colour by. ------
  # Per workspace, then merged on the RELATIVE path (the only key comparable
  # across workspaces whose ancestors are named differently).
  cands <- lapply(names(gs_list), function(w) {
    pop_paths_union <- unique(unlist(paths_by_sample[[w]]))
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
  # Gates layout only: per sample, the kept cells' calls, cuts and values plus
  # the call-pattern counts over ALL cells, and the gate decisions.
  per_sample_gl     <- list()
  per_sample_dec    <- list()
  # Progress is reported against the sample count the GatingSets actually hold,
  # which is the only total that is known before the loop and true after it.
  n_total_samples <- sum(vapply(gs_list, function(g)
    length(flowWorkspace::sampleNames(g)), integer(1)))
  if (!all(sample_keep)) n_total_samples <- sum(sample_keep)
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
      if (!sample_keep[[display]]) next   # left out (filter, or no gate)
      if (gate_labels) {
        res <- gate_layout_extract_sample(gs[[s]], pop_by_wsp[[w]],
                                          dictionary = channel_dict$dictionary,
                                          aliases = marker_alias_map,
                                          exclude_markers = viability_pattern)
      } else {
        res <- extract_population_sample(gs[[s]], pop_by_wsp[[w]],
                                         subpop_paths_by_wsp[[w]],
                                         subpop_depths, subpop_labels,
                                         assignment = subpop_assignment,
                                         dictionary = channel_dict$dictionary,
                                         aliases = marker_alias_map,
                                         levels = isTRUE(subpop_levels),
                                         level_spec = level_spec)
      }
      # No usable channel on this sample. Collected rather than fatal — see
      # extract_population_sample() — and reported together below.
      if (is.null(res)) {
        no_marker_samples <- c(no_marker_samples, display)
        cat(sprintf("  %s: SKIPPED — no named marker channel\n", display))
        next
      }
      mtx <- res$mtx
      sp  <- res$subpop

      # Gates layout: keep up to cells_per_combination cells of every gate-call
      # pattern (rare ones whole), each carrying the weight that restores its
      # sample's true count, and count every pattern over ALL cells first so
      # the frequencies reported later never depend on the sampling.
      # Markers layout coloured by combinations: every cell (or the uniform
      # --downsample of them), each weighing 1, so the embedding sees the
      # cohort as it is.
      if (gate_labels) {
        ss <- if (gates_layout) {
          stratified_sample(res$key, cells_per_combination)
        } else {
          n_res <- length(res$key)
          keep_rows <- if (downsample > 0 && n_res > downsample) {
            sort(sample(n_res, downsample))
          } else {
            seq_len(n_res)
          }
          list(idx = keep_rows, weight = rep(1, length(keep_rows)))
        }
        kc <- table(res$key)
        per_sample_gl[[display]] <- list(
          calls  = res$calls[ss$idx, , drop = FALSE],
          cuts   = res$cuts[ss$idx, , drop = FALSE],
          anchors = res$anchors[ss$idx, , drop = FALSE],
          values = res$values[ss$idx, , drop = FALSE],
          mtx    = res$mtx[ss$idx, , drop = FALSE],
          weight = ss$weight, key = res$key[ss$idx],
          key_counts = stats::setNames(as.integer(kc), names(kc)),
          population_cells = nrow(res$mtx))
        if (nrow(res$decisions)) {
          d <- res$decisions
          d$sample <- display
          d$workspace <- gs_source[[w]]
          per_sample_dec[[display]] <- d
        }
        mtx <- per_sample_gl[[display]]$mtx
        # Placeholder: labels need the whole cohort's tiers; set after pooling.
        sp  <- rep(NA_character_, nrow(mtx))
        gl_note <- sprintf("    gate calls: %d cells in %d pattern(s), %d kept; decided: %s\n",
                           nrow(res$mtx), length(kc), nrow(mtx),
                           if (ncol(res$calls)) paste(colnames(res$calls), collapse = ", ")
                           else "(none)")
      }

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
      if (gate_labels) cat(gl_note)
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

  # --- 4b. Gates layout: pool the samples, label every cell. -----------------
  # Tiers, the fill of missing decisions and the rare-combination pooling are
  # properties of the cohort, so they are decided here, over every sample, and
  # not inside the loop. See pool_gate_layout().
  gl <- NULL
  if (gate_labels) {
    decisions <- do.call(rbind, unname(per_sample_dec[names(per_sample_mtx)]))
    if (is.null(decisions) || !nrow(decisions)) {
      stop("no gate below '", population, "' decides a marker in any sample, so ",
           if (gates_layout) "there is nothing for the gates layout to separate. "
           else "there are no gate combinations to colour by. ",
           "Pick a population with gates under it",
           if (gates_layout) ", or use --layout markers." else
             ", or another --subpop_mode.", call. = FALSE)
    }
    rownames(decisions) <- NULL
    panel_order <- unique(unlist(lapply(per_sample_mtx, colnames),
                                 use.names = FALSE))
    gl <- pool_gate_layout(per_sample_gl[names(per_sample_mtx)], decisions,
                           panel_order = panel_order,
                           min_combination_cells = min_combination_cells)
    dec_markers <- names(gl$tiers)
    subpop_mode_tag <- "gate combinations"

    cat("\nGate decisions under ", population, " (marker: tier, samples decided, median cut):\n",
        sep = "")
    for (m in dec_markers) {
      dm <- decisions[decisions$marker == m, , drop = FALSE]
      cat(sprintf("  %-10s tier %d, %d of %d sample(s), cut %.1f\n", m,
                  gl$tiers[[m]], length(unique(dm$sample)),
                  length(per_sample_mtx), stats::median(dm$cut)))
    }
    if (!is.null(gl$filled) && nrow(gl$filled)) {
      cat(sprintf(paste0("  %d sample(s) lack a decision most of their lineage has; ",
                         "those cells are treated as negative for it:\n"),
                  length(unique(gl$filled$sample))))
      fl <- stats::aggregate(cells ~ marker + group, data = gl$filled, FUN = sum)
      # The group is the first-tier call pattern ("10"), spelled out by marker.
      t1 <- names(gl$tiers)[gl$tiers == min(gl$tiers)]
      sign_of <- c("1" = "+", "0" = "-", "x" = "?")
      for (i in seq_len(nrow(fl))) {
        grp <- paste0(t1, sign_of[strsplit(fl$group[i], "")[[1]]], collapse = " ")
        cat(sprintf("    %s in %s cells: %d kept cell(s)\n", fl$marker[i], grp,
                    as.integer(fl$cells[i])))
      }
    }

    if (gates_layout && is.null(markers)) {
      # The features are the markers a gate decides, as long as every sample
      # carries the channel: a sample without a marker's DECISION still has its
      # values (pool_gate_layout() takes them from the channel), a sample
      # without the CHANNEL has nothing to embed.
      use_markers <- dec_markers[dec_markers %in% common_markers]
      not_shared <- setdiff(dec_markers, common_markers)
      if (length(not_shared)) {
        cat("  Decided but not measured in every sample, so labels only, not features: ",
            paste(not_shared, collapse = ", "), "\n", sep = "")
      }
      if (!length(use_markers)) {
        stop("none of the decided markers (", paste(dec_markers, collapse = ", "),
             ") is measured in every sample.", call. = FALSE)
      }
    }
    extra <- setdiff(use_markers, dec_markers)
    if (gates_layout) {
      cat(sprintf("  Gate-anchored features: %s%s\n",
                  paste(intersect(use_markers, dec_markers), collapse = ", "),
                  if (length(extra)) paste0("; without a gate: ",
                                            paste(extra, collapse = ", ")) else ""))
    } else {
      cat("  Gate combinations colour the plots only; the embedding uses the ",
          "marker values.\n", sep = "")
    }

    tot <- tapply(gl$weight, gl$combination, sum)
    tot <- sort(tot, decreasing = TRUE)
    cat(sprintf("  %d combination(s) of gate calls (combinations under %d cells pooled as '%s'):\n",
                length(tot), as.integer(min_combination_cells), .GL_RARE))
    for (k in names(tot)) {
      cat(sprintf("    %-40s %10.0f cells, %d kept\n", k, tot[[k]],
                  sum(gl$combination == k)))
    }
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
  # Gates layout (and the markers layout coloured by combinations): the label
  # of a cell is its combination of gate calls, pooled in the same sample
  # order by pool_gate_layout().
  if (gate_labels) {
    stopifnot(length(gl$combination) == nrow(combined_matrix))
    combined_subpop <- gl$combination
  }

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
  # The FlowJo values themselves, kept for cells.csv.gz: a re-plot can then
  # colour by any scale exactly. The same object as combined_matrix unless the
  # embedding is normalised.
  raw_matrix <- combined_matrix
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
    raw_matrix      = raw_matrix,
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
    inspected       = FALSE,
    layout          = layout,
    # Per sample (display label): kept by the sample filter or not.
    sample_keep     = sample_keep,
    # Gates layout only (NULL otherwise): calls, values, cuts, weight,
    # combination per kept cell (row-aligned to combined_matrix), decisions,
    # counts over all cells, tiers, filled. See pool_gate_layout().
    gate_layout     = gl
  )
}


# -----------------------------------------------------------------------------
# Helpers of the gates layout and the metadata splits (umap_from_flowjo()).
# -----------------------------------------------------------------------------

# Annotation column names a metadata table or a sample-name pattern will give,
# read without loading a workspace: the table's header minus its id column, or
# the pattern's named groups. This is what lets --split_by and --filter be
# checked before the slow part of a run.
annotation_columns_offered <- function(sample_metadata = NULL,
                                       sample_pattern = NULL,
                                       sample_id_column = NULL) {
  if (!is.null(sample_metadata)) {
    md <- read_sample_metadata(sample_metadata)
    return(setdiff(colnames(md), metadata_id_column(md, sample_id_column)))
  }
  if (!is.null(sample_pattern)) {
    nms <- attr(regexpr(sample_pattern, "", perl = TRUE), "capture.names")
    return(nms[nzchar(nms)])
  }
  character(0)
}

# A per-cell table as gzipped CSV. data.table::fwrite is many times faster
# than write.csv on millions of rows; base R is the fallback when data.table
# (or its zlib support) is missing.
write_table_gz <- function(df, path) {
  done <- requireNamespace("data.table", quietly = TRUE) &&
    isTRUE(tryCatch({
      data.table::fwrite(df, path, compress = "gzip")
      TRUE
    }, error = function(e) FALSE))
  if (!done) {
    con <- gzfile(path, "w")
    on.exit(close(con))
    utils::write.csv(df, con, row.names = FALSE)
  }
  invisible(path)
}

# run_info.json of a run: how to reproduce it (argv, parameters, seed,
# threads, code version, package versions) and what it found (cell counts,
# labels; for the gates layout the decisions and feature scales, for a
# supercell fit set the supercell count). Vectors that must stay JSON arrays
# even with one element are wrapped in I(), because auto_unbox would turn them
# into scalars. The gate-only arguments are NULL in the markers layout, where
# gate_markers holds the embedded markers.
write_gate_run_info <- function(path, parameters, seed, embed_threads,
                                sigma = NULL, gate_markers, extra_markers = NULL,
                                tiers = NULL, decisions = NULL, filled = NULL,
                                n_cells, n_fit, n_samples,
                                samples_left_out, combinations,
                                supercells = NULL, embedding = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("run_info.json not written: install jsonlite.")
    return(invisible(NULL))
  }
  repo <- if (exists("src_dir")) dirname(normalizePath(src_dir, mustWork = FALSE))
          else getwd()
  git <- function(args) tryCatch(
    suppressWarnings(system2("git", c("-C", repo, args), stdout = TRUE,
                             stderr = FALSE)),
    error = function(e) structure(character(0), status = 1L))
  commit <- git(c("rev-parse", "HEAD"))
  in_git <- length(commit) == 1L && is.null(attr(commit, "status"))
  pkg <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                              error = function(e) NA_character_)
  info <- list(
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    argv = I(commandArgs()),
    R = R.version.string,
    packages = lapply(stats::setNames(nm = c("uwot", "ggplot2", "CytoML",
                                             "flowWorkspace", "ggrastr", "RANN",
                                             "data.table", "jsonlite",
                                             "SuperCellCyto")), pkg),
    git = list(repo = repo,
               commit = if (in_git) commit else NA_character_,
               dirty = if (in_git) length(git(c("status", "--porcelain"))) > 0
                       else NA),
    seed = seed,
    embed_threads = embed_threads,
    parameters = parameters,
    features = if (is.null(tiers)) list(markers = I(gate_markers))
               else list(gate_anchored = I(gate_markers),
                         without_gate = I(extra_markers)),
    sigma = if (is.null(sigma)) NULL else as.list(sigma),
    decisions = if (is.null(tiers) || is.null(decisions)) NULL
                else lapply(names(tiers), function(m) {
      d <- decisions[decisions$marker == m, , drop = FALSE]
      list(marker = m, tier = tiers[[m]], samples = length(unique(d$sample)),
           median_cut = stats::median(d$cut),
           cut_range = I(range(d$cut)))
    }),
    filled = if (is.null(filled)) list() else filled,
    cells = n_cells, fit_cells = n_fit, samples = n_samples,
    samples_left_out = I(samples_left_out),
    combinations = as.list(combinations),
    supercells = supercells,
    # embedding_check(): whether the map is fit to look at, and why not.
    embedding = embedding)
  # The full settings (every option, inputs, host) live in the run's settings
  # record, written when the run started; see R_shared/run_record.R.
  settings <- file.path(dirname(path), "_run", "settings.json")
  if (file.exists(settings)) info$settings <- "_run/settings.json"
  info <- info[!vapply(info, is.null, logical(1))]
  jsonlite::write_json(info, path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", na = "null", digits = NA)
  invisible(path)
}

# The plots plot_flowjo_umap_facets() wrote in this run, as index entries:
# files modified since `since` (whole seconds, so a coarse file system clock
# cannot drop one), minus the gate-layout plots already indexed.
facet_plot_entries <- function(output_dir, since, exclude = NULL) {
  f <- list.files(output_dir, pattern = "\\.(png|pdf|svg)$", recursive = TRUE)
  f <- f[!startsWith(f, "checkpoints/")]
  f <- f[as.numeric(file.mtime(file.path(output_dir, f))) >= floor(as.numeric(since))]
  rel <- unique(sub("\\.(png|pdf|svg)$", "", f))
  # `exclude` (the styled plots' index) gives paths; compare by file stem.
  rel <- rel[!basename(rel) %in% basename(as.character(exclude))]
  stems <- rel
  if (!length(stems)) return(NULL)
  kind <- ifelse(grepl("_marker_", stems, fixed = TRUE), "marker",
                 ifelse(grepl("_annot_", stems, fixed = TRUE), "annot",
                        "overview"))
  data.frame(file = stems, kind = kind, stringsAsFactors = FALSE)
}

# plot_index.json: one entry per plot (file stem, kind, column, level,
# reference), which the app uses to sort the gallery instead of guessing from
# file names.
write_plot_index <- function(path, gate_index, facet_entries) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("plot_index.json not written: install jsonlite.")
    return(invisible(NULL))
  }
  cols <- c("file", "kind", "column", "level", "reference")
  norm <- function(df) {
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    for (cl in setdiff(cols, names(df))) df[[cl]] <- NA_character_
    df[, cols, drop = FALSE]
  }
  idx <- rbind(norm(gate_index), norm(facet_entries))
  jsonlite::write_json(if (is.null(idx)) list() else idx, path,
                       dataframe = "rows", pretty = TRUE, na = "null")
  invisible(path)
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
                             normalize_data = "none", normalize_umap = "same",
                             subpop_levels = TRUE,
                             subpop_levels_select = NULL,
                             downsample = 0L, n_neighbors = 15L, min_dist = NULL,
                             metric = "euclidean", seed = 42L,
                             point_size = NULL, point_alpha = NULL,
                             plot_width = 10, plot_height = 8,
                             output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
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
                             umap_subsample = -1L,
                             layout = "markers",
                             gate_gap = 2,
                             gate_tier_ratio = 2,
                             gate_softness = 4,
                             supervise = 0,
                             cells_per_combination = 2000L,
                             min_combination_cells = 50L,
                             viability_pattern = NULL,
                             filter = NULL,
                             split_by = NULL,
                             level_order = NULL,
                             split_reference = NULL,
                             fit_subset = "random",
                             supercell_gamma = 20,
                             supercell_k_knn = 5L,
                             supercell_max = 100000L,
                             resume = TRUE,
                             background = NULL,
                             pair_by = NULL,
                             umap_model = NULL,
                             nn_method = NULL,
                             gate_separation = 0) {
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

  # --- Layout, metadata splits and sample filter: checked BEFORE the load. ----
  # A misspelt column should cost seconds, not the hour a cohort takes to load.
  # Column names are read from the metadata header (or the pattern's group
  # names), which is exactly what resolve_sample_annotations() will offer.
  gates_layout <- identical(layout, "gates")
  # min_dist depends on the layout when it is not given: the gate layout's
  # islands are rounder and closer to the look asked for at 0.25 (tuned on the
  # MTB cohort, docs/plans/umap_gate_layout_plan.md section 11); the markers
  # layout keeps the value it always had.
  if (is.null(min_dist) || is.na(min_dist)) {
    min_dist <- if (gates_layout) 0.25 else 0.1
  }
  if (gates_layout && (!is.finite(gate_tier_ratio) || gate_tier_ratio < 1)) {
    stop("gate_tier_ratio must be 1 or more: below 1 a lineage gate's jump ",
         "would be smaller than the cytokine gates' and could fall under the ",
         "within-side spread, so its two sides could overlap.", call. = FALSE)
  }
  if (gates_layout && gate_gap <= 1) {
    warning("gate_gap <= 1: the two sides of a gate are no longer guaranteed ",
            "to be apart (the within-side part of a feature spans -1 to 1).",
            call. = FALSE)
  }
  if (gates_layout && !identical(tolower(method), "umap")) {
    stop("layout = \"gates\" embeds with UMAP only (it fits on a set balanced ",
         "across gate combinations and projects the rest); method = '", method,
         "' is not supported there.", call. = FALSE)
  }
  # What the markers layout fits UMAP on before projecting every cell. The
  # gates layout has its own balanced fit set, so a supercell choice there is
  # said and dropped rather than silently ignored.
  if (!fit_subset %in% c("random", "supercell")) {
    stop("fit_subset must be 'random' or 'supercell' (got '", fit_subset, "').",
         call. = FALSE)
  }
  if (gates_layout && identical(fit_subset, "supercell")) {
    warning("--fit_subset supercell is ignored in the gates layout, which fits ",
            "on a set balanced across gate combinations.", call. = FALSE)
    fit_subset <- "random"
  }
  gate_separation <- suppressWarnings(as.numeric(gate_separation))
  if (!length(gate_separation) || is.na(gate_separation) || gate_separation < 0) {
    stop("--gate_separation must be a number of 0 or more.", call. = FALSE)
  }
  if (gate_separation > 0 && gates_layout) {
    warning("--gate_separation is ignored with --layout gates, which separates ",
            "the gates already.", call. = FALSE)
    gate_separation <- 0
  }
  if (gate_separation > 0 && !is.null(umap_model) && nzchar(umap_model)) {
    stop("--gate_separation cannot be combined with --umap_model: the saved ",
         "map was fitted on other columns.", call. = FALSE)
  }
  # An existing map: no fit at all, every cell projected onto the saved model.
  # The gates layout's features depend on the cohort they were built from, so
  # its model cannot take another cohort's cells.
  model_file <- NULL
  if (!is.null(umap_model) && nzchar(umap_model)) {
    if (gates_layout || !identical(tolower(method), "umap")) {
      stop("--umap_model projects onto a saved markers-layout UMAP model; it ",
           "does not apply to the gates layout or to method '", method, "'.",
           call. = FALSE)
    }
    if (!file.exists(umap_model)) {
      stop("--umap_model: no such file: ", umap_model, call. = FALSE)
    }
    model_file <- normalizePath(umap_model)
    if (identical(fit_subset, "supercell")) {
      message("--fit_subset supercell is not needed with --umap_model: ",
              "nothing is fitted.")
      fit_subset <- "random"
    }
  }
  supercell_fit <- identical(fit_subset, "supercell")
  if (supercell_fit) {
    if (!identical(tolower(method), "umap")) {
      stop("--fit_subset supercell projects cells onto a UMAP model; method = '",
           method, "' cannot project.", call. = FALSE)
    }
    if (!exists("require_supercellcyto", mode = "function")) {
      stop("--fit_subset supercell needs src/visualization/umap_supercells.R to ",
           "be sourced first.", call. = FALSE)
    }
    if (!isTRUE(list_populations)) require_supercellcyto()
    if (!is.finite(supercell_gamma) || supercell_gamma <= 1) {
      stop("supercell_gamma must be above 1 (cells per supercell).", call. = FALSE)
    }
    if (is.na(as.integer(supercell_k_knn)) || as.integer(supercell_k_knn) < 1L) {
      stop("supercell_k_knn must be a whole number of at least 1.", call. = FALSE)
    }
    if (is.na(as.integer(supercell_max)) || as.integer(supercell_max) < 1000L) {
      stop("supercell_max must be a whole number of at least 1000.", call. = FALSE)
    }
    if (umap_subsample > 0) {
      warning("--umap_subsample is ignored with --fit_subset supercell: the ",
              "supercells are the fit set.", call. = FALSE)
    }
  }
  split_cols <- unique(trimws(as.character(split_by)))
  split_cols <- split_cols[nzchar(split_cols)]
  level_order_list <- list()
  split_ref_list   <- list()
  bg_spec          <- NULL
  filter_clauses   <- filter
  if (!isTRUE(list_populations)) {
    uses_gl <- gates_layout || length(split_cols) || !is.null(filter) ||
      !is.null(level_order) || !is.null(split_reference)
    if (uses_gl && !exists("parse_level_order", mode = "function")) {
      stop("--layout gates, --split_by, --filter, --level_order and ",
           "--split_reference need src/visualization/umap_gate_layout.R to be ",
           "sourced first.", call. = FALSE)
    }
    if (is.character(filter_clauses)) {
      filter_clauses <- parse_sample_filter(filter_clauses)
    }
    level_order_list <- if (uses_gl) parse_level_order(level_order) else list()
    split_ref_list   <- if (uses_gl) parse_split_reference(split_reference) else list()
    if (length(split_cols) || length(filter_clauses)) {
      if (is.null(sample_metadata) && is.null(sample_pattern)) {
        stop("--split_by and --filter read the sample metadata, but neither ",
             "--sample_metadata nor --sample_pattern was given.", call. = FALSE)
      }
      offered <- annotation_columns_offered(sample_metadata, sample_pattern,
                                            sample_id_column)
      bad <- setdiff(split_cols, offered)
      if (length(bad)) {
        stop("--split_by names column(s) ", paste(bad, collapse = ", "),
             ", which the sample annotations do not have. Available: ",
             if (length(offered)) paste(offered, collapse = ", ") else "(none)",
             ".", call. = FALSE)
      }
      # A split column needs its values per cell, so it joins the --annotate
      # selection rather than being refused by it: --annotate chooses which
      # columns get their own annotation plots, not which ones exist.
      if (length(annotate)) annotate <- union(annotate, split_cols)
      bad <- setdiff(vapply(filter_clauses, `[[`, character(1), "column"),
                     offered)
      if (length(bad)) {
        stop("--filter names column(s) ", paste(unique(bad), collapse = ", "),
             ", which the sample metadata does not have. Available: ",
             if (length(offered)) paste(offered, collapse = ", ") else "(none)",
             ".", call. = FALSE)
      }
    }
    bg_spec <- parse_background(background)
    if (!is.null(bg_spec)) {
      if (is.null(sample_metadata) && is.null(sample_pattern)) {
        stop("--background reads the sample metadata, but neither ",
             "--sample_metadata nor --sample_pattern was given.", call. = FALSE)
      }
      pair_by <- unique(trimws(as.character(pair_by)))
      pair_by <- pair_by[nzchar(pair_by)]
      if (!length(pair_by)) {
        stop("--background needs --pair_by: the metadata column(s) a sample ",
             "shares with its background sample, e.g. patient,visit.",
             call. = FALSE)
      }
      offered <- annotation_columns_offered(sample_metadata, sample_pattern,
                                            sample_id_column)
      bad <- setdiff(c(bg_spec$column, pair_by), offered)
      if (length(bad)) {
        stop("--background/--pair_by name column(s) ", paste(bad, collapse = ", "),
             ", which the sample metadata does not have. Available: ",
             paste(offered, collapse = ", "), ".", call. = FALSE)
      }
    }
    for (nm in c("level_order", "split_reference")) {
      given <- names(if (nm == "level_order") level_order_list else split_ref_list)
      stray <- setdiff(given, split_cols)
      if (length(stray)) {
        warning("--", nm, " names column(s) not in --split_by, ignored: ",
                paste(stray, collapse = ", "), call. = FALSE)
      }
    }
  }

  # --- Checkpoints: a cohort run is an hour of extraction, supercells, fit and
  # projection, and each step is saved under output_dir/checkpoints/ as soon as
  # it is done, with a key of the settings it depends on (each key includes
  # the one before it). A rerun into the same folder reuses every step whose
  # key still matches and recomputes from the first one that does not, so a
  # crash, a cancel or a changed plot setting costs only what has to be redone.
  extract_args <- list(
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
    workspaces = workspaces, report_json = report_json,
    layout = layout, sample_filter = filter_clauses,
    viability_pattern = if (is.null(viability_pattern)) VIABILITY_MARKER_PATTERN
                        else viability_pattern,
    cells_per_combination = cells_per_combination,
    min_combination_cells = min_combination_cells)
  use_ckpt <- isTRUE(resume) && !isTRUE(list_populations)
  key_extract <- ckpt_key(list(
    "extract", extract_args[setdiff(names(extract_args), "report_json")],
    inputs = ckpt_input_stamp(wsp_path, fcs_dir, sample_metadata)))
  ext <- if (use_ckpt) ckpt_load(output_dir, "extract", key_extract) else NULL
  if (is.null(ext)) {
    ext <- do.call(extract_flowjo_population, extract_args)
    # Inspect mode already printed and returned everything it has.
    if (isTRUE(ext$inspected)) return(invisible(ext))
    if (use_ckpt) ckpt_save(output_dir, "extract", key_extract, ext)
  }

  combined_matrix     <- ext$combined_matrix
  # The value_<marker> columns of cells.csv.gz: the FlowJo values, so that a
  # re-plot can colour by any scale exactly. An extract checkpoint from before
  # 2026-11-03 has no raw copy; then the values on the scale they were kept on.
  saved_values <- if (!is.null(ext$raw_matrix)) {
    list(matrix = ext$raw_matrix, scale = "none")
  } else if (identical(layout, "gates")) {
    list(matrix = ext$colour_matrix, scale = ext$normalize_umap)
  } else {
    list(matrix = combined_matrix, scale = ext$normalize_data)
  }
  # Each sample's background sample, from every metadata column (not only the
  # --annotate subset), as the filter reads them.
  bg_map <- NULL
  if (!is.null(bg_spec)) {
    full_ann <- suppressWarnings(resolve_sample_annotations(
      ext$sample_index, metadata_file = sample_metadata, pattern = sample_pattern,
      annotate = NULL, id_column = sample_id_column))
    rows <- match(ext$samples, ext$sample_index$display)
    bg_map <- background_pairs(ext$samples, full_ann[rows, , drop = FALSE],
                               bg_spec$column, bg_spec$level, pair_by)
    make_outdir(output_dir)
    utils::write.csv(data.frame(sample = names(bg_map), background = unname(bg_map),
                                stringsAsFactors = FALSE),
                     run_path(output_dir, "data", "background_pairs.csv"), row.names = FALSE)
    cat(sprintf(paste0("Background: %d sample(s) paired with their %s = %s ",
                       "sample by %s (background_pairs.csv)\n"),
                length(bg_map), bg_spec$column, bg_spec$level,
                paste(pair_by, collapse = ", ")))
  }
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

  # Gates layout: the embedding is computed from gate-anchored features, not
  # from combined_matrix (which keeps the raw values for the marker colours).
  # Each decided marker becomes tanh((x - anchor) / (softness * sigma)) plus a
  # jump of +-gap_m on the side of its call (anchor: the median of the
  # sample's cells of the same combination; gap_m larger one tier up), so
  # cells on different sides of any gate cannot overlap; markers no gate
  # decides enter without a jump. The fit set is
  # balanced across combinations, so a rare cytokine combination gets as much
  # say in the layout as the bulk CD4 cells, and the rest are projected.
  gl <- ext$gate_layout
  if (gates_layout) {
    gate_markers  <- intersect(use_markers, colnames(gl$calls))
    extra_markers <- setdiff(use_markers, gate_markers)
    sigma <- gate_scales(gl$values, gl$anchors, gl$calls, gl$weight)
    gaps <- tier_gaps(gl$tiers[gate_markers], gate_gap, gate_tier_ratio)
    features <- gate_features(
      gl$values[, gate_markers, drop = FALSE],
      gl$anchors[, gate_markers, drop = FALSE],
      gl$calls[, gate_markers, drop = FALSE],
      sigma[gate_markers], gap = gaps, softness = gate_softness,
      extra = if (length(extra_markers))
        combined_matrix[, extra_markers, drop = FALSE] else NULL)
    # Cells no gate decides on any marker (with a population above the live
    # gate: the cells that gate excluded) would otherwise sit on the negative
    # side of every marker, on top of the all-negative combination. One more
    # coordinate with the largest jump gives them their own island.
    undecided <- if (ncol(gl$calls)) rowSums(!is.na(gl$calls)) == 0 else
      rep(TRUE, nrow(features))
    if (any(undecided) && !all(undecided)) {
      jump <- if (length(gaps)) max(gaps) else gate_gap
      features <- cbind(features,
                        no_gate_decision = ifelse(undecided, jump, -jump))
      cat(sprintf(paste0("NOTE: %s of %s embedded cells (%.1f%% of the ",
                         "population, weighted) fall outside every gate that ",
                         "decides a marker; they are labelled '(no decision)' ",
                         "and kept apart. With a population above the live ",
                         "gate these are the cells that gate excluded.\n"),
                  format(sum(undecided), big.mark = ","),
                  format(length(undecided), big.mark = ","),
                  100 * sum(gl$weight[undecided]) / sum(gl$weight)))
    }
    if (!all(is.finite(features))) {
      bad_cols <- colnames(features)[
        apply(features, 2, function(x) any(!is.finite(x)))]
      stop("gate-anchored features are not finite for marker(s): ",
           paste(bad_cols, collapse = ", "), ". Refusing to embed.",
           call. = FALSE)
    }
    # The fit set is drawn from its own seed so it does not depend on how many
    # random numbers the extraction happened to use.
    fit_budget <- if (umap_subsample > 0) as.integer(umap_subsample) else 100000L
    set.seed(seed)
    fit_idx <- balanced_fit_index(gl$combination, fit_budget)
    cat(sprintf(paste0("Gate layout: gap %g (x%g per tier up), softness %g, ",
                       "supervision %g; ",
                       "sigma %s\n  fit set: %d of %d cells, balanced over %d ",
                       "combination(s)\n"),
                gate_gap, gate_tier_ratio, gate_softness, supervise,
                paste(sprintf("%s=%.2f", gate_markers, sigma[gate_markers]),
                      collapse = ", "),
                length(fit_idx), nrow(features),
                length(unique(gl$combination))))
  }

  # Through run_dimred(), not run_umap(): the dispatcher is what the rest of the
  # repo embeds with, it re-checks that the coordinate count still matches the
  # matrix, and going through it means `method` can be tsne or diffusion_map
  # here without another line of code. draw_plots = FALSE because the figures
  # below are this path's own — the standard by-sample plot is not one of them.
  # The gates layout calls run_umap() directly: it needs the fit set and the
  # supervision, which only run_umap() takes.
  # The default fit set. A full fit is superlinear in the cell count, and past
  # SPECTRAL_LIMIT (run_umap) it also loses the spectral start; fitting on
  # AUTO_FIT_CELLS random cells and projecting the rest gets the spectral start
  # and finishes in minutes. Only the markers layout's random fit set is
  # affected: the gates layout reads a negative value as its own default
  # budget (above), and supercells are their own fit set. 0 still means every
  # cell, whatever the count.
  AUTO_FIT_CELLS <- 300000L
  if (is.na(umap_subsample) || umap_subsample < 0) {
    umap_subsample <- if (!gates_layout && !supercell_fit &&
                          nrow(combined_matrix) > AUTO_FIT_CELLS) AUTO_FIT_CELLS
                      else 0L
    if (umap_subsample > 0) {
      cat("Fit set: ", format(AUTO_FIT_CELLS, big.mark = ","), " random cells of ",
          format(nrow(combined_matrix), big.mark = ","), "; the rest are ",
          "projected onto them. Pass --umap_subsample 0 to fit on every cell.\n",
          sep = "")
    }
  }
  # Vector formats store every point. Said before the hour of embedding, so the
  # choice can be undone while it is still cheap: the 2026-11-02 J445 run wrote
  # 8.5 GB, a 417 MB SVG and a 100 MB PDF per plot for 2.86M cells.
  vector_formats <- intersect(c("svg", "pdf"),
                              names(Filter(isTRUE, output_formats)))
  if (length(vector_formats) && nrow(combined_matrix) > 500000L) {
    cat("NOTE: ", paste(toupper(vector_formats), collapse = " and "),
        " store one vector point per cell; at ",
        format(nrow(combined_matrix), big.mark = ","), " cells each plot will ",
        "be hundreds of megabytes. PNG alone is usually what is wanted.\n",
        sep = "")
  }
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
  # Said BEFORE the call: the neighbour search runs before the first epoch and
  # reports nothing, so without this the bar would sit empty through a phase
  # that can take minutes on a few million cells.
  # uwot's own stage messages go to the log. Without them the pre-epoch phases —
  # neighbour search, smooth-kNN calibration, and the spectral initialisation,
  # which is a single-threaded eigendecomposition and on millions of cells is
  # the longest of the three — are indistinguishable from a hang. There is no
  # epoch to count yet, so the only honest signal is what uwot says it is doing.
  # Markers layout, supercell fit set: SuperCellCyto groups each sample's cells
  # into supercells of about supercell_gamma cells; UMAP is fitted on their
  # means and every cell is projected onto that model.
  # What the markers layout embeds: the markers, plus the gate columns when
  # --gate_separation is above 0 (see gate_separation_columns()). Colours and
  # saved values are unaffected; only the distances UMAP sees change.
  embed_matrix <- combined_matrix
  if (!gates_layout && gate_separation > 0) {
    gcols <- gate_separation_columns(combined_subpop, combined_matrix,
                                     gate_separation)
    if (is.null(gcols)) {
      cat("Gate separation: only one gate label among the cells; nothing to ",
          "separate.\n", sep = "")
    } else {
      embed_matrix <- cbind(combined_matrix, gcols)
      cat(sprintf(paste0("Gate separation %g: %d gate label(s) (%s colouring) ",
                         "added to the embedding as columns; cells in ",
                         "different gates are pushed apart.\n"),
                  gate_separation, ncol(gcols), subpop_mode))
    }
    rm(gcols)
  }
  sc_fit <- NULL
  key_fitset <- key_extract
  if (supercell_fit) {
    key_fitset <- ckpt_key(list("supercells", key_extract, supercell_gamma,
                                as.integer(supercell_k_knn),
                                as.integer(supercell_max), seed,
                                if (gate_separation > 0) gate_separation))
    if (use_ckpt) sc_fit <- ckpt_load(output_dir, "supercells", key_fitset)
  }
  if (supercell_fit && is.null(sc_fit)) {
    emit_progress("embed", 0L, 100L, "SuperCellCyto supercells")
    sc_started <- Sys.time()
    set.seed(seed)
    # A cohort of millions at gamma 20 gives hundreds of thousands of
    # supercells: a slow fit, and past SPECTRAL_LIMIT (run_umap) a PCA start
    # instead of the spectral one. Raising gamma keeps the fit at most
    # supercell_max supercells, every sample still grouped in proportion.
    gamma_used <- max(supercell_gamma,
                      ceiling(nrow(combined_matrix) / as.integer(supercell_max)))
    if (gamma_used > supercell_gamma) {
      cat(sprintf(paste0("SuperCellCyto: gamma raised from %g to %g so that %s ",
                         "cells give at most %s supercells (--supercell_max).\n"),
                  supercell_gamma, gamma_used,
                  format(nrow(combined_matrix), big.mark = ","),
                  format(as.integer(supercell_max), big.mark = ",")))
    }
    sc_fit <- supercell_fit_matrix(embed_matrix, file_labels,
                                   gamma = gamma_used,
                                   k_knn = supercell_k_knn,
                                   n_workers = n_cores)
    cat(sprintf(paste0("SuperCellCyto: %s cells -> %s supercells (gamma %g, ",
                       "k_knn %d)%s in %s; fitting on them, projecting every ",
                       "cell.\n"),
                format(sc_fit$n_cells, big.mark = ","),
                format(sc_fit$n_supercells, big.mark = ","), gamma_used,
                as.integer(supercell_k_knn),
                if (sc_fit$passed_through) sprintf(
                  ", plus %d cell(s) of samples too small to group",
                  sc_fit$passed_through) else "",
                format(round(difftime(Sys.time(), sc_started), 1))))
    sc_fit$seconds <- as.numeric(difftime(Sys.time(), sc_started, units = "secs"))
    sc_fit$gamma_used <- gamma_used
    # Saved at once: an hour of a cohort run goes into these, and nothing after
    # this point (the embedding, the plots) may cost them.
    make_outdir(output_dir)
    write_table_gz(cbind(sc_fit$supercells, as.data.frame(sc_fit$matrix)),
                   run_path(output_dir, "data", "supercells.csv.gz"))
    cat("Wrote supercells.csv.gz (the fit set: one row per supercell)\n")
    if (use_ckpt) ckpt_save(output_dir, "supercells", key_fitset, sc_fit)
  }
  # The fitted model is saved the moment it exists, before the projection
  # (most of an hour on a cohort), so a failure there cannot cost the fit:
  # uwot::load_uwot() gives it back and uwot::umap_transform() places cells on
  # this same map without refitting. The path through run_dimred() has no such
  # hook; its model is saved after the embedding instead (below).
  # The embedding's key: the fit set's, plus everything that shapes the fit.
  key_embed <- ckpt_key(list(
    "embedding", key_fitset, layout, method, fit_subset,
    if (!is.null(model_file)) ckpt_input_stamp(model_file, "", NULL),
    if (supercell_fit) NULL else umap_subsample, n_neighbors, min_dist, metric,
    seed, umap_init, embed_sgd_cores, nn_method,
    if (gate_separation > 0) gate_separation,
    if (gates_layout) list(gate_gap, gate_tier_ratio, gate_softness, supervise)))
  model_saved <- FALSE
  save_model <- function(model) {
    make_outdir(output_dir)
    # save_uwot() will not overwrite, and a rerun into this folder with other
    # settings must not leave the old model beside new coordinates: written
    # under a temporary name, then renamed over the old one.
    model_saved <<- tryCatch({
      final <- run_path(output_dir, "model", "umap_model.uwot")
      tmp <- paste0(final, ".part")
      unlink(tmp)
      uwot::save_uwot(model, tmp)
      if (!file.rename(tmp, final)) stop("could not rename ", tmp)
      cat("Wrote umap_model.uwot (load with uwot::load_uwot())\n")
      TRUE
    }, error = function(e) {
      cat("NOTE: umap_model.uwot not written: ", conditionMessage(e), "\n",
          sep = "")
      FALSE
    })
    if (model_saved && use_ckpt) {
      ckpt_save(output_dir, "model", key_embed, "model/umap_model.uwot")
    }
  }
  # Resume: the finished coordinates, or else the fitted model (then only the
  # projection runs).
  umap_coords <- if (use_ckpt) ckpt_load(output_dir, "embedding", key_embed)
                 else NULL
  saved_model <- NULL
  if (is.null(umap_coords) && use_ckpt && !is.null(ckpt_load(output_dir, "model", key_embed)) &&
      file.exists(run_path(output_dir, "model", "umap_model.uwot"))) {
    saved_model <- tryCatch(
      uwot::load_uwot(run_path(output_dir, "model", "umap_model.uwot")),
      error = function(e) NULL)
    model_saved <- !is.null(saved_model)
  }
  if (is.null(umap_coords) && !is.null(model_file)) {
    saved_model <- uwot::load_uwot(model_file)
    # The model's run: model/ inside it (2026-10-15 on), or the run folder.
    info_file <- c(file.path(dirname(dirname(model_file)), "run_info.json"),
                   file.path(dirname(model_file), "run_info.json"))
    info_file <- info_file[file.exists(info_file)][1]
    if (!is.na(info_file) && requireNamespace("jsonlite", quietly = TRUE)) {
      mk <- unlist(jsonlite::fromJSON(info_file)$features$markers)
      if (length(mk) && !identical(mk, colnames(combined_matrix))) {
        stop("--umap_model was fitted on markers ", paste(mk, collapse = ", "),
             "; this run embeds ", paste(colnames(combined_matrix), collapse = ", "),
             ". Pass the same --markers, in the same order.", call. = FALSE)
      }
    }
    cat("Projecting every cell onto the saved model ", model_file, "\n", sep = "")
    if (!identical(normalizePath(run_path(output_dir, "model", "umap_model.uwot"),
                                 mustWork = FALSE), model_file)) {
      make_outdir(output_dir)
      file.copy(model_file, run_path(output_dir, "model", "umap_model.uwot"),
                overwrite = TRUE)
    }
    model_saved <- TRUE
  }
  embedding_resumed <- !is.null(umap_coords)
  # From here on the sub-phases report themselves: emit_uwot_stage()
  # (R_shared/stats_utils.R) translates uwot's stage announcements into
  # embed_nn / embed_init, and the epoch callback carries the optimisation.
  emit_progress("embed_nn", 0L, 3L, "starting the neighbour search")
  if (embedding_resumed) {
    # Coordinates from the checkpoint; nothing to compute.
  } else if (!is.null(model_file)) {
    umap_coords <- run_umap(
      embed_matrix, n_neighbors = n_neighbors, min_dist = min_dist,
      metric = metric, seed = seed, n_threads = embed_cores,
      sgd_threads = embed_sgd_cores, nn_method = nn_method, progress = TRUE,
      verbose = TRUE, init = umap_init, model = saved_model,
      transform_threads = max(embed_cores, n_cores))
  } else if (supercell_fit) {
    umap_coords <- run_umap(
      embed_matrix, n_neighbors = n_neighbors, min_dist = min_dist,
      metric = metric, seed = seed, n_threads = embed_cores,
      sgd_threads = embed_sgd_cores, nn_method = nn_method, progress = TRUE,
      verbose = TRUE, init = umap_init, fit_matrix = sc_fit$matrix,
      transform_threads = max(embed_cores, n_cores), on_fit = save_model,
      model = saved_model)
  } else if (gates_layout) {
    umap_coords <- run_umap(
      features, n_neighbors = n_neighbors, min_dist = min_dist,
      metric = metric, seed = seed, n_threads = embed_cores,
      sgd_threads = embed_sgd_cores, nn_method = nn_method, progress = TRUE,
      verbose = TRUE, init = umap_init,
      y = if (supervise > 0) factor(gl$combination[fit_idx]) else NULL,
      target_weight = supervise, fit_idx = fit_idx,
      # Projection's neighbour search on every plotting core: same neighbours
      # at any thread count, so the coordinates do not depend on it.
      transform_threads = max(embed_cores, n_cores), on_fit = save_model,
      model = saved_model)
  } else {
    dimred <- run_dimred(
      expr_matrix = embed_matrix,
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
  }
  stopifnot(nrow(umap_coords) == nrow(combined_matrix))
  rm(embed_matrix)
  emb_check <- embedding_check(umap_coords)
  if (isTRUE(emb_check$suspect)) {
    cat("WARNING: this embedding looks wrong. ",
        paste(emb_check$reasons, collapse = " "), "\n", sep = "")
  }
  # The random fit set (--umap_subsample) goes through run_dimred(), which
  # has no hook to save the model early; it comes back attached to the
  # coordinates and is saved now. Any model already saved is not written twice.
  umap_model <- attr(umap_coords, "model")
  attr(umap_coords, "model") <- NULL
  if (!is.null(umap_model) && !model_saved) save_model(umap_model)
  rm(umap_model, saved_model)
  if (use_ckpt && !embedding_resumed) {
    ckpt_save(output_dir, "embedding", key_embed, umap_coords)
  }

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

  # One row per plotted cell in the form plot_gate_layout() takes. In the
  # markers layout the "combination" is the subgate label, every cell weighs 1
  # and which cells a --umap_subsample fit used is not tracked.
  cells <- NULL
  if (!is.null(gl) || length(split_cols)) {
    cells <- data.frame(
      UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2],
      sample = file_labels, workspace = wsp_labels_per_cell,
      combination = if (gates_layout) gl$combination else combined_subpop,
      weight = if (gates_layout) gl$weight else 1,
      fit = if (gates_layout) seq_len(nrow(umap_df)) %in% fit_idx else TRUE,
      stringsAsFactors = FALSE)
    for (cl in annot_cols) {
      cells[[paste0("annot_", cl)]] <- umap_df[[paste0("annot_", cl)]]
    }
  }

  # --- 6b. Gates layout: the record of the run, written BEFORE any plot. ------
  # Everything needed to redraw, re-score or audit the embedding without
  # re-reading a workspace: coordinates and calls per cell, every gate
  # decision, the true counts per sample, how well combinations separate, and
  # the exact inputs. A plotting failure must not cost the run these.
  if (gates_layout) {
    make_outdir(output_dir)
    calls_df <- as.data.frame(gl$calls)
    names(calls_df) <- paste0("call_", colnames(gl$calls))
    first <- c("UMAP1", "UMAP2", "sample", "workspace", "combination")
    vals <- as.data.frame(saved_values$matrix)
    names(vals) <- paste0("value_", colnames(saved_values$matrix))
    cells_out <- cbind(cells[, first], calls_df,
                       cells[, setdiff(names(cells), first), drop = FALSE], vals)
    write_table_gz(cells_out, run_path(output_dir, "data", "cells.csv.gz"))
    utils::write.csv(gl$decisions, run_path(output_dir, "data", "gate_decisions.csv"),
                     row.names = FALSE)
    utils::write.csv(gl$counts, run_path(output_dir, "data", "combination_counts.csv"),
                     row.names = FALSE)
    sep_metrics <- separation_metrics(umap_coords, gl$combination,
                                      groups = wsp_labels_per_cell, seed = seed,
                                      samples = file_labels)
    if (!is.null(sep_metrics)) {
      utils::write.csv(sep_metrics,
                       run_path(output_dir, "data", "separation_metrics.csv"),
                       row.names = FALSE)
      macro <- sep_metrics[1, ]
      cat(sprintf(paste0("Separation (macro average over %d combination(s)): ",
                         "15-NN purity %.3f, isolated %.3f, silhouette %.3f, ",
                         "workspace mixing %.3f, sample mixing %.3f\n"),
                  nrow(sep_metrics) - 1L, macro$knn_purity, macro$isolated,
                  macro$silhouette, macro$group_mixing, macro$sample_mixing))
      # Islands whose samples sit apart: worth a look before reading a split
      # or fold-change plot of them as biology.
      low <- sep_metrics[-1, , drop = FALSE]
      low_mix <- function(v) is.finite(v) & v < 0.5
      low <- low[(low_mix(low$sample_mixing) | low_mix(low$group_mixing)) &
                   low$cells_scored >= 200, , drop = FALSE]
      if (nrow(low)) {
        cat("  Low mixing (< 0.5; samples / workspaces) in: ",
            paste(sprintf("%s (%.2f / %.2f)", low$label, low$sample_mixing,
                          low$group_mixing), collapse = "; "), "\n", sep = "")
      }
      # Islands that are apart in feature space can still be drawn touching:
      # UMAP places small disconnected pieces of the graph in 2-D with little
      # to keep them off each other. Named so a reader checks before reading
      # two touching colours as one population.
      touch <- sep_metrics[-1, , drop = FALSE]
      touch <- touch[is.finite(touch$isolated) & touch$isolated < 0.9, ,
                     drop = FALSE]
      if (nrow(touch)) {
        cat("  Islands drawn touching another (isolated < 0.9): ",
            paste(sprintf("%s (%.2f)", touch$label, touch$isolated),
                  collapse = "; "), "\n", sep = "")
      }
    }
    if (!is.null(sample_metadata) && file.exists(sample_metadata)) {
      file.copy(sample_metadata, run_path(output_dir, "data", "sample_metadata.csv"),
                overwrite = TRUE)
    }
    write_gate_run_info(
      file.path(output_dir, "run_info.json"),
      parameters = list(
        layout = layout, population = population,
        population_by_workspace = as.list(ext$pop_by_wsp),
        gate_gap = gate_gap, gate_tier_ratio = gate_tier_ratio,
        viability_pattern = if (is.null(viability_pattern))
          VIABILITY_MARKER_PATTERN else viability_pattern,
        gate_softness = gate_softness,
        supervise = supervise, cells_per_combination = cells_per_combination,
        min_combination_cells = min_combination_cells,
        filter = filter,
        split_by = I(split_cols), level_order = level_order,
        split_reference = split_reference, n_neighbors = n_neighbors,
        min_dist = min_dist, metric = metric, umap_init = umap_init,
        fit_budget = fit_budget,
        markers = if (length(markers)) I(markers) else NULL,
        marker_aliases = marker_aliases, sample_metadata = sample_metadata,
        sample_pattern = sample_pattern,
        annotate = if (length(annotate)) I(annotate) else NULL,
        workspaces = workspaces, group = group,
        normalize_umap = ext$normalize_umap,
        # The scale of the value_<marker> columns in cells.csv.gz.
        values_scale = saved_values$scale,
        gate_separation = gate_separation),
      seed = seed, embed_threads = embed_cores, sigma = sigma,
      gate_markers = gate_markers, extra_markers = extra_markers,
      tiers = gl$tiers, decisions = gl$decisions, filled = gl$filled,
      n_cells = nrow(umap_df), n_fit = length(fit_idx),
      n_samples = length(sample_labels),
      embedding = emb_check,
      samples_left_out = names(ext$sample_keep)[!ext$sample_keep],
      combinations = sort(tapply(gl$weight, gl$combination, sum),
                          decreasing = TRUE))
    written <- c("data/cells.csv.gz", "data/gate_decisions.csv",
                 "data/combination_counts.csv", "data/separation_metrics.csv",
                 "run_info.json", "model/umap_model.uwot")
    written <- written[file.exists(file.path(output_dir, written))]
    cat("Wrote ", paste(written, collapse = ", "), "\n", sep = "")
  }

  # --- 6c. Markers layout: the record of the run, written BEFORE any plot. ----
  # Coordinates of every cell with its sample, workspace, label, metadata, the
  # values it was embedded from and (with supercells) its supercell, so the
  # figures can be redrawn or re-analysed without re-reading a workspace or
  # recomputing the embedding.
  if (!gates_layout) {
    make_outdir(output_dir)
    cells_out <- data.frame(
      UMAP1 = umap_df$UMAP1, UMAP2 = umap_df$UMAP2,
      sample = umap_df$file, workspace = umap_df$workspace,
      label = combined_subpop, stringsAsFactors = FALSE)
    if (!is.null(sc_fit)) cells_out$supercell <- sc_fit$cell_supercell
    for (cl in annot_cols) {
      cells_out[[paste0("annot_", cl)]] <- umap_df[[paste0("annot_", cl)]]
    }
    vals <- as.data.frame(saved_values$matrix)
    names(vals) <- paste0("value_", colnames(saved_values$matrix))
    write_table_gz(cbind(cells_out, vals), run_path(output_dir, "data", "cells.csv.gz"))
    rm(cells_out, vals)
    if (!is.null(gl)) {
      utils::write.csv(gl$decisions, run_path(output_dir, "data", "gate_decisions.csv"),
                       row.names = FALSE)
      utils::write.csv(gl$counts, run_path(output_dir, "data", "combination_counts.csv"),
                       row.names = FALSE)
    }
    if (!is.null(sample_metadata) && file.exists(sample_metadata)) {
      file.copy(sample_metadata, run_path(output_dir, "data", "sample_metadata.csv"),
                overwrite = TRUE)
    }
    write_gate_run_info(
      file.path(output_dir, "run_info.json"),
      parameters = list(
        layout = layout, population = population,
        population_by_workspace = as.list(ext$pop_by_wsp),
        subpop_mode = subpop_mode,
        fit_subset = fit_subset,
        umap_model = model_file,
        background = background,
        pair_by = if (length(pair_by)) I(pair_by) else NULL,
        supercell_gamma = if (supercell_fit) supercell_gamma else NULL,
        supercell_k_knn = if (supercell_fit) as.integer(supercell_k_knn) else NULL,
        supercell_max = if (supercell_fit) as.integer(supercell_max) else NULL,
        umap_subsample = if (supercell_fit) NULL else umap_subsample,
        filter = filter, split_by = I(split_cols), level_order = level_order,
        split_reference = split_reference, n_neighbors = n_neighbors,
        min_dist = min_dist, metric = metric, umap_init = umap_init,
        normalize_data = ext$normalize_data, downsample = downsample,
        markers = if (length(markers)) I(markers) else NULL,
        marker_aliases = marker_aliases, sample_metadata = sample_metadata,
        sample_pattern = sample_pattern,
        annotate = if (length(annotate)) I(annotate) else NULL,
        workspaces = workspaces, group = group,
        normalize_umap = ext$normalize_umap,
        # The scale of the value_<marker> columns in cells.csv.gz.
        values_scale = saved_values$scale,
        gate_separation = gate_separation),
      seed = seed, embed_threads = embed_cores,
      gate_markers = use_markers,
      n_cells = nrow(umap_df),
      n_fit = if (!is.null(sc_fit)) nrow(sc_fit$matrix)
              else if (umap_subsample > 0) min(umap_subsample, nrow(umap_df))
              else nrow(umap_df),
      n_samples = length(sample_labels),
      embedding = emb_check,
      samples_left_out = names(ext$sample_keep)[!ext$sample_keep],
      combinations = sort(table(combined_subpop), decreasing = TRUE),
      supercells = if (is.null(sc_fit)) NULL else list(
        supercells = sc_fit$n_supercells, gamma_used = sc_fit$gamma_used,
        passed_through = sc_fit$passed_through,
        seconds = round(sc_fit$seconds)))
    written <- c("data/supercells.csv.gz", "model/umap_model.uwot", "data/cells.csv.gz",
                 "data/gate_decisions.csv", "data/combination_counts.csv", "run_info.json")
    written <- written[file.exists(file.path(output_dir, written))]
    cat("Wrote ", paste(written, collapse = ", "), "\n", sep = "")
  }

  # --- 7. Plotting. ----------------------------------------------------------
  # The existing figures first, the styled combination, split, fold-change and
  # marker plots last: plot_gate_layout() draws its marker plots under the same
  # umap_<pop>_marker_<m>_ALL names the app already files, so it has to be the
  # one that writes them last.
  facets_started <- Sys.time()
  plot_flowjo_umap_facets(
    umap_df = umap_df, combined_matrix = combined_matrix,
    use_markers = use_markers,
    # The gates layout colours by combination in its own plots; the subgate
    # plots of the markers layout would only repeat them in another palette.
    subpop_labels = if (gates_layout) character(0) else subpop_labels,
    annot_cols = annot_cols, population = population,
    output_dir = output_dir, n_workspaces = length(gs_list),
    subpop_mode_tag = subpop_mode_tag, n_cores = n_cores,
    colour_matrix = ext$colour_matrix,
    subpop_by_level = if (gates_layout) list() else ext$subpop_by_level %||% list(),
    normalize = ext$normalize_data, normalize_umap = ext$normalize_umap,
    point_size = point_size %||% 0.4,
    point_alpha = point_alpha %||% 0.6, plot_width = plot_width,
    plot_height = plot_height, output_formats = output_formats)

  gate_index <- NULL
  if (!is.null(cells)) {
    if (exists("plot_gate_layout", mode = "function")) {
      gate_index <- plot_gate_layout(
        cells, colour_matrix = ext$colour_matrix, use_markers = use_markers,
        population = population, output_dir = output_dir,
        split_cols = split_cols, level_order = level_order_list,
        split_reference = split_ref_list,
        counts = if (!is.null(gl)) gl$counts else NULL,
        background = bg_map,
        point_size = point_size, point_alpha = point_alpha,
        normalize_umap = ext$normalize_umap,
        plot_width = plot_width, plot_height = plot_height,
        output_formats = output_formats, n_cores = n_cores)
    } else {
      cat("NOTE: plot_gate_layout() is not available (src/visualization/",
          "umap_split_plots.R was not sourced), so the combination, split and ",
          "fold-change plots were not drawn.\n", sep = "")
    }
  }

  # What was drawn, for the app: the gate-layout index plus the plots the
  # function above wrote in this run (found by modification time).
  if (!is.null(cells)) {
    write_plot_index(file.path(output_dir, "plot_index.json"), gate_index,
                     facet_plot_entries(output_dir, facets_started,
                                        exclude = gate_index$file))
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
                 output_dir = normalizePath(output_dir),
                 layout = layout, cells = cells, gate_layout = gl,
                 plot_index = gate_index))
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

  # Capped so a deep gate path cannot crowd out the sample tag that follows it
  # in the per-file plot names. A full boolean gate name runs past 60 characters
  # easily ("Viable cells/CD8+/CD107a-IFNg-IL2-TNF-" and friends).
  pop_tag <- create_safe_filename(sub("^/", "", population), max_length = 60)

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

  # Every figure goes into its kind's folder of the run (run_subdir()).
  in_subdir <- function(subdir, basename_no_ext) {
    d <- file.path(output_dir, subdir)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.path(d, basename_no_ext)
  }
  save_one <- function(p, basename_no_ext, subdir = "") {
    save_plot_formats(
      p,
      in_subdir(subdir, basename_no_ext),
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
  save_categorical <- function(p, basename_no_ext, legend_strings, ncol = 1L,
                               subdir = run_subdir("overview")) {
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
      in_subdir(subdir, basename_no_ext),
      width   = width,
      height  = plot_height,
      formats = output_formats
    )
  }

  # The overview plots, each drawing every cell, are collected as jobs and
  # drawn in parallel by run_overview_jobs() below. Every job reads only
  # shared, read-only objects and writes its own file.
  overview_jobs <- list()
  add_job <- function(name, fn) {
    overview_jobs[[length(overview_jobs) + 1L]] <<- list(name = name, fn = fn)
    invisible(NULL)
  }
  run_overview_jobs <- function() {
    n_job <- length(overview_jobs)
    if (!n_job) return(invisible(NULL))
    started <- Sys.time()
    invisible(parallel_lapply(seq_len(n_job), function(i) {
      overview_jobs[[i]]$fn()
      TRUE
    }, n_cores = min(max(1L, as.integer(n_cores)), 12L)))
    cat(sprintf("Saved %d overview plot(s) in %s\n", n_job,
                format(round(difftime(Sys.time(), started), 1))))
  }

  add_job("by file", function() {
    # 7a. Categorical plot coloured by file/sample (all cells).
    file_levels <- sort(unique(umap_df$file))
    file_ncol <- legend_ncol(length(file_levels), plot_height)
    p_file <- ggplot(umap_df[draw_order_shuffled(nrow(umap_df)), , drop = FALSE],
                     aes(x = UMAP1, y = UMAP2, color = file)) +
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
  })

  add_job("by workspace", function() {
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
      df_ws <- df_ws[draw_order_shuffled(nrow(df_ws)), , drop = FALSE]

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
  })

  add_job("by subpopulation", function() {
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
      df_sp <- df_sp[draw_order_by_category(
        df_sp$subpop, bottom = c(.SUBPOP_MULTI, .SUBPOP_NONE)), , drop = FALSE]

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
        # Then the rarest gate last, so it is not buried under the common ones.
        df_lv <- df_lv[draw_order_by_category(df_lv$level), , drop = FALSE]

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
  })

  # 7a-quater. Sample annotations: for EACH column one plot coloured by it, and
  # then one plot per level with that level drawn over the whole embedding in
  # grey. One embedding, many views — the grey context is the point, because an
  # isolated subset gives you no idea where it sits relative to everything else,
  # which is exactly what a shared embedding is for.
  for (cl in annot_cols) local({
    cl <- cl
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

    add_job(paste0("annot ", cl), function() {
      df_a <- umap_df
      df_a$.annot <- factor(values, levels = levels_all)
      df_a <- df_a[draw_order_by_category(df_a$.annot,
                                          bottom = .ANNOT_UNASSIGNED), , drop = FALSE]
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
        legend_strings = c(lab, cl), ncol = a_ncol,
        subdir = run_subdir("annot", column = cl))
    })

    # One plot per level: the whole embedding in grey, this level on top.
    for (lv in levels_all) local({
      lv <- lv
      rows <- which(values == lv)
      if (length(rows)) add_job(paste0("annot ", cl, ": ", lv), function() {
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
          paste0("umap_", pop_tag, "_annot_", cl, "_level_", lv)),
          subdir = run_subdir("annot", column = cl))
      })
    })
  })

  # 7b. Conjoint marker plots (all cells), one job each. ONE colour scale per
  # marker, shared by the conjoint plot and every per-file plot below, so a
  # given colour means the same expression value in every plot of this marker.
  # The limits span the full min-max range of the marker across ALL pooled
  # cells. Per-file subsets only ever cover part of that range.
  marker_limits_all <- lapply(stats::setNames(nm = use_markers), function(m)
    color_limits(colour_matrix[, m]))
  for (marker in use_markers) local({
    marker <- marker
    marker_limits <- marker_limits_all[[marker]]
    add_job(paste0("marker ", marker), function() {
      df_all <- umap_df
      df_all$expr <- colour_matrix[, marker]
      df_all <- df_all[draw_order_by_value(df_all$expr), , drop = FALSE]
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
        create_safe_filename(paste0("umap_", pop_tag, "_marker_", marker, "_ALL")),
        subdir = run_subdir("marker")
      )
    })
  })

  # The overview plots above, all at once. Each draws every cell on one core
  # (about 1.3 min per plot at 8 million cells), and they are independent, so
  # in a row they were half an hour of a cohort run. Forked, so umap_df and the
  # matrices are shared copy-on-write; each worker still holds one plot's
  # working copy, hence at most 12 of them.
  run_overview_jobs()

  files_unique <- unique(file_labels)

  # One tag per sample, computed ONCE and reused for every marker. Building the
  # tag inside the marker loop made its length depend on the marker's, so the
  # same sample appeared under several names and samples sharing a prefix
  # overwrote each other's plots. See safe_label_tags() for the full story.
  file_tags <- safe_label_tags(files_unique, max_length = 40)

  # 7c. For each marker, one plot per file on the conjoint plot's colour scale.
  n_marker_done <- 0L
  for (marker in use_markers) {
    n_marker_done <- n_marker_done + 1L
    emit_progress("plot", n_marker_done, length(use_markers), marker)
    marker_limits <- marker_limits_all[[marker]]

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
    # The OUTER loop over markers stays serial: progress reads sensibly per
    # marker, and 311 items
    # is already far more parallelism than any core count here will absorb.
    plot_one_file <- function(f) {
      rows <- which(file_labels == f)               # index-aligned subset
      if (length(rows) == 0) return(invisible(NULL)) # nothing for this file
      df_f <- umap_df[rows, , drop = FALSE]
      df_f$expr <- colour_matrix[rows, marker]      # aligned colour vector
      df_f <- df_f[draw_order_by_value(df_f$expr), , drop = FALSE]

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
        # The budget is per component (pop_tag <= 60, marker capped below, tag
        # <= 40), so the composed name lands well inside any filesystem limit
        # and this outer cap never has to eat the sample tag.
        create_safe_filename(
          paste0("umap_", pop_tag, "_marker_", substr(marker, 1, 20),
                 "_file_", file_tags[[f]]),
          max_length = 180),
        subdir = run_subdir("marker_file", marker = marker)
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



# =============================================================================
# replot_from_saved() — redraw every plot of a finished run from what it saved
# (cells.csv.gz, run_info.json, combination_counts.csv), without reading a
# workspace or recomputing anything. A cohort run spends most of an hour on the
# embedding; changing a point size, a plot dimension or the core count used for
# drawing should not cost that again. The plots are written into `output_dir`
# (default: the run's own folder), replacing the old ones.
#
# Markers colour by the saved value_<marker> columns. Those are the values the
# run plotted (the gates layout) or embedded (the markers layout, after
# --normalize_data); --normalize_umap is re-applied to them.
# =============================================================================
replot_from_saved <- function(run_dir, output_dir = run_dir, n_cores = 1L,
                              point_size = NULL, point_alpha = NULL,
                              normalize_umap = NULL,
                              plot_width = 10, plot_height = 8,
                              output_formats = list(pdf = FALSE, png = TRUE,
                                                    svg = FALSE)) {
  # The layout's folders, or flat for a run before 2026-10-15 (run_find()).
  need <- c(run_find(run_dir, "data", "cells.csv.gz"),
            run_find(run_dir, "", "run_info.json"))
  if (anyNA(need)) {
    stop("replot needs cells.csv.gz and run_info.json in ", run_dir,
         " (written by runs since 2026-10-14).", call. = FALSE)
  }
  n_cores <- resolve_cores(n_cores)
  info <- jsonlite::fromJSON(need[2], simplifyVector = TRUE)
  par <- info$parameters
  cat("Replotting ", run_dir, " (", par$layout, " layout, ", info$cells,
      " cells) with ", n_cores, " core(s)\n", sep = "")

  t0 <- Sys.time()
  cells <- if (requireNamespace("data.table", quietly = TRUE) &&
               nzchar(Sys.which("zcat"))) {
    as.data.frame(data.table::fread(cmd = paste("zcat", shQuote(need[1]))))
  } else {
    utils::read.csv(gzfile(need[1]), check.names = FALSE,
                    stringsAsFactors = FALSE)
  }
  cat(sprintf("Read %s cells in %s\n", format(nrow(cells), big.mark = ","),
              format(round(difftime(Sys.time(), t0), 1))))
  label_col <- if ("combination" %in% names(cells)) "combination" else "label"
  value_cols <- grep("^value_", names(cells), value = TRUE)
  if (!length(value_cols)) {
    stop("cells.csv.gz has no value_<marker> columns, so the marker plots ",
         "cannot be coloured; rerun to save them.", call. = FALSE)
  }
  use_markers <- sub("^value_", "", value_cols)
  values <- as.matrix(cells[, value_cols, drop = FALSE])
  colnames(values) <- use_markers
  # The colour scale: the one asked for, "same" = the embedding's, and by
  # default the one the run used. The saved values are FlowJo values for runs
  # since 2026-11-03 (values_scale "none"); older runs saved them on the scale
  # they were embedded (markers) or coloured (gates) on.
  embed_scale <- par$normalize_data %||% "none"
  values_scale <- par$values_scale %||%
    (if (identical(par$layout, "gates")) par$normalize_umap %||% "none"
     else embed_scale)
  target <- if (is.null(normalize_umap)) par$normalize_umap %||% embed_scale
            else resolve_colour_scale(normalize_umap, embed_scale)
  coloured <- colour_values_for(values, values_scale, target)
  colour_matrix <- coloured$matrix
  normalize_umap <- coloured$scale
  cat("Marker colours: ", switch(normalize_umap, none = "FlowJo values",
                                 zscore = "z-scores", minmax = "min-max"),
      " (embedding: ", embed_scale, ")\n", sep = "")
  annot_cols <- sub("^annot_", "", grep("^annot_", names(cells), value = TRUE))

  umap_df <- data.frame(UMAP1 = cells$UMAP1, UMAP2 = cells$UMAP2,
                        file = as.character(cells$sample),
                        workspace = as.character(cells$workspace),
                        subpop = as.character(cells[[label_col]]),
                        stringsAsFactors = FALSE)
  for (cl in annot_cols) {
    umap_df[[paste0("annot_", cl)]] <- as.character(cells[[paste0("annot_", cl)]])
  }
  # Combinations get the styled plots below, so the plain subgate plot is drawn
  # only for subgate labels, as in the run itself.
  combos <- identical(par$layout, "gates") ||
    identical(par$subpop_mode, "combinations")
  subpop_labels <- if (combos) character(0) else
    sort(setdiff(unique(umap_df$subpop), c(.SUBPOP_MULTI, .SUBPOP_NONE)))
  population <- par$population

  facets_started <- Sys.time()
  plot_flowjo_umap_facets(
    umap_df = umap_df, combined_matrix = values, use_markers = use_markers,
    subpop_labels = subpop_labels, annot_cols = annot_cols,
    population = population, output_dir = output_dir,
    n_workspaces = length(unique(umap_df$workspace)),
    subpop_mode_tag = par$subpop_mode %||% "", n_cores = n_cores,
    colour_matrix = colour_matrix, subpop_by_level = list(),
    normalize = "none", normalize_umap = normalize_umap,
    point_size = point_size %||% 0.4, point_alpha = point_alpha %||% 0.6,
    plot_width = plot_width, plot_height = plot_height,
    output_formats = output_formats)

  split_cols <- unlist(par$split_by)
  split_cols <- split_cols[nzchar(split_cols)]
  gate_index <- NULL
  if ((combos || length(split_cols)) &&
      exists("plot_gate_layout", mode = "function")) {
    counts_file <- run_find(run_dir, "data", "combination_counts.csv")
    counts <- if (combos && !is.na(counts_file))
      utils::read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
      else NULL
    gl_cells <- data.frame(
      UMAP1 = cells$UMAP1, UMAP2 = cells$UMAP2,
      sample = as.character(cells$sample),
      workspace = as.character(cells$workspace),
      combination = as.character(cells[[label_col]]),
      weight = if ("weight" %in% names(cells)) cells$weight else 1,
      fit = if ("fit" %in% names(cells)) as.logical(cells$fit) else TRUE,
      stringsAsFactors = FALSE)
    for (cl in annot_cols) {
      gl_cells[[paste0("annot_", cl)]] <- umap_df[[paste0("annot_", cl)]]
    }
    bg_file <- run_find(run_dir, "data", "background_pairs.csv")
    bg_map <- NULL
    if (!is.null(par$background) && !is.na(bg_file)) {
      bp <- utils::read.csv(bg_file, stringsAsFactors = FALSE)
      bg_map <- stats::setNames(as.character(bp$background), as.character(bp$sample))
    }
    gate_index <- plot_gate_layout(
      gl_cells, colour_matrix = colour_matrix, use_markers = use_markers,
      population = population, output_dir = output_dir,
      split_cols = split_cols,
      level_order = parse_level_order(par$level_order),
      split_reference = parse_split_reference(par$split_reference),
      background = bg_map,
      counts = counts, point_size = point_size, point_alpha = point_alpha,
      normalize_umap = normalize_umap, plot_width = plot_width,
      plot_height = plot_height, output_formats = output_formats,
      n_cores = n_cores)
  }
  write_plot_index(file.path(output_dir, "plot_index.json"), gate_index,
                   facet_plot_entries(output_dir, facets_started,
                                      exclude = gate_index$file))
  cat("Done. Plots written to: ", output_dir, " in ",
      format(round(difftime(Sys.time(), t0), 1)), "\n", sep = "")
  invisible(output_dir)
}


# -----------------------------------------------------------------------------
# Checkpoints of umap_from_flowjo(): <output_dir>/checkpoints/<stage>.rds, each
# list(key, saved, value). The key is the md5 of the settings the stage depends
# on (base R only: tools::md5sum of the serialised settings), so a checkpoint
# made with other settings is recognised and recomputed, never reused.
# -----------------------------------------------------------------------------
ckpt_key <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(x, f, compress = FALSE, version = 2)
  unname(tools::md5sum(f))
}

# What the inputs looked like: the workspace files, the FCS files and the
# metadata table, by size and modification time. A re-exported workspace or an
# edited metadata table changes the key; reading the files would take as long
# as the extraction the checkpoint saves.
ckpt_input_stamp <- function(wsp_path, fcs_dir, sample_metadata) {
  stamp <- function(paths) {
    paths <- sort(unique(paths[nzchar(paths) & file.exists(paths)]))
    if (!length(paths)) return(NULL)
    fi <- file.info(paths, extra_cols = FALSE)
    data.frame(path = paths, size = fi$size,
               mtime = as.numeric(fi$mtime), stringsAsFactors = FALSE)
  }
  wsp <- if (length(wsp_path) && dir.exists(wsp_path[1]))
    list.files(wsp_path, pattern = "\\.wsp$", recursive = TRUE, full.names = TRUE,
               ignore.case = TRUE) else wsp_path
  fcs <- if (dir.exists(fcs_dir))
    list.files(fcs_dir, pattern = "\\.fcs$", recursive = TRUE, full.names = TRUE,
               ignore.case = TRUE) else character(0)
  list(wsp = stamp(wsp), fcs = stamp(fcs),
       metadata = stamp(if (is.null(sample_metadata)) character(0) else sample_metadata))
}

ckpt_load <- function(output_dir, stage, key) {
  f <- file.path(output_dir, "checkpoints", paste0(stage, ".rds"))
  if (!file.exists(f)) return(NULL)
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj) || !identical(obj$key, key)) {
    cat("Checkpoint '", stage, "' was saved with other settings; recomputing ",
        "it.\n", sep = "")
    return(NULL)
  }
  cat("Resuming: '", stage, "' reused from checkpoints/", stage, ".rds (saved ",
      obj$saved, ")\n", sep = "")
  obj$value
}

# Written to a temporary name and renamed, so a run killed mid-write leaves
# the previous checkpoint (or none), never a truncated one.
ckpt_save <- function(output_dir, stage, key, value) {
  d <- file.path(output_dir, "checkpoints")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(d, paste0(stage, ".rds"))
  tmp <- paste0(f, ".part")
  ok <- tryCatch({
    saveRDS(list(key = key, saved = format(Sys.time(), "%Y-%m-%d %H:%M"),
                 value = value), tmp, compress = FALSE)
    file.rename(tmp, f)
  }, error = function(e) {
    cat("NOTE: checkpoint '", stage, "' not saved: ", conditionMessage(e), "\n",
        sep = "")
    FALSE
  })
  if (isTRUE(ok)) cat("Checkpoint saved: checkpoints/", stage, ".rds\n", sep = "")
  invisible(ok)
}
