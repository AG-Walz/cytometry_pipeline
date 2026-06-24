#!/usr/bin/env Rscript
# =============================================================================
# visualization/umap_from_flowjo.R — Shared UMAP from a FlowJo workspace population
# =============================================================================
# Core logic for the FlowJo-population UMAP analysis. Opens a FlowJo .wsp, builds
# a GatingSet over a directory of FCS files, extracts the cells of ONE selected
# population from each sample (data is already compensated + transformed in the
# FlowJo working space, so NO further transform is applied), concatenates the
# per-sample expression matrices, OPTIONALLY normalizes each marker using
# statistics computed from those concatenated cells, computes a SINGLE shared
# UMAP embedding, and plots that embedding coloured by sample, by the FlowJo
# subgates found under the selected population (one colour per subgate; see
# `subpop_mode`), and by every marker (conjoint + per-file).
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
#   combined_matrix   : rows = cells (sample 1 stacked on top of sample 2)
#   file_labels       : length == nrow(combined_matrix); file_labels[i] is the
#                       source file of combined_matrix row i
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
# -----------------------------------------------------------------------------
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
       "list_populations = TRUE to see the available paths.")
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
                                      mode = "cd4cd8", select = NULL) {
  prefix    <- paste0(population, "/")
  desc      <- pop_paths_union[startsWith(pop_paths_union, prefix)]
  rel       <- substring(desc, nchar(prefix) + 1L)
  depth     <- lengths(strsplit(rel, "/", fixed = TRUE))
  base_name <- basename(rel)  # leaf gate name (last path component)

  if (!is.null(select)) {
    keep <- rel %in% select | base_name %in% select
    matched <- unique(c(rel[keep], base_name[keep]))
    missing_sel <- setdiff(select, matched)
    if (length(missing_sel)) {
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
# -----------------------------------------------------------------------------
get_subpop_membership <- function(gh, viable_idx, subpop_paths) {
  nv <- sum(viable_idx)
  M  <- matrix(FALSE, nrow = nv, ncol = length(subpop_paths))
  if (length(subpop_paths) == 0L || nv == 0L) return(M)
  gh_nodes <- flowWorkspace::gh_get_pop_paths(gh, path = "full")
  present  <- subpop_paths %in% gh_nodes
  if (any(present)) {
    mat <- flowWorkspace::gh_pop_get_indices_mat(gh, y = subpop_paths[present])
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = sum(present))
    mat <- mat[viable_idx, , drop = FALSE]     # root order -> population order
    M[, present] <- mat
  }
  M
}

# -----------------------------------------------------------------------------
# Extract one sample: returns list(mtx = cells x markers, subpop = per-cell
# subgate label). Both are row-aligned to gh_pop_get_data(gh, population); the
# subpop vector is built from the SAME population indices so it stays aligned
# through downsampling and rbind.
#
# We use explore_data(ff)$markers (per-column marker names, "" when none) for
# consistency with the rest of the repo. A channel is kept as a biological
# marker only when it has a non-empty marker name AND its CHANNEL name does not
# match the scatter/Time exclusion regex (applied to the CHANNEL name only, so a
# legitimate antibody whose marker string contains "FSC"/"Time" is not dropped).
# Kept columns are renamed to their MARKER name.
# -----------------------------------------------------------------------------
extract_population_sample <- function(gh, population, subpop_paths, subpop_depths,
                                      subpop_labels, assignment = "depth",
                                      exclude_pattern = "Time|FSC|SSC|scatter") {
  # NOTE: data is already compensated + transformed in the FlowJo working space.
  # We take exprs() AS-IS, no further transform.
  ff   <- flowWorkspace::gh_pop_get_data(gh, y = population)
  expr <- flowCore::exprs(ff)
  info <- explore_data(ff)          # dyes (channels), markers, mcols

  channels <- colnames(expr)
  markers  <- info$markers          # length == ncol(expr)

  keep <- vapply(seq_along(channels), function(i) {
    m <- markers[i]
    if (is.na(m) || m == "") return(FALSE)                    # empty marker
    if (grepl(exclude_pattern, channels[i],
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
# Main entry point: compute and plot a shared UMAP over one FlowJo population.
#
# @param wsp_path Path to the FlowJo .wsp workspace file.
# @param fcs_dir Directory containing the FCS files referenced by the wsp.
# @param population Population path to extract (leading '/' optional; a suffix
#   such as 'Single Cells/Viable cells' is matched against the gating tree).
# @param output_dir Output directory for plots.
# @param group FlowJo workspace group index to load (flowjo_to_gatingset name=).
# @param expected_samples If > 0, warn when the loaded group lacks this many.
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
                             group = 1L, expected_samples = 0L,
                             markers = NULL, subpop_mode = "furthest",
                             subpop_select = NULL, normalize = "none",
                             downsample = 0L, n_neighbors = 15L, min_dist = 0.1,
                             metric = "euclidean", seed = 42L,
                             point_size = 0.4, point_alpha = 0.6,
                             plot_width = 10, plot_height = 8,
                             output_formats = list(pdf = FALSE, png = FALSE, svg = TRUE),
                             list_populations = FALSE) {
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

  if (!file.exists(wsp_path)) stop("wsp file not found: ", wsp_path)
  if (!dir.exists(fcs_dir))   stop("fcs_dir not found: ", fcs_dir)
  wsp_path <- normalizePath(wsp_path, mustWork = TRUE)
  fcs_dir  <- normalizePath(fcs_dir,  mustWork = TRUE)

  # --- 1. Open workspace and build the GatingSet over the FCS dir. -----------
  cat("Opening FlowJo workspace:", wsp_path, "\n")
  ws <- CytoML::open_flowjo_xml(wsp_path)
  cat("Building GatingSet from FCS dir:", fcs_dir,
      "(group index =", group, ")\n")
  gs <- CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir)

  sample_names <- flowWorkspace::sampleNames(gs)
  cat("Samples in loaded group (", length(sample_names), "):\n", sep = "")
  for (s in sample_names) cat("  -", s, "\n")
  if (expected_samples > 0 && length(sample_names) != expected_samples) {
    warning(sprintf(
      "Loaded group %d has %d sample(s) but expected_samples=%d. The wrong ",
      group, length(sample_names), expected_samples),
      "workspace group may be selected; re-check group or list_populations.")
  }

  all_pops <- flowWorkspace::gs_get_pop_paths(gs, path = "full")

  # --- 2. list_populations mode: print the tree and return. ------------------
  if (isTRUE(list_populations)) {
    cat("\nPopulation tree (", length(all_pops), " nodes):\n", sep = "")
    for (p in all_pops) cat("  ", p, "\n", sep = "")
    cat("\nSelect one with population = (leading '/' optional), e.g.:\n")
    cat("  population = \"", sub("^/", "", all_pops[length(all_pops)]),
        "\"\n", sep = "")
    return(invisible(all_pops))
  }

  if (is.null(population) || !nzchar(population)) {
    stop("population is required (or use list_populations = TRUE to see options).")
  }

  population <- resolve_population(population, all_pops)
  cat("\nTarget population resolved to:", population, "\n")

  # --- 2b. Discover the FlowJo subgates UNDER `population` to colour by. ------
  pop_paths_union <- unique(unlist(lapply(sample_names, function(s)
    flowWorkspace::gh_get_pop_paths(gs[[s]], path = "full"))))
  cand <- flowjo_subgate_candidates(pop_paths_union, all_pops, population,
                                    mode = subpop_mode, select = subpop_select)
  subpop_paths  <- cand$paths
  subpop_labels <- cand$labels
  subpop_depths <- cand$depths
  subpop_mode_tag <- cand$tag

  # 'furthest' labels each cell by the most specific (smallest) gate it reaches;
  # the other modes assign by path depth (ties -> "(multiple)").
  subpop_assignment <- if (subpop_mode == "furthest" && is.null(subpop_select)) {
    "specificity"
  } else {
    "depth"
  }

  do_subpop <- length(subpop_paths) > 0
  subpop_requested <- !is.null(subpop_select) || subpop_mode != "none"
  if (subpop_requested) {
    if (do_subpop) {
      cat(sprintf("\nSubpopulation colouring (%s): %d candidate subgate(s) under %s\n",
                  subpop_mode_tag, length(subpop_paths), population))
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
  for (s in sample_names) {
    res <- extract_population_sample(gs[[s]], population, subpop_paths,
                                     subpop_depths, subpop_labels,
                                     assignment = subpop_assignment)
    mtx <- res$mtx
    sp  <- res$subpop

    # Optional per-sample downsampling AFTER gating. The SAME row index subsets
    # the subpopulation labels so they stay aligned to mtx.
    if (downsample > 0 && nrow(mtx) > downsample) {
      idx <- sort(sample(nrow(mtx), downsample))
      mtx <- mtx[idx, , drop = FALSE]
      sp  <- sp[idx]
    }

    cat(sprintf("  %s: %d cells x %d markers\n", s, nrow(mtx), ncol(mtx)))
    if (nrow(mtx) == 0) {
      warning("Sample '", s, "' has 0 cells in population '", population,
              "' after gating/downsampling; it will not contribute rows.")
    }
    per_sample_mtx[[s]]    <- mtx
    per_sample_subpop[[s]] <- sp
  }

  # Drop samples that contributed no cells (keeps alignment + plotting sane).
  nonempty <- vapply(per_sample_mtx, nrow, integer(1)) > 0
  if (!any(nonempty)) {
    stop("No sample has any cells in population '", population, "'.")
  }
  per_sample_mtx    <- per_sample_mtx[nonempty]
  per_sample_subpop <- per_sample_subpop[nonempty]

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
  raw_names   <- names(per_sample_mtx)
  clean_names <- clean_sample_label(raw_names)
  if (anyDuplicated(clean_names)) {
    warning("Cleaned sample labels collide; keeping raw sample names for display.")
    clean_names <- raw_names
  }
  names(clean_names) <- raw_names  # raw -> display lookup

  combined_matrix <- do.call(rbind, per_sample_mtx)
  file_labels     <- rep(clean_names[raw_names],
                         vapply(per_sample_mtx, nrow, integer(1)))
  stopifnot(length(file_labels) == nrow(combined_matrix))

  # Subpopulation labels, concatenated in the SAME sample order as
  # combined_matrix (the two lists are filtered/ordered identically).
  combined_subpop <- unlist(per_sample_subpop[raw_names], use.names = FALSE)
  stopifnot(length(combined_subpop) == nrow(combined_matrix))

  cat(sprintf("\nConcatenated matrix: %d cells x %d markers (%d samples)\n",
              nrow(combined_matrix), ncol(combined_matrix),
              length(per_sample_mtx)))

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
    UMAP1  = umap_coords[, 1],
    UMAP2  = umap_coords[, 2],
    file   = file_labels,
    subpop = combined_subpop,
    stringsAsFactors = FALSE
  )

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
  invisible(list(umap_coords = umap_coords, umap_df = umap_df,
                 output_dir = normalizePath(output_dir)))
}
