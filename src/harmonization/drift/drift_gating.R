# =============================================================================
# src/harmonization/drift/drift_gating.R — FlowJo workspaces to per-tube data
# =============================================================================
# Needs CytoML, flowWorkspace, flowCore (highmem only). Every function takes a
# GatingHierarchy or GatingSet and returns plain data frames so that the base-R
# modules never see Bioconductor objects.
# =============================================================================

#' Build GatingSets for one workspace. Tries the whole sample group first; if
#' CytoML refuses (mixed channel sets, e.g. a patient acquired on both
#' cytometers), falls back to one GatingSet per sample. Returns a list of
#' GatingSets, names are workspace label plus part index.
load_workspace_gatingsets <- function(wsp_path, fcs_dir, group = NULL) {
  ws <- CytoML::open_flowjo_xml(wsp_path)
  groups <- CytoML::fj_ws_get_sample_groups(ws)
  if (is.null(group)) {
    g <- setdiff(unique(groups$groupName), "Compensation")
    group <- if ("All Samples" %in% g) "All Samples" else g[1]
  }
  label <- sub("\\.wsp$", "", basename(wsp_path), ignore.case = TRUE)
  out <- list()
  whole <- tryCatch(
    CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir),
    error = function(e) e)
  if (!inherits(whole, "error")) { out[[label]] <- whole; return(out) }
  cat(sprintf("  %s: whole-group build failed (%s); building per sample\n",
              label, conditionMessage(whole)))
  n <- sum(groups$groupName == group)
  for (i in seq_len(n)) {
    part <- tryCatch(
      CytoML::flowjo_to_gatingset(ws, name = group, subset = i, path = fcs_dir),
      error = function(e) { cat(sprintf("    sample %d skipped: %s\n", i,
                                        conditionMessage(e))); NULL })
    if (!is.null(part)) out[[sprintf("%s_part%02d", label, i)]] <- part
  }
  out
}

#' The FCS file name a GatingHierarchy was built from.
gh_filename <- function(gh) {
  kw <- flowWorkspace::keyword(gh)
  f <- kw[["$FIL"]]
  if (is.null(f) || !nzchar(f)) f <- kw[["FILENAME"]]
  if (is.null(f) || !nzchar(f)) f <- sub("_\\d+$", "", flowWorkspace::sampleNames(gh))
  basename(f)
}

#' Full path of the first node whose path ends with `suffix` (case-insensitive).
find_population <- function(gh, suffix) {
  paths <- flowWorkspace::gh_get_pop_paths(gh, path = "full")
  hit <- paths[grepl(paste0("(^|/)", gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", suffix), "$"),
                     paths, ignore.case = TRUE)]
  if (!length(hit)) return(NA_character_)
  hit[which.min(nchar(hit))]
}

#' Events of one population on the COMPENSATED LINEAR scale, as a matrix whose
#' columns are canonical marker names. Channels not in `channel_map` (a named
#' character: canonical -> channel) are dropped.
population_matrix <- function(gh, node, channel_map, n_max = 20000L, seed = 1L) {
  cf <- flowWorkspace::gh_pop_get_data(gh, y = node, inverse.transform = TRUE)
  m <- flowCore::exprs(cf)
  if (nrow(m) > n_max) { set.seed(seed); m <- m[sample.int(nrow(m), n_max), , drop = FALSE] }
  keep <- channel_map[channel_map %in% colnames(m)]
  out <- m[, keep, drop = FALSE]
  colnames(out) <- names(keep)
  out
}

#' Per-dimension bounds of a gate, mapped back to the compensated linear scale
#' via the hierarchy's inverse transformations, so gates from different
#' workspaces (different FlowJo transforms) are comparable.
gate_bounds <- function(gh, node) {
  g <- flowWorkspace::gh_pop_get_gate(gh, node)
  inv <- flowWorkspace::gh_get_transformations(gh, inverse = TRUE)
  dims <- names(flowCore::parameters(g))
  if (methods::is(g, "rectangleGate")) {
    lo <- g@min; hi <- g@max
  } else if (methods::is(g, "polygonGate")) {
    lo <- apply(g@boundaries, 2, min); hi <- apply(g@boundaries, 2, max)
  } else if (methods::is(g, "ellipsoidGate")) {
    r <- sqrt(diag(g@cov)) * sqrt(g@distance); lo <- g@mean - r; hi <- g@mean + r
  } else return(NULL)
  back <- function(v, d) if (!is.null(inv[[d]])) inv[[d]](v) else v
  data.frame(node = node, dimension = dims,
             lower = vapply(seq_along(dims), function(i) back(unname(lo[i]), dims[i]), numeric(1)),
             upper = vapply(seq_along(dims), function(i) back(unname(hi[i]), dims[i]), numeric(1)),
             gate_type = class(g)[1], stringsAsFactors = FALSE)
}

#' Counts of every node of every sample in a GatingSet, one row per node.
gatingset_counts <- function(gs) {
  ct <- flowWorkspace::gs_pop_get_count_fast(gs, statistic = "count", path = "full")
  ct <- as.data.frame(ct, stringsAsFactors = FALSE)
  fn <- vapply(flowWorkspace::sampleNames(gs), function(s) gh_filename(gs[[s]]), character(1))
  ct$filename <- unname(fn[ct$name])
  ct[, c("filename", "Population", "Parent", "Count", "ParentCount")]
}
