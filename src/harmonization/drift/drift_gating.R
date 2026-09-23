# src/harmonization/drift/drift_gating.R — FlowJo workspaces to per-tube data
# Needs CytoML, flowWorkspace, flowCore (highmem only). Every function takes a
# GatingHierarchy or GatingSet and returns plain data frames so that the base-R
# modules never see Bioconductor objects.

#' Escape a literal string for use inside a regular expression.
regex_escape <- function(x) gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", x)

#' Build GatingSets for one workspace. Tries the whole sample group first; if
#' CytoML refuses because the samples do not share one channel set (a patient
#' acquired on both cytometers, or with and without FSC-W), the group is split
#' into channel-homogeneous parts and each part is built with `subset` = the
#' sample names of that part, which is the form this CytoML accepts.
#' `signatures` is the step-2 table (filename, panel_signature); files not in it
#' are read with flowCore::read.FCSheader. Returns a named list of GatingSets.
load_workspace_gatingsets <- function(wsp_path, fcs_dir, signatures = NULL) {
  ws <- CytoML::open_flowjo_xml(wsp_path)
  groups <- CytoML::fj_ws_get_sample_groups(ws)
  g <- setdiff(unique(groups$groupName), "Compensation")
  group <- if ("All Samples" %in% g) "All Samples" else g[1]
  label <- sub("\\.wsp$", "", basename(wsp_path), ignore.case = TRUE)
  out <- list()
  whole <- tryCatch(CytoML::flowjo_to_gatingset(ws, name = group, path = fcs_dir), error = function(e) e)
  if (!inherits(whole, "error")) { out[[label]] <- whole; return(out) }
  if (!grepl("channel inconsistency", conditionMessage(whole), ignore.case = TRUE)) stop(whole)

  # sample names of the group: join the group table (groupName, sampleID) to the
  # full sample table (sampleID, name); group_id numbering differs between CytoML
  # versions, the sampleID join does not
  all_samples <- CytoML::fj_ws_get_samples(ws)
  ids <- groups$sampleID[groups$groupName == group]
  names_ws <- as.character(all_samples$name[all_samples$sampleID %in% ids])
  if (!length(names_ws)) names_ws <- as.character(all_samples$name)
  if (!length(names_ws)) stop(label, ": workspace lists no samples")
  sig <- vapply(names_ws, function(nm) {
    fn <- basename(nm)
    hit <- if (!is.null(signatures)) signatures$panel_signature[signatures$filename == fn] else character(0)
    if (length(hit)) return(paste(sort(strsplit(hit[1], ";", fixed = TRUE)[[1]]), collapse = "|"))
    f <- list.files(fcs_dir, pattern = paste0("^", regex_escape(fn), "$"), recursive = TRUE, full.names = TRUE)
    if (!length(f)) return(NA_character_)
    kw <- tryCatch(flowCore::read.FCSheader(f[1])[[1]], error = function(e) NULL)
    if (is.null(kw)) return(NA_character_)
    n <- as.integer(kw[["$PAR"]])
    paste(sort(vapply(seq_len(n), function(k) kw[[paste0("$P", k, "N")]], character(1))), collapse = "|")
  }, character(1))
  known <- !is.na(sig)
  parts <- if (any(known)) split(names_ws[known], sig[known]) else list(names_ws)
  parts <- parts[order(-lengths(parts))]
  if (any(!known)) parts[[1]] <- c(parts[[1]], names_ws[!known])
  cat(sprintf("  %s: %d channel sets, building %d part(s)\n", label, length(parts), length(parts)))
  for (i in seq_along(parts)) {
    # do.call splices the VALUE in: flowjo_to_gatingset() runs substitute(subset)
    # and would otherwise see the symbol `parts[[i]]`, which it cannot evaluate
    part <- tryCatch(
      do.call(CytoML::flowjo_to_gatingset,
              list(ws, name = group, subset = parts[[i]], path = fcs_dir)),
      error = function(e) { cat(sprintf("    part %d (%d samples) failed: %s\n", i, length(parts[[i]]),
                                        conditionMessage(e))); NULL })
    if (!is.null(part)) out[[sprintf("%s_part%d", label, i)]] <- part
  }
  out
}

#' GatingSets of one workspace for the step 4/5 loops: announces the workspace
#' and turns a load failure into a message and an empty list.
open_workspace <- function(wsp_path, fcs_dir, signatures) {
  cat("\nWorkspace", basename(wsp_path), "\n")
  tryCatch(load_workspace_gatingsets(wsp_path, fcs_dir, signatures = signatures),
           error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); list() })
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
  hit <- paths[grepl(paste0("(^|/)", regex_escape(suffix), "$"), paths, ignore.case = TRUE)]
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
  # FlowJo-derived cytoframes name compensated channels "Comp-<channel>"
  bare <- sub("^Comp-", "", colnames(m))
  idx <- match(channel_map, bare)
  keep <- !is.na(idx)
  out <- m[, idx[keep], drop = FALSE]
  colnames(out) <- names(channel_map)[keep]
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
  back <- function(v, d) if (!is.finite(v)) NA_real_ else if (!is.null(inv[[d]])) inv[[d]](v) else v
  data.frame(node = node, dimension = sub("^Comp-", "", dims),
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
