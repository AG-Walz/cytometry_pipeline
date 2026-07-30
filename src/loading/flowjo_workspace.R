#!/usr/bin/env Rscript
# =============================================================================
# loading/flowjo_workspace.R — THE FlowJo workspace layer
# =============================================================================
# Every call into CytoML/flowWorkspace for a .wsp workspace lives here and
# nowhere else: opening a workspace, listing its sample groups, listing a
# group's samples, and building a GatingSet. `load_flowjo_wsp.R` and
# `clustering/flowjo_import.R` delegate to these rather than calling CytoML
# themselves, so an API change is a fix in ONE file.
#
# That is not hypothetical: `fj_ws_get_sample_groups` is exported by CytoML on
# some installs and flowWorkspace on others, and the hard-coded
# `flowWorkspace::` form that used to appear in two files failed outright on a
# real analysis host. See find_wsp_accessor().
#
# On top of that it turns two FOLDER ROOTS — one holding .wsp workspaces, one
# holding .fcs files, both searched recursively — into a list of ready-to-load
# workspace specs.
#
# Why this exists: a .wsp records a full `uri` per sample, but it is the path on
# the (usually Windows) machine that built the workspace, so it does not exist
# on the analysis host. What the workspace reliably gives us is the LIST OF FCS
# FILE NAMES. So we read those names out of the XML first, index the FCS root
# recursively by file name, and pair the two. That imposes no directory
# convention on the data: a flat dump, mirrored subfolders and a messy mix all
# work, and a workspace's files may even be split across directories.
#
# Everything is resolved BEFORE any GatingSet is built, and every problem across
# every workspace is reported together (see format_resolution_report), because a
# GatingSet load is slow and discovering a missing file per 20-minute load is
# not a workflow.
#
# This file holds only function definitions (no top-level side effects). The
# caller is responsible for having CytoML/flowWorkspace available; the pure
# discovery, indexing and labelling helpers below deliberately need neither, so
# they can be tested with base R alone (test_data/test_umap_multi_wsp.R).
# =============================================================================

# -----------------------------------------------------------------------------
# The four CytoML entry points. Everything else in the repo goes through these.
# -----------------------------------------------------------------------------

#' Open a FlowJo .wsp workspace.
wsp_open <- function(wsp_path) {
  if (!file.exists(wsp_path)) stop("wsp file not found: ", wsp_path)
  CytoML::open_flowjo_xml(wsp_path)
}

#' Build a GatingSet from an open workspace.
#'
#' @param ws An object from wsp_open().
#' @param group Sample group: an integer index, or a group name.
#' @param path Either a directory holding the FCS files, or a data.frame with
#'   `sampleID` and `file` columns mapping each sample to its own file.
wsp_to_gatingset <- function(ws, group = 1L, path, ...) {
  CytoML::flowjo_to_gatingset(ws, name = group, path = path, ...)
}

# -----------------------------------------------------------------------------
# Locate the .wsp workspaces under `root`.
#
# `root` may be a single .wsp FILE (the original single-workspace usage, kept
# working unchanged) or a DIRECTORY, which is searched recursively.
# -----------------------------------------------------------------------------
find_workspace_files <- function(root) {
  if (!nzchar(root)) stop("No workspace path given.")
  if (!file.exists(root)) stop("Workspace path not found: ", root)

  info <- file.info(root)
  if (!isTRUE(info$isdir)) {
    if (!grepl("\\.wsp$", root, ignore.case = TRUE)) {
      stop("Not a .wsp file: ", root)
    }
    return(normalizePath(root, mustWork = TRUE))
  }

  hits <- list.files(root, pattern = "\\.wsp$", recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) {
    stop("No .wsp workspace found under: ", normalizePath(root),
         "\n(the directory is searched recursively)")
  }
  sort(normalizePath(hits, mustWork = FALSE))
}

# -----------------------------------------------------------------------------
# Short display label per workspace: the file name without its extension.
#
# Two workspaces in different subfolders can share a base name, and the label
# ends up in plot legends and — via create_safe_filename() — in file names, so a
# collision would silently merge two workspaces. Colliding labels are therefore
# qualified with their parent directory, and if that still collides, with the
# full path. Non-colliding labels are never touched.
# -----------------------------------------------------------------------------
workspace_labels <- function(wsp_paths) {
  base <- sub("\\.wsp$", "", basename(wsp_paths), ignore.case = TRUE)
  out  <- base
  dup  <- base %in% base[duplicated(base)]
  if (any(dup)) {
    out[dup] <- file.path(basename(dirname(wsp_paths[dup])), base[dup])
  }
  # Still colliding (same parent name in two different trees) -> full path.
  dup2 <- out %in% out[duplicated(out)]
  if (any(dup2)) out[dup2] <- wsp_paths[dup2]
  out
}

# -----------------------------------------------------------------------------
# Recursive index of the FCS root: lower-cased file name -> every full path
# carrying that name. Built ONCE and shared by every workspace, because the
# scan is the expensive part and the workspaces all search the same tree.
# -----------------------------------------------------------------------------
index_fcs_dir <- function(root) {
  if (!nzchar(root)) stop("No FCS directory given.")
  if (!dir.exists(root)) stop("FCS directory not found: ", root)

  files <- list.files(root, pattern = "\\.fcs$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop("No .fcs file found under: ", normalizePath(root),
         "\n(the directory is searched recursively)")
  }
  files <- normalizePath(files, mustWork = FALSE)
  split(files, tolower(basename(files)))
}

# -----------------------------------------------------------------------------
# Resolve the FCS file names a workspace asks for against the index.
#
# Matching is deliberately tolerant, in this order, because FlowJo sample names
# are not always the bare file name:
#   1. the name as given (case-insensitively — the index key is lower-cased)
#   2. the name with ".fcs" appended (FlowJo sometimes drops the extension)
#   3. the name's base name, in case the workspace stored a full path. Windows
#      separators are handled explicitly: base R's basename() only splits on
#      "/", so a stored "C:\Data\Pool_001.fcs" would otherwise never match and
#      be reported as missing.
#
# Returns list(file = character vector aligned to `names`, NA where unresolved,
# missing = the unmatched names, ambiguous = named list of candidate vectors).
# Pure: no CytoML, no file system access. This is the core the tests exercise.
# -----------------------------------------------------------------------------
path_basename <- function(x) {
  sub(".*[\\\\/]", "", x)
}

resolve_fcs_names <- function(names, fcs_index) {
  file      <- rep(NA_character_, length(names))
  missing   <- character(0)
  ambiguous <- list()

  for (i in seq_along(names)) {
    nm   <- names[i]
    base <- path_basename(nm)
    keys <- unique(tolower(c(nm, paste0(nm, ".fcs"), base,
                             paste0(base, ".fcs"))))
    hit <- NULL
    for (k in keys) {
      if (!is.null(fcs_index[[k]])) { hit <- fcs_index[[k]]; break }
    }
    if (is.null(hit)) {
      missing <- c(missing, nm)
    } else if (length(hit) > 1L) {
      ambiguous[[nm]] <- hit
    } else {
      file[i] <- hit
    }
  }
  list(file = file, missing = missing, ambiguous = ambiguous)
}

# -----------------------------------------------------------------------------
# Find a workspace accessor wherever the installed stack happens to export it.
#
# These functions have moved between CytoML and flowWorkspace and been renamed
# across versions (fj_ws_get_samples was getSamples; fj_ws_get_sample_groups was
# getSampleGroups), and a hard-coded `flowWorkspace::fj_ws_get_sample_groups`
# fails outright on installs where CytoML owns it. Both packages are searched,
# CytoML first, then the legacy names, so one accessor moving is not fatal.
#
# Returns NULL rather than stopping when nothing matches, so a caller that can
# carry on without the accessor (see resolve_group) is free to.
# -----------------------------------------------------------------------------
find_wsp_accessor <- function(candidates) {
  for (pkg in c("CytoML", "flowWorkspace")) {
    if (!requireNamespace(pkg, quietly = TRUE)) next
    for (nm in candidates) {
      f <- tryCatch(getExportedValue(pkg, nm), error = function(e) NULL)
      if (is.function(f)) return(f)
    }
  }
  NULL
}

# The column holding a file/group name differs across versions, so take the
# first plausible one rather than hard-coding it.
.first_col <- function(df, candidates, what) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop("could not find the ", what, " column (saw: ",
         paste(colnames(df), collapse = ", "), "). Update ",
         "src/loading/flowjo_workspace.R.")
  }
  hit[1]
}

# -----------------------------------------------------------------------------
# The workspace's sample groups, as a table we control the shape of.
# `index` is the 1-based position used by fj_ws_get_samples(group_id=) and by
# flowjo_to_gatingset(name=<integer>).
# -----------------------------------------------------------------------------
wsp_group_table <- function(ws) {
  f <- find_wsp_accessor(c("fj_ws_get_sample_groups", "getSampleGroups"))
  if (is.null(f)) {
    stop("neither CytoML nor flowWorkspace exports fj_ws_get_sample_groups() ",
         "or getSampleGroups().")
  }
  groups <- f(ws)
  gcol <- .first_col(groups, c("groupName", "groupname", "name"), "group-name")
  nm <- unique(as.character(groups[[gcol]]))
  data.frame(
    index     = seq_along(nm),
    groupName = nm,
    n_samples = as.integer(vapply(
      nm, function(g) sum(as.character(groups[[gcol]]) == g), integer(1))),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# The FCS file names of one group, in workspace order. This one is required:
# the whole pairing scheme is built on the names a workspace asks for.
# -----------------------------------------------------------------------------
wsp_sample_files <- function(ws, group_index) {
  f <- find_wsp_accessor(c("fj_ws_get_samples", "getSamples"))
  if (is.null(f)) {
    stop("neither CytoML nor flowWorkspace exports fj_ws_get_samples() or ",
         "getSamples(); the workspace's sample list cannot be read.")
  }
  # group_id= is the modern argument name; the legacy accessor took it
  # positionally.
  samples <- tryCatch(f(ws, group_id = group_index),
                      error = function(e) f(ws, group_index))
  ncol_ <- .first_col(samples, c("name", "sampleName", "file", "filename"),
                      "file-name")
  icol  <- .first_col(samples, c("sampleID", "sampleId", "id"), "sample-ID")
  data.frame(
    sampleID = as.character(samples[[icol]]),
    name     = as.character(samples[[ncol_]]),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Resolve --group against one workspace. `group` is a character string holding
# EITHER an integer index (the original meaning, unchanged) OR a group name.
#
# A name is what makes one --group work across a folder of workspaces, where the
# group of interest can sit at a different index in each. Returns
# list(index, name, n_samples).
# -----------------------------------------------------------------------------
resolve_group <- function(group, group_table, where = "") {
  group <- trimws(as.character(group))
  ctx   <- if (nzchar(where)) paste0(" in workspace '", where, "'") else ""

  # group_table = NULL means the group accessor was unavailable. An INDEX still
  # works — flowjo_to_gatingset(name=<integer>) does its own lookup — so the
  # ordinary run is not held hostage to an accessor we only need for names and
  # for the informational listing.
  if (is.null(group_table)) {
    if (grepl("^[0-9]+$", group)) {
      return(list(index = as.integer(group), name = NA_character_,
                  n_samples = NA_integer_))
    }
    stop("--group \"", group, "\" is a group NAME, but this workspace's group ",
         "list could not be read", ctx, ", so the name cannot be resolved. ",
         "Use the group index instead (e.g. --group 1).")
  }

  avail <- paste0("  ", group_table$index, "  \"", group_table$groupName,
                  "\"  (", group_table$n_samples, " sample(s))",
                  collapse = "\n")

  if (grepl("^[0-9]+$", group)) {
    idx <- as.integer(group)
    if (idx < 1L || idx > nrow(group_table)) {
      stop("--group ", idx, " is out of range", ctx, " (it has ",
           nrow(group_table), " group(s)). Available:\n", avail)
    }
  } else {
    idx <- match(group, group_table$groupName)
    if (is.na(idx)) {
      stop("--group \"", group, "\" is not a sample group", ctx,
           ". Available:\n", avail)
    }
  }

  if (group_table$n_samples[idx] == 0L) {
    stop("--group \"", group_table$groupName[idx], "\" is empty", ctx, ".")
  }
  list(index = as.integer(group_table$index[idx]),
       name  = group_table$groupName[idx],
       n_samples = group_table$n_samples[idx])
}

# How to render a spec's group in a log line. The name is NA when the group
# accessor was unavailable and an index was used, so it is dropped rather than
# printed as the string "NA".
describe_group <- function(spec) {
  if (is.null(spec$group_name) || is.na(spec$group_name)) {
    sprintf("group %d", spec$group_index)
  } else {
    sprintf("group %d \"%s\"", spec$group_index, spec$group_name)
  }
}

# -----------------------------------------------------------------------------
# Build one spec per workspace: which group, which files, where they live.
#
# Deliberately does NOT stop on a resolution failure — it collects everything so
# the caller can report every problem in every workspace at once (see
# format_resolution_report). `open_ws` and `read_samples` are injectable so the
# orchestration is testable without CytoML.
#
# `allow_missing` = TRUE relaxes ONE failure mode: a sample whose FCS file is
# nowhere under the root is dropped and the workspace loads without it, rather
# than failing the whole run. Real workspace collections routinely reference
# negative controls or re-exports that were never copied to the analysis host,
# and one such file should not cost the other 32 workspaces.
#
# It deliberately does NOT relax the others:
#   * AMBIGUOUS names stay fatal. Missing is a known unknown — the sample is
#     absent and the report says so. Ambiguous means two files could be the
#     sample, and picking one silently puts the WRONG cells in the embedding.
#   * A workspace that fails to OPEN, or a bad --group, stays fatal: that is a
#     broken input, not a gap in the data.
#
# Returns list(specs = per-workspace specs, problems = per-workspace problem
# records, dropped = per-workspace records of what `allow_missing` discarded).
# A spec carries:
#   label, wsp, group_index, group_name, files (data.frame sampleID/name/file),
#   dirs (the unique directories its files live in)
# A dropped record carries: label, wsp, samples (the dropped names), kept, total.
# A workspace whose samples ALL go missing yields a dropped record with kept = 0
# and no spec — it is reported as dropped, not silently absent.
# -----------------------------------------------------------------------------
resolve_workspace_specs <- function(wsp_paths, labels, fcs_index, group,
                                    open_ws = NULL, read_samples = NULL,
                                    allow_missing = FALSE) {
  if (is.null(open_ws))      open_ws      <- wsp_open
  if (is.null(read_samples)) read_samples <- function(ws, idx) wsp_sample_files(ws, idx)

  specs    <- list()
  problems <- list()
  dropped  <- list()

  # One workspace -> a spec and/or a problem/dropped record. Kept as an explicit
  # return value (rather than assigning into `specs`/`problems` from inside
  # tryCatch) so the error handler and the happy path write to the same place.
  one <- function(path, label) {
    tryCatch({
      ws  <- open_ws(path)
      # Optional: only needed to resolve a group NAME and to report the group
      # list. An index-based run carries on without it (see resolve_group).
      tbl <- tryCatch(wsp_group_table(ws), error = function(e) NULL)
      grp <- resolve_group(group, tbl, where = label)
      smp <- read_samples(ws, grp$index)
      res <- resolve_fcs_names(smp$name, fcs_index)

      fatal <- length(res$ambiguous) > 0 ||
               (length(res$missing) > 0 && !allow_missing)
      if (fatal) {
        return(list(problem = list(label = label, wsp = path,
                                   missing = res$missing,
                                   ambiguous = res$ambiguous)))
      }

      out <- list()
      if (length(res$missing)) {
        keep <- !is.na(res$file)
        out$dropped <- list(label = label, wsp = path,
                            samples = res$missing,
                            kept = sum(keep), total = length(keep))
        smp <- smp[keep, , drop = FALSE]
        res$file <- res$file[keep]
        # Nothing survived: report the workspace as dropped and emit no spec.
        if (!nrow(smp)) return(out)
      }
      smp$file <- res$file
      out$spec <- list(label = label, wsp = path,
                       group_index = grp$index, group_name = grp$name,
                       files = smp, dirs = unique(dirname(res$file)),
                       # `partial` is what build_workspace_gatingset() keys on to
                       # name the samples it wants. Without it a partial spec
                       # would be handed a DIRECTORY and CytoML would go looking
                       # for the very files that were just dropped.
                       partial = !is.null(out$dropped),
                       groups = tbl)
      out
    }, error = function(e) {
      list(problem = list(label = label, wsp = path,
                          error = conditionMessage(e)))
    })
  }

  for (i in seq_along(wsp_paths)) {
    got <- one(wsp_paths[i], labels[i])
    if (!is.null(got$spec))    specs[[labels[i]]]    <- got$spec
    if (!is.null(got$problem)) problems[[labels[i]]] <- got$problem
    if (!is.null(got$dropped)) dropped[[labels[i]]]  <- got$dropped
  }

  list(specs = specs, problems = problems, dropped = dropped)
}

# -----------------------------------------------------------------------------
# One consolidated, actionable error message from the problem records. Returns
# "" when there is nothing wrong, so the caller can `if (nzchar(msg)) stop(msg)`.
# -----------------------------------------------------------------------------
format_resolution_report <- function(problems, fcs_root = "") {
  if (length(problems) == 0) return("")

  lines <- c("could not build the workspace set.", "")
  for (p in problems) {
    lines <- c(lines, sprintf("  %s (%s)", p$label, p$wsp))
    if (!is.null(p$error)) {
      lines <- c(lines, sprintf("    %s", p$error))
    }
    if (length(p$missing)) {
      lines <- c(lines, sprintf("    missing:   %s",
                                paste(p$missing, collapse = ", ")))
    }
    for (nm in names(p$ambiguous)) {
      lines <- c(lines, sprintf("    ambiguous: %s matches", nm),
                 sprintf("                 %s", p$ambiguous[[nm]]))
    }
  }
  # Only offered when EVERY problem is a missing file: the flag does not help
  # with an ambiguous name or a workspace that would not open, and suggesting
  # it there would send you down a dead end.
  only_missing <- all(vapply(problems, function(p) {
    length(p$missing) > 0 && length(p$ambiguous) == 0 && is.null(p$error)
  }, logical(1)))

  lines <- c(lines, "",
             if (nzchar(fcs_root))
               paste0("FCS root searched (recursively): ", fcs_root),
             "Nothing was loaded. Fix the FCS root (or --group) and re-run.",
             if (only_missing)
               paste0("To load the samples that ARE present and drop the rest, ",
                      "re-run with --allow_missing_samples."))
  paste(lines, collapse = "\n")
}

# -----------------------------------------------------------------------------
# Report what --allow_missing_samples discarded. Returns "" when nothing was.
#
# This is a WARNING the run continues past, so it has to be readable at a
# glance: the embedding it produces is missing cells that the workspaces asked
# for, and that fact belongs in the log next to the results, not only in the
# scrollback of whoever started the run. A whole workspace dropped (kept = 0) is
# called out separately, because losing one subject entirely is a different
# thing from losing a negative control.
# -----------------------------------------------------------------------------
format_dropped_report <- function(dropped, max_show = 8L) {
  if (length(dropped) == 0) return("")

  n_samples <- sum(vapply(dropped, function(d) length(d$samples), integer(1)))
  gone      <- Filter(function(d) d$kept == 0L, dropped)

  lines <- c(sprintf(
    "%d sample(s) across %d workspace(s) were NOT FOUND and have been dropped:",
    n_samples, length(dropped)), "")
  for (d in dropped) {
    shown <- d$samples
    extra <- length(shown) - max_show
    if (extra > 0) shown <- c(shown[seq_len(max_show)],
                              sprintf("... and %d more", extra))
    lines <- c(lines,
               sprintf("  %s  (%d of %d sample(s) kept)%s",
                       d$label, d$kept, d$total,
                       if (d$kept == 0L) "  <- WORKSPACE DROPPED ENTIRELY" else ""),
               sprintf("    %s", paste(shown, collapse = ", ")))
  }
  lines <- c(lines, "",
             if (length(gone))
               sprintf("%d workspace(s) contributed NO cells at all: %s",
                       length(gone),
                       paste(vapply(gone, function(d) d$label, character(1)),
                             collapse = ", ")),
             paste0("The embedding below is built WITHOUT those samples. ",
                    "Drop --allow_missing_samples to make this fatal again."))
  paste(lines, collapse = "\n")
}

# -----------------------------------------------------------------------------
# Warn about FCS files claimed by more than one workspace.
#
# This checks FILE IDENTITY, not names: a repeated patient or subject ID across
# workspaces is normal and is NOT what this reports. In a longitudinal design
# the same subject legitimately appears many times, at different timepoints, in
# different files — none of that trips this.
#
# What it reports is the same FCS FILE being pulled in by two workspaces, which
# does happen: a workspace saved twice under slightly different names, or one
# sample re-gated in a combined workspace as well as its own. Each such file
# contributes its cells to the pooled embedding ONCE PER WORKSPACE, silently
# over-weighting that sample.
#
# Deliberately a warning, not an error: re-gating one file in two workspaces is
# a legitimate thing to compare. The caller decides.
# Returns "" when every file belongs to exactly one workspace.
# -----------------------------------------------------------------------------
format_shared_file_report <- function(specs, max_show = 10L) {
  if (length(specs) < 2L) return("")

  owners <- list()
  for (spec in specs) {
    for (f in spec$files$file) {
      owners[[f]] <- c(owners[[f]], spec$label)
    }
  }
  shared <- owners[vapply(owners, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(shared) == 0L) return("")

  lines <- c(sprintf(
    paste0("%d FCS file(s) are referenced by more than one workspace, so their ",
           "cells enter the shared embedding ONCE PER WORKSPACE and those ",
           "samples are over-weighted. This is about the same FILE being used ",
           "twice, not about a subject appearing at several timepoints:"),
    length(shared)))
  for (f in utils::head(names(shared), max_show)) {
    lines <- c(lines, sprintf("    %s  <- %s", basename(f),
                              paste(unique(shared[[f]]), collapse = ", ")))
  }
  if (length(shared) > max_show) {
    lines <- c(lines, sprintf("    ... and %d more",
                              length(shared) - max_show))
  }
  paste(lines, collapse = "\n")
}

# -----------------------------------------------------------------------------
# Build the GatingSet for one spec.
#
# When every one of the workspace's files lives in ONE directory we pass that
# directory, which is the flowjo_to_gatingset() form the rest of this repo
# already uses and is therefore certain to work. Only a workspace whose files
# are split across directories needs the per-sample data.frame form; if the
# installed CytoML does not support it we say so explicitly rather than
# returning a quietly wrong GatingSet.
#
# A PARTIAL spec (--allow_missing_samples dropped some of its samples) must not
# take the directory form even when one directory would do. That form leaves
# CytoML to pair the group's samples with files itself, and the samples we
# dropped are exactly the ones it cannot pair — it would fail on them, which is
# the failure the flag exists to avoid. Such a spec is named sample by sample
# instead, via `subset`, so CytoML is asked only for what we know is there.
#
# The sample count is checked against what we resolved, because
# fj_ws_get_samples(group_id=) and flowjo_to_gatingset(name=) index groups
# independently — a disagreement would otherwise load the wrong group silently.
# -----------------------------------------------------------------------------
build_workspace_gatingset <- function(spec) {
  ws <- wsp_open(spec$wsp)

  if (isTRUE(spec$partial)) {
    cat(sprintf("    %s: loading %d named sample(s); the rest have no FCS file.\n",
                spec$label, nrow(spec$files)))
    map <- data.frame(sampleID = spec$files$sampleID,
                      file     = spec$files$file,
                      stringsAsFactors = FALSE)
    gs <- tryCatch(
      wsp_to_gatingset(ws, group = spec$group_index, path = map,
                       subset = spec$files$name),
      error = function(e) {
        stop("Workspace '", spec$label, "' is missing ",
             "some of its FCS files and this CytoML rejected loading the ",
             "remaining samples by name:\n  ", conditionMessage(e),
             "\nDrop --allow_missing_samples and supply the missing files ",
             "instead, or exclude this workspace from --wsp.")
      })
  } else if (length(spec$dirs) == 1L) {
    gs <- wsp_to_gatingset(ws, group = spec$group_index, path = spec$dirs)
  } else {
    cat(sprintf(
      "    %s: FCS files span %d directories; using a per-sample path map.\n",
      spec$label, length(spec$dirs)))
    map <- data.frame(sampleID = spec$files$sampleID,
                      file     = spec$files$file,
                      stringsAsFactors = FALSE)
    gs <- tryCatch(
      wsp_to_gatingset(ws, group = spec$group_index, path = map),
      error = function(e) {
        stop("Workspace '", spec$label, "' has FCS files in ",
             length(spec$dirs), " different directories, and this CytoML ",
             "rejected the per-sample path map:\n  ", conditionMessage(e),
             "\nPut this workspace's FCS files under one directory, or ",
             "point --fcs_dir at a root where they are not split.\n",
             "Directories involved:\n  ",
             paste(spec$dirs, collapse = "\n  "))
      })
  }

  loaded <- length(flowWorkspace::sampleNames(gs))
  if (loaded != nrow(spec$files)) {
    warning(sprintf(
      paste0("Workspace '%s': %s lists %d sample(s) but the GatingSet loaded ",
             "%d. The workspace's group indexing may differ from the sample ",
             "table; check the group with --list_populations."),
      spec$label, describe_group(spec), nrow(spec$files), loaded))
  }
  gs
}

# -----------------------------------------------------------------------------
# Display label per sample across the WHOLE run.
#
# A sample name is kept exactly as it is when it is unique across every
# workspace — so a single-workspace run produces the same legends and the same
# file names as before this feature existed. Only names that genuinely collide
# are qualified with their workspace, because two same-named samples would
# otherwise merge into one legend entry and one set of per-file plots, with the
# second silently overwriting the first on disk.
#
# Pure function over character vectors. `clean_names` is clean_sample_label()
# output; `wsp_labels` is the owning workspace of each.
# -----------------------------------------------------------------------------
disambiguate_sample_labels <- function(clean_names, wsp_labels) {
  stopifnot(length(clean_names) == length(wsp_labels))
  out <- clean_names
  dup <- clean_names %in% clean_names[duplicated(clean_names)]
  if (any(dup)) {
    out[dup] <- paste0(wsp_labels[dup], ": ", clean_names[dup])
  }
  out
}
