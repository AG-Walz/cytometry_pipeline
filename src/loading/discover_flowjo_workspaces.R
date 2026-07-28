#!/usr/bin/env Rscript
# =============================================================================
# loading/discover_flowjo_workspaces.R — Find and pair many FlowJo workspaces
# =============================================================================
# Turns two FOLDER ROOTS — one holding .wsp workspaces, one holding .fcs files,
# both searched recursively — into a list of ready-to-load workspace specs.
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
# The workspace's sample groups, as a table we control the shape of.
#
# Isolates the CytoML/flowWorkspace accessor so an API difference on the
# analysis host is a one-function fix rather than a hunt. `index` is the 1-based
# position used by fj_ws_get_samples(group_id=) and by
# flowjo_to_gatingset(name=<integer>).
# -----------------------------------------------------------------------------
wsp_group_table <- function(ws) {
  groups <- flowWorkspace::fj_ws_get_sample_groups(ws)
  nm <- unique(as.character(groups$groupName))
  data.frame(
    index     = seq_along(nm),
    groupName = nm,
    n_samples = as.integer(vapply(
      nm, function(g) sum(as.character(groups$groupName) == g), integer(1))),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# The FCS file names of one group, in workspace order.
#
# fj_ws_get_samples() names its file-name column differently across versions, so
# the first plausible column is used rather than a hard-coded one.
# -----------------------------------------------------------------------------
wsp_sample_files <- function(ws, group_index) {
  samples <- flowWorkspace::fj_ws_get_samples(ws, group_id = group_index)
  cols <- c("name", "sampleName", "file", "filename")
  hit  <- cols[cols %in% colnames(samples)]
  if (length(hit) == 0) {
    stop("fj_ws_get_samples() returned no recognisable file-name column ",
         "(saw: ", paste(colnames(samples), collapse = ", "), "). ",
         "Update wsp_sample_files() in ",
         "src/loading/discover_flowjo_workspaces.R.")
  }
  data.frame(
    sampleID = as.character(samples$sampleID),
    name     = as.character(samples[[hit[1]]]),
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

# -----------------------------------------------------------------------------
# Build one spec per workspace: which group, which files, where they live.
#
# Deliberately does NOT stop on a resolution failure — it collects everything so
# the caller can report every problem in every workspace at once (see
# format_resolution_report). `open_ws` and `read_samples` are injectable so the
# orchestration is testable without CytoML.
#
# Returns list(specs = list of per-workspace specs, problems = list of
# per-workspace problem records). A spec carries:
#   label, wsp, group_index, group_name, files (data.frame sampleID/name/file),
#   dirs (the unique directories its files live in)
# -----------------------------------------------------------------------------
resolve_workspace_specs <- function(wsp_paths, labels, fcs_index, group,
                                    open_ws = NULL, read_samples = NULL) {
  if (is.null(open_ws))      open_ws      <- function(p) CytoML::open_flowjo_xml(p)
  if (is.null(read_samples)) read_samples <- function(ws, idx) wsp_sample_files(ws, idx)

  specs    <- list()
  problems <- list()

  # One workspace -> either a spec or a problem record. Kept as an explicit
  # return value (rather than assigning into `specs`/`problems` from inside
  # tryCatch) so the error handler and the happy path write to the same place.
  one <- function(path, label) {
    tryCatch({
      ws  <- open_ws(path)
      tbl <- wsp_group_table(ws)
      grp <- resolve_group(group, tbl, where = label)
      smp <- read_samples(ws, grp$index)
      res <- resolve_fcs_names(smp$name, fcs_index)

      if (length(res$missing) || length(res$ambiguous)) {
        return(list(problem = list(label = label, wsp = path,
                                   missing = res$missing,
                                   ambiguous = res$ambiguous)))
      }
      smp$file <- res$file
      list(spec = list(label = label, wsp = path,
                       group_index = grp$index, group_name = grp$name,
                       files = smp, dirs = unique(dirname(res$file)),
                       groups = tbl))
    }, error = function(e) {
      list(problem = list(label = label, wsp = path,
                          error = conditionMessage(e)))
    })
  }

  for (i in seq_along(wsp_paths)) {
    got <- one(wsp_paths[i], labels[i])
    if (!is.null(got$spec))    specs[[labels[i]]]    <- got$spec
    if (!is.null(got$problem)) problems[[labels[i]]] <- got$problem
  }

  list(specs = specs, problems = problems)
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
  lines <- c(lines, "",
             if (nzchar(fcs_root))
               paste0("FCS root searched (recursively): ", fcs_root),
             "Nothing was loaded. Fix the FCS root (or --group) and re-run.")
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
# The sample count is checked against what we resolved, because
# fj_ws_get_samples(group_id=) and flowjo_to_gatingset(name=) index groups
# independently — a disagreement would otherwise load the wrong group silently.
# -----------------------------------------------------------------------------
build_workspace_gatingset <- function(spec) {
  ws <- CytoML::open_flowjo_xml(spec$wsp)

  if (length(spec$dirs) == 1L) {
    gs <- CytoML::flowjo_to_gatingset(ws, name = spec$group_index,
                                      path = spec$dirs)
  } else {
    cat(sprintf(
      "    %s: FCS files span %d directories; using a per-sample path map.\n",
      spec$label, length(spec$dirs)))
    map <- data.frame(sampleID = spec$files$sampleID,
                      file     = spec$files$file,
                      stringsAsFactors = FALSE)
    gs <- tryCatch(
      CytoML::flowjo_to_gatingset(ws, name = spec$group_index, path = map),
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
      paste0("Workspace '%s': group %d (\"%s\") lists %d sample(s) but the ",
             "GatingSet loaded %d. The workspace's group indexing may differ ",
             "from the sample table; check the group with --list_populations."),
      spec$label, spec$group_index, spec$group_name,
      nrow(spec$files), loaded))
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
