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
#' @param subset Which samples to load: a character vector of sample names, an
#'   integer index, or `list()` (the default) for all of them.
#'
#' `subset` is spliced in by value rather than forwarded, because
#' flowjo_to_gatingset() resolves it with `eval(substitute(subset))` — evaluated
#' in ITS OWN frame, not the caller's. An expression naming any local variable
#' therefore fails to resolve there, silently falls through to a keyword-filter
#' fallback, and dies inside it with "arguments imply differing number of rows".
#' do.call() puts the evaluated vector into the call, so `substitute()` returns a
#' literal and the scope of the caller stops mattering.
wsp_to_gatingset <- function(ws, group = 1L, path, subset = list(), ...) {
  do.call(CytoML::flowjo_to_gatingset,
          c(list(ws, name = group, path = path, subset = subset), list(...)))
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
# `recorded` is the file name the workspace itself recorded for each sample (the
# basename of its <DataSet uri>, see wsp_sample_uris). It is NOT used to match —
# that would change which cells land in the embedding — but when a name fails to
# resolve it is the single most useful thing to be able to say, because the
# usual cause is that FlowJo's sample LABEL was edited and no longer equals the
# file name. It goes into the diagnosis so the log can name the actual file.
#
# Returns list(file = character vector aligned to `names`, NA where unresolved,
# missing = the unmatched names, ambiguous = named list of candidate vectors,
# why = named list of per-missing-name diagnoses; see diagnose_unresolved).
# Pure: no CytoML, no file system access. This is the core the tests exercise.
# -----------------------------------------------------------------------------
path_basename <- function(x) {
  sub(".*[\\\\/]", "", x)
}

# A deliberately blunt key for "is this the same name, give or take punctuation":
# case, spaces, underscores, dashes, dots and parentheses all disappear. Used
# ONLY to explain a failure, never to match — "MTB048_2" and "MTB048 2" collapse
# together here, which is exactly what you want in a hint and exactly what you
# do not want when deciding whose cells go in the plot.
loose_key <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", sub("\\.fcs$", "", x, ignore.case = TRUE)))
}

# fcs_index keyed by loose_key() instead of by exact lower-cased name.
loose_fcs_index <- function(fcs_index) {
  files <- unlist(fcs_index, use.names = FALSE)
  if (!length(files)) return(list())
  split(files, loose_key(basename(files)))
}

# Why one sample name did not resolve, in the terms someone reading the log
# needs: is the file there under another name, is something like it there, or is
# it simply not present? Returns list(kind, candidates), kind being one of:
#   "recorded" - the workspace's own recorded file name IS in the FCS folder.
#                Almost always an edited FlowJo label. candidates = that file.
#   "near"     - a file differing only in punctuation/spacing exists.
#   "absent"   - nothing of that name or anything like it.
diagnose_unresolved <- function(name, recorded, fcs_index, loose_index) {
  if (!is.null(recorded) && !is.na(recorded) && nzchar(recorded)) {
    hit <- fcs_index[[tolower(recorded)]]
    if (!is.null(hit)) {
      return(list(kind = "recorded", candidates = hit, recorded = recorded))
    }
  }
  near <- loose_index[[loose_key(name)]]
  if (!is.null(near)) {
    return(list(kind = "near", candidates = near, recorded = recorded))
  }
  list(kind = "absent", candidates = character(0), recorded = recorded)
}

resolve_fcs_names <- function(names, fcs_index, recorded = NULL,
                              loose_index = NULL) {
  file      <- rep(NA_character_, length(names))
  missing   <- character(0)
  ambiguous <- list()
  why       <- list()

  # The RECORDED file name is tried before the sample name. A FlowJo sample name
  # is an editable label; the recorded name is the file FlowJo actually opened
  # when the gates were drawn. Where they disagree AND both exist on disk, the
  # label is the wrong one to follow: in a real workspace here that meant
  # pairing a 2025 export with a 2024 file of the same name and a different
  # event count — a different acquisition entirely.
  via <- rep(NA_character_, length(names))

  lookup <- function(cand) {
    if (is.null(cand) || is.na(cand) || !nzchar(cand)) return(NULL)
    base <- path_basename(cand)
    keys <- unique(tolower(c(cand, paste0(cand, ".fcs"), base,
                             paste0(base, ".fcs"))))
    for (k in keys) if (!is.null(fcs_index[[k]])) return(fcs_index[[k]])
    NULL
  }

  for (i in seq_along(names)) {
    nm  <- names[i]
    rec <- if (is.null(recorded)) NA_character_ else recorded[i]
    hit <- lookup(rec)
    src <- "recorded"
    if (is.null(hit)) { hit <- lookup(nm); src <- "name" }

    if (is.null(hit)) {
      missing <- c(missing, nm)
      if (is.null(loose_index)) loose_index <- loose_fcs_index(fcs_index)
      why[[nm]] <- diagnose_unresolved(nm, rec, fcs_index, loose_index)
    } else if (length(hit) > 1L) {
      ambiguous[[nm]] <- hit
    } else {
      file[i] <- hit
      via[i]  <- src
    }
  }
  list(file = file, missing = missing, ambiguous = ambiguous, why = why,
       via = via)
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
# The file name the WORKSPACE recorded for each sample, keyed by sampleID.
#
# fj_ws_get_samples() returns sampleID/name/count only — the file is not among
# them — so this reads the workspace XML directly: every <DataSet> carries both
# a uri and its sampleID, which is an exact pairing and not a guess from
# document order. The uri itself is the path on the Windows machine that built
# the workspace and is dead here; only its BASE NAME is of any use.
#
# Purely diagnostic. Nothing pairs samples to FCS files on this — see
# resolve_fcs_names — because a FlowJo sample name and its recorded file name
# can legitimately differ and only the name is the documented key. It exists so
# a failure can say "the file it wants is called X", which is the difference
# between a biologist fixing the data and filing a ticket.
#
# Returns a named character vector (names = sampleID), or an empty one when the
# workspace cannot be read as XML — a diagnosis is a nicety and must never be
# the reason a run fails.
# -----------------------------------------------------------------------------
wsp_recorded_files <- function(wsp_path) {
  empty <- character(0)
  if (!requireNamespace("xml2", quietly = TRUE)) return(empty)
  out <- tryCatch({
    doc  <- xml2::read_xml(wsp_path)
    sets <- xml2::xml_find_all(doc, "//DataSet")
    if (!length(sets)) return(empty)
    ids  <- xml2::xml_attr(sets, "sampleID")
    uris <- xml2::xml_attr(sets, "uri")
    keep <- !is.na(ids) & !is.na(uris)
    if (!any(keep)) return(empty)
    # A uri is percent-encoded ("MTB%20Manuskript") and may use either
    # separator. URLdecode() is vectorised on current R but was not always, so
    # it is applied one at a time rather than silently decoding only uris[1].
    base <- path_basename(vapply(uris[keep], utils::URLdecode, character(1),
                                 USE.NAMES = FALSE))
    stats::setNames(base, ids[keep])
  }, error = function(e) empty)
  if (is.null(out)) empty else out
}

# -----------------------------------------------------------------------------
# Which of the discovered workspaces this run should use.
#
# `spec` is an INCLUDE list: workspace labels, either as a character vector or
# comma-separated in one string. Empty or NULL means all of them, so the option
# costs nothing when nobody sets it.
#
# Why this exists rather than reusing --population_map, which already excludes a
# workspace that has no entry: the map requires a GATE per workspace. Wanting to
# leave one workspace out of a run is not the same as wanting to name a
# different gate in each of the other twenty-two, and forcing the second in
# order to express the first is how a one-click exclusion turns into a
# twenty-two-line configuration file. The two compose instead — select the
# workspaces here, then choose their gates however you like.
#
# Applied BEFORE any GatingSet is built, so deselecting a workspace makes the
# load faster rather than merely hiding it at the end.
#
# An unknown label is fatal and lists what was available: silently ignoring it
# would run on more workspaces than the caller asked for, which is the one
# outcome a selection must never produce.
#
# Returns a logical vector over `labels`.
# -----------------------------------------------------------------------------
select_workspaces <- function(labels, spec = NULL) {
  if (is.null(spec)) return(rep(TRUE, length(labels)))
  wanted <- if (length(spec) == 1L && is.character(spec)) {
    trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  } else {
    trimws(as.character(spec))
  }
  wanted <- wanted[nzchar(wanted)]
  if (!length(wanted)) return(rep(TRUE, length(labels)))

  unknown <- setdiff(wanted, labels)
  if (length(unknown)) {
    stop("--workspaces names ", length(unknown),
         " workspace(s) that were not found: ",
         paste(unknown, collapse = ", "),
         "\nAvailable: ", paste(labels, collapse = ", "), call. = FALSE)
  }
  labels %in% wanted
}

# -----------------------------------------------------------------------------
# A readable label for one gate path.
#
# FlowJo lets a gate be named with nothing but spaces, and real workspaces
# contain them — two siblings called " " and "  " sit under "Live cells" in the
# study this was written for. Rendered verbatim they are invisible rows in a
# tree: impossible to see, impossible to tell apart, and impossible to click on
# with any confidence about which one you got.
#
# The PATH is never changed — it is what --population has to resolve against, so
# it stays byte-for-byte what the workspace said. Only the label shown to a
# human is substituted, and it states the width so two blank siblings remain
# distinguishable.
# -----------------------------------------------------------------------------
gate_display_label <- function(path) {
  vapply(path, function(p) {
    leaf <- sub(".*/", "", p)
    if (nzchar(trimws(leaf))) return(leaf)
    n <- nchar(leaf)
    if (n == 0L) "(unnamed gate)"
    else sprintf("(blank name: %d space%s)", n, if (n == 1L) "" else "s")
  }, character(1), USE.NAMES = FALSE)
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
# `match_mode` decides what a FILE NAME is worth:
#
#   "strict"  — the original contract, above. A name that does not resolve is
#               the end of the matter for that sample.
#   "content" — a name that does not resolve is not a verdict. The workspace is
#               handed to CytoML anyway, with the FCS ROOT as its search path,
#               and greedy_match pairs each sample by the $FIL keyword stored
#               INSIDE the file plus its event count. That is a CONTENT match:
#               it survives a rename on disk, which is the single most common
#               reason a real cohort will not load. It is also stricter than
#               matching names, not looser — a name can be edited to anything,
#               while ($FIL, $TOT) has to agree with what FlowJo gated.
#
#               Names stay the fast path: anything they resolve is used as-is
#               and recorded as such, so the usual case does not change. Only
#               the leftovers are deferred to CytoML. Ambiguous names are also
#               deferred rather than fatal here, because ambiguity is a property
#               of the NAME and content is what decides.
#
# Returns list(specs = per-workspace specs, problems = per-workspace problem
# records, dropped = per-workspace records of what `allow_missing` discarded,
# deferred = per-workspace records of what "content" mode left to CytoML).
# A spec carries:
#   label, wsp, group_index, group_name, files (data.frame sampleID/name/file/
#   via), dirs (the directories to search)
# A dropped record carries: label, wsp, samples (the dropped names), kept, total.
# A workspace whose samples ALL go missing yields a dropped record with kept = 0
# and no spec — it is reported as dropped, not silently absent. In "content"
# mode it yields a spec instead, because "no name matched" and "no cells are
# there" are different claims and only CytoML can tell them apart.
#
# `fcs_root` is only consulted in "content" mode, where it is the search path a
# workspace with unresolved names needs. It defaults to "" so every existing
# caller and test keeps working unchanged.
# -----------------------------------------------------------------------------
resolve_workspace_specs <- function(wsp_paths, labels, fcs_index, group,
                                    open_ws = NULL, read_samples = NULL,
                                    allow_missing = FALSE,
                                    read_recorded = NULL,
                                    match_mode = "strict", fcs_root = "") {
  match_mode <- match.arg(match_mode, c("strict", "content"))
  by_content <- identical(match_mode, "content")
  if (by_content && !nzchar(fcs_root)) {
    stop("match_mode = \"content\" needs fcs_root: a workspace whose names do ",
         "not resolve has no directory of its own to search.")
  }
  if (is.null(open_ws))      open_ws      <- wsp_open
  if (is.null(read_samples)) read_samples <- function(ws, idx) wsp_sample_files(ws, idx)
  if (is.null(read_recorded)) read_recorded <- wsp_recorded_files

  specs    <- list()
  problems <- list()
  dropped  <- list()
  deferred <- list()
  # Built once: every workspace diagnoses against the same FCS folder.
  loose <- loose_fcs_index(fcs_index)

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
      # Diagnostic only, and never fatal: a workspace that will not parse as XML
      # simply yields no recorded names and the log says less.
      rec_map  <- tryCatch(read_recorded(path), error = function(e) character(0))
      recorded <- unname(rec_map[smp$sampleID])
      res <- resolve_fcs_names(smp$name, fcs_index, recorded = recorded,
                               loose_index = loose)

      # In "content" mode nothing about a NAME is fatal and nothing is dropped:
      # the unresolved samples are handed to CytoML, which pairs them on $FIL +
      # event count. Only a workspace that would not open or a bad --group
      # (handled by the error branch below) can still fail here.
      fatal <- !by_content &&
               (length(res$ambiguous) > 0 ||
                (length(res$missing) > 0 && !allow_missing))
      if (fatal) {
        return(list(problem = list(label = label, wsp = path,
                                   missing = res$missing,
                                   ambiguous = res$ambiguous,
                                   why = res$why,
                                   total = nrow(smp))))
      }

      # The group's FULL sample listing, before anything is dropped. Kept so the
      # load record can compare what arrived against what the workspace actually
      # asked for, rather than against the trimmed set we chose to look for.
      listed <- smp$name

      out <- list()
      if (by_content) {
        unresolved <- unique(c(res$missing, names(res$ambiguous)))
        if (length(unresolved)) {
          out$deferred <- list(label = label, wsp = path,
                               samples = unresolved, why = res$why,
                               ambiguous = res$ambiguous,
                               resolved = sum(!is.na(res$file)),
                               total = nrow(smp))
        }
      } else if (length(res$missing)) {
        keep <- !is.na(res$file)
        out$dropped <- list(label = label, wsp = path,
                            samples = res$missing, why = res$why,
                            kept = sum(keep), total = length(keep))
        smp <- smp[keep, , drop = FALSE]
        res$file <- res$file[keep]
        res$via  <- res$via[keep]
        # Nothing survived: report the workspace as dropped and emit no spec.
        if (!nrow(smp)) return(out)
      }
      smp$file <- res$file
      # How each sample was paired, so the run can say so afterwards instead of
      # leaving "it loaded" as the only fact. NA here means no name matched;
      # whether that sample arrives at all is then CytoML's answer to give.
      smp$via  <- res$via

      # Which directories CytoML is pointed at. Resolved files supply their own;
      # "content" mode adds the FCS root, because the samples whose names did
      # not resolve could be anywhere under it — and a workspace where NOTHING
      # resolved would otherwise have no search path at all.
      dirs <- unique(dirname(res$file[!is.na(res$file)]))
      if (by_content) dirs <- unique(c(dirs, fcs_root))

      out$spec <- list(label = label, wsp = path,
                       group_index = grp$index, group_name = grp$name,
                       files = smp, listed = listed, dirs = dirs,
                       # `partial` is what build_workspace_gatingset() keys on to
                       # name the samples it wants. Without it a partial spec
                       # would be handed a DIRECTORY and CytoML would go looking
                       # for the very files that were just dropped.
                       partial = !is.null(out$dropped),
                       match_mode = match_mode,
                       groups = tbl)
      out
    }, error = function(e) {
      list(problem = list(label = label, wsp = path,
                          error = conditionMessage(e)))
    })
  }

  for (i in seq_along(wsp_paths)) {
    got <- one(wsp_paths[i], labels[i])
    if (!is.null(got$spec))     specs[[labels[i]]]    <- got$spec
    if (!is.null(got$problem))  problems[[labels[i]]] <- got$problem
    if (!is.null(got$dropped))  dropped[[labels[i]]]  <- got$dropped
    if (!is.null(got$deferred)) deferred[[labels[i]]] <- got$deferred
  }

  list(specs = specs, problems = problems, dropped = dropped,
       deferred = deferred)
}

# -----------------------------------------------------------------------------
# Explain, sample by sample, why a set of names did not resolve.
#
# The audience is the person who owns the data, not the person who wrote this:
# they need to know which samples are affected, whether the cells still exist,
# and what to do — so the samples are grouped by CAUSE and every one is named.
# Nothing is ever elided with "and N more": a biologist cannot act on a sample
# they have not been told about, and a list too long to read is itself the
# finding.
#
# `why` is the diagnosis map from resolve_fcs_names(); a name missing from it
# degrades to "absent", which is the safe thing to claim.
# -----------------------------------------------------------------------------
.explain_unresolved <- function(samples, why, ind = "    ") {
  kinds <- vapply(samples, function(nm) {
    w <- why[[nm]]; if (is.null(w$kind)) "absent" else w$kind
  }, character(1))
  fmt <- function(nm) {
    w <- why[[nm]]
    if (length(w$candidates)) basename(w$candidates[1]) else NA_character_
  }
  out <- character(0)
  pad <- max(c(0L, nchar(samples))) + 2L

  rec <- samples[kinds == "recorded"]
  if (length(rec)) {
    out <- c(out,
      sprintf("%s%d sample(s): the FCS file IS in your folder, under the name this",
              ind, length(rec)),
      sprintf("%sworkspace recorded for it. The sample's LABEL in FlowJo was edited and",
              ind),
      sprintf("%sno longer equals the file name, and samples are matched by label.", ind),
      sprintf("%s  %-*s %s", ind, pad, "FlowJo label", "file sitting in your FCS folder"),
      vapply(rec, function(nm)
        sprintf("%s  %-*s %s", ind, pad, nm, fmt(nm)), character(1)),
      sprintf("%s  -> These cells are NOT lost. Rename the file to the label (or fix the", ind),
      sprintf("%s     label in FlowJo) and they will load. Check first whether another", ind),
      sprintf("%s     workspace already uses the file under its current name.", ind))
  }

  near <- samples[kinds == "near"]
  if (length(near)) {
    out <- c(out,
      sprintf("%s%d sample(s): no exact match, but a file differing only in spacing or",
              ind, length(near)),
      sprintf("%spunctuation is present. Probably the same file, renamed.", ind),
      sprintf("%s  %-*s %s", ind, pad, "FlowJo label", "closest file in your FCS folder"),
      vapply(near, function(nm)
        sprintf("%s  %-*s %s", ind, pad, nm, fmt(nm)), character(1)),
      sprintf("%s  -> Confirm they are the same acquisition, then make the names agree.", ind))
  }

  abs_ <- samples[kinds == "absent"]
  if (length(abs_)) {
    out <- c(out,
      sprintf("%s%d sample(s): NOT in your FCS folder under this name, and nothing",
              ind, length(abs_)),
      sprintf("%sresembling it is either. These cells are genuinely unavailable:", ind),
      vapply(abs_, function(nm) {
        w <- why[[nm]]
        r <- w$recorded
        if (!is.null(r) && !is.na(r) && nzchar(r) && !identical(r, nm))
          sprintf("%s  %-*s (workspace records the file as %s)", ind, pad, nm, r)
        else sprintf("%s  %s", ind, nm)
      }, character(1)),
      sprintf("%s  -> Copy these files into the FCS folder, or accept that these", ind),
      sprintf("%s     samples are not in the analysis.", ind))
  }
  out
}

# -----------------------------------------------------------------------------
# One consolidated, actionable error message from the problem records. Returns
# "" when there is nothing wrong, so the caller can `if (nzchar(msg)) stop(msg)`.
# -----------------------------------------------------------------------------
format_resolution_report <- function(problems, fcs_root = "") {
  if (length(problems) == 0) return("")

  lines <- c("COULD NOT LOAD THE WORKSPACES — nothing was analysed.", "")
  for (p in problems) {
    n_miss <- length(p$missing) + length(p$ambiguous)
    head <- if (!is.null(p$total) && p$total > 0)
      sprintf("  %s — %d of %d sample(s) could not be used", p$label, n_miss, p$total)
    else sprintf("  %s", p$label)
    lines <- c(lines, head, sprintf("    workspace file: %s", p$wsp))
    if (!is.null(p$error)) {
      lines <- c(lines, sprintf("    this workspace could not be read at all:"),
                 sprintf("      %s", p$error))
    }
    if (length(p$missing)) {
      lines <- c(lines, .explain_unresolved(p$missing, p$why))
    }
    if (length(p$ambiguous)) {
      lines <- c(lines, sprintf(
        "    %d sample(s): the name matches MORE THAN ONE file, so which cells are",
        length(p$ambiguous)),
        "    meant is undecidable. Nothing is guessed — remove or rename the copies:")
      for (nm in names(p$ambiguous)) {
        lines <- c(lines, sprintf("      %s matches", nm),
                   sprintf("        %s", p$ambiguous[[nm]]))
      }
    }
    lines <- c(lines, "")
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
format_dropped_report <- function(dropped) {
  if (length(dropped) == 0) return("")

  n_samples <- sum(vapply(dropped, function(d) length(d$samples), integer(1)))
  n_total   <- sum(vapply(dropped, function(d) d$total, integer(1)))
  n_kept    <- sum(vapply(dropped, function(d) d$kept, integer(1)))
  gone      <- Filter(function(d) d$kept == 0L, dropped)

  # Counts are scoped to the AFFECTED workspaces and say so. Summing kept across
  # only these and printing it as "loaded" reads as the run total, which it is
  # not — the untouched workspaces are not in this report at all.
  lines <- c(
    sprintf("SAMPLES MISSING FROM THIS ANALYSIS — %d sample(s) dropped.", n_samples),
    sprintf("In the %d affected workspace(s), %d of %d sample(s) loaded; the other %d could not be",
            length(dropped), n_kept, n_total, n_samples),
    "used and are listed below, each with the reason.", "")

  for (d in dropped) {
    flag <- if (d$kept == 0L)
      "   <<< NOTHING FROM THIS WORKSPACE IS IN THE UMAP" else ""
    lines <- c(lines,
               sprintf("  %s — %d of %d sample(s) loaded, %d dropped%s",
                       d$label, d$kept, d$total, length(d$samples), flag),
               .explain_unresolved(d$samples, d$why))
    # A workspace where nothing resolved and nothing even resembles its files is
    # a different diagnosis from scattered losses, and worth saying outright.
    kinds <- vapply(d$samples, function(nm) {
      w <- d$why[[nm]]; if (is.null(w$kind)) "absent" else w$kind
    }, character(1))
    if (d$kept == 0L && all(kinds == "absent")) {
      lines <- c(lines, sprintf(
        "    -> No file of this workspace is present under any name: the whole set"),
        "       is missing, rather than individual files having been lost.")
    }
    lines <- c(lines, "")
  }

  lines <- c(lines,
             if (length(gone))
               sprintf("%d workspace(s) contributed no cells at all: %s",
                       length(gone),
                       paste(vapply(gone, function(d) d$label, character(1)),
                             collapse = ", ")),
             paste0("The plots below were built WITHOUT the dropped samples. ",
                    "Turn off \"Load anyway when FCS"),
             paste0("files are missing\" (--allow_missing_samples) to make this ",
                    "stop the run instead."))
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
    # NA in "content" mode means no NAME matched, so there is no file to
    # attribute yet — CytoML has not been asked. Counting those would group
    # every such sample under one "NA" key and invent an overlap that is really
    # just "not looked up".
    for (f in spec$files$file[!is.na(spec$files$file)]) {
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
# `path` MUST be a character scalar on CytoML 2.14: it calls normalizePath() on
# it, so the per-sample data.frame form this used to pass is rejected outright
# ("invalid 'path' argument"). It searches RECURSIVELY, so a workspace whose
# files are split across directories is served by their common parent.
#
# Restricting WHICH samples load is left to CytoML, which already skips a sample
# whose file is not under `path` — verified on CytoML 2.14. The previous attempt
# to name them via `subset=` could not work: flowjo_to_gatingset() evaluates
# that argument with eval(substitute(subset)) in its OWN frame, so any variable
# of ours is out of scope and it falls through to "invalid 'subset' argument!".
#
# greedy_match=TRUE lets CytoML accept a file whose NAME differs from the sample
# label when the file's $FIL keyword matches it AND the event count is right.
# That is content-verified, not a guess on names, and it is what makes a
# workspace load whose analyst renamed the files after gating.
#
# The sample count is checked against what we resolved, because
# fj_ws_get_samples(group_id=) and flowjo_to_gatingset(name=) index groups
# independently — a disagreement would otherwise load the wrong group silently.
# -----------------------------------------------------------------------------

# The deepest directory containing all of `dirs`. Used because `path` takes one
# scalar and searches under it.
common_parent <- function(dirs) {
  if (length(dirs) == 1L) return(dirs)
  parts <- strsplit(dirs, .Platform$file.sep, fixed = TRUE)
  n <- min(vapply(parts, length, integer(1)))
  keep <- 0L
  for (i in seq_len(n)) {
    seg <- unique(vapply(parts, function(p) p[i], character(1)))
    if (length(seg) > 1L) break
    keep <- i
  }
  if (keep == 0L) return(.Platform$file.sep)
  out <- paste(parts[[1]][seq_len(keep)], collapse = .Platform$file.sep)
  # An absolute path splits with an empty first segment, so sharing only the
  # root collapses to "" — which is not a directory anyone can search.
  if (!nzchar(out)) .Platform$file.sep else out
}

build_workspace_gatingset <- function(spec) {
  ws   <- wsp_open(spec$wsp)
  root <- common_parent(spec$dirs)
  if (length(spec$dirs) > 1L) {
    cat(sprintf("    %s: FCS files span %d directories; searching under %s\n",
                spec$label, length(spec$dirs), root))
  }
  wsp_to_gatingset(ws, group = spec$group_index, path = root,
                   greedy_match = TRUE)
}

# -----------------------------------------------------------------------------
# One workspace, one channel set — and what to do when that is not true.
#
# CytoML materialises a sample group as a single cytoset: one HDF5-backed
# collection with ONE column schema. flowjo_to_gatingset() therefore refuses a
# group whose FCS files do not all carry exactly the same $PnN set, with
# "Found channel inconsistency across samples", and the whole workspace is lost.
#
# FlowJo never needs that. It holds every sample as its own object and resolves
# each gate's dimension against that sample's own parameter list, so a workspace
# that opens perfectly in FlowJo can be unloadable here. The divergence is not
# corruption and not a mistake — a cohort acquired over two years on two
# instruments really does have:
#
#   * the viability detector named "Aqua LiveDead-A" on the FACSCanto II and
#     "BV480-A" on the FACSymphony A1 — same dye, same $PnS, different name in
#     the cytometer configuration;
#   * the scatter-width parameters FSC-W and SSC-W recorded under some
#     acquisition templates and not others ($PAR 14 vs 12);
#   * occasional extra detectors ("SP SSC-A") from a third template.
#
# So the samples are split into channel-homogeneous groups and one GatingSet is
# built per group. Nothing is renamed and nothing is dropped: each sample keeps
# the channel names its gates were drawn against, which is the only way the
# gates stay the gates FlowJo applied. The parts rejoin downstream, where cells
# are pooled on MARKER — the channel dictionary already maps both viability
# channel names onto the same marker — and not on channel.
# -----------------------------------------------------------------------------

#' Is this CytoML's cross-sample channel check talking?
is_channel_inconsistency <- function(message) {
  grepl("channel inconsistency", message, ignore.case = TRUE)
}

#' The set of channels an FCS file carries, as one comparable string.
#'
#' Reads the TEXT segment only — no events — so this costs a stat and a short
#' read per file even on a 400 MB acquisition. NA when the header cannot be
#' read or does not describe its own parameters, which the caller treats as
#' "unknown" rather than as a group of its own.
fcs_channel_signature <- function(path) {
  tryCatch({
    kw <- flowCore::read.FCSheader(path)[[1]]
    n  <- suppressWarnings(as.integer(kw[["$PAR"]]))
    if (is.na(n) || n < 1L) return(NA_character_)
    nms <- vapply(seq_len(n), function(i) {
      v <- kw[[paste0("$P", i, "N")]]
      if (is.null(v) || is.na(v)) NA_character_ else as.character(v)
    }, character(1))
    if (anyNA(nms)) return(NA_character_)
    # Sorted: parameter ORDER also differs between templates and is irrelevant
    # here — flowWorkspace matches channels by name, not by index.
    paste(sort(nms), collapse = "|")
  }, error = function(e) NA_character_)
}

#' Channel signature of every sample a spec lists.
#'
#' Keyed on the LISTED sample names, not on the resolved ones: the subset handed
#' to CytoML has to name every sample the group contains, or content matching
#' loses exactly the samples it exists to rescue. Samples whose file was never
#' resolved by name have no signature to read and come back NA.
workspace_channel_signatures <- function(spec) {
  listed <- if (is.null(spec$listed)) spec$files$name else spec$listed
  paths  <- stats::setNames(as.character(spec$files$file), spec$files$name)
  vapply(listed, function(nm) {
    # A listed sample the spec never resolved is not IN `paths` at all, so the
    # lookup has to tolerate a miss rather than index out of bounds.
    p <- unname(paths[nm])
    if (is.na(p) || !nzchar(p) || !file.exists(p)) return(NA_character_)
    fcs_channel_signature(p)
  }, character(1), USE.NAMES = TRUE)
}

#' Split sample names into channel-homogeneous groups, largest first.
#'
#' Pure over a named character vector, so it is testable without the cytometry
#' stack. Samples of unknown signature (NA) join the largest known group: they
#' have to go somewhere, and the largest group is both the likeliest and the one
#' whose failure costs the least, since a part that still refuses to build only
#' loses its own samples.
#'
#' @return list(groups = <list of character vectors>, signature = <character>).
channel_signature_groups <- function(signatures) {
  nms   <- names(signatures)
  known <- !is.na(signatures)
  if (!any(known)) {
    return(list(groups = list(nms), signature = NA_character_))
  }
  grouped <- split(nms[known], signatures[known])
  # Largest first; ties broken on the signature text so the split is stable
  # across runs and the "part 1" everything is compared against never moves.
  ord     <- order(-lengths(grouped), names(grouped))
  grouped <- grouped[ord]
  sigs    <- names(grouped)
  if (any(!known)) grouped[[1]] <- c(grouped[[1]], nms[!known])
  list(groups = unname(grouped), signature = sigs)
}

#' Say what the split is, in terms of the channels that differ.
#'
#' Pure; returns "" when there is nothing to split, so the caller can print it
#' unconditionally. Every later part is described against part 1 rather than
#' listing four full panels, because what a person needs is the difference.
format_panel_split <- function(label, split) {
  if (length(split$groups) < 2L) return("")
  ref <- strsplit(split$signature[1], "|", fixed = TRUE)[[1]]
  n   <- length(unlist(split$groups, use.names = FALSE))
  out <- c(sprintf(
    paste0("%s: its %d sample(s) were acquired with %d different channel sets, ",
           "so CytoML\n  refused the workspace whole. Loading it as %d parts, ",
           "each with its own\n  channels; they rejoin in the UMAP on marker ",
           "names."),
    label, n, length(split$groups), length(split$groups)))
  for (i in seq_along(split$groups)) {
    g <- split$groups[[i]]
    out <- c(out, sprintf("    part %d — %d sample(s), e.g. %s",
                          i, length(g), g[1]))
    if (i > 1L && !is.na(split$signature[i])) {
      cur   <- strsplit(split$signature[i], "|", fixed = TRUE)[[1]]
      miss  <- setdiff(ref, cur)
      extra <- setdiff(cur, ref)
      if (length(miss)) {
        out <- c(out, paste0("        lacks: ", paste(miss, collapse = ", ")))
      }
      if (length(extra)) {
        out <- c(out, paste0("        adds:  ", paste(extra, collapse = ", ")))
      }
    }
  }
  paste(out, collapse = "\n")
}

#' Build the GatingSet(s) a workspace needs — one normally, more when its
#' samples do not share a channel set.
#'
#' Always returns a NAMED list, so the caller has one shape to handle. The
#' single-part name is the workspace label itself, which keeps every existing
#' run byte-identical; parts are suffixed "[panel k]" and are the only case
#' where a workspace appears more than once downstream.
#'
#' Failures of individual parts are attached as `attr(, "failed_parts")` rather
#' than raised: losing one acquisition template is not a reason to lose the
#' other three.
build_workspace_gatingsets <- function(spec) {
  ws   <- wsp_open(spec$wsp)
  root <- common_parent(spec$dirs)
  if (length(spec$dirs) > 1L) {
    cat(sprintf("    %s: FCS files span %d directories; searching under %s\n",
                spec$label, length(spec$dirs), root))
  }

  whole <- tryCatch(
    list(gs = wsp_to_gatingset(ws, group = spec$group_index, path = root,
                               greedy_match = TRUE)),
    error = function(e) list(error = e))
  if (!is.null(whole$gs)) {
    return(stats::setNames(list(whole$gs), spec$label))
  }
  # Any other failure is the caller's to report unchanged — a missing file or a
  # faulty gate is not fixed by splitting anything.
  if (!is_channel_inconsistency(conditionMessage(whole$error))) stop(whole$error)

  split <- channel_signature_groups(workspace_channel_signatures(spec))
  if (length(split$groups) < 2L) stop(whole$error)
  cat("\n  ", format_panel_split(spec$label, split), "\n", sep = "")

  parts  <- list()
  failed <- list()
  for (i in seq_along(split$groups)) {
    nm <- sprintf("%s [panel %d]", spec$label, i)
    gs <- tryCatch(
      wsp_to_gatingset(ws, group = spec$group_index, path = root,
                       subset = split$groups[[i]], greedy_match = TRUE),
      error = function(e) {
        failed[[nm]] <<- conditionMessage(e)
        NULL
      })
    if (!is.null(gs)) parts[[nm]] <- gs
  }
  if (!length(parts)) stop(whole$error)
  attr(parts, "failed_parts") <- failed
  parts
}

# -----------------------------------------------------------------------------
# What a build actually produced, as data.
#
# Everything upstream of this describes what we ASKED for. This is the only
# place that knows what ARRIVED, and the two are routinely different — so the
# run's headline counts are taken from here and nowhere else.
#
# Three outcomes are recorded, because each says something different:
#
#   rejected — the workspace listed the sample and it did not arrive. CytoML
#              keys a sample on its name AND its event count, so a file of the
#              right name with the wrong count is a different acquisition and is
#              refused rather than used. Those cells are NOT in the analysis.
#   rescued  — we could not find a file of that NAME and had written the sample
#              off, and CytoML loaded it anyway by content. The cells ARE in the
#              analysis. This is the measure of what content matching is worth,
#              and in "strict" mode it is the measure of how wrong the run's own
#              dropped-sample report was.
#   unlisted — CytoML loaded a sample the group listing never mentioned at all.
#              Rare, and a real signal: fj_ws_get_samples(group_id=) and
#              flowjo_to_gatingset(name=) index groups independently, so this is
#              what a disagreement about which group is "group N" looks like.
#
# `routes` counts how each loaded sample was paired, so "content" mode can show
# what it bought: name / recorded (the file name in the workspace XML) /
# content (CytoML's $FIL + event-count match).
#
# Pure over character vectors — `loaded_raw` is flowWorkspace::sampleNames(gs) —
# so it is testable without the cytometry stack.
# -----------------------------------------------------------------------------
workspace_load_record <- function(spec, loaded_raw) {
  # Sample names come back as "<name>_<$TOT>"; compare on the name part.
  loaded <- sub("_[0-9]+$", "", loaded_raw)
  sought <- spec$files$name
  # `listed` is the group's full sample list; older specs predate it and can
  # only offer the sought set, which makes "rescued" empty rather than wrong.
  listed <- if (is.null(spec$listed)) sought else spec$listed
  via    <- if (is.null(spec$files$via)) rep(NA_character_, length(sought))
            else spec$files$via

  # A sample no name matched can only have been paired by content; a sample we
  # did look up is credited to whichever lookup won.
  route_of <- function(nm) {
    i <- match(nm, sought)
    if (is.na(i) || is.na(via[i])) "content" else via[i]
  }
  routes <- vapply(loaded, route_of, character(1), USE.NAMES = FALSE)

  list(
    label       = spec$label,
    wsp         = spec$wsp,
    group_index = spec$group_index,
    group_name  = spec$group_name,
    match_mode  = if (is.null(spec$match_mode)) "strict" else spec$match_mode,
    n_listed    = length(listed),
    n_sought    = length(sought),
    n_loaded    = length(loaded),
    loaded      = loaded,
    rejected    = setdiff(listed, loaded),
    # A sample no file NAME accounted for, that arrived anyway. Taken from the
    # route rather than from set arithmetic on `sought`: in "content" mode the
    # unresolved samples stay in the spec (with file = NA), so `sought` and
    # `listed` are the same set there and the difference would always be empty.
    # `unlisted` is subtracted so a sample is never reported twice under two
    # different headings — it is the more specific finding of the two.
    rescued     = setdiff(loaded[routes == "content"], setdiff(loaded, listed)),
    unlisted    = setdiff(loaded, listed),
    routes      = c(name     = sum(routes == "name"),
                    recorded = sum(routes == "recorded"),
                    content  = sum(routes == "content"))
  )
}

# -----------------------------------------------------------------------------
# The load records as a block of prose. Returns "" when every workspace loaded
# exactly what it listed, so the caller can stay quiet in the ordinary case.
# -----------------------------------------------------------------------------
format_load_records <- function(records, max_show = 20L) {
  odd <- Filter(function(r)
    length(r$rejected) || length(r$rescued) || length(r$unlisted), records)
  if (!length(odd)) return("")

  show <- function(x) {
    c(paste0("      ", utils::head(x, max_show)),
      if (length(x) > max_show)
        sprintf("      ... and %d more", length(x) - max_show))
  }

  lines <- c("WHAT ACTUALLY LOADED, sample by sample, against what each workspace listed.")
  for (r in odd) {
    lines <- c(lines, sprintf("  %s — listed %d, loaded %d",
                              r$label, r$n_listed, r$n_loaded))
    if (length(r$rescued)) {
      lines <- c(lines, sprintf(
        "    %d sample(s) had NO file of a matching name and loaded anyway, paired by",
        length(r$rescued)),
        "    the $FIL keyword inside the file plus its event count. These cells ARE in",
        "    the analysis:",
        show(r$rescued),
        if (identical(r$match_mode, "strict"))
          "    -> Strict matching reported these as lost. They were not." else NULL)
    }
    if (length(r$rejected)) {
      lines <- c(lines, sprintf(
        "    %d listed sample(s) did NOT arrive and are not in the analysis. CytoML",
        length(r$rejected)),
        "    matches on the sample name AND its event count, so a file of the right",
        "    name with a different count is a different acquisition:",
        show(r$rejected))
    }
    if (length(r$unlisted)) {
      lines <- c(lines, sprintf(
        "    %d sample(s) ARE in the analysis that this workspace's group listing never",
        length(r$unlisted)),
        "    mentioned, so nothing checked them against your FCS folder. The group",
        sprintf("    accessor and the loader disagree about which group is number %d —",
                r$group_index),
        "    naming the group instead of numbering it usually settles it:",
        show(r$unlisted))
    }
  }
  paste(lines, collapse = "\n")
}

# -----------------------------------------------------------------------------
# How the loaded samples were paired, across the whole run. One line, because in
# the ordinary case it is reassurance rather than news — but on a renamed cohort
# it is the difference between "26 samples" and "311".
# -----------------------------------------------------------------------------
format_match_routes <- function(records) {
  if (!length(records)) return("")
  total <- c(name = 0L, recorded = 0L, content = 0L)
  for (r in records) total <- total + r$routes
  if (sum(total) == 0L) return("")
  parts <- c(
    if (total[["name"]])     sprintf("%d by file name", total[["name"]]),
    if (total[["recorded"]]) sprintf("%d by the name the workspace recorded",
                                     total[["recorded"]]),
    if (total[["content"]])  sprintf("%d by content ($FIL + event count)",
                                     total[["content"]]))
  sprintf("Samples paired: %s.", paste(parts, collapse = ", "))
}

# -----------------------------------------------------------------------------
# Report what "content" mode left to CytoML. Returns "" when nothing was.
#
# This is NOT a list of losses — that is the whole point of the mode, and saying
# "42 samples could not be found" about samples that then loaded fine would be
# worse than saying nothing. It is stated as pending work, before the build, and
# the load records above are what settles it afterwards.
# -----------------------------------------------------------------------------
format_deferred_report <- function(deferred, max_show = 10L) {
  if (length(deferred) == 0) return("")
  n <- sum(vapply(deferred, function(d) length(d$samples), integer(1)))
  lines <- c(sprintf(
    paste0("%d sample(s) in %d workspace(s) have no file of a matching NAME in your ",
           "FCS folder."), n, length(deferred)),
    "They are not being dropped: CytoML will try to pair each one by the $FIL keyword",
    "stored inside the files plus the event count, which survives a rename on disk.",
    "What actually loaded is reported once the workspaces are built.", "")
  for (d in deferred) {
    lines <- c(lines, sprintf("  %s — %d of %d matched by name, %d left to content",
                              d$label, d$resolved, d$total, length(d$samples)))
    lines <- c(lines, paste0("      ", utils::head(d$samples, max_show)))
    if (length(d$samples) > max_show) {
      lines <- c(lines, sprintf("      ... and %d more",
                                length(d$samples) - max_show))
    }
  }
  paste(lines, collapse = "\n")
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
