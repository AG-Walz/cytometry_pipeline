#!/usr/bin/env Rscript
# =============================================================================
# visualization/umap_load_report.R — the workspace load, as data
# =============================================================================
# One JSON document describing everything a caller needs to know about a
# --list_populations load: which workspaces resolved, how each sample was
# paired, what each workspace's gating tree and marker panel look like, what
# diverges between them, and what could not be used.
#
# WHY THIS EXISTS
#
# The Dash front-end used to recover all of that by regular expression from the
# run's stdout — four block headers, three indentation levels, and a comment on
# both sides begging the next person to keep them in step. That is a parser for
# a format that was never designed to be one, and it broke the way such things
# break: a population list that came back empty because the first line was
# "root" and the pattern required a leading slash.
#
# A 33-workspace load also produces an 862-line report. Prose is the right shape
# for one workspace and the wrong shape for thirty; the GUI needs a table, and a
# table needs data.
#
# So R — which is the only side that can compute any of this — writes it out as
# data, and the front-end reads it. The text report stays exactly as it was: it
# is what someone running the CLI reads, and it is not going away.
#
# CONTRACT
#
# `schema` is versioned. A consumer should check the major version and refuse a
# document it does not understand rather than silently reading a field that has
# moved. Every list-valued field is always present, empty rather than absent, so
# a reader never has to distinguish "none" from "missing".
#
# Pure except for the final write: build_load_report() takes values and returns
# a list, so it is testable without CytoML, a workspace, or a file system.
# =============================================================================

#: Bumped on any change to the document's shape. Minor = additions only.
UMAP_REPORT_SCHEMA <- "1.0"

# -----------------------------------------------------------------------------
# Assemble the report.
#
# Everything here has already been computed by the caller — this arranges it,
# and deliberately performs no analysis of its own, so the numbers in the JSON
# and the numbers in the log can never disagree.
#
# @param inputs list(wsp, fcs_dir, group, match_mode, allow_missing)
# @param workspaces_found Character vector of every workspace label discovered.
# @param load_records Named list from workspace_load_record(), one per loaded
#   workspace.
# @param specs The resolved specs (for the workspace file path and group).
# @param problems,dropped,deferred The three failure/deferral maps from
#   resolve_workspace_specs().
# @param build_failed Named list of label -> CytoML's own error message.
# @param pops_by_wsp Named list of each workspace's full gate paths.
# @param shared_pops The gates every workspace has.
# @param markers_by_wsp Named list of each workspace's marker panel.
# @param common_markers The markers every workspace has.
# @param panel_cmp compare_marker_sets() over markers_by_wsp, or NULL.
# @param suggestions suggest_marker_aliases() output, or NULL.
# @param channel_dict channel_marker_dictionary() output, or NULL.
# @param samples data.frame with workspace/raw/display, or NULL.
# @param annotations Named list of column -> levels, or NULL.
# @param unnamed_samples How many samples carried no $PnS at all.
# -----------------------------------------------------------------------------
build_load_report <- function(inputs,
                              workspaces_found = character(0),
                              load_records = list(),
                              specs = list(),
                              problems = list(), dropped = list(),
                              deferred = list(), build_failed = list(),
                              pops_by_wsp = list(), shared_pops = character(0),
                              markers_by_wsp = list(),
                              common_markers = character(0),
                              panel_cmp = NULL, suggestions = NULL,
                              channel_dict = NULL, samples = NULL,
                              annotations = NULL, unnamed_samples = 0L) {

  # A gate is reported as {path, label, depth}: the path is what --population
  # needs, the label is what a human can read (see gate_display_label — a gate
  # named " " is invisible otherwise), and the depth saves every consumer from
  # re-deriving the nesting by counting slashes.
  gate_entries <- function(paths) {
    if (!length(paths)) return(list())
    lapply(seq_along(paths), function(i) {
      p <- paths[i]
      list(path = p,
           label = gate_display_label(p),
           depth = if (identical(p, "root")) 0L
                   else length(strsplit(sub("^/", "", p), "/", fixed = TRUE)[[1]]))
    })
  }

  # One row per workspace DISCOVERED — loaded, failed and dropped alike. A
  # front-end that only ever hears about successes cannot show you why the other
  # thirty-one are not in your analysis, which is the question actually being
  # asked of it.
  ws_rows <- lapply(workspaces_found, function(label) {
    rec  <- load_records[[label]]
    spec <- specs[[label]]
    prob <- problems[[label]]
    drop <- dropped[[label]]
    defr <- deferred[[label]]

    status <- if (!is.null(rec)) "loaded"
              else if (!is.null(build_failed[[label]])) "build_failed"
              else if (!is.null(prob)) "unusable"
              else if (!is.null(drop) && identical(drop$kept, 0L)) "no_files"
              else "not_loaded"

    reason <- if (!is.null(build_failed[[label]])) build_failed[[label]]
              else if (!is.null(prob$error)) prob$error
              else if (identical(status, "no_files"))
                "no FCS file of any name this workspace asks for was found"
              else NA_character_

    list(
      label      = label,
      wsp        = if (!is.null(spec)) spec$wsp else
                   if (!is.null(prob)) prob$wsp else
                   if (!is.null(drop)) drop$wsp else NA_character_,
      status     = status,
      reason     = reason,
      group      = if (!is.null(spec)) describe_group(spec) else NA_character_,
      n_listed   = if (!is.null(rec)) rec$n_listed else
                   if (!is.null(drop)) drop$total else
                   if (!is.null(prob)) prob$total else NA_integer_,
      n_loaded   = if (!is.null(rec)) rec$n_loaded else 0L,
      # How the loaded samples were paired. The whole point of content matching
      # is visible here and nowhere else.
      routes     = if (!is.null(rec)) as.list(rec$routes) else
                   list(name = 0L, recorded = 0L, content = 0L),
      rescued    = if (!is.null(rec)) as.character(rec$rescued) else character(0),
      rejected   = if (!is.null(rec)) as.character(rec$rejected) else character(0),
      unlisted   = if (!is.null(rec)) as.character(rec$unlisted) else character(0),
      deferred   = if (!is.null(defr)) as.character(defr$samples) else character(0),
      gates      = gate_entries(pops_by_wsp[[label]]),
      n_gates    = length(pops_by_wsp[[label]]),
      markers    = as.character(markers_by_wsp[[label]]),
      # The gates and markers this workspace does NOT share with the others.
      # Precomputed because it is what decides whether a joint run is possible,
      # and every consumer would otherwise compute it the same way.
      gates_not_shared   = setdiff(pops_by_wsp[[label]], shared_pops),
      markers_not_shared = setdiff(markers_by_wsp[[label]], common_markers)
    )
  })
  names(ws_rows) <- NULL

  routes_total <- list(name = 0L, recorded = 0L, content = 0L)
  for (r in load_records) {
    for (k in names(routes_total)) {
      routes_total[[k]] <- routes_total[[k]] + as.integer(r$routes[[k]])
    }
  }

  not_common <- if (!is.null(panel_cmp) && length(panel_cmp$not_common)) {
    lapply(panel_cmp$not_common, function(m) list(
      marker = m,
      present_in = colnames(panel_cmp$presence)[panel_cmp$presence[m, ]]))
  } else list()

  alias_rows <- if (!is.null(suggestions) && nrow(suggestions)) {
    lapply(seq_len(nrow(suggestions)), function(i) list(
      names      = strsplit(suggestions$names[i], ", ", fixed = TRUE)[[1]],
      suggested  = suggestions$suggested[i],
      confidence = suggestions$confidence[i]))
  } else list()

  conflict_rows <- if (!is.null(channel_dict) && nrow(channel_dict$conflicts)) {
    lapply(seq_len(nrow(channel_dict$conflicts)), function(i) list(
      channel = channel_dict$conflicts$channel[i],
      names   = strsplit(channel_dict$conflicts$names[i], ", ", fixed = TRUE)[[1]],
      chosen  = channel_dict$conflicts$chosen[i],
      counts  = channel_dict$conflicts$counts[i]))
  } else list()

  sample_rows <- if (!is.null(samples) && nrow(samples)) {
    lapply(seq_len(nrow(samples)), function(i) list(
      display = samples$display[i], workspace = samples$workspace[i],
      raw = samples$raw[i]))
  } else list()

  list(
    schema  = UMAP_REPORT_SCHEMA,
    inputs  = inputs,
    totals  = list(
      workspaces_found  = length(workspaces_found),
      workspaces_loaded = length(load_records),
      samples           = sum(vapply(load_records,
                                     function(r) r$n_loaded, integer(1))),
      samples_listed    = sum(vapply(load_records,
                                     function(r) r$n_listed, integer(1))),
      routes            = routes_total,
      samples_without_marker_names = as.integer(unnamed_samples)
    ),
    workspaces  = ws_rows,
    populations = gate_entries(shared_pops),
    markers = list(
      common      = as.character(common_markers),
      not_common  = not_common,
      suggestions = alias_rows,
      dictionary  = if (!is.null(channel_dict) && length(channel_dict$dictionary))
                      as.list(channel_dict$dictionary) else
                      structure(list(), names = character(0)),
      conflicts   = conflict_rows
    ),
    samples     = sample_rows,
    annotations = if (is.null(annotations)) structure(list(), names = character(0))
                  else annotations
  )
}

# -----------------------------------------------------------------------------
# Write the report. Never fatal: a run that produced its figures must not fail
# because a side-car file could not be written, so a failure is a warning and
# the caller carries on.
#
# auto_unbox = TRUE renders a length-1 vector as a scalar rather than a
# one-element array, which is what makes the document natural to read on the
# other side. The always-a-list fields above are built as lists precisely so
# that unboxing cannot turn a one-element collection into a bare value.
# -----------------------------------------------------------------------------
write_load_report <- function(report, path) {
  if (is.null(path) || !nzchar(path)) return(invisible(NULL))
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("jsonlite is not installed; no JSON report was written to ", path,
            call. = FALSE)
    return(invisible(NULL))
  }
  ok <- tryCatch({
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(report, path, auto_unbox = TRUE, null = "null",
                         na = "null", pretty = TRUE)
    TRUE
  }, error = function(e) {
    warning("could not write the JSON report to ", path, ": ",
            conditionMessage(e), call. = FALSE)
    FALSE
  })
  if (ok) cat("\nLoad report: ", path, "\n", sep = "")
  invisible(ok)
}

# -----------------------------------------------------------------------------
# What --sample_metadata / --sample_pattern would resolve to, for a given list
# of sample names. Written as JSON; no workspace, no FCS file, no GatingSet.
#
# This exists so a form can populate its annotation-column picker and show a
# match rate the moment a file is chosen. That used to be a Python
# reimplementation of the matching rules living in the Dash app, whose own
# docstring conceded it was "a mirror" and that "if the two ever disagree, R
# wins and this is the bug". There is now nothing to disagree with: the preview
# calls resolve_sample_annotations(), which is what a run calls.
#
# Errors become an "error" field rather than a stop(): the caller is a form
# reacting to half-typed input, and a bad path is an ordinary state there.
# -----------------------------------------------------------------------------
preview_annotations <- function(out_path, samples = NULL, metadata = NULL,
                                pattern = NULL, id_column = NULL) {
  fail <- function(msg) list(columns = list(), levels = structure(list(), names = character(0)),
                             total = 0L, matched = 0L, unmatched = list(),
                             source = NA_character_, error = msg)

  result <- tryCatch({
    names_v <- if (!is.null(samples) && nzchar(samples) && file.exists(samples))
      readLines(samples, warn = FALSE) else character(0)
    names_v <- names_v[nzchar(trimws(names_v))]

    # resolve_sample_annotations() keys on a sample's several identities; a
    # preview only has the display name, which is the one a GUI shows.
    idx <- data.frame(workspace = rep("", length(names_v)),
                      raw = names_v, clean = names_v, display = names_v,
                      stringsAsFactors = FALSE)

    md  <- if (!is.null(metadata) && nzchar(trimws(metadata))) metadata else NULL
    pat <- if (!is.null(pattern) && nzchar(trimws(pattern))) pattern else NULL
    if (is.null(md) && is.null(pat)) {
      list(columns = list(), levels = structure(list(), names = character(0)),
           total = length(names_v), matched = 0L, unmatched = list(),
           source = NA_character_, error = NULL)
    } else {
      ann <- suppressWarnings(resolve_sample_annotations(
        idx, metadata_file = md, pattern = pat, id_column = id_column))
      cols <- colnames(ann)
      # A sample is "matched" when at least one column resolved to something
      # other than the unassigned placeholder.
      unmatched <- if (length(cols)) {
        hit <- Reduce(`|`, lapply(cols, function(cl)
          !is.na(ann[[cl]]) & ann[[cl]] != .ANNOT_UNASSIGNED))
        names_v[!hit]
      } else character(0)
      list(
        columns = as.list(cols),
        levels  = stats::setNames(
          lapply(cols, function(cl) as.list(unique(ann[[cl]]))), cols),
        total = length(names_v),
        matched = length(names_v) - length(unmatched),
        unmatched = as.list(unmatched),
        source = if (!is.null(md)) basename(md) else "the sample name",
        error = NULL)
    }
  }, error = function(e) fail(conditionMessage(e)))

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(result, out_path, auto_unbox = TRUE, null = "null",
                         na = "null")
  }
  invisible(result)
}
