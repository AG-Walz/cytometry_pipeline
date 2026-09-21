#!/usr/bin/env Rscript
# =============================================================================
# analysis/cytokine_combinations.R — Boolean cytokine combinations per sample
# =============================================================================
# One stacked bar per FCS file, one figure per patient. For every cell of a
# parent gate (CD4, CD8) the script asks which of the cytokine gates beneath it
# the cell falls in, and counts the 2^k Boolean combinations ("IFNy+ IL-2+",
# "none", ...). Bars are grouped by timepoint (weeks after the first
# vaccination, in chronological order but on a DISCRETE axis) and faceted by
# stimulation cocktail, so the drift of the T-cell response across vaccinations
# is read left to right within each cocktail column.
#
# The workspace loading reuses src/loading/flowjo_workspace.R and the per-cell
# gate membership reuses get_subpop_membership() from
# src/visualization/umap_from_flowjo.R; everything that does not need CytoML is
# pure base R and covered by test_data/test_cytokine_combinations.R.
#
# Sample names follow the MTB study convention
#   MTB<patient>_<batch>_V<vaccine>W<week>_<cocktail>.fcs
# see scripts/fcs_overview.py for the same pattern in Python.
# =============================================================================

# ---- cytometry src bootstrap (Nextflow-safe); see src/R_shared/paths.R ----
if (!exists(".cyto_paths_loaded", inherits = TRUE)) {
  local({
    root <- Sys.getenv("CYTOMETRY_SRC", "")
    if (!nzchar(root) || !dir.exists(file.path(root, "R_shared"))) {
      of <- NULL
      frames <- sys.frames()
      for (i in rev(seq_along(frames))) {
        f <- frames[[i]]$ofile
        if (!is.null(f) && nzchar(f)) { of <- f; break }
      }
      if (is.null(of)) {
        a <- commandArgs(FALSE); m <- grep("^--file=", a)
        if (length(m)) of <- sub("^--file=", "", a[m[[1]]])
      }
      d <- if (!is.null(of)) dirname(normalizePath(of, mustWork = FALSE)) else getwd()
      while (!dir.exists(file.path(d, "R_shared")) && d != dirname(d)) d <- dirname(d)
      root <- d
    }
    source(file.path(root, "R_shared", "paths.R"), local = FALSE)
  })
}
if (!exists("shared_dir")) shared_dir <- cyto_shared_dir()
source(file.path(shared_dir, "io_utils.R"))    # make_outdir
source(file.path(shared_dir, "plot_utils.R"))  # save_plot_formats

# =============================================================================
# 1. Sample names -> patient / timepoint / cocktail  (pure)
# =============================================================================

MTB_SAMPLE_PATTERN <- paste0(
  "^MTB(?<patient>[0-9]+)_(?<batch>[0-9]+)_V(?<vaccine>[0-9]+)W(?<week>[0-9]+)",
  "_(?<cocktail>.+)$")

# Gate names that mean the same thing but are spelled differently across
# workspaces. Applied to gate BASENAMES only.
GATE_ALIASES <- c(INFy = "IFNy")

gate_alias <- function(x) {
  hit <- x %in% names(GATE_ALIASES)
  x[hit] <- unname(GATE_ALIASES[x[hit]])
  x
}

# The cocktail token of a file name, normalised:
#   "neg", "1_neg"        -> "neg"   (unstimulated control)
#   "194335_018"          -> "194335" (trailing _NNN is a re-acquisition index)
# A leading "1" on a repeated pool ("1215008" for pool 215008 at a later
# timepoint) is left alone on purpose: it cannot be told apart from a genuine
# seven-digit ID by the name alone.
normalise_cocktail <- function(x) {
  x <- as.character(x)
  neg <- !is.na(x) & grepl("neg$", x)
  x[neg] <- "neg"
  x[!neg] <- sub("_[0-9]{3}$", "", x[!neg])
  x
}

# One row per sample name. Patients are zero-padded to three digits so
# "MTB59_1_..." (the FCS files) and "MTB059_1.wsp" (the workspace) agree.
parse_mtb_sample <- function(names) {
  sample <- clean_sample_label(names)
  cap <- annotations_from_pattern(sample, MTB_SAMPLE_PATTERN)
  bad <- is.na(cap$patient)
  if (any(bad)) {
    warning(sum(bad), " sample name(s) do not follow MTB<p>_<b>_V<v>W<w>_<cocktail>",
            " and carry no patient/timepoint: ",
            paste(utils::head(sample[bad], 5), collapse = ", "), call. = FALSE)
  }
  num <- function(v) suppressWarnings(as.integer(v))
  data.frame(
    sample    = sample,
    patient   = ifelse(bad, NA_character_, sprintf("MTB%03d", num(cap$patient))),
    batch     = num(cap$batch),
    vaccine   = num(cap$vaccine),
    week      = num(cap$week),
    timepoint = ifelse(bad, NA_character_,
                       sprintf("V%d W%d", num(cap$vaccine), num(cap$week))),
    cocktail  = normalise_cocktail(cap$cocktail),
    stringsAsFactors = FALSE)
}

# =============================================================================
# 2. Boolean combinations  (pure)
# =============================================================================

# All 2^k combinations of `markers`, in a fixed order: by number of positive
# markers, then by marker order. The order is what keeps colours stable across
# patients and what the legend follows.
combination_levels <- function(markers) {
  k <- length(markers)
  codes <- 0:(2^k - 1)
  labs <- vapply(codes, function(code) {
    pos <- as.logical(bitwAnd(code, 2^(seq_len(k) - 1)))
    if (!any(pos)) "none" else paste0(markers[pos], "+", collapse = " ")
  }, character(1))
  npos <- vapply(codes, function(code) sum(as.logical(bitwAnd(code, 2^(seq_len(k) - 1)))),
                 integer(1))
  labs[order(npos, codes)]
}

# Per-cell combination label from an n x k logical membership matrix.
combination_labels <- function(mem, markers) {
  mem <- matrix(as.logical(mem), ncol = length(markers))
  mem[is.na(mem)] <- FALSE
  k <- length(markers)
  lookup <- vapply(0:(2^k - 1), function(code) {
    pos <- as.logical(bitwAnd(code, 2^(seq_len(k) - 1)))
    if (!any(pos)) "none" else paste0(markers[pos], "+", collapse = " ")
  }, character(1))
  code <- as.integer(mem %*% 2^(seq_len(k) - 1))
  factor(lookup[code + 1L], levels = combination_levels(markers))
}

# Counts per combination, every level present (zero counts kept so a bar always
# closes at 100 %).
combination_counts <- function(mem, markers) {
  lab <- combination_labels(mem, markers)
  tab <- table(lab)
  data.frame(combination = names(tab), count = as.integer(tab),
             stringsAsFactors = FALSE)
}

# percent of the parent gate, per (workspace, sample, parent).
combination_percent <- function(tbl) {
  key <- paste(tbl$workspace, tbl$sample, tbl$parent, sep = "\r")
  tot <- ave(tbl$count, key, FUN = sum)
  tbl$percent <- ifelse(tot > 0, 100 * tbl$count / tot, 0)
  tbl
}

# Cocktail minus the unstimulated control of the same patient x timepoint x
# parent, floored at 0. "none" is recomputed so the stack still closes at 100.
# Samples without a `neg` at their timepoint cannot be corrected and are
# dropped with a warning.
subtract_background <- function(tbl) {
  key <- function(d) paste(d$patient, d$timepoint, d$parent, sep = "\r")
  neg <- tbl[tbl$cocktail == "neg", , drop = FALSE]
  base <- stats::setNames(neg$percent, paste(key(neg), neg$combination, sep = "\r"))
  out <- tbl[tbl$cocktail != "neg", , drop = FALSE]
  has_neg <- key(out) %in% key(neg)
  if (any(!has_neg)) {
    warning("no `neg` control for ", length(unique(out$sample[!has_neg])),
            " sample(s); they are left out of the background-subtracted plot: ",
            paste(utils::head(unique(out$sample[!has_neg]), 5), collapse = ", "),
            call. = FALSE)
    out <- out[has_neg, , drop = FALSE]
  }
  if (!nrow(out)) return(out)
  b <- base[paste(key(out), out$combination, sep = "\r")]
  b[is.na(b)] <- 0
  pos <- out$combination != "none"
  out$percent[pos] <- pmax(out$percent[pos] - b[pos], 0)
  k <- paste(out$workspace, out$sample, out$parent, sep = "\r")
  rest <- ave(ifelse(pos, out$percent, 0), k, FUN = sum)
  out$percent[!pos] <- pmax(100 - rest[!pos], 0)
  out$count <- NA_integer_
  out
}

# =============================================================================
# 3. Gate paths  (pure over the character vector gs_get_pop_paths() returns)
# =============================================================================

# Full path of `parent` and of its direct children named like `markers`
# (after aliasing). A marker without a gate is NA and reported once.
resolve_gate_paths <- function(pop_paths, parent, markers) {
  base <- sub(".*/", "", pop_paths)
  hit <- which(base == parent)
  if (length(hit) != 1L) {
    stop("gate '", parent, "' matches ", length(hit), " population(s) in the ",
         "gating tree; expected exactly one. Tree: ",
         paste(pop_paths, collapse = ", "), call. = FALSE)
  }
  parent_path <- pop_paths[hit]
  prefix <- paste0(parent_path, "/")
  kids <- pop_paths[startsWith(pop_paths, prefix)]
  kids <- kids[!grepl("/", substring(kids, nchar(prefix) + 1L))]
  kid_names <- gate_alias(sub(".*/", "", kids))
  marker_paths <- kids[match(markers, kid_names)]
  names(marker_paths) <- markers
  miss <- markers[is.na(marker_paths)]
  if (length(miss)) {
    warning("under ", parent_path, " no gate is named ",
            paste(miss, collapse = ", "), "; those markers count as negative.",
            call. = FALSE)
  }
  list(parent = parent_path, markers = marker_paths)
}

# =============================================================================
# 4. CytoML-dependent part (analysis hosts only)
# =============================================================================

# Discover workspaces, pair them with FCS files by name and build one GatingSet
# per workspace — the same chain umap_from_flowjo() uses. Returns a named list
# of GatingSets (a workspace acquired with two channel sets yields two parts).
load_workspaces <- function(wsp, fcs_dir, group = "All Samples") {
  wsp_paths <- find_workspace_files(wsp)
  wsp_names <- workspace_labels(wsp_paths)
  cat("Workspaces found (", length(wsp_paths), "): ",
      paste(wsp_names, collapse = ", "), "\n", sep = "")
  fcs_index <- index_fcs_dir(fcs_dir)
  cat("FCS files indexed: ", length(fcs_index), " under ", fcs_dir, "\n", sep = "")

  resolved <- resolve_workspace_specs(wsp_paths, wsp_names, fcs_index, group)
  report <- format_resolution_report(resolved$problems, fcs_root = fcs_dir)
  if (nzchar(report)) {
    cat(report, "\n", sep = "", file = stderr())
    stop("could not build the workspace set: ", length(resolved$problems),
         " of ", length(wsp_paths), " workspace(s) unusable. See the report above.",
         call. = FALSE)
  }
  gs_list <- list()
  for (label in names(resolved$specs)) {
    spec <- resolved$specs[[label]]
    cat(sprintf("Building GatingSet for '%s' — %s, %d sample(s)\n",
                label, describe_group(spec), length(spec$listed)))
    parts <- build_workspace_gatingsets(spec)
    for (part in names(parts)) gs_list[[part]] <- parts[[part]]
    failed <- attr(parts, "failed_parts")
    for (part in names(failed)) {
      warning("part '", part, "' of workspace '", label, "' did not load: ",
              failed[[part]], call. = FALSE)
    }
  }
  gs_list
}

# One sample, one parent gate: counts of every cytokine combination among the
# parent's cells. `marker_paths` may contain NA (gate absent -> all negative).
sample_combination_table <- function(gh, parent_path, marker_paths, markers) {
  viable <- flowWorkspace::gh_pop_get_indices(gh, parent_path)
  mem <- get_subpop_membership(gh, viable, unname(marker_paths))
  out <- combination_counts(mem, markers)
  out$parent_count <- sum(viable)
  out
}

# The long table behind every plot: one row per workspace x sample x parent x
# combination, with the sample's patient/timepoint/cocktail alongside.
cytokine_combination_table <- function(gs_list, parents = c("CD4", "CD8"),
                                       markers = c("CD107a", "IFNy", "IL-2", "TNF")) {
  rows <- list()
  n_total <- sum(vapply(gs_list, function(gs) length(flowWorkspace::sampleNames(gs)),
                        integer(1)))
  done <- 0L
  for (w in names(gs_list)) {
    gs <- gs_list[[w]]
    pops <- flowWorkspace::gs_get_pop_paths(gs, path = "full")
    paths <- lapply(parents, function(p) resolve_gate_paths(pops, p, markers))
    names(paths) <- parents
    for (s in flowWorkspace::sampleNames(gs)) {
      gh <- gs[[s]]
      for (p in parents) {
        one <- sample_combination_table(gh, paths[[p]]$parent, paths[[p]]$markers,
                                        markers)
        one$workspace <- w
        one$sample <- clean_sample_label(s)
        one$parent <- p
        rows[[length(rows) + 1L]] <- one
      }
      done <- done + 1L
      emit_progress("extract", done, n_total, s)
    }
  }
  tbl <- do.call(rbind, rows)
  meta <- parse_mtb_sample(tbl$sample)
  tbl <- cbind(tbl[c("workspace", "sample")],
               meta[c("patient", "batch", "vaccine", "week", "timepoint", "cocktail")],
               tbl[c("parent", "parent_count", "combination", "count")])
  combination_percent(tbl)
}

# =============================================================================
# 5. Plots and table output
# =============================================================================

# Colours: "none" is a neutral grey; the positive combinations are grouped by
# how many cytokines are positive (one hue family per count), so
# polyfunctional cells stand out and the 16 legend entries read as four groups.
combination_palette <- function(levels) {
  npos <- ifelse(levels == "none", 0L, lengths(strsplit(levels, " ")))
  families <- list(
    grDevices::colorRampPalette(c("#9ecae1", "#08519c")),   # 1 positive: blues
    grDevices::colorRampPalette(c("#a1d99b", "#006d2c")),   # 2: greens
    grDevices::colorRampPalette(c("#fdae6b", "#a63603")),   # 3: oranges
    grDevices::colorRampPalette(c("#e7298a", "#67001f")))   # 4+: magenta/red
  cols <- rep("grey85", length(levels))
  for (n in sort(unique(npos[npos > 0L]))) {
    idx <- which(npos == n)
    fam <- families[[min(n, length(families))]]
    cols[idx] <- fam(length(idx))
  }
  stats::setNames(cols, levels)
}

# One patient: x = timepoint (chronological, discrete), stack = combinations,
# facet rows = parent gate, facet columns = cocktail (`neg` first).
#
# `drop_none` leaves the all-negative cells out and frees the y axis per row:
# the cytokine-positive cells are a few percent of the parent, and on a 0-100
# scale they are a sliver nobody can read. That figure shows the drift; the
# full stack shows how small the response is.
plot_patient_combinations <- function(tbl, patient, output_dir,
                                      formats = list(png = TRUE),
                                      subtract = FALSE, drop_none = FALSE,
                                      png_dpi = 200) {
  require(ggplot2)
  d <- tbl[!is.na(tbl$patient) & tbl$patient == patient, , drop = FALSE]
  if (subtract) d <- subtract_background(d)
  if (drop_none) d <- d[d$combination != "none", , drop = FALSE]
  if (!nrow(d)) {
    warning("nothing to plot for ", patient, call. = FALSE)
    return(invisible(NULL))
  }
  tp <- unique(d[c("timepoint", "week", "vaccine")])
  tp <- tp[order(tp$week, tp$vaccine), ]
  d$timepoint <- factor(d$timepoint, levels = tp$timepoint)
  cocktails <- sort(unique(d$cocktail))
  cocktails <- c(intersect("neg", cocktails), setdiff(cocktails, "neg"))
  d$cocktail <- factor(d$cocktail, levels = cocktails)
  lv <- combination_levels_from(d$combination)
  d$combination <- factor(d$combination, levels = rev(lv))  # "none" at the bottom

  p <- ggplot(d, aes(x = timepoint, y = percent, fill = combination)) +
    geom_col(position = "stack", width = 0.8, colour = "white", linewidth = 0.2) +
    facet_grid(parent ~ cocktail, space = "free_x",
               scales = if (drop_none) "free" else "free_x") +
    scale_fill_manual(values = combination_palette(lv), breaks = lv,
                      name = "Cytokine combination") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste0(patient, ": cytokine combinations per timepoint",
                     if (subtract) ", cocktail minus neg" else "",
                     if (drop_none) ", cytokine-positive cells only" else ""),
      subtitle = "Columns: stimulation cocktail; rows: parent gate; x: vaccination and week after first vaccination",
      x = "Timepoint", y = "% of parent gate") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")

  base <- file.path(output_dir, paste0(patient, "_cytokine_combinations",
                                       if (subtract) "_minus_neg" else "",
                                       if (drop_none) "_positive_only" else ""))
  width <- max(8, 2.2 * length(cocktails) + 3.5)
  height <- 3.2 * length(unique(d$parent)) + 1.5
  saved <- save_plot_formats(p, base, width = width, height = height,
                             formats = formats, png_dpi = png_dpi)
  cat("  wrote ", base, " [", paste(saved, collapse = ", "), "]\n", sep = "")
  invisible(p)
}

# The level order for the labels actually present, recovered from the labels
# themselves (so a table read back from CSV plots the same way).
combination_levels_from <- function(labels) {
  labels <- unique(as.character(labels))
  npos <- ifelse(labels == "none", 0L, lengths(strsplit(labels, " ")))
  labels[order(npos, labels)]
}

write_combination_table <- function(tbl, output_dir) {
  path <- file.path(output_dir, "cytokine_combinations.csv")
  utils::write.table(tbl, path, sep = ";", row.names = FALSE, quote = FALSE)
  cat("  wrote ", path, "\n", sep = "")
  invisible(path)
}

read_combination_table <- function(path) {
  utils::read.table(path, sep = ";", header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, colClasses = c(cocktail = "character"))
}

# =============================================================================
# 6. Entry point
# =============================================================================

cytokine_combinations <- function(wsp, fcs_dir, output_dir,
                                  patients = NULL,
                                  parents = c("CD4", "CD8"),
                                  markers = c("CD107a", "IFNy", "IL-2", "TNF"),
                                  group = "All Samples",
                                  output_formats = list(png = TRUE),
                                  from_table = NULL) {
  make_outdir(output_dir)
  if (is.null(from_table)) {
    gs_list <- load_workspaces(wsp, fcs_dir, group)
    tbl <- cytokine_combination_table(gs_list, parents, markers)
    write_combination_table(tbl, output_dir)
  } else {
    cat("Re-plotting from ", from_table, "\n", sep = "")
    tbl <- read_combination_table(from_table)
  }
  found <- sort(unique(tbl$patient[!is.na(tbl$patient)]))
  if (!is.null(patients)) {
    miss <- setdiff(patients, found)
    if (length(miss)) {
      warning("no samples for patient(s) ", paste(miss, collapse = ", "),
              "; available: ", paste(found, collapse = ", "), call. = FALSE)
    }
    found <- intersect(patients, found)
  }
  cat("Plotting ", length(found), " patient(s): ", paste(found, collapse = ", "),
      "\n", sep = "")
  for (i in seq_along(found)) {
    plot_patient_combinations(tbl, found[i], output_dir, output_formats)
    plot_patient_combinations(tbl, found[i], output_dir, output_formats,
                              subtract = TRUE)
    plot_patient_combinations(tbl, found[i], output_dir, output_formats,
                              subtract = TRUE, drop_none = TRUE)
    emit_progress("plot", i, length(found), found[i])
  }
  invisible(tbl)
}
