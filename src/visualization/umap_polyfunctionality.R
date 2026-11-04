# =============================================================================
# visualization/umap_polyfunctionality.R — polyfunctionality maps of a
# gate-combination run (umap_from_flowjo(), --layout gates or --subpop_mode
# combinations)
# =============================================================================
# Every cell is coloured by how many of the functional gates it is in, read from
# its combination label: "CD4+ | IFNy+ TNFa+" -> 2, "CD4+ | none" -> 0. Cells
# with more positive gates are drawn larger and on top, in the five blues of
# the MTB cohort maps (scripts/plot_polyfunctionality.R --size_by_polyfunctionality).
# Cells whose label has no functional tier ("CD4- CD8-", "CD4+ CD8+", "other
# (rare)") have no count and are left out.
#
# One map of the whole embedding, and one per level of every --split_by column,
# all on the frame of the whole embedding so that they line up; and a table of
# the counts per sample.
#
# Needs umap_style.R (theme_umap_clean, umap_axis_arrows, umap_point_layer,
# auto_point_size) and R_shared/plot_utils.R (save_plot_formats,
# create_safe_filename) sourced first; uses run_subdir() of umap_from_flowjo.R,
# ordered_levels() of umap_gate_layout.R and parallel_lapply() of
# R_shared/memory_utils.R when they are there.
# =============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#: The MTB cohort colours for 0 to 4 positive gates (light to dark blue).
.PF_COLOURS <- c("#DCE9F6", "#9CC3E5", "#4A8FC9", "#0B4F8A", "#062F52")
#: Labels that are not a combination of gate calls.
.PF_NOT_COUNTED <- c("other (rare)", "(no decision)", "(none)", "(multiple)")

# Positive markers of one label part ("IFNy+ TNFa+" -> c("IFNy", "TNFa");
# "none" and "IFNy- TNFa-" -> none). A call ends at its sign, so the part is
# split after each "+" or "-" followed by space, not at every space: marker
# names can contain spaces ("TNF a+ IFN y+" -> c("TNF a", "IFN y")). Split at
# every space, the J445 run counted gates "a" and "y" besides "TNF a" and
# "IFN y", and drew a 0-7 scale for five gates.
.pf_positive <- function(part) {
  tok <- strsplit(trimws(part), "(?<=[+-])[[:space:]]+", perl = TRUE)[[1]]
  sub("\\+$", "", tok[endsWith(tok, "+")])
}

# -----------------------------------------------------------------------------
# Number of positive functional gates per cell, from its combination label.
#
# The labels are tiers joined by " | " (combination_labels() in
# umap_gate_layout.R). The functional gates are the last tier, counted on the
# labels that have every tier; a label with fewer (a lineage without gates
# below it, e.g. "CD4- CD8-") is NA, as are .PF_NOT_COUNTED. A population with a
# single tier ("IFNy+ TNFa+", no " | " anywhere) counts the whole label.
#
# Returns list(degree = integer per label, markers = the gates counted, in the
# order the labels spell them), or NULL when no label has a positive gate.
# -----------------------------------------------------------------------------
polyfunctional_degree <- function(labels) {
  lab <- as.character(labels)
  u <- unique(lab[!is.na(lab) & !lab %in% .PF_NOT_COUNTED])
  if (!length(u)) return(NULL)
  parts <- strsplit(u, " | ", fixed = TRUE)
  n_parts <- lengths(parts)
  full <- n_parts == max(n_parts)
  last <- vapply(parts, function(p) p[length(p)], character(1))
  pos <- lapply(last, .pf_positive)
  pos[!full] <- list(NULL)
  # Marker order: as the most polyfunctional label spells them, then the rest.
  by_n <- order(-lengths(pos))
  markers <- unique(unlist(pos[by_n], use.names = FALSE))
  if (!length(markers)) return(NULL)
  deg_u <- ifelse(full, lengths(pos), NA_integer_)
  list(degree = as.integer(deg_u[match(lab, u)]), markers = markers)
}

# The functional gates of a run: the markers of the last tier of its gate
# decisions (marker_tiers() of umap_gate_layout.R gives marker -> tier).
polyfunctional_gates <- function(tiers) {
  if (!length(tiers)) return(NULL)
  names(tiers)[tiers == max(tiers)]
}

# Colours for 0 to k positive gates: the MTB blues for k = 4, the same ramp
# stretched or squeezed otherwise.
polyfunctionality_colours <- function(k) {
  if (k == length(.PF_COLOURS) - 1L) return(.PF_COLOURS)
  grDevices::colorRampPalette(.PF_COLOURS)(k + 1L)
}

# Point sizes for 0 to k positive gates: `base` times 1 to 3, evenly spaced
# (1, 1.5, 2, 2.5, 3 for four gates, as in the MTB maps).
polyfunctionality_sizes <- function(k, base) {
  base * (1 + 2 * (0:k) / max(k, 1L))
}

# Per sample: cells (weighted, as the gates layout samples them) and percent
# of the counted cells with 0 to k positive gates.
polyfunctionality_table <- function(cells, degree, k, annot_cols = character(0)) {
  w <- if (is.null(cells$weight)) rep(1, nrow(cells)) else cells$weight
  keep <- !is.na(degree)
  s <- as.character(cells$sample)[keep]
  d <- factor(degree[keep], levels = 0:k)
  tot <- tapply(w[keep], s, sum)
  tab <- tapply(w[keep], list(s, d), sum)
  tab[is.na(tab)] <- 0
  out <- data.frame(
    sample = rep(rownames(tab), times = k + 1L),
    n_positive = rep(0:k, each = nrow(tab)),
    cells = round(as.vector(tab), 2),
    stringsAsFactors = FALSE)
  out$percent <- round(100 * out$cells / as.numeric(tot[out$sample]), 4)
  first <- match(out$sample, as.character(cells$sample))
  for (cl in annot_cols) {
    out[[cl]] <- as.character(cells[[paste0("annot_", cl)]])[first]
  }
  out[order(out$sample, out$n_positive), , drop = FALSE]
}

# -----------------------------------------------------------------------------
# Entry point. `cells`: one row per plotted cell with UMAP1, UMAP2, sample,
# combination (the gate-combination label), optionally weight and annot_<col>
# columns, as plot_gate_layout() takes it. `gates`: the functional gates the
# run decided (its last tier, polyfunctional_gates()); without it they are read
# off the labels, which miss a gate no kept combination is positive for, and
# the scale then stops short (0 to 3 when IL-2+ cells were all pooled as rare).
# Returns (invisibly) the plot index
# entries (file, kind, column, level, reference), or NULL when the labels have
# no functional gates.
# -----------------------------------------------------------------------------
plot_polyfunctionality <- function(cells, population, output_dir,
                                   split_cols = character(0),
                                   level_order = list(),
                                   point_size = NULL, point_alpha = NULL,
                                   plot_width = 10, plot_height = 8,
                                   output_formats = list(pdf = FALSE, png = TRUE,
                                                         svg = FALSE),
                                   n_cores = 1L, gates = NULL) {
  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  pf <- polyfunctional_degree(cells$combination)
  if (!is.null(pf) && length(gates)) {
    pf$markers <- c(gates, setdiff(pf$markers, gates))
  }
  if (is.null(pf)) {
    cat("Polyfunctionality maps skipped: no combination label has a positive ",
        "gate below the lineage.\n", sep = "")
    return(invisible(NULL))
  }
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rng_state <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", rng_state, envir = globalenv()), add = TRUE)
  }
  deg <- pf$degree
  k <- length(pf$markers)
  w <- if (is.null(cells$weight)) rep(1, nrow(cells)) else cells$weight
  cols <- polyfunctionality_colours(k)
  pop_tag  <- create_safe_filename(sub("^/", "", population))
  pop_name <- basename(sub("/$", "", population))
  if (!nzchar(pop_name)) pop_name <- population
  subdir <- function(column = NULL) {
    if (!exists("run_subdir", mode = "function")) return("")
    if (is.null(column)) run_subdir("poly") else run_subdir("poly_split", column = column)
  }
  lims_x <- range(cells$UMAP1, na.rm = TRUE)
  lims_y <- range(cells$UMAP2, na.rm = TRUE)
  legend_title <- paste0("Positive gates per cell\n(",
                         paste(pf$markers, collapse = ", "), ")")
  alpha <- point_alpha %||% 0.85

  draw <- function(rows, title, subdir_rel, stem) {
    counted <- rows[!is.na(deg[rows])]
    d_c <- deg[counted]
    tab <- vapply(0:k, function(i) sum(w[counted][d_c == i]), 0)
    n_t <- sum(tab)
    pct <- if (n_t > 0) sprintf("%.2f %%", 100 * tab / n_t) else rep("-", k + 1L)
    leg <- sprintf("%d  (%s, %s)", 0:k,
                   format(round(tab), big.mark = ",", trim = TRUE), pct)
    # MTB look: points no larger than 0.4 at the least polyfunctional, smaller
    # for a dense map (auto_point_size), so the largest stay dots, not blots.
    base <- point_size %||% min(0.4, auto_point_size(length(counted)))
    sizes <- polyfunctionality_sizes(k, base)
    df <- data.frame(UMAP1 = cells$UMAP1[counted], UMAP2 = cells$UMAP2[counted],
                     lvl = factor(d_c, levels = 0:k, labels = leg),
                     sz = sizes[d_c + 1L])
    df <- df[order(d_c), , drop = FALSE]  # most polyfunctional drawn last, on top
    leg_size <- 2 + 1.5 * (sizes - min(sizes)) / max(diff(range(sizes)), 1e-9)
    n_other <- length(rows) - length(counted)
    # One invisible point per count, so a count no cell has keeps its key.
    keys <- data.frame(UMAP1 = lims_x[1], UMAP2 = lims_y[1],
                       lvl = factor(leg, levels = leg), sz = sizes)
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = UMAP1, y = UMAP2, colour = lvl, size = sz),
                          data = keys, alpha = 0) +
      umap_point_layer(ggplot2::aes(x = UMAP1, y = UMAP2, colour = lvl, size = sz),
                       data = df, alpha = alpha) +
      ggplot2::scale_size_identity() +
      ggplot2::scale_colour_manual(
        values = stats::setNames(cols, leg), drop = FALSE, name = legend_title,
        guide = ggplot2::guide_legend(override.aes = list(size = leg_size, alpha = 1))) +
      umap_axis_arrows(lims_x, lims_y) +
      ggplot2::coord_fixed(xlim = lims_x, ylim = lims_y) +
      theme_umap_clean(11) +
      ggplot2::labs(
        title = title,
        subtitle = sprintf("%s cells from %d sample(s)%s",
                           format(round(n_t), big.mark = ","),
                           length(unique(cells$sample[rows])),
                           if (n_other > 0) sprintf(
                             "; %s without gates below the lineage not shown",
                             format(n_other, big.mark = ",")) else ""),
        caption = paste0("Same embedding as the whole run. Percentages: of the ",
                         "cells shown. Cells with more positive gates are drawn ",
                         "on top and larger.")) +
      ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0),
                     plot.caption.position = "plot")
    d <- file.path(output_dir, subdir_rel)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    save_plot_formats(p, file.path(d, stem), width = plot_width,
                      height = plot_height, formats = output_formats)
    stem
  }

  jobs <- list()
  index <- list()
  add_job <- function(rows, title, column, level, stem) {
    sub <- subdir(if (is.na(column)) NULL else column)
    stem <- create_safe_filename(stem)
    index[[length(index) + 1L]] <<- data.frame(
      file = if (nzchar(sub)) file.path(sub, stem) else stem,
      kind = if (is.na(column)) "polyfunctionality" else "polyfunctionality_split",
      column = column, level = level, reference = NA_character_,
      stringsAsFactors = FALSE)
    jobs[[length(jobs) + 1L]] <<- list(rows = rows, title = title, sub = sub,
                                       stem = stem)
  }
  add_job(seq_len(nrow(cells)), pop_name, NA_character_, NA_character_,
          paste0("umap_", pop_tag, "_polyfunctionality"))
  for (col in split_cols) {
    acol <- paste0("annot_", col)
    if (!acol %in% names(cells)) next
    vals <- as.character(cells[[acol]])
    levels <- if (exists("ordered_levels", mode = "function"))
      ordered_levels(vals, level_order[[col]]) else sort(unique(vals[!is.na(vals)]))
    for (lv in levels) {
      add_job(which(vals == lv), sprintf("%s, %s: %s", pop_name, col, lv), col, lv,
              paste0("umap_", pop_tag, "_polyfunctionality_", col, "_", lv))
    }
  }

  run_job <- function(j) draw(j$rows, j$title, j$sub, j$stem)
  if (n_cores > 1L && length(jobs) > 1L && exists("parallel_lapply", mode = "function")) {
    parallel_lapply(jobs, run_job, n_cores = n_cores)
  } else {
    lapply(jobs, run_job)
  }

  annot_cols <- sub("^annot_", "", grep("^annot_", names(cells), value = TRUE))
  tab_path <- file.path(output_dir, subdir(), "polyfunctionality_per_sample.csv")
  dir.create(dirname(tab_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(polyfunctionality_table(cells, deg, k, annot_cols), tab_path,
                   row.names = FALSE)
  cat(sprintf("Polyfunctionality: %d map(s) of 0 to %d positive gates (%s), and %s\n",
              length(jobs), k, paste(pf$markers, collapse = ", "),
              file.path(subdir(), "polyfunctionality_per_sample.csv")))
  invisible(do.call(rbind, index))
}
