# =============================================================================
# visualization/umap_split_plots.R: styled figures of the gate-layout UMAP,
# the whole embedding by combination, one figure per metadata level, a panel
# of all levels, maps of each combination's change in frequency against a
# reference level, and one figure per marker
# =============================================================================
# The embedding is computed once over every sample and drawn many times. Every
# figure uses the same limits and an equal aspect ratio, so an island sits at
# the same place on every page and two levels can be compared by eye. A level
# figure shows that level's cells coloured by combination over a pale grey
# silhouette of all cells, so an island that is empty at that level is still
# visibly there.
#
# Frequencies shown in legends and written to combination_frequencies_<col>.csv
# come from the counts over ALL cells (taken before sampling), as a mean over
# samples of each sample's percentage, so every sample counts equally whatever
# its size. Without counts they are rebuilt from the sampling weights of the
# plotted cells, which carry the same information for the sampled cells.
#
# Plan: docs/plans/umap_gate_layout_plan.md, sections 6 and 7. Look and colours:
# umap_style.R. Density grid of the outline: umap_gate_layout.R.
#
# Needs, already sourced: umap_style.R, umap_gate_layout.R (ordered_levels,
# balanced_fit_index, weighted_density_grid),
# R_shared/plot_utils.R (save_plot_formats, create_safe_filename) and
# R_shared/io_utils.R (make_outdir). R_shared/memory_utils.R (parallel_lapply)
# and emit_progress (io_utils.R) are used when present. Nothing here needs the
# %||% polyfill, so it runs on R 4.3 as sourced.
# =============================================================================

# Mean, sd and number of samples of each combination's percentage per level.
# `pct` is samples x combinations (rows named by sample), `level_of` names
# each sample's level. Returns one row per (level, combination).
combination_frequency_table <- function(pct, level_of, levels, column) {
  rows <- lapply(levels, function(lv) {
    s <- intersect(names(level_of)[level_of == lv], rownames(pct))
    m <- pct[s, , drop = FALSE]
    data.frame(column = column, level = lv, combination = colnames(pct),
               mean_pct = if (length(s)) colMeans(m) else NA_real_,
               sd_pct = if (length(s) > 1L) apply(m, 2, stats::sd) else NA_real_,
               n_samples = length(s), stringsAsFactors = FALSE,
               row.names = NULL)
  })
  do.call(rbind, rows)
}

# Per-sample percentage matrix (samples x combinations) over the samples in
# `cells`. With `counts`, combinations the cells do not know (pooled as rare in
# the plot) are folded into "other (rare)" and missing (sample, combination)
# pairs count as 0. Without counts, each sample's weighted plotted cells are
# used: a cell's weight is the number of cells it stands for in its sample.
# The cells each sample's percentages are taken of are kept in the attribute
# "population" (named by sample), for the pseudo-count of the change maps.
.sp_percent_matrix <- function(cells, counts, combos) {
  samples <- unique(as.character(cells$sample))
  if (is.null(counts) || !nrow(counts)) {
    w <- if (is.null(cells$weight)) rep(1, nrow(cells)) else cells$weight
    tab <- tapply(w, list(factor(cells$sample, levels = samples),
                          factor(cells$combination, levels = combos)), sum)
    tab[is.na(tab)] <- 0
    pop <- rowSums(tab)
    pct <- 100 * tab / ifelse(pop > 0, pop, NA_real_)
  } else {
    cb <- as.character(counts$combination)
    cb[!cb %in% combos] <- .US_RARE
    keep <- as.character(counts$sample) %in% samples
    tab <- tapply(as.numeric(counts$cells[keep]),
                  list(factor(counts$sample[keep], levels = samples),
                       factor(cb[keep], levels = union(combos, .US_RARE))), sum)
    tab[is.na(tab)] <- 0
    tab <- tab[, combos, drop = FALSE]
    pop <- as.numeric(tapply(as.numeric(counts$population_cells[keep]),
                             factor(counts$sample[keep], levels = samples), max))
    pct <- 100 * tab / ifelse(!is.na(pop) & pop > 0, pop, NA_real_)
  }
  keep <- stats::complete.cases(pct)
  pct <- pct[keep, , drop = FALSE]
  dimnames(pct) <- list(rownames(pct), combos)
  attr(pct, "population") <- stats::setNames(as.numeric(pop)[keep], rownames(pct))
  pct
}

# Change of every combination's frequency at `level` against `reference`,
# from the per-level table of combination_frequency_table(): the log2 ratio of
# the mean % per sample, with `pseudo_pct` added to both means. The pseudo-
# count gives a combination absent at one level a finite value and pulls
# changes of well under one cell per sample towards 0; a combination absent
# at both levels has no change to show and gets NA.
combination_change_table <- function(freq, level, reference, pseudo_pct) {
  a <- freq[freq$level == level, , drop = FALSE]
  b <- freq[freq$level == reference, , drop = FALSE]
  b <- b[match(a$combination, b$combination), , drop = FALSE]
  lfc <- log2((a$mean_pct + pseudo_pct) / (b$mean_pct + pseudo_pct))
  lfc[which(a$mean_pct == 0 & b$mean_pct == 0)] <- NA_real_
  data.frame(column = a$column, level = rep(level, nrow(a)),
             reference = rep(reference, nrow(a)), combination = a$combination,
             mean_pct_level = a$mean_pct, mean_pct_reference = b$mean_pct,
             log2_fold_change = lfc, n_samples_level = a$n_samples,
             n_samples_reference = b$n_samples,
             pseudo_pct = rep(pseudo_pct, nrow(a)), stringsAsFactors = FALSE,
             row.names = NULL)
}

# The run folder's subfolder for a figure or table (run_subdir() of
# umap_from_flowjo.R); "" (flat) when that file is not sourced.
.sp_subdir <- function(kind, column = NULL) {
  if (!exists("run_subdir", mode = "function")) return("")
  run_subdir(kind, column = if (is.null(column) || is.na(column)) NULL else column)
}

.sp_fmt_pct <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p >= 10, sprintf("%.0f%%", p),
  ifelse(p >= 1, sprintf("%.1f%%", p),
  ifelse(p >= 0.01, sprintf("%.2f%%", p), "<0.01%"))))
}

# Evaluate `expr` under a fixed seed without disturbing the caller's RNG
# stream (the embedding and sampling upstream depend on it).
.sp_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

# Row order for drawing: the largest combination first and the rarest last,
# so small islands are not buried under big ones; within a combination the
# rows are shuffled (deterministically) so that no sample is drawn on top of
# the others.
.sp_draw_order <- function(comb, seed = 42L) {
  sizes <- table(comb)
  ranked <- names(sizes)[order(-as.numeric(sizes), names(sizes), method = "radix")]
  shuffle <- .sp_with_seed(seed, stats::runif(length(comb)))
  order(match(comb, ranked), shuffle)
}

# Shared limits: the range of ALL cells plus a margin, wider at the lower left
# where the axis arrows go.
.sp_limits <- function(x, y) {
  rx <- range(x, na.rm = TRUE); ry <- range(y, na.rm = TRUE)
  dx <- max(diff(rx), 1e-6); dy <- max(diff(ry), 1e-6)
  list(x = c(rx[1] - 0.07 * dx, rx[2] + 0.03 * dx),
       y = c(ry[1] - 0.07 * dy, ry[2] + 0.03 * dy))
}

# Legend geometry: about 0.17 in per key row, wrapping into more columns
# instead of growing the figure, and the extra width those columns need.
.sp_legend_ncol <- function(n_keys, height) {
  per_col <- max(4L, floor((height - 1.2) / 0.17))
  max(1L, as.integer(ceiling(n_keys / per_col)))
}
.sp_legend_width <- function(labels, ncol) {
  longest <- max(nchar(as.character(labels)), 8L)
  ncol * (0.35 + 0.055 * longest) + 0.3
}

# Per-combination point size and alpha. The balanced fit gives every
# combination comparable room, so island densities on the page differ by
# orders of magnitude (200,000 cells of one combination in about the area of
# 50 of another). One size and alpha for all would either bury the rare
# islands or turn the abundant ones into solid blots, so each combination gets
# the automatic size and alpha of a figure as dense as its own island: its
# cells over the share of the panel it occupies (`occ`, bins of a 64 x 64 grid
# holding any of its cells), relative to the 30 % the automatic curves assume.
# `scale` corrects for a panel smaller than a full figure. Explicit
# point_size / point_alpha apply to every combination.
.sp_point_style <- function(n_by, occ, scale = 1, point_size = NULL,
                            point_alpha = NULL) {
  cb <- names(n_by)
  n_eq <- pmax(as.numeric(n_by), 1) * 0.3 / pmax(occ[cb], 1 / 4096) * scale
  data.frame(
    size  = if (is.null(point_size)) vapply(n_eq, auto_point_size, 0) else rep(point_size, length(cb)),
    alpha = if (is.null(point_alpha)) vapply(n_eq, auto_point_alpha, 0) else rep(point_alpha, length(cb)),
    row.names = cb)
}

# Share of the panel each combination occupies (named by combination).
.sp_occupancy <- function(x, y, comb, lims, g = 64L) {
  ix <- pmin(pmax(floor((x - lims$x[1]) / diff(lims$x) * g), 0), g - 1)
  iy <- pmin(pmax(floor((y - lims$y[1]) / diff(lims$y) * g), 0), g - 1)
  bins <- unique(data.frame(comb = comb, b = ix * g + iy, stringsAsFactors = FALSE))
  occ <- table(bins$comb) / g^2
  stats::setNames(as.numeric(occ), names(occ))
}

# The combination layer and its colour scale. `df` carries .size and .alpha
# per row; `labels` maps every combination in `breaks` to its legend text.
.sp_combination_layers <- function(df, cmap, breaks, labels, legend_title,
                                   ncol) {
  list(
    umap_point_layer(ggplot2::aes(x = UMAP1, y = UMAP2, colour = combination,
                                  size = .size, alpha = .alpha), data = df),
    ggplot2::scale_size_identity(),
    ggplot2::scale_alpha_identity(),
    ggplot2::scale_colour_manual(
      values = cmap, breaks = breaks, labels = unname(labels[breaks]),
      name = legend_title,
      guide = ggplot2::guide_legend(ncol = ncol,
                                    override.aes = list(size = 2.4, alpha = 1)))
  )
}

# Attach the per-combination style to a frame of cells.
.sp_styled <- function(df, style) {
  df$.size  <- style$size[match(df$combination, rownames(style))]
  df$.alpha <- style$alpha[match(df$combination, rownames(style))]
  df
}

# Arrow length as a fraction of the shorter axis: `base`, but never shorter
# on the page than `min_in` inches, so the "UMAP2" label (written along the
# vertical arrow) does not overrun it in a small or flat panel. `width_in`
# and `height_in` are the room the panel has; coord_fixed fits the limits
# into it at the smaller of the two scales.
.sp_arrow_frac <- function(lims, width_in, height_in, base = 0.07, min_in = 0.45) {
  dx <- diff(lims$x); dy <- diff(lims$y)
  scale <- min(width_in / dx, height_in / dy)
  max(base, min_in / (scale * min(dx, dy)))
}

# Smoothing bandwidth (fraction of the axis range) for a density over n
# cells: the given floor for a dense embedding, wider when the cells are few,
# so that a small run gives a smooth field rather than one blob per cell.
# sqrt(0.3 / n) is the typical spacing of n cells spread over 30 % of the
# panel.
.sp_bandwidth <- function(n, floor, mult) {
  max(floor, mult * sqrt(0.3 / max(n, 1)))
}

# Colour limits of a marker: its full finite range, or NULL (ggplot chooses)
# when that is empty or a single value. The same rule as color_limits() in
# umap_from_flowjo.R, which sets the scale of the per-file marker plots; it is
# repeated here so that this file does not have to source that one.
.sp_colour_limits <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NULL)
  lims <- range(x)
  if (!all(is.finite(lims)) || lims[1] >= lims[2]) return(NULL)
  lims
}

.sp_frame <- function(lims, clip = "off") {
  ggplot2::coord_fixed(ratio = 1, xlim = lims$x, ylim = lims$y, expand = FALSE,
                       clip = clip)
}

# -----------------------------------------------------------------------------
# Entry point. See the "Plotting entry point" section of the spec in the plan
# for the contract; `cells` holds one row per embedded cell with UMAP1, UMAP2,
# sample, workspace, combination, weight, fit and annot_<col> columns, and
# `colour_matrix` the FlowJo values of `use_markers` the marker colours are
# taken from, row-aligned to it: raw, or scaled per marker as named by
# `normalize_umap` ("none", "zscore" or "minmax", as --normalize_umap).
# Returns (invisibly) a data.frame of the figures written: file (stem), kind,
# column, level, reference.
# -----------------------------------------------------------------------------
plot_gate_layout <- function(cells, colour_matrix, use_markers, population,
                             output_dir, split_cols = character(0),
                             level_order = list(), split_reference = list(),
                             counts = NULL, point_size = NULL,
                             point_alpha = NULL, plot_width = 10,
                             plot_height = 8,
                             output_formats = list(pdf = FALSE, png = TRUE,
                                                   svg = FALSE),
                             n_cores = 1L, normalize_umap = "none",
                             background = NULL) {
  normalize_umap <- match.arg(normalize_umap, c("none", "zscore", "minmax"))
  # ggrastr draws random numbers while it rasterises a layer, which would move
  # the caller's RNG stream; drawing figures must not change what a later
  # random step of the run does, so the state is put back on exit.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rng_state <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", rng_state, envir = globalenv()), add = TRUE)
  }
  make_outdir(output_dir)
  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  cells$combination <- as.character(cells$combination)
  cells$sample <- as.character(cells$sample)
  if (is.null(cells$weight)) cells$weight <- 1
  n_all <- nrow(cells)
  pop_tag  <- create_safe_filename(sub("^/", "", population))
  pop_name <- basename(sub("/$", "", population))
  if (!nzchar(pop_name)) pop_name <- population

  lims <- .sp_limits(cells$UMAP1, cells$UMAP2)
  cmap <- combination_colours(cells$combination)
  combos <- names(cmap)
  ord_all <- .sp_draw_order(cells$combination)
  arrows <- umap_axis_arrows(lims$x, lims$y, frac = .sp_arrow_frac(
    lims, plot_width - 0.5, plot_height - 1.4))
  theme_u <- theme_umap_clean()
  pct_all <- .sp_percent_matrix(cells, counts, combos)
  # Background subtraction (background_pairs(): sample -> its unstimulated
  # tube): every frequency below is then the stimulated tube's minus its
  # background's, negatives 0, and the background tubes drop out of the
  # frequencies and of the split plots. The whole-embedding plot still shows
  # every cell.
  bg_on <- length(background) > 0L
  if (bg_on) {
    pct_all <- subtract_background(pct_all, background,
                                   subtract = response_combinations(combos))
    if (!nrow(pct_all)) {
      stop("background subtraction left no sample: none of the paired ",
           "samples has both its own and its background's counts.", call. = FALSE)
    }
  }
  in_freq <- if (bg_on) cells$sample %in% rownames(pct_all) else rep(TRUE, nrow(cells))

  occ <- .sp_occupancy(cells$UMAP1, cells$UMAP2, cells$combination, lims)
  style_of <- function(n_by, scale = 1) {
    .sp_point_style(n_by, occ, scale, point_size, point_alpha)
  }

  # The pale grey silhouette of the whole embedding drawn under every level.
  # A capped subsample, balanced across combinations so that rare islands get
  # a silhouette too, drawn opaque so that overlap does not show.
  sil_idx <- .sp_with_seed(7L, balanced_fit_index(cells$combination, 150000L))
  sil <- cells[sil_idx, c("UMAP1", "UMAP2", "combination")]
  sil_style <- style_of(table(sil$combination))
  sil$.size <- 1.5 * sil_style$size[match(sil$combination, rownames(sil_style))]
  # Its size scale is the identity scale the combination layer adds.
  sil_layer <- function(data) {
    umap_point_layer(ggplot2::aes(x = UMAP1, y = UMAP2, size = .size),
                     data = data, alpha = 1, colour = "#EBEBEB")
  }

  # File stems are claimed up front, in order, so that names truncated to the
  # same 100 characters get _2, _3 rather than overwrite each other, and so
  # that the parallel jobs below never race for a name.
  used <- character(0)
  claim <- function(stem) {
    s <- create_safe_filename(stem)
    base <- s; i <- 1L
    while (s %in% used) {
      i <- i + 1L
      suffix <- paste0("_", i)
      s <- paste0(substr(base, 1L, 100L - nchar(suffix)), suffix)
    }
    used <<- c(used, s)
    s
  }
  index <- list()
  jobs  <- list()
  # Each figure goes into its kind's folder of the run (run_subdir() in
  # umap_from_flowjo.R; flat when that is not sourced). The index keeps the
  # path relative to the run folder.
  sub_of <- stats::setNames(character(0), character(0))
  add_job <- function(stem, kind, column, level, reference, fn) {
    sub <- .sp_subdir(if (kind == "combination") "overview" else kind, column)
    sub_of[[stem]] <<- sub
    index[[length(index) + 1L]] <<- data.frame(
      file = if (nzchar(sub)) file.path(sub, stem) else stem,
      kind = kind, column = column, level = level,
      reference = reference, stringsAsFactors = FALSE)
    jobs[[length(jobs) + 1L]] <<- list(stem = stem, fn = fn)
  }
  save_fig <- function(p, stem, width, height) {
    d <- file.path(output_dir, sub_of[[stem]])
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    save_plot_formats(p, file.path(d, stem), width = width,
                      height = height, formats = output_formats)
  }

  # --- Whole embedding, coloured by combination. ------------------------------
  freq_all <- if (nrow(pct_all)) colMeans(pct_all, na.rm = FALSE) else stats::setNames(rep(NA_real_, length(combos)), combos)
  lab_all <- stats::setNames(paste0(combos, "  ", .sp_fmt_pct(freq_all[combos])), combos)
  legend_title <- paste0("Mean % of ", pop_name, " per sample",
                         if (bg_on) " minus background" else "")
  stem_comb <- claim(paste0("umap_", pop_tag, "_by_combination"))
  add_job(stem_comb, "combination", NA_character_, NA_character_, NA_character_,
    function() {
      df <- .sp_styled(cells[ord_all, c("UMAP1", "UMAP2", "combination")],
                       style_of(table(cells$combination)))
      ncol <- .sp_legend_ncol(length(combos), plot_height)
      p <- ggplot2::ggplot() +
        .sp_combination_layers(df, cmap, combos, lab_all, legend_title, ncol) +
        arrows + .sp_frame(lims) + theme_u +
        ggplot2::labs(title = pop_name, subtitle = sprintf(
          "%s cells from %d samples, coloured by gate combination",
          format(n_all, big.mark = ","), length(unique(cells$sample))))
      save_fig(p, stem_comb, plot_width + .sp_legend_width(lab_all, ncol),
               plot_height)
    })

  # --- Per metadata column: levels, panel, changes, frequency tables. -------
  # Occupancy of the embedding (where any cell lies, whatever its weight),
  # shared by every change map, where its contour at 0.15 is the grey outline
  # of the embedding; rare islands are outlined as well as abundant ones. The
  # change maps draw the silhouette's cells, so an island has the same shape
  # on every map, with the rarest islands drawn last, in points a little
  # larger than the silhouette's (up to its largest) so that each island
  # reads as one filled shape of one colour.
  outline <- NULL
  if (length(split_cols)) {
    gn <- 256L
    ix <- pmin(pmax(round((cells$UMAP1 - lims$x[1]) / diff(lims$x) * (gn - 1)) + 1, 1), gn)
    iy <- pmin(pmax(round((cells$UMAP2 - lims$y[1]) / diff(lims$y) * (gn - 1)) + 1, 1), gn)
    occ_bins <- unique(cbind(ix, iy))
    gx <- seq(lims$x[1], lims$x[2], length.out = gn)
    gy <- seq(lims$y[1], lims$y[2], length.out = gn)
    og <- weighted_density_grid(gx[occ_bins[, 1]], gy[occ_bins[, 2]], NULL,
                                lims$x, lims$y, n = gn,
                                bandwidth = .sp_bandwidth(n_all, 0.006, 1.5))
    outline <- data.frame(x = rep(og$x, times = length(og$y)),
                          y = rep(og$y, each = length(og$x)),
                          occ = as.vector(og$z) / max(og$z))
    chg_cells <- sil[.sp_draw_order(sil$combination), ]
    chg_cells$.size <- pmin(1.4 * chg_cells$.size, 1.65)
  }

  for (col in split_cols) {
    acol <- paste0("annot_", col)
    if (!acol %in% names(cells)) {
      warning("split column '", col, "' is not among the cell annotations; skipped.",
              call. = FALSE)
      next
    }
    vals <- as.character(cells[[acol]])
    # Background tubes are not drawn per level: they were subtracted.
    vals[!in_freq] <- NA_character_
    levels <- ordered_levels(vals, level_order[[col]])
    if (!length(levels)) {
      warning("split column '", col, "' has no assigned levels; skipped.", call. = FALSE)
      next
    }
    ref <- split_reference[[col]]
    if (is.null(ref) || !ref %in% levels) {
      if (!is.null(ref)) {
        warning("reference level '", ref, "' of '", col, "' is not present; using '",
                levels[1], "'.", call. = FALSE)
      }
      ref <- levels[1]
    }

    # Level of each sample, and the frequency table for this column.
    s_lv <- tapply(vals, cells$sample, function(v) v[1])
    level_of <- stats::setNames(as.character(s_lv), names(s_lv))
    freq <- combination_frequency_table(pct_all, level_of, levels, col)
    freq$background_subtracted <- rep(bg_on, nrow(freq))
    tab_dir <- file.path(output_dir, .sp_subdir("table", col))
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(freq, file.path(tab_dir, paste0(
      "combination_frequencies_", create_safe_filename(col), ".csv")),
      row.names = FALSE)

    col_tag <- create_safe_filename(col)
    # One size and alpha per combination for all levels of the column, from
    # its mean cells per level, so that an island thinner at one level than
    # at another also looks thinner.
    in_lv <- vals %in% levels
    n_typ <- table(factor(cells$combination[in_lv], levels = combos)) / length(levels)
    n_typ <- n_typ[n_typ > 0]
    style_col <- style_of(n_typ)

    # One figure per level.
    for (lv in levels) {
      stem <- claim(paste0("umap_", pop_tag, "_split_", col_tag, "_", lv))
      add_job(stem, "split", col, lv, ref, local({
        stem <- stem; lv <- lv; col <- col; vals <- vals; freq <- freq
        style_col <- style_col
        function() {
          ord <- ord_all[which(vals[ord_all] == lv)]
          df <- .sp_styled(cells[ord, c("UMAP1", "UMAP2", "combination")], style_col)
          f <- freq[freq$level == lv, ]
          present <- combos[combos %in% unique(df$combination)]
          lab <- stats::setNames(paste0(present, "  ",
                                        .sp_fmt_pct(f$mean_pct[match(present, f$combination)])),
                                 present)
          ncol <- .sp_legend_ncol(length(present), plot_height)
          n_s <- f$n_samples[1]
          p <- ggplot2::ggplot() + sil_layer(sil) +
            .sp_combination_layers(df, cmap, present, lab, legend_title, ncol) +
            arrows + .sp_frame(lims) + theme_u +
            ggplot2::labs(title = lv, subtitle = sprintf(
              "%s: %d sample%s, %s cells shown (grey: all cells)", col, n_s,
              if (n_s == 1L) "" else "s", format(nrow(df), big.mark = ",")))
          save_fig(p, stem, plot_width + .sp_legend_width(lab, ncol), plot_height)
        }
      }))
    }

    # All levels side by side, one shared legend.
    stem_panel <- claim(paste0("umap_", pop_tag, "_split_", col_tag, "_panel"))
    add_job(stem_panel, "split_panel", col, NA_character_, ref, local({
      stem_panel <- stem_panel; col <- col; vals <- vals; levels <- levels
      n_typ <- n_typ
      function() {
        nl <- length(levels)
        # Near-square grid (5 levels: 3 x 2), so each facet stays large enough
        # to read an island; a single row of five was 4800 px wide and flat.
        ncol_f <- if (nl <= 3L) nl else as.integer(ceiling(sqrt(nl)))
        nrow_f <- as.integer(ceiling(nl / ncol_f))
        facet_in <- 4.2
        asp <- diff(lims$y) / diff(lims$x)
        # The same cells per facet would crowd a panel this much smaller than
        # a full figure, so size and alpha are chosen as for a denser figure.
        shrink <- (8 / facet_in)^2
        ord <- ord_all[vals[ord_all] %in% levels]
        df <- .sp_styled(cells[ord, c("UMAP1", "UMAP2", "combination")],
                         style_of(n_typ, shrink))
        df$.level <- factor(vals[ord], levels = levels)
        sil_s <- sil[seq(1L, nrow(sil), by = max(1L, nrow(sil) %/% 60000L)), ]
        sil_s$.size <- sil_s$.size / sqrt(shrink) * 1.6
        sil_f <- do.call(rbind, lapply(levels, function(lv)
          data.frame(sil_s, .level = lv, stringsAsFactors = FALSE)))
        sil_f$.level <- factor(sil_f$.level, levels = levels)
        present <- combos[combos %in% unique(df$combination)]
        ncol <- .sp_legend_ncol(length(present), nrow_f * facet_in * asp + 1)
        p <- ggplot2::ggplot() + sil_layer(sil_f) +
          .sp_combination_layers(df, cmap, present,
                                 stats::setNames(present, present),
                                 "Combination", ncol) +
          ggplot2::facet_wrap(~ .level, ncol = ncol_f) +
          # A facet is a third of a figure wide: a longer arrow (relative to
          # the panel) and smaller labels keep "UMAP2" within its arrow.
          umap_axis_arrows(lims$x, lims$y, base_size = 8, frac = .sp_arrow_frac(
            lims, facet_in, facet_in * asp, base = 0.1, min_in = 0.4)) +
          .sp_frame(lims) + theme_u +
          ggplot2::labs(title = paste0(pop_name, " by ", col))
        save_fig(p, stem_panel,
                 ncol_f * facet_in + .sp_legend_width(present, ncol) + 0.4,
                 max(3, nrow_f * (facet_in * asp + 0.35) + 0.9))
      }
    }))

    # How the frequency of every combination changes at each level against
    # the reference. Each island is coloured by the log2 fold change of its
    # combination's mean % per sample, from the same per-sample percentages
    # as the frequency table (the exact counts when given). A map of cell
    # density differences would be set by the abundant lineage islands, and
    # every cytokine island would be near-white however much it changed. The
    # pseudo-count added to both means is one cell in the median sample of
    # the two levels, so changes of well under a cell per sample stay pale.
    pop_s <- attr(pct_all, "population")
    changes <- list()
    for (lv in setdiff(levels, ref)) {
      s_two <- intersect(names(level_of)[level_of %in% c(lv, ref)], names(pop_s))
      pseudo <- if (length(s_two)) 100 / stats::median(pop_s[s_two]) else NA_real_
      if (!is.finite(pseudo) || pseudo <= 0) pseudo <- 0.01
      chg <- combination_change_table(freq, lv, ref, pseudo)
      changes[[length(changes) + 1L]] <- chg
      stem <- claim(paste0("umap_", pop_tag, "_diff_", col_tag, "_", lv))
      add_job(stem, "diff", col, lv, ref, local({
        stem <- stem; lv <- lv; ref <- ref; chg <- chg; pseudo <- pseudo
        function() {
          df <- chg_cells
          df$lfc <- chg$log2_fold_change[match(df$combination, chg$combination)]
          p <- ggplot2::ggplot() +
            umap_point_layer(ggplot2::aes(x = UMAP1, y = UMAP2, colour = lfc,
                                          size = .size), data = df, alpha = 1) +
            ggplot2::scale_size_identity() +
            ggplot2::geom_contour(data = outline, ggplot2::aes(x = x, y = y, z = occ),
                                  breaks = 0.15, colour = "grey70", linewidth = 0.3) +
            umap_diverging_colour(3, name = "Fold change", breaks = -3:3,
                                  labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8")) +
            arrows + .sp_frame(lims) + theme_u +
            ggplot2::labs(
              title = sprintf("Change in frequency, %s vs %s", lv, ref),
              subtitle = sprintf("log2 fold change of the mean %%%s per sample; %d / %d samples",
                                 if (bg_on) " (background-subtracted)" else "",
                                 chg$n_samples_level[1], chg$n_samples_reference[1]),
              caption = sprintf(paste0(
                if (bg_on) "Grey: absent at both levels, or not a cytokine response (not subtracted). "
                else "Grey: absent at both levels. ",
                "Colours stop at 8-fold. ",
                "%s %% (one cell in the median sample) is added to both means.%s"),
                format(signif(pseudo, 2), scientific = FALSE),
                if (bg_on) paste0("\nEach stimulated tube minus the unstimulated ",
                                  "tube of the same patient and visit, negatives 0.")
                else ""))
          save_fig(p, stem, plot_width + 1.2, plot_height)
        }
      }))
    }
    if (length(changes)) {
      utils::write.csv(do.call(rbind, changes), file.path(tab_dir, paste0(
        "combination_change_", col_tag, ".csv")), row.names = FALSE)
    }
  }

  # --- One figure per marker, on the per-file plots' colour scale. ----------
  # plot_flowjo_umap_facets() has drawn these under the same names already,
  # in viridis over the full range of each marker across all cells; this
  # redraw keeps that scale, so the styled conjoint plot and the per-file
  # plots of a marker read against one colour bar.
  norm_tag <- switch(normalize_umap, none = "", zscore = " (z-score)",
                     minmax = " (min-max)")
  norm_sub <- switch(normalize_umap,
                     none   = "FlowJo value (transformed scale)",
                     zscore = "FlowJo value (transformed scale), z-scored over all cells",
                     minmax = "FlowJo value (transformed scale), min-max scaled over all cells")
  if (!is.null(colour_matrix) && length(use_markers)) {
    for (m in intersect(use_markers, colnames(colour_matrix))) {
      stem <- claim(paste0("umap_", pop_tag, "_marker_", m, "_ALL"))
      add_job(stem, "marker", NA_character_, NA_character_, NA_character_, local({
        stem <- stem; m <- m
        function() {
          v <- colour_matrix[, m]
          # High values on top, so the positive cells of a marker show.
          o <- order(v, na.last = NA)
          df <- .sp_styled(data.frame(UMAP1 = cells$UMAP1[o], UMAP2 = cells$UMAP2[o],
                                      combination = cells$combination[o], expr = v[o],
                                      stringsAsFactors = FALSE),
                           style_of(table(cells$combination)))
          p <- ggplot2::ggplot() +
            umap_point_layer(ggplot2::aes(x = UMAP1, y = UMAP2, colour = expr,
                                          size = .size, alpha = .alpha), data = df) +
            ggplot2::scale_size_identity() + ggplot2::scale_alpha_identity() +
            umap_marker_colour(.sp_colour_limits(v), name = paste0(m, norm_tag)) +
            arrows + .sp_frame(lims) + theme_u +
            ggplot2::labs(title = m, subtitle = norm_sub)
          save_fig(p, stem, plot_width + 1.2, plot_height)
        }
      }))
    }
  }

  # --- Draw. Every job writes one figure under a name claimed above. --------
  n_jobs <- length(jobs)
  run_job <- function(i) {
    jobs[[i]]$fn()
    if (exists("emit_progress", mode = "function")) {
      emit_progress("plot", i, n_jobs, jobs[[i]]$stem)
    }
    jobs[[i]]$stem
  }
  if (n_cores > 1L && n_jobs > 1L && exists("parallel_lapply", mode = "function")) {
    parallel_lapply(seq_len(n_jobs), run_job, n_cores = n_cores)
  } else {
    lapply(seq_len(n_jobs), run_job)
  }

  idx <- if (length(index)) do.call(rbind, index) else
    data.frame(file = character(0), kind = character(0), column = character(0),
               level = character(0), reference = character(0),
               stringsAsFactors = FALSE)
  invisible(idx)
}
