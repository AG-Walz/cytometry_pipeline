# =============================================================================
# visualization/umap_style.R: the look of the gate-layout UMAP figures
# =============================================================================
# One place for the palette, theme, axis arrows and point drawing used by
# umap_split_plots.R, so every figure of a run (overall, per metadata level,
# panel, change in frequency, marker) reads as one set.
#
# The target is the figure the user supplied (/mnt/volume/data/pasted-image.png):
# muted pastel colours, small translucent points that read as soft clouds, no
# axes, grid or box, a short L-shaped UMAP1/UMAP2 arrow at the lower left, a
# small sans title and a legend with small round keys. There are deliberately
# no leader lines, centroid dots or text labels on the clusters: the user asked
# for the palette and the island layout, not the annotation.
#
# Plan: docs/plans/umap_gate_layout_plan.md, section 7.
#
# Everything here is pure (base R plus ggplot2 for the layer constructors).
# ggrastr and ragg are optional: points are rasterised when ggrastr is
# installed, and plain vector points are drawn otherwise.
# =============================================================================

.US_RARE      <- "other (rare)"
.US_UNDEFINED <- "#D9D9D9"   # light grey of "Undefined" in the reference

# Twelve colours sampled from the reference figure (its clusters, read where
# the points are dense enough to show the full colour), ordered so that
# neighbouring entries differ in hue AND lightness.
.US_BASE_PALETTE <- c(
  "#D4A19C",  # dusty rose
  "#9DB8A3",  # sage
  "#B3A1C9",  # lavender
  "#D1A961",  # ochre / mustard
  "#7FAAAD",  # teal
  "#A85B55",  # brick red
  "#8FA5C8",  # slate blue
  "#5F8A6C",  # dark green
  "#DCC9A6",  # sand
  "#4A5F8C",  # navy
  "#E3A68A",  # salmon
  "#957240"   # dark ochre
)

# CIE-Lab coordinates of colours, one row each. Rounded so that the greedy
# choices below do not hang on the last bits of the conversion, which can
# differ between BLAS builds: the same labels must give the same colours on
# every host.
.us_lab <- function(cols) {
  round(grDevices::convertColor(t(grDevices::col2rgb(cols)) / 255,
                                from = "sRGB", to = "Lab"), 6)
}

# Muted HCL candidates: every combination of the given hues, chromas and
# luminances that lies inside the sRGB gamut (hcl() returns NA otherwise) and
# is no more saturated (HSV) than the reference's most saturated colour, the
# dark ochre. Without that cap a dark teal or olive at the same chroma comes
# out close to a pure colour, far from the reference's muted register.
.us_candidates <- function(h, c, l) {
  g <- expand.grid(h = h %% 360, c = c, l = l)
  cand <- grDevices::hcl(h = g$h, c = g$c, l = g$l, fixup = FALSE)
  cand <- unique(cand[!is.na(cand)])
  rgb <- grDevices::col2rgb(cand)
  sat <- 1 - apply(rgb, 2, min) / pmax(apply(rgb, 2, max), 1)
  cand[sat <= 0.57]
}

# Pick `n` of the candidate colours one at a time, each time the one whose
# smallest CIE-Lab distance to every colour chosen so far (`chosen` included)
# is largest. Ties go to the earlier candidate, so the result is
# deterministic, and a longer pick extends a shorter one.
.us_maximin <- function(cand, n, chosen = character(0)) {
  cl <- t(.us_lab(cand))
  d <- rep(Inf, ncol(cl))
  for (h in chosen) d <- pmin(d, sqrt(colSums((cl - .us_lab(h)[1, ])^2)))
  out <- character(0)
  for (k in seq_len(n)) {
    j <- which.max(d)
    out <- c(out, cand[j])
    d <- pmin(d, sqrt(colSums((cl - cl[, j])^2)))
  }
  out
}

# n muted categorical colours. The first twelve are the reference palette.
# Beyond that each colour is the muted candidate (hue every 5 degrees, chroma
# 26 to 40, luminance 45 to 82) furthest in CIE-Lab from every colour already
# in the palette, so a late colour never lands next to an early one.
umap_muted_palette <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n <= 0L) return(character(0))
  if (n <= length(.US_BASE_PALETTE)) return(.US_BASE_PALETTE[seq_len(n)])
  cand <- .us_candidates(h = seq(0, 355, by = 5), c = c(26, 33, 40),
                         l = c(45, 52, 59, 66, 74, 82))
  c(.US_BASE_PALETTE,
    .us_maximin(cand, n - length(.US_BASE_PALETTE), .US_BASE_PALETTE))
}

# -----------------------------------------------------------------------------
# Colours for gate-combination labels such as "CD4+ | IFNy+ TNFa+".
#
# The part before the first " | " is the first-tier group (the lineage). A
# group with several labels is a family and gets its own sector of the hue
# wheel, so lineage reads at a glance: its "| none" label is the lightest,
# least saturated tint, and its positive combinations are mid to dark shades
# within the sector, darker the more markers are positive. A group whose first
# tier is all negative ("CD4- CD8-") is drawn in neutral greys, since it is
# defined by what it lacks. A group with a single label (e.g. "CD4+ CD8+") is
# an accent colour placed in a gap between the family sectors. Labels without
# " | " are single-label groups; when there are no families at all the plain
# reference palette is used. "other (rare)" is light grey.
#
# The result depends only on the SET of labels (they are sorted in a
# locale-independent way first), so a run gives the same colours on every
# host. Names are in display order (family by family, light to dark, then
# accents, greys, and "other (rare)" last), which is the legend order.
# -----------------------------------------------------------------------------
combination_colours <- function(labels) {
  lab <- unique(as.character(labels[!is.na(labels)]))
  if (!length(lab)) return(stats::setNames(character(0), character(0)))
  lab <- sort(lab, method = "radix")
  rare <- lab[lab == .US_RARE]
  lab  <- lab[lab != .US_RARE]

  bar <- regexpr(" | ", lab, fixed = TRUE)
  has_bar <- bar > 0L
  group <- lab
  rest  <- rep("", length(lab))
  group[has_bar] <- substr(lab[has_bar], 1L, bar[has_bar] - 1L)
  rest[has_bar]  <- substring(lab[has_bar], bar[has_bar] + 3L)
  # Positive markers below the first tier: tokens that end in "+".
  n_pos <- vapply(strsplit(rest, "[[:space:]|]+"), function(t)
    sum(grepl("\\+$", t)), integer(1))
  is_none <- has_bar & grepl("^none$", rest)
  # A first tier spelled entirely as negatives ("CD4- CD8-").
  tok1 <- strsplit(group, "[[:space:]]+")
  all_neg <- vapply(tok1, function(t) length(t) > 0L && all(grepl("-$", t)),
                    logical(1))
  # Cells no gate decides (with a population above the live gate: the cells
  # it excluded) get one fixed near-black slate, darker than any of the
  # negative greys and far from "other (rare)".
  undecided <- lab[group == "(no decision)"]
  keep <- group != "(no decision)"
  lab <- lab[keep]; group <- group[keep]; rest <- rest[keep]
  n_pos <- n_pos[keep]; is_none <- is_none[keep]; all_neg <- all_neg[keep]

  groups <- unique(group)
  size <- table(group)[groups]
  fam_groups    <- groups[!groups %in% group[all_neg] & size >= 2L]
  single_groups <- groups[!groups %in% group[all_neg] & size < 2L]
  grey_groups   <- unique(group[all_neg])

  # Within a group: zero-positive labels first ("| none" before a bare group
  # label), then by number of positive markers, then alphabetically.
  within_order <- function(idx) {
    idx[order(n_pos[idx], !is_none[idx], lab[idx], method = "radix")]
  }

  out <- character(0)
  k <- length(fam_groups)
  step <- if (k) 360 / k else 360
  width <- min(0.75 * step, 100)
  centre0 <- 30   # first family: rose, salmon, ochre and brick, the reference's warm end

  for (j in seq_len(k)) {
    idx <- within_order(which(group == fam_groups[j]))
    centre <- centre0 + (j - 1) * step
    zero <- idx[n_pos[idx] == 0L]
    pos  <- idx[n_pos[idx] > 0L]
    cols <- character(length(idx)); names(cols) <- lab[idx]
    # Zero-positive labels: pale tints of the family centre, the first one
    # (the "| none" label when present) the palest.
    if (length(zero)) {
      lz <- 84 - 10 * (seq_along(zero) - 1L)
      cols[lab[zero]] <- grDevices::hcl(h = centre %% 360, c = 18, l = pmax(lz, 60))
    }
    if (length(pos)) {
      max_p <- max(n_pos[pos])
      hues <- centre + width * seq(-0.5, 0.5, length.out = max(2, round(width / 2.5) + 1))
      for (p in sort(unique(n_pos[pos]))) {
        ip <- pos[n_pos[pos] == p]
        # Each number of positive markers has its own band of lightness, from
        # L 76 for one marker down to L 42 (the reference's navy and brick)
        # for the most, so darker still means more markers. Within the band
        # and the family's sector of hues, every label takes the muted
        # candidate furthest in CIE-Lab from every colour chosen so far, so
        # even a family of 16 labels stays tellable apart.
        top <- 76 - (p - 1) * 34 / max_p
        bottom <- 76 - p * 34 / max_p
        cand <- .us_candidates(h = hues, c = c(22, 28, 34, 40),
                               l = seq(bottom, top, length.out = max(2, round((top - bottom) / 2.5) + 1)))
        cols[lab[ip]] <- .us_maximin(cand, length(ip), c(out, cols[nzchar(cols)]))
      }
    }
    out <- c(out, cols)
  }

  if (length(single_groups)) {
    idx <- unlist(lapply(single_groups, function(g) within_order(which(group == g))))
    if (!k) {
      cols <- umap_muted_palette(length(idx))
    } else {
      # Accents sit in the gaps between family sectors, starting with the gap
      # after the last family (lavender when there are two families), each
      # the candidate in its gap furthest from every colour chosen so far.
      gaps <- (centre0 + (seq_len(k) - 0.5) * step) %% 360
      gaps <- gaps[c(k, seq_len(k - 1L))]
      cols <- character(0)
      for (i in seq_along(idx)) {
        g <- gaps[((i - 1L) %% k) + 1L]
        cand <- .us_candidates(h = g + (step - width) * seq(-0.3, 0.3, length.out = 7),
                               c = c(26, 30, 34), l = seq(48, 78, by = 3))
        cols <- c(cols, .us_maximin(cand, 1L, c(out, cols)))
      }
    }
    out <- c(out, stats::setNames(cols, lab[idx]))
  }

  if (length(grey_groups)) {
    idx <- unlist(lapply(grey_groups, function(g) within_order(which(group == g))))
    m <- length(idx)
    # Warm neutral greys, light to dark, kept clearly darker than the
    # "other (rare)" grey.
    l <- if (m == 1L) 70 else seq(72, 38, length.out = m)
    out <- c(out, stats::setNames(grDevices::hcl(h = 60, c = 4, l = l), lab[idx]))
  }
  if (length(undecided)) out <- c(out, stats::setNames("#2E3440", undecided))
  if (length(rare)) out <- c(out, stats::setNames(.US_UNDEFINED, .US_RARE))
  out
}

# -----------------------------------------------------------------------------
# Theme: theme_void with a small, centred sans title, a quiet subtitle and a
# compact legend. The font family itself is stamped by save_plot_formats().
# -----------------------------------------------------------------------------
theme_umap_clean <- function(base_size = 9) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title       = ggplot2::element_text(size = base_size + 1.5, colour = "grey15",
                                               hjust = 0.5, margin = ggplot2::margin(b = 4)),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, colour = "grey45",
                                               hjust = 0.5, margin = ggplot2::margin(b = 6)),
      plot.caption     = ggplot2::element_text(size = base_size - 2, colour = "grey50"),
      plot.margin      = ggplot2::margin(12, 12, 12, 12),
      legend.title     = ggplot2::element_text(size = base_size - 0.5, colour = "grey25",
                                               margin = ggplot2::margin(b = 4)),
      legend.text      = ggplot2::element_text(size = base_size - 1.5, colour = "grey25"),
      legend.key.size  = grid::unit(0.34, "cm"),
      legend.key.spacing.y = grid::unit(1, "pt"),
      legend.margin    = ggplot2::margin(0, 0, 0, 8),
      strip.text       = ggplot2::element_text(size = base_size, colour = "grey20",
                                               margin = ggplot2::margin(2, 0, 4, 0)),
      panel.spacing    = grid::unit(10, "pt")
    )
}

# -----------------------------------------------------------------------------
# The short L-shaped arrow pair at the lower left of the panel, labelled UMAP1
# (under the horizontal arrow) and UMAP2 (left of the vertical one). Both
# arrows have the same length in data units, which with coord_fixed means the
# same length on the page. `xlim`, `ylim` are the panel limits; the origin is
# their lower-left corner, so the limits should carry a margin there.
# -----------------------------------------------------------------------------
umap_axis_arrows <- function(xlim, ylim, frac = 0.12, base_size = 9,
                             colour = "grey20") {
  len <- frac * min(diff(xlim), diff(ylim))
  x0 <- xlim[1] + 0.02 * diff(xlim)
  y0 <- ylim[1] + 0.03 * diff(ylim)
  arr <- grid::arrow(length = grid::unit(0.06, "in"), type = "closed", angle = 22)
  txt <- (base_size - 2) / ggplot2::.pt
  list(
    ggplot2::annotate("segment", x = x0, y = y0, xend = x0 + len, yend = y0,
                      colour = colour, linewidth = 0.35, arrow = arr),
    ggplot2::annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + len,
                      colour = colour, linewidth = 0.35, arrow = arr),
    ggplot2::annotate("text", x = x0 + len / 2, y = y0, label = "UMAP1",
                      vjust = 1.6, size = txt, colour = "grey35"),
    ggplot2::annotate("text", x = x0, y = y0 + len / 2, label = "UMAP2",
                      vjust = -0.8, angle = 90, size = txt, colour = "grey35")
  )
}

# -----------------------------------------------------------------------------
# A point layer drawn as soft dots: filled circles without an outline (the
# default stroke would make every tiny point a hard-edged ring). Rasterised
# through ggrastr when it is installed, so SVG and PDF output stay small; PNG
# output is raster anyway, so nothing is lost without it. `size` and `alpha`
# may be NULL when the mapping carries them (with identity scales); extra
# arguments (e.g. a constant colour) go to geom_point.
# -----------------------------------------------------------------------------
umap_point_layer <- function(mapping = NULL, data = NULL, size = NULL,
                             alpha = NULL, ...) {
  args <- list(mapping = mapping, data = data, shape = 16, stroke = 0, ...)
  if (!is.null(size))  args$size  <- size
  if (!is.null(alpha)) args$alpha <- alpha
  layer <- do.call(ggplot2::geom_point, args)
  if (!requireNamespace("ggrastr", quietly = TRUE)) return(layer)
  dev <- if (requireNamespace("ragg", quietly = TRUE)) "ragg"
         else if (isTRUE(capabilities("cairo"))) "cairo" else NA_character_
  if (is.na(dev)) return(layer)
  ggrastr::rasterise(layer, dpi = 300, dev = dev)
}

# Point size and alpha for n points spread over a typical share (about 30 %)
# of an ~8-inch panel. The model: a ggplot point of size s (shape 16, no
# stroke) covers about 27.6 s^2 px^2 at 200 dpi, and the colour a region
# reaches is 1 - exp(-density * area * alpha). Alpha is set so that this
# exposure is about 1.5 (close to the full colour in the dense core, fading at
# the edges), and size shrinks with n so that a dense cloud is made of many
# faint specks rather than a few solid dots: that is what reads as a soft
# cloud. Few points get larger, nearly opaque dots. Tuned by eye against the
# reference figure.
auto_point_size <- function(n) {
  n <- max(1, n)
  min(1.1, max(0.18, 0.3 * (1e6 / n)^0.33))
}
auto_point_alpha <- function(n) {
  n <- max(1, n)
  min(0.85, max(0.04, 31300 / (n * auto_point_size(n)^2)))
}

# Colour bars are sized through the guide's own theme: there the length is
# taken as given, while a plot-level legend.key.height is multiplied by five
# for colour bars (ggplot2 3.5 and 4 alike).
.us_colourbar <- function(length_in = 1.5) {
  ggplot2::guide_colourbar(ticks = FALSE, theme = ggplot2::theme(
    legend.key.height = grid::unit(length_in, "in"),
    legend.key.width  = grid::unit(0.12, "in")))
}

# Colour scale for marker expression: viridis, the scale of the per-file
# marker plots of plot_flowjo_umap_facets() (viridis::scale_color_viridis()
# is this same gradient through 256 viridis colours). With the same limits a
# colour means the same value in the styled conjoint plot and in every
# per-file plot. Values beyond the limits are squished, not dropped; NULL
# limits leave the range to ggplot. viridisLite comes with ggplot2 (through
# scales), so this does not need the viridis package.
umap_marker_colour <- function(limits, name) {
  ggplot2::scale_colour_gradientn(
    colours = viridisLite::viridis(256), limits = limits,
    oob = scales::squish, name = name, guide = .us_colourbar())
}

# Diverging colour for the change maps, after panel E of the reference: muted
# slate blue below zero, a near-white centre and muted brick red above, over
# -limit to limit, with values beyond it squished to the end colours. NA
# (nothing to compare) is a mid grey, dark enough not to be read as the
# near-white of "no change".
umap_diverging_colour <- function(limit, name, breaks = ggplot2::waiver(),
                                  labels = ggplot2::waiver()) {
  if (!is.finite(limit) || limit <= 0) limit <- 1
  ggplot2::scale_colour_gradientn(
    colours = c("#3F5E8C", "#7D97BD", "#C3CEDF", "#F6F4F1",
                "#E5C4BE", "#C7837A", "#9A3F3D"),
    values = c(0, 0.2, 0.38, 0.5, 0.62, 0.8, 1),
    limits = c(-limit, limit), breaks = breaks, labels = labels,
    oob = scales::squish, na.value = "#BDBDBD", name = name,
    guide = .us_colourbar())
}
