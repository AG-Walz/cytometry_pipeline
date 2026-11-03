# =============================================================================
# visualization/umap_gate_layout.R — gate-anchored features for a UMAP in which
# every gater-defined population lands in its own island
# =============================================================================
# The plain UMAP embeds FlowJo values as they are, so a cytokine-positive cell
# sits a few units away from its negative neighbours and is drawn among them.
# Here the gates themselves are read: for every gate below the selected
# population we find the marker it decides and where its cut lies in that
# sample, record the gater's call per cell and marker, and rebuild each decision
# marker as
#
#     f = tanh((x - anchor) / (softness * sigma)) + gap_m * (call ? +1 : -1)
#
# The tanh part is bounded by 1, so with gap > 1 the two sides of every gate can
# never overlap, whatever the voltages of that visit were. Cells whose calls
# differ as positive against negative (or undecided) differ by at least
# 2 * (gap - 1) in that coordinate, so every combination of calls is its own
# region of feature space, while the bounded continuous part keeps a natural
# cloud inside each island. gap_m grows by `ratio` per tier up the tree
# (tier_gaps()), so the lineage split is the largest distance. The continuous
# part is centred on the median of the sample's own cells of
# the same combination (the anchor, combination_anchors()): that removes
# per-visit voltage shifts without depending on where the gater drew the line
# in that tube, or on which other populations share a side.
#
# Plan and rationale: docs/plans/umap_gate_layout_plan.md
#
# Everything here except gate_layout_extract_sample() is pure base R, so
# test_data/test_umap_gate_layout.R exercises it with no CytoML installed. The
# flowWorkspace touch points are confined to that one function.
#
# Conventions: calls are integer 1 (positive), 0 (negative), NA (not decided).
# A "decision" is one (gate, marker) pair with a polarity (+1 for a gate that
# keeps cells above its cut, -1 for one that keeps them below) and a cut.
# =============================================================================

.GL_RARE <- "other (rare)"

# Marker names of common viability dyes. Case-insensitive, matched against the
# marker (not the channel) name after dictionary recovery and aliasing.
VIABILITY_MARKER_PATTERN <- paste(
  "live|dead|viab|zombie|fixable|^ald$|7-?aad|dapi|^pi$|propidium|fvd|fvs",
  "sytox|draq ?7|^aqua$|^l/?-?d$|^ld[- ]|^ld$|l/d", sep = "|")

is_viability_marker <- function(markers, pattern = VIABILITY_MARKER_PATTERN) {
  if (is.null(pattern) || !nzchar(pattern)) return(rep(FALSE, length(markers)))
  !is.na(markers) & grepl(pattern, markers, ignore.case = TRUE)
}

# Parent path of a FlowJo full path. gh paths are "/A/B/C"; a gate directly
# under the root has parent "root".
.gl_parent_path <- function(paths) {
  p <- sub("/[^/]*$", "", paths)
  p[!nzchar(p)] <- "root"
  p
}

# Gates strictly below `population`, with their parent and depth relative to
# it. `parents` (named by path) overrides the string rule when the caller has
# the real tree (a gate name may itself contain "/").
gate_descendants <- function(pop_paths, population, parents = NULL) {
  pref <- if (identical(population, "root")) "/" else paste0(sub("/$", "", population), "/")
  below <- pop_paths[startsWith(pop_paths, pref)]
  if (!length(below)) {
    return(data.frame(path = character(0), parent = character(0),
                      depth = integer(0), stringsAsFactors = FALSE))
  }
  par <- if (!is.null(parents)) unname(parents[below]) else .gl_parent_path(below)
  par[is.na(par)] <- .gl_parent_path(below)[is.na(par)]
  depth <- integer(length(below))
  for (i in seq_along(below)) {
    d <- 1L; p <- par[i]
    while (!identical(p, population) && p %in% below) {
      d <- d + 1L; p <- par[match(p, below)]
    }
    depth[i] <- d
  }
  data.frame(path = below, parent = par, depth = depth,
             stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Axis rows of one gate. `x` is cells x channels (population cells, channel
# names as columns), `inside` and `in_parent` are logical vectors over the same
# cells, `geom` is list(type, dims, min, max) from gate_geometry().
# `channel_marker` maps channel -> marker (NA for scatter/Time/unnamed).
# Returns one row per marker axis the gate binds on: `type` "decision" (with a
# polarity and a cut) or "band" (both sides bind, no single cut). Axes that
# exclude nothing, scatter and unnamed channels give no row.
# -----------------------------------------------------------------------------
gate_decisions_one <- function(x, inside, in_parent, geom, channel_marker,
                               min_excluded = 0.005, min_auc = 0.8,
                               min_class = 10L, max_clip = 0.05) {
  empty <- data.frame(channel = character(0), marker = character(0),
                      polarity = integer(0), cut = numeric(0),
                      method = character(0), excluded = numeric(0),
                      type = character(0), stringsAsFactors = FALSE)
  if (is.null(geom) || identical(geom$type, "boolean") || !length(geom$dims)) {
    return(empty)
  }
  n_par <- sum(in_parent)
  if (n_par == 0L) return(empty)
  row <- function(ch, mk, pol, cut, method, ex, type) {
    data.frame(channel = ch, marker = mk, polarity = pol, cut = cut,
               method = method, excluded = ex, type = type,
               stringsAsFactors = FALSE)
  }
  rows <- list()
  for (ch in geom$dims) {
    mk <- unname(channel_marker[ch])
    if (is.na(mk) || !nzchar(mk) || !(ch %in% colnames(x))) next
    v <- x[in_parent, ch]
    if (identical(geom$type, "rectangle")) {
      lo <- unname(geom$min[ch]); hi <- unname(geom$max[ch])
      ex_lo <- if (is.finite(lo)) mean(v < lo) else 0
      ex_hi <- if (is.finite(hi)) mean(v > hi) else 0
      b_lo <- ex_lo >= min_excluded
      b_hi <- ex_hi >= min_excluded
      # A far edge that clips only a small tail (a box not dragged to the end
      # of the axis) does not make the gate a band: the other edge is the
      # decision. A real band (a "CD4 dim" gate) excludes a share on both
      # sides.
      if (b_lo && b_hi && min(ex_lo, ex_hi) < max_clip) {
        if (ex_lo >= ex_hi) b_hi <- FALSE else b_lo <- FALSE
      }
      if (b_lo && !b_hi) {
        rows[[length(rows) + 1L]] <- row(ch, mk, 1L, lo, "rectangle", ex_lo,
                                         "decision")
      } else if (b_hi && !b_lo) {
        rows[[length(rows) + 1L]] <- row(ch, mk, -1L, hi, "rectangle", ex_hi,
                                         "decision")
      } else if (b_lo && b_hi) {
        rows[[length(rows) + 1L]] <- row(ch, mk, NA_integer_, NA_real_,
                                         "rectangle", ex_lo + ex_hi, "band")
      }
    } else {
      mem <- inside[in_parent]
      n1 <- sum(mem); n0 <- length(mem) - n1
      if (n1 < min_class || n0 < min_class) next
      auc <- .gl_auc(v[mem], v[!mem])
      tol <- 1e-9   # 1 - 0.8 is not exactly 0.2 in floating point
      if (!is.finite(auc) || (auc < min_auc - tol && auc > 1 - min_auc + tol)) next
      pol <- if (auc >= 0.5) 1L else -1L
      rows[[length(rows) + 1L]] <- row(ch, mk, pol, .gl_youden_cut(v, mem, pol),
                                       "membership", mean(!mem), "decision")
    }
  }
  if (!length(rows)) return(empty)
  do.call(rbind, rows)
}

# Mann-Whitney AUC: P(member > non-member), ties counted half.
.gl_auc <- function(a, b) {
  n1 <- length(a); n0 <- length(b)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(c(a, b))
  # Doubles, not integers: n1 * n0 passes 2^31 at about 46000 x 46000 cells.
  n1 <- as.numeric(n1); n0 <- as.numeric(n0)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# The 1-D cut that best reproduces membership (maximum Youden J). For a
# negative-polarity gate the members are BELOW the cut.
.gl_youden_cut <- function(v, mem, polarity) {
  if (polarity < 0) v <- -v
  o <- order(v); vs <- v[o]; ms <- mem[o]
  n1 <- sum(ms); n0 <- length(ms) - n1
  # Cut between position k and k+1: members above it are TP, non-members above
  # it are FP. Evaluated at every split in one pass.
  tp <- n1 - cumsum(ms); fp <- n0 - cumsum(!ms)
  j <- tp / n1 - fp / n0
  k <- which.max(j)
  cut <- if (k < length(vs)) (vs[k] + vs[k + 1L]) / 2 else vs[k]
  if (polarity < 0) -cut else cut
}

# -----------------------------------------------------------------------------
# All decisions of one sample. `gates` is the gate_descendants() frame with a
# `geom` list column; `mem` is cells x gates (population cells) named by path.
#
# Two rules keep a gate from deciding what it does not decide:
#   - A gate never re-decides a marker that one of its ANCESTOR gates already
#     decides. A cytokine box drawn on CD4 x IFNg under "CD4+" may clip the
#     CD4-dim tail; that edge trims cells inside a side the CD4+ gate already
#     chose, it is not a second CD4 threshold. (Sibling subtrees are not
#     ancestors: an IFNg gate under "Viable" and another under "CD4+" are both
#     kept, and gate_marker_calls() lets the deeper one speak for CD4+ cells.)
#   - `n_axes` counts the marker axes a gate binds on after that filter
#     (decisions and bands). gate_marker_calls() only infers a call for cells
#     OUTSIDE a gate from its polarity when the gate binds on that marker
#     alone; a cell can be outside a CD4+CD8- box because of its CD8.
# Returns a frame with gate, parent, depth, channel, marker, polarity, cut,
# method, excluded, n_parent, n_member, n_axes (decision rows only).
# -----------------------------------------------------------------------------
gate_decisions_sample <- function(x, mem, gates, channel_marker, population,
                                  min_excluded = 0.005, min_auc = 0.8,
                                  min_class = 10L) {
  cols <- c("gate", "parent", "depth", "channel", "marker", "polarity", "cut",
            "method", "excluded", "n_parent", "n_member", "n_axes")
  out <- list()
  decided <- list()            # gate path -> markers it decides (after filter)
  parent_of <- stats::setNames(gates$parent, gates$path)
  for (i in order(gates$depth)) {
    g <- gates$path[i]
    if (!(g %in% colnames(mem))) next
    in_parent <- gate_parent_membership(mem, gates$parent[i], population)
    d <- gate_decisions_one(x, mem[, g], in_parent, gates$geom[[i]],
                            channel_marker, min_excluded, min_auc, min_class)
    if (!nrow(d)) next
    anc_markers <- character(0)
    p <- gates$parent[i]
    while (!identical(p, population) && p %in% names(parent_of)) {
      anc_markers <- c(anc_markers, decided[[p]])
      p <- parent_of[[p]]
    }
    d <- d[!(d$marker %in% anc_markers), , drop = FALSE]
    if (!nrow(d)) next
    d$n_axes <- nrow(d)
    d <- d[d$type == "decision", , drop = FALSE]
    if (!nrow(d)) next
    decided[[g]] <- unique(d$marker)
    d$gate <- g; d$parent <- gates$parent[i]; d$depth <- gates$depth[i]
    d$n_parent <- sum(in_parent); d$n_member <- sum(mem[, g] & in_parent)
    out[[length(out) + 1L]] <- d[, cols]
  }
  if (!length(out)) {
    return(data.frame(gate = character(0), parent = character(0),
                      depth = integer(0), channel = character(0),
                      marker = character(0), polarity = integer(0),
                      cut = numeric(0), method = character(0),
                      excluded = numeric(0), n_parent = integer(0),
                      n_member = integer(0), n_axes = integer(0),
                      stringsAsFactors = FALSE))
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

gate_parent_membership <- function(mem, parent, population) {
  if (identical(parent, population)) return(rep(TRUE, nrow(mem)))
  if (parent %in% colnames(mem)) return(mem[, parent])
  # A parent outside the population's subtree cannot happen for a descendant;
  # treat it as "everyone" rather than dropping the decision.
  rep(TRUE, nrow(mem))
}

# -----------------------------------------------------------------------------
# Per-cell calls and cuts for every decided marker of one sample.
# `x` cells x channels, `mem` cells x gates, `dec` from gate_decisions_sample().
# Returns list(calls  = integer matrix cells x markers (1/0/NA),
#              cuts   = numeric matrix cells x markers (the gate cut that
#                       decided the cell; reported, not used for features),
#              values = numeric matrix cells x markers (the raw value of the
#                       channel each marker was decided on)).
# -----------------------------------------------------------------------------
gate_marker_calls <- function(x, mem, dec, population) {
  markers <- unique(dec$marker)
  n <- nrow(x)
  calls   <- matrix(NA_integer_, n, length(markers), dimnames = list(NULL, markers))
  cuts    <- matrix(NA_real_,    n, length(markers), dimnames = list(NULL, markers))
  values  <- matrix(NA_real_,    n, length(markers), dimnames = list(NULL, markers))
  n_axes_all <- if ("n_axes" %in% names(dec)) dec$n_axes else rep(1L, nrow(dec))
  for (m in markers) {
    keep_m <- dec$marker == m
    d <- dec[keep_m, , drop = FALSE]
    d$n_axes <- n_axes_all[keep_m]
    # One channel per marker: the most common one among its decisions (a
    # marker read from two detectors in one sample would be a panel error).
    ch <- names(sort(table(d$channel), decreasing = TRUE))[1]
    d <- d[d$channel == ch, , drop = FALSE]
    v <- x[, ch]
    values[, m] <- v
    k <- nrow(d)
    app <- matrix(FALSE, n, k)
    for (r in seq_len(k)) {
      app[, r] <- gate_parent_membership(mem, d$parent[r], population)
    }

    # Cells share a handful of applicability patterns (which decisions' parents
    # contain them), so everything is resolved per pattern. Within a pattern
    # only the DEEPEST applicable decisions count: when a marker is gated both
    # under "Viable" and again under "CD4+", the gater's CD4+-level gate is the
    # decision for a CD4+ cell, and the shallower gate only speaks for cells no
    # deeper gate covers. Same-depth decisions (the four quadrant gates) are
    # all used together.
    pat <- if (k == 1L) as.character(app[, 1]) else
      do.call(paste0, lapply(seq_len(k), function(r) as.integer(app[, r])))
    cut_cell <- rep(stats::median(d$cut), n)
    cl <- rep(NA_integer_, n)
    for (p in unique(pat)) {
      sel <- which(pat == p)
      a <- app[sel[1], ]
      if (!any(a)) next                        # no decision covers these cells
      rows <- which(a)
      rows <- rows[d$depth[rows] == max(d$depth[rows])]
      cut_cell[sel] <- stats::median(d$cut[rows])
      ins <- mem[sel, d$gate[rows], drop = FALSE]
      pos <- d$polarity[rows] > 0
      any_in <- rowSums(ins) > 0
      pos_in <- if (any(pos)) rowSums(ins[, pos, drop = FALSE]) > 0 else rep(FALSE, length(sel))
      c_p <- rep(NA_integer_, length(sel))
      c_p[any_in] <- as.integer(pos_in[any_in])
      # Outside every deciding gate. When each of them binds on this marker
      # alone, being outside IS the decision: negative when every one keeps
      # the cells above its cut (outside the IFNg+ gate, whatever its FSC-A
      # edge did), positive when every one keeps them below. When a gate also
      # binds on another marker, a cell may be outside it because of that
      # other axis, so the cell's own value against the cut decides.
      out <- !any_in
      if (any(out)) {
        single <- all(d$n_axes[rows] == 1L)
        if (single && all(pos))       c_p[out] <- 0L
        else if (single && !any(pos)) c_p[out] <- 1L
        else c_p[out] <- as.integer(v[sel][out] >= cut_cell[sel][out])
      }
      cl[sel] <- c_p
    }
    cuts[, m] <- cut_cell
    calls[, m] <- cl
  }
  list(calls = calls, cuts = cuts, values = values)
}

# -----------------------------------------------------------------------------
# Anchors: per sample, the median of every marker over the cells that share a
# combination of calls. The continuous part of a feature is measured from it
# (gate_features()), so a cell sits inside its island by where it lies among
# its own sample's cells OF THAT COMBINATION. Two things that differ between
# tubes therefore cannot move an island's cells apart: where the gater drew
# the line (a cut-centred feature splits identical cells by tube), and which
# other populations share a side (CD8 T cells and double negatives are both
# CD4-negative, in proportions that differ per tube, so a side median puts CD8
# T cells of different tubes at different CD4 positions). Voltage shifts are
# removed along the way. What it gives up is between-tube location within an
# island, which is mostly instrument setting here; the marker plots keep the
# raw values.
# -----------------------------------------------------------------------------
combination_anchors <- function(values, calls) {
  anchors <- values
  if (!ncol(values)) return(anchors)
  key <- calls_key(calls)
  for (m in colnames(values)) {
    anchors[, m] <- stats::ave(values[, m], key, FUN = stats::median)
  }
  anchors
}

# -----------------------------------------------------------------------------
# Stratification key of a calls matrix: one string per cell, equal for cells
# with identical calls. NA is spelled "x". Column order is the caller's.
# -----------------------------------------------------------------------------
calls_key <- function(calls) {
  if (!ncol(calls)) return(rep("", nrow(calls)))
  parts <- lapply(seq_len(ncol(calls)), function(j) {
    v <- as.character(calls[, j]); v[is.na(v)] <- "x"; v
  })
  do.call(paste0, parts)
}

# -----------------------------------------------------------------------------
# Keep up to `cap` cells per key; rare keys are kept whole. Returns
# list(idx = sorted row indices, weight = cells of that key / cells kept).
# cap <= 0 keeps everything with weight 1.
# -----------------------------------------------------------------------------
stratified_sample <- function(key, cap) {
  n <- length(key)
  if (cap <= 0 || n == 0L) return(list(idx = seq_len(n), weight = rep(1, n)))
  by_key <- split(seq_len(n), key)
  idx <- integer(0); w <- numeric(0)
  # Iterated by position, not looked up by name: a sample with no decision has
  # the key "" for every cell, and by_key[[""]] is NULL in R.
  for (rows in by_key) {
    take <- if (length(rows) > cap) rows[sample.int(length(rows), cap)] else rows
    idx <- c(idx, take)
    w   <- c(w, rep(length(rows) / length(take), length(take)))
  }
  o <- order(idx)
  list(idx = idx[o], weight = w[o])
}

# -----------------------------------------------------------------------------
# Row indices of a fit set in which no label has more than budget / (number of
# labels) cells; smaller labels are kept whole and what they leave unused is
# NOT handed to the large ones. Island area in a UMAP follows the number of
# fitted cells, so a fixed per-label cap is what gives the abundant and the
# rare combinations comparable room, and it does so whether or not the run has
# more cells than the budget (a small run is balanced too). The set holds at
# most `budget` cells, or one per label when the budget is smaller than that.
# -----------------------------------------------------------------------------
balanced_fit_index <- function(labels, budget) {
  n <- length(labels)
  if (budget <= 0 || n == 0L) return(seq_len(n))
  by <- split(seq_len(n), labels)
  cap <- max(1L, floor(budget / length(by)))
  idx <- unlist(lapply(by, function(rows) {
    if (length(rows) > cap) rows[sample.int(length(rows), cap)] else rows
  }), use.names = FALSE)
  sort(idx)
}

# -----------------------------------------------------------------------------
# Marker tiers and order from the pooled decisions: a marker's tier is the
# (rounded median) depth below the population of the gates that decide it, and
# markers within a tier follow `panel_order` (first appearance otherwise).
# -----------------------------------------------------------------------------
marker_tiers <- function(decisions, panel_order = NULL) {
  if (!nrow(decisions)) return(integer(0))
  t <- tapply(decisions$depth, decisions$marker,
              function(d) as.integer(round(stats::median(d))))
  mk <- names(t)
  ord_ref <- c(panel_order[panel_order %in% mk],
               setdiff(unique(decisions$marker), panel_order))
  mk <- mk[order(as.integer(t[mk]), match(mk, ord_ref))]
  stats::setNames(as.integer(t[mk]), mk)
}

# -----------------------------------------------------------------------------
# Human-readable combination label per cell from a calls matrix whose columns
# include names(tiers). Per tier: the positive markers ("IFNy+ TNFa+"); none
# positive reads "none", except the first tier, which spells its markers as
# negative ("CD4- CD8-"). Tiers in which every marker is NA for a cell are left
# out. Tiers are joined with " | ". Computed once per distinct call pattern.
# -----------------------------------------------------------------------------
combination_labels <- function(calls, tiers) {
  mk <- names(tiers)
  calls <- calls[, mk, drop = FALSE]
  key <- calls_key(calls)
  first <- match(unique(key), key)
  tier_levels <- sort(unique(tiers))
  lab_u <- vapply(first, function(i) {
    cl <- calls[i, ]
    parts <- character(0)
    for (t in tier_levels) {
      m_t <- mk[tiers == t]
      c_t <- cl[m_t]
      if (all(is.na(c_t))) next
      pos <- m_t[!is.na(c_t) & c_t == 1L]
      parts <- c(parts,
                 if (length(pos)) paste0(pos, "+", collapse = " ")
                 else if (t == tier_levels[1]) paste0(m_t[!is.na(c_t)], "-", collapse = " ")
                 else "none")
    }
    if (!length(parts)) "(no decision)" else paste(parts, collapse = " | ")
  }, character(1))
  lab_u[match(key, unique(key))]
}

# Labels with fewer than `min_cells` (weighted) cells become .GL_RARE.
pool_rare_labels <- function(labels, weights = NULL, min_cells = 0) {
  if (min_cells <= 0) return(labels)
  if (is.null(weights)) weights <- rep(1, length(labels))
  tot <- tapply(weights, labels, sum)
  rare <- names(tot)[tot < min_cells]
  labels[labels %in% rare] <- .GL_RARE
  labels
}

# -----------------------------------------------------------------------------
# Fill decisions a sample is missing. Within each first-tier group (e.g. all
# "CD4+" cells), a deeper marker decided for most of the group's cells is
# treated as negative where a sample has no decision for it, and reported.
# Returns list(calls, filled = data.frame(sample, group, marker, cells)).
# -----------------------------------------------------------------------------
fill_missing_decisions <- function(calls, tiers, sample, majority = 0.5,
                                   weights = NULL) {
  # The majority is a share of the group's CELLS, so the sampling weights
  # count: a share of kept cells would count tubes, and flip with the cap.
  if (is.null(weights)) weights <- rep(1, nrow(calls))
  filled <- list()
  if (length(unique(tiers)) < 2L) return(list(calls = calls, filled = NULL))
  t1 <- names(tiers)[tiers == min(tiers)]
  deeper <- names(tiers)[tiers > min(tiers)]
  grp <- calls_key(calls[, t1, drop = FALSE])
  for (g in unique(grp)) {
    rows <- which(grp == g)
    for (m in deeper) {
      na <- is.na(calls[rows, m])
      if (!any(na) || sum(weights[rows][!na]) / sum(weights[rows]) <= majority) next
      hit <- rows[na]
      for (s in unique(sample[hit])) {
        filled[[length(filled) + 1L]] <- data.frame(
          sample = s, group = g, marker = m, cells = sum(sample[hit] == s),
          stringsAsFactors = FALSE)
      }
      calls[hit, m] <- 0L
    }
  }
  list(calls = calls,
       filled = if (length(filled)) do.call(rbind, filled) else NULL)
}

# -----------------------------------------------------------------------------
# Pooled scale per marker: weighted MAD of (value - anchor) over negative cells,
# floored at `floor_frac` of the channel range (FlowJo channels span 0-256).
# -----------------------------------------------------------------------------
gate_scales <- function(values, anchors, calls, weights = NULL,
                        channel_range = 256, floor_frac = 0.01) {
  if (is.null(weights)) weights <- rep(1, nrow(values))
  s <- vapply(colnames(values), function(m) {
    neg <- !is.na(calls[, m]) & calls[, m] == 0L & is.finite(values[, m])
    if (sum(neg) < 10L) neg <- is.finite(values[, m])
    d <- values[neg, m] - anchors[neg, m]
    w <- weights[neg]
    med <- .gl_wquantile(d, w, 0.5)
    1.4826 * .gl_wquantile(abs(d - med), w, 0.5)
  }, numeric(1))
  pmax(s, floor_frac * channel_range)
}

.gl_wquantile <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= p)[1]]
}

# -----------------------------------------------------------------------------
# One jump per marker, larger for shallower tiers: a marker of tier t gets
# gap * ratio^(deepest tier - t). With ratio > 1 the first split of the gating
# tree (CD4 vs CD8) is the largest distance in feature space, so the islands of
# one lineage sit together and the cytokine combinations arrange themselves
# inside that region instead of being scattered across the whole embedding.
# -----------------------------------------------------------------------------
tier_gaps <- function(tiers, gap, ratio = 1) {
  if (!length(tiers)) return(numeric(0))
  stats::setNames(gap * ratio^(max(tiers) - tiers), names(tiers))
}

# -----------------------------------------------------------------------------
# The gate-anchored feature matrix: tanh((value - anchor) / (softness *
# sigma)) plus the marker's jump on the side of its call. NA calls are placed
# on the negative side. Extra (undecided) markers in `extra` enter as
# tanh(robust z / softness), without a jump.
# -----------------------------------------------------------------------------
gate_features <- function(values, anchors, calls, sigma, gap = 2, softness = 2,
                          extra = NULL) {
  # `gap` is one number, or one per marker (named), e.g. from tier_gaps().
  gap_of <- function(m) if (length(gap) > 1L) gap[[m]] else gap
  f <- values
  for (m in colnames(values)) {
    z <- (values[, m] - anchors[, m]) / (softness * sigma[[m]])
    side <- ifelse(!is.na(calls[, m]) & calls[, m] == 1L, 1, -1)
    f[, m] <- tanh(z) + gap_of(m) * side
  }
  if (!is.null(extra) && ncol(extra)) {
    ex <- extra
    for (m in colnames(extra)) {
      v <- extra[, m]
      med <- stats::median(v); s <- stats::mad(v)
      if (!is.finite(s) || s <= 0) s <- stats::sd(v)
      if (!is.finite(s) || s <= 0) s <- 1
      ex[, m] <- tanh((v - med) / (softness * s))
    }
    f <- cbind(f, ex)
  }
  f
}

# -----------------------------------------------------------------------------
# Sample selection and level ordering from metadata.
# -----------------------------------------------------------------------------

# "col=a,b;col2!=c" -> list of list(column, negate, values).
parse_sample_filter <- function(spec) {
  if (is.null(spec) || !nzchar(trimws(spec))) return(list())
  clauses <- trimws(strsplit(spec, ";", fixed = TRUE)[[1]])
  clauses <- clauses[nzchar(clauses)]
  lapply(clauses, function(cl) {
    m <- regmatches(cl, regexec("^([^!=]+?)\\s*(!=|=)\\s*(.*)$", cl))[[1]]
    if (length(m) != 4L) {
      stop("--filter clause '", cl, "' is not of the form column=value[,value] ",
           "or column!=value[,value].", call. = FALSE)
    }
    vals <- trimws(strsplit(m[4], ",", fixed = TRUE)[[1]])
    list(column = trimws(m[2]), negate = identical(m[3], "!="),
         values = vals[nzchar(vals)])
  })
}

# Logical keep-vector over the rows of `annotations` (one row per sample).
apply_sample_filter <- function(annotations, clauses) {
  keep <- rep(TRUE, nrow(annotations))
  for (cl in clauses) {
    if (!(cl$column %in% names(annotations))) {
      stop("--filter names column '", cl$column, "', which the metadata does ",
           "not have. Columns: ", paste(names(annotations), collapse = ", "),
           call. = FALSE)
    }
    hit <- annotations[[cl$column]] %in% cl$values
    keep <- keep & if (cl$negate) !hit else hit
  }
  keep
}

# Natural order: digit runs compare as numbers ("W5" < "W12" < "W110").
natural_order <- function(x) {
  x <- as.character(x)
  if (!length(x)) return(integer(0))
  toks <- regmatches(x, gregexpr("[0-9]+|[^0-9]+", x))
  width <- max(1L, vapply(toks, length, integer(1)))
  keys <- lapply(seq_len(width), function(j) {
    t <- vapply(toks, function(v) if (length(v) >= j) v[j] else "", character(1))
    num <- grepl("^[0-9]+$", t)
    list(num = ifelse(num, suppressWarnings(as.numeric(t)), -Inf),
         chr = ifelse(num, "", tolower(t)))
  })
  args <- unlist(lapply(keys, function(k) list(k$chr, k$num)), recursive = FALSE)
  do.call(order, c(args, list(method = "radix")))
}

# "col=lvA,lvB;col2=..." -> named list of character vectors.
parse_level_order <- function(spec) {
  if (is.null(spec) || !nzchar(trimws(spec))) return(list())
  parts <- trimws(strsplit(spec, ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  out <- list()
  for (p in parts) {
    kv <- regmatches(p, regexpr("=", p), invert = TRUE)[[1]]
    if (length(kv) != 2L) {
      stop("--level_order entry '", p, "' is not of the form column=a,b,c.",
           call. = FALSE)
    }
    v <- trimws(strsplit(kv[2], ",", fixed = TRUE)[[1]])
    out[[trimws(kv[1])]] <- v[nzchar(v)]
  }
  out
}

# Levels of `values` in display order: an explicit order first (levels it does
# not mention follow in natural order), "Unassigned" never included.
ordered_levels <- function(values, explicit = NULL, drop = "Unassigned") {
  lv <- setdiff(unique(as.character(values)), drop)
  lv <- lv[!is.na(lv)]
  lv <- lv[natural_order(lv)]
  if (length(explicit)) {
    unknown <- setdiff(explicit, lv)
    if (length(unknown)) {
      warning("level order names level(s) not present: ",
              paste(unknown, collapse = ", "), call. = FALSE)
    }
    lv <- c(explicit[explicit %in% lv], setdiff(lv, explicit))
  }
  lv
}

# -----------------------------------------------------------------------------
# Separation metrics of a 2-D embedding against labels. Uses an exact kNN from
# RANN (or FNN); returns NULL with a message when neither is installed.
# `groups` (e.g. patient) adds a mixing score per label: the mean, over a
# label's cells, of the share of its neighbours that come from a DIFFERENT
# group, divided by what random mixing within the label would give (1 = as
# mixed as chance, 0 = every neighbour from the cell's own group).
# -----------------------------------------------------------------------------
separation_metrics <- function(coords, labels, groups = NULL, k = 15L,
                               max_cells = 50000L, seed = 42L,
                               samples = NULL) {
  knn <- .gl_knn_fun()
  if (is.null(knn)) {
    message("separation metrics skipped: install RANN or FNN.")
    return(NULL)
  }
  set.seed(seed)
  labels <- as.character(labels)
  idx <- balanced_fit_index(labels, max_cells)
  xy <- as.matrix(coords[idx, 1:2]); lab <- labels[idx]
  grp <- if (!is.null(groups)) as.character(groups[idx]) else NULL
  smp <- if (!is.null(samples)) as.character(samples[idx]) else NULL
  nn <- knn(xy, k + 1L)[, -1, drop = FALSE]
  same <- matrix(lab[nn] == lab, ncol = k)
  purity <- rowMeans(same)
  isolated <- purity == 1
  # Simplified silhouette: own centroid vs nearest other centroid.
  cent <- do.call(rbind, lapply(split(seq_along(lab), lab), function(r)
    colMeans(xy[r, , drop = FALSE])))
  dmat <- sapply(seq_len(nrow(cent)), function(j)
    sqrt((xy[, 1] - cent[j, 1])^2 + (xy[, 2] - cent[j, 2])^2))
  if (is.null(dim(dmat))) dmat <- matrix(dmat, ncol = nrow(cent))
  own <- match(lab, rownames(cent))
  a <- dmat[cbind(seq_along(lab), own)]
  dmat[cbind(seq_along(lab), own)] <- Inf
  b <- apply(dmat, 1, min)
  sil <- ifelse(is.finite(b), (b - a) / pmax(a, b), NA_real_)
  # Mixing within a label, relative to chance: the share of same-label
  # neighbours from another group over what random placement among that
  # label's groups would give. Scored for workspaces (`groups`) and for
  # samples (`samples`): a low sample mixing in one island means its samples
  # sit apart inside it, which is how a cut placed differently per tube, or a
  # real per-sample shift, shows up.
  mixing <- function(g) {
    mix <- rep(NA_real_, length(lab))
    if (is.null(g)) return(mix)
    other <- matrix(g[nn] != g, ncol = k)
    obs <- rowSums(same & other) / pmax(rowSums(same), 1)
    for (l in unique(lab)) {
      r <- which(lab == l)
      sh <- table(g[r]) / length(r)
      expct <- 1 - sum(sh^2)
      mix[r] <- if (expct > 0) obs[r] / expct else NA_real_
    }
    mix
  }
  mix <- mixing(grp)
  smix <- mixing(smp)
  per <- do.call(rbind, lapply(split(seq_along(lab), lab), function(r)
    data.frame(label = lab[r[1]], cells_scored = length(r),
               knn_purity = mean(purity[r]), isolated = mean(isolated[r]),
               silhouette = mean(sil[r], na.rm = TRUE),
               group_mixing = if (is.null(grp)) NA_real_ else mean(mix[r], na.rm = TRUE),
               sample_mixing = if (is.null(smp)) NA_real_ else mean(smix[r], na.rm = TRUE),
               stringsAsFactors = FALSE)))
  overall <- data.frame(label = "(macro average)", cells_scored = sum(per$cells_scored),
                        knn_purity = mean(per$knn_purity),
                        isolated = mean(per$isolated),
                        silhouette = mean(per$silhouette, na.rm = TRUE),
                        group_mixing = if (is.null(grp)) NA_real_ else mean(per$group_mixing, na.rm = TRUE),
                        sample_mixing = if (is.null(smp)) NA_real_ else mean(per$sample_mixing, na.rm = TRUE),
                        stringsAsFactors = FALSE)
  rbind(overall, per[order(-per$cells_scored), , drop = FALSE])
}

.gl_knn_fun <- function() {
  if (requireNamespace("RANN", quietly = TRUE)) {
    return(function(x, k) RANN::nn2(x, k = k)$nn.idx)
  }
  if (requireNamespace("FNN", quietly = TRUE)) {
    return(function(x, k) cbind(seq_len(nrow(x)),
                                FNN::get.knn(x, k = k - 1L)$nn.index))
  }
  NULL
}

# -----------------------------------------------------------------------------
# Weighted 2-D density on a grid: weighted binning, then a separable Gaussian
# kernel. `bandwidth` is a fraction of each axis range. Returns
# list(x, y, z) with z normalised to sum to 1 (or all zero when w sums to 0).
# -----------------------------------------------------------------------------
weighted_density_grid <- function(x, y, w = NULL, xlim, ylim, n = 256L,
                                  bandwidth = 0.012) {
  if (is.null(w)) w <- rep(1, length(x))
  gx <- seq(xlim[1], xlim[2], length.out = n)
  gy <- seq(ylim[1], ylim[2], length.out = n)
  ix <- pmin(pmax(round((x - xlim[1]) / diff(xlim) * (n - 1)) + 1, 1), n)
  iy <- pmin(pmax(round((y - ylim[1]) / diff(ylim) * (n - 1)) + 1, 1), n)
  h <- matrix(0, n, n)
  agg <- tapply(w, list(factor(ix, levels = seq_len(n)),
                        factor(iy, levels = seq_len(n))), sum)
  agg[is.na(agg)] <- 0
  h[] <- agg
  kern <- function(len, bw) {
    s <- bw * (len - 1)
    d <- outer(seq_len(len), seq_len(len), "-")
    k <- exp(-0.5 * (d / max(s, 1e-9))^2)
    k / rowSums(k)
  }
  z <- kern(n, bandwidth) %*% h %*% t(kern(n, bandwidth))
  tot <- sum(z)
  list(x = gx, y = gy, z = if (tot > 0) z / tot else z)
}

# -----------------------------------------------------------------------------
# CytoML-facing extraction of one sample. The only function in this file that
# touches flowWorkspace. Returns NULL when the sample has no named marker
# channel, else list(
#   mtx       = cells x markers raw values for every named marker channel
#               (same contract as extract_population_sample()),
#   decisions = gate_decisions_sample() frame,
#   calls, cuts, values = per decision marker (gate_marker_calls()),
#   key       = stratification key per cell).
# -----------------------------------------------------------------------------
gate_layout_extract_sample <- function(gh, population, dictionary = NULL,
                                       aliases = NULL,
                                       exclude_pattern = "Time|FSC|SSC|scatter",
                                       exclude_markers = VIABILITY_MARKER_PATTERN,
                                       min_excluded = 0.005) {
  ff   <- flowWorkspace::gh_pop_get_data(gh, y = population)
  expr <- flowCore::exprs(ff)
  sel  <- marker_columns(ff, exclude_pattern, dictionary = dictionary,
                         aliases = aliases)
  if (!any(sel$keep)) return(NULL)
  chans <- colnames(expr)
  # A viability gate is a gate like any other, so a population chosen ABOVE the
  # live gate would make the live/dead dye a decision marker and give dead cells
  # their own island. Nobody wants that axis in this layout; marker names that
  # look like a viability dye are never decisions (the raw values still reach
  # the marker plots through mtx).
  decidable <- sel$keep & !is_viability_marker(sel$markers, exclude_markers)
  channel_marker <- stats::setNames(ifelse(decidable, sel$markers, NA_character_),
                                    chans)

  mtx <- expr[, sel$keep, drop = FALSE]
  colnames(mtx) <- sel$markers[sel$keep]
  mtx <- mtx[, !duplicated(colnames(mtx)), drop = FALSE]

  pop_idx <- flowWorkspace::gh_pop_get_indices(gh, population)
  stopifnot(sum(pop_idx) == nrow(expr))

  all_paths <- flowWorkspace::gh_get_pop_paths(gh, path = "full")
  pref <- if (identical(population, "root")) "/" else paste0(sub("/$", "", population), "/")
  below <- all_paths[startsWith(all_paths, pref)]
  parents <- stats::setNames(vapply(below, function(p)
    flowWorkspace::gh_pop_get_parent(gh, p, path = "full"), character(1)), below)
  gates <- gate_descendants(all_paths, population, parents)
  gates$geom <- lapply(gates$path, function(p) gate_geometry(gh, p))

  mem <- matrix(FALSE, nrow(expr), 0)
  if (nrow(gates)) {
    m <- flowWorkspace::gh_pop_get_indices_mat(gh, y = gates$path)
    if (is.null(dim(m))) m <- matrix(m, ncol = nrow(gates))
    mem <- m[pop_idx, , drop = FALSE]
    colnames(mem) <- gates$path
  }

  dec <- gate_decisions_sample(expr, mem, gates, channel_marker, population,
                               min_excluded = min_excluded)
  cm <- if (nrow(dec)) gate_marker_calls(expr, mem, dec, population) else
    list(calls = matrix(NA_integer_, nrow(expr), 0),
         cuts = matrix(NA_real_, nrow(expr), 0),
         values = matrix(NA_real_, nrow(expr), 0))
  cm$calls <- cm$calls[, sort(colnames(cm$calls)), drop = FALSE]
  cm$cuts <- cm$cuts[, colnames(cm$calls), drop = FALSE]
  cm$values <- cm$values[, colnames(cm$calls), drop = FALSE]
  cm$anchors <- combination_anchors(cm$values, cm$calls)
  list(mtx = mtx, decisions = dec, calls = cm$calls, cuts = cm$cuts,
       anchors = cm$anchors, values = cm$values, key = calls_key(cm$calls))
}

# Geometry of one gate in the form gate_decisions_one() reads.
gate_geometry <- function(gh, path) {
  g <- tryCatch(flowWorkspace::gh_pop_get_gate(gh, path), error = function(e) NULL)
  if (is.null(g) || methods::is(g, "booleanFilter")) {
    return(list(type = "boolean", dims = character(0)))
  }
  if (methods::is(g, "rectangleGate")) {
    return(list(type = "rectangle", dims = names(g@min),
                min = g@min, max = g@max))
  }
  dims <- tryCatch(as.character(flowCore::parameters(g)),
                   error = function(e) character(0))
  list(type = tolower(sub("Gate$", "", class(g)[1])), dims = dims)
}

# -----------------------------------------------------------------------------
# Pool the per-sample gate-layout extractions of a run into one set of cells.
#
# `per_sample` is a named list (sample label -> list(calls, cuts, values, mtx,
# weight, key, key_counts, population_cells)), every matrix and vector
# row-aligned to the cells KEPT for that sample; `key_counts` is the named
# count of every stratification key over ALL of the sample's cells, taken
# before sampling. `decisions` holds every sample's decisions with a `sample`
# column. Rows of the result follow the order of `per_sample`.
#
# Labels are computed once, on the pooled calls, because tiers and the
# missing-decision fill are properties of the cohort, not of one sample. The
# counts are then read off the labels: cells sharing a sample and a key have
# identical calls, so they receive identical labels, and every key a sample
# has is represented among its kept cells (stratified_sample keeps at least
# one cell per key). That is what lets the counts cover every cell while the
# labels are computed on the kept ones only.
#
# A decision marker that a sample has no decision for is NA in that sample's
# calls (fill_missing_decisions() may then set it to negative); its values are
# taken from the sample's own channel when it has one and its cut is the
# pooled median cut of that marker, so the feature stays finite.
# -----------------------------------------------------------------------------
pool_gate_layout <- function(per_sample, decisions, panel_order = NULL,
                             min_combination_cells = 50) {
  samples <- names(per_sample)
  tiers <- marker_tiers(decisions, panel_order)
  mk <- names(tiers)
  pooled_cut <- if (nrow(decisions)) {
    tapply(decisions$cut, decisions$marker, stats::median)
  } else numeric(0)

  parts <- lapply(samples, function(s) {
    p <- per_sample[[s]]
    n <- length(p$weight)
    calls  <- matrix(NA_integer_, n, length(mk), dimnames = list(NULL, mk))
    cuts   <- matrix(NA_real_,    n, length(mk), dimnames = list(NULL, mk))
    anchors <- matrix(NA_real_,   n, length(mk), dimnames = list(NULL, mk))
    values <- matrix(NA_real_,    n, length(mk), dimnames = list(NULL, mk))
    have <- intersect(mk, colnames(p$calls))
    calls[, have]  <- p$calls[, have, drop = FALSE]
    cuts[, have]   <- p$cuts[, have, drop = FALSE]
    anchors[, have] <- p$anchors[, have, drop = FALSE]
    values[, have] <- p$values[, have, drop = FALSE]
    for (m in setdiff(mk, have)) {
      cuts[, m] <- pooled_cut[[m]]
      if (m %in% colnames(p$mtx)) values[, m] <- p$mtx[, m]
      # A marker this sample never decided: centred per combination like the
      # others (sampling is uniform within a combination, so the median of the
      # kept cells is the combination's median).
      vv <- values[, m]
      anchors[, m] <- if (all(is.finite(vv))) {
        stats::ave(vv, p$key, FUN = stats::median)
      } else pooled_cut[[m]]
    }
    list(calls = calls, cuts = cuts, anchors = anchors, values = values)
  })
  stack <- function(what) {
    out <- do.call(rbind, lapply(parts, `[[`, what))
    if (is.null(out)) out <- matrix(if (what == "calls") NA_integer_ else NA_real_,
                                    0, length(mk), dimnames = list(NULL, mk))
    out
  }
  calls  <- stack("calls")
  cuts   <- stack("cuts")
  anchors <- stack("anchors")
  values <- stack("values")
  n_per  <- vapply(per_sample, function(p) length(p$weight), integer(1))
  sample <- rep(samples, n_per)
  weight <- unlist(lapply(per_sample, `[[`, "weight"), use.names = FALSE)
  key    <- unlist(lapply(per_sample, `[[`, "key"), use.names = FALSE)

  fill <- fill_missing_decisions(calls, tiers, sample, weights = weight)
  calls <- fill$calls
  labels <- if (length(mk)) combination_labels(calls, tiers)
            else rep("(no decision)", nrow(calls))
  combination <- pool_rare_labels(labels, weight, min_combination_cells)

  # Counts over ALL cells, per sample and combination, zero-filled so that a
  # mean over samples counts a sample without the combination as 0 %.
  sk <- paste(sample, key, sep = "\r")
  per <- lapply(samples, function(s) {
    kc <- per_sample[[s]]$key_counts
    comb <- combination[match(paste(s, names(kc), sep = "\r"), sk)]
    if (anyNA(comb)) {
      stop("internal: sample '", s, "' has a call pattern with no kept cell.",
           call. = FALSE)
    }
    tapply(as.numeric(kc), comb, sum)
  })
  tot <- tapply(weight, combination, sum)
  comb_levels <- names(sort(tot, decreasing = TRUE))
  counts <- do.call(rbind, lapply(seq_along(samples), function(i) {
    v <- per[[i]]
    data.frame(sample = samples[i], combination = comb_levels,
               cells = ifelse(is.na(v[comb_levels]), 0, v[comb_levels]),
               population_cells = per_sample[[samples[i]]]$population_cells,
               stringsAsFactors = FALSE, row.names = NULL)
  }))

  list(calls = calls, values = values, cuts = cuts, anchors = anchors,
       weight = weight,
       combination = combination, decisions = decisions, counts = counts,
       tiers = tiers, filled = fill$filled)
}

# "col=level;col2=level" -> named list of single reference levels.
parse_split_reference <- function(spec) {
  if (is.null(spec) || !nzchar(trimws(spec))) return(list())
  parts <- trimws(strsplit(spec, ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  out <- list()
  for (p in parts) {
    kv <- regmatches(p, regexpr("=", p), invert = TRUE)[[1]]
    if (length(kv) != 2L || !nzchar(trimws(kv[1])) || !nzchar(trimws(kv[2]))) {
      stop("--split_reference entry '", p, "' is not of the form column=level.",
           call. = FALSE)
    }
    out[[trimws(kv[1])]] <- trimws(kv[2])
  }
  out
}

# -----------------------------------------------------------------------------
# Background subtraction (--background, --pair_by): the frequency of a
# combination in a stimulated tube minus its frequency in the unstimulated
# tube of the same patient and visit. Cells cannot be subtracted, only their
# shares per tube, so the embedding is untouched; the per-sample percentages
# behind the frequency tables, legends and fold-change maps are what change.
# -----------------------------------------------------------------------------

# "stimulation=unstimulated" -> list(column = "stimulation", level = "unstimulated").
parse_background <- function(spec) {
  if (is.null(spec) || !nzchar(trimws(spec))) return(NULL)
  parts <- strsplit(trimws(spec), "=", fixed = TRUE)[[1]]
  if (length(parts) != 2L || !all(nzchar(trimws(parts)))) {
    stop("--background takes column=level, e.g. stimulation=unstimulated (got '",
         spec, "').", call. = FALSE)
  }
  list(column = trimws(parts[1]), level = trimws(parts[2]))
}

# Background sample of every other sample: the sample whose `column` is `level`
# and whose `pair_by` values are the same. Named character (sample ->
# background); samples without one are left out, with a warning.
#
# @param samples Sample labels, one per row of `ann`.
# @param ann Data frame of sample annotations (every metadata column).
background_pairs <- function(samples, ann, column, level, pair_by) {
  stopifnot(nrow(ann) == length(samples), length(pair_by) > 0L)
  miss <- setdiff(c(column, pair_by), names(ann))
  if (length(miss)) {
    stop("background pairing needs metadata column(s) ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  is_bg <- as.character(ann[[column]]) == level
  is_bg[is.na(is_bg)] <- FALSE
  if (!any(is_bg)) {
    stop("no sample has ", column, " = ", level, ", so there is no background ",
         "to subtract (a --filter may have left the background tubes out).",
         call. = FALSE)
  }
  key <- do.call(paste, c(lapply(pair_by, function(cl) as.character(ann[[cl]])),
                          sep = "\r"))
  bgs <- split(samples[is_bg], key[is_bg])
  multi <- lengths(bgs) > 1L
  if (any(multi)) {
    warning(sum(multi), " ", paste(pair_by, collapse = "/"), " group(s) have ",
            "more than one background sample; the first of each is used.",
            call. = FALSE)
  }
  first <- vapply(bgs, function(x) x[1], character(1))
  out <- stats::setNames(unname(first[key[!is_bg]]), samples[!is_bg])
  if (anyNA(out)) {
    warning(sum(is.na(out)), " sample(s) have no background sample with the ",
            "same ", paste(pair_by, collapse = ", "), "; they are left out of ",
            "the background-subtracted frequencies: ",
            paste(utils::head(names(out)[is.na(out)], 5), collapse = ", "),
            if (sum(is.na(out)) > 5) ", ..." else "", call. = FALSE)
  }
  out[!is.na(out)]
}

# The combinations background subtraction applies to: those with a positive
# call below the lineage, e.g. "CD4+ | IFNy+ TNFa+", not "CD4+ | none",
# "CD4- CD8-" or "other (rare)". Labels without the " | " of
# combination_labels() (subgate names) are all kept.
response_combinations <- function(combos) {
  combos <- as.character(combos)
  if (!any(grepl(" | ", combos, fixed = TRUE))) return(combos)
  combos[grepl(" | ", combos, fixed = TRUE) & !grepl("\\| none$", combos)]
}

# Per-sample percentages (samples x combinations, .sp_percent_matrix) minus
# each sample's background, negatives set to 0 (a combination no more frequent
# than in the unstimulated tube has no response). Only samples with a
# background are kept; the background samples themselves drop out.
#
# Only `subtract` combinations are background-subtracted: the responses. For
# a bulk population (CD4+ with no cytokine, CD4- CD8-) stimulated minus
# unstimulated is a near-zero difference of two large numbers, and its fold
# change is noise; those columns become NA ("not a response").
subtract_background <- function(pct, background, subtract = colnames(pct)) {
  pop <- attr(pct, "population")
  s <- names(background)[names(background) %in% rownames(pct) &
                           background %in% rownames(pct)]
  out <- pct[s, , drop = FALSE] - pct[background[s], , drop = FALSE]
  out[out < 0] <- 0
  out[, !colnames(out) %in% subtract] <- NA_real_
  rownames(out) <- s
  attr(out, "population") <- pop[s]
  attr(out, "background") <- background[s]
  out
}
