#!/usr/bin/env Rscript
# =============================================================================
# scripts/plot_polyfunctionality.R: polyfunctionality maps of one patient on a
# finished umap_from_flowjo.R run's embedding.
#
# Every CD4 or CD8 T cell is coloured by how many of the cytokine gates it is
# in (0 to 4 of IFNy, IL-2, TNFa, CD107a), read from its gate-combination label
# ("CD4+ | IFNy+ TNFa+" -> 2). CD4-CD8- and CD4+CD8+ cells have no cytokine
# gates, so they are left out (--other_cells grey draws them in light grey
# behind the counted cells). Nothing is
# recomputed: the coordinates are the run's (data/cells.csv.gz), so the maps
# are on the same embedding as all its other plots. The run must colour by gate
# combination (--subpop_mode combinations or --layout gates).
#
# One map per timepoint x stimulation, and one per timepoint with all its
# cells, written to <run>/polyfunctionality/<patient>/. Every map has the
# frame of the whole embedding, so maps of different patients line up.
#
# Usage:
#   Rscript scripts/plot_polyfunctionality.R --run DIR --patient MTB073 \
#     [--formats svg] [--colours "#DCE9F6,#9CC3E5,#4A8FC9,#0B4F8A,#062F52"]
#   Rscript scripts/plot_polyfunctionality.R --run DIR --patient all --pooled --n_cores 12
#
# --pooled also draws the chosen patients together, into <output>/all_patients/:
# per timepoint x stimulation, per timepoint, and over all timepoints.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run", type = "character", help = "Run folder (with data/cells.csv.gz)"),
  make_option("--patient", type = "character",
              help = "Patient(s), comma-separated, as in the metadata column 'patient'; 'all' for every one"),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Where to write; one subfolder per patient [default <run>/polyfunctionality]"),
  make_option("--other_cells", type = "character", default = "drop",
              help = paste0("CD4-CD8- and CD4+CD8+ cells, which have no cytokine gates: ",
                            "'drop' (not drawn) or 'grey' (light grey behind) [default %default]")),
  make_option("--point_size", type = "double", default = 0.4,
              help = "Size of every point, the same in every map [default %default]"),
  make_option("--emphasise", action = "store_true", default = FALSE,
              help = paste0("Draw cells with more positive gates larger, and every point of a ",
                            "map with few cells larger (the drawing before 2026-10-21)")),
  make_option("--size_by_polyfunctionality", action = "store_true", default = FALSE,
              help = paste0("Draw cells with more positive gates larger: --point_size times ",
                            "1, 1.5, 2, 2.5, 3 for 0 to 4 gates (or --sizes), the same in every map")),
  make_option("--sizes", type = "character", default = NULL,
              help = "Five point sizes for 0 to 4 gates, with --size_by_polyfunctionality"),
  make_option("--alpha", type = "double", default = NULL,
              help = "Opacity of the points [default 0.85, 1 with --projector]"),
  make_option("--equal_cells", action = "store_true", default = FALSE,
              help = paste0("Draw the same number of randomly chosen CD4 and CD8 T cells in every ",
                            "timepoint's map of a patient (per stimulation, and for all tubes): ",
                            "the smallest timepoint's. Percentages stay those of all cells")),
  make_option("--projector", action = "store_true", default = FALSE,
              help = paste0("Readable from afar and on a weak projector: larger dark text, ",
                            "opaque points, thicker axis arrows")),
  make_option("--centre_labels", action = "store_true", default = FALSE,
              help = paste0("Title and legend inside the map, near its middle and in larger text, ",
                            "so they stay readable when the figure is shown small")),
  make_option("--pooled", action = "store_true", default = FALSE,
              help = "Also draw all chosen patients together, in all_patients/"),
  make_option("--n_cores", type = "integer", default = 1L,
              help = "Patients drawn in parallel [default %default]"),
  make_option("--colours", type = "character",
              default = "#DCE9F6,#9CC3E5,#4A8FC9,#0B4F8A,#062F52",
              help = "Five colours for 0 to 4 positive gates [default %default]"),
  make_option("--formats", type = "character", default = "svg",
              help = "Comma-separated: svg, png, pdf [default %default]"),
  make_option("--timepoint_column", type = "character", default = "timepoint"),
  make_option("--stimulation_column", type = "character", default = "stimulation")
)))
for (nm in c("run", "patient")) {
  if (is.null(opt[[nm]])) stop("--", nm, " is required.", call. = FALSE)
}
script_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
src_dir <- normalizePath(file.path(script_dir, "..", "src"), mustWork = TRUE)
source(file.path(src_dir, "R_shared", "plot_utils.R"))
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))
source(file.path(src_dir, "visualization", "umap_style.R"))

if (!opt$other_cells %in% c("drop", "grey")) {
  stop("--other_cells must be 'drop' or 'grey'.", call. = FALSE)
}
cols <- trimws(strsplit(opt$colours, ",", fixed = TRUE)[[1]])
if (length(cols) != 5L) stop("--colours needs five colours (0 to 4).", call. = FALSE)
formats <- trimws(strsplit(opt$formats, ",", fixed = TRUE)[[1]])
out_root <- if (is.null(opt$output_dir)) file.path(opt$run, "polyfunctionality") else opt$output_dir

cells_file <- run_find(opt$run, "data", "cells.csv.gz")
meta_file <- run_find(opt$run, "data", "sample_metadata.csv")
if (is.na(cells_file) || is.na(meta_file)) {
  stop("the run needs data/cells.csv.gz and data/sample_metadata.csv.", call. = FALSE)
}
md <- fread(meta_file)
md[, sample := sub("\\.fcs$", "", basename(as.character(md[[1]])))]
patients <- if (identical(tolower(opt$patient), "all")) sort(unique(md$patient)) else
  trimws(strsplit(opt$patient, ",", fixed = TRUE)[[1]])
missing <- setdiff(patients, md$patient)
if (length(missing)) stop("no sample of patient(s) ", paste(missing, collapse = ", "), ".",
                          call. = FALSE)

tcol <- paste0("annot_", opt$timepoint_column)
scol <- paste0("annot_", opt$stimulation_column)
all_cells <- fread(cmd = paste("zcat", shQuote(cells_file)),
                   select = c("UMAP1", "UMAP2", "sample", "label", tcol, scol))
setnames(all_cells, c(tcol, scol), c("timepoint", "stimulation"))
# The frame of the whole embedding, shared by every map.
lims_x <- range(all_cells$UMAP1); lims_y <- range(all_cells$UMAP2)
all_cells[, patient := md$patient[match(sample, md$sample)]]
all_cells <- all_cells[patient %in% patients]
x <- all_cells
if (!any(grepl(" | ", x$label, fixed = TRUE))) {
  stop("the run's labels are not gate combinations; rerun with ",
       "--subpop_mode combinations.", call. = FALSE)
}

# Positive gates per cell from its label: the "+" markers after " | ".
cyto <- sub("^[^|]*\\| ?", "", x$label)
x[, n_pos := ifelse(!grepl(" | ", label, fixed = TRUE), NA_integer_,
                    ifelse(cyto == "none", 0L,
                           lengths(regmatches(cyto, gregexpr("\\+", cyto)))))]
cat(sprintf("%d patient(s): %s cells in %d tubes\n", length(patients),
            format(nrow(x), big.mark = ","), length(unique(x$sample))))

size_of <- c(0.25, 0.45, 0.6, 0.75, 0.9)
poly_size <- if (is.null(opt$sizes)) opt$point_size * c(1, 1.5, 2, 2.5, 3) else
  as.numeric(strsplit(opt$sizes, ",", fixed = TRUE)[[1]])
if (length(poly_size) != 5L || anyNA(poly_size)) stop("--sizes needs five numbers.", call. = FALSE)
# Legend dots grow in the same order as the map's, so the
# smallest stays visible and the largest does not dominate (2 to 3.5).
leg_size <- 2 + 1.5 * (poly_size - min(poly_size)) / max(diff(range(poly_size)), 1e-9)
lab_of <- c("0", "1", "2", "3", "4")
raster <- requireNamespace("ggrastr", quietly = TRUE)
base_size <- if (opt$projector) 17 else 11
pt_alpha <- if (!is.null(opt$alpha)) opt$alpha else if (opt$projector) 1 else 0.85

# Smallest number of counted cells over the maps in `sets` (--equal_cells).
n_equal <- function(sets) {
  if (!opt$equal_cells) return(NULL)
  min(vapply(sets, function(d) sum(!is.na(d$n_pos)), integer(1)))
}

draw <- function(d, title, stem, out, n_draw = NULL) {
  counted <- d[!is.na(n_pos)]
  n_t <- nrow(counted)
  tab <- tabulate(counted$n_pos + 1L, nbins = 5L)
  pct <- if (n_t > 0) sprintf("%.2f %%", 100 * tab / n_t) else rep("-", 5L)
  leg <- sprintf("%s  (%s, %s)", lab_of, format(tab, big.mark = ","), pct)
  counted[, lvl := factor(n_pos, levels = 0:4, labels = leg)]
  drawn <- if (is.null(n_draw)) "" else
    sprintf(", %s drawn at random (as many as in the other timepoints)",
            format(min(n_draw, n_t), big.mark = ","))
  keep_frac <- 1
  if (!is.null(n_draw) && n_draw < n_t) {
    set.seed(1L)
    keep_frac <- n_draw / n_t
    counted <- counted[sample.int(n_t, n_draw)]
  }
  counted <- counted[order(n_pos)]          # most polyfunctional drawn last, on top
  # One size for every point unless --emphasise: then larger with more positive
  # gates, and larger in a map with few cells, so a single tube stays readable.
  # --size_by_polyfunctionality: larger with more positive gates only.
  grow <- if (opt$emphasise) min(3, max(1, sqrt(150000 / max(nrow(d), 1)))) else 1
  counted[, sz := if (opt$emphasise) size_of[n_pos + 1L] * grow else
                  if (opt$size_by_polyfunctionality) poly_size[n_pos + 1L] else opt$point_size]
  pt <- function(...) {
    g <- geom_point(..., stroke = 0)
    if (raster) ggrastr::rasterise(g, dpi = 300) else g
  }
  others <- if (opt$other_cells == "grey") d[is.na(n_pos)] else d[0]
  if (keep_frac < 1) others <- others[sample.int(nrow(others), round(keep_frac * nrow(others)))]
  p <- ggplot() +
    pt(data = others, aes(UMAP1, UMAP2), colour = "#E4E4E4",
       size = if (opt$emphasise) 0.25 * grow else opt$point_size, alpha = 0.7) +
    pt(data = counted, aes(UMAP1, UMAP2, colour = lvl, size = sz), alpha = pt_alpha) +
    scale_size_identity() +
    scale_colour_manual(values = stats::setNames(cols, leg), drop = FALSE,
                        name = "Positive gates per cell\n(IFNy, IL-2, TNFa, CD107a)",
                        guide = guide_legend(override.aes = list(
                          size = (if (opt$size_by_polyfunctionality) leg_size
                                  else 3) * (if (opt$projector) 2 else 1), alpha = 1))) +
    umap_axis_arrows(lims_x, lims_y, base_size = if (opt$projector) 15 else 9) +
    coord_fixed(xlim = lims_x, ylim = lims_y) + theme_umap_clean(base_size) +
    labs(title = title,
         subtitle = sprintf(paste0("%s CD4 and CD8 T cells from %d tube(s)%s\n%s CD4-CD8- ",
                                   "and CD4+CD8+ cells %s (no cytokine gates)"),
                            format(n_t, big.mark = ","), length(unique(d$sample)), drawn,
                            format(sum(is.na(d$n_pos)), big.mark = ","),
                            if (opt$other_cells == "grey") "in grey" else "not shown"),
         caption = paste0("Same embedding as the whole cohort. Percentages: of ",
                          if (is.null(n_draw)) "the CD4 and CD8 T cells shown" else
                            "all CD4 and CD8 T cells of the tubes, not only those drawn",
                          ".\nCells with more positive gates are drawn on top",
                          if (opt$emphasise || opt$size_by_polyfunctionality) " and larger."
                          else "; every point has the same size.")) +
    theme(plot.caption = element_text(hjust = 0), plot.caption.position = "plot")
  if (opt$projector) {
    p <- p + theme(plot.title = element_text(colour = "black", face = "bold", hjust = 0),
                   plot.subtitle = element_text(colour = "grey15", hjust = 0),
                   plot.title.position = "plot",
                   plot.caption = element_text(colour = "grey20", hjust = 0),
                   legend.title = element_text(colour = "black", face = "bold"),
                   legend.text = element_text(colour = "black"),
                   legend.key.size = grid::unit(0.7, "cm"))
  }
  if (opt$centre_labels) {
    # The title goes in the empty upper middle of the panel and the legend just
    # below it, so that it may overlap the upper clusters slightly.
    p <- p + labs(title = NULL) +
      annotate("text", x = mean(lims_x), y = lims_y[1] + 0.9 * diff(lims_y), label = title,
               size = 6, colour = "black", hjust = 0.5) +
      theme(legend.position = "inside", legend.position.inside = c(0.5, 0.62),
            legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA),
            legend.title = element_text(size = 14, colour = "black"),
            legend.text = element_text(size = 13, colour = "black"),
            legend.key.size = grid::unit(0.6, "cm")) +
      coord_fixed(xlim = lims_x, ylim = lims_y, clip = "off")
  }
  for (fmt in formats) {
    ggsave(file.path(out, paste0(stem, ".", fmt)), p, width = 10, height = 7.5,
           dpi = 200, bg = "white")
  }
  cat("  ", stem, ": ", format(nrow(d), big.mark = ","), " cells\n", sep = "")
}

short <- function(tp) sub(" .*", "", tp)
one_patient <- function(pt) {
  xp <- x[patient == pt]
  out <- file.path(out_root, create_safe_filename(pt))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  n <- 0L
  tps <- sort(unique(xp$timepoint))
  n_all <- n_equal(lapply(tps, function(tp) xp[timepoint == tp]))
  for (tp in tps) {
    for (st in sort(unique(xp[timepoint == tp, stimulation]))) {
      n_st <- n_equal(lapply(tps, function(t) xp[timepoint == t & stimulation == st])
                      [vapply(tps, function(t) any(xp$timepoint == t & xp$stimulation == st), TRUE)])
      draw(xp[timepoint == tp & stimulation == st], sprintf("%s, %s, %s", pt, tp, st),
           create_safe_filename(sprintf("polyfunctionality_%s_%s_%s", pt, short(tp), st)), out,
           n_st)
      n <- n + 1L
    }
    draw(xp[timepoint == tp], sprintf("%s, %s, all tubes", pt, tp),
         create_safe_filename(sprintf("polyfunctionality_%s_%s_all", pt, short(tp))), out, n_all)
    n <- n + 1L
  }
  n
}
# All chosen patients together: per timepoint (x stimulation) and over all
# timepoints, the same maps as for one patient.
pooled <- function() {
  out <- file.path(out_root, "all_patients")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  who <- sprintf("All %d patients", length(unique(x$patient)))
  sets <- c(lapply(sort(unique(x$timepoint)), function(tp) list(tp = tp, label = tp,
                                                                  tag = short(tp))),
            list(list(tp = NULL, label = "all timepoints", tag = "all_timepoints")))
  n <- 0L
  for (st_set in sets) {
    xs <- if (is.null(st_set$tp)) x else x[timepoint == st_set$tp]
    for (st in sort(unique(xs$stimulation))) {
      draw(xs[stimulation == st], sprintf("%s, %s, %s", who, st_set$label, st),
           create_safe_filename(sprintf("polyfunctionality_all_patients_%s_%s",
                                        st_set$tag, st)), out)
      n <- n + 1L
    }
    draw(xs, sprintf("%s, %s, all tubes", who, st_set$label),
         create_safe_filename(sprintf("polyfunctionality_all_patients_%s_all", st_set$tag)),
         out)
    n <- n + 1L
  }
  n
}

setDTthreads(1L)
n_made <- if (opt$n_cores > 1L && .Platform$OS.type == "unix") {
  unlist(parallel::mclapply(patients, one_patient, mc.cores = opt$n_cores,
                            mc.preschedule = FALSE))
} else {
  vapply(patients, one_patient, integer(1))
}
if (length(n_made) != length(patients) || !all(vapply(n_made, is.numeric, logical(1)))) {
  stop("a worker failed; rerun with fewer --n_cores.", call. = FALSE)
}
if (opt$pooled) n_made <- c(n_made, pooled = pooled())
cat(sprintf("Written %d map(s) for %d patient(s)%s to %s\n", sum(n_made),
            length(patients), if (opt$pooled) " and the patients pooled" else "",
            out_root))
