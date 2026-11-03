#!/usr/bin/env Rscript
# =============================================================================
# scripts/plot_polyfunctionality_timepoints.R: compare the timepoints of a
# patient on a finished umap_from_flowjo.R run's embedding, coloured by how
# many cytokine gates each CD4 or CD8 T cell is in (0 to 4 of IFNy, IL-2, TNFa,
# CD107a; see scripts/plot_polyfunctionality.R).
#
# Maps of different timepoints only compare when they hold the same number of
# cells: a timepoint with more tubes otherwise looks more polyfunctional. So
# every map of a patient draws the same number of randomly chosen cells (the
# smallest stimulated timepoint of that patient), and the percentages, from
# all cells, are plotted beside them.
#
# Written to --output_dir:
#   <patient>/polyfunctionality_timepoints_<patient>   maps of the stimulated
#       tubes, one per timepoint, over the percentage of cells with 1 to 4
#       positive gates per timepoint (stimulated and unstimulated)
#   all_patients/polyfunctionality_timepoints_grid     the maps, one row per
#       patient (equal cells within a row)
#   all_patients/polyfunctionality_timepoints_percent  the percentages, one
#       row per patient
#   polyfunctionality_timepoints.csv                   the numbers
#
# With --single_maps, instead of the figures: one map per file,
#   <patient>/polyfunctionality_<patient>_<timepoint>  (same cells per map
#       within a patient, as above)
#
# Usage:
#   Rscript scripts/plot_polyfunctionality_timepoints.R --run DIR \
#     --patient MTB059,MTB074,MTB079,MTB087 --output_dir DIR
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run", type = "character", help = "Run folder (with data/cells.csv.gz)"),
  make_option("--patient", type = "character", help = "Patient(s), comma-separated"),
  make_option("--output_dir", type = "character", help = "Where to write the figures"),
  make_option("--colours", type = "character",
              default = "#C6DBEF,#4FA3E0,#1F6FC5,#1A3A9A,#050B30",
              help = "Colours of 0 to 4 positive gates [default %default]"),
  make_option("--point_size", type = "double", default = 1,
              help = "Scales every point (0 to 4 gates: 0.2 to 1.4 times this) [default %default]"),
  make_option("--formats", type = "character", default = "svg,png"),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--timepoint_column", type = "character", default = "timepoint"),
  make_option("--stimulation_column", type = "character", default = "stimulation"),
  make_option("--single_maps", action = "store_true", default = FALSE,
              help = "Write every map to its own file instead of the combined figures"),
  make_option("--stimulated", type = "character", default = "stimulated",
              help = "Stimulation level drawn on the maps [default %default]")
)))
for (nm in c("run", "patient", "output_dir")) {
  if (is.null(opt[[nm]])) stop("--", nm, " is required.", call. = FALSE)
}
script_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
src_dir <- normalizePath(file.path(script_dir, "..", "src"), mustWork = TRUE)
source(file.path(src_dir, "R_shared", "plot_utils.R"))
source(file.path(src_dir, "visualization", "umap_from_flowjo.R"))
source(file.path(src_dir, "visualization", "umap_style.R"))

cols <- trimws(strsplit(opt$colours, ",", fixed = TRUE)[[1]])
if (length(cols) != 5L) stop("--colours needs five colours (0 to 4).", call. = FALSE)
names(cols) <- as.character(0:4)
sizes <- opt$point_size * c(0.2, 0.5, 0.8, 1.1, 1.4)
formats <- trimws(strsplit(opt$formats, ",", fixed = TRUE)[[1]])
patients <- trimws(strsplit(opt$patient, ",", fixed = TRUE)[[1]])

cells_file <- run_find(opt$run, "data", "cells.csv.gz")
meta_file <- run_find(opt$run, "data", "sample_metadata.csv")
if (is.na(cells_file) || is.na(meta_file)) {
  stop("the run needs data/cells.csv.gz and data/sample_metadata.csv.", call. = FALSE)
}
md <- fread(meta_file)
md[, sample := sub("\\.fcs$", "", basename(as.character(md[[1]])))]
missing <- setdiff(patients, md$patient)
if (length(missing)) stop("no sample of patient(s) ", paste(missing, collapse = ", "), ".",
                          call. = FALSE)
tcol <- paste0("annot_", opt$timepoint_column)
scol <- paste0("annot_", opt$stimulation_column)
x <- fread(cmd = paste("zcat", shQuote(cells_file)),
           select = c("UMAP1", "UMAP2", "sample", "label", tcol, scol))
setnames(x, c(tcol, scol), c("timepoint", "stimulation"))
lims_x <- range(x$UMAP1); lims_y <- range(x$UMAP2)
x[, patient := md$patient[match(sample, md$sample)]]
x <- x[patient %in% patients & grepl(" | ", label, fixed = TRUE)]
cyto <- sub("^[^|]*\\| ?", "", x$label)
x[, n_pos := ifelse(cyto == "none", 0L, lengths(regmatches(cyto, gregexpr("\\+", cyto))))]
x[, tp := sub(" .*", "", timepoint)]
# Visit weeks per timepoint, e.g. "V2 (W17, W42)" where a timepoint has two.
if ("week" %in% names(md)) {
  wk <- unique(data.table(patient = md$patient, tp = sub(" .*", "", md[[opt$timepoint_column]]),
                          week = sub("^W0*(?=[0-9])", "W", md$week, perl = TRUE)))
  wk <- wk[order(as.integer(sub("^W", "", week)))]
  wk <- wk[, .(weeks = paste(unique(week), collapse = ", ")), by = .(patient, tp)]
  x[wk, weeks := i.weeks, on = .(patient, tp)]
} else x[, weeks := NA_character_]
x[, panel := ifelse(is.na(weeks), tp, sprintf("%s (%s)", tp, weeks))]

# --- Percentages, from all cells. -------------------------------------------
pct <- x[, .(n_cells = .N, n0 = sum(n_pos == 0), n1 = sum(n_pos == 1), n2 = sum(n_pos == 2),
             n3 = sum(n_pos == 3), n4 = sum(n_pos == 4)),
         by = .(patient, tp, stimulation)]
long <- melt(pct, id.vars = c("patient", "tp", "stimulation", "n_cells"),
             measure.vars = paste0("n", 1:4), variable.name = "gates", value.name = "n")
long[, gates := sub("^n", "", gates)][, percent := 100 * n / n_cells]
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(long[order(patient, stimulation, tp, gates)],
       file.path(opt$output_dir, "polyfunctionality_timepoints.csv"))

projector_theme <- function(base = 17) {
  theme(plot.title = element_text(size = base + 3, face = "bold", colour = "black"),
        plot.subtitle = element_text(size = base - 2, colour = "grey15"),
        plot.title.position = "plot")
}

percent_plot <- function(d, base = 17) {
  d <- copy(d)
  d[, gates := factor(gates, levels = 1:4,
                      labels = paste0(1:4, " positive gate", c("", "s", "s", "s")))]
  d[, stimulation := factor(stimulation, levels = c(opt$stimulated,
                                                    setdiff(unique(stimulation), opt$stimulated)))]
  ggplot(d, aes(tp, percent, group = stimulation)) +
    geom_line(aes(linetype = stimulation, colour = gates), linewidth = 1.2) +
    geom_point(aes(shape = stimulation, colour = gates), size = 4, stroke = 1.4,
               fill = "white") +
    facet_wrap(~ gates, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = stats::setNames(cols[-1], levels(d$gates)), guide = "none") +
    scale_shape_manual(values = c(16, 21), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12)),
                       labels = function(v) paste0(format(v, trim = TRUE, drop0trailing = TRUE), " %")) +
    labs(x = NULL, y = "% of CD4 and CD8 T cells") +
    theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold", colour = "black"),
          axis.text = element_text(colour = "black"),
          legend.position = "top", legend.justification = "left",
          legend.key.width = grid::unit(1.6, "cm"))
}

map_panel <- function(d, n_draw, base = 17) {
  share <- d[, .(lab = sprintf("2+ gates: %.2f %%", 100 * mean(n_pos >= 2))), by = panel]
  set.seed(opt$seed)
  d <- d[, .SD[sample.int(.N, min(.N, n_draw))], by = panel]
  d <- d[order(n_pos)]
  g <- geom_point(aes(colour = factor(n_pos), size = factor(n_pos)), stroke = 0)
  if (requireNamespace("ggrastr", quietly = TRUE)) g <- ggrastr::rasterise(g, dpi = 300)
  ggplot(d, aes(UMAP1, UMAP2)) + g +
    geom_text(data = share, aes(label = lab), x = lims_x[2], y = lims_y[2], hjust = 1,
              vjust = 1, size = (base - 2) / .pt, fontface = "bold", colour = "black",
              inherit.aes = FALSE) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = cols, name = "Positive gates\nper cell", drop = FALSE) +
    scale_size_manual(values = stats::setNames(sizes, names(cols)), guide = "none") +
    guides(colour = guide_legend(override.aes = list(size = 6))) +
    coord_fixed(xlim = lims_x, ylim = lims_y) +
    theme_void(base_size = base) +
    theme(strip.text = element_text(face = "bold", size = base + 1, margin = margin(b = 6)),
          legend.title = element_text(face = "bold"), panel.spacing = grid::unit(18, "pt"),
          panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5))
}

save_fig <- function(p, dir, stem, w, h) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (fmt in formats) ggsave(file.path(dir, paste0(stem, ".", fmt)), p, width = w, height = h,
                              dpi = 200, bg = "white")
  cat("  ", file.path(basename(dir), stem), "\n", sep = "")
}

stim <- x[stimulation == opt$stimulated]
n_draw <- stim[, .N, by = .(patient, panel)][, .(n_draw = min(N)), by = patient]
n_tp <- max(stim[, uniqueN(panel), by = patient]$V1)

if (opt$single_maps) {
  for (pt in patients) {
    nd <- n_draw[patient == pt, n_draw]
    sp <- stim[patient == pt]
    for (pn in sort(unique(sp$panel))) {
      tp <- sp[panel == pn, tp[1]]
      p <- map_panel(sp[panel == pn], nd) +
        labs(title = sprintf("%s, %s, %s", pt, pn, opt$stimulated),
             subtitle = sprintf("%s randomly chosen CD4 and CD8 T cells
(the same in every map of %s)",
                                format(nd, big.mark = ","), pt)) +
        projector_theme() + theme(strip.text = element_blank())
      save_fig(p, file.path(opt$output_dir, create_safe_filename(pt)),
               create_safe_filename(sprintf("polyfunctionality_%s_%s", pt, tp)), 8.5, 6.5)
    }
  }
  cat("Written to ", opt$output_dir, "\n", sep = "")
  quit(status = 0)
}

for (pt in patients) {
  nd <- n_draw[patient == pt, n_draw]
  p <- map_panel(stim[patient == pt], nd) / percent_plot(long[patient == pt]) +
    plot_layout(heights = c(1.5, 1)) +
    plot_annotation(title = sprintf("%s: polyfunctionality over time", pt),
                    subtitle = sprintf(paste0("Maps: %s tubes, %s randomly chosen CD4 and CD8 ",
                                              "T cells per timepoint, so the maps compare. ",
                                              "Percentages: all cells."),
                                       opt$stimulated, format(nd, big.mark = ",")),
                    theme = projector_theme())
  save_fig(p, file.path(opt$output_dir, create_safe_filename(pt)),
           create_safe_filename(paste0("polyfunctionality_timepoints_", pt)), 18, 10)
}

# All patients: one row each, every row with its own (equal) number of cells.
rows <- lapply(patients, function(pt) {
  map_panel(stim[patient == pt], n_draw[patient == pt, n_draw]) +
    labs(tag = sprintf("%s\n%s cells\nper map", pt,
                       format(n_draw[patient == pt, n_draw], big.mark = ","))) +
    theme(plot.tag.position = "left", plot.tag = element_text(size = 15, face = "bold",
                                                              hjust = 1))
})
widths <- vapply(patients, function(pt) stim[patient == pt, uniqueN(panel)], integer(1))
rows <- Map(function(p, k) if (k < n_tp) p + plot_spacer() +
              plot_layout(widths = c(k, n_tp - k)) else p, rows, widths)
grid <- wrap_plots(rows, ncol = 1) + plot_layout(guides = "collect") +
  plot_annotation(title = "Polyfunctionality over time, patients with a V0 visit",
                  subtitle = sprintf(paste0("Maps of the %s tubes. Within a row every map has ",
                                            "the same number of cells; percentages from all cells."),
                                     opt$stimulated),
                  theme = projector_theme())
save_fig(grid, file.path(opt$output_dir, "all_patients"), "polyfunctionality_timepoints_grid",
         4.2 * n_tp + 3, 3.4 * length(patients) + 1.5)

pp <- percent_plot(long) +
  facet_wrap(~ patient + gates, ncol = 4, scales = "free_y",
             labeller = function(l) list(paste0(l$patient, ", ", l$gates))) +
  plot_annotation(title = "Cells with 1 to 4 positive gates per timepoint",
                  subtitle = "Percent of all CD4 and CD8 T cells of the timepoint's tubes.",
                  theme = projector_theme())
save_fig(pp, file.path(opt$output_dir, "all_patients"), "polyfunctionality_timepoints_percent",
         18, 3 * length(patients) + 2)
cat("Written to ", opt$output_dir, "\n", sep = "")
