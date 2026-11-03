#!/usr/bin/env Rscript
# =============================================================================
# scripts/compare_umap_runs.R: put two umap_from_flowjo.R runs side by side.
#
# Made for comparing a run with its background-subtracted twin (--background,
# projected onto the same map with --umap_model), but works for any two runs
# that split by the same metadata column. For that column it writes, into
# --output_dir:
#
#   compare_<col>_<level>.png   the two runs' fold-change maps of a level,
#                               side by side (needs the PNGs of both runs)
#   compare_<col>_scatter.png   every combination's log2 fold change in run A
#                               against run B, one facet per level
#   compare_<col>.csv           the numbers behind it, one row per level and
#                               combination, both runs' frequencies and changes
#
# Usage:
#   Rscript scripts/compare_umap_runs.R --run_a DIR --run_b DIR \
#     --label_a "Stimulated" --label_b "Stimulated minus unstimulated" \
#     --column timepoint --output_dir DIR
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run_a", type = "character", help = "First run folder"),
  make_option("--run_b", type = "character", help = "Second run folder"),
  make_option("--label_a", type = "character", default = "Run A",
              help = "Name of the first run in the figures [default %default]"),
  make_option("--label_b", type = "character", default = "Run B",
              help = "Name of the second run in the figures [default %default]"),
  make_option("--column", type = "character",
              help = "The --split_by column both runs share, e.g. timepoint"),
  make_option("--output_dir", type = "character",
              help = "Where to write the comparison"),
  make_option("--label_fold", type = "double", default = 2,
              help = paste0("Label combinations whose fold change in either ",
                            "run is at least this (or 1/this) [default %default]"))
)))
for (nm in c("run_a", "run_b", "column", "output_dir")) {
  if (is.null(opt[[nm]])) stop("--", nm, " is required.", call. = FALSE)
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

safe <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  substr(gsub("^[_-]+|[_-]+$", "", x), 1L, 100L)
}
col_tag <- safe(opt$column)

# A run keeps its tables in changes/<column>/ (or flat, before 2026-10-15).
find_in_run <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  f[!grepl("/checkpoints/", f, fixed = TRUE)]
}
read_change <- function(dir) {
  name <- paste0("combination_change_", col_tag, ".csv")
  f <- find_in_run(dir, paste0("^", gsub(".", "\\.", name, fixed = TRUE), "$"))[1]
  if (is.na(f)) {
    stop("no ", name, " in ", dir, ": was it split by '", opt$column,
         "'?", call. = FALSE)
  }
  utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}
a <- read_change(opt$run_a)
b <- read_change(opt$run_b)
keep <- c("level", "reference", "combination", "mean_pct_level",
          "mean_pct_reference", "log2_fold_change", "n_samples_level",
          "n_samples_reference")
m <- merge(a[, keep], b[, keep], by = c("level", "reference", "combination"),
           suffixes = c("_a", "_b"), all = TRUE)
m$change_of_log2_fold_change <- m$log2_fold_change_b - m$log2_fold_change_a
utils::write.csv(m, file.path(opt$output_dir, paste0("compare_", col_tag, ".csv")),
                 row.names = FALSE)

# --- Scatter: fold change in A against B, per level. -------------------------
lineage <- ifelse(grepl("^CD4\\+ \\|", m$combination), "CD4+",
           ifelse(grepl("^CD8\\+ \\|", m$combination), "CD8+", "other"))
d <- data.frame(m, lineage = lineage, stringsAsFactors = FALSE)
d <- d[is.finite(d$log2_fold_change_a) | is.finite(d$log2_fold_change_b), ]
cap <- 3
d$x <- pmax(pmin(d$log2_fold_change_a, cap), -cap)
d$y <- pmax(pmin(d$log2_fold_change_b, cap), -cap)
lim <- log2(opt$label_fold)
d$lab <- ifelse(pmax(abs(d$x), abs(d$y), na.rm = TRUE) >= lim, d$combination, "")
lvl <- unique(a$level)
d$level <- factor(d$level, levels = lvl)
ticks <- -3:3
tick_lab <- c("1/8", "1/4", "1/2", "1", "2", "4", "8")
p <- ggplot(d, aes(x, y, colour = lineage)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
  geom_point(size = 1.8, alpha = 0.85, na.rm = TRUE) +
  facet_wrap(~ level) +
  scale_colour_manual(values = c("CD4+" = "#B5838D", "CD8+" = "#6D97A8",
                                 other = "#9A9A9A"), name = NULL) +
  scale_x_continuous(breaks = ticks, labels = tick_lab, limits = c(-cap, cap)) +
  scale_y_continuous(breaks = ticks, labels = tick_lab, limits = c(-cap, cap)) +
  coord_fixed() +
  labs(x = paste0("Fold change vs ", a$reference[1], ": ", opt$label_a),
       y = paste0("Fold change vs ", b$reference[1], ": ", opt$label_b),
       title = paste0("Change of every gate combination, by ", opt$column),
       subtitle = paste0("On the diagonal: the same change in both. Capped at ",
                         "8-fold. Labelled: at least ", opt$label_fold,
                         "-fold in either.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(aes(label = lab), size = 2.2,
                                    show.legend = FALSE, max.overlaps = 30,
                                    segment.colour = "grey60", segment.size = 0.2,
                                    min.segment.length = 0, na.rm = TRUE)
}
nf <- length(lvl)
nc <- min(nf, 2L)
ggsave(file.path(opt$output_dir, paste0("compare_", col_tag, "_scatter.png")), p,
       width = 4.6 * nc + 0.6, height = 4.8 * ceiling(nf / nc) + 1, dpi = 170,
       bg = "white")

# --- Side by side: the two runs' fold-change maps of each level. -------------
if (requireNamespace("png", quietly = TRUE)) {
  find_map <- function(dir, lv) {
    f <- find_in_run(dir, paste0("_diff_", col_tag, "_.*\\.png$"))
    f[endsWith(f, paste0("_diff_", col_tag, "_", safe(lv), ".png"))][1]
  }
  for (lv in lvl) {
    fa <- find_map(opt$run_a, lv)
    fb <- find_map(opt$run_b, lv)
    if (is.na(fa) || is.na(fb)) next
    ia <- png::readPNG(fa)
    ib <- png::readPNG(fb)
    panel <- function(img, title) {
      grid::gTree(children = grid::gList(
        grid::textGrob(title, y = grid::unit(1, "npc") - grid::unit(0.6, "lines"),
                       gp = grid::gpar(fontsize = 15, fontface = "bold")),
        grid::rasterGrob(img, y = 0.47, height = grid::unit(0.92, "npc"))))
    }
    w <- (ncol(ia) + ncol(ib)) / 170
    h <- max(nrow(ia), nrow(ib)) / 170 * 1.08
    out <- file.path(opt$output_dir, paste0("compare_", col_tag, "_", safe(lv), ".png"))
    grDevices::png(out, width = w, height = h, units = "in", res = 170, bg = "white")
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
    grid::pushViewport(grid::viewport(layout.pos.col = 1))
    grid::grid.draw(panel(ia, opt$label_a))
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.col = 2))
    grid::grid.draw(panel(ib, opt$label_b))
    grid::popViewport(2)
    invisible(grDevices::dev.off())
  }
} else {
  message("Package 'png' is not installed: the side-by-side maps are skipped.")
}
cat("Comparison written to ", opt$output_dir, "\n", sep = "")
