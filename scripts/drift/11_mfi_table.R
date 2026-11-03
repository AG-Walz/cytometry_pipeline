#!/usr/bin/env Rscript
# Step 11: per-tube intensities at several gating levels, for the variance partition
# (docs/plans/drift_three_questions_plan.md, part D, and the gating comparison of 2026-10-13).
# For every tube, stimulated and unstimulated:
#   - the median of each technical marker within the neg population ("live"), as before;
#   - the median of each cytokine marker on all events ("all", no gate), within the neg population ("live")
#     and within each parent population ("CD4", "CD8"), so the effect of gating on the cytokine MFI shows.
# Medians on the compensated linear scale; the arcsinh median is asinh(median / cofactor) and the FlowJo median
# the median through that sample's own FlowJo display transformation as stored in the workspace (biexponential
# here), a median being unchanged by a monotone transform. Writes mfi.csv (filename, population, marker, n,
# median_linear, median_asinh, median_flowjo).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "11_mfi_table.R --wsp_dir DIR --fcs_dir DIR --out DIR --config FILE"
a <- parse_args(list(wsp_dir = "", fcs_dir = "", out = "", config = ""), usage = usage)
cfg <- read_drift_config(a$config)
need_packages(c("CytoML", "flowWorkspace", "flowCore"))

tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
pm <- read.csv(file.path(a$out, "panel_map.csv"), stringsAsFactors = FALSE)
sigs <- read.csv(file.path(a$out, "panel_signatures.csv"), stringsAsFactors = FALSE)
excl <- excluded_tubes(a$out)
fluor <- setdiff(c(cfg$technical_markers, cfg$cytokine_markers), c("FSC_A", "SSC_A"))
cyto <- unlist(cfg$cytokine_markers)
# gating level -> FlowJo node suffix, label, markers summarised
plan <- c(list(list(node = "root", label = "all", markers = cyto),
               list(node = cfg$neg_population, label = "live", markers = c(unlist(cfg$technical_markers), cyto))),
          lapply(names(cfg$parent_populations), function(p) list(node = cfg$parent_populations[[p]], label = p, markers = cyto)))

rows <- list()
for (w in list_workspaces(a$wsp_dir)) {
  for (gs in open_workspace(w, a$fcs_dir, sigs)) for (s in flowWorkspace::sampleNames(gs)) {
    gh <- gs[[s]]; fn <- gh_filename(gh)
    if (!fn %in% tt$filename || fn %in% excl) next
    cm <- pm[pm$filename == fn & !is.na(pm$canonical), ]
    channel_map <- setNames(cm$channel, cm$canonical)
    fj <- flowWorkspace::gh_get_transformations(gh)   # this sample's FlowJo transformation per channel
    names(fj) <- sub("^Comp-", "", names(fj))
    for (pl in plan) {
      node <- find_population(gh, pl$node)
      if (is.na(node)) { cat("  population", pl$node, "not found in", fn, "\n"); next }
      m <- population_matrix(gh, node, channel_map, n_max = .Machine$integer.max)
      for (mk in intersect(pl$markers, colnames(m))) {
        x <- m[, mk]; x <- x[is.finite(x)]
        med <- if (length(x)) median(x) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(filename = fn, population = pl$label, marker = mk, n = length(x),
          median_linear = med, median_asinh = if (mk %in% fluor) asinh(med / cfg$arcsinh_cofactor) else med,
          median_flowjo = if (channel_map[[mk]] %in% names(fj) && is.finite(med)) fj[[channel_map[[mk]]]](med) else med,
          stringsAsFactors = FALSE)
      }
    }
  }
}
mfi <- do.call(rbind, rows)
write_table(mfi, file.path(a$out, "mfi.csv"))
cat("\nTubes per population:\n"); print(tapply(mfi$filename, mfi$population, function(x) length(unique(x))))
