#!/usr/bin/env Rscript
# Step 4, layer A: per-neg-tube marker summaries on the compensated arcsinh scale, EMD to the
# reference distribution, and CD4/CD8 gate bounds read from the workspaces.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "04_neg_summaries.R --wsp_dir DIR --fcs_dir DIR --out DIR --config FILE [--list_populations] [--all_tubes]"
a <- parse_args(list(wsp_dir = "", fcs_dir = "", out = "",
                     config = "",
                     list_populations = FALSE, all_tubes = FALSE), usage = usage)
cfg <- read_drift_config(a$config)
need_packages(c("CytoML", "flowWorkspace", "flowCore"))

tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
pm <- read.csv(file.path(a$out, "panel_map.csv"), stringsAsFactors = FALSE)
excl <- excluded_tubes(a$out)
wsps <- list_workspaces(a$wsp_dir)
sigs <- read.csv(file.path(a$out, "panel_signatures.csv"), stringsAsFactors = FALSE)

summ <- list(); emds <- list(); gates <- list(); pops_seen <- character(0)
fluor <- setdiff(c(cfg$technical_markers, cfg$cytokine_markers), c("FSC_A", "SSC_A"))

for (w in wsps) {
  for (gs in open_workspace(w, a$fcs_dir, sigs)) for (s in flowWorkspace::sampleNames(gs)) {
    gh <- gs[[s]]
    fn <- gh_filename(gh)
    row <- tt[tt$filename == fn, ]
    if (!nrow(row)) { cat("  not in tube table:", fn, "\n"); next }
    if (fn %in% excl) next
    if (a$list_populations) { pops_seen <- union(pops_seen, flowWorkspace::gh_get_pop_paths(gh, path = "full")); next }
    if (!a$all_tubes && !isTRUE(row$is_neg)) next
    node <- find_population(gh, cfg$neg_population)
    if (is.na(node)) { cat("  population", cfg$neg_population, "not found in", fn, "\n"); next }
    cm <- pm[pm$filename == fn & !is.na(pm$canonical), ]
    channel_map <- setNames(cm$channel, cm$canonical)
    m <- population_matrix(gh, node, channel_map, n_max = cfg$n_cells_summary, seed = cfg$seed)
    for (mk in intersect(fluor, colnames(m))) m[, mk] <- asinh_cof(m[, mk], cfg$arcsinh_cofactor)
    for (mk in colnames(m)) {
      s1 <- marker_summary(m[, mk], probs = unlist(cfg$quantiles), valley = mk %in% cfg$bimodal_markers)
      summ[[length(summ) + 1L]] <- cbind(filename = fn, marker = mk, s1, stringsAsFactors = FALSE)
    }
    # keep the (subsampled, arcsinh) cells per tube for EMD and for BayesFlow variant 2
    dir.create(file.path(a$out, "neg_cells"), showWarnings = FALSE)
    con <- gzfile(file.path(a$out, "neg_cells", paste0(fn, ".csv.gz")), "w")
    write.csv(round(m, 4), con, row.names = FALSE); close(con)
    # gate bounds of the parent populations
    for (pp in names(cfg$parent_populations)) {
      gnode <- find_population(gh, cfg$parent_populations[[pp]])
      if (is.na(gnode)) next
      gb <- tryCatch(gate_bounds(gh, gnode), error = function(e) NULL)
      if (!is.null(gb)) {
        gb$parent <- pp; gb$filename <- fn
        # map dimension channel -> canonical and arcsinh the bounds of fluorescence dims
        gb$marker <- names(channel_map)[match(gb$dimension, channel_map)]
        isf <- gb$marker %in% fluor
        gb$lower[isf] <- asinh_cof(gb$lower[isf], cfg$arcsinh_cofactor)
        gb$upper[isf] <- asinh_cof(gb$upper[isf], cfg$arcsinh_cofactor)
        gates[[length(gates) + 1L]] <- gb
      }
    }
  }
}

if (a$list_populations) { cat("\nPopulations seen:\n"); cat(sort(unique(pops_seen)), sep = "\n"); quit(status = 0) }

summary_df <- do.call(rbind, summ)
gates_df <- if (length(gates)) do.call(rbind, gates) else data.frame()
write_table(summary_df, file.path(a$out, "neg_summaries.csv"))
write_table(gates_df, file.path(a$out, "gate_bounds.csv"))

# EMD of every neg tube to pooled reference-instrument negs: those of the same patient
# (excluding itself) when available, else all of them.
cells <- list.files(file.path(a$out, "neg_cells"), "\\.csv\\.gz$", full.names = TRUE)
names(cells) <- sub("\\.csv\\.gz$", "", basename(cells))
mats <- lapply(cells, function(f) as.matrix(read.csv(gzfile(f), check.names = FALSE)))
ref_files <- tt$filename[tt$instrument == cfg$reference_instrument & tt$filename %in% names(mats)]
pool <- function(fns, mk) unlist(lapply(fns, function(f) if (mk %in% colnames(mats[[f]])) mats[[f]][, mk]))
for (fn in names(mats)) {
  pat <- tt$patient[tt$filename == fn]
  same_pat <- setdiff(ref_files[tt$patient[match(ref_files, tt$filename)] == pat], fn)
  ref_kind <- if (length(same_pat)) "same_patient_reference_instrument" else "pooled_reference_instrument"
  ref <- if (length(same_pat)) same_pat else setdiff(ref_files, fn)
  for (mk in colnames(mats[[fn]])) {
    emds[[length(emds) + 1L]] <- data.frame(filename = fn, marker = mk, ref_kind = ref_kind,
      emd = emd_1d(mats[[fn]][, mk], pool(ref, mk)), stringsAsFactors = FALSE)
  }
}
write_table(do.call(rbind, emds), file.path(a$out, "neg_emd.csv"))
