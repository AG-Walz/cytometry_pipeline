#!/usr/bin/env Rscript
# Step 11b: acquisition audit. What the cytometer and FlowJo did to every tube, and how these settings overlap
# with the design factors (src/harmonization/drift/drift_audit.R). Writes:
#   fcs_keywords.csv.gz          every FCS keyword of every tube (long format)
#   acquisition_settings.csv     detector voltage per marker, voltage set, compensation matrix id, CS&T setup and
#                                bead lot, CS&T performance expired, cytometer configuration, experiment name
#   sample_quality.csv           viability (live-cell gate / its parent) and live-cell count per tube
#   flowjo_transformations.csv   FlowJo display transformation per sample and parameter, from the workspaces
#   gate_positions.csv           every gate of every sample, on FlowJo's scale and in raw fluorescence units
#   factor_overlap.csv           how much each factor determines each other one (1 = fully, 0 = independent)
# Settings that one acquisition day fully determines cannot be separate terms of the variance partition; they
# explain why days differ and are reported, not modelled. Sample quality enters step 12 at the visit level.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
usage <- "11b_acquisition_audit.R --out DIR --config FILE [--wsp_dir DIR --fcs_dir DIR]"
a <- parse_args(list(out = "", config = "", wsp_dir = "none", fcs_dir = "none"), usage = usage)
cfg <- read_drift_config(a$config)
need_packages("flowCore")
tt <- read.csv(file.path(a$out, "tubes.csv"), stringsAsFactors = FALSE)
pm <- read.csv(file.path(a$out, "panel_map.csv"), stringsAsFactors = FALSE)

# FCS keywords and acquisition settings
kw <- do.call(rbind, lapply(seq_len(nrow(tt)), function(i) {
  k <- flowCore::read.FCSheader(tt$path[i])[[1]]
  data.frame(filename = tt$filename[i], key = toupper(names(k)), value = unname(k), stringsAsFactors = FALSE)
}))
con <- gzfile(file.path(a$out, "fcs_keywords.csv.gz"), "w"); write.csv(kw, con, row.names = FALSE); close(con)
acq <- acquisition_settings(kw, pm)
write_table(acq, file.path(a$out, "acquisition_settings.csv"))

# sample quality from the FlowJo gate counts of step 5
gc <- file.path(a$out, "all_gate_counts.csv")
if (file.exists(gc)) {
  q <- sample_quality(read.csv(gc, stringsAsFactors = FALSE), cfg$neg_population)
  write_table(q, file.path(a$out, "sample_quality.csv"))
} else cat("no all_gate_counts.csv (step 5 not run): no sample quality\n")

# FlowJo transformations and gate positions
if (a$wsp_dir != "none") {
  wsps <- list_workspaces(a$wsp_dir)
  write_table(do.call(rbind, lapply(wsps, flowjo_transformations)), file.path(a$out, "flowjo_transformations.csv"))
  if (a$fcs_dir != "none" && requireNamespace("CytoML", quietly = TRUE)) {
    sigs <- read.csv(file.path(a$out, "panel_signatures.csv"), stringsAsFactors = FALSE)
    rows <- list()
    for (w in wsps) for (gs in open_workspace(w, a$fcs_dir, sigs, quiet = TRUE)) for (s in flowWorkspace::sampleNames(gs)) {
      gh <- gs[[s]]; fn <- gh_filename(gh); inv <- flowWorkspace::gh_get_transformations(gh, inverse = TRUE)
      for (node in setdiff(flowWorkspace::gh_get_pop_paths(gh, path = "full"), "root")) {
        g <- tryCatch(flowWorkspace::gh_pop_get_gate(gh, node), error = function(e) NULL)
        if (is.null(g) || methods::is(g, "booleanFilter")) next
        dims <- names(flowCore::parameters(g))
        if (methods::is(g, "rectangleGate")) { lo <- g@min; hi <- g@max }
        else if (methods::is(g, "polygonGate")) { lo <- apply(g@boundaries, 2, min); hi <- apply(g@boundaries, 2, max) }
        else next
        back <- function(v, d) if (!is.finite(v)) NA_real_ else if (!is.null(inv[[d]])) inv[[d]](v) else v
        rows[[length(rows) + 1L]] <- data.frame(filename = fn, gate = node, gate_type = class(g)[1], dimension = sub("^Comp-", "", dims),
          lower_flowjo = unname(lo), upper_flowjo = unname(hi),
          lower_raw = vapply(seq_along(dims), function(i) back(unname(lo[i]), dims[i]), numeric(1)),
          upper_raw = vapply(seq_along(dims), function(i) back(unname(hi[i]), dims[i]), numeric(1)), stringsAsFactors = FALSE)
      }
    }
    write_table(do.call(rbind, rows), file.path(a$out, "gate_positions.csv"))
  }
}

# overlap of design factors and acquisition settings, on the stimulated tubes
d <- merge(tt, acq[, intersect(c("filename", "voltage_set", "compensation_id", "cst_setup", "experiment"), names(acq))], by = "filename")
ds <- cfg$design
if (!is.null(ds$time)) d$visit <- paste(d[[ds$unit]], d[[ds$time]], sep = "_")
cols <- c(day = design_day_term(cfg), setNames(as.list(design_technical_factors(cfg)), design_technical_factors(cfg)),
          voltages = "voltage_set", compensation = "compensation_id", "CS&T setup" = "cst_setup",
          setNames(ds$unit, ds$unit), visit = if (!is.null(ds$time)) "visit", time = ds$time, pool = ds$stimulation)
cols <- cols[vapply(cols, function(x) length(x) == 1 && x %in% names(d) && length(unique(stats::na.omit(d[[x]]))) > 1, logical(1))]
M <- determination_matrix(d[!d$is_neg %in% TRUE, ], unlist(cols))
write.csv(round(M, 3), file.path(a$out, "factor_overlap.csv"))
cat("\nHow much the row factor determines the column factor (stimulated tubes):\n"); print(round(M, 2))
if (exists("q")) {
  dq <- merge(d, q, by = "filename")
  E <- sapply(unlist(cols), function(f) eta_squared(qlogis(pmin(pmax(dq$viability, .01), .99)), dq[[f]]))
  cat("\nViability explained by each factor alone (eta squared):\n"); print(round(E, 2))
}
