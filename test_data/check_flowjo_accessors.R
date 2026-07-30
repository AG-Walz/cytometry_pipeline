#!/usr/bin/env Rscript
# Which FlowJo workspace accessors does THIS install actually export?
suppressPackageStartupMessages({library(CytoML); library(flowWorkspace)})
cat("CytoML ", as.character(packageVersion("CytoML")),
    " | flowWorkspace ", as.character(packageVersion("flowWorkspace")), "\n\n", sep = "")
for (n in c("fj_ws_get_sample_groups", "getSampleGroups",
            "fj_ws_get_samples", "getSamples", "open_flowjo_xml",
            "flowjo_to_gatingset")) {
  where <- character(0)
  for (p in c("CytoML", "flowWorkspace")) {
    ok <- !is.null(tryCatch(getExportedValue(p, n), error = function(e) NULL))
    if (ok) where <- c(where, p)
  }
  cat(sprintf("  %-26s %s\n", n,
              if (length(where)) paste(where, collapse = " + ") else "-- NOT FOUND --"))
}
cat("\nflowjo_to_gatingset() arguments:\n  ",
    paste(names(formals(CytoML::flowjo_to_gatingset)), collapse = ", "), "\n", sep = "")
