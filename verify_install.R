#!/usr/bin/env Rscript
#
# verify_install.R — Check that all required packages are installed
#
# Usage:
#   Rscript verify_install.R
#
# Prints a summary of installed/missing packages, grouped by category.
# Exit code 0 if all required packages are present, 1 otherwise.

cat("=== Cytometry Pipeline — Installation Verification ===\n\n")

# --- R version check ---
r_ver <- paste0(R.version$major, ".", R.version$minor)
cat(sprintf("R version: %s", r_ver))
if (numeric_version(r_ver) < "4.3") {
  cat("  [WARNING: R >= 4.3 recommended]\n")
} else {
  cat("  [OK]\n")
}
cat("\n")

# --- Package lists ---
required_bioc <- c(
  "flowCore", "flowWorkspace", "CytoML", "openCyto",
  "flowAI", "PeacoQC", "FlowSOM", "flowVS",
  "preprocessCore", "flowTrans", "Biobase", "BiocGenerics",
  "BiocParallel", "cytolib", "limma", "diffcyt",
  "flowDensity", "sva", "CytoNorm"
)

required_cran <- c(
  "data.table", "ggplot2", "patchwork", "viridis", "ggrepel",
  "pheatmap", "uwot", "Rtsne", "leiden", "igraph", "Matrix",
  "lme4", "lmerTest", "emmeans", "MASS", "dbscan", "mclust",
  "survival", "survminer", "harmony",
  "tidyverse", "future", "future.apply", "progress",
  "ggridges", "ggpubr", "scales", "openxlsx", "jsonlite",
  "optparse", "reshape2", "diptest", "ggcorrplot",
  "RColorBrewer", "UpSetR", "VennDiagram", "eulerr",
  "nnls", "devtools", "reticulate"
)

optional_cran <- c(
  "ggforce", "performance", "destiny"
)

optional_github <- c(
  "Rphenograph"
)

# --- Check function ---
check_packages <- function(pkgs, label) {
  cat(sprintf("--- %s ---\n", label))
  installed <- vapply(pkgs, function(p) {
    requireNamespace(p, quietly = TRUE)
  }, logical(1))

  for (i in seq_along(pkgs)) {
    status <- if (installed[i]) "[OK]" else "[MISSING]"
    cat(sprintf("  %-25s %s\n", pkgs[i], status))
  }
  cat("\n")
  return(installed)
}

# --- Run checks ---
bioc_ok   <- check_packages(required_bioc, "Bioconductor (required)")
cran_ok   <- check_packages(required_cran, "CRAN (required)")
opt_cran  <- check_packages(optional_cran, "CRAN (optional)")
opt_gh    <- check_packages(optional_github, "GitHub (optional)")

# --- Python / leidenalg check ---
cat("--- Python (for Leiden clustering) ---\n")
python_ok <- FALSE
tryCatch({
  if (requireNamespace("reticulate", quietly = TRUE)) {
    py_config <- reticulate::py_config()
    cat(sprintf("  Python:     %s  [OK]\n", py_config$python))
    leidenalg_ok <- reticulate::py_module_available("leidenalg")
    cat(sprintf("  leidenalg:  %s\n", if (leidenalg_ok) "[OK]" else "[MISSING]"))
    igraph_ok <- reticulate::py_module_available("igraph")
    cat(sprintf("  igraph:     %s\n", if (igraph_ok) "[OK]" else "[MISSING]"))
    python_ok <- leidenalg_ok && igraph_ok
  } else {
    cat("  reticulate package not installed  [MISSING]\n")
  }
}, error = function(e) {
  cat(sprintf("  Python check failed: %s\n", conditionMessage(e)))
})
cat("\n")

# --- Summary ---
n_required_missing <- sum(!bioc_ok) + sum(!cran_ok)
n_optional_missing <- sum(!opt_cran) + sum(!opt_gh)

cat("=== Summary ===\n")
cat(sprintf("Required packages: %d / %d installed\n",
            length(bioc_ok) + length(cran_ok) - n_required_missing,
            length(bioc_ok) + length(cran_ok)))
cat(sprintf("Optional packages: %d / %d installed\n",
            length(opt_cran) + length(opt_gh) - n_optional_missing,
            length(opt_cran) + length(opt_gh)))
cat(sprintf("Python/leidenalg:  %s\n", if (python_ok) "OK" else "INCOMPLETE"))

if (n_required_missing > 0) {
  cat(sprintf("\n[!] %d required package(s) missing. Install them before running pipelines.\n", n_required_missing))
  missing_bioc <- required_bioc[!bioc_ok]
  missing_cran <- required_cran[!cran_ok]
  if (length(missing_bioc) > 0) {
    cat(sprintf("\n  BiocManager::install(c(%s))\n",
                paste0("'", missing_bioc, "'", collapse = ", ")))
  }
  if (length(missing_cran) > 0) {
    cat(sprintf("\n  install.packages(c(%s))\n",
                paste0("'", missing_cran, "'", collapse = ", ")))
  }
  quit(status = 1)
} else {
  cat("\nAll required packages are installed.\n")
  if (!python_ok) {
    cat("\nNote: Python/leidenalg is missing. Leiden clustering will not work.\n")
    cat("  Fix: pip install leidenalg igraph numpy\n")
  }
}
