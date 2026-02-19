#!/usr/bin/env Rscript
# =============================================================================
# test_data/generate_test_data.R — Generate Synthetic Test Data for Pipelines
# =============================================================================
# Creates structured test fixtures for Pipelines A, B, C, E.
# Run from the project root: Rscript test_data/generate_test_data.R
# =============================================================================

cat("=== Generating Synthetic Test Data ===\n")

# Determine script location
.script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) {
    # Fallback: search the call stack
    for (i in rev(seq_len(sys.nframe()))) {
      f <- sys.frame(i)$ofile
      if (!is.null(f)) return(normalizePath(f))
    }
    # Last resort: assume run from project root
    normalizePath("test_data/generate_test_data.R")
  }
)
base_dir <- dirname(.script_path)
if (basename(base_dir) != "test_data") {
  base_dir <- file.path(base_dir, "test_data")
}

markers <- c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR")
sample_names <- c("sample_001", "sample_002", "sample_003")
n_cells <- 50
set.seed(42)

# =============================================================================
# Pipeline A & E: Synthetic FCS files
# =============================================================================
fcs_dir <- file.path(base_dir, "fcs_samples")
dir.create(fcs_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n--- Generating FCS files ---\n")

for (i in seq_along(sample_names)) {
  # Generate expression data with some structure
  # Sample 1-2: "Control" pattern, Sample 3: "Treatment" pattern
  expr_data <- matrix(rnorm(n_cells * length(markers), mean = 3, sd = 1),
                      nrow = n_cells, ncol = length(markers))
  colnames(expr_data) <- markers

  # Add some biological signal: CD3+ cells in treatment are shifted
  if (i == 3) {
    expr_data[1:25, "CD3"] <- expr_data[1:25, "CD3"] + 2
    expr_data[1:25, "CD4"] <- expr_data[1:25, "CD4"] + 1.5
  }

  # Ensure non-negative values (flow data)
  expr_data <- abs(expr_data)

  # Create flowFrame
  params <- Biobase::AnnotatedDataFrame(
    data.frame(
      name = markers,
      desc = markers,
      range = rep(1024, length(markers)),
      minRange = rep(0, length(markers)),
      maxRange = rep(1024, length(markers)),
      row.names = paste0("$P", seq_along(markers)),
      stringsAsFactors = FALSE
    )
  )
  ff <- flowCore::flowFrame(exprs = expr_data, parameters = params)

  fcs_path <- file.path(fcs_dir, paste0(sample_names[i], ".fcs"))
  flowCore::write.FCS(ff, fcs_path)
  cat("  Created:", fcs_path, "(", n_cells, "cells x", length(markers), "markers)\n")
}

# Conditions CSV for FCS pipelines
conditions_fcs <- data.frame(
  file_name = sample_names,
  condition = c("Control", "Control", "Treatment"),
  stringsAsFactors = FALSE
)
cond_path <- file.path(fcs_dir, "conditions.csv")
data.table::fwrite(conditions_fcs, cond_path, sep = ";")
cat("  Created:", cond_path, "\n")

# =============================================================================
# Pipeline B: Per-sample expression CSVs
# =============================================================================
csv_expr_dir <- file.path(base_dir, "csv_expression")
dir.create(csv_expr_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n--- Generating Expression CSVs ---\n")

# Pipeline B marker_statistics expects filenames like {sample}_{population}.csv
populations <- c("CD4_T_cells", "CD8_T_cells")

for (sname in sample_names) {
  for (pop in populations) {
    expr_data <- matrix(abs(rnorm(n_cells * length(markers), mean = 3, sd = 1)),
                        nrow = n_cells, ncol = length(markers))
    colnames(expr_data) <- markers

    csv_path <- file.path(csv_expr_dir, paste0(sname, "_", pop, ".csv"))
    data.table::fwrite(as.data.frame(expr_data), csv_path, sep = ";")
    cat("  Created:", csv_path, "\n")
  }
}

# =============================================================================
# Pipeline C: CSV gating + dummy FCS files
# =============================================================================
csv_gating_dir <- file.path(base_dir, "csv_gating")
dir.create(csv_gating_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n--- Generating Gating CSVs + FCS pairs ---\n")

for (i in seq_along(sample_names)) {
  sname <- sample_names[i]

  # Expression data
  expr_data <- matrix(abs(rnorm(n_cells * length(markers), mean = 3, sd = 1)),
                      nrow = n_cells, ncol = length(markers))
  colnames(expr_data) <- markers

  # Binary gating columns
  gate1 <- rep(0L, n_cells)
  gate1[sample(n_cells, 20)] <- 1L  # 20 cells in cluster 1
  gate2 <- rep(0L, n_cells)
  gate2[sample(n_cells, 15)] <- 1L  # 15 cells in cluster 2
  raw_cluster <- rep(0L, n_cells)
  raw_cluster[sample(n_cells, 25)] <- 1L

  csv_data <- cbind(
    as.data.frame(expr_data),
    `CD4+ (Clustering annotated)` = gate1,
    `CD8+ (Clustering annotated)` = gate2,
    `Clustering raw` = raw_cluster
  )

  csv_path <- file.path(csv_gating_dir, paste0(sname, ".csv"))
  data.table::fwrite(csv_data, csv_path, sep = ";")
  cat("  Created:", csv_path, "\n")

  # Matching dummy FCS file (needed by load_csv_gating for validation)
  params <- Biobase::AnnotatedDataFrame(
    data.frame(
      name = markers,
      desc = markers,
      range = rep(1024, length(markers)),
      minRange = rep(0, length(markers)),
      maxRange = rep(1024, length(markers)),
      row.names = paste0("$P", seq_along(markers)),
      stringsAsFactors = FALSE
    )
  )
  ff <- flowCore::flowFrame(
    exprs = as.matrix(expr_data),
    parameters = params
  )
  fcs_path <- file.path(csv_gating_dir, paste0(sname, ".fcs"))
  flowCore::write.FCS(ff, fcs_path)
  cat("  Created:", fcs_path, "\n")
}

# Conditions CSV for gating pipeline
conditions_gating <- data.frame(
  file_name = sample_names,
  condition = c("Control", "Control", "Treatment"),
  stringsAsFactors = FALSE
)
cond_gating_path <- file.path(base_dir, "csv_gating_conditions.csv")
data.table::fwrite(conditions_gating, cond_gating_path, sep = ";")
cat("  Created:", cond_gating_path, "\n")

cat("\n=== Test Data Generation Complete ===\n")
cat("  FCS samples:        ", fcs_dir, "\n")
cat("  Expression CSVs:    ", csv_expr_dir, "\n")
cat("  Gating CSVs + FCS:  ", csv_gating_dir, "\n")
