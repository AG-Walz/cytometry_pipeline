
library(bigstatsr)
library(ClusterR)
results_path <- "/data/Maerklin/ICK/ICK_T1_PR/Results/"
target_rows  <- 500000      # desired number of representative rows
input_rds   <- file.path(results_path, "concatenated_data.rds")

# --- 1) Attach FBM ---
fbm <- bigstatsr::big_attach(input_rds)
n <- nrow(fbm); p <- ncol(fbm)
cat("Concatenated data dimensions:", n, "cells x", p, "markers\n")

