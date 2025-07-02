#!/usr/bin/env Rscript

# Test script to create expression heatmaps with proper cluster naming
library(data.table)
library(gplots)
library(RColorBrewer)

# Source the main script to get functions
source("/home/ngr24/Desktop/Maerklin/Cytometry/R_scripts/05_analysis/Create_labelled_exprmat.R")

# Set directory paths
output_dir <- "/home/ngr24/Desktop/Maerklin/Cytometry/Output/Split_by_cluster"
matrices_dir <- file.path(output_dir, "group_expression_matrices")
heatmaps_dir <- file.path(output_dir, "group_expression_heatmaps")

# Create directories
dir.create(heatmaps_dir, recursive = TRUE, showWarnings = FALSE)

cat("Testing expression heatmap generation with proper cluster naming...\n")

# Call the function directly to test cluster naming
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
create_group_expression_heatmaps(timestamp)

cat("Heatmap generation complete.\n")