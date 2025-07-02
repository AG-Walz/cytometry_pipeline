#!/usr/bin/env Rscript

# Test script for expression heatmap creation
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(gplots)
})

# Configuration
CONFIG <- list(
  output_dir = "/home/ngr24/Desktop/Maerklin/Cytometry/Output/Split_by_cluster/",
  create_group_matrices = TRUE
)

# Test the heatmap creation function
source("Create_labelled_exprmat.R")

# Define timestamp
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

# Test heatmap creation
cat("Testing expression heatmap creation...\n")
create_group_expression_heatmaps(timestamp)
cat("=== Done ===\n")