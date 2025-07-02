#!/usr/bin/env Rscript

# Test script to verify cluster naming logic
library(data.table)

# Source the functions
source("/home/ngr24/Desktop/Maerklin/Cytometry/R_scripts/05_analysis/Create_labelled_exprmat.R")

# Test cluster names from user examples
test_names <- c(
  "AI_Cluster_BCs_128 (Clustering annotated)",
  "AI_Cluster_NKCs_12 (Clustering annotated)", 
  "AI_Cluster_CD3-_121 (Clustering annotated)",
  "AI_Cluster_NKTCs_120 (Clustering annotated)",
  "AI_Cluster_NKCs_123 (Clustering annotated)",
  "AI_Cluster_CD4TCs_82 (Clustering annotated)",
  "AI_Cluster_CD8TCs_76 (Clustering annotated)",
  "AI_Cluster_gdTCR-CD56-TCells_120 (Clustering annotated)",
  "AI_Cluster_γδTCells_60 (Clustering annotated)"
)

cat("Testing cluster name processing:\n")
cat("================================\n")

for (test_name in test_names) {
  clean_name <- clean_cluster_name(test_name)
  group_name <- extract_cluster_group(test_name)
  
  # Apply the display name logic from the script
  display_name <- clean_name
  if (grepl("^([A-Za-z0-9+-]+)_[0-9]+$", clean_name)) {
    # Already in correct format
    display_name <- clean_name
  } else if (grepl("_[0-9]+$", clean_name)) {
    # Has number but might need group prefix adjustment
    display_name <- clean_name
  } else {
    # No number, might need special handling
    display_name <- paste0(group_name, "s_", gsub("[^0-9]", "", test_name))
    if (display_name == paste0(group_name, "s_")) {
      display_name <- clean_name
    }
  }
  
  cat(sprintf("Original: %s\n", test_name))
  cat(sprintf("Cleaned:  %s\n", clean_name))
  cat(sprintf("Group:    %s\n", group_name))
  cat(sprintf("Display:  %s\n", display_name))
  cat("---\n")
}