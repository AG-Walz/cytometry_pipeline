#!/usr/bin/env Rscript

# Test expression matrix function
cat("=== Testing Expression Matrix Function ===\n")

# Source the main script to get functions
source("Create_labelled_exprmat.R")

# Check configuration
cat("CONFIG create_group_matrices:", CONFIG$create_group_matrices, "\n")

# Find cluster files
output_dir <- "/data/Cytometry/ICK/ICK_T1_P+_poR/Output/Split_by_cluster"
files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("cluster_tracking_matrix", files)]

cat("Found", length(files), "cluster files\n")

# Test if function exists
if (exists("create_group_expression_matrices")) {
  cat("Function create_group_expression_matrices exists\n")
  
  # Test the function directly with some debug output
  cat("Calling create_group_expression_matrices...\n")
  
  # Create a simplified version with debug output
  create_group_expression_matrices_debug <- function() {
    cat("Starting create_group_expression_matrices_debug...\n")
    
    # Check if we should create matrices
    if (!CONFIG$create_group_matrices) {
      cat("CONFIG$create_group_matrices is FALSE, skipping\n")
      return()
    }
    
    cat("CONFIG check passed\n")
    
    # Create matrices directory
    matrices_dir <- file.path(dirname(output_dir), "group_expression_matrices")
    cat("Creating matrices directory:", matrices_dir, "\n")
    
    if (!dir.exists(matrices_dir)) {
      dir.create(matrices_dir, recursive = TRUE)
      cat("Directory created successfully\n")
    } else {
      cat("Directory already exists\n")
    }
    
    # Get cluster files
    cluster_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
    cluster_files <- cluster_files[!grepl("cluster_tracking_matrix", cluster_files)]
    
    cat("Found", length(cluster_files), "cluster files to process\n")
    
    if (length(cluster_files) == 0) {
      cat("No cluster files found!\n")
      return()
    }
    
    # Show first few files
    cat("Sample files:\n")
    cat(head(basename(cluster_files), 5), sep = "\n")
    
    # Try to read one file to test
    test_file <- cluster_files[1]
    cat("Testing read of file:", basename(test_file), "\n")
    
    tryCatch({
      test_data <- fread(test_file, nrows = 5)
      cat("Successfully read test file with", ncol(test_data), "columns\n")
      cat("Column names:", paste(head(colnames(test_data), 10), collapse = ", "), "\n")
    }, error = function(e) {
      cat("Error reading test file:", e$message, "\n")
    })
    
    cat("Debug function completed\n")
  }
  
  # Run debug version
  create_group_expression_matrices_debug()
  
} else {
  cat("Function create_group_expression_matrices does NOT exist!\n")
}

cat("=== Test Complete ===\n")