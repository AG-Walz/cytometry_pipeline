#!/usr/bin/env Rscript

# Load the main script functions
source("Create_labelled_exprmat.R")

cat("=== Running Expression Matrix Creation ===\n")
cat("CONFIG create_group_matrices:", CONFIG$create_group_matrices, "\n")

if (CONFIG$create_group_matrices) {
  cat("Creating group expression matrices...\n")
  
  # Create timestamp like in main script
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  cat("Using timestamp:", timestamp, "\n")
  
  # Call the actual function from the script with timestamp
  result <- tryCatch({
    create_group_expression_matrices(timestamp)
    cat("Function completed successfully\n")
  }, error = function(e) {
    cat("Error in function:", e$message, "\n")
    cat("Traceback:\n")
    print(traceback())
  })
  
  # Check if the output directory was created and has files
  matrices_dir <- file.path(dirname("/home/ngr24/Desktop/Maerklin/Cytometry/Output/Split_by_cluster/"), "group_expression_matrices")
  
  if (dir.exists(matrices_dir)) {
    cat("Matrices directory exists\n")
    files <- list.files(matrices_dir, recursive = TRUE)
    cat("Files in matrices directory:", length(files), "\n")
    if (length(files) > 0) {
      cat("Sample files:\n")
      cat(head(files, 10), sep = "\n")
    }
  } else {
    cat("Matrices directory does NOT exist\n")
  }
} else {
  cat("create_group_matrices is FALSE, skipping\n")
}

cat("=== Done ===\n")