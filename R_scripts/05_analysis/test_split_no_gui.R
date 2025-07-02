#!/usr/bin/env Rscript

# Test version of split clusters with hardcoded marker selection (no GUI)
cat("=== Testing Cluster Splitting (No GUI) ===\n")

# Source the main script
source("01_split_clusters.R")

# Override the GUI function with a hardcoded selection
select_columns_gui <- function(csv_file) {
  cat("Using hardcoded column selection for testing...\n")
  
  # Read sample to get available columns
  sample_data <- fread(csv_file, nrows = 10)
  all_columns <- colnames(sample_data)
  
  # Select common marker patterns
  marker_patterns <- c("CD3", "CD4", "CD8", "CD19", "CD56", "NK", "TCR", "HLA", "Ki67", "PD1", "LAG3", "TIM3")
  selected_markers <- all_columns[grepl(paste(marker_patterns, collapse = "|"), all_columns, ignore.case = TRUE)]
  
  # Find UMAP columns
  umap_x <- all_columns[grepl("umap.*1|umap.*x", all_columns, ignore.case = TRUE)][1]
  umap_y <- all_columns[grepl("umap.*2|umap.*y", all_columns, ignore.case = TRUE)][1]
  
  cat(sprintf("Auto-selected %d markers: %s\n", length(selected_markers), paste(head(selected_markers, 5), collapse = ", ")))
  if (!is.na(umap_x)) cat(sprintf("UMAP X: %s\n", umap_x))
  if (!is.na(umap_y)) cat(sprintf("UMAP Y: %s\n", umap_y))
  
  return(list(
    markers = selected_markers,
    umap_x = if(!is.na(umap_x)) umap_x else NULL,
    umap_y = if(!is.na(umap_y)) umap_y else NULL
  ))
}

# Test with a small subset if too many files
test_main <- function() {
  cat("=== Testing Cluster Extraction (Limited) ===\n")
  
  start_time <- Sys.time()
  
  # Validate configuration
  validate_config()
  
  # Get input files
  csv_files <- get_csv_files()
  
  # Limit to first 3 files for testing
  if (length(csv_files) > 3) {
    cat(sprintf("Limiting test to first 3 of %d files\n", length(csv_files)))
    csv_files <- csv_files[1:3]
  }
  
  # Auto column selection using first file
  cat("\n=== Automatic Column Selection ===\n")
  selected_columns <- select_columns_gui(csv_files[1])
  
  if (length(selected_columns$markers) == 0) {
    cat("No markers auto-selected. Using fallback selection.\n")
    sample_data <- fread(csv_files[1], nrows = 10)
    all_columns <- colnames(sample_data)
    selected_columns$markers <- all_columns[!grepl("^(File|Time|Event|Original_ID|cluster|umap)", all_columns, ignore.case = TRUE)][1:10]
  }
  
  cat(sprintf("Will process with %d markers\n", length(selected_columns$markers)))
  
  # Process files
  cat(sprintf("\n=== Processing %d test files ===\n", length(csv_files)))
  
  total_clusters_processed <- 0
  
  for (i in seq_along(csv_files)) {
    csv_file <- csv_files[i]
    cat(sprintf("Processing file %d/%d: %s\n", i, length(csv_files), basename(csv_file)))
    
    clusters_from_file <- process_file_by_cluster(csv_file, selected_columns)
    if (!is.null(clusters_from_file) && clusters_from_file > 0) {
      total_clusters_processed <- total_clusters_processed + clusters_from_file
      cat(sprintf("  -> Created %d cluster files\n", clusters_from_file))
    } else {
      cat("  -> No clusters processed\n")
    }
  }
  
  # Final summary
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("\n=== TEST COMPLETE ===\n")
  cat(sprintf("Files processed: %d\n", length(csv_files)))
  cat(sprintf("Total cluster files created: %d\n", total_clusters_processed))
  cat(sprintf("Processing time: %.2f minutes\n", total_time))
  cat(sprintf("Output directory: %s\n", CONFIG$output_dir))
  
  # Verify output
  output_files <- if (dir.exists(CONFIG$output_dir)) {
    list.files(CONFIG$output_dir, pattern = "\\.csv$", full.names = FALSE)
  } else {
    character(0)
  }
  
  cat(sprintf("Verification: %d CSV files in output\n", length(output_files)))
  
  if (total_clusters_processed > 0) {
    cat("✅ SUCCESS: Test cluster extraction completed!\n")
  } else {
    cat("⚠️  WARNING: No cluster files created in test\n")
  }
  
  return(invisible(total_clusters_processed))
}

# Run the test
if (!interactive()) {
  test_main()
}