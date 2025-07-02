#!/usr/bin/env Rscript

# Test script for expression heatmap creation only
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(gplots)
  library(stringr)
})

# Configuration
CONFIG <- list(
  output_dir = "/home/ngr24/Desktop/Maerklin/Cytometry/Output/Split_by_cluster/",
  create_group_matrices = TRUE
)

# Source required functions only
cat("Loading required functions...\n")

# Clean cluster name function
clean_cluster_name <- function(cluster_name) {
  cleaned <- gsub("AI_Cluster_", "", cluster_name)
  cleaned <- gsub(" \\(Clustering annotated\\)", "", cleaned)
  return(cleaned)
}

# Extract cluster group function
extract_cluster_group <- function(cluster_name) {
  # First clean the name
  clean_name <- clean_cluster_name(cluster_name)
  
  # Extract group patterns
  if (grepl("BCs", clean_name)) return("BCs")
  if (grepl("CD3-", clean_name)) return("CD3-")
  if (grepl("CD4TCs", clean_name)) return("CD4TCs")
  if (grepl("CD8TCs", clean_name)) return("CD8TCs")
  if (grepl("NKCs", clean_name)) return("NKCs")
  if (grepl("NKTCells", clean_name)) return("NKTCells")
  if (grepl("γδTCells|gdTCR-CD56-TCells", clean_name)) return("gdTCR-CD56-TCells")
  
  # Default: return the name itself if no pattern matches
  return(clean_name)
}

# Test the heatmap creation function
create_group_expression_heatmaps <- function(timestamp) {
  if (!CONFIG$create_group_matrices) {
    return(invisible(NULL))
  }
  
  cat("\n=== Creating expression heatmaps per cluster group ===\n")
  
  # Create subdirectory for group heatmaps
  heatmaps_dir <- file.path(CONFIG$output_dir, "group_expression_heatmaps")
  if (!dir.exists(heatmaps_dir)) {
    dir.create(heatmaps_dir, recursive = TRUE)
  }
  
  # Get all individual cluster output files (not the aggregated matrices)
  output_files <- list.files(CONFIG$output_dir, pattern = "\\.csv$", full.names = TRUE)
  output_files <- output_files[!grepl("cluster_tracking_matrix|group_matrices_summary", output_files)]
  
  if (length(output_files) == 0) {
    warning("No cluster output files found")
    return(invisible(NULL))
  }
  
  cat(sprintf("Processing %d cluster files for heatmap creation...\n", length(output_files)))
  
  # Initialize list to collect data by group
  group_cluster_data <- list()
  
  # Process each cluster file to collect mean expressions by cluster and group
  for (file_path in output_files) {  # Process all files
    file_basename <- basename(file_path)
    
    # Extract cluster name from filename
    parts <- strsplit(gsub("\\.csv$", "", file_basename), "_")[[1]]
    if (length(parts) < 2) next
    
    cluster_name <- paste(parts[2:length(parts)], collapse = "_")
    
    # Determine group
    group_name <- extract_cluster_group(cluster_name)
    
    # Clean cluster name for display
    clean_cluster <- clean_cluster_name(cluster_name)
    
    tryCatch({
      # Read cluster data
      cluster_data <- fread(file_path)
      
      if (nrow(cluster_data) == 0) next
      
      # Get marker columns (exclude UMAP coordinates if present)
      marker_cols <- setdiff(names(cluster_data), c("UMAP_1", "UMAP_2", "umap_x", "umap_y"))
      
      if (length(marker_cols) == 0) next
      
      # Calculate mean expression for this cluster
      cluster_means <- cluster_data[, lapply(.SD, mean, na.rm = TRUE), .SDcols = marker_cols]
      cluster_means$cluster_name <- clean_cluster
      cluster_means$original_cluster <- cluster_name
      cluster_means$file_origin <- gsub("_[^_]+\\.csv$", "", file_basename)
      
      # Add to group collection
      if (group_name %in% names(group_cluster_data)) {
        group_cluster_data[[group_name]] <- rbind(group_cluster_data[[group_name]], cluster_means, fill = TRUE)
      } else {
        group_cluster_data[[group_name]] <- cluster_means
      }
      
    }, error = function(e) {
      warning(sprintf("Error processing file %s: %s", file_basename, e$message))
    })
  }
  
  if (length(group_cluster_data) == 0) {
    warning("No valid group data found for heatmaps")
    return(invisible(NULL))
  }
  
  cat(sprintf("Creating heatmaps for %d cluster groups...\n", length(group_cluster_data)))
  
  # Create heatmaps for each group
  for (group_name in names(group_cluster_data)) {
    group_data <- group_cluster_data[[group_name]]
    
    # Skip if too few clusters
    if (nrow(group_data) < 2) {
      cat(sprintf("  Skipping %s (only %d cluster)\n", group_name, nrow(group_data)))
      next
    }
    
    tryCatch({
      # Get marker columns for heatmap
      marker_cols <- setdiff(names(group_data), c("cluster_name", "original_cluster", "file_origin"))
      
      if (length(marker_cols) == 0) {
        cat(sprintf("  Skipping %s (no marker columns)\n", group_name))
        next
      }
      
      # Create expression matrix
      expr_matrix <- as.matrix(group_data[, marker_cols, with = FALSE])
      rownames(expr_matrix) <- group_data$cluster_name
      
      # Handle duplicate cluster names by making them unique
      if (any(duplicated(rownames(expr_matrix)))) {
        rownames(expr_matrix) <- make.unique(rownames(expr_matrix))
      }
      
      # Create heatmap
      safe_group_name <- gsub("[^A-Za-z0-9_-]", "_", group_name)
      heatmap_filename <- sprintf("expression_heatmap_%s_%s.png", safe_group_name, timestamp)
      heatmap_path <- file.path(heatmaps_dir, heatmap_filename)
      
      # Calculate appropriate dimensions
      n_clusters <- nrow(expr_matrix)
      n_markers <- ncol(expr_matrix)
      
      # Dynamic sizing
      base_width <- max(10, min(24, n_markers * 0.6 + 5))
      base_height <- max(8, min(20, n_clusters * 0.4 + 4))
      
      png(heatmap_path, width = base_width, height = base_height, units = "in", res = 300)
      
      # Create color palette
      colors <- colorRampPalette(c("blue", "white", "red"))(100)
      
      # Create heatmap with clustering
      tryCatch({
        heatmap.2(expr_matrix,
                 main = sprintf("Marker Expression - %s\n(%d clusters across files)", 
                               group_name, nrow(expr_matrix)),
                 col = colors,
                 scale = "column",  # Scale by marker (column)
                 trace = "none",
                 density.info = "none",
                 margins = c(15, 18),
                 cexRow = max(0.5, min(1.0, 15/n_clusters)),
                 cexCol = max(0.5, min(1.0, 12/n_markers)),
                 key = TRUE,
                 keysize = 1.0,
                 key.title = "Scaled\nExpression",
                 key.xlab = "",
                 symkey = FALSE,
                 symm = FALSE,
                 distfun = function(x) dist(x, method = "euclidean"),
                 hclustfun = function(x) hclust(x, method = "ward.D2"),
                 Rowv = TRUE,
                 Colv = TRUE,
                 lwid = c(1.5, 4),
                 lhei = c(1.5, 4))
      }, error = function(e) {
        # Fallback to simple heatmap if clustering fails
        cat(sprintf("    Clustering failed for %s, using simple heatmap: %s\n", group_name, e$message))
        image(t(scale(expr_matrix)), 
              main = sprintf("Marker Expression - %s\n(%d clusters across files)", 
                            group_name, nrow(expr_matrix)),
              col = colors,
              axes = FALSE,
              xlab = "Markers", ylab = "Clusters")
        axis(1, at = seq(0, 1, length.out = ncol(expr_matrix)), 
             labels = colnames(expr_matrix), las = 2, cex.axis = max(0.5, min(1.0, 12/n_markers)))
        axis(2, at = seq(0, 1, length.out = nrow(expr_matrix)), 
             labels = rownames(expr_matrix), las = 2, cex.axis = max(0.5, min(1.0, 15/n_clusters)))
      })
      
      dev.off()
      
      cat(sprintf("  %s: %d clusters -> %s\n", 
                  group_name, nrow(expr_matrix), heatmap_filename))
      
    }, error = function(e) {
      warning(sprintf("Error creating heatmap for %s: %s", group_name, e$message))
    })
  }
  
  cat(sprintf("Expression heatmaps saved in: %s\n", heatmaps_dir))
  
  return(invisible(NULL))
}

# Define timestamp and test
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

# Test heatmap creation
cat("Testing expression heatmap creation...\n")
create_group_expression_heatmaps(timestamp)
cat("=== Done ===\n")