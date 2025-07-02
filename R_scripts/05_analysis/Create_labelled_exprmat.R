#!/usr/bin/env Rscript

# =========================
# Cluster-by-Cluster Cell Extractor
# Extract ALL cells by cluster into separate files + create cluster tracking heatmap
# Output: Many files per cluster + comprehensive heatmap of file×cluster cell counts
# =========================

# =========================
# CONFIG — EDIT THIS BLOCK
# =========================
CONFIG <- list(
  # 1) Input folder with per-cell CSVs
  input_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR",
  
  # 2) Output folder for cluster-split files
  output_dir = "/home/ngr24/Desktop/Maerklin/Cytometry/Output/Split_by_cluster/",
  
  # 3) Interactive column selection
  markers_of_interest = NULL,      # Will be selected via GUI
  umap_x = NULL,                   # Will be selected via GUI  
  umap_y = NULL,                   # Will be selected via GUI
  
  # 4) Cluster filtering
  cluster_regex = "Clustering annotated",   # Regex to identify cluster columns
  min_cells_per_cluster = 5,                # Skip clusters with fewer cells
  
  # 5) Cluster name cleaning
  clean_cluster_names = TRUE,      # Remove AI_Cluster_ prefix and (Clustering annotated) suffix
  create_group_heatmaps = TRUE,    # Create separate heatmaps per cluster group
  create_group_matrices = TRUE,    # Create expression matrices per cluster group
  
  # 6) Processing options
  overwrite_existing = FALSE,      # Skip files that already exist in output
  create_heatmap = TRUE,           # Create cluster×file heatmap
  chunk_size = 50000,              # For reading very large files only
  
  # 7) Random seed
  seed = 123
)

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(data.table)   # Fast CSV operations
  library(progress)     # Progress bars
  library(stringr)      # String manipulation
  library(tcltk)        # GUI components
  library(ggplot2)      # For heatmap
  library(reshape2)     # For data reshaping
  library(RColorBrewer) # Color palettes
  library(gplots)       # For heatmap.2
})

# Set seed for reproducibility
set.seed(CONFIG$seed)

# =========================
# GUI Selection Functions
# =========================

# Interactive column selection via GUI
select_columns_gui <- function(csv_file) {
  cat("Loading sample for column selection...\n")
  sample_data <- fread(csv_file, nrows = 10)
  all_columns <- colnames(sample_data)
  
  cat(sprintf("Total columns available: %d\n", length(all_columns)))
  
  # Create main window
  tt <- tktoplevel()
  tkwm.title(tt, "Column Selection")
  tkwm.geometry(tt, "800x600")
  
  # Variables to store selections
  selected_markers <- character(0)
  umap_x_var <- character(0)
  umap_y_var <- character(0)
  
  # Create notebook for tabs
  nb <- tkwidget(tt, "ttk::notebook")
  
  # Tab 1: Markers of interest
  tab1 <- tkframe(nb)
  tkadd(nb, tab1, text = "Markers of Interest")
  
  # Scrollable listbox for markers
  frame1 <- tkframe(tab1)
  tkpack(frame1, fill = "both", expand = TRUE)
  
  scrollbar1 <- tkscrollbar(frame1, repeatinterval = 5, command = function(...) tkyview(listbox1, ...))
  listbox1 <- tklistbox(frame1, height = 20, selectmode = "multiple", yscrollcommand = function(...) tkset(scrollbar1, ...))
  
  tkpack(scrollbar1, side = "right", fill = "y")
  tkpack(listbox1, side = "left", fill = "both", expand = TRUE)
  
  # Populate markers listbox (exclude obvious metadata columns)
  marker_columns <- all_columns[!grepl("^(File|Time|Event|Original_ID)", all_columns, ignore.case = TRUE)]
  for (col in marker_columns) {
    tkinsert(listbox1, "end", col)
  }
  
  # Tab 2: UMAP coordinates
  tab2 <- tkframe(nb)
  tkadd(nb, tab2, text = "UMAP Coordinates")
  
  # UMAP X selection
  tkpack(tklabel(tab2, text = "Select UMAP X coordinate:"), anchor = "w")
  umap_x_var <- tclVar("")
  umap_x_combo <- tkwidget(tab2, "ttk::combobox", textvariable = umap_x_var, values = all_columns)
  tkpack(umap_x_combo, fill = "x", padx = 5, pady = 5)
  
  # UMAP Y selection  
  tkpack(tklabel(tab2, text = "Select UMAP Y coordinate:"), anchor = "w")
  umap_y_var <- tclVar("")
  umap_y_combo <- tkwidget(tab2, "ttk::combobox", textvariable = umap_y_var, values = all_columns)
  tkpack(umap_y_combo, fill = "x", padx = 5, pady = 5)
  
  tkpack(nb, fill = "both", expand = TRUE)
  
  # Button frame
  button_frame <- tkframe(tt)
  tkpack(button_frame, fill = "x", side = "bottom")
  
  # OK button
  ok_button <- tkbutton(button_frame, text = "OK", command = function() {
    # Get selected markers
    sel_indices <- as.integer(tkcurselection(listbox1))
    selected_markers <<- marker_columns[sel_indices + 1]  # R is 1-indexed
    
    # Get UMAP coordinates
    umap_x_val <- tclvalue(umap_x_var)
    umap_y_val <- tclvalue(umap_y_var)
    
    if (umap_x_val != "") umap_x_var <<- umap_x_val
    if (umap_y_val != "") umap_y_var <<- umap_y_val
    
    tkdestroy(tt)
  })
  
  # Cancel button
  cancel_button <- tkbutton(button_frame, text = "Cancel", command = function() {
    selected_markers <<- character(0)
    tkdestroy(tt)
  })
  
  tkpack(ok_button, cancel_button, side = "left", padx = 5, pady = 5)
  
  # Wait for window to close
  tkwait.window(tt)
  
  return(list(
    markers = selected_markers,
    umap_x = if(length(umap_x_var) > 0) umap_x_var else NULL,
    umap_y = if(length(umap_y_var) > 0) umap_y_var else NULL
  ))
}

# =========================
# Utility Functions
# =========================

# Clean cluster names by removing prefixes and suffixes
clean_cluster_name <- function(cluster_name) {
  if (!CONFIG$clean_cluster_names) {
    return(cluster_name)
  }
  
  # Remove "AI_Cluster_" prefix
  cleaned <- gsub("^AI_Cluster_", "", cluster_name)
  
  # Remove " (Clustering annotated)" suffix
  cleaned <- gsub(" \\(Clustering annotated\\)$", "", cleaned)
  
  return(cleaned)
}

# Extract cluster group from cluster name using regex patterns
extract_cluster_group <- function(cluster_name) {
  cleaned_name <- clean_cluster_name(cluster_name)
  
  # Define patterns to extract groups
  # Pattern 1: GroupName_Number (e.g., CD3-_227 -> CD3-)
  if (grepl("^([A-Za-z0-9+-]+)_[0-9]+$", cleaned_name)) {
    group <- gsub("^([A-Za-z0-9+-]+)_[0-9]+$", "\\1", cleaned_name)
    return(group)
  }
  
  # Pattern 2: Simple names without numbers (e.g., LymphoidCells -> LymphoidCells)
  if (!grepl("_[0-9]+$", cleaned_name)) {
    return(cleaned_name)
  }
  
  # Default: return the cleaned name as its own group
  return(cleaned_name)
}

# Group clusters by their base pattern
group_clusters <- function(cluster_names) {
  cleaned_names <- sapply(cluster_names, clean_cluster_name)
  groups <- sapply(cluster_names, extract_cluster_group)
  
  # Create mapping
  cluster_mapping <- data.frame(
    original = cluster_names,
    cleaned = cleaned_names,
    group = groups,
    stringsAsFactors = FALSE
  )
  
  return(cluster_mapping)
}

# Validate input parameters
validate_config <- function() {
  if (!dir.exists(CONFIG$input_dir)) {
    stop("Input directory does not exist: ", CONFIG$input_dir)
  }
  
  if (!dir.exists(CONFIG$output_dir)) {
    cat("Creating output directory:", CONFIG$output_dir, "\n")
    dir.create(CONFIG$output_dir, recursive = TRUE)
  }
  
  invisible(TRUE)
}

# Get list of CSV files to process
get_csv_files <- function() {
  csv_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    stop("No CSV files found in input directory: ", CONFIG$input_dir)
  }
  
  cat(sprintf("Found %d CSV files to process\n", length(csv_files)))
  return(csv_files)
}

# Extract cluster columns from data
extract_cluster_columns <- function(data) {
  cluster_cols <- colnames(data)[grepl(CONFIG$cluster_regex, colnames(data), ignore.case = TRUE)]
  
  if (length(cluster_cols) == 0) {
    warning("No cluster columns found matching pattern: ", CONFIG$cluster_regex)
    return(character(0))
  }
  
  cat(sprintf("Found %d cluster columns: %s\n", 
              length(cluster_cols), 
              paste(head(cluster_cols, 3), collapse = ", ")))
  
  return(cluster_cols)
}

# Extract and save cells by cluster for a single file
process_file_by_cluster <- function(csv_file, selected_columns, tracking_matrix) {
  file_basename <- tools::file_path_sans_ext(basename(csv_file))
  cat(sprintf("\n=== Processing: %s ===\n", file_basename))
  
  # Load data (chunked if very large)
  file_size <- file.info(csv_file)$size / 1024^2  # MB
  
  if (file_size > 100) {  # If larger than 100MB, read in chunks
    cat(sprintf("Large file detected (%.1f MB). Reading in chunks...\n", file_size))
    
    # Get total rows first
    total_rows <- length(readLines(csv_file)) - 1  # Subtract header
    
    # Read in chunks
    data_chunks <- list()
    chunk_start <- 1
    chunk_num <- 1
    
    while (chunk_start <= total_rows) {
      chunk_end <- min(chunk_start + CONFIG$chunk_size - 1, total_rows)
      
      chunk_data <- fread(csv_file, 
                         skip = ifelse(chunk_start == 1, 0, chunk_start), 
                         nrows = chunk_end - chunk_start + 1,
                         header = ifelse(chunk_start == 1, TRUE, FALSE))
      
      if (chunk_start > 1) {
        # Add column names for subsequent chunks
        colnames(chunk_data) <- colnames(data_chunks[[1]])
      }
      
      data_chunks[[chunk_num]] <- chunk_data
      chunk_start <- chunk_end + 1
      chunk_num <- chunk_num + 1
    }
    
    # Combine chunks
    data <- rbindlist(data_chunks)
    rm(data_chunks)
    gc()
    
  } else {
    # Read entire file at once
    data <- fread(csv_file)
  }
  
  cat(sprintf("Loaded %d rows and %d columns\n", nrow(data), ncol(data)))
  
  # Extract cluster columns
  cluster_cols <- extract_cluster_columns(data)
  
  if (length(cluster_cols) == 0) {
    warning(sprintf("No cluster columns found in %s\n", file_basename))
    return(tracking_matrix)
  }
  
  # Subset data to selected columns + cluster columns
  keep_cols <- c(selected_columns$markers, selected_columns$umap_x, selected_columns$umap_y, cluster_cols)
  keep_cols <- keep_cols[!is.null(keep_cols)]  # Remove NULL values first
  keep_cols <- keep_cols[keep_cols %in% colnames(data)]  # Only keep columns that exist
  keep_cols <- as.character(keep_cols)  # Ensure character vector
  
  if (length(keep_cols) == 0) {
    warning(sprintf("No selected columns found in %s\n", file_basename))
    return(tracking_matrix)
  }
  
  cat(sprintf("Keeping %d columns: %s\n", length(keep_cols), paste(head(keep_cols, 5), collapse = ", ")))
  
  # Subset data using standard data.table column selection
  data_subset <- data[, keep_cols, with = FALSE]
  
  # Process each cluster
  cluster_progress <- progress_bar$new(
    format = paste0("  [:bar] :percent :current/:total clusters processed in ",
                   sprintf("%.1f min", as.numeric(difftime(Sys.time(), start_time, units = "mins")))),
    total = length(cluster_cols), 
    clear = FALSE, 
    width = 80
  )
  
  for (cluster_col in cluster_cols) {
    cluster_progress$tick()
    
    # Extract cells belonging to this cluster (where column value is 1)
    cluster_cells <- data_subset[get(cluster_col) == 1]
    
    # Skip if too few cells
    if (nrow(cluster_cells) < CONFIG$min_cells_per_cluster) {
      next
    }
    
    # Remove cluster columns from output (keep only markers and UMAP)
    output_cols <- setdiff(colnames(cluster_cells), cluster_cols)
    cluster_output <- cluster_cells[, output_cols, with = FALSE]
    
    # Create output filename using cleaned cluster name
    cleaned_cluster_name <- clean_cluster_name(cluster_col)
    output_filename <- sprintf("%s_%s.csv", file_basename, cleaned_cluster_name)
    output_path <- file.path(CONFIG$output_dir, output_filename)
    
    # Skip if file already exists (unless overwrite is enabled)
    if (file.exists(output_path) && !CONFIG$overwrite_existing) {
      cat(sprintf("    Skipping %s (already exists)\n", cleaned_cluster_name))
      next
    }
    
    # Save cluster data
    fwrite(cluster_output, output_path)
    cat(sprintf("    Saved %d cells for cluster: %s\n", nrow(cluster_output), cleaned_cluster_name))
    
    # Update tracking matrix (use cleaned name as column key)
    cleaned_col_name <- clean_cluster_name(cluster_col)
    if (file_basename %in% rownames(tracking_matrix)) {
      if (cleaned_col_name %in% colnames(tracking_matrix)) {
        tracking_matrix[file_basename, cleaned_col_name] <- nrow(cluster_output)
      } else {
        # Add new column if cluster not in matrix yet
        tracking_matrix <- cbind(tracking_matrix, 0)
        colnames(tracking_matrix)[ncol(tracking_matrix)] <- cleaned_col_name
        tracking_matrix[file_basename, cleaned_col_name] <- nrow(cluster_output)
      }
    } else {
      # Add new row if file not in matrix yet
      new_row <- rep(0, ncol(tracking_matrix))
      names(new_row) <- colnames(tracking_matrix)
      if (cleaned_col_name %in% names(new_row)) {
        new_row[cleaned_col_name] <- nrow(cluster_output)
      } else {
        # Expand matrix to include new cluster
        tracking_matrix <- cbind(tracking_matrix, 0)
        colnames(tracking_matrix)[ncol(tracking_matrix)] <- cleaned_col_name
        new_row <- c(new_row, nrow(cluster_output))
        names(new_row)[length(new_row)] <- cleaned_col_name
      }
      tracking_matrix <- rbind(tracking_matrix, new_row)
      rownames(tracking_matrix)[nrow(tracking_matrix)] <- file_basename
    }
  }
  
  return(tracking_matrix)
}

# Create heatmap of file×cluster combinations
create_cluster_heatmap <- function(tracking_matrix) {
  if (!CONFIG$create_heatmap) {
    return(invisible(NULL))
  }
  
  cat("\n=== Creating cluster tracking heatmap ===\n")
  
  # Convert matrix to long format for ggplot
  tracking_df <- melt(as.matrix(tracking_matrix), varnames = c("File", "Cluster"), value.name = "Cell_Count")
  
  # Remove combinations with zero cells
  tracking_df <- tracking_df[tracking_df$Cell_Count > 0, ]
  
  if (nrow(tracking_df) == 0) {
    warning("No data to plot in heatmap")
    return(invisible(NULL))
  }
  
  # Add cluster groups
  tracking_df$Group <- sapply(tracking_df$Cluster, extract_cluster_group)
  
  # Create main heatmap
  p <- ggplot(tracking_df, aes(x = Cluster, y = File, fill = Cell_Count)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradientn(
      colors = brewer.pal(9, "YlOrRd"),
      name = "Cell Count",
      trans = "log10",
      labels = function(x) format(x, scientific = FALSE, big.mark = ",")
    ) +
    labs(
      title = "Cell Count Distribution Across Files and Clusters",
      subtitle = sprintf("Total files: %d, Total clusters: %d, Groups: %d", 
                        length(unique(tracking_df$File)), 
                        length(unique(tracking_df$Cluster)),
                        length(unique(tracking_df$Group))),
      x = "Cluster",
      y = "File"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )
  
  # Save main heatmap
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  heatmap_filename <- sprintf("cluster_distribution_heatmap_%s.png", timestamp)
  heatmap_path <- file.path(CONFIG$output_dir, heatmap_filename)
  
  ggsave(heatmap_path, p, width = 12, height = 8, dpi = 300, bg = "white")
  cat(sprintf("Main heatmap saved: %s\n", heatmap_path))
  
  # Create group-specific heatmaps if requested
  if (CONFIG$create_group_heatmaps) {
    create_group_heatmaps(tracking_df, timestamp)
  }
  
  # Create group expression matrices if requested
  if (CONFIG$create_group_matrices) {
    create_group_expression_matrices(timestamp)
    create_group_expression_heatmaps(timestamp)
  }
  
  # Save tracking matrix as CSV
  matrix_filename <- sprintf("cluster_tracking_matrix_%s.csv", timestamp)
  matrix_path <- file.path(CONFIG$output_dir, matrix_filename)
  write.csv(tracking_matrix, matrix_path, row.names = TRUE)
  cat(sprintf("Tracking matrix saved: %s\n", matrix_path))
  
  return(invisible(p))
}

# Create separate heatmaps for each cluster group
create_group_heatmaps <- function(tracking_df, timestamp) {
  cat("Creating group-specific heatmaps...\n")
  
  # Get unique groups
  groups <- unique(tracking_df$Group)
  groups <- groups[!is.na(groups)]
  
  cat(sprintf("Found %d cluster groups: %s\n", 
              length(groups), 
              paste(head(groups, 5), collapse = ", ")))
  
  # Create subdirectory for group heatmaps
  group_dir <- file.path(CONFIG$output_dir, "group_heatmaps")
  if (!dir.exists(group_dir)) {
    dir.create(group_dir, recursive = TRUE)
  }
  
  for (group in groups) {
    # Filter data for this group
    group_data <- tracking_df[tracking_df$Group == group, ]
    
    if (nrow(group_data) == 0) {
      next
    }
    
    # Skip groups with only one cluster (not very informative)
    if (length(unique(group_data$Cluster)) < 2) {
      next
    }
    
    # Create group-specific heatmap
    p_group <- ggplot(group_data, aes(x = Cluster, y = File, fill = Cell_Count)) +
      geom_tile(color = "white", linewidth = 0.1) +
      scale_fill_gradientn(
        colors = brewer.pal(9, "YlOrRd"),
        name = "Cell Count",
        trans = "log10",
        labels = function(x) format(x, scientific = FALSE, big.mark = ",")
      ) +
      labs(
        title = sprintf("Cluster Group: %s", group),
        subtitle = sprintf("Files: %d, Clusters in group: %d", 
                          length(unique(group_data$File)), 
                          length(unique(group_data$Cluster))),
        x = "Cluster",
        y = "File"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
      )
    
    # Save group heatmap
    safe_group_name <- gsub("[^A-Za-z0-9_-]", "_", group)  # Make filename safe
    group_filename <- sprintf("group_%s_heatmap_%s.png", safe_group_name, timestamp)
    group_path <- file.path(group_dir, group_filename)
    
    # Adjust plot size based on number of clusters
    n_clusters <- length(unique(group_data$Cluster))
    plot_width <- max(8, min(16, n_clusters * 0.5 + 4))
    plot_height <- max(6, min(12, length(unique(group_data$File)) * 0.3 + 4))
    
    ggsave(group_path, p_group, width = plot_width, height = plot_height, dpi = 300, bg = "white")
    cat(sprintf("  Group '%s' heatmap saved: %s\n", group, basename(group_path)))
  }
  
  cat(sprintf("Group heatmaps saved in: %s\n", group_dir))
}

# Create expression matrices per cluster group
create_group_expression_matrices <- function(timestamp) {
  if (!CONFIG$create_group_matrices) {
    return(invisible(NULL))
  }
  
  cat("\n=== Creating expression matrices per cluster group ===\n")
  
  # Create subdirectory for group matrices
  matrices_dir <- file.path(CONFIG$output_dir, "group_expression_matrices")
  if (!dir.exists(matrices_dir)) {
    dir.create(matrices_dir, recursive = TRUE)
  }
  
  # Get all cluster output files
  output_files <- list.files(CONFIG$output_dir, pattern = "\\.csv$", full.names = TRUE)
  output_files <- output_files[!grepl("cluster_tracking_matrix", output_files)]  # Exclude tracking matrix
  
  if (length(output_files) == 0) {
    warning("No cluster output files found")
    return(invisible(NULL))
  }
  
  cat(sprintf("Processing %d cluster files...\n", length(output_files)))
  
  # Initialize list to collect data by group
  group_data_list <- list()
  
  # Process each cluster file
  for (file_path in output_files) {
    file_basename <- basename(file_path)
    
    # Extract cluster name from filename (format: filename_clustername.csv)
    parts <- strsplit(gsub("\\.csv$", "", file_basename), "_")[[1]]
    if (length(parts) < 2) next
    
    # Reconstruct cluster name (everything after first underscore)
    cluster_name <- paste(parts[2:length(parts)], collapse = "_")
    
    # Determine group
    group_name <- extract_cluster_group(cluster_name)
    
    # Skip if group is just the cluster itself (single-cluster groups)
    if (group_name == cluster_name) next
    
    # Read cluster data
    tryCatch({
      cluster_data <- fread(file_path)
      
      if (nrow(cluster_data) == 0) next
      
      # Add metadata columns
      cluster_data$file_origin <- gsub("_[^_]+\\.csv$", "", file_basename)
      cluster_data$cluster_name <- cluster_name
      cluster_data$cluster_group <- group_name
      
      # Add to group collection
      if (group_name %in% names(group_data_list)) {
        group_data_list[[group_name]] <- rbind(group_data_list[[group_name]], cluster_data)
      } else {
        group_data_list[[group_name]] <- cluster_data
      }
      
    }, error = function(e) {
      warning(sprintf("Error reading file %s: %s", file_basename, e$message))
    })
  }
  
  if (length(group_data_list) == 0) {
    warning("No valid group data found")
    return(invisible(NULL))
  }
  
  cat(sprintf("Found %d cluster groups with multiple clusters\n", length(group_data_list)))
  
  # Save combined matrices per group
  summary_stats <- data.frame(
    Group = character(),
    Total_Cells = numeric(),
    Unique_Clusters = numeric(),
    Files_Represented = numeric(),
    Matrix_File = character(),
    stringsAsFactors = FALSE
  )
  
  for (group_name in names(group_data_list)) {
    group_data <- group_data_list[[group_name]]
    
    # Skip groups with too few cells
    if (nrow(group_data) < 100) {
      cat(sprintf("  Skipping group '%s' (only %d cells)\n", group_name, nrow(group_data)))
      next
    }
    
    # Create safe filename
    safe_group_name <- gsub("[^A-Za-z0-9_-]", "_", group_name)
    matrix_filename <- sprintf("expression_matrix_%s_%s.csv", safe_group_name, timestamp)
    matrix_path <- file.path(matrices_dir, matrix_filename)
    
    # Sort data by cluster and file origin for better organization
    group_data <- group_data[order(group_data$cluster_name, group_data$file_origin), ]
    
    # Save matrix
    fwrite(group_data, matrix_path)
    
    # Collect summary statistics
    n_cells <- nrow(group_data)
    n_clusters <- length(unique(group_data$cluster_name))
    n_files <- length(unique(group_data$file_origin))
    
    summary_stats <- rbind(summary_stats, data.frame(
      Group = group_name,
      Total_Cells = n_cells,
      Unique_Clusters = n_clusters,
      Files_Represented = n_files,
      Matrix_File = matrix_filename,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  Group '%s': %s cells from %d clusters across %d files -> %s\n", 
                group_name, format(n_cells, big.mark = ","), n_clusters, n_files, matrix_filename))
  }
  
  # Save summary statistics
  summary_path <- file.path(matrices_dir, sprintf("group_matrices_summary_%s.csv", timestamp))
  fwrite(summary_stats, summary_path)
  
  cat(sprintf("Expression matrices saved in: %s\n", matrices_dir))
  cat(sprintf("Summary statistics saved: %s\n", basename(summary_path)))
  
  return(invisible(summary_stats))
}

# Create expression heatmaps per cluster group
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
  for (file_path in output_files) {
    file_basename <- basename(file_path)
    
    # Extract cluster name from filename
    parts <- strsplit(gsub("\\.csv$", "", file_basename), "_")[[1]]
    if (length(parts) < 2) next
    
    cluster_name <- paste(parts[2:length(parts)], collapse = "_")
    
    # Determine group
    group_name <- extract_cluster_group(cluster_name)
    
    # Clean cluster name for display
    clean_cluster_name <- clean_cluster_name(cluster_name)
    
    tryCatch({
      # Read cluster data
      cluster_data <- fread(file_path)
      
      if (nrow(cluster_data) == 0) next
      
      # Get marker columns (exclude UMAP coordinates if present)
      marker_cols <- setdiff(names(cluster_data), c("UMAP_1", "UMAP_2", "umap_x", "umap_y"))
      
      if (length(marker_cols) == 0) next
      
      # Calculate mean expression for this cluster
      cluster_means <- cluster_data[, lapply(.SD, mean, na.rm = TRUE), .SDcols = marker_cols]
      cluster_means$cluster_name <- clean_cluster_name
      cluster_means$original_cluster <- cluster_name
      cluster_means$file_origin <- gsub("_[^_]+\\.csv$", "", file_basename)
      
      # Create proper display name for heatmap (GroupName_Number format)
      display_name <- clean_cluster_name
      if (grepl("^([A-Za-z0-9+-]+)_[0-9]+$", clean_cluster_name)) {
        # Already in correct format
        display_name <- clean_cluster_name
      } else if (grepl("_[0-9]+$", clean_cluster_name)) {
        # Has number but might need group prefix adjustment
        display_name <- clean_cluster_name
      } else {
        # No number, might need special handling
        display_name <- paste0(group_name, "s_", gsub("[^0-9]", "", cluster_name))
        if (display_name == paste0(group_name, "s_")) {
          display_name <- clean_cluster_name
        }
      }
      cluster_means$display_name <- display_name
      
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
      marker_cols <- setdiff(names(group_data), c("cluster_name", "original_cluster", "file_origin", "display_name"))
      
      if (length(marker_cols) == 0) {
        cat(sprintf("  Skipping %s (no marker columns)\n", group_name))
        next
      }
      
      # Create expression matrix
      expr_matrix <- as.matrix(group_data[, marker_cols, with = FALSE])
      # Use display_name for proper cluster naming in heatmap
      rownames(expr_matrix) <- group_data$display_name
      
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

# =========================
# Main Execution
# =========================

main <- function() {
  cat("=== Cluster-by-Cluster Cell Extractor ===\n")
  cat("Extracting ALL cells by cluster into separate files\n\n")
  
  # Record start time
  start_time <<- Sys.time()
  
  # Validate configuration
  validate_config()
  
  # Get input files
  csv_files <- get_csv_files()
  
  # Interactive column selection using first file as sample
  cat("\n=== Column Selection ===\n")
  selected_columns <- select_columns_gui(csv_files[1])
  
  if (length(selected_columns$markers) == 0) {
    stop("No markers selected. Exiting.")
  }
  
  cat(sprintf("Selected %d markers of interest\n", length(selected_columns$markers)))
  if (!is.null(selected_columns$umap_x)) cat(sprintf("UMAP X: %s\n", selected_columns$umap_x))
  if (!is.null(selected_columns$umap_y)) cat(sprintf("UMAP Y: %s\n", selected_columns$umap_y))
  
  # Initialize tracking matrix
  # First, scan all files to get all possible cluster names
  cat("\n=== Scanning for cluster columns ===\n")
  all_clusters <- character(0)
  
  for (csv_file in csv_files) {
    sample_data <- fread(csv_file, nrows = 1)
    file_clusters <- extract_cluster_columns(sample_data)
    all_clusters <- unique(c(all_clusters, file_clusters))
  }
  
  cat(sprintf("Found %d unique clusters across all files\n", length(all_clusters)))
  
  # Create cluster mapping with cleaned names
  if (length(all_clusters) > 0) {
    cluster_mapping <- group_clusters(all_clusters)
    cleaned_cluster_names <- unique(cluster_mapping$cleaned)
    
    cat(sprintf("After cleaning: %d unique cluster names\n", length(cleaned_cluster_names)))
    cat(sprintf("Cluster groups identified: %s\n", 
                paste(unique(cluster_mapping$group), collapse = ", ")))
  } else {
    stop("No cluster columns found in any files")
  }
  
  # Initialize tracking matrix (files × cleaned clusters)
  tracking_matrix <- matrix(0, nrow = 0, ncol = length(cleaned_cluster_names))
  colnames(tracking_matrix) <- cleaned_cluster_names
  
  # Process each file
  cat(sprintf("\n=== Processing %d files ===\n", length(csv_files)))
  
  file_progress <- progress_bar$new(
    format = "[:bar] :percent :current/:total files processed (:elapsed elapsed, :eta remaining)",
    total = length(csv_files), 
    clear = FALSE, 
    width = 80
  )
  
  for (i in seq_along(csv_files)) {
    file_progress$tick()
    csv_file <- csv_files[i]
    
    tracking_matrix <- process_file_by_cluster(csv_file, selected_columns, tracking_matrix)
    
    # Garbage collection after each file to manage memory
    gc()
  }
  
  # Create summary heatmap
  create_cluster_heatmap(tracking_matrix)
  
  # Final summary
  end_time <- Sys.time()
  processing_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat(sprintf("\n=== Processing Complete ===\n"))
  cat(sprintf("Total processing time: %.1f minutes\n", processing_time))
  cat(sprintf("Output directory: %s\n", CONFIG$output_dir))
  
  # Count output files
  output_files <- list.files(CONFIG$output_dir, pattern = "\\.csv$")
  cat(sprintf("Total output files created: %d\n", length(output_files)))
  
  # Summary statistics
  total_combinations <- sum(tracking_matrix > 0)
  cat(sprintf("File×Cluster combinations with data: %d\n", total_combinations))
  cat(sprintf("Average cells per combination: %.0f\n", mean(tracking_matrix[tracking_matrix > 0])))
  
  return(invisible(tracking_matrix))
}

# Run main function
if (!interactive()) {
  main()
}