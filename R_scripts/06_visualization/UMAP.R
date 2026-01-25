# =============================================================================
# UMAP Visualization Script for Cytometry Data
# =============================================================================
# Description: Creates UMAP dimensionality reduction plots for FCS flow 
#              cytometry data with marker expression overlays
#
# Dependencies: Packages.R, Functions.R, 06_visualization/visualization_functions.R
#
# TODOs:
#   - [ ] Add support for reading pre-computed UMAP coordinates from CSV
#   - [ ] Implement t-SNE as alternative dimensionality reduction
#   - [ ] Add clustering overlay options (FlowSOM, Leiden, etc.)
#   - [ ] Support for batch processing multiple experiments
#   - [ ] Add interactive HTML plots using plotly
#   - [ ] Implement UMAP parameter optimization (grid search)
# =============================================================================

# ========================= PARAMETERS =========================

fallback_dir <- "/home/ngr24/Desktop/Maerklin/Cytometry/R_scripts"
input_file_path <- "/data/Cytometry/Conventional FC/raw_MTB081_SARC/Cytolution/"
results_path <- "/data/Cytometry/Conventional FC/raw_MTB081_SARC/Results"
project_name <- "ICSAML"

# Marker selection options:
# - "popup": Open GUI window for interactive marker selection
# - "terminal": Choose from all filtered markers in terminal
# - "all": Use all filtered markers
# - vector of marker names: Use specific markers (e.g. c("CD4", "CD8", "TNF"))
# - vector of column numbers: Use specific columns (e.g. c(5, 7, 9))
MARKER_SELECTION_MODE <- "popup"  # Change this to control marker selection

# Patterns to use to filter out markers
EXCLUDE_PATTERNS <- c("^$", "ALD", "Column", "index", "Time", "FSC", "SSC")

# Output format options: Choose which formats to generate
OUTPUT_FORMATS <- list(
  pdf = FALSE,    # Vector format for publications (scalable, small file size)
  png = TRUE,    # Raster format for presentations (high quality, larger file size)
  svg = FALSE     # Vector format for web (scalable, editable)
)

# Plot dimensions (inches)
PLOT_WIDTH <- 8
PLOT_HEIGHT <- 6
PNG_DPI <- 300  # Resolution for PNG files

# File naming options
COMPACT_NAMING <- TRUE  # Remove repetitive parts from file names
MAX_SAMPLES <- NULL     # NULL = process all samples, or set number to limit

# ========================= HELPER FUNCTIONS =========================

# Function to create compact sample names by removing repetitive parts
create_compact_names <- function(file_names) {
  # Simple fallback if there's any issue
  if (is.null(file_names) || length(file_names) == 0) {
    return(character(0))
  }
  
  if (length(file_names) <= 1) {
    return(file_names)
  }
  
  # Direct approach for your specific filename pattern
  tryCatch({
    cat("DEBUG: Processing file names for compacting...\n")
    cat("DEBUG: Input names:", paste(file_names, collapse=" | "), "\n")
    
    # Your pattern: "export_MTB081SARC_1 V1W5_216080_Live"
    # We want to extract just the part after the space: "V1W5_216080_Live"
    
    # Split each name on spaces
    space_split <- strsplit(file_names, "\\s+")
    cat("DEBUG: After space split, lengths:", sapply(space_split, length), "\n")
    
    # Check if all names have at least 2 parts (before and after space)
    all_have_space <- all(sapply(space_split, length) >= 2)
    
    if (all_have_space) {
      # Extract the part after the first space (everything from second part onwards)
      compact_names <- sapply(space_split, function(parts) {
        if (length(parts) >= 2) {
          # Join everything after the first part
          paste(parts[2:length(parts)], collapse = " ")
        } else {
          paste(parts, collapse = " ")
        }
      })
      
      cat("DEBUG: Extracted compact names:", paste(compact_names, collapse=" | "), "\n")
      
      # Check if we got unique names
      if (length(unique(compact_names)) == length(compact_names)) {
        cat("DEBUG: Success - all names are unique\n")
        return(compact_names)
      } else {
        cat("DEBUG: Names not unique, trying alternative approach\n")
      }
    }
    
    # Alternative approach: remove common prefix from full names
    # Find the longest common prefix
    if (length(file_names) > 1) {
      # Find common prefix by comparing characters
      min_length <- min(nchar(file_names))
      common_prefix_length <- 0
      
      for (i in 1:min_length) {
        chars_at_i <- substr(file_names, i, i)
        if (length(unique(chars_at_i)) == 1) {
          common_prefix_length <- i
        } else {
          break
        }
      }
      
      cat("DEBUG: Common prefix length:", common_prefix_length, "\n")
      
      if (common_prefix_length > 0) {
        # Remove common prefix and any leading/trailing spaces or underscores
        compact_names <- substr(file_names, common_prefix_length + 1, nchar(file_names))
        compact_names <- gsub("^[\\s_]+|[\\s_]+$", "", compact_names)  # trim spaces/underscores
        
        cat("DEBUG: After prefix removal:", paste(compact_names, collapse=" | "), "\n")
        
        # Check if we got meaningful unique names
        if (length(unique(compact_names)) == length(compact_names) && 
            all(nchar(compact_names) > 0)) {
          return(compact_names)
        }
      }
    }
    
    # Final fallback: use original names (don't create Sample_1, Sample_2...)
    cat("DEBUG: Using original names as fallback\n")
    return(file_names)
    
  }, error = function(e) {
    cat("Warning: Error in create_compact_names:", e$message, "\n")
    cat("Using original names\n")
    return(file_names)
  })
}

# Function to create safe filename
create_safe_filename <- function(base_name, marker_name, sample_name) {
  # Clean marker name
  clean_marker <- gsub("[^A-Za-z0-9._-]", "_", marker_name)
  clean_marker <- gsub("_{2,}", "_", clean_marker)  # Remove multiple underscores
  clean_marker <- gsub("^_|_$", "", clean_marker)  # Remove leading/trailing underscores
  
  # Clean sample name
  clean_sample <- gsub("[^A-Za-z0-9._-]", "_", sample_name)
  clean_sample <- gsub("_{2,}", "_", clean_sample)
  clean_sample <- gsub("^_|_$", "", clean_sample)
  
  # Combine
  filename <- paste(base_name, clean_marker, clean_sample, sep = "_")
  return(filename)
}

# ========================= SETUP =========================
# Set working directory to parent R_scripts folder for consistent sourcing
setwd(fallback_dir)

# Source dependencies from proper locations
source("Packages.R")
source("Functions.R")
source("06_visualization/visualization_functions.R")

# Setup output directory
umap_dir <- paste0(results_path, "/umap_plots/")
dir_create(c(results_path, umap_dir))

# ========================= DIRECT FCS LOADING =========================

# Get FCS file paths and extract original names from filenames
input_fcs_files <- list.files(path = input_file_path, pattern = "\\.fcs$", full.names = TRUE)

if (length(input_fcs_files) == 0) {
  stop("No FCS files found in ", input_file_path)
}

original_names <- input_fcs_files %>%
  basename() %>%
  tools::file_path_sans_ext()

cat("Found", length(input_fcs_files), "FCS files\n")
cat("Full paths:", paste(head(input_fcs_files, 3), collapse = ", "), if(length(input_fcs_files) > 3) "..." else "", "\n")
cat("Original filenames:", paste(original_names, collapse = ", "), "\n")

# Load FCS files into flow frames
flow_frames <- load_fcs_list(input_fcs_files)

# Ensure we have the same number of names as flow frames
if (length(original_names) != length(flow_frames)) {
  cat("Warning: Mismatch between filenames and loaded frames, using sequential names\n")
  original_names <- paste0("Sample_", seq_along(flow_frames))
}

# Apply sample limit if specified
if (!is.null(MAX_SAMPLES) && length(flow_frames) > MAX_SAMPLES) {
  cat("Limiting to", MAX_SAMPLES, "samples out of", length(flow_frames), "\n")
  flow_frames <- flow_frames[1:MAX_SAMPLES]
  original_names <- original_names[1:MAX_SAMPLES]
}

# Create compact sample names
if (COMPACT_NAMING && length(original_names) > 1) {
  compact_names <- create_compact_names(original_names)
  if (length(compact_names) == length(original_names)) {
    cat("Original names:", paste(original_names, collapse = ", "), "\n")
    cat("Compact names:", paste(compact_names, collapse = ", "), "\n")
  } else {
    cat("Warning: compact names length mismatch, using original names\n")
    compact_names <- original_names
  }
} else {
  compact_names <- original_names
}

cat("Processing", length(flow_frames), "samples:", paste(compact_names, collapse = ", "), "\n")

# ========================= MARKER SELECTION =========================
# Get marker information and filter markers
marker_info <- explore_data(flow_frames[[1]])

# Filter out empty/unnamed markers and uninteresting ones
all_markers <- marker_info$markers[marker_info$mcols]
all_indices <- marker_info$mcols

# Apply filtering using base R to avoid namespace conflicts
exclude_regex <- paste(EXCLUDE_PATTERNS, collapse = "|")
keep_mask <- !grepl(exclude_regex, all_markers, ignore.case = TRUE)
filtered_markers <- all_markers[keep_mask]
filtered_indices <- all_indices[keep_mask]

# Simplified marker selection
marker_indices_of_interest <- select_markers_by_mode(MARKER_SELECTION_MODE,
                                                     filtered_markers,
                                                     filtered_indices)

cat("Final markers to use for UMAP:", length(marker_indices_of_interest),
    "markers | Indices:", paste(marker_indices_of_interest, collapse = ", "),
    "\n")
# ========================= DATA SUBSAMPLING =========================
cat("Extracting expression data using selected markers...\n")

expression_data <- map(seq_along(flow_frames), ~ {
  # Extract full expression matrix first
  full_matrix <- flowCore::exprs(flow_frames[[.x]])
  
  # Debug: Check if indices are accidentally included
  if (.x == 1) {
    cat("DEBUG: Full matrix dimensions:", nrow(full_matrix), "x", ncol(full_matrix), "\n")
    cat("DEBUG: Selected marker indices:", paste(marker_indices_of_interest, collapse=", "), "\n")
    cat("DEBUG: Max column index:", max(marker_indices_of_interest), "vs matrix cols:", ncol(full_matrix), "\n")
  }
  
  # Extract expression matrix for selected markers only - ensure numeric conversion
  expr_matrix <- full_matrix[, marker_indices_of_interest, drop = FALSE]
  
  # Ensure it's purely numeric and remove any row names/indices
  expr_matrix <- as.matrix(expr_matrix)
  rownames(expr_matrix) <- NULL
  
  # Create sample labels for each cell using compact names
  n_cells <- nrow(expr_matrix)
  sample_name <- compact_names[.x]
  labels <- rep(sample_name, n_cells)
  
  cat("Sample", .x, "(", sample_name, "):", n_cells, "cells x", ncol(expr_matrix), "markers\n")
  
  # Return both matrix and labels
  list(matrix = expr_matrix, labels = labels, sample_name = sample_name)
})

# Combine data and clean
df_combined <- map(expression_data, "matrix") %>% do.call(rbind, .)
sample_labels <- map(expression_data, "labels") %>% unlist()

# Remove row names that might contain indices
rownames(df_combined) <- NULL

# Set proper marker names as column names
selected_marker_names <- marker_info$markers[marker_indices_of_interest] %>%
  map_chr(~ if(.x == "" || is.na(.x)) paste0("Marker_", which(marker_info$markers == .x)) else .x)

if (ncol(df_combined) == length(selected_marker_names)) {
  colnames(df_combined) <- selected_marker_names
  cat("Set column names to marker names:", paste(head(selected_marker_names, 5), collapse = ", "), "\n")
} else {
  cat("Warning: Column count mismatch - cannot set marker names as column names\n")
}

# Ensure we only have numeric columns (no index columns)
if (ncol(df_combined) != length(marker_indices_of_interest)) {
  cat("Warning: Combined data has", ncol(df_combined), "columns but expected", 
      length(marker_indices_of_interest), "markers\n")
  # If there are extra columns, assume they are indices and remove them
  if (ncol(df_combined) > length(marker_indices_of_interest)) {
    cat("Removing", ncol(df_combined) - length(marker_indices_of_interest), "extra columns (likely indices)\n")
    df_combined <- df_combined[, 1:length(marker_indices_of_interest), drop = FALSE]
  }
}

cat("Combined data:", nrow(df_combined), "cells x",
    ncol(df_combined), "markers from", length(flow_frames), "samples\n")

# Remove NA values efficiently
complete_rows <- complete.cases(df_combined)
df_pre_clean <- df_combined[complete_rows, ]
sample_labels_clean <- sample_labels[complete_rows]

# Ensure we have only numeric expression data for UMAP
# Convert to numeric matrix and verify all values are numeric
df_clean <- as.matrix(df_pre_clean)
if (!is.numeric(df_clean)) {
  cat("Warning: Converting non-numeric data to numeric for UMAP\n")
  df_clean <- apply(df_clean, 2, as.numeric)
}

# Remove any rows with NAs that might have been introduced during conversion
final_complete_rows <- complete.cases(df_clean)
if (sum(!final_complete_rows) > 0) {
  cat("Removing", sum(!final_complete_rows), "rows with non-numeric values\n")
  df_clean <- df_clean[final_complete_rows, ]
  sample_labels_clean <- sample_labels_clean[final_complete_rows]
}

# Verify dimensions and data types
cat("After cleaning:", nrow(df_clean), "cells x", ncol(df_clean), "markers\n")
cat("Data type:", class(df_clean), "| All numeric:", all(sapply(df_clean, is.numeric)), "\n")
cat("Column names:", if(is.null(colnames(df_clean))) "None" else paste(head(colnames(df_clean), 5), collapse=", "), "\n")

# Debug: Check first few values of first column to see if it contains indices
if (nrow(df_clean) > 0) {
  first_col_sample <- head(df_clean[, 1], 10)
  cat("First column sample values:", paste(first_col_sample, collapse=", "), "\n")
  cat("First column range:", min(df_clean[, 1], na.rm=TRUE), "to", max(df_clean[, 1], na.rm=TRUE), "\n")
}

# ========================= DATA NORMALIZATION =========================
cat("Applying quantile normalization and scaling...\n")

# Store original data for comparison
df_raw <- df_clean

# Apply quantile normalization
# This makes the distribution of values identical across all markers
cat("Before normalization - marker ranges:\n")
for (i in seq_len(min(5, ncol(df_clean)))) {
  marker_range <- range(df_clean[, i], na.rm = TRUE)
  cat("  ", colnames(df_clean)[i], ": ", marker_range[1], " to ", marker_range[2], "\n")
}

# Quantile normalization function
quantile_normalize <- function(data_matrix) {
  # Calculate quantiles for each column
  ranks <- apply(data_matrix, 2, rank, ties.method = "average")
  
  # Calculate the empirical distribution (average quantile values across columns)
  sorted_means <- rowMeans(apply(data_matrix, 2, sort, na.last = TRUE), na.rm = TRUE)
  
  # Apply the common distribution to all columns
  normalized_matrix <- apply(ranks, 2, function(r) {
    # Handle NAs
    if (any(is.na(r))) {
      result <- rep(NA, length(r))
      result[!is.na(r)] <- sorted_means[r[!is.na(r)]]
      return(result)
    } else {
      return(sorted_means[r])
    }
  })

  # Preserve original matrix properties
  dimnames(normalized_matrix) <- dimnames(data_matrix)
  
  return(normalized_matrix)
}

# Apply quantile normalization
df_quantile_norm <- quantile_normalize(df_clean)

# Apply scaling (z-score normalization) after quantile normalization
df_scaled <- scale(df_quantile_norm, center = TRUE, scale = TRUE)

# Handle any NAs that might have been introduced
if (any(is.na(df_scaled))) {
  na_cols <- apply(df_scaled, 2, function(x) any(is.na(x)))
  if (any(na_cols)) {
    cat("Warning: NAs introduced during scaling for markers:", paste(colnames(df_scaled)[na_cols], collapse=", "), "\n")
    cat("Replacing NAs with 0 (mean after centering)\n")
    df_scaled[is.na(df_scaled)] <- 0
  }
}

# Update df_clean with normalized and scaled data
df_clean <- as.matrix(df_scaled)

cat("After quantile normalization and scaling - marker ranges:\n")
for (i in seq_len(min(5, ncol(df_clean)))) {
  marker_range <- range(df_clean[, i], na.rm = TRUE)
  cat("  ", colnames(df_clean)[i], ": ", sprintf("%.3f", marker_range[1]), " to ", sprintf("%.3f", marker_range[2]), "\n")
}

# Save the combined data (both raw and normalized versions)
data.table::fwrite(as.data.frame(df_raw), paste0(umap_dir, "/combined_expression_raw.csv"))
data.table::fwrite(as.data.frame(df_clean), paste0(umap_dir, "/combined_expression_normalized.csv"))

# ========================= UMAP COMPUTATION =========================
cat("Computing UMAP embedding on normalized data...\n")
umap_result <- uwot::umap(X = df_clean,
                          n_neighbors = 15,# 30
                          min_dist = 0.1,  # 0.5
                          verbose = TRUE)

umap_result2 <- uwot::umap(X = df_clean, n_neighbors = 30, min_dist = 0.5)

# Create final dataframe with UMAP coordinates and labels
umap_df_final <- data.frame(
  UMAP1 = umap_result2[, 1],
  UMAP2 = umap_result2[, 2],
  Sample = sample_labels_clean
)

# Save UMAP results
data.table::fwrite(umap_df_final, paste0(umap_dir, "/umap_results.csv"))

# Create and save sample overview plot
unique_samples <- unique(sample_labels_clean)
colors <- rainbow(length(unique_samples)) %>% setNames(unique_samples)

umap_plot <- umap_df_final %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Sample)) +
  geom_point(size = 0.1, alpha = 0.6) +
  scale_color_manual(values = colors) +
  labs(title = "UMAP - Colored by Sample", x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 2)))

cat("Saving sample overview plot...\n")
sample_formats <- save_plot_formats(umap_plot, paste0(umap_dir, "/umap_by_sample"))
cat("Sample plot saved in formats:", paste(sample_formats, collapse = ", "), "\n")

# ========================= MARKER-SPECIFIC PLOTS =========================
cat("Creating individual plots for each marker and sample combination...\n")

# Add expression data to UMAP dataframe and calculate global ranges for consistency
umap_with_markers <- umap_df_final
marker_ranges <- list()

iwalk(marker_indices_of_interest, ~ {
  marker_name <- selected_marker_names[.y]
  clean_marker_name <- make.names(marker_name)
  
  umap_with_markers[[clean_marker_name]] <<- df_clean[, .y]
  marker_ranges[[clean_marker_name]] <<- range(df_clean[, .y], na.rm = TRUE)
  
  cat("Marker", clean_marker_name, "range:",
      marker_ranges[[clean_marker_name]][1], "to",
      marker_ranges[[clean_marker_name]][2], "\n")
})

# Get unique samples
unique_samples <- unique(sample_labels_clean)
n_samples <- length(unique_samples)
n_markers <- length(marker_indices_of_interest)

cat("Will create", n_markers, "markers ×", n_samples, "samples =", n_markers * n_samples, "individual plots\n")

# Create individual plots for each marker and sample combination
plot_counter <- 0
total_plots <- n_markers * n_samples

# Loop through each marker
for (marker_idx in seq_along(marker_indices_of_interest)) {
  marker_name <- selected_marker_names[marker_idx]
  clean_marker_name <- make.names(marker_name)
  
  cat("Processing marker:", marker_name, "(", marker_idx, "/", n_markers, ")\n")
  
  # Loop through each sample
  for (sample_name in unique_samples) {
    plot_counter <- plot_counter + 1
    
    # Filter data for this sample
    sample_data <- umap_with_markers[umap_with_markers$Sample == sample_name, ]
    
    if (nrow(sample_data) == 0) {
      cat("  Warning: No data for sample", sample_name, "- skipping\n")
      next
    }
    
    # Create plot
    marker_plot <- sample_data %>%
      ggplot(aes(x = UMAP1, y = UMAP2, color = .data[[clean_marker_name]])) +
      geom_point(size = 0.1, alpha = 0.7) +
      scale_color_viridis_c(name = marker_name, 
                           limits = marker_ranges[[clean_marker_name]], 
                           option = "viridis") +
      labs(title = paste(marker_name, "expression"),
           subtitle = paste("Sample:", sample_name, 
                           paste0("(", nrow(sample_data), " cells)")),
           x = "UMAP1", y = "UMAP2") +
      theme_minimal() +
      theme(legend.position = "right",
            plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 10))

    # Create safe filename
    filename_base <- create_safe_filename("umap", marker_name, sample_name)
    full_path <- paste0(umap_dir, "/", filename_base)

    # Save plot in specified formats
    marker_formats <- save_plot_formats(marker_plot, full_path)

    cat("  Plot", plot_counter, "/", total_plots, ":", filename_base, 
        "(", paste(marker_formats, collapse = ", "), ")\n")
  }
}

# ========================= COMPLETION SUMMARY =========================

# Get enabled output formats
enabled_formats <- names(OUTPUT_FORMATS)[unlist(OUTPUT_FORMATS)]
sample_overview_plots <- 1
individual_plots <- n_markers * n_samples
total_plots <- sample_overview_plots + individual_plots

cat("✓ UMAP analysis complete!\n")
cat("Results saved in:", umap_dir, "\n")
cat("Processed", n_samples, "samples with compact naming:", paste(compact_names, collapse = ", "), "\n")
cat("Files created:\n")
cat("  - combined_expression.csv (", nrow(df_clean), "cells ×", ncol(df_clean), "markers)\n")
cat("  - umap_results.csv (UMAP coordinates + sample labels)\n")
cat("  - Sample overview plot (", paste(enabled_formats, collapse = ", "), ")\n")
cat("  -", individual_plots, "individual marker-sample plots (", paste(enabled_formats, collapse = ", "), ")\n")
cat("Total plots generated:", total_plots, "plots in", length(enabled_formats), "format(s)\n")

# Summary by sample
cat("\nPlots per sample:\n")
for (i in seq_along(compact_names)) {
  cat("  ", compact_names[i], ":", n_markers, "marker plots\n")
}

# Summary by marker
cat("\nPlots per marker:\n")
for (i in seq_along(selected_marker_names)) {
  cat("  ", selected_marker_names[i], ":", n_samples, "sample plots\n")
}