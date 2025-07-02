# ========================= PARAMETERS =========================

fallback_dir <- "/home/ngr24/Desktop/Maerklin/Cytometry/R_scripts"
input_file_path <- "/data/Cytometry/Conventional FC/files_ICS AML_pre-gating_9_16_2025/"
results_path <- "/data/Cytometry/Conventional FC/files_ICS AML_pre-gating_9_16_2025/Results"
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
# Set to TRUE for formats you want, FALSE to skip
# You can choose any combination: none, one, two, or all three formats
OUTPUT_FORMATS <- list(
  pdf = TRUE,    # Vector format for publications (scalable, small file size)
  png = TRUE,    # Raster format for pwp (high quality, larger file size)
  svg = TRUE     # Vector format for web (scalable, editable)
)

# Plot dimensions (inches)
PLOT_WIDTH <- 8
PLOT_HEIGHT <- 6
PNG_DPI <- 300  # Resolution for PNG files

# ========================= SETUP =========================
# Set working directory to script location
setwd(fallback_dir)

source_files <- c("Packages.R", "Functions.R", "UMAP_functions.R") # UMAP.R

lapply(source_files, source, chdir = TRUE)

# Setup output directory
umap_dir <- paste0(results_path, "/umap_plots/")
dir_create(c(results_path, umap_dir))

# ========================= DIRECT FCS LOADING =========================

# Load FCS files into flow frames and extract sample names
flow_frames <- load_fcs_folder(input_file_path)

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
ffr <- flow_frames[1:2]
cat("Extracting expression data using selected markers...\n")

expression_data <- map(seq_along(ffr), ~ {
  # Extract expression matrix for selected markers only
  expr_matrix <- flowCore::exprs(flow_frames[[.x]])[,marker_indices_of_interest,
                                                    drop = FALSE]

  # Create sample labels for each cell
  n_cells <- nrow(expr_matrix)
  sample_name <- paste0("Sample_", .x)  # or extract from file names if available
  labels <- rep(sample_name, n_cells)
  
  cat("Sample", .x, ":", n_cells, "cells x", ncol(expr_matrix), "markers\n")
  
  # Return both matrix and labels
  list(matrix = expr_matrix, labels = labels)
})

# Combine data and clean
df_combined <- map(expression_data, "matrix") %>% do.call(rbind, .)
sample_labels <- map(expression_data, "labels") %>% unlist()

cat("Combined data:", nrow(df_combined), "cells x",
    ncol(df_combined), "markers\n")

# Remove NA values efficiently
complete_rows <- complete.cases(df_combined)
df_clean <- df_combined[complete_rows, ]
sample_labels_clean <- sample_labels[complete_rows]

cat("After removing NAs:", nrow(df_clean), "cells\n")

# Save the combined data
data.table::fwrite(df_clean, paste0(umap_dir, "/combined_expression.csv"))

# ========================= UMAP COMPUTATION =========================
cat("Computing UMAP embedding...\n")
umap_result <- uwot::umap(X = df_clean,
                          n_neighbors = 15,
                          min_dist = 0.1,
                          verbose = TRUE)

# Create final dataframe with UMAP coordinates and labels
umap_df_final <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
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
cat("Creating marker-specific plots for each sample...\n")

# Get marker names and prepare data
selected_marker_names <- marker_info$markers[marker_indices_of_interest] %>%
  map_chr(~ if(.x == "" || is.na(.x)) paste0("Marker_", which(marker_info$markers == .x)) else .x)

# Add expression data to UMAP dataframe and calculate ranges
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

# Generate plots efficiently using base R
plot_combinations <- expand_grid(
  marker_idx = seq_along(marker_indices_of_interest),
  sample = unique_samples
)

plot_combinations$marker_name <- selected_marker_names[plot_combinations$marker_idx]
plot_combinations$clean_marker_name <- make.names(plot_combinations$marker_name)
plot_combinations$plot_num <- seq_len(nrow(plot_combinations))

total_plots <- nrow(plot_combinations)

pwalk(plot_combinations, function(marker_idx, sample, marker_name, clean_marker_name, plot_num) {
  
  cat("Creating plots for marker:", marker_name, "\n")
  
  # Filter data for this sample using base R
  sample_data <- umap_with_markers[umap_with_markers$Sample == sample, ]
  
  # Create plot
  marker_plot <- sample_data %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = .data[[clean_marker_name]])) +
    geom_point(size = 0.1, alpha = 0.7) +
    scale_color_viridis_c(name = marker_name, 
                         limits = marker_ranges[[clean_marker_name]], 
                         option = "viridis") +
    labs(title = paste("UMAP -", marker_name, "expression"),
         subtitle = paste("Sample:", sample),
         x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(legend.position = "right")

  # Create safe filename and save
  filename_base <- paste0(umap_dir, "/umap_",
                          gsub("[^A-Za-z0-9_-]", "_", marker_name), "_",
                          gsub("[^A-Za-z0-9_-]", "_", sample))

  marker_formats <- save_plot_formats(marker_plot, filename_base)

  cat("  Plot", plot_num, "/", total_plots, ":",
      gsub("[^A-Za-z0-9_-]", "_", marker_name), "in",
      gsub("[^A-Za-z0-9_-]", "_", sample),
      "(", paste(marker_formats, collapse=", "), ")\n")
})

# ========================= COMPLETION SUMMARY =========================
enabled_formats <- names(OUTPUT_FORMATS)[unlist(OUTPUT_FORMATS)]
total_plots <- 1 + (length(marker_indices_of_interest) * length(unique_samples))

cat("✓ UMAP analysis complete!\n")
cat("Results saved in:", umap_dir, "\n")
cat("Files created:\n  - combined_expression.csv\n  - umap_results.csv\n")
cat("  - Sample overview plot (", paste(enabled_formats, collapse = ", "), ")\n")
cat("  -", length(marker_indices_of_interest) * length(unique_samples), 
    "marker-specific plots (", paste(enabled_formats, collapse = ", "), ")\n")
cat("Total plots generated:", total_plots, "plots in", length(enabled_formats), "format(s)\n")
