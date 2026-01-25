# =============================================================================
# UMAP Visualization Script with Memory-Efficient Processing
# =============================================================================
# Description: Memory-optimized UMAP dimensionality reduction for large FCS 
#              datasets. Supports subsampling and chunked processing.
#
# Dependencies: Packages.R, Functions.R, visualization_functions.R,
#               02_preprocessing/Downsampling_functions.R, 01_loading/Memory_config.R
#
# TODOs:
#   - [ ] Add automatic memory estimation based on file sizes
#   - [ ] Implement GPU acceleration for UMAP (via rapids/cuml)
#   - [ ] Add support for incremental/streaming UMAP
#   - [ ] Implement checkpoint saving for long computations
#   - [ ] Add multi-threaded file loading
#   - [ ] Support for distributed processing across multiple machines
# =============================================================================

# ========================= PARAMETERS =========================

fallback_dir <- "/home/ngr24/Desktop/Maerklin/Cytometry/R_scripts"
input_file_path <- "/data/Maerklin/ICK/ICK_T1_PR/"
output_dir <- "./Results"
project_name <- "ICK_T1_PR"

# Marker selection options:
# - "popup": Open GUI window for interactive marker selection
# - "auto": Use popup in interactive mode, all filtered markers in non-interactive mode
# - "all": Use all filtered markers (no popup)
# - vector of marker names: Use specific markers (e.g. c("CD4", "CD8", "TNF"))
# - vector of column numbers: Use specific columns (e.g. c(5, 7, 9))
MARKER_SELECTION <- "popup"  # Change this to control marker selection

# Patterns to use to filter out markers
EXCLUDE_PATTERNS <- c("^$", "ALD", "Column", "index", "Time", "FSC", "SSC")

# Memory-efficient processing settings for large FCS files

# Maximum cells to keep per sample (for memory management)
# Set to NULL to keep all cells, or a number to downsample
MAX_CELLS_PER_SAMPLE <- 50000  # Adjust based on available memory

# Maximum total cells for UMAP computation
# UMAP becomes very slow/memory intensive with >100k cells
MAX_TOTAL_CELLS_UMAP <- 100000

# Chunk size for processing large files (cells per chunk)
CHUNK_SIZE <- 10000

# Memory-safe plotting: limit number of points in plots
MAX_PLOT_POINTS <- 20000  # For visualization, more points don't improve quality

# Enable/disable memory monitoring
MONITOR_MEMORY <- TRUE

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
if (interactive()) {
  # Running in RStudio
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  # Running from command line
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) > 0) {
    script_dir <- dirname(script_path)
  } else {
    script_dir <- fallback_dir
  }
}
setwd(script_dir)

# Source dependencies from proper locations
source("../Packages.R")
source("../Functions.R")
source("visualization_functions.R")
source("../02_preprocessing/Downsampling_functions.R")
source("../01_loading/Memory_config.R")

# ========================= MEMORY MANAGEMENT FUNCTIONS =========================

# Function to check available system memory
check_system_memory <- function() {
  if (MONITOR_MEMORY) {
    tryCatch({
      if (Sys.info()["sysname"] == "Linux") {
        # Linux memory check
        mem_info <- system("cat /proc/meminfo", intern = TRUE)
        available_line <- grep("MemAvailable", mem_info, value = TRUE)
        if (length(available_line) > 0) {
          available_kb <- as.numeric(gsub(".*MemAvailable:\\s*(\\d+).*", "\\1", available_line))
          available_mb <- available_kb / 1024
          cat("Available system memory:", round(available_mb, 1), "MB\n")
          return(available_mb)
        }
      }
      # Fallback: estimate from R's memory usage
      mem_info <- gc()
      used_mb <- sum(mem_info[, 2]) * 8 / 1024
      cat("Estimated memory usage by R:", round(used_mb, 1), "MB\n")
      return(NULL)
    }, error = function(e) {
      cat("Could not determine system memory\n")
      return(NULL)
    })
  }
}

# Function to monitor memory usage
print_memory_usage <- function(label = "") {
  if (MONITOR_MEMORY) {
    mem_info <- gc()
    used_mb <- sum(mem_info[, 2]) * 8 / 1024  # Convert to MB (rough estimate)
    cat(sprintf("Memory %s: %.1f MB used\n", label, used_mb))
  }
}

# Function to safely load and downsample a single FCS file
load_fcs_safe <- function(file_path, max_cells = NULL, markers_only = NULL) {
  print_memory_usage("before loading FCS")
  
  # Load FCS file
  ff <- flowWorkspace::load_cytoframe_from_fcs(file_path, 
                                                transformation = FALSE, 
                                                truncate_max_range = FALSE)
  
  n_cells_original <- nrow(ff)
  n_markers <- ncol(ff)
  
  # Extract expression matrix
  exprs_data <- flowCore::exprs(ff)
  
  # Select only needed markers if specified
  if (!is.null(markers_only)) {
    valid_markers <- markers_only[markers_only <= ncol(exprs_data)]
    if (length(valid_markers) > 0) {
      exprs_data <- exprs_data[, valid_markers, drop = FALSE]
    }
  }
  
  # Downsample if needed
  if (!is.null(max_cells) && n_cells_original > max_cells) {
    sample_indices <- sample(n_cells_original, max_cells)
    exprs_data <- exprs_data[sample_indices, , drop = FALSE]
    cat("  Downsampled from", n_cells_original, "to", max_cells, "cells\n")
  }
  
  # Remove rows with any NA values
  complete_rows <- complete.cases(exprs_data)
  exprs_data <- exprs_data[complete_rows, , drop = FALSE]
  
  print_memory_usage("after processing FCS")
  
  list(
    data = exprs_data,
    n_cells_final = nrow(exprs_data),
    n_cells_original = n_cells_original,
    n_markers = ncol(exprs_data)
  )
}

# Function to create balanced sample from multiple files
create_balanced_sample <- function(file_paths, sample_names, markers, total_target_cells) {
  n_files <- length(file_paths)
  cells_per_file <- ceiling(total_target_cells / n_files)
  
  cat("Creating balanced sample:", cells_per_file, "cells per file for", n_files, "files\n")
  
  combined_data <- list()
  combined_labels <- c()
  
  for (i in seq_along(file_paths)) {
    cat("Processing file", i, "/", n_files, ":", basename(file_paths[i]), "\n")
    
    file_data <- load_fcs_safe(file_paths[i], 
                               max_cells = cells_per_file, 
                               markers_only = markers)
    
    combined_data[[i]] <- file_data$data
    combined_labels <- c(combined_labels, rep(sample_names[i], file_data$n_cells_final))
    
    cat("  Loaded:", file_data$n_cells_final, "cells x", file_data$n_markers, "markers\n")
    
    # Force garbage collection after each file
    gc()
  }
  
  # Combine all data
  final_data <- do.call(rbind, combined_data)
  
  print_memory_usage("after combining all files")
  
  list(
    data = final_data,
    labels = combined_labels,
    total_cells = nrow(final_data)
  )
}

# Helper function to save plots in selected formats
save_plot_formats <- function(plot, filename_base,
                width = PLOT_WIDTH, height = PLOT_HEIGHT) {
  formats_saved <- c()

  format_configs <- list(
    pdf = list(ext = ".pdf", device = "pdf", extra_args = list()),
    png = list(ext = ".png", device = "png", extra_args = list(dpi = PNG_DPI)),
    svg = list(ext = ".svg", device = "svg", extra_args = list())
  )

  for (format_name in names(OUTPUT_FORMATS)) {
    if (OUTPUT_FORMATS[[format_name]]) {
      config <- format_configs[[format_name]]
      args <- c(
        list(
          filename = paste0(filename_base, config$ext),
          plot = plot,
          width = width,
          height = height,
          device = config$device
        ),
        config$extra_args
      )
      do.call(ggplot2::ggsave, args)
      formats_saved <- c(formats_saved, toupper(format_name))
    }
  }
  return(formats_saved)
}

# ========================= MEMORY-SAFE DATA LOADING =========================

# Get file list and sample names first
input_fcs_files <- list.files(path = input_file_path, pattern = "\\.fcs$", full.names = TRUE)
sample_names <- input_fcs_files %>%
  basename() %>%
  gsub("\\.fcs$", "", .) %>%
  sapply(function(x) strsplit(x, "_")[[1]][1:2] %>% paste(collapse = "_"))

cat("Found", length(input_fcs_files), "FCS files in", input_file_path, "\n")
cat("Sample names:", paste(sample_names, collapse = ", "), "\n")

print_memory_usage("initial")
check_system_memory()

# Estimate dataset size and provide recommendations
file_sizes_mb <- file.size(input_fcs_files) / (1024^2)
total_size_mb <- sum(file_sizes_mb)
cat("Total FCS file size on disk:", round(total_size_mb, 1), "MB\n")

if (total_size_mb > MAX_AUTO_FILE_SIZE_MB) {
  cat("⚠️  LARGE DATASET WARNING: Total file size >", MAX_AUTO_FILE_SIZE_MB, "MB\n")
  cat("Consider using smaller MAX_CELLS_PER_SAMPLE and MAX_TOTAL_CELLS_UMAP values\n")
  cat("Current settings: MAX_CELLS_PER_SAMPLE =", MAX_CELLS_PER_SAMPLE %||% "unlimited", 
      ", MAX_TOTAL_CELLS_UMAP =", MAX_TOTAL_CELLS_UMAP %||% "unlimited", "\n")
}

# Load first file to get marker information (memory-efficient)
cat("Loading first file to analyze markers...\n")
first_file_data <- load_fcs_safe(input_fcs_files[1], max_cells = 1000)  # Just sample for markers
temp_ff <- flowWorkspace::load_cytoframe_from_fcs(input_fcs_files[1], 
                                                  transformation = FALSE, 
                                                  truncate_max_range = FALSE)

print_memory_usage("after loading first file")

# ========================= DATA PROCESSING =========================
# Get marker information and filter markers  
marker_info <- explore_data(temp_ff)
cat("Available markers:", length(marker_info$markers), "| Marker columns:", paste(marker_info$mcols, collapse = ", "), "\n")

# Clean up temporary objects
rm(temp_ff, first_file_data)
gc()

print_memory_usage("after marker analysis")

# Filter out empty/unnamed markers and uninteresting ones
exclude_regex <- paste(EXCLUDE_PATTERNS, collapse="|")
interesting_markers <- marker_info$markers[marker_info$mcols]
interesting_indices <- marker_info$mcols

# Apply filtering using base R to avoid namespace conflicts
keep_mask <- !grepl(exclude_regex, interesting_markers, ignore.case = TRUE)
filtered_markers <- interesting_markers[keep_mask]
filtered_indices <- interesting_indices[keep_mask]

cat("Filtered markers available for selection:\n")
walk2(filtered_indices, filtered_markers, ~ cat(.x, ":", .y, "\n"))

# Interactive marker selection with popup window
marker_indices_of_interest <- {
  if (is.character(MARKER_SELECTION) && length(MARKER_SELECTION) == 1) {
    if (MARKER_SELECTION == "all") {
      cat("Using all filtered markers\n")
      filtered_indices
    } else if (MARKER_SELECTION %in% c("popup", "auto")) {
      if (length(filtered_markers) == 0) {
        cat("No suitable markers found! Using all marker columns > 2\n")
        marker_info$mcols[marker_info$mcols > 2]
      } else {
        # Display available markers
        cat("\nFiltered markers available for UMAP analysis:\n")
        iwalk(filtered_markers, ~ cat(sprintf("%d. Column %d: %s\n", .y, filtered_indices[.y], .x)))
        
        # Try popup selection for "popup" mode or interactive "auto" mode
        if (MARKER_SELECTION == "popup" || (MARKER_SELECTION == "auto" && interactive())) {
          tryCatch({
            cat("\nOpening enhanced marker selection window...\n")
            selected_indices <- select_markers(list(indices = filtered_indices, names = filtered_markers))
            
            if (!is.null(selected_indices) && length(selected_indices) > 0) {
              cat("Selected markers from popup:\n")
              walk(selected_indices, ~ cat("  Column", .x, ":", marker_info$markers[.x], "\n"))
              selected_indices
            } else {
              cat("No markers selected or popup cancelled. Using all filtered markers.\n")
              filtered_indices
            }
          }, error = function(e) {
            cat("Popup selection failed:", e$message, "\nUsing all filtered markers.\n")
            filtered_indices
          })
        } else {
          cat("Non-interactive mode. Using all filtered markers.\n")
          filtered_indices
        }
      }
    } else {
      # Single marker name
      if (MARKER_SELECTION %in% marker_info$markers) {
        marker_idx <- which(marker_info$markers == MARKER_SELECTION)
        cat("Using specified marker:", MARKER_SELECTION, "(column", marker_idx, ")\n")
        marker_idx
      } else {
        cat("Warning: Marker", MARKER_SELECTION, "not found. Using all filtered markers.\n")
        filtered_indices
      }
    }
  } else if (is.character(MARKER_SELECTION) && length(MARKER_SELECTION) > 1) {
    # Vector of marker names
    cat("Using specified markers:", paste(MARKER_SELECTION, collapse = ", "), "\n")
    found_indices <- MARKER_SELECTION %>%
      map_int(~ if(.x %in% marker_info$markers) which(marker_info$markers == .x) else NA) %>%
      .[!is.na(.)]
    
    if (length(found_indices) == 0) {
      cat("No valid markers found. Using all filtered markers.\n")
      filtered_indices
    } else {
      found_indices
    }
  } else if (is.numeric(MARKER_SELECTION)) {
    # Vector of column numbers
    cat("Using specified columns:", paste(MARKER_SELECTION, collapse = ", "), "\n")
    valid_cols <- MARKER_SELECTION[MARKER_SELECTION %in% marker_info$mcols]
    if (length(valid_cols) > 0) {
      valid_cols
    } else {
      cat("Warning: No valid column numbers found. Using all filtered markers.\n")
      filtered_indices
    }
  } else {
    cat("Warning: Invalid MARKER_SELECTION parameter. Using all filtered markers.\n")
    filtered_indices
  }
}

# Final safety check and validation
marker_indices_of_interest <- marker_indices_of_interest %>%
  {if(length(.) == 0) {
    cat("Warning: No markers selected. Using all available marker columns.\n")
    marker_info$mcols
  } else .} %>%
  .[. <= ncol(flowCore::exprs(flowWorkspace::load_cytoframe_from_fcs(input_fcs_files[1])))]  # Check against first file

cat("Final markers to use for UMAP:", length(marker_indices_of_interest), "markers | Indices:", paste(marker_indices_of_interest, collapse=", "), "\n")

# Setup output directory
results_path <- "./Results"
umap_dir <- paste0(results_path, "/umap_plots/")
dir_create(c(results_path, umap_dir))

# ========================= MEMORY-EFFICIENT DATA LOADING AND CONCATENATION =========================
print_memory_usage("before loading all data")

# Calculate optimal sampling strategy
cat("Planning memory-efficient data loading...\n")

# Quick check of file sizes to estimate memory needs
file_info <- input_fcs_files %>%
  map_dfr(~ {
    ff_temp <- flowWorkspace::load_cytoframe_from_fcs(.x, transformation = FALSE, truncate_max_range = FALSE)
    n_cells <- nrow(ff_temp)
    n_markers <- ncol(ff_temp)
    rm(ff_temp)
    gc()
    data.frame(file = basename(.x), cells = n_cells, markers = n_markers)
  })

total_cells_available <- sum(file_info$cells)
cat("Total cells available across all files:", total_cells_available, "\n")

# Determine sampling strategy
if (!is.null(MAX_TOTAL_CELLS_UMAP) && total_cells_available > MAX_TOTAL_CELLS_UMAP) {
  cat("Large dataset detected. Using balanced sampling to", MAX_TOTAL_CELLS_UMAP, "cells\n")
  use_cells_per_file <- min(MAX_CELLS_PER_SAMPLE %||% Inf, 
                           floor(MAX_TOTAL_CELLS_UMAP / length(input_fcs_files)))
} else if (!is.null(MAX_CELLS_PER_SAMPLE)) {
  cat("Using per-file limit of", MAX_CELLS_PER_SAMPLE, "cells\n")
  use_cells_per_file <- MAX_CELLS_PER_SAMPLE
} else {
  cat("Loading all available cells\n")
  use_cells_per_file <- NULL
}

# Load data with memory-efficient approach
cat("Loading expression data with memory management...\n")
combined_result <- create_balanced_sample(
  input_fcs_files, 
  sample_names, 
  marker_indices_of_interest,
  MAX_TOTAL_CELLS_UMAP %||% total_cells_available
)

df_clean <- combined_result$data
sample_labels_clean <- combined_result$labels

cat("Final dataset:", nrow(df_clean), "cells x", ncol(df_clean), "markers\n")
print_memory_usage("after loading and combining data")

# Save the combined data
data.table::fwrite(df_clean, paste0(umap_dir, "/combined_expression.csv"))
print_memory_usage("after saving combined data")

# ========================= UMAP COMPUTATION =========================
cat("Computing UMAP embedding...\n")
print_memory_usage("before UMAP")

# Determine optimal number of threads
n_threads <- if (is.null(N_CORES)) {
  min(8, parallel::detectCores())
} else {
  min(N_CORES, parallel::detectCores())
}

cat("Using", n_threads, "CPU cores for UMAP computation\n")

# For very large datasets, consider using approximate UMAP
umap_result <- uwot::umap(
  X = df_clean, 
  n_neighbors = 15, 
  min_dist = 0.1, 
  verbose = TRUE,
  n_threads = n_threads
)

print_memory_usage("after UMAP computation")

# Create final dataframe with UMAP coordinates and labels
umap_df_final <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2],
  Sample = sample_labels_clean
)

# Clean up large objects we no longer need
rm(df_clean, umap_result)
gc()
print_memory_usage("after cleaning up UMAP objects")

# Save UMAP results
data.table::fwrite(umap_df_final, paste0(umap_dir, "/umap_results.csv"))

# Create and save sample overview plot
unique_samples <- unique(sample_labels_clean)
colors <- rainbow(length(unique_samples)) %>% setNames(unique_samples)

# Downsample for plotting if needed
plot_data_sample <- umap_df_final
if (nrow(umap_df_final) > MAX_PLOT_POINTS) {
  cat("Downsampling for visualization:", MAX_PLOT_POINTS, "points\n")
  sample_indices <- sample(nrow(umap_df_final), MAX_PLOT_POINTS)
  plot_data_sample <- umap_df_final[sample_indices, ]
}

umap_plot <- plot_data_sample %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Sample)) +
  geom_point(size = 0.1, alpha = 0.6) +
  scale_color_manual(values = colors) +
  labs(title = paste("UMAP - Colored by Sample", 
                     if(nrow(plot_data_sample) < nrow(umap_df_final)) 
                       paste("(", nrow(plot_data_sample), "of", nrow(umap_df_final), "points)") 
                     else ""),
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 2)))

cat("Saving sample overview plot...\n")
sample_formats <- save_plot_formats(umap_plot, paste0(umap_dir, "/umap_by_sample"))
cat("Sample plot saved in formats:", paste(sample_formats, collapse = ", "), "\n")

print_memory_usage("after sample plot")

# ========================= MEMORY-EFFICIENT MARKER-SPECIFIC PLOTS =========================
cat("Creating marker-specific plots for each sample...\n")

# Load expression data back efficiently for plotting
cat("Loading expression data for marker plots...\n")
if (file.exists(paste0(umap_dir, "/combined_expression.csv"))) {
  # Load only what we need for plotting
  expression_for_plots <- data.table::fread(paste0(umap_dir, "/combined_expression.csv"), 
                                           nThread = min(4, parallel::detectCores()))
  
  # Downsample expression data if needed for memory
  if (nrow(expression_for_plots) > MAX_PLOT_POINTS) {
    cat("Downsampling expression data for plotting...\n")
    sample_indices <- sample(nrow(expression_for_plots), MAX_PLOT_POINTS)
    expression_for_plots <- expression_for_plots[sample_indices, ]
    umap_for_markers <- umap_df_final[sample_indices, ]
  } else {
    umap_for_markers <- umap_df_final
  }
} else {
  cat("Warning: Could not find combined expression data. Skipping marker plots.\n")
  return()
}

print_memory_usage("after loading expression data for plots")

# Get marker names and prepare data
selected_marker_names <- marker_info$markers[marker_indices_of_interest] %>%
  map_chr(~ if(.x == "" || is.na(.x)) paste0("Marker_", which(marker_info$markers == .x)) else .x)

# Add expression data to UMAP dataframe and calculate ranges
umap_with_markers <- umap_for_markers
marker_ranges <- list()

# Process markers in chunks to manage memory
cat("Processing", length(marker_indices_of_interest), "markers...\n")

iwalk(marker_indices_of_interest, ~ {
  marker_name <- selected_marker_names[.y]
  clean_marker_name <- make.names(marker_name)
  
  # Add marker expression data
  umap_with_markers[[clean_marker_name]] <<- expression_for_plots[, .y]
  marker_ranges[[clean_marker_name]] <<- range(expression_for_plots[, .y], na.rm = TRUE)
  
  cat("Marker", clean_marker_name, "range:", 
      marker_ranges[[clean_marker_name]][1], "to", 
      marker_ranges[[clean_marker_name]][2], "\n")
})

# Clean up expression data to free memory
rm(expression_for_plots)
gc()
print_memory_usage("after processing markers")

# Generate plots efficiently using base R
plot_combinations <- expand_grid(
  marker_idx = seq_along(marker_indices_of_interest),
  sample = unique_samples
)

plot_combinations$marker_name <- selected_marker_names[plot_combinations$marker_idx]
plot_combinations$clean_marker_name <- make.names(plot_combinations$marker_name)
plot_combinations$plot_num <- seq_len(nrow(plot_combinations))

total_plots <- nrow(plot_combinations)

cat("Generating", total_plots, "marker-specific plots...\n")

pwalk(plot_combinations, function(marker_idx, sample, marker_name, clean_marker_name, plot_num) {
  
  if (plot_num %% 5 == 0) {  # Progress every 5 plots
    cat("Creating plot", plot_num, "/", total_plots, "for marker:", marker_name, "\n")
    print_memory_usage(paste("plot", plot_num))
  }
  
  # Filter data for this sample using base R
  sample_data <- umap_with_markers[umap_with_markers$Sample == sample, ]
  
  # Further downsample if this sample is still too large
  if (nrow(sample_data) > MAX_PLOT_POINTS / length(unique_samples)) {
    max_points_per_sample <- ceiling(MAX_PLOT_POINTS / length(unique_samples))
    sample_indices <- sample(nrow(sample_data), max_points_per_sample)
    sample_data <- sample_data[sample_indices, ]
  }
  
  # Create plot
  marker_plot <- sample_data %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = .data[[clean_marker_name]])) +
    geom_point(size = 0.1, alpha = 0.7) +
    scale_color_viridis_c(name = marker_name, 
                         limits = marker_ranges[[clean_marker_name]], 
                         option = "viridis") +
    labs(title = paste("UMAP -", marker_name, "expression"),
         subtitle = paste("Sample:", sample, 
                         if(nrow(sample_data) < nrow(umap_with_markers[umap_with_markers$Sample == sample, ]))
                           paste("(", nrow(sample_data), "points)")
                         else ""),
         x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Create safe filename and save
  filename_base <- paste0(umap_dir, "/umap_", 
                         gsub("[^A-Za-z0-9_-]", "_", marker_name), "_",
                         gsub("[^A-Za-z0-9_-]", "_", sample))
  
  marker_formats <- save_plot_formats(marker_plot, filename_base)
  
  if (plot_num %% 5 == 0) {  # Progress every 5 plots
    cat("  Plot", plot_num, "/", total_plots, ":", 
        gsub("[^A-Za-z0-9_-]", "_", marker_name), "in", 
        gsub("[^A-Za-z0-9_-]", "_", sample),
        "(", paste(marker_formats, collapse=", "), ")\n")
  }
  
  # Clean up plot object
  rm(marker_plot, sample_data)
  
  # Force garbage collection every few plots
  if (plot_num %% 10 == 0) {
    gc()
  }
})

print_memory_usage("after all marker plots")

# ========================= COMPLETION SUMMARY =========================
enabled_formats <- names(OUTPUT_FORMATS)[unlist(OUTPUT_FORMATS)]
total_plots <- 1 + (length(marker_indices_of_interest) * length(unique_samples))

print_memory_usage("final")

cat("✓ MEMORY-EFFICIENT UMAP analysis complete!\n")
cat("Results saved in:", umap_dir, "\n")
cat("Files created:\n  - combined_expression.csv\n  - umap_results.csv\n")
cat("  - Sample overview plot (", paste(enabled_formats, collapse = ", "), ")\n")
cat("  -", length(marker_indices_of_interest) * length(unique_samples), 
    "marker-specific plots (", paste(enabled_formats, collapse = ", "), ")\n")
cat("Total plots generated:", total_plots, "plots in", length(enabled_formats), "format(s)\n")

# Final memory report
if (MONITOR_MEMORY) {
  cat("\n=== MEMORY MANAGEMENT REPORT ===\n")
  cat("Memory limits applied:\n")
  cat("  - Max cells per sample:", MAX_CELLS_PER_SAMPLE %||% "unlimited", "\n")
  cat("  - Max total cells for UMAP:", MAX_TOTAL_CELLS_UMAP %||% "unlimited", "\n") 
  cat("  - Max points per plot:", MAX_PLOT_POINTS, "\n")
  cat("  - Final dataset size:", nrow(umap_df_final), "cells\n")
  cat("Memory optimizations used:\n")
  cat("  - ✓ Progressive file loading\n")
  cat("  - ✓ Automatic downsampling\n") 
  cat("  - ✓ Garbage collection\n")
  cat("  - ✓ Plot-level data reduction\n")
  print_memory_usage("final cleanup")
}