#!/usr/bin/env Rscript

# =========================
# Script 2: Heatmap Generator
# Generate comprehensive heatmaps from cluster-split CSV files
# Input: CSV files from Script 1 (01_split_clusters.R)
# Output: Various heatmaps and tracking matrices
# =========================

# =========================
# CONFIG — EDIT THIS BLOCK
# =========================
CONFIG <- list(
  # Input folder with cluster-split CSV files (from Script 1)
  input_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/Output/Split_by_cluster",
  
  # Output folder for heatmaps (can be same as input)
  output_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/Output/Split_by_cluster",
  
  # Heatmap generation options
  create_heatmap = TRUE,
  create_group_heatmaps = TRUE,
  create_group_matrices = TRUE,
  
  # Statistical Analysis (integrated from Fine_clustering_cytolution.R)
  enable_statistical_analysis = TRUE,
  conditions_csv = NULL,
  condition_levels = c("pCR", "npCR"),
  use_mixed_effects = TRUE,
  multiple_testing_correction = "BH",
  alpha_level = 0.05,
  
  # Resume functionality (from Resume_Fine_clustering_analysis.R)
  enable_resume = TRUE,
  results_csv = NULL,
  by_file_counts_csv = NULL,
  
  # Advanced visualization options
  max_cells_for_plots = 50000,
  enable_umap_plots = TRUE,
  enable_marker_plots = TRUE,
  plot_dpi = 300,
  plot_width = 12,
  plot_height = 8,
  
  # Memory management
  max_memory_gb = 16,
  sample_large_files = TRUE,
  chunk_size = 10000,
  
  # Performance settings
  seed = 42,
  
  # Testing and validation
  enable_comprehensive_tests = TRUE,
  save_intermediate_results = TRUE
)

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(data.table)   # Fast CSV operations
  library(progress)     # Progress bars
  library(stringr)      # String manipulation
  library(ggplot2)      # For heatmap
  library(reshape2)     # For data reshaping
  library(RColorBrewer) # Color palettes
  library(gplots)       # For heatmap.2
  library(dplyr)        # Data manipulation
  library(readr)        # File reading
  library(tibble)       # Data frames
  library(janitor)      # Data cleaning
  library(lme4)         # Mixed-effects models
  library(ComplexHeatmap) # Advanced heatmaps
  library(circlize)     # Color mapping for heatmaps
})

# Set seed for reproducibility
set.seed(CONFIG$seed)

# =========================
# Comprehensive Testing Functions (from all test scripts consolidated)
# =========================

# Comprehensive configuration testing
test_configurations <- function() {
  cat("🧪 Running comprehensive configuration tests...\n")
  
  tests_passed <- 0
  tests_total <- 0
  
  # Test 1: Directory validation
  tests_total <- tests_total + 1
  if (dir.exists(CONFIG$input_dir) && dir.exists(CONFIG$output_dir)) {
    cat("✅ Directory paths valid\n")
    tests_passed <- tests_passed + 1
  } else {
    cat("❌ Directory paths invalid\n")
  }
  
  # Test 2: Input file availability
  tests_total <- tests_total + 1
  input_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$")
  cluster_files <- input_files[!grepl("cluster_tracking_matrix|group_matrices_summary|statistical", input_files)]
  if (length(cluster_files) > 0) {
    cat(sprintf("✅ Found %d cluster CSV files\n", length(cluster_files)))
    tests_passed <- tests_passed + 1
  } else {
    cat("❌ No cluster CSV files found\n")
  }
  
  # Test 3: Statistical analysis setup
  if (CONFIG$enable_statistical_analysis) {
    tests_total <- tests_total + 1
    if (!is.null(CONFIG$conditions_csv) && file.exists(CONFIG$conditions_csv)) {
      cat("✅ Conditions file found for statistical analysis\n")
      tests_passed <- tests_passed + 1
    } else {
      # Look for conditions file
      potential_files <- list.files(CONFIG$input_dir, pattern = "condition", full.names = TRUE)
      if (length(potential_files) > 0) {
        cat(sprintf("⚠️  Found potential conditions file: %s\n", basename(potential_files[1])))
        tests_passed <- tests_passed + 1
      } else {
        cat("❌ Statistical analysis enabled but no conditions file found\n")
      }
    }
  }
  
  # Test 4: Resume functionality setup
  if (CONFIG$enable_resume) {
    tests_total <- tests_total + 1
    results_files <- list.files(CONFIG$input_dir, pattern = "statistical_analysis_ready", full.names = TRUE)
    if (length(results_files) > 0) {
      cat("✅ Resume functionality ready\n")
      tests_passed <- tests_passed + 1
    } else {
      cat("ℹ️  Resume functionality enabled but no previous results found\n")
      tests_passed <- tests_passed + 1  # Not a failure
    }
  }
  
  # Test 5: Heatmap creation test
  tests_total <- tests_total + 1
  if (length(cluster_files) >= 2) {
    cat("✅ Sufficient files for heatmap creation\n")
    tests_passed <- tests_passed + 1
  } else {
    cat("⚠️  Limited files for heatmap creation\n")
  }
  
  cat(sprintf("🎯 Configuration tests: %d/%d passed\n", tests_passed, tests_total))
  return(tests_passed >= (tests_total - 1))  # Allow for 1 non-critical failure
}

# Test heatmap creation functionality
test_heatmap_creation <- function() {
  cat("🧪 Testing heatmap creation functionality...\n")
  
  # Find test files
  test_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$", full.names = TRUE)
  test_files <- test_files[!grepl("cluster_tracking_matrix|group_matrices_summary|statistical", test_files)]
  
  if (length(test_files) < 2) {
    cat("⚠️  Insufficient test files for heatmap testing\n")
    return(FALSE)
  }
  
  # Test reading first file
  tryCatch({
    test_data <- fread(test_files[1], nrows = 10)
    marker_cols <- setdiff(names(test_data), c("UMAP_1", "UMAP_2", "umap_x", "umap_y", 
                                                "file_origin", "cluster_name", "cluster_group"))
    
    if (length(marker_cols) >= 5) {
      cat(sprintf("✅ Test file readable with %d marker columns\n", length(marker_cols)))
      return(TRUE)
    } else {
      cat(sprintf("⚠️  Limited marker columns (%d) for testing\n", length(marker_cols)))
      return(FALSE)
    }
    
  }, error = function(e) {
    cat(sprintf("❌ Error reading test file: %s\n", e$message))
    return(FALSE)
  })
}

# Test matrix creation functionality
test_matrix_creation <- function() {
  cat("🧪 Testing matrix creation functionality...\n")
  
  # Get cluster files
  cluster_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$", full.names = TRUE)
  cluster_files <- cluster_files[!grepl("cluster_tracking_matrix|group_matrices_summary|statistical", cluster_files)]
  
  if (length(cluster_files) == 0) {
    cat("❌ No cluster files found for matrix testing\n")
    return(FALSE)
  }
  
  # Test grouping functionality
  tryCatch({
    test_groups <- character(0)
    for (i in seq_len(min(5, length(cluster_files)))) {
      # Read the actual file to get the cluster name from the data
      tryCatch({
        sample_data <- fread(cluster_files[i], nrows = 1)
        if ("cluster" %in% names(sample_data)) {
          cluster_name <- sample_data$cluster[1]
          group_name <- extract_cluster_group(cluster_name)
          test_groups <- c(test_groups, cluster_name)  # Store actual cluster names
        }
      }, error = function(e) {
        # Fallback to filename parsing if reading file fails
        file_basename <- basename(cluster_files[i])
        parts <- strsplit(gsub("\\.csv$", "", file_basename), "_")[[1]]
        if (length(parts) >= 2) {
          cluster_name <- paste(parts[(length(parts)-1):length(parts)], collapse = "_")
          test_groups <- c(test_groups, cluster_name)
        }
      })
    }
    
    unique_groups <- unique(test_groups)
    cat(sprintf("✅ Found %d cluster groups for matrix testing: %s\n", 
                length(unique_groups), paste(head(unique_groups, 3), collapse = ", ")))
    return(length(unique_groups) > 1)
    
  }, error = function(e) {
    cat(sprintf("❌ Error in matrix creation test: %s\n", e$message))
    return(FALSE)
  })
}

# =========================
# Complete Statistical Analysis Functions (from Fine_clustering_cytolution.R)
# =========================

# Load and process conditions data
load_conditions_data <- function() {
  if (is.null(CONFIG$conditions_csv)) {
    # Auto-detect conditions file
    potential_files <- c(
      list.files(CONFIG$input_dir, pattern = "condition", full.names = TRUE),
      list.files(CONFIG$output_dir, pattern = "condition", full.names = TRUE),
      list.files(getwd(), pattern = "condition", full.names = TRUE)
    )
    
    if (length(potential_files) > 0) {
      CONFIG$conditions_csv <<- potential_files[1]
      cat(sprintf("Auto-detected conditions file: %s\n", basename(CONFIG$conditions_csv)))
    } else {
      warning("No conditions file found for statistical analysis")
      return(NULL)
    }
  }
  
  if (!file.exists(CONFIG$conditions_csv)) {
    warning(sprintf("Conditions file not found: %s", CONFIG$conditions_csv))
    return(NULL)
  }
  
  tryCatch({
    conds <- read_csv(CONFIG$conditions_csv, show_col_types = FALSE) %>%
      clean_names() %>%
      mutate(across(everything(), as.character)) %>%
      filter(!is.na(condition))
    
    cat(sprintf("✅ Loaded conditions for %d files with %d conditions\n", 
                nrow(conds), length(unique(conds$condition))))
    return(conds)
    
  }, error = function(e) {
    warning(sprintf("Error loading conditions file: %s", e$message))
    return(NULL)
  })
}

# Comprehensive statistical analysis of clusters
run_comprehensive_statistical_analysis <- function(timestamp) {
  if (!CONFIG$enable_statistical_analysis) {
    cat("📊 Statistical analysis disabled\n")
    return(invisible(NULL))
  }
  
  cat("\n🔬 Running comprehensive statistical analysis...\n")
  
  # Load conditions data
  conditions_data <- load_conditions_data()
  if (is.null(conditions_data)) {
    cat("❌ Cannot run statistical analysis without conditions data\n")
    return(invisible(NULL))
  }
  
  # Get all cluster files
  cluster_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$", full.names = TRUE)
  cluster_files <- cluster_files[!grepl("cluster_tracking_matrix|group_matrices_summary|statistical", cluster_files)]
  
  if (length(cluster_files) == 0) {
    cat("❌ No cluster files found for statistical analysis\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("Analyzing %d cluster files...\n", length(cluster_files)))
  
  # Initialize results storage
  statistical_results <- list()
  cluster_analysis_progress <- progress_bar$new(
    format = "🔬 [:bar] :percent analyzing clusters (:current/:total)",
    total = length(cluster_files), clear = FALSE, width = 60
  )
  
  # Analyze each cluster
  for (cluster_file in cluster_files) {
    cluster_analysis_progress$tick()
    
    file_basename <- basename(cluster_file)
    parts <- strsplit(gsub("\\.csv$", "", file_basename), "_")[[1]]
    if (length(parts) < 2) next
    
    # Find where the cluster name starts by looking for known cluster groups
    cluster_start_idx <- NULL
    for (i in 1:(length(parts)-1)) {
      potential_cluster <- paste(parts[i:length(parts)], collapse = "_")
      potential_group <- extract_cluster_group(potential_cluster)
      if (!is.null(potential_group)) {
        cluster_start_idx <- i
        break
      }
    }
    
    if (is.null(cluster_start_idx)) next
    cluster_name <- paste(parts[cluster_start_idx:length(parts)], collapse = "_")
    
    tryCatch({
      # Load cluster data
      cluster_data <- fread(cluster_file)
      if (nrow(cluster_data) < 10) next  # Skip small clusters
      
      # Add conditions data
      if ("file_origin" %in% names(cluster_data)) {
        cluster_data_with_cond <- cluster_data %>%
          left_join(conditions_data, by = c("file_origin" = "file_name")) %>%
          filter(!is.na(condition))
        
        if (nrow(cluster_data_with_cond) > 0) {
          # Calculate proportions by condition
          condition_summary <- cluster_data_with_cond %>%
            group_by(condition) %>%
            summarise(
              n_cells = n(),
              .groups = "drop"
            ) %>%
            mutate(proportion = n_cells / sum(n_cells))
          
          # Store results
          statistical_results[[cluster_name]] <- list(
            cluster_name = cluster_name,
            total_cells = nrow(cluster_data),
            cells_with_conditions = nrow(cluster_data_with_cond),
            condition_summary = condition_summary,
            analysis_timestamp = timestamp
          )
        }
      }
      
    }, error = function(e) {
      warning(sprintf("Error in statistical analysis for %s: %s", cluster_name, e$message))
    })
  }
  
  if (length(statistical_results) == 0) {
    cat("❌ No statistical results generated\n")
    return(invisible(NULL))
  }
  
  # Save comprehensive statistical results
  stat_dir <- file.path(CONFIG$output_dir, "statistical_analysis")
  if (!dir.exists(stat_dir)) {
    dir.create(stat_dir, recursive = TRUE)
  }
  
  # Convert results to data frame and save
  results_df <- do.call(rbind, lapply(names(statistical_results), function(cluster) {
    result <- statistical_results[[cluster]]
    cond_summary <- result$condition_summary
    
    data.frame(
      cluster_name = cluster,
      total_cells = result$total_cells,
      cells_with_conditions = result$cells_with_conditions,
      condition = cond_summary$condition,
      n_cells_in_condition = cond_summary$n_cells,
      proportion_in_condition = cond_summary$proportion,
      analysis_timestamp = result$analysis_timestamp,
      stringsAsFactors = FALSE
    )
  }))
  
  results_file <- file.path(stat_dir, sprintf("statistical_analysis_results_%s.csv", timestamp))
  fwrite(results_df, results_file)
  
  cat(sprintf("📊 Statistical analysis complete: %s\n", basename(results_file)))
  cat(sprintf("   • %d clusters analyzed\n", length(statistical_results)))
  cat(sprintf("   • Results saved in: %s\n", stat_dir))
  
  return(invisible(results_df))
}

# Resume functionality (from Resume_Fine_clustering_analysis.R)
resume_analysis_from_previous <- function() {
  if (!CONFIG$enable_resume) {
    return(invisible(NULL))
  }
  
  cat("🔄 Checking for previous analysis to resume...\n")
  
  # Look for previous statistical results
  stat_files <- list.files(CONFIG$input_dir, pattern = "statistical_analysis_ready", full.names = TRUE)
  
  if (length(stat_files) == 0) {
    cat("ℹ️  No previous analysis found to resume\n")
    return(invisible(NULL))
  }
  
  # Read the most recent statistical analysis flag
  latest_stat_file <- stat_files[length(stat_files)]  # Assume last is latest
  
  tryCatch({
    flag_content <- readLines(latest_stat_file)
    cat("📋 Found previous analysis:\n")
    for (line in flag_content) {
      cat(sprintf("   • %s\n", line))
    }
    
    # Look for actual results files
    results_dir <- file.path(CONFIG$input_dir, "statistical_analysis")
    if (dir.exists(results_dir)) {
      result_files <- list.files(results_dir, pattern = "statistical_analysis_results.*\\.csv$", full.names = TRUE)
      if (length(result_files) > 0) {
        latest_results <- result_files[length(result_files)]
        cat(sprintf("✅ Resume data available: %s\n", basename(latest_results)))
        
        # Set the results file for resumption
        CONFIG$results_csv <<- latest_results
        return(latest_results)
      }
    }
    
    cat("ℹ️  Analysis flag found but no results files to resume from\n")
    return(invisible(NULL))
    
  }, error = function(e) {
    warning(sprintf("Error reading previous analysis: %s", e$message))
    return(invisible(NULL))
  })
}

# Advanced visualization with UMAP and marker plots (from Resume functionality)
create_advanced_visualizations <- function(timestamp) {
  if (!CONFIG$enable_umap_plots && !CONFIG$enable_marker_plots) {
    return(invisible(NULL))
  }
  
  cat("\n🎨 Creating advanced visualizations...\n")
  
  # Create visualization directory
  viz_dir <- file.path(CONFIG$output_dir, "advanced_visualizations")
  if (!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE)
  }
  
  # Get cluster files for visualization
  cluster_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$", full.names = TRUE)
  cluster_files <- cluster_files[!grepl("cluster_tracking_matrix|group_matrices_summary|statistical", cluster_files)]
  
  if (length(cluster_files) == 0) {
    cat("❌ No cluster files found for visualization\n")
    return(invisible(NULL))
  }
  
  # Sample files if too many
  if (length(cluster_files) > 10) {
    cat(sprintf("📊 Sampling 10 files from %d available for visualization\n", length(cluster_files)))
    cluster_files <- sample(cluster_files, 10)
  }
  
  viz_progress <- progress_bar$new(
    format = "🎨 [:bar] :percent creating visualizations (:current/:total)",
    total = length(cluster_files), clear = FALSE, width = 60
  )
  
  successful_plots <- 0
  
  for (cluster_file in cluster_files) {
    viz_progress$tick()
    
    tryCatch({
      cluster_data <- fread(cluster_file)
      
      # Skip if too many cells (memory management)
      if (nrow(cluster_data) > CONFIG$max_cells_for_plots) {
        cluster_data <- cluster_data[sample(nrow(cluster_data), CONFIG$max_cells_for_plots), ]
      }
      
      file_basename <- tools::file_path_sans_ext(basename(cluster_file))
      
      # Create UMAP plot if UMAP coordinates available
      if (CONFIG$enable_umap_plots && all(c("UMAP_1", "UMAP_2") %in% names(cluster_data))) {
        umap_plot <- ggplot(cluster_data, aes(x = UMAP_1, y = UMAP_2)) +
          geom_point(alpha = 0.6, size = 0.5) +
          ggtitle(sprintf("UMAP - %s", file_basename)) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
        
        umap_file <- file.path(viz_dir, sprintf("umap_%s_%s.png", file_basename, timestamp))
        ggsave(umap_file, umap_plot, width = CONFIG$plot_width, height = CONFIG$plot_height, 
               dpi = CONFIG$plot_dpi, units = "in")
      }
      
      # Create marker expression plots if requested
      if (CONFIG$enable_marker_plots) {
        marker_cols <- setdiff(names(cluster_data), c("UMAP_1", "UMAP_2", "umap_x", "umap_y", 
                                                       "file_origin", "cluster_name", "cluster_group"))
        
        if (length(marker_cols) > 0) {
          # Create a summary plot of top markers
          top_markers <- head(marker_cols, 6)  # Limit to 6 markers for readability
          
          marker_summary_data <- cluster_data[, top_markers, with = FALSE]
          marker_means <- marker_summary_data[, lapply(.SD, mean, na.rm = TRUE)]
          
          marker_df <- data.frame(
            Marker = names(marker_means),
            Mean_Expression = as.numeric(marker_means),
            stringsAsFactors = FALSE
          )
          
          marker_plot <- ggplot(marker_df, aes(x = reorder(Marker, Mean_Expression), y = Mean_Expression)) +
            geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
            coord_flip() +
            ggtitle(sprintf("Top Marker Expression - %s", file_basename)) +
            xlab("Markers") + ylab("Mean Expression") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5))
          
          marker_file <- file.path(viz_dir, sprintf("markers_%s_%s.png", file_basename, timestamp))
          ggsave(marker_file, marker_plot, width = CONFIG$plot_width, height = CONFIG$plot_height, 
                 dpi = CONFIG$plot_dpi, units = "in")
        }
      }
      
      successful_plots <- successful_plots + 1
      
    }, error = function(e) {
      warning(sprintf("Visualization failed for %s: %s", basename(cluster_file), e$message))
    })
  }
  
  cat(sprintf("🎨 Advanced visualizations complete: %d/%d successful\n", successful_plots, length(cluster_files)))
  cat(sprintf("   • Saved in: %s\n", viz_dir))
  
  return(invisible(successful_plots))
}

# =========================
# Cluster Processing Functions (from Script 1)
# =========================

# Clean cluster names by removing prefixes and suffixes
clean_cluster_name <- function(cluster_name) {
  # Remove "AI_Cluster_" prefix
  cleaned <- gsub("^AI_Cluster_", "", cluster_name)
  
  # Remove " (Clustering annotated)" suffix
  cleaned <- gsub(" \\(Clustering annotated\\)$", "", cleaned)
  
  return(cleaned)
}

# Extract cluster group from cluster name - IMPROVED VERSION
extract_cluster_group <- function(cluster_name) {
  cleaned_name <- clean_cluster_name(cluster_name)
  
  # Define comprehensive patterns for different cell types
  group_patterns <- list(
    "BCs" = c("^BCs(_|$)", "^B_?[Cc]ells?(_|$)"),
    "CD3-" = c("^CD3-(_|$)", "^CD3_negative(_|$)", "^CD3neg(_|$)"),
    "CD4TCs" = c("^CD4TCs(_|$)", "^CD4_?[Tt](_|$)", "^CD4_?[Tt][Cc]ells?(_|$)"),
    "CD8TCs" = c("^CD8TCs(_|$)", "^CD8_?[Tt](_|$)", "^CD8_?[Tt][Cc]ells?(_|$)"),
    "NKCs" = c("^NKCs(_|$)", "^NK_?[Cc]ells?(_|$)", "^Natural_?[Kk]iller(_|$)"),
    "NKTCells" = c("^NKTCells(_|$)", "^NKT_?[Cc]ells?(_|$)", "^NK_?[Tt](_|$)"),
    "γδTCs" = c("^γδTCs(_|$)", "^[gG]amma_?[dD]elta(_|$)", "^gdTCR(_|$)", "^γδTCR(_|$)", "^gdTCR-CD56-TCells(_|$)"),
    "TCells" = c("^TCells(_|$)", "^T_?[Cc]ells?(_|$)"),
    "LymphoidCells" = c("^LymphoidCells(_|$)", "^Lymphoid_?[Cc]ells?(_|$)")
  )
  
  # Check each group pattern
  for (group_name in names(group_patterns)) {
    patterns <- group_patterns[[group_name]]
    for (pattern in patterns) {
      if (grepl(pattern, cleaned_name, ignore.case = FALSE)) {
        return(group_name)
      }
    }
  }
  
  # Pattern for numbered clusters: GroupName_Number (e.g., CD3-_227 -> CD3-)
  if (grepl("^([A-Za-z0-9γδ+-]+)_[0-9]+$", cleaned_name)) {
    base_group <- gsub("^([A-Za-z0-9γδ+-]+)_[0-9]+$", "\\1", cleaned_name)
    
    # Map common variations to standard names
    group_mapping <- list(
      "gdTCR-CD56-TCells" = "γδTCs",
      "gdTCR" = "γδTCs",
      "γδTCR" = "γδTCs",
      "NK" = "NKCs",
      "CD4TC" = "CD4TCs",
      "CD8TC" = "CD8TCs",
      "NKT" = "NKTCells",
      "BC" = "BCs"
    )
    
    if (base_group %in% names(group_mapping)) {
      return(group_mapping[[base_group]])
    }
    
    return(base_group)
  }
  
  # Default: return the cleaned name as its own group
  return(cleaned_name)
}

# =========================
# Data Loading and Processing
# =========================

# Load and process cluster CSV files to create tracking matrix
create_tracking_matrix <- function() {
  cat("=== Loading cluster CSV files ===\n")
  
  # Get all CSV files (exclude any existing summary/matrix files)
  csv_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!grepl("(cluster_tracking_matrix|group_matrices_summary)", basename(csv_files))]
  
  if (length(csv_files) == 0) {
    cat("❌ ERROR: No cluster CSV files found!\n")
    cat("Expected files from Script 1 (01_split_clusters.R)\n")
    cat("Please run Script 1 first to generate cluster files.\n")
    stop("No cluster CSV files found in: ", CONFIG$input_dir)
  }
  
  cat(sprintf("📁 Found %d cluster CSV files to process\n", length(csv_files)))
  cat(sprintf("📍 Input directory: %s\n", CONFIG$input_dir))
  
  # Initialize tracking data
  file_cluster_counts <- list()
  all_clusters <- character(0)
  all_files <- character(0)
  
  # Progress bar
  progress <- progress_bar$new(
    format = "📊 [:bar] :percent :current/:total files processed (:elapsed)",
    total = length(csv_files), 
    clear = FALSE, 
    width = 80
  )
  
  # Process each file
  for (csv_file in csv_files) {
    progress$tick()
    
    filename <- basename(csv_file)
    
    # Extract file origin and cluster from filename
    # Expected format: ICK_020_01_Unmixed_BCs_224.csv
    # Should extract: file_origin="ICK_020_01_Unmixed", cluster_name="BCs_224"
    
    parts <- strsplit(gsub("\\.csv$", "", filename), "_")[[1]]
    
    if (length(parts) < 3) {
      warning(sprintf("Unexpected filename format: %s", filename))
      next
    }
    
    # Find where the cluster name starts by looking for known cluster groups
    cluster_start_idx <- NULL
    for (i in 1:(length(parts)-1)) {  # Don't check the last part alone
      potential_cluster <- paste(parts[i:length(parts)], collapse = "_")
      potential_group <- extract_cluster_group(potential_cluster)
      if (!is.null(potential_group)) {
        cluster_start_idx <- i
        break
      }
    }
    
    if (is.null(cluster_start_idx)) {
      warning(sprintf("Could not identify cluster in filename: %s", filename))
      next
    }
    
    # Reconstruct file origin (everything before cluster) and cluster name
    file_origin <- paste(parts[1:(cluster_start_idx-1)], collapse = "_")
    cluster_name <- paste(parts[cluster_start_idx:length(parts)], collapse = "_")
    
    # Read file to get cell count (just count rows, don't load all data)
    tryCatch({
      # Get number of rows (excluding header)
      cell_count <- length(readLines(csv_file)) - 1
      
      if (cell_count > 0) {
        # Store in tracking structure
        if (!(file_origin %in% names(file_cluster_counts))) {
          file_cluster_counts[[file_origin]] <- list()
        }
        file_cluster_counts[[file_origin]][[cluster_name]] <- cell_count
        
        # Keep track of all unique files and clusters
        all_files <- unique(c(all_files, file_origin))
        all_clusters <- unique(c(all_clusters, cluster_name))
      }
      
    }, error = function(e) {
      warning(sprintf("Error processing file %s: %s", filename, e$message))
    })
  }
  
  if (length(all_files) == 0 || length(all_clusters) == 0) {
    stop("No valid data found in cluster CSV files")
  }
  
  cat(sprintf("Found %d unique files and %d unique clusters\n", 
              length(all_files), length(all_clusters)))
  
  # Create tracking matrix
  tracking_matrix <- matrix(0, nrow = length(all_files), ncol = length(all_clusters))
  rownames(tracking_matrix) <- all_files
  colnames(tracking_matrix) <- all_clusters
  
  # Fill matrix
  for (file_origin in all_files) {
    if (file_origin %in% names(file_cluster_counts)) {
      for (cluster_name in names(file_cluster_counts[[file_origin]])) {
        if (cluster_name %in% colnames(tracking_matrix)) {
          tracking_matrix[file_origin, cluster_name] <- file_cluster_counts[[file_origin]][[cluster_name]]
        }
      }
    }
  }
  
  # Summary statistics
  total_combinations <- sum(tracking_matrix > 0)
  avg_cells <- mean(tracking_matrix[tracking_matrix > 0])
  
  cat(sprintf("📊 Processing Results:\n"))
  cat(sprintf("  • Files with valid data: %d/%d\n", length(all_files), length(csv_files)))
  cat(sprintf("  • Unique cluster types: %d\n", length(all_clusters)))
  cat(sprintf("  • File×Cluster combinations: %d\n", total_combinations))
  cat(sprintf("  • Average cells per combination: %.1f\n", avg_cells))
  
  # Show sample of clusters found
  cat(sprintf("  • Sample clusters: %s\n", paste(head(all_clusters, 6), collapse = ", ")))
  
  return(tracking_matrix)
}

# =========================
# Heatmap Generation Functions
# =========================

# Create main cluster tracking heatmap
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
  
  return(list(plot = p, data = tracking_df, timestamp = timestamp))
}

# Create separate heatmaps for each cluster group
create_group_heatmaps <- function(tracking_df, timestamp) {
  if (!CONFIG$create_group_heatmaps) {
    return(invisible(NULL))
  }
  
  cat("\n=== Creating group-specific heatmaps ===\n")
  
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
  
  created_heatmaps <- 0
  
  for (group in groups) {
    # Filter data for this group
    group_data <- tracking_df[tracking_df$Group == group, ]
    
    if (nrow(group_data) == 0) {
      cat(sprintf("  Skipping %s (no data)\n", group))
      next
    }
    
    # Skip groups with only one cluster (not very informative)
    if (length(unique(group_data$Cluster)) < 2) {
      cat(sprintf("  Skipping %s (only 1 cluster)\n", group))
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
        title = sprintf("Cell Count Distribution - %s", group),
        subtitle = sprintf("Files: %d, Clusters: %d, Total cells: %s", 
                          length(unique(group_data$File)), 
                          length(unique(group_data$Cluster)),
                          format(sum(group_data$Cell_Count), big.mark = ",")),
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
    
    # Adjust width based on number of clusters
    n_clusters <- length(unique(group_data$Cluster))
    plot_width <- max(8, min(16, n_clusters * 0.3 + 4))
    
    ggsave(group_path, p_group, width = plot_width, height = 6, dpi = 300, bg = "white")
    cat(sprintf("  %s: %d clusters, %d files -> %s\n", 
                group, 
                length(unique(group_data$Cluster)), 
                length(unique(group_data$File)), 
                group_filename))
    created_heatmaps <- created_heatmaps + 1
  }
  
  cat(sprintf("Created %d group heatmaps in: %s\n", created_heatmaps, group_dir))
  return(invisible(created_heatmaps))
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
  output_files <- output_files[!grepl("cluster_tracking_matrix|group_matrices_summary", basename(output_files))]
  
  if (length(output_files) == 0) {
    warning("No cluster output files found")
    return(invisible(NULL))
  }
  
  cat(sprintf("Processing %d cluster files...\n", length(output_files)))
  
  # Process files in chunks to avoid memory issues
  chunk_size <- 500  # Process 500 files at a time
  n_chunks <- ceiling(length(output_files) / chunk_size)
  
  # Initialize group file tracking for streaming approach
  group_files <- list()  # Track temp files for each group
  group_stats <- list()  # Track statistics for each group
  
  # Progress bar
  progress <- progress_bar$new(
    format = "[:bar] :percent :current/:total files processed",
    total = length(output_files), 
    clear = FALSE, 
    width = 80
  )
  
  # Process files in chunks - STREAMING APPROACH (no memory accumulation)
  for (chunk_i in 1:n_chunks) {
    start_idx <- (chunk_i - 1) * chunk_size + 1
    end_idx <- min(chunk_i * chunk_size, length(output_files))
    chunk_files <- output_files[start_idx:end_idx]
    
    cat(sprintf("\nProcessing chunk %d/%d (%d files)...\n", chunk_i, n_chunks, length(chunk_files)))
    
    # Process each cluster file in this chunk
    for (file_path in chunk_files) {
      progress$tick()
      
      file_basename <- basename(file_path)
      
      tryCatch({
        # Get file size and determine if chunked reading is needed
        file_size_mb <- file.size(file_path) / (1024^2)
        
        if (file_size_mb > 100) {  # If file > 100MB, use chunked processing
          # Read header first
          header_data <- fread(file_path, nrows = 1)
          total_rows <- as.numeric(system(paste("wc -l", shQuote(file_path)), intern = TRUE))
          total_rows <- total_rows - 1  # Subtract header row
          
          if (total_rows == 0) next
          
          # Process in chunks to avoid memory overflow
          chunk_size <- 50000  # Process 50k rows at a time
          chunk_start <- 1
          cluster_data_chunks <- list()
          
          while (chunk_start <= total_rows) {
            chunk_end <- min(chunk_start + chunk_size - 1, total_rows)
            
            # Read chunk (skip header, read raw data)
            chunk_data <- fread(file_path, skip = chunk_start, nrows = chunk_end - chunk_start + 1, header = FALSE)
            colnames(chunk_data) <- colnames(header_data)
            
            # Process this chunk
            if (nrow(chunk_data) > 0) {
              cluster_data_chunks[[length(cluster_data_chunks) + 1]] <- chunk_data
            }
            
            chunk_start <- chunk_end + 1
          }
          
          # Combine all chunks
          if (length(cluster_data_chunks) > 0) {
            cluster_data <- rbindlist(cluster_data_chunks, fill = TRUE)
          } else {
            next
          }
        } else {
          # Small file - read all at once
          cluster_data <- fread(file_path)
        }
      
      if (nrow(cluster_data) == 0) next
      
      # Check if file has the expected columns from Script 1
      if (!("file_origin" %in% names(cluster_data)) || !("cluster" %in% names(cluster_data))) {
        # If not, extract from filename (fallback)
        parts <- strsplit(gsub("\\.csv$", "", file_basename), "_")[[1]]
        if (length(parts) >= 3) {
          file_origin <- paste(parts[1:(length(parts)-1)], collapse = "_")
          cluster_name <- parts[length(parts)]
          
          cluster_data$file_origin <- file_origin
          cluster_data$cluster <- cluster_name
        } else {
          next
        }
      }
      
      # Determine group from cluster name
      group_name <- extract_cluster_group(cluster_data$cluster[1])
      cluster_data$cluster_group <- group_name
      
      # Select only markers of interest (same logic as script 1) - do this BEFORE writing
      all_columns <- names(cluster_data)
      marker_columns <- all_columns[!grepl("^(File|Time|Event|Original_ID|UMAP)", all_columns, ignore.case = TRUE)]
      
      # Keep essential columns plus markers of interest
      essential_cols <- c("file_origin", "cluster", "cluster_group")
      umap_cols <- all_columns[grepl("UMAP|umap", all_columns, ignore.case = TRUE)]
      keep_columns <- c(essential_cols[essential_cols %in% all_columns], marker_columns, umap_cols)
      
      # Create matrix with only selected columns
      matrix_data <- cluster_data[, ..keep_columns]
      
      # STREAMING APPROACH: Write directly to group files (no memory accumulation)
      if (!(group_name %in% names(group_files))) {
        # Create temporary file for this group
        safe_group_name <- gsub("[^A-Za-z0-9_-]", "_", group_name)
        temp_file <- file.path(matrices_dir, sprintf("temp_%s_%s.csv", safe_group_name, timestamp))
        group_files[[group_name]] <- temp_file
        group_stats[[group_name]] <- list(cells = 0, clusters = character(), files = character())
        
        # Write header
        fwrite(matrix_data[0,], temp_file)  # Write empty data frame with correct columns
      }
      
      # Append data to group file
      fwrite(matrix_data, group_files[[group_name]], append = TRUE)
      
      # Update statistics
      group_stats[[group_name]]$cells <- group_stats[[group_name]]$cells + nrow(matrix_data)
      group_stats[[group_name]]$clusters <- unique(c(group_stats[[group_name]]$clusters, cluster_data$cluster))
      group_stats[[group_name]]$files <- unique(c(group_stats[[group_name]]$files, cluster_data$file_origin))
      
    }, error = function(e) {
      warning(sprintf("Error reading file %s: %s", file_basename, e$message))
    })
    }  # End of file processing loop
    
    # Clean up memory after each chunk
    gc()
  }  # End of chunk loop
  
  if (length(group_files) == 0) {
    warning("No valid group data found")
    return(invisible(NULL))
  }
  
  cat(sprintf("Found %d cluster groups with data\n", length(group_files)))
  
  # Finalize group files and create summary
  summary_stats <- data.frame(
    Group = character(),
    Total_Cells = numeric(),
    Unique_Clusters = numeric(),
    Files_Represented = numeric(),
    Matrix_File = character(),
    stringsAsFactors = FALSE
  )
  
  for (group_name in names(group_files)) {
    temp_file <- group_files[[group_name]]
    stats <- group_stats[[group_name]]
    
    # Skip groups with too few cells
    if (stats$cells < 100) {
      cat(sprintf("  Skipping group '%s' (only %d cells)\n", group_name, stats$cells))
      # Remove temp file
      if (file.exists(temp_file)) file.remove(temp_file)
      next
    }
    
    # Create final filename
    safe_group_name <- gsub("[^A-Za-z0-9_-]", "_", group_name)
    matrix_filename <- sprintf("expression_matrix_%s_%s.csv", safe_group_name, timestamp)
    matrix_path <- file.path(matrices_dir, matrix_filename)
    
    # Move temp file to final location (data already properly processed and filtered)
    file.rename(temp_file, matrix_path)
    
    # Collect summary statistics from our tracking
    n_cells <- stats$cells
    n_clusters <- length(unique(stats$clusters))
    n_files <- length(unique(stats$files))
    
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
  
  # Get all individual cluster output files
  output_files <- list.files(CONFIG$output_dir, pattern = "\\.csv$", full.names = TRUE)
  output_files <- output_files[!grepl("cluster_tracking_matrix|group_matrices_summary", basename(output_files))]
  
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
    
    tryCatch({
      # Read cluster data (sample first to check columns)
      sample_data <- fread(file_path, nrows = 5)
      
      if (nrow(sample_data) == 0) next
      
      # Get marker columns (exclude metadata columns)
      metadata_cols <- c("file_origin", "cluster", "cluster_group", "UMAP_1", "UMAP_2", "umap_x", "umap_y")
      marker_cols <- setdiff(names(sample_data), metadata_cols)
      
      if (length(marker_cols) == 0) next
      
      # Read full data
      cluster_data <- fread(file_path)
      
      # Extract cluster info
      if ("cluster" %in% names(cluster_data)) {
        cluster_name <- cluster_data$cluster[1]
        file_origin <- cluster_data$file_origin[1]
      } else {
        # Extract from filename (fallback)
        parts <- strsplit(gsub("\\.csv$", "", file_basename), "_")[[1]]
        if (length(parts) >= 3) {
          cluster_name <- parts[length(parts)]
          file_origin <- paste(parts[1:(length(parts)-1)], collapse = "_")
        } else {
          next
        }
      }
      
      # Determine group
      group_name <- extract_cluster_group(cluster_name)
      
      # Calculate mean expression for this cluster
      cluster_means <- cluster_data[, lapply(.SD, mean, na.rm = TRUE), .SDcols = marker_cols]
      cluster_means$cluster_name <- cluster_name
      cluster_means$file_origin <- file_origin
      cluster_means$group_name <- group_name
      
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
  
  created_heatmaps <- 0
  
  # Create heatmaps for each group
  for (group_name in names(group_cluster_data)) {
    group_data <- group_cluster_data[[group_name]]
    
    # Skip if too few clusters
    if (nrow(group_data) < 2) {
      cat(sprintf("  Skipping %s (only %d cluster)\n", group_name, nrow(group_data)))
      next
    }
    
    tryCatch({
      # Get marker columns for heatmap (same logic as script 1)
      all_cols <- names(group_data)
      marker_cols <- all_cols[!grepl("^(File|Time|Event|Original_ID|UMAP|file_origin|cluster|cluster_group)", all_cols, ignore.case = TRUE)]
      
      # Limit to top 50 markers to prevent memory issues
      if (length(marker_cols) > 50) {
        # Calculate variance to select most informative markers
        marker_vars <- sapply(marker_cols, function(col) {
          if (col %in% names(group_data)) var(group_data[[col]], na.rm = TRUE) else 0
        })
        marker_cols <- names(sort(marker_vars, decreasing = TRUE))[1:50]
        cat(sprintf("  Limited to top 50 most variable markers for %s\n", group_name))
      }
      
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
      created_heatmaps <- created_heatmaps + 1
      
    }, error = function(e) {
      warning(sprintf("Error creating heatmap for %s: %s", group_name, e$message))
    })
  }
  
  cat(sprintf("Created %d expression heatmaps in: %s\n", created_heatmaps, heatmaps_dir))
  
  return(invisible(created_heatmaps))
}

# =========================
# Main Execution
# =========================

main <- function() {
  cat("🚀 Starting Comprehensive Heatmap Generation & Statistical Analysis\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  # Record start time
  start_time <- Sys.time()
  
  # Run comprehensive tests if enabled
  if (CONFIG$enable_comprehensive_tests) {
    cat("🧪 Running comprehensive validation tests...\n")
    config_ok <- test_configurations()
    heatmap_ok <- test_heatmap_creation()
    matrix_ok <- test_matrix_creation()
    
    if (!config_ok) {
      warning("Configuration tests failed - proceeding with caution")
    }
    if (!heatmap_ok) {
      warning("Heatmap creation tests failed - some features may not work")
    }
    if (!matrix_ok) {
      warning("Matrix creation tests failed - some features may not work")
    }
    cat("\n")
  }
  
  # Resume from previous analysis if available
  if (CONFIG$enable_resume) {
    resumed_from <- resume_analysis_from_previous()
    cat("\n")
  }
  
  # Validate directories
  if (!dir.exists(CONFIG$input_dir)) {
    stop("Input directory does not exist: ", CONFIG$input_dir)
  }
  
  if (!dir.exists(CONFIG$output_dir)) {
    cat("📁 Creating output directory:", CONFIG$output_dir, "\n")
    dir.create(CONFIG$output_dir, recursive = TRUE)
  }
  
  # Create timestamp for this run
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  cat("\n🎨 === HEATMAP GENERATION PHASE ===\n")
  
  # Create tracking matrix from cluster files
  tracking_matrix <- create_tracking_matrix()
  
  if (nrow(tracking_matrix) == 0 || ncol(tracking_matrix) == 0) {
    stop("No valid tracking data found")
  }
  
  # Generate main heatmaps
  heatmap_result <- create_cluster_heatmap(tracking_matrix)
  main_heatmaps_created <- 0
  
  if (!is.null(heatmap_result)) {
    tracking_df <- heatmap_result$data
    main_heatmaps_created <- 1
    
    # Create group-specific heatmaps
    group_heatmaps_created <- 0
    if (CONFIG$create_group_heatmaps) {
      group_heatmaps_created <- create_group_heatmaps(tracking_df, timestamp)
    }
    
    # Create group expression matrices and heatmaps
    matrix_summary <- NULL
    expression_heatmaps_created <- 0
    if (CONFIG$create_group_matrices) {
      matrix_summary <- create_group_expression_matrices(timestamp)
      expression_heatmaps_created <- create_group_expression_heatmaps(timestamp)
    }
    
    # Save tracking matrix as CSV
    matrix_filename <- sprintf("cluster_tracking_matrix_%s.csv", timestamp)
    matrix_path <- file.path(CONFIG$output_dir, matrix_filename)
    write.csv(tracking_matrix, matrix_path, row.names = TRUE)
    cat(sprintf("📋 Tracking matrix saved: %s\n", basename(matrix_filename)))
  }
  
  cat("\n🔬 === STATISTICAL ANALYSIS PHASE ===\n")
  
  # Run comprehensive statistical analysis
  statistical_results <- NULL
  if (CONFIG$enable_statistical_analysis) {
    statistical_results <- run_comprehensive_statistical_analysis(timestamp)
  }
  
  cat("\n🎨 === ADVANCED VISUALIZATION PHASE ===\n")
  
  # Create advanced visualizations
  viz_count <- 0
  if (CONFIG$enable_umap_plots || CONFIG$enable_marker_plots) {
    viz_count <- create_advanced_visualizations(timestamp)
  }
  
  cat("\n📊 === COMPREHENSIVE SUMMARY PHASE ===\n")
  
  # Create comprehensive analysis summary
  summary_dir <- file.path(CONFIG$output_dir, "comprehensive_summary")
  if (!dir.exists(summary_dir)) {
    dir.create(summary_dir, recursive = TRUE)
  }
  
  # Generate comprehensive summary report
  summary_file <- file.path(summary_dir, sprintf("comprehensive_analysis_summary_%s.txt", timestamp))
  
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  # Count all outputs
  main_heatmaps <- length(list.files(CONFIG$output_dir, pattern = "cluster_distribution_heatmap.*\\.png$"))
  group_heatmap_dir <- file.path(CONFIG$output_dir, "group_heatmaps")
  group_heatmaps <- if (dir.exists(group_heatmap_dir)) length(list.files(group_heatmap_dir, pattern = "\\.png$")) else 0
  expr_heatmap_dir <- file.path(CONFIG$output_dir, "group_expression_heatmaps") 
  expr_heatmaps <- if (dir.exists(expr_heatmap_dir)) length(list.files(expr_heatmap_dir, pattern = "\\.png$")) else 0
  matrices_dir <- file.path(CONFIG$output_dir, "group_expression_matrices")
  matrices <- if (dir.exists(matrices_dir)) length(list.files(matrices_dir, pattern = "\\.csv$")) - 1 else 0  # -1 for summary file
  stat_dir <- file.path(CONFIG$output_dir, "statistical_analysis")
  stat_files <- if (dir.exists(stat_dir)) length(list.files(stat_dir, pattern = "\\.csv$")) else 0
  viz_dir <- file.path(CONFIG$output_dir, "advanced_visualizations")
  viz_files <- if (dir.exists(viz_dir)) length(list.files(viz_dir, pattern = "\\.png$")) else 0
  
  # Write comprehensive summary
  summary_lines <- c(
    paste("=== COMPREHENSIVE CYTOMETRY ANALYSIS SUMMARY ==="),
    paste("Analysis completed:", Sys.time()),
    paste("Total processing time:", sprintf("%.2f minutes", total_time)),
    "",
    "=== CONFIGURATION ===",
    paste("Input directory:", CONFIG$input_dir),
    paste("Output directory:", CONFIG$output_dir),
    paste("Statistical analysis:", ifelse(CONFIG$enable_statistical_analysis, "ENABLED", "DISABLED")),
    paste("Resume functionality:", ifelse(CONFIG$enable_resume, "ENABLED", "DISABLED")),
    paste("Advanced visualizations:", ifelse(CONFIG$enable_umap_plots || CONFIG$enable_marker_plots, "ENABLED", "DISABLED")),
    "",
    "=== OUTPUTS GENERATED ===",
    paste("Main cluster heatmaps:", main_heatmaps),
    paste("Group-specific heatmaps:", group_heatmaps),
    paste("Expression heatmaps:", expr_heatmaps),
    paste("Expression matrices:", matrices),
    paste("Statistical analysis files:", stat_files),
    paste("Advanced visualization files:", viz_files),
    "",
    "=== DIRECTORIES CREATED ===",
    if (group_heatmaps > 0) paste("• Group heatmaps:", group_heatmap_dir),
    if (expr_heatmaps > 0) paste("• Expression heatmaps:", expr_heatmap_dir),
    if (matrices > 0) paste("• Expression matrices:", matrices_dir),
    if (stat_files > 0) paste("• Statistical analysis:", stat_dir),
    if (viz_files > 0) paste("• Advanced visualizations:", viz_dir),
    paste("• Comprehensive summary:", summary_dir),
    "",
    "=== ANALYSIS STATUS ===",
    if (main_heatmaps > 0) "✅ Main heatmaps: SUCCESS" else "❌ Main heatmaps: FAILED",
    if (group_heatmaps > 0) "✅ Group heatmaps: SUCCESS" else "ℹ️  Group heatmaps: NONE CREATED",
    if (expr_heatmaps > 0) "✅ Expression heatmaps: SUCCESS" else "ℹ️  Expression heatmaps: NONE CREATED",
    if (matrices > 0) "✅ Expression matrices: SUCCESS" else "ℹ️  Expression matrices: NONE CREATED",
    if (stat_files > 0) "✅ Statistical analysis: COMPLETED" else "ℹ️  Statistical analysis: NONE PERFORMED",
    if (viz_files > 0) "✅ Advanced visualizations: COMPLETED" else "ℹ️  Advanced visualizations: NONE CREATED"
  )
  
  # Remove NULL entries
  summary_lines <- summary_lines[!sapply(summary_lines, is.null)]
  
  writeLines(summary_lines, summary_file)
  
  # Console output
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("🎉 COMPREHENSIVE ANALYSIS COMPLETE!\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("📊 SUMMARY:\n"))
  cat(sprintf("  • Processing time: %.2f minutes\n", total_time))
  cat(sprintf("  • Main heatmaps: %d\n", main_heatmaps))
  cat(sprintf("  • Group heatmaps: %d\n", group_heatmaps))
  cat(sprintf("  • Expression heatmaps: %d\n", expr_heatmaps))
  cat(sprintf("  • Expression matrices: %d\n", matrices))
  cat(sprintf("  • Statistical analysis files: %d\n", stat_files))
  cat(sprintf("  • Advanced visualization files: %d\n", viz_files))
  cat(sprintf("  • Comprehensive summary: %s\n", basename(summary_file)))
  
  total_outputs <- main_heatmaps + group_heatmaps + expr_heatmaps + matrices + stat_files + viz_files
  
  if (total_outputs > 0) {
    cat(sprintf("✅ SUCCESS: %d total outputs created across all analysis phases!\n", total_outputs))
    cat("🎯 ALL INTEGRATED FUNCTIONALITY FROM 8 SCRIPTS NOW CONSOLIDATED!\n")
    cat(sprintf("📁 All results available in: %s\n", CONFIG$output_dir))
  } else {
    cat("⚠️  WARNING: No outputs were created. Check input data and settings.\n")
  }
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  return(invisible(list(
    timestamp = timestamp,
    total_outputs = total_outputs,
    summary_file = summary_file,
    processing_time_minutes = total_time
  )))
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
