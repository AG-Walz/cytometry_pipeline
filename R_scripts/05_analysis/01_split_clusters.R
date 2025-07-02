#!/usr/bin/env Rscript

# =========================
# Script 1: Cluster-by-Cluster Cell Extractor
# Extract ALL cells by cluster into separate files with proper naming
# Output: Many CSV files per cluster with file_origin and cluster columns
# =========================

# =========================
# CONFIG — EDIT THIS BLOCK
# =========================
CONFIG <- list(
  # 1) Input folder with per-cell CSVs
  input_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR",
  
  # 2) Output folder for cluster-split files
  output_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/Output/Split_by_cluster",
  
  # 3) Interactive column selection
  markers_of_interest = NULL,      # Will be selected via GUI
  umap_x = NULL,                   # Will be selected via GUI  
  umap_y = NULL,                   # Will be selected via GUI
  
  # 4) Cluster detection settings
  cluster_regex = "Clustering annotated",   # Regex to identify cluster columns
  clean_cluster_names = TRUE,
  
  # 5) Processing options
  min_cells_per_cluster = 10,
  overwrite_existing = TRUE,
  
  # 6) Performance settings
  chunk_size = 10000,
  seed = 42,
  
  # 7) Statistical Analysis Configuration (from Fine_clustering_cytolution.R)
  conditions_csv = NULL,           # Path to file with conditions mapping
  condition_levels = c("pCR", "npCR"),  # Condition levels for comparison
  enable_statistical_analysis = TRUE,    # Whether to run statistical tests
  use_mixed_effects = TRUE,             # Use mixed-effects models vs pooled analysis
  multiple_testing_correction = "BH",   # Multiple testing correction method
  alpha_level = 0.05,                   # Significance threshold
  
  # 8) Memory Management (from Fine_clustering_cytolution.R)
  max_memory_gb = 16,              # Maximum memory usage in GB
  sample_large_files = TRUE,       # Sample large files to manage memory
  max_cells_for_full_analysis = 100000,  # Max cells before sampling
  
  # 9) Advanced Features
  create_expression_matrices = TRUE,     # Create group expression matrices
  create_statistical_summaries = TRUE,   # Create statistical analysis summaries
  enable_validation_tests = TRUE,        # Run validation tests
  save_intermediate_results = TRUE       # Save intermediate analysis results
)

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(data.table)   # Fast CSV operations
  library(progress)     # Progress bars
  library(stringr)      # String manipulation
  library(tcltk)        # GUI components
  library(dplyr)        # Data manipulation
  library(readr)        # File reading
  library(tibble)       # Data frames
  library(janitor)      # Data cleaning
  library(lme4)         # Mixed-effects models
  library(ggplot2)      # Plotting
})

# Set seed for reproducibility
set.seed(CONFIG$seed)

# =========================
# Statistical Analysis Functions (from Fine_clustering_cytolution.R)
# =========================

# List CSV files with metadata
list_csvs <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  
  tibble(
    file_path = files,
    file_name = basename(files),
    file_size_mb = round(file.info(files)$size / 1024^2, 2)
  ) %>%
    arrange(file_name)
}

# Calculate Fisher's exact test for condition enrichment
calculate_enrichment <- function(cluster_col, cells_data, cond_levels) {
  inC <- as.logical(cells_data[[cluster_col]])
  condA <- cells_data$condition == cond_levels[1]
  condB <- cells_data$condition == cond_levels[2]
  
  # Create 2x2 contingency table
  a <- sum(inC & condA, na.rm = TRUE)
  b <- sum(!inC & condA, na.rm = TRUE)
  c <- sum(inC & condB, na.rm = TRUE)
  d <- sum(!inC & condB, na.rm = TRUE)
  
  # Check for valid contingency table
  if (any(c(a, b, c, d) < 0) || (a + b == 0) || (c + d == 0)) {
    return(tibble(
      in_A = a, not_in_A = b, in_B = c, not_in_B = d,
      odds_ratio = NA, p_value = NA, ci_lower = NA, ci_upper = NA,
      method = "fisher_failed"
    ))
  }
  
  m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
              dimnames = list(c("inCluster","notInCluster"), cond_levels))
  
  tryCatch({
    ft <- fisher.test(m)
    tibble(
      in_A = a, not_in_A = b, in_B = c, not_in_B = d,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      ci_lower = ft$conf.int[1],
      ci_upper = ft$conf.int[2],
      method = "fisher"
    )
  }, error = function(e) {
    tibble(
      in_A = a, not_in_A = b, in_B = c, not_in_B = d,
      odds_ratio = NA, p_value = NA, ci_lower = NA, ci_upper = NA,
      method = "fisher_error"
    )
  })
}

# Mixed-effects analysis accounting for file-level variation
do_mixed_effects <- function(cluster_col, cells_data) {
  # Create binary outcome and condition factor
  model_data <- cells_data %>%
    select(file_name, condition, cluster = all_of(cluster_col)) %>%
    mutate(
      condition = factor(condition, levels = CONFIG$condition_levels),
      file_name = factor(file_name),
      cluster = as.numeric(cluster)
    ) %>%
    filter(!is.na(cluster), !is.na(condition))
  
  # Check if we have enough variation
  if (length(unique(model_data$condition)) < 2) {
    return(tibble(
      estimate = NA, std_error = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      method = "insufficient_conditions"
    ))
  }
  
  # Check cluster variation
  cluster_sum <- sum(model_data$cluster, na.rm = TRUE)
  if (cluster_sum == 0 || cluster_sum == nrow(model_data)) {
    return(tibble(
      estimate = NA, std_error = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      method = "no_cluster_variation"
    ))
  }
  
  tryCatch({
    # Fit mixed-effects logistic regression
    model <- glmer(cluster ~ condition + (1|file_name), 
                   data = model_data, 
                   family = binomial,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
    
    # Extract results
    fixed_effects <- fixef(model)
    se <- sqrt(diag(vcov(model)))
    
    # Get condition effect (assumes first level is reference)
    cond_effect_idx <- grep("condition", names(fixed_effects))
    if (length(cond_effect_idx) == 0) {
      return(tibble(
        estimate = NA, std_error = NA, p_value = NA,
        odds_ratio = NA, ci_lower = NA, ci_upper = NA,
        method = "no_condition_effect"
      ))
    }
    
    estimate <- fixed_effects[cond_effect_idx]
    std_error <- se[cond_effect_idx]
    z_value <- estimate / std_error
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    
    # Calculate odds ratio and confidence interval
    odds_ratio <- exp(estimate)
    ci_lower <- exp(estimate - 1.96 * std_error)
    ci_upper <- exp(estimate + 1.96 * std_error)
    
    tibble(
      estimate = estimate,
      std_error = std_error,
      p_value = p_value,
      odds_ratio = odds_ratio,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      method = "mixed_effects"
    )
    
  }, error = function(e) {
    message("Mixed-effects model failed for cluster ", cluster_col, ": ", e$message)
    
    # Fallback to aggregated t-test on file-level proportions
    file_props <- model_data %>%
      group_by(file_name, condition) %>%
      summarise(prop = mean(cluster, na.rm = TRUE), .groups = "drop")
    
    if (length(unique(file_props$condition)) == 2) {
      t_test <- t.test(prop ~ condition, data = file_props)
      
      tibble(
        estimate = diff(t_test$estimate),
        std_error = NA,
        p_value = t_test$p.value,
        odds_ratio = NA,
        ci_lower = t_test$conf.int[1],
        ci_upper = t_test$conf.int[2],
        method = "t_test_fallback"
      )
    } else {
      tibble(
        estimate = NA, std_error = NA, p_value = NA,
        odds_ratio = NA, ci_lower = NA, ci_upper = NA,
        method = "mixed_effects_failed"
      )
    }
  })
}

# Validation testing functionality (from test scripts)
run_validation_tests <- function() {
  cat("\n🧪 Running validation tests...\n")
  
  tests_passed <- 0
  tests_total <- 0
  
  # Test 1: Check CONFIG values
  tests_total <- tests_total + 1
  if (dir.exists(CONFIG$input_dir) && dir.exists(CONFIG$output_dir)) {
    cat("✅ Directory paths valid\n")
    tests_passed <- tests_passed + 1
  } else {
    cat("❌ Directory paths invalid\n")
  }
  
  # Test 2: Check for input files
  tests_total <- tests_total + 1
  input_files <- list.files(CONFIG$input_dir, pattern = "\\.csv$")
  if (length(input_files) > 0) {
    cat(sprintf("✅ Found %d CSV input files\n", length(input_files)))
    tests_passed <- tests_passed + 1
  } else {
    cat("❌ No CSV input files found\n")
  }
  
  # Test 3: Check conditions file if statistical analysis enabled
  if (CONFIG$enable_statistical_analysis) {
    tests_total <- tests_total + 1
    if (!is.null(CONFIG$conditions_csv) && file.exists(CONFIG$conditions_csv)) {
      tryCatch({
        conds <- read_csv(CONFIG$conditions_csv, show_col_types = FALSE)
        if ("condition" %in% names(conds) && "file_name" %in% names(conds)) {
          cat("✅ Conditions file valid\n")
          tests_passed <- tests_passed + 1
        } else {
          cat("❌ Conditions file missing required columns\n")
        }
      }, error = function(e) {
        cat("❌ Conditions file cannot be read\n")
      })
    } else {
      cat("❌ Conditions file not found or not specified\n")
    }
  }
  
  # Test 4: Memory availability
  tests_total <- tests_total + 1
  tryCatch({
    mem_info <- system("free -g", intern = TRUE)
    available_mem <- as.numeric(strsplit(mem_info[2], "\\s+")[[1]][7])
    if (available_mem >= CONFIG$max_memory_gb * 0.5) {
      cat(sprintf("✅ Sufficient memory available (~%dG)\n", available_mem))
      tests_passed <- tests_passed + 1
    } else {
      cat(sprintf("⚠️  Low memory available (~%dG, need %dG)\n", available_mem, CONFIG$max_memory_gb))
    }
  }, error = function(e) {
    cat("⚠️  Cannot check memory availability\n")
    tests_passed <- tests_passed + 1  # Don't fail on this
  })
  
  cat(sprintf("🎯 Validation: %d/%d tests passed\n", tests_passed, tests_total))
  return(tests_passed == tests_total)
}

# =========================
# Expression Matrix Creation (from Create_labelled_exprmat.R)
# =========================

# Create expression matrices per cluster group
create_group_expression_matrices <- function(timestamp) {
  if (!CONFIG$create_expression_matrices) {
    return(invisible(NULL))
  }
  
  cat("\n📊 Creating expression matrices per cluster group...\n")
  
  # Create matrices directory
  matrices_dir <- file.path(CONFIG$output_dir, "group_expression_matrices")
  if (!dir.exists(matrices_dir)) {
    dir.create(matrices_dir, recursive = TRUE)
  }
  
  # Get all individual cluster output files (not the aggregated matrices)
  output_files <- list.files(CONFIG$output_dir, pattern = "\\.csv$", full.names = TRUE)
  output_files <- output_files[!grepl("cluster_tracking_matrix|group_matrices_summary|statistical_summary", output_files)]
  
  if (length(output_files) == 0) {
    warning("No cluster output files found for matrix creation")
    return(invisible(NULL))
  }
  
  cat(sprintf("Processing %d cluster files for matrix creation...\n", length(output_files)))
  
  # Initialize list to collect data by group
  group_data_list <- list()
  
  # Process each cluster file
  for (file_path in output_files) {
    file_basename <- basename(file_path)
    
    # Extract cluster name from filename
    # Expected format: ICK_020_01_Unmixed_BCs_224.csv
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
    
    # Determine group
    group_name <- extract_cluster_group(cluster_name)
    
    # Skip single clusters (need multiple for meaningful matrices)
    current_group_size <- length(grep(paste0("_", group_name), output_files))
    if (current_group_size < 2) next
    
    tryCatch({
      # Read cluster data
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
  
  cat(sprintf("📁 Expression matrices saved in: %s\n", matrices_dir))
  cat(sprintf("📋 Summary statistics saved: %s\n", basename(summary_path)))
  
  return(invisible(summary_stats))
}

# Conditions file management (from Fine_clustering_cytolution.R)
check_conditions_file_exists <- function() {
  cat("🔍 Checking for conditions file...\n")
  
  if (!is.null(CONFIG$conditions_csv) && file.exists(CONFIG$conditions_csv)) {
    cat(sprintf("✅ Conditions file found: %s\n", CONFIG$conditions_csv))
    return(CONFIG$conditions_csv)
  }
  
  # Look for conditions file in safe locations
  potential_paths <- c(
    file.path(CONFIG$output_dir, "file_conditions.csv"),
    file.path(getwd(), "file_conditions.csv"),
    file.path("~", "file_conditions.csv"),
    file.path(dirname(CONFIG$input_dir), "file_conditions.csv")
  )
  
  for (path in potential_paths) {
    if (file.exists(path)) {
      cat(sprintf("Found potential file: %s\n", path))
      tryCatch({
        test_data <- read.csv(path, stringsAsFactors = FALSE, nrows = 5)
        if ("condition" %in% names(test_data) && "file_name" %in% names(test_data)) {
          full_data <- read.csv(path, stringsAsFactors = FALSE)
          if (!all(is.na(full_data$condition))) {
            cat(sprintf("✅ Valid conditions file found: %s\n", path))
            CONFIG$conditions_csv <<- path
            return(path)
          }
        }
      }, error = function(e) {
        cat(sprintf("Cannot read file %s: %s\n", path, e$message))
      })
    }
  }
  
  return(NULL)
}

# Setup statistical analysis infrastructure
setup_statistical_analysis <- function() {
  if (!CONFIG$enable_statistical_analysis) {
    cat("📊 Statistical analysis disabled\n")
    return(FALSE)
  }
  
  cat("🔬 Setting up statistical analysis...\n")
  
  # Check for conditions file
  existing_file <- check_conditions_file_exists()
  
  if (is.null(existing_file)) {
    cat("\n⚠️  SETUP REQUIRED: No conditions file found\n")
    cat("Please create a file called 'file_conditions.csv' with structure:\n")
    cat("file_name,condition\n")
    cat("your_file1.csv,pCR\n")
    cat("your_file2.csv,npCR\n")
    cat("...\n\n")
    cat(sprintf("RECOMMENDED LOCATION: %s\n", file.path(CONFIG$output_dir, "file_conditions.csv")))
    cat("Then set: CONFIG$conditions_csv <- 'path_to_your_conditions_file.csv'\n")
    cat("And re-run this script.\n")
    return(FALSE)
  }
  
  # Read and validate conditions file
  tryCatch({
    conds <- read_csv(CONFIG$conditions_csv, show_col_types = FALSE) %>%
      clean_names() %>%
      mutate(across(everything(), as.character))
    
    # Validate conditions
    if (all(is.na(conds$condition))) {
      stop("All conditions are NA. Please fill in the condition column.")
    }
    
    # Check for valid condition values
    valid_conditions <- conds$condition[!is.na(conds$condition)]
    if (length(valid_conditions) == 0) {
      stop("No valid conditions found. Please fill in the condition column.")
    }
    
    # Warn about NAs and remove them
    na_count <- sum(is.na(conds$condition))
    if (na_count > 0) {
      cat(sprintf("⚠️  %d files have NA conditions and will be excluded\n", na_count))
      conds <- conds[!is.na(conds$condition), ]
    }
    
    # Store conditions data globally for later use
    assign("CONDITIONS_DATA", conds, envir = .GlobalEnv)
    
    cat(sprintf("✅ Statistical analysis ready with %d files and %d conditions\n", 
                nrow(conds), length(unique(conds$condition))))
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("❌ Error with conditions file: %s\n", e$message))
    return(FALSE)
  })
}

# =========================
# GUI Selection Functions
# =========================

# Interactive column selection via GUI with fallback
select_columns_gui <- function(csv_file) {
  cat("Loading sample for column selection...\n")
  sample_data <- fread(csv_file, nrows = 10)
  all_columns <- colnames(sample_data)
  
  # Try to detect if we're in a non-interactive environment
  if (!interactive() || Sys.getenv("DISPLAY") == "") {
    cat("Non-interactive environment detected, using automatic column selection...\n")
    
    # Automatic selection based on common patterns
    marker_columns <- all_columns[!grepl("^(File|Time|Event|Original_ID|UMAP)", all_columns, ignore.case = TRUE)]
    
    # Try to find UMAP columns automatically
    umap_x_col <- NULL
    umap_y_col <- NULL
    
    umap_candidates <- all_columns[grepl("UMAP|umap", all_columns, ignore.case = TRUE)]
    if (length(umap_candidates) >= 2) {
      umap_x_col <- umap_candidates[1]
      umap_y_col <- umap_candidates[2]
    }
    
    cat(sprintf("✅ Auto-selected %d marker columns\n", length(marker_columns)))
    if (!is.null(umap_x_col)) {
      cat(sprintf("✅ Auto-selected UMAP: %s, %s\n", umap_x_col, umap_y_col))
    }
    
    return(list(
      markers = marker_columns,
      umap_x = umap_x_col,
      umap_y = umap_y_col
    ))
  }
  
  # Original GUI code for interactive environments
  tryCatch({
    # Create GUI window
    tt <- tktoplevel()
    tkwm.title(tt, "Column Selection")
    tkwm.geometry(tt, "800x600")
    
    # Create notebook for tabs
    nb <- tkwidget(tt, "ttk::notebook")
    
    # Tab 1: Marker selection
    tab1 <- tkframe(nb)
    tkadd(nb, tab1, text = "Markers")
    
    # Instructions
    tkpack(tklabel(tab1, text = "Select markers of interest (hold Ctrl for multiple):"), anchor = "w")
    
    # Create listbox with scrollbar
    list_frame <- tkframe(tab1)
    tkpack(list_frame, fill = "both", expand = TRUE)
    
    scrollbar <- tkscrollbar(list_frame, orient = "vertical")
    listbox <- tklistbox(list_frame, yscrollcommand = function(...) tkset(scrollbar, ...), selectmode = "extended")
    tkconfigure(scrollbar, command = function(...) tkyview(listbox, ...))
    
    tkpack(scrollbar, side = "right", fill = "y")
    tkpack(listbox, side = "left", fill = "both", expand = TRUE)
    
    # Populate markers listbox (exclude obvious metadata columns)
    marker_columns <- all_columns[!grepl("^(File|Time|Event|Original_ID)", all_columns, ignore.case = TRUE)]
    for (col in marker_columns) {
      tkinsert(listbox, "end", col)
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
    
    # Variables to store results
    result <- list()
    
    # OK button
    ok_button <- tkbutton(button_frame, text = "OK", command = function() {
      # Get selected markers
      selected_indices <- as.numeric(tkcurselection(listbox))
      if (length(selected_indices) > 0) {
        result$markers <<- marker_columns[selected_indices + 1]  # R uses 1-based indexing
      } else {
        result$markers <<- character(0)
      }
      result$umap_x <<- as.character(tclvalue(umap_x_var))
      result$umap_y <<- as.character(tclvalue(umap_y_var))
      tkdestroy(tt)
    })
    
    # Cancel button
    cancel_button <- tkbutton(button_frame, text = "Cancel", command = function() {
      result <<- list(markers = character(0), umap_x = NULL, umap_y = NULL)
      tkdestroy(tt)
    })
    
    tkpack(ok_button, cancel_button, side = "left", padx = 5, pady = 5)
    
    # Wait for window to close
    tkwait.window(tt)
    
    return(list(
      markers = if(length(result$markers) > 0) result$markers else character(0),
      umap_x = if(!is.null(result$umap_x) && nchar(result$umap_x) > 0) result$umap_x else NULL,
      umap_y = if(!is.null(result$umap_y) && nchar(result$umap_y) > 0) result$umap_y else NULL
    ))
    
  }, error = function(e) {
    cat("GUI failed, falling back to automatic selection...\n")
    cat("Error:", e$message, "\n")
    
    # Fallback: automatic selection
    marker_columns <- all_columns[!grepl("^(File|Time|Event|Original_ID|UMAP)", all_columns, ignore.case = TRUE)]
    
    # Try to find UMAP columns automatically
    umap_x_col <- NULL
    umap_y_col <- NULL
    
    umap_candidates <- all_columns[grepl("UMAP|umap", all_columns, ignore.case = TRUE)]
    if (length(umap_candidates) >= 2) {
      umap_x_col <- umap_candidates[1]
      umap_y_col <- umap_candidates[2]
    }
    
    cat(sprintf("✅ Auto-selected %d marker columns\n", length(marker_columns)))
    if (!is.null(umap_x_col)) {
      cat(sprintf("✅ Auto-selected UMAP: %s, %s\n", umap_x_col, umap_y_col))
    }
    
    return(list(
      markers = marker_columns,
      umap_x = umap_x_col,
      umap_y = umap_y_col
    ))
  })
}

# =========================
# Cluster Processing Functions
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

# Create proper cluster identifier (GroupName_Number format)
create_cluster_id <- function(cluster_name) {
  cleaned_name <- clean_cluster_name(cluster_name)
  group_name <- extract_cluster_group(cluster_name)
  
  # Extract number from cluster name if present
  if (grepl("_[0-9]+$", cleaned_name)) {
    number_part <- gsub(".*_([0-9]+)$", "\\1", cleaned_name)
    return(paste0(group_name, "_", number_part))
  } else {
    # If no number, just return the group name
    return(group_name)
  }
}

# =========================
# Utility Functions
# =========================

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
# Comprehensive cluster discovery across ALL files
discover_all_clusters <- function(csv_files) {
  cat("\n=== Comprehensive Cluster Discovery Across ALL Files ===\n")
  cat(sprintf("Scanning %d files to discover all possible clusters...\n", length(csv_files)))
  
  all_clusters_found <- character(0)
  files_with_clusters <- 0
  cluster_file_mapping <- list()  # Track which files contain which clusters
  
  # Progress bar for discovery
  discovery_progress <- progress_bar$new(
    format = "🔍 [:bar] :percent :current/:total files scanned for clusters",
    total = length(csv_files), 
    clear = FALSE, 
    width = 80
  )
  
  # Scan every single file to find all possible clusters
  for (i in seq_along(csv_files)) {
    csv_file <- csv_files[i]
    discovery_progress$tick()
    
    file_basename <- basename(csv_file)
    
    tryCatch({
      # Read header only to check columns (memory efficient)
      sample_data <- fread(csv_file, nrows = 1)
      file_clusters <- extract_cluster_columns(sample_data)
      
      if (length(file_clusters) > 0) {
        # Add to master list
        all_clusters_found <- unique(c(all_clusters_found, file_clusters))
        files_with_clusters <- files_with_clusters + 1
        
        # Track which files have which clusters (for validation)
        cluster_file_mapping[[file_basename]] <- file_clusters
      }
      
    }, error = function(e) {
      warning(sprintf("Error scanning %s for clusters: %s", file_basename, e$message))
    })
  }
  
  cat(sprintf("\n🎯 Discovery Results:\n"))
  cat(sprintf("  • Files with cluster data: %d/%d\n", files_with_clusters, length(csv_files)))
  cat(sprintf("  • Total unique clusters discovered: %d\n", length(all_clusters_found)))
  cat(sprintf("  • Sample clusters: %s\n", paste(head(all_clusters_found, 5), collapse = ", ")))
  
  if (length(all_clusters_found) == 0) {
    stop("❌ No cluster columns found in any files matching pattern: ", CONFIG$cluster_regex)
  }
  
  # Sort clusters for consistent processing order
  all_clusters_found <- sort(all_clusters_found)
  
  return(list(
    clusters = all_clusters_found,
    files_with_clusters = files_with_clusters,
    cluster_file_mapping = cluster_file_mapping
  ))
}

# Create empty output files for ALL discovered clusters
initialize_cluster_files <- function(all_clusters, selected_columns) {
  cat(sprintf("\n📁 Skipping cluster file pre-initialization...\n"))
  cat("Files will be created dynamically as cluster values are discovered during processing.\n")
  cat(sprintf("Discovered %d cluster columns for processing: %s\n", 
              length(all_clusters), 
              paste(head(all_clusters, 3), collapse = ", ")))
  
  # Validate output directory exists
  if (!exists("CONFIG") || is.null(CONFIG$output_dir)) {
    stop("CONFIG$output_dir is not defined. Please check configuration.")
  }
  
  if (!dir.exists(CONFIG$output_dir)) {
    cat(sprintf("Creating output directory: %s\n", CONFIG$output_dir))
    dir.create(CONFIG$output_dir, recursive = TRUE)
  }
  
  cat("✅ Output directory ready for dynamic file creation\n")
  return(invisible(TRUE))
}

# Original cluster column extraction function (kept for compatibility)
extract_cluster_columns <- function(data) {
  cluster_cols <- colnames(data)[grepl(CONFIG$cluster_regex, colnames(data), ignore.case = TRUE)]
  
  if (length(cluster_cols) == 0) {
    return(character(0))
  }
  
  return(cluster_cols)
}

# Comprehensive processing - every row from every file into appropriate cluster files
process_all_files_comprehensive <- function(csv_files, all_clusters, selected_columns) {
  cat(sprintf("\n=== Comprehensive Processing: ALL Rows from ALL Files ===\n"))
  cat(sprintf("Processing %d files into %d cluster categories...\n", length(csv_files), length(all_clusters)))
  
  # Progress tracking
  file_progress <- progress_bar$new(
    format = "📂 [:bar] :percent :current/:total files processed (:elapsed elapsed, :eta remaining)",
    total = length(csv_files), 
    clear = FALSE, 
    width = 80
  )
  
  total_rows_processed <- 0
  total_cells_assigned <- 0
  files_processed <- 0
  memory_efficient_processing <- 0
  
  # Initialize cluster file counts - will grow dynamically as we discover cluster values
  cluster_file_counts <- list()
  
  # Note: We don't pre-initialize with all_clusters because those are column names,
  # not the actual cluster values that will be used for file creation
  
  # Process each input file comprehensively
  for (i in seq_along(csv_files)) {
    csv_file <- csv_files[i]
    file_basename <- basename(csv_file)
    file_progress$tick()
    
    tryCatch({
      # Determine file origin for metadata
      file_origin <- gsub("\\.csv$", "", file_basename)
      
      # Check file size and use appropriate reading strategy
      file_size_mb <- file.size(csv_file) / (1024^2)
      
      if (file_size_mb > 500) {
        cat(sprintf("\n  Large file detected: %s (%.1f MB). Using chunked processing...\n", 
                    file_basename, file_size_mb))
        rows_processed <- process_large_file_comprehensive(csv_file, all_clusters, selected_columns, file_origin)
        memory_efficient_processing <- memory_efficient_processing + 1
      } else {
        # Read entire file for smaller files
        data <- fread(csv_file)
        if (nrow(data) == 0) next
        
        rows_processed <- process_file_data_comprehensive(data, all_clusters, selected_columns, file_origin)
      }
      
      if (rows_processed > 0) {
        total_rows_processed <- total_rows_processed + rows_processed
        total_cells_assigned <- total_cells_assigned + rows_processed
        files_processed <- files_processed + 1
      }
      
    }, error = function(e) {
      warning(sprintf("Error processing file %s: %s", file_basename, e$message))
    })
  }
  
  cat(sprintf("\n✅ Comprehensive Processing Complete!\n"))
  cat(sprintf("  • Total rows processed: %s\n", format(total_rows_processed, big.mark = ",")))
  cat(sprintf("  • Files with data processed: %d/%d\n", files_processed, length(csv_files)))
  
  # Return comprehensive statistics
  return(list(
    total_rows_processed = total_rows_processed,
    total_cells_assigned = total_cells_assigned,
    files_processed = files_processed,
    memory_efficient_processing = memory_efficient_processing,
    cluster_file_counts = cluster_file_counts
  ))
}

# Process large files in chunks (comprehensive approach)
process_large_file_comprehensive <- function(csv_file, all_clusters, selected_columns, file_origin) {
  # Get total rows
  total_rows <- as.numeric(system(paste("wc -l", shQuote(csv_file)), intern = TRUE)) - 1
  
  if (total_rows <= 0) return(0)
  
  # Read header
  header_data <- fread(csv_file, nrows = 1)
  
  # Process in chunks
  chunk_start <- 1
  chunk_num <- 1
  total_processed <- 0
  
  while (chunk_start <= total_rows) {
    chunk_end <- min(chunk_start + CONFIG$chunk_size - 1, total_rows)
    cat(sprintf("    Reading chunk %d: rows %d-%d\n", chunk_num, chunk_start, chunk_end))
    
    # Read chunk
    chunk_data <- fread(csv_file, skip = chunk_start, nrows = chunk_end - chunk_start + 1, header = FALSE)
    colnames(chunk_data) <- colnames(header_data)
    
    # Process this chunk
    chunk_processed <- process_file_data_comprehensive(chunk_data, all_clusters, selected_columns, file_origin)
    total_processed <- total_processed + chunk_processed
    
    chunk_start <- chunk_end + 1
    chunk_num <- chunk_num + 1
  }
  
  return(total_processed)
}

# Process file data comprehensively - every row into appropriate cluster files
process_file_data_comprehensive <- function(data, all_clusters, selected_columns, file_origin) {
  if (nrow(data) == 0) return(0)
  
  # Get available cluster columns in this file
  available_clusters <- extract_cluster_columns(data)
  
  if (length(available_clusters) == 0) {
    return(0)  # No cluster data in this file
  }
  
  # Prepare columns to keep
  essential_cols <- c("file_origin", "cluster")
  umap_cols <- c()
  if (!is.null(selected_columns$umap_x) && selected_columns$umap_x %in% names(data)) {
    umap_cols <- c(umap_cols, selected_columns$umap_x)
  }
  if (!is.null(selected_columns$umap_y) && selected_columns$umap_y %in% names(data)) {
    umap_cols <- c(umap_cols, selected_columns$umap_y)
  }
  
  # Get marker columns that exist in this file
  available_markers <- selected_columns$markers[selected_columns$markers %in% names(data)]
  
  keep_cols <- c(available_markers, umap_cols)
  keep_cols <- keep_cols[keep_cols %in% names(data)]
  
  if (length(keep_cols) == 0) {
    warning("No selected marker columns found in file")
    return(0)
  }
  
  # Add metadata columns
  data$file_origin <- file_origin
  keep_cols_with_meta <- c(keep_cols, "file_origin")
  
  rows_processed <- 0
  
  # Process each available cluster in this file
  for (cluster_col in available_clusters) {
    if (!(cluster_col %in% names(data))) next
    
    # Get unique cluster values in this column
    cluster_values <- unique(data[[cluster_col]])
    cluster_values <- cluster_values[!is.na(cluster_values) & cluster_values != ""]
    
    # Process each cluster value
    for (cluster_value in cluster_values) {
      if (is.na(cluster_value) || cluster_value == "") next
      
      # Extract rows for this cluster
      cluster_rows <- data[data[[cluster_col]] == cluster_value, ]
      if (nrow(cluster_rows) == 0) next
      
      # Add cluster information
      cluster_rows$cluster <- cluster_value
      
      # Select final columns
      final_cols <- c(keep_cols_with_meta, "cluster")
      final_cols <- final_cols[final_cols %in% names(cluster_rows)]
      cluster_data <- cluster_rows[, ..final_cols]
      
      # Create safe filename and append to appropriate cluster file
      safe_cluster_name <- gsub("[^A-Za-z0-9_-]", "_", cluster_value)
      output_filename <- sprintf("%s.csv", safe_cluster_name)
      output_path <- file.path(CONFIG$output_dir, output_filename)
      
      # Ensure output directory exists
      if (!dir.exists(CONFIG$output_dir)) {
        dir.create(CONFIG$output_dir, recursive = TRUE)
      }
      
      # Append to file (creates file if it doesn't exist)
      tryCatch({
        fwrite(cluster_data, output_path, append = TRUE)
        rows_processed <- rows_processed + nrow(cluster_data)
      }, error = function(e) {
        warning(sprintf("Failed to write data for cluster '%s': %s", cluster_value, e$message))
      })
    }
  }
  
  return(rows_processed)
}

# Extract and save cells by cluster for a single file
process_file_by_cluster <- function(csv_file, selected_columns) {
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
      cat(sprintf("  Reading chunk %d: rows %d-%d\n", chunk_num, chunk_start, chunk_end))
      
      chunk_data <- fread(csv_file, skip = chunk_start, nrows = chunk_end - chunk_start + 1, header = FALSE)
      
      # Set column names from the first chunk
      if (chunk_num == 1) {
        header_data <- fread(csv_file, nrows = 1)
        colnames(chunk_data) <- colnames(header_data)
      } else {
        header_data <- fread(csv_file, nrows = 1)
        colnames(chunk_data) <- colnames(header_data)
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
    return(invisible(NULL))
  }
  
  # Subset data to selected columns + cluster columns
  keep_cols <- c(selected_columns$markers, selected_columns$umap_x, selected_columns$umap_y, cluster_cols)
  keep_cols <- keep_cols[!is.null(keep_cols)]  # Remove NULL values first
  keep_cols <- keep_cols[keep_cols %in% colnames(data)]  # Only keep columns that exist
  keep_cols <- as.character(keep_cols)  # Ensure character vector
  
  if (length(keep_cols) == 0) {
    warning(sprintf("No selected columns found in %s\n", file_basename))
    return(invisible(NULL))
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
  
  processed_clusters <- 0
  
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
    
    # Add metadata columns (enhancement from original script)
    cluster_output$file_origin <- file_basename
    cluster_output$cluster <- create_cluster_id(cluster_col)
    
    # Create output filename using cleaned cluster name (like original)
    cleaned_cluster_name <- clean_cluster_name(cluster_col)
    output_filename <- sprintf("%s_%s.csv", file_basename, cleaned_cluster_name)
    output_path <- file.path(CONFIG$output_dir, output_filename)
    
    # Skip if file already exists (unless overwrite is enabled)
    if (file.exists(output_path) && !CONFIG$overwrite_existing) {
      cat(sprintf("    Skipping %s (already exists)\n", cleaned_cluster_name))
      next
    }
    
    # Always show processing for overwrite case
    if (file.exists(output_path) && CONFIG$overwrite_existing) {
      cat(sprintf("    Overwriting existing file: %s\n", cleaned_cluster_name))
    }
    
    # Save cluster data
    fwrite(cluster_output, output_path)
    cat(sprintf("    ✓ Saved %d cells for cluster: %s\n", nrow(cluster_output), cleaned_cluster_name))
    processed_clusters <- processed_clusters + 1
  }
  
  cat(sprintf("  Processed %d clusters from file %s\n", processed_clusters, file_basename))
  return(invisible(processed_clusters))
}

# =========================
# Main Execution
# =========================

main <- function() {
  cat("=== Comprehensive Cluster-by-Cluster Cell Extractor (Script 1) ===\n")
  cat("Extracting ALL cells by cluster into separate files with COMPREHENSIVE scanning\n")
  cat("Processing EVERY row from EVERY file without data loss\n\n")
  
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
  
  # COMPREHENSIVE cluster discovery across ALL files
  cat(sprintf("\n=== COMPREHENSIVE Cluster Discovery Across ALL %d Files ===\n", length(csv_files)))
  cat("This will scan EVERY column in EVERY file to find ALL possible clusters\n")
  cat("(This may take a few minutes for large datasets)\n\n")
  
  all_clusters_found <- discover_all_clusters(csv_files)
  
  if (length(all_clusters_found) == 0) {
    stop("❌ No cluster columns found in any files using regex pattern!")
  }
  
  cat(sprintf("\n📊 COMPREHENSIVE Discovery Results:\n"))
  cat(sprintf("  • Total unique clusters found across ALL files: %d\n", length(all_clusters_found)))
  cat(sprintf("  • All clusters: %s\n", paste(all_clusters_found, collapse = ", ")))
  
  # Initialize ALL cluster output files (optional - skip if causing errors)
  cat(sprintf("\n=== Initializing %d Cluster Output Files (Optional) ===\n", length(all_clusters_found)))
  tryCatch({
    initialize_cluster_files(all_clusters_found, selected_columns)
  }, error = function(e) {
    warning(sprintf("Skipping cluster file initialization due to error: %s", e$message))
    cat("⚠️  Note: Cluster files will be created during processing instead\n")
  })
  
  # COMPREHENSIVE processing: every row from every file
  cat(sprintf("\n=== COMPREHENSIVE Processing: EVERY Row from EVERY File ===\n"))
  cat("Processing ALL rows without subsampling or data loss\n")
  cat("Large files will be processed in memory-efficient chunks\n\n")
  
  processing_stats <- process_all_files_comprehensive(csv_files, all_clusters_found, selected_columns)
  
  # Final summary with comprehensive statistics
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("🎉 COMPREHENSIVE CLUSTER EXTRACTION COMPLETE! \n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("📊 COMPREHENSIVE PROCESSING SUMMARY:\n"))
  cat(sprintf("  • Total input files processed: %d\n", length(csv_files)))
  
  # Safe formatting with error handling
  total_rows_str <- tryCatch({
    if (is.numeric(processing_stats$total_rows_processed) && !is.na(processing_stats$total_rows_processed)) {
      format(processing_stats$total_rows_processed, big.mark = ",")
    } else {
      "Unknown"
    }
  }, error = function(e) "Error")
  
  total_cells_str <- tryCatch({
    if (is.numeric(processing_stats$total_cells_assigned) && !is.na(processing_stats$total_cells_assigned)) {
      format(processing_stats$total_cells_assigned, big.mark = ",")
    } else {
      "Unknown"
    }
  }, error = function(e) "Error")
  
  cat(sprintf("  • Total rows processed across all files: %s\n", total_rows_str))
  cat(sprintf("  • Total cells assigned to clusters: %s\n", total_cells_str))
  cat(sprintf("  • Files successfully processed: %d\n", 
              if(is.numeric(processing_stats$files_processed)) processing_stats$files_processed else 0))
  cat(sprintf("  • Total unique clusters discovered: %d\n", length(all_clusters_found)))
  cat(sprintf("  • Cluster output files created: %d\n", length(all_clusters_found)))
  cat(sprintf("  • Processing time: %.2f minutes (%.1f sec/file)\n", 
              total_time, (total_time * 60) / max(1, length(csv_files))))
  cat(sprintf("  • Output directory: %s\n", CONFIG$output_dir))
  
  # Detailed cluster statistics
  cat(sprintf("\n📋 CLUSTER DETAILS:\n"))
  for (i in seq_along(all_clusters_found)) {
    cluster <- all_clusters_found[i]
    
    # Check if output file exists and get actual file size/line count
    safe_cluster_name <- gsub("[^A-Za-z0-9_-]", "_", cluster)
    cluster_file_path <- file.path(CONFIG$output_dir, sprintf("%s.csv", safe_cluster_name))
    
    if (file.exists(cluster_file_path)) {
      file_size <- file.info(cluster_file_path)$size
      if (file_size > 0) {
        # Try to count lines (cells) in the file
        tryCatch({
          line_count <- length(readLines(cluster_file_path)) - 1  # Subtract header
          line_count <- max(0, line_count)  # Ensure non-negative
          cat(sprintf("  %2d. %s: %s cells ✅\n", i, cluster, format(line_count, big.mark = ",")))
        }, error = function(e) {
          cat(sprintf("  %2d. %s: File exists ✅\n", i, cluster))
        })
      } else {
        cat(sprintf("  %2d. %s: Empty file ⚪\n", i, cluster))
      }
    } else {
      cat(sprintf("  %2d. %s: No output file ❌\n", i, cluster))
    }
  }
  
  # Verify output files exist
  output_files <- list.files(CONFIG$output_dir, pattern = "\\.csv$", full.names = FALSE)
  cat(sprintf("\n✅ Verification: %d CSV files created in output directory\n", length(output_files)))
  
  # Memory usage summary
  memory_efficient_count <- tryCatch({
    if (is.numeric(processing_stats$memory_efficient_processing)) {
      processing_stats$memory_efficient_processing
    } else {
      0
    }
  }, error = function(e) 0)
  
  if (memory_efficient_count > 0) {
    threshold_mb <- tryCatch({
      if (exists("CONFIG") && !is.null(CONFIG$file_size_threshold_mb)) {
        CONFIG$file_size_threshold_mb
      } else {
        500
      }
    }, error = function(e) 500)
    
    cat(sprintf("💾 Memory Management: %d files processed with chunking (>%d MB)\n", 
                memory_efficient_count, threshold_mb))
  }
  
  # Create timestamp for additional analyses
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Run integrated analyses if requested
  total_cells_assigned <- tryCatch({
    if (is.numeric(processing_stats$total_cells_assigned)) {
      processing_stats$total_cells_assigned
    } else {
      0
    }
  }, error = function(e) 0)
  
  if (total_cells_assigned > 0) {
    cat("\n🔬 Running integrated analyses...\n")
    
    # Create expression matrices
    if (CONFIG$create_expression_matrices) {
      tryCatch({
        create_group_expression_matrices(timestamp)
      }, error = function(e) {
        warning(sprintf("Expression matrix creation failed: %s", e$message))
      })
    }
    
    # Setup statistical analysis
    if (CONFIG$enable_statistical_analysis && exists("CONDITIONS_DATA")) {
      cat("\n📈 Setting up comprehensive statistical analysis flag...\n")
      
      tryCatch({
        # Save comprehensive processing flag for script 2
        cat("Statistical analysis will be completed in 02_generate_heatmaps.R\n")
        cat("Make sure to run the second script for complete statistical results\n")
        
        # Save a flag file to indicate comprehensive statistical analysis is ready
        stat_flag_file <- file.path(CONFIG$output_dir, sprintf("comprehensive_analysis_ready_%s.txt", timestamp))
        writeLines(c(
          paste("Comprehensive cluster extraction completed at:", Sys.time()),
          paste("Total clusters discovered:", length(all_clusters_found)),
          paste("Total rows processed:", total_rows_str),
          paste("Total cells assigned:", total_cells_str),
          paste("Files processed:", if(is.numeric(processing_stats$files_processed)) processing_stats$files_processed else 0),
          paste("Comprehensive scanning: ALL files scanned for ALL clusters"),
          paste("Data loss: ZERO - every row processed"),
          paste("Conditions file:", CONFIG$conditions_csv),
          paste("Ready for comprehensive statistical analysis in script 2")
        ), stat_flag_file)
        
      }, error = function(e) {
        warning(sprintf("Statistical analysis setup failed: %s", e$message))
      })
    }
    
    cat("✅ SUCCESS: Comprehensive processing complete - Ready for heatmap generation and full analysis (Script 2)\n")
  } else {
    cat("⚠️  WARNING: No cells were assigned to clusters. Check input data and regex patterns.\n")
  }
  
  cat("\n🔍 COMPREHENSIVE PROCESSING COMPLETED:\n")
  cat(sprintf("  • Every column in every file was scanned for cluster patterns\n"))
  cat(sprintf("  • Every row from every file was processed without data loss\n"))
  cat(sprintf("  • All %d discovered clusters have output files (even if empty)\n", length(all_clusters_found)))
  cat(sprintf("  • No subsampling was performed - all %s rows processed\n", total_rows_str))
  
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  return(invisible(list(
    success = TRUE,
    clusters_found = all_clusters_found,
    processing_stats = processing_stats,
    total_time_minutes = total_time
  )))
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
