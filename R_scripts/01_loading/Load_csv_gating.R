# ============================================================================
# Load CSV Gating Data - Unified Version
# ============================================================================
# This script combines functionality from all versions to provide:
# - Automatic CSV-to-FCS matching
# - Interactive marker/gate selection
# - Memory-efficient data.table processing
# - Flexible separator detection
# - Chunked processing for large files
# - Comprehensive error handling and validation
# ============================================================================
#
# VERSION HISTORY & MERGED FEATURES:
# - Load_csv_gating.R (original): Basic CSV loading with manual label dict
# - Load_csv_gating_2.R: Added filter_valid_csv_files, process_single_csv
# - Load_csv_gating_3.R: Added future.apply parallelization, lme4/emmeans stats
# - This unified version: Combines all features with CAFE format support
#
# TODO (from older versions to consider integrating):
# - [ ] Load_csv_gating_3.R uses future::plan(multisession) for parallelization
#       which may be more stable than foreach/doParallel for some systems
# - [ ] Load_csv_gating_3.R includes lme4, lmerTest, emmeans for statistics
#       Consider adding optional statistical analysis after gating
# - [ ] Add option to output in long format (population column) as in _3.R
# ============================================================================

# Required Libraries: flowWorkspace, flowCore, dplyr, openxlsx, data.table, tcltk
# !For RStudio only: Set working directory to the same path as the current file

# Optional: Set custom library path if needed
# .libPaths("Q:/IM/AG_Märklin/FSFC_Pipeline/library/")

# Auto-install missing packages
pkgs <- c("data.table", "flowCore", "flowWorkspace", "dplyr", "openxlsx", 
          "tcltk", "stringr", "progress", "parallel", "foreach", "doParallel")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) {
  cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE))

# Detect available cores for parallel processing
# Conservative setting: use max 4 cores to prevent memory issues
detected_cores <- parallel::detectCores()
n_cores <- min(4, max(1, detected_cores - 1))  # Max 4 cores, leave one free
cat("Detected", detected_cores, "cores, will use max", n_cores, "for parallel processing\n")
cat("NOTE: Parallel processing disabled for stability - using sequential mode\n")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# =============================== Parameters ===================================

# Data source mode:
# "csv_with_fcs" - CSV files reference separate FCS files (traditional workflow)
# "csv_only" - CSV files contain all data (expression + optional gating)
# "auto" - Automatically detect based on file presence
data_source_mode <- "auto"

# Processing mode:
# "automatic" - Interactive selection of markers/gates
# "manual" - Use predefined label dictionaries
# "expression_only" - No gating data, treat all cells as single population
processing_mode <- "automatic"

# Input paths
input_fcs_path <- "/data/Cytometry/ICK/files_ICK18_80_pre-gating_11_18_2025/"
input_csv_path <- input_fcs_path  # Can be different if CSVs are in separate folder

# Pattern to identify gating columns (for automatic mode)
gating_pattern <- "\\(Manual gating\\)$"  # Adjust based on your CSV column naming

# Chunked processing settings (for large files)
use_chunked_processing <- FALSE  # Set TRUE for memory-efficient processing
chunk_size <- 100000  # Number of rows per chunk
auto_chunk_threshold_mb <- 100  # Auto-enable chunking if file > this size (MB)
max_combined_file_rows <- 5000000  # Skip combined file if total rows exceed this

# Parallel processing settings
use_parallel <- TRUE  # Set FALSE to disable parallel processing
max_cores <- NULL  # NULL = auto-detect, or specify number (e.g., 4)

# Output path
output_path <- paste0(input_fcs_path, "Expression_output/")

# Export format for expression matrices
export_individual_files <- TRUE  # Export one file per input CSV
export_combined_file <- FALSE    # Disabled - combined files cause memory crashes

# CAFE-format renaming settings
enable_cafe_renaming <- TRUE     # Enable two-stage CAFE format renaming
group_mapping_file <- "sample_group_mapping.xlsx"  # Excel file for user input

# ========== Manual Mode Settings (used when processing_mode = "manual") ==========
# Label translation dictionary - maps CSV column names to human-readable names
csv_label_dict <- c(
  "CD3n..Manual.gating." = "ncCD3-",
  "NKCs..Manual.gating." = "NK cells",
  "gdTCRnTCs..Manual.gating." = "ncT cells",
  "CD4TCs..Manual.gating." = "CD4T cells",
  "CD8TCs..Manual.gating." = "CD8T cells",
  "BCs..Manual.gating." = "B cells",
  "gdTCRnCD56nTCs..Manual.gating." = "ncT cells",
  "gdTCs..Manual.gating." = "gdT cells",
  "NKTCs..Manual.gating." = "NKT cells",
  "root" = "Removed"
)

# Gating hierarchy tree - defines parent-child relationships
# Higher numbers = deeper in hierarchy
gating_tree <- c(
  "Other" = 0,
  "T cells" = 1, "CD3-" = 1, "gdT cells" = 1,
  "NK cells" = 2, "NKT cells" = 2, "B cells" = 2,
  "CD8T cells" = 3, "CD4T cells" = 3
)

# ================================ Load Functions ==============================

source("Packages.R")
source("Functions.R")
dir_create(output_path)

# ======================== Validation & Setup ==================================

# Validate paths
if (!dir.exists(input_csv_path)) {
  stop("CSV input path does not exist: ", input_csv_path)
}

if (data_source_mode == "csv_with_fcs" && !dir.exists(input_fcs_path)) {
  stop("FCS input path does not exist: ", input_fcs_path)
}

cat("============================================================================\n")
cat("Load CSV Data - Unified Processing Pipeline\n")
cat("============================================================================\n")
cat("Data source mode:", data_source_mode, "\n")
cat("Processing mode:", processing_mode, "\n")
cat("CSV path:", input_csv_path, "\n")
if (data_source_mode == "csv_with_fcs") {
  cat("FCS path:", input_fcs_path, "\n")
}
cat("Output path:", output_path, "\n")
cat("Chunked processing:", ifelse(use_chunked_processing, 
                                   paste("YES (", chunk_size, "rows)"), "NO"), "\n")
if (!is.null(max_cores)) {
  n_cores <- min(max_cores, parallel::detectCores() - 1)
}
cat("Parallel processing:", ifelse(use_parallel, 
                                    paste("ENABLED (", n_cores, "cores)"), "DISABLED"), "\n")
cat("============================================================================\n\n")

#================================= Load Files ==================================

cat("[1/6] Discovering files...\n")

# Get list of CSV files
csv_files <- list.files(input_csv_path, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No CSV files found in: ", input_csv_path)
}

cat("  Found", length(csv_files), "CSV files\n")

# Check total file size and enable auto-chunking if needed
total_size_mb <- sum(file.info(csv_files)$size) / (1024^2)
cat("  Total CSV file size:", round(total_size_mb, 1), "MB\n")

if (!use_chunked_processing && total_size_mb > auto_chunk_threshold_mb) {
  cat("  ⚠ Large dataset detected - enabling chunked processing\n")
  use_chunked_processing <- TRUE
}

# Get list of FCS files if applicable
fcs_files <- NULL
if (data_source_mode == "auto") {
  fcs_files <- list.files(input_fcs_path, pattern = "\\.fcs$", full.names = TRUE)
  
  if (length(fcs_files) > 0) {
    cat("  Found", length(fcs_files), "FCS files - switching to csv_with_fcs mode\n")
    data_source_mode <- "csv_with_fcs"
  } else {
    cat("  No FCS files found - switching to csv_only mode\n")
    data_source_mode <- "csv_only"
  }
} else if (data_source_mode == "csv_with_fcs") {
  fcs_files <- list.files(input_fcs_path, pattern = "\\.fcs$", full.names = TRUE)
  if (length(fcs_files) == 0) {
    stop("csv_with_fcs mode selected but no FCS files found in: ", input_fcs_path)
  }
  cat("  Found", length(fcs_files), "FCS files\n")
}

#===================== CAFE Sample ID Mapping (MUST BE FIRST) =================

cafe_mapping <- NULL

if (enable_cafe_renaming && export_individual_files) {
  cat("\n[2/7] Preparing CAFE format mapping...\n")
  
  mapping_file_path <- file.path(output_path, group_mapping_file)
  
  # Extract sample IDs from input CSV files
  input_sample_ids <- sapply(csv_files, function(f) {
    tools::file_path_sans_ext(basename(f))
  })
  
  # Compact sample IDs: remove repeated substrings
  find_common_substrings <- function(strings) {
    if (length(strings) <= 1) return(character(0))
    
    # Split by underscores
    parts_list <- strsplit(strings, "_")
    
    # Find parts that appear in ALL strings
    all_parts <- unlist(parts_list)
    unique_parts <- unique(all_parts)
    
    common_parts <- character(0)
    for (part in unique_parts) {
      if (all(sapply(parts_list, function(p) part %in% p))) {
        common_parts <- c(common_parts, part)
      }
    }
    
    return(common_parts)
  }
  
  compact_sample_id <- function(sample_id, common_parts) {
    # Remove common substrings
    parts <- strsplit(sample_id, "_")[[1]]
    parts <- parts[!parts %in% common_parts]
    
    # Remove spaces and collapse
    compacted <- paste(parts, collapse="")
    compacted <- gsub("[[:space:]]", "", compacted)
    compacted <- gsub("-+", "-", compacted)  # Collapse multiple dashes
    
    # Remove leading zeros from each numeric segment
    # Split into character-by-character and remove leading zeros before digits
    chars <- strsplit(compacted, "")[[1]]
    result <- character(0)
    skip_zeros <- TRUE
    
    for (i in seq_along(chars)) {
      ch <- chars[i]
      
      # If it's a zero and we're in "skip zeros" mode
      if (ch == "0" && skip_zeros) {
        # Check if next char exists and is a digit
        if (i < length(chars) && grepl("[0-9]", chars[i + 1])) {
          next  # Skip this leading zero
        }
      }
      
      # Add the character
      result <- c(result, ch)
      
      # Reset skip_zeros when we hit a non-digit
      if (!grepl("[0-9]", ch)) {
        skip_zeros <- TRUE
      } else {
        skip_zeros <- FALSE
      }
    }
    
    compacted <- paste(result, collapse="")
    
    return(compacted)
  }
  
  common_parts <- find_common_substrings(input_sample_ids)
  if (length(common_parts) > 0) {
    cat("  Common substrings found:", paste(common_parts, collapse=", "), "\n")
  }
  
  compacted_ids <- sapply(input_sample_ids, function(id) compact_sample_id(id, common_parts))
  
  # Check if mapping file exists
  if (!file.exists(mapping_file_path)) {
    # Create new mapping file
    cat("  Creating group mapping file:", group_mapping_file, "\n")
    
    mapping_df <- data.frame(
      SampleID = compacted_ids,
      OriginalID = input_sample_ids,
      Group = "",
      stringsAsFactors = FALSE
    )
    
    openxlsx::write.xlsx(mapping_df, mapping_file_path, overwrite = TRUE)
    
    cat("\n")
    cat("  ╔════════════════════════════════════════════════════════════════╗\n")
    cat("  ║  ACTION REQUIRED: Fill in the 'Group' column                       ║\n")
    cat("  ╠════════════════════════════════════════════════════════════════╣\n")
    cat("  ║  1. Open:", format(basename(mapping_file_path), width=50), "║\n")
    cat("  ║  2. Fill the 'Group' column for each sample                        ║\n")
    cat("  ║  3. Save the file                                                  ║\n")
    cat("  ║  4. Re-run this script to generate files with CAFE format          ║\n")
    cat("  ╚════════════════════════════════════════════════════════════════╝\n")
    cat("\n")
    cat("  File location:", mapping_file_path, "\n")
    cat("\n")
    
    stop("\n\n  ⚠️  STOPPED: Mapping file created but Group column is EMPTY.\n  ⚠️  Please fill the 'Group' column and re-run the script.\n  ⚠️  NO OUTPUT FILES will be generated until mapping is complete.\n\n")
    
  } else {
    # Read existing mapping file
    cat("  Reading group mapping file:", group_mapping_file, "\n")
    
    mapping_df <- openxlsx::read.xlsx(mapping_file_path)
    
    # Check if Group column is filled
    if (!"Group" %in% colnames(mapping_df)) {
      stop("Mapping file is corrupted or invalid - please delete and re-run to recreate")
    }
    
    if (all(is.na(mapping_df$Group) | mapping_df$Group == "")) {
      cat("  File location:", mapping_file_path, "\n\n")
      stop("Group column is empty - script execution halted")
    }
    
    # Filter to samples with valid groups
    empty_groups <- sum(is.na(mapping_df$Group) | mapping_df$Group == "")
    if (empty_groups > 0) {
      cat("  ⚠ Warning:", empty_groups, "samples have no group defined - will skip those\n")
      mapping_df <- mapping_df[!is.na(mapping_df$Group) & mapping_df$Group != "", ]
    }
    
    # Convert to named list for easy lookup: cafe_mapping[[file_id]] = list(sample_id, group)
    cafe_mapping <- list()
    for (i in seq_len(nrow(mapping_df))) {
      original_id <- mapping_df$OriginalID[i]
      cafe_mapping[[original_id]] <- list(
        sample_id = mapping_df$SampleID[i],
        group = mapping_df$Group[i]
      )
    }
    
    cat("  ✓ Loaded", length(cafe_mapping), "sample-to-group mappings\n")
    cat("  Files will be written directly in CAFE format: {SampleID}_GroupX_{Population}.csv\n")
  }
}

# Detect if CSVs contain gating data
first_csv <- csv_files[1]
first_csv_cols <- names(data.table::fread(first_csv, nrows = 0))
has_gating_cols <- any(grepl(gating_pattern, first_csv_cols))

if (!has_gating_cols && processing_mode == "automatic") {
  cat("  No gating columns detected - switching to expression_only mode\n")
  processing_mode <- "expression_only"
} else if (has_gating_cols) {
  cat("  Detected gating columns in CSV files\n")
}

#===================== Create CSV to FCS Mapping ==============================

csv_to_fcs <- NULL
fcs_to_csvs <- NULL

if (data_source_mode == "csv_with_fcs") {
  if (enable_cafe_renaming && export_individual_files) {
    cat("\n[3/7] Creating CSV-to-FCS mapping...\n")
  } else {
    cat("\n[2/6] Creating CSV-to-FCS mapping...\n")
  }
  
  # Use the more robust mapping function from version 2
  csv_to_fcs <- create_csv_fcs_mapping(csv_files, fcs_files)
  
  # Group CSVs by FCS for summary
  fcs_to_csvs <- split(names(csv_to_fcs), csv_to_fcs)
  
  cat("  Mapping created for", length(csv_to_fcs), "CSV files\n")
  cat("  Mapped to", length(fcs_to_csvs), "unique FCS files\n\n")
  
  # Print detailed mapping
  cat("CSV to FCS mapping:\n")
  for (fcs in names(fcs_to_csvs)) {
    cat("  ", basename(fcs), ":\n")
    for (csv in fcs_to_csvs[[fcs]]) {
      cat("    -", basename(csv), "\n")
    }
  }
} else {
  if (enable_cafe_renaming && export_individual_files) {
    cat("\n[3/7] CSV-only mode - no FCS mapping needed\n")
  } else {
    cat("\n[2/6] CSV-only mode - no FCS mapping needed\n")
  }
  cat("  Each CSV will be processed independently\n")
}

#===================== Column Discovery & Selection ===========================

if (enable_cafe_renaming && export_individual_files) {
  cat("\n[4/7] Discovering columns...\n")
} else {
  cat("\n[3/6] Discovering columns...\n")
}

# Read first CSV to discover columns
first_csv <- csv_files[1]
separator <- detect_csv_sep(first_csv)
colnames_all <- names(data.table::fread(first_csv, nrows = 0))

cat("  Total columns in CSV:", length(colnames_all), "\n")

# Identify gating columns and markers based on processing mode
if (processing_mode == "expression_only") {
  # No gating data - all columns are markers or metadata
  cat("  Expression-only mode: treating all data as single population\n")
  
  # Exclude common metadata columns
  metadata_patterns <- c("^index$", "^Index$", "^ID$", "^id$", "^File", "^file", 
                        "^Time$", "^time$", "^Event", "^event")
  metadata_cols <- character(0)
  for (pattern in metadata_patterns) {
    metadata_cols <- c(metadata_cols, grep(pattern, colnames_all, value = TRUE))
  }
  metadata_cols <- unique(metadata_cols)
  
  candidate_marker_cols <- setdiff(colnames_all, metadata_cols)
  candidate_marker_cols <- candidate_marker_cols[!is.na(candidate_marker_cols) & 
                                                  candidate_marker_cols != ""]
  
  cat("  Candidate marker columns:", length(candidate_marker_cols), "\n")
  cat("  Metadata columns excluded:", length(metadata_cols), "\n")
  
  # Interactive marker selection
  cat("\n  Opening interactive marker selection dialog...\n")
  sel_markers <- tk_multiselect(candidate_marker_cols, "Select markers for expression matrix")
  if (!length(sel_markers)) stop("No markers selected – aborting.")
  
  cat("  Selected markers:", length(sel_markers), "\n")
  
  sel_gates <- NULL
  sel_gates_display <- "All_cells"
  
} else if (processing_mode == "automatic") {
  gate_cols <- grep(gating_pattern, colnames_all, value = TRUE)
  candidate_marker_cols <- setdiff(colnames_all, gate_cols)
  
  # Clean candidate columns
  candidate_marker_cols <- candidate_marker_cols[!is.na(candidate_marker_cols) & 
                                                  candidate_marker_cols != ""]
  
  cat("  Gating columns:", length(gate_cols), "\n")
  cat("  Candidate marker columns:", length(candidate_marker_cols), "\n")
  
  # Interactive selection
  cat("\n  Opening interactive marker selection dialog...\n")
  sel_markers <- tk_multiselect(candidate_marker_cols, "Select markers to analyse")
  if (!length(sel_markers)) stop("No markers selected – aborting.")
  
  cat("  Selected markers:", length(sel_markers), "\n")
  
  gate_display <- gsub(gating_pattern, "", gate_cols)
  cat("\n  Opening interactive gate selection dialog...\n")
  sel_gates_display <- tk_multiselect(gate_display, "Select cell populations to analyse")
  if (!length(sel_gates_display)) stop("No gates selected – aborting.")
  
  # Map back to original gate names
  sel_gates <- gate_cols[gate_display %in% sel_gates_display]
  cat("  Selected gates:", length(sel_gates), "\n")
  
} else {
  # Manual mode - use predefined dictionary
  cat("  Using manual mode with predefined labels\n")
  check_gates <- colnames(read.csv(first_csv, sep = separator, header = TRUE))
  sel_gates <- names(csv_label_dict)[names(csv_label_dict) %in% check_gates]
  sel_markers <- NULL  # Will extract all non-gating columns
}

#===================== Main Processing =========================================

if (enable_cafe_renaming && export_individual_files) {
  cat("\n[5/7] Processing CSV files...\n")
} else {
  cat("\n[4/6] Processing CSV files...\n")
}

if (processing_mode == "expression_only") {
  # ========== Expression-only mode: no gating data ===========
  cat("  Processing expression matrices without gating data\n")
  
  # Function to process a single file
  process_expr_only_file <- function(csv_file) {
    tryCatch({
      separator <- detect_csv_sep(csv_file)
      
      # Read full data to check original dimensions
      csv_full <- data.table::fread(csv_file, showProgress = FALSE)
      original_rows <- nrow(csv_full)
      original_cols <- ncol(csv_full)
      rm(csv_full)  # Free memory immediately
      gc(verbose = FALSE)
      
      # Read selected markers only
      csv_data <- data.table::fread(csv_file, select = sel_markers, showProgress = FALSE)
      
      # Get file identifier for validation tracking
      file_id <- tools::file_path_sans_ext(basename(csv_file))
      
      # Validate row counts
      output_rows <- nrow(csv_data)
      row_check <- (original_rows == output_rows)
      col_check <- (length(sel_markers) == original_cols)
      
      validation <- list(
        original_rows = original_rows,
        output_rows = output_rows,
        original_cols = original_cols,
        selected_cols = length(sel_markers),
        row_match = row_check,
        col_match = col_check
      )
      
      # Export individual file if requested
      if (export_individual_files) {
        if (!is.null(cafe_mapping) && file_id %in% names(cafe_mapping)) {
          # Use CAFE format: {SampleID}_{Group}.csv
          cafe_info <- cafe_mapping[[file_id]]
          output_filename <- paste0(cafe_info$sample_id, "_", cafe_info$group, ".csv")
        } else {
          # Use standard format
          output_filename <- paste0(file_id, "_All_cells_expression.csv")
        }
        output_file <- file.path(output_path, output_filename)
        data.table::fwrite(csv_data, output_file)
      }
      
      return(list(data = csv_data, validation = validation, file_id = file_id))
      
    }, error = function(e) {
      cat("\n  ERROR processing", basename(csv_file), ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Process files sequentially (safer for large files)
  # Parallel processing disabled in this section due to memory issues with large datasets
  cat("  Using sequential processing (memory-safe mode)\n")
  
  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent :current/:total ETA: :eta",
    total = length(csv_files), clear = FALSE, width = 80
  )
  
  results <- list()
  for (i in seq_along(csv_files)) {
    pb$tick()
    results[[i]] <- process_expr_only_file(csv_files[i])
    
    # Force garbage collection after each file
    gc(verbose = FALSE)
    
    # Optional: brief pause to allow system to release memory
    Sys.sleep(0.1)
  }
  
  # Extract results
  results <- results[!sapply(results, is.null)]
  file_validation <- lapply(results, function(x) {
    v <- x$validation
    names(v) <- NULL
    v
  })
  names(file_validation) <- sapply(results, function(x) x$file_id)
  
  cat("  Processed", length(results), "files successfully\n")
  
  # Clean up
  rm(results)
  gc(verbose = FALSE)
  
  # Validation report
  cat("\n  Row/Column Validation:\n")
  for (file_id in names(file_validation)) {
    v <- file_validation[[file_id]]
    row_status <- if(v$row_match) "✓" else "✗"
    col_status <- if(v$col_match) "✓" else "✗"
    cat(sprintf("    %s %s: %d/%d rows %s, %d/%d cols %s\n",
                row_status, file_id, v$output_rows, v$original_rows,
                if(v$row_match) "OK" else "MISMATCH",
                v$selected_cols, v$original_cols,
                if(v$col_match) "OK" else "SUBSET"))
  }
  
} else if (processing_mode == "automatic" && use_chunked_processing) {
  # ========== Automatic mode with data.table (memory efficient) ===========
  cat("  Using data.table with chunked processing\n")
  
  need_cols <- c(sel_markers, sel_gates)
  need_cols <- need_cols[!is.na(need_cols) & need_cols != ""]
  
  # Function to read a single file
  read_file_with_validation <- function(f) {
    tryCatch({
      # Read full data to track original dimensions
      csv_full <- data.table::fread(f, showProgress = FALSE)
      original_rows <- nrow(csv_full)
      original_cols <- ncol(csv_full)
      rm(csv_full)  # Free memory immediately
      gc(verbose = FALSE)
      
      dt <- data.table::fread(f, select = need_cols, showProgress = FALSE)
      file_id <- tools::file_path_sans_ext(basename(f))
      dt[, file := file_id]
      
      return(list(
        file_id = file_id,
        data = dt,
        original_rows = original_rows,
        original_cols = original_cols
      ))
      
    }, error = function(e) {
      cat("  ERROR reading", basename(f), ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Read files sequentially (memory-safe for large datasets)
  cat("  Reading files sequentially (memory-safe mode)\n")
  
  file_results <- list()
  for (i in seq_along(csv_files)) {
    file_results[[i]] <- read_file_with_validation(csv_files[i])
    
    # Force garbage collection after each file
    gc(verbose = FALSE)
  }
  
  # Process results
  file_results <- file_results[!sapply(file_results, is.null)]
  file_data_list <- list()
  file_validation <- list()
  
  for (result in file_results) {
    file_data_list[[result$file_id]] <- list(
      data = result$data,
      original_rows = result$original_rows,
      original_cols = result$original_cols
    )
  }
  
  # Export individual files per gate/population
  if (export_individual_files) {
    cat("\n  Exporting individual expression matrices per file and population...\n")
    
    # Function to export all gates for one file
    export_file_gates <- function(file_id) {
      file_info <- file_data_list[[file_id]]
      dt <- file_info$data
      validations <- list()
      
      for (gate_idx in seq_along(sel_gates)) {
        gate_col <- sel_gates[gate_idx]
        gate_name <- sel_gates_display[gate_idx]
        
        # Extract data for this gate
        gate_data <- dt[get(gate_col) == 1, sel_markers, with = FALSE]
        
        # Validate row counts
        gate_rows <- nrow(gate_data)
        gate_cells <- sum(dt[[gate_col]] == 1)
        
        if (gate_rows != gate_cells) {
          cat(sprintf("\n  Warning: %s - %s: row mismatch %d != %d\n",
                     file_id, gate_name, gate_rows, gate_cells))
        }
        
        # Save individual file with CAFE format
        if (!is.null(cafe_mapping) && file_id %in% names(cafe_mapping)) {
          # Use CAFE format: {SampleID}_Group{Group}_{Population}.csv
          cafe_info <- cafe_mapping[[file_id]]
          safe_gate_name <- gsub("[^A-Za-z0-9_-]", "_", gate_name)
          output_filename <- paste0(cafe_info$sample_id, "_Group", cafe_info$group, "_", safe_gate_name, ".csv")
        } else {
          # Use standard format
          safe_gate_name <- gsub("[^A-Za-z0-9_-]", "_", gate_name)
          output_filename <- paste0(file_id, "_", safe_gate_name, "_expression.csv")
        }
        output_file <- file.path(output_path, output_filename)
        data.table::fwrite(gate_data, output_file)
        
        # Track validation
        validation_key <- paste0(file_id, "_", gate_name)
        validations[[validation_key]] <- list(
          file = file_id,
          population = gate_name,
          original_total_rows = file_info$original_rows,
          gated_rows = gate_rows,
          selected_cols = length(sel_markers),
          original_cols = file_info$original_cols
        )
      }
      return(validations)
    }
    
    # Process sequentially (memory-safe for large datasets)
    cat("  Sequential export (memory-safe mode)\n")
    individual_pb <- progress::progress_bar$new(
      format = "  [:bar] :percent :current/:total",
      total = length(file_data_list), clear = FALSE, width = 80
    )
    
    file_validation <- list()
    for (file_id in names(file_data_list)) {
      individual_pb$tick()
      validations <- export_file_gates(file_id)
      file_validation <- c(file_validation, validations)
      
      # Force garbage collection after each file
      gc(verbose = FALSE)
    }
  }
  
  # Clean up file data to free memory
  cat("\n  Cleaning up memory...\n")
  rm(file_data_list)
  gc(verbose = FALSE)
  cat("  Individual files per population have been saved\n")
  
  # Validation report
  if (export_individual_files && length(file_validation) > 0) {
    cat("\n  Individual File Validation:\n")
    for (v_key in names(file_validation)) {
      v <- file_validation[[v_key]]
      cat(sprintf("    %s [%s]: %d/%d rows (%.1f%%), %d cols\n",
                  v$file, v$population, v$gated_rows, v$original_total_rows,
                  100 * v$gated_rows / v$original_total_rows, v$selected_cols))
    }
  }
  
} else if (processing_mode == "manual") {
  # ========== Manual mode with traditional processing ===========
  cat("  Using manual mode with label dictionaries\n")
  
  csv_data <- list()
  label_lists <- list()
  flat_list <- list()
  idx_list <- list()
  idx_labels <- list()
  
  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent :current/:total ETA: :eta",
    total = length(csv_files), clear = FALSE, width = 80
  )
  
  for (csv_file in csv_files) {
    pb$tick()
    
    tryCatch({
      separator <- detect_csv_sep(csv_file)
      csv_data[[csv_file]] <- read.csv(csv_file, sep = separator, header = TRUE)
      
      # Validate columns
      if (!identical(colnames(csv_data[[csv_file]]), check_gates)) {
        diff1 <- setdiff(check_gates, colnames(csv_data[[csv_file]]))
        diff2 <- setdiff(colnames(csv_data[[csv_file]]), check_gates)
        
        if (length(diff1) > 0 || length(diff2) > 0) {
          cat("\n  Warning: Column mismatch in", basename(csv_file), "\n")
          if (length(diff1) > 0) cat("    Missing:", paste(diff1, collapse = ", "), "\n")
          if (length(diff2) > 0) cat("    Extra:", paste(diff2, collapse = ", "), "\n")
        }
      }
      
      idx_list[[csv_file]] <- csv_data[[csv_file]][, 1]
      colnames(csv_data[[csv_file]]) <- translate_labels(
        colnames(csv_data[[csv_file]]), csv_label_dict)
      label_lists[[csv_file]] <- labels_from_binary_df(csv_data[[csv_file]])
      flat_list[[csv_file]] <- flatten_label_list(label_lists[[csv_file]],
                                                   gating_tree)
      
      # Extract the correct CSV-associated FCS file
      fcs_file <- csv_to_fcs[[csv_file]]
      num_events <- nrow(flowCore::read.FCS(fcs_file, truncate_max_range = FALSE))
      idx_labels_vec <- rep("", num_events)
      idx_labels_vec[idx_list[[csv_file]] + 1] <- translate_labels(
        flat_list[[csv_file]], csv_label_dict)
      
      idx_labels[[csv_file]] <- data.frame(idx_labels_vec, stringsAsFactors = FALSE)
      
    }, error = function(e) {
      cat("\n  ERROR processing", basename(csv_file), ":", e$message, "\n")
    })
  }
  
  cat("\n  Processed", length(idx_labels), "files successfully\n")
  
  # Save results for each file
  for (csv_file in names(idx_labels)) {
    output_file <- file.path(output_path, paste0(
      tools::file_path_sans_ext(basename(csv_file)), "_labels.csv"))
    write.csv(idx_labels[[csv_file]], output_file, row.names = FALSE)
  }
  
} else {
  # ========== Automatic mode without chunking ===========
  cat("  Using automatic mode with standard processing\n")
  
  # Function to process a single file
  process_auto_file <- function(csv_file) {
    file_id <- tools::file_path_sans_ext(basename(csv_file))
    
    tryCatch({
      # Check file size first
      file_size_mb <- file.info(csv_file)$size / (1024^2)
      if (file_size_mb > 500) {
        cat("\n  Warning: Large file (", round(file_size_mb, 1), "MB):", basename(csv_file), "\n")
      }
      
      # Read original data for validation
      csv_full <- data.table::fread(csv_file, showProgress = FALSE)
      original_rows <- nrow(csv_full)
      original_cols <- ncol(csv_full)
      rm(csv_full)  # Free memory immediately
      gc(verbose = FALSE)
      
      # Read with selected columns and gates
      separator <- detect_csv_sep(csv_file)
      csv_data <- data.table::fread(csv_file, select = c(sel_markers, sel_gates), 
                                    showProgress = FALSE)
      
      # Export individual expression matrices per gate
      validations <- list()
      if (export_individual_files) {
        for (gate_idx in seq_along(sel_gates)) {
          gate_col <- sel_gates[gate_idx]
          gate_name <- sel_gates_display[gate_idx]
          
          # Extract gated data
          gate_data <- csv_data[get(gate_col) == 1, sel_markers, with = FALSE]
          
          # Validate
          gate_rows <- nrow(gate_data)
          validation_key <- paste0(file_id, "_", gate_name)
          validations[[validation_key]] <- list(
            file = file_id,
            population = gate_name,
            original_total_rows = original_rows,
            gated_rows = gate_rows,
            selected_cols = length(sel_markers),
            original_cols = original_cols
          )
          
          # Save with CAFE format
          if (!is.null(cafe_mapping) && file_id %in% names(cafe_mapping)) {
            # Use CAFE format: {SampleID}_Group{Group}_{Population}.csv
            cafe_info <- cafe_mapping[[file_id]]
            safe_gate_name <- gsub("[^A-Za-z0-9_-]", "_", gate_name)
            output_filename <- paste0(cafe_info$sample_id, "_Group", cafe_info$group, "_", safe_gate_name, ".csv")
          } else {
            # Use standard format
            safe_gate_name <- gsub("[^A-Za-z0-9_-]", "_", gate_name)
            output_filename <- paste0(file_id, "_", safe_gate_name, "_expression.csv")
          }
          output_file <- file.path(output_path, output_filename)
          data.table::fwrite(gate_data, output_file)
        }
      }
      
      # Process using original function for label-based output
      result_file <- process_single_csv(
        csv_files = csv_file,
        pattern = gating_pattern,
        csv_to_fcs = csv_to_fcs,
        csv_label_dict = csv_label_dict,
        gating_tree = gating_tree,
        output_path = output_path,
        chunk_size = chunk_size
      )
      
      return(list(
        result_file = result_file,
        validations = validations,
        file_name = basename(csv_file)
      ))
      
    }, error = function(e) {
      cat("\n  ERROR:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Process files sequentially (memory-safe for large datasets)
  cat("  Sequential processing (memory-safe mode)\n")
  
  pb <- progress::progress_bar$new(
    format = "  [:bar] :percent :current/:total ETA: :eta",
    total = length(csv_files), clear = FALSE, width = 80
  )
  
  file_results <- list()
  for (i in seq_along(csv_files)) {
    pb$tick()
    file_results[[i]] <- process_auto_file(csv_files[i])
    
    # Force garbage collection after each file
    gc(verbose = FALSE)
    
    # Optional: brief pause to allow system to release memory
    Sys.sleep(0.1)
  }
  
  # Extract results
  file_results <- file_results[!sapply(file_results, is.null)]
  results <- lapply(file_results, function(x) x$result_file)
  names(results) <- sapply(file_results, function(x) x$file_name)
  
  file_validation <- list()
  for (result in file_results) {
    if (!is.null(result$validations)) {
      file_validation <- c(file_validation, result$validations)
    }
  }
  
  successful_results <- results[!sapply(results, is.null)]
  cat("\n  Successfully processed", length(successful_results), "files\n")
  
  # Validation report
  if (export_individual_files && length(file_validation) > 0) {
    cat("\n  Individual File Validation:\n")
    for (v_key in names(file_validation)) {
      v <- file_validation[[v_key]]
      cat(sprintf("    %s [%s]: %d/%d rows (%.1f%%), %d cols\n",
                  v$file, v$population, v$gated_rows, v$original_total_rows,
                  100 * v$gated_rows / v$original_total_rows, v$selected_cols))
    }
  }
}

#===================== Validation & Summary ====================================

if (enable_cafe_renaming && export_individual_files) {
  cat("\n[6/7] Validating results...\n")
} else {
  cat("\n[5/6] Validating results...\n")
}

if (processing_mode == "expression_only") {
  if (export_individual_files && exists("file_validation")) {
    cat("  Individual expression files:", length(file_validation), "\n")
    
    # Summary validation
    all_rows_match <- all(sapply(file_validation, function(v) v$row_match))
    if (all_rows_match) {
      cat("  ✓ All files: input rows = output rows (no filtering applied)\n")
    } else {
      mismatches <- sum(!sapply(file_validation, function(v) v$row_match))
      cat("  ✗ Warning:", mismatches, "files have row count mismatches\n")
    }
  }
  
} else if (processing_mode == "automatic" && use_chunked_processing) {
  if (export_individual_files && exists("file_validation")) {
    total_individual_files <- length(file_validation)
    cat("  Individual expression files created:", total_individual_files, "\n")
  }
  
} else if (processing_mode == "manual") {
  total_labels <- sum(sapply(idx_labels, function(x) sum(x[[1]] != "")))
  cat("  Total labeled events:", total_labels, "\n")
  cat("  Files with labels:", length(idx_labels), "\n")
  
} else {
  if (exists("successful_results")) {
    validate_csv_results(successful_results)
  }
  if (export_individual_files && exists("file_validation")) {
    total_individual_files <- length(file_validation)
    cat("  Individual expression files created:", total_individual_files, "\n")
  }
}

#===================== Final Summary ===========================================

if (enable_cafe_renaming && export_individual_files) {
  cat("\n[7/7] Complete!\n")
} else {
  cat("\n[6/6] Complete!\n")
}
cat("============================================================================\n")
cat("Processing Summary:\n")
cat("  Data source mode:", data_source_mode, "\n")
cat("  Processing mode:", processing_mode, "\n")
cat("  CSV files processed:", length(csv_files), "\n")
cat("  Output directory:", output_path, "\n")

if (processing_mode == "expression_only") {
  if (export_individual_files && exists("file_validation")) {
    cat("  Individual expression files:", length(file_validation), "\n")
    if (enable_cafe_renaming && !is.null(cafe_mapping)) {
      cat("  Format: {SampleID}_GroupX_AllCells.csv\n")
    } else {
      cat("  Format: {filename}_All_cells_expression.csv\n")
    }
  }
} else if (processing_mode == "automatic") {
  if (use_chunked_processing) {
    if (export_individual_files && exists("file_validation")) {
      cat("  Individual expression files:", length(file_validation), "\n")
      if (enable_cafe_renaming && !is.null(cafe_mapping)) {
        cat("  Format: {SampleID}_GroupX_{Population}.csv\n")
      } else {
        cat("  Format: {filename}_{population}_expression.csv\n")
      }
    }
  } else {
    cat("  Label files:", length(successful_results), "\n")
    if (export_individual_files && exists("file_validation")) {
      cat("  Individual expression files:", length(file_validation), "\n")
      if (enable_cafe_renaming && !is.null(cafe_mapping)) {
        cat("  Format: {SampleID}_GroupX_{Population}.csv\n")
      } else {
        cat("  Format: {filename}_{population}_expression.csv\n")
      }
    }
  }
} else {
  cat("  Label files created:", length(idx_labels), "\n")
}

cat("============================================================================\n")
cat("\nAll processing complete. Results saved to:\n")
cat("  ", output_path, "\n\n")