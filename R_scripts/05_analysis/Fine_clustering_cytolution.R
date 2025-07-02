#!/usr/bin/env Rscript

# =========================
# Fine Clustering Cytolution Analysis
# Statistical analysis of cluster enrichment between conditions
# =========================

# =========================
# CONFIG — EDIT THIS BLOCK
# =========================
CONFIG <- list(
  # 1) Folder with per-cell CSVs
  data_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR",

  # 2) File condition mapping.
  #    Option A: Create a CSV "file_conditions.csv" in data_dir with columns: file_name,condition
  #              (file_name must match the CSV filename without extension).
  #    Option B: Leave as NULL and the script will scaffold a template for you then stop.
  conditions_csv = NULL,  # e.g., "/data/.../file_conditions.csv"

  # 3) UMAP columns - will be selected via GUI
  umap_x = NULL,  # Selected interactively
  umap_y = NULL,  # Selected interactively

  # 4) Regex to pick cluster columns (binary 0/1 membership).
  #    Example: "^cluster_" or "(?i)meta.*cluster" or "^(?:C[0-9]+)$"
  cluster_regex = "Clustering annotated",

  # 5) Markers of interest - will be selected via GUI
  markers_of_interest = NULL,  # Selected interactively

  # 6) Two condition labels (order used for log2 enrichment sign / plots)
  condition_levels = c("pCR", "npCR"),

  # 7) Output directory (will be created)
  out_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/",

  # 8) Minimum fraction of cells in a cluster to keep (filters ultra-rare clusters)
  min_cluster_frac = 0.001,  # 0.1%

  # 9) Multiple-testing correction method ("BH" recommended)
  p_adjust_method = "BH",

  # 10) Statistical method: "pooled" (faster) or "mixed_effects" (accounts for file-level variation)
  stat_method = "mixed_effects",

  # 11) Random seed (for reproducible sampling in plots if needed)
  seed = 123,
  
  # 12) Memory management for large files
  chunk_size = 50000,      # Reduced chunk size
  max_memory_gb = 4,       # Conservative memory limit
  max_cells_per_file = 500000,  # Skip files larger than this
  sample_large_files = TRUE     # Sample large files instead of loading entirely
)

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)   # fast CSV read
  library(janitor)
  library(stringr)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(lme4)         # mixed effects models
  library(emmeans)      # estimated marginal means
  library(ggrepel)      # text repel for plots
  library(viridis)      # color scales
  library(progress)     # progress bars
  library(tcltk)        # GUI components
})

set.seed(CONFIG$seed)

# =========================
# Helper Functions
# =========================
die <- function(...) { message(...); quit(status = 1) }

make_outdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

list_csvs <- function(dir_path) {
  csvs <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) die("No .csv files found in: ", dir_path)
  tibble(file = csvs, file_name = tools::file_path_sans_ext(basename(csvs)))
}

# Memory-efficient CSV reading with chunking for large files
read_csv_chunked <- function(file_path, chunk_size = CONFIG$chunk_size) {
  file_size_mb <- file.info(file_path)$size / (1024^2)
  
  if (file_size_mb > CONFIG$max_memory_gb * 1000) {
    message("Large file detected (", round(file_size_mb, 1), " MB). Using chunked reading...")
    
    # Read header first
    header <- fread(file_path, nrows = 0)
    
    # Safe row counting with fallbacks
    total_rows <- tryCatch({
      # Method 1: Try wc command (Linux/Mac)
      if (file.exists(file_path)) {
        wc_result <- system2("wc", args = c("-l", shQuote(file_path)), stdout = TRUE, stderr = FALSE)
        if (length(wc_result) > 0 && !is.na(wc_result)) {
          as.numeric(str_extract(wc_result, "^\\d+"))
        } else {
          stop("wc command failed")
        }
      } else {
        stop("File does not exist")
      }
    }, error = function(e) {
      message("wc command failed, using R fallback method...")
      # Method 2: R fallback - count lines by reading file
      tryCatch({
        length(readLines(file_path))
      }, error = function(e2) {
        message("Line counting failed, estimating from file size...")
        # Method 3: Rough estimate based on file size
        round(file.size(file_path) / 1000)  # Assume ~1000 bytes per row
      })
    })
    
    chunks <- list()
    n_chunks <- ceiling(total_rows / chunk_size)
    
    pb <- progress_bar$new(
      format = "Reading [:bar] :percent :current/:total chunks",
      total = n_chunks, clear = FALSE
    )
    
    for (i in seq_len(n_chunks)) {
      skip_rows <- (i - 1) * chunk_size
      chunk <- fread(file_path, skip = skip_rows + 1, nrows = chunk_size, header = FALSE)
      
      if (nrow(chunk) > 0) {
        names(chunk) <- names(header)
        chunks[[i]] <- chunk
      }
      pb$tick()
    }
    
    return(rbindlist(chunks, fill = TRUE))
  } else {
    return(fread(file_path))
  }
}

# Check memory usage
check_memory <- function() {
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -m", intern = TRUE)
    mem_line <- mem_info[2]
    mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
    used_pct <- mem_values[2] / mem_values[1]
    
    if (used_pct > 0.85) {
      message("Warning: High memory usage (", round(used_pct * 100, 1), "%). Running garbage collection...")
      gc()
      return(FALSE)  # Signal memory pressure
    }
  }
  return(TRUE)  # Memory OK
}

# GUI for column selection
select_columns_gui <- function(column_names, title = "Select Columns", multiple = TRUE) {
  require(tcltk)
  
  if (length(column_names) == 0) {
    message("No columns available for selection")
    return(character(0))
  }
  
  # Create main window
  tt <- tktoplevel()
  tkwm.title(tt, title)
  tkwm.geometry(tt, "500x400")
  
  # Main frame
  main_frame <- tkframe(tt)
  tkpack(main_frame, fill = "both", expand = TRUE, padx = 10, pady = 10)
  
  # Instructions
  instr_text <- if (multiple) "Select multiple columns (Ctrl+click):" else "Select one column:"
  instr_label <- tklabel(main_frame, text = instr_text, 
                         font = tkfont.create(size = 10, weight = "bold"))
  tkpack(instr_label, pady = 5)
  
  # Listbox with scrollbar
  list_frame <- tkframe(main_frame)
  tkpack(list_frame, fill = "both", expand = TRUE)
  
  scrollbar <- tkscrollbar(list_frame, repeatinterval = 5,
                          command = function(...) tkyview(lb, ...))
  select_mode <- if (multiple) "extended" else "single"
  lb <- tklistbox(list_frame, yscrollcommand = function(...) tkset(scrollbar, ...),
                  selectmode = select_mode, height = 15)
  
  tkpack(scrollbar, side = "right", fill = "y")
  tkpack(lb, side = "left", fill = "both", expand = TRUE)
  
  # Add columns to listbox
  for (col in column_names) {
    tkinsert(lb, "end", col)
  }
  
  # Buttons
  button_frame <- tkframe(main_frame)
  tkpack(button_frame, pady = 10)
  
  selected_items <- character(0)
  
  ok_function <- function() {
    selections <- as.integer(tkcurselection(lb))
    if (length(selections) > 0) {
      selected_items <<- column_names[selections + 1]
    }
    tkdestroy(tt)
  }
  
  cancel_function <- function() {
    selected_items <<- character(0)
    tkdestroy(tt)
  }
  
  ok_button <- tkbutton(button_frame, text = "OK", command = ok_function)
  cancel_button <- tkbutton(button_frame, text = "Cancel", command = cancel_function)
  
  tkpack(ok_button, side = "left", padx = 5)
  tkpack(cancel_button, side = "left", padx = 5)
  
  # Focus and wait
  tkfocus(lb)
  tkwait.window(tt)
  
  return(selected_items)
}

# Get file info without loading entire file
get_file_info <- function(file_path) {
  tryCatch({
    # Read just the header
    header <- fread(file_path, nrows = 0)
    
    # Get approximate row count with safe fallbacks
    row_count <- tryCatch({
      if (file.exists(file_path)) {
        # Try wc command with proper quoting
        wc_result <- system2("wc", args = c("-l", shQuote(file_path)), stdout = TRUE, stderr = FALSE)
        if (length(wc_result) > 0 && !is.na(wc_result)) {
          as.numeric(str_extract(wc_result, "^\\d+")) - 1  # Subtract header
        } else {
          stop("wc command failed")
        }
      } else {
        stop("File does not exist")
      }
    }, error = function(e) {
      message("wc command failed for ", basename(file_path), ", using R fallback...")
      # Fallback method using data.table
      tryCatch({
        nrow(fread(file_path, select = 1))
      }, error = function(e2) {
        message("Row counting failed, estimating...")
        # Rough estimate
        max(1000, round(file.size(file_path) / 800))
      })
    })
    
    file_size_mb <- file.info(file_path)$size / (1024^2)
    
    list(
      columns = names(header),
      n_rows = row_count,
      size_mb = file_size_mb,
      too_large = row_count > CONFIG$max_cells_per_file
    )
  }, error = function(e) {
    message("Error reading file info for ", basename(file_path), ": ", e$message)
    list(columns = character(0), n_rows = 0, size_mb = 0, too_large = TRUE)
  })
}

# =========================
# Statistical Functions
# =========================

# Fisher's exact test for pooled analysis
do_fisher <- function(vec01, cond_factor, cond_levels) {
  stopifnot(length(vec01) == length(cond_factor))
  inC <- vec01 == 1
  condA <- cond_factor == cond_levels[1]
  condB <- cond_factor == cond_levels[2]
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

# =========================
# Data Loading and Processing
# =========================

# Load files + metadata with error handling
tryCatch({
  files_tbl <- list_csvs(CONFIG$data_dir)
  message("Found ", nrow(files_tbl), " CSV files")
}, error = function(e) {
  message("Error listing CSV files: ", e$message)
  message("Please check:")
  message("1. CONFIG$data_dir path is correct: ", CONFIG$data_dir)
  message("2. Directory exists and contains .csv files")
  message("3. You have read permissions to the directory")
  quit(save = "no", status = 1)
})

# Function to safely check for existing conditions file
check_conditions_file_exists <- function() {
  message("Checking for existing conditions file...")
  
  if (!is.null(CONFIG$conditions_csv)) {
    message("✓ Conditions CSV already specified: ", CONFIG$conditions_csv)
    return(CONFIG$conditions_csv)
  }
  
  # Look for conditions file in safe locations (NOT in data directory to avoid reading as input)
  potential_paths <- c(
    file.path(CONFIG$out_dir, "file_conditions.csv"),    # Output directory (preferred)
    file.path(getwd(), "file_conditions.csv"),          # Working directory  
    file.path("~", "file_conditions.csv"),              # Home directory
    file.path(dirname(CONFIG$data_dir), "file_conditions.csv")  # Parent of data dir
  )
  
  for (path in potential_paths) {
    if (file.exists(path)) {
      message("Found file: ", path)
      tryCatch({
        test_data <- read.csv(path, stringsAsFactors = FALSE, nrows = 5)  # Only read first few rows
        if ("condition" %in% names(test_data) && "file_name" %in% names(test_data)) {
          # Check if file has some non-NA conditions
          full_data <- read.csv(path, stringsAsFactors = FALSE)
          if (!all(is.na(full_data$condition))) {
            message("✓ Found valid conditions file: ", path)
            CONFIG$conditions_csv <<- path
            return(path)
          }
        }
      }, error = function(e) {
        message("Cannot read file: ", e$message)
      })
    }
  }
  
  return(NULL)
}

# Safe scaffolding that doesn't crash the session
safe_scaffold_check <- function() {
  existing_file <- check_conditions_file_exists()
  
  if (!is.null(existing_file)) {
    message("Using existing conditions file: ", existing_file)
    return(TRUE)  # Success - can continue
  }
  
  # If no file exists, provide instructions but don't crash
  message("\n" , "="*50)
  message("SETUP REQUIRED: No conditions file found")
  message("="*50)
  message("Please create a file called 'file_conditions.csv' with this structure:")
  message("")
  message("file_name,condition")
  message("your_file1.csv,pCR")
  message("your_file2.csv,npCR")
  message("...")
  message("")
  message("RECOMMENDED LOCATION (to avoid interference with input data):")
  message("Save to: ", file.path(CONFIG$out_dir, "file_conditions.csv"))
  message("Or alternatively to: ", file.path(getwd(), "file_conditions.csv"))
  message("")
  message("Then set: CONFIG$conditions_csv <- '", file.path(CONFIG$out_dir, "file_conditions.csv"), "'")
  message("And re-run this script.")
  message("="*50)
  
  return(FALSE)  # Cannot continue, but don't crash
}

# Check if we can proceed
message("Step: Checking conditions file setup...")
can_proceed <- safe_scaffold_check()

if (!can_proceed) {
  message("Script stopping safely - please create conditions file as instructed above.")
  # Don't use quit() - just stop execution here with a clear message
  stop("Please create conditions file and re-run script. See instructions above.")
}

# Read and validate conditions file
conds <- tryCatch({
  temp_conds <- read_csv(CONFIG$conditions_csv, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(across(everything(), as.character))
  
  # Validate that conditions are not all NAs
  if (all(is.na(temp_conds$condition))) {
    stop("All conditions are NA. Please fill in the condition column with valid values.")
  }
  
  # Check for valid condition values
  valid_conditions <- temp_conds$condition[!is.na(temp_conds$condition)]
  if (length(valid_conditions) == 0) {
    stop("No valid conditions found. Please fill in the condition column.")
  }
  
  # Warn about any remaining NAs
  na_count <- sum(is.na(temp_conds$condition))
  if (na_count > 0) {
    message("Warning: ", na_count, " files have NA conditions and will be excluded:")
    na_files <- temp_conds$file_name[is.na(temp_conds$condition)]
    message(paste(head(na_files, 5), collapse = ", "))
    if (length(na_files) > 5) message("... and ", length(na_files) - 5, " more")
    
    # Remove rows with NA conditions
    temp_conds <- temp_conds[!is.na(temp_conds$condition), ]
  }
  
  temp_conds
}, error = function(e) {
  message("Error reading conditions CSV file: ", e$message)
  message("File: ", CONFIG$conditions_csv)
  
  # Try to provide more specific help
  if (file.exists(CONFIG$conditions_csv)) {
    tryCatch({
      # Check if file is readable
      test_read <- read.csv(CONFIG$conditions_csv, nrows = 3)
      message("File exists and is readable. First few rows:")
      print(test_read)
      
      if (!"condition" %in% names(test_read)) {
        message("ERROR: Missing 'condition' column. Required columns: file_name, condition")
      }
      if (!"file_name" %in% names(test_read)) {
        message("ERROR: Missing 'file_name' column. Required columns: file_name, condition")
      }
      if (all(is.na(test_read$condition))) {
        message("ERROR: All conditions are NA. Please fill in actual condition values.")
      }
    }, error = function(e2) {
      message("File exists but cannot be read: ", e2$message)
    })
  } else {
    message("File does not exist: ", CONFIG$conditions_csv)
  }
  
  message("\nPlease check:")
  message("1. File exists: ", CONFIG$conditions_csv)
  message("2. File is a valid CSV with columns: file_name, condition")
  message("3. The 'condition' column is filled with '", CONFIG$condition_levels[1], "' or '", CONFIG$condition_levels[2], "'")
  message("4. You have read permissions to the file")
  quit(save = "no", status = 1)
})

if (!all(c("file_name","condition") %in% names(conds))) {
  die("`conditions_csv` must have columns: file_name, condition")
}

# ensure mapping covers all files
missing <- setdiff(files_tbl$file_name, conds$file_name)
if (length(missing) > 0) {
  die("These file_names are missing from the condition map: ",
      paste(missing, collapse = ", "))
}

# enforce two conditions
unique_conditions <- sort(unique(conds$condition))
if (length(unique_conditions) != 2) {
  die("Condition mapping must contain exactly two unique labels; found: ",
      paste(unique_conditions, collapse = ", "))
}

# set desired condition order
conds <- conds %>%
  mutate(condition = factor(condition, levels = CONFIG$condition_levels)) %>%
  arrange(file_name)

# =========================
# Interactive Column Selection
# =========================
message("Analyzing available columns...")

# Get column information from first file
first_file_info <- get_file_info(files_tbl$file[1])
all_columns <- first_file_info$columns

if (length(all_columns) == 0) {
  die("Could not read columns from first file: ", files_tbl$file[1])
}

message("Found ", length(all_columns), " columns in data files")

# Interactive selection of UMAP columns
message("Please select UMAP X coordinate column...")
umap_x_candidates <- grep("umap|UMAP|tSNE|tsne|PC1|pc1", all_columns, value = TRUE, ignore.case = TRUE)
if (length(umap_x_candidates) > 0) {
  message("Suggested UMAP X candidates: ", paste(head(umap_x_candidates, 5), collapse = ", "))
}

CONFIG$umap_x <- select_columns_gui(all_columns, "Select UMAP X coordinate", multiple = FALSE)
if (length(CONFIG$umap_x) == 0) die("No UMAP X column selected")

message("Please select UMAP Y coordinate column...")
umap_y_candidates <- grep("umap|UMAP|tSNE|tsne|PC2|pc2", all_columns, value = TRUE, ignore.case = TRUE)
if (length(umap_y_candidates) > 0) {
  message("Suggested UMAP Y candidates: ", paste(head(umap_y_candidates, 5), collapse = ", "))
}

CONFIG$umap_y <- select_columns_gui(all_columns, "Select UMAP Y coordinate", multiple = FALSE)
if (length(CONFIG$umap_y) == 0) die("No UMAP Y column selected")

# Interactive selection of markers
message("Please select markers of interest...")
marker_candidates <- grep("CD[0-9]|HLA|Ki67|IL|TNF|IFN|Foxp3", all_columns, value = TRUE, ignore.case = TRUE)
if (length(marker_candidates) > 0) {
  message("Suggested marker candidates: ", paste(head(marker_candidates, 10), collapse = ", "))
}

CONFIG$markers_of_interest <- select_columns_gui(all_columns, "Select Markers of Interest", multiple = TRUE)
if (length(CONFIG$markers_of_interest) == 0) {
  message("No markers selected - proceeding without marker analysis")
}

message("Selected UMAP X: ", CONFIG$umap_x)
message("Selected UMAP Y: ", CONFIG$umap_y)
message("Selected markers: ", paste(CONFIG$markers_of_interest, collapse = ", "))

# =========================
# Read + combine per-cell data with memory management
# =========================
message("Reading CSVs with aggressive memory management...")

read_one_optimized <- function(path, file_name) {
  file_info <- get_file_info(path)
  
  if (file_info$too_large && CONFIG$sample_large_files) {
    message("Large file detected (", file_info$n_rows, " rows, ", 
            round(file_info$size_mb, 1), " MB). Sampling...")
    
    # Sample rows instead of loading everything
    sample_size <- min(CONFIG$max_cells_per_file, file_info$n_rows)
    skip_rows <- sort(sample(2:(file_info$n_rows + 1), sample_size))
    
    # Read header first
    header <- fread(path, nrows = 0)
    
    # Read sampled rows
    dt <- fread(path, skip = min(skip_rows) - 1, nrows = length(skip_rows))
    names(dt) <- names(header)
    
    message("Sampled ", nrow(dt), " of ", file_info$n_rows, " cells")
    
  } else if (file_info$size_mb > CONFIG$max_memory_gb * 200) {
    # File is extremely large, use chunked reading
    dt <- read_csv_chunked(path)
  } else {
    # Normal reading
    dt <- fread(path)
  }
  
  # Convert to data.frame and add file identifier
  if (is.data.table(dt)) dt <- as.data.frame(dt)
  dt[["file_name"]] <- file_name
  
  # Select only essential columns early to save memory
  essential_cols <- c("file_name", CONFIG$umap_x, CONFIG$umap_y, CONFIG$markers_of_interest)
  cluster_cols_in_file <- grep(CONFIG$cluster_regex, names(dt), value = TRUE)
  essential_cols <- c(essential_cols, cluster_cols_in_file)
  
  # Keep only columns that exist
  essential_cols <- intersect(essential_cols, names(dt))
  dt <- dt[, essential_cols, drop = FALSE]
  
  return(dt)
}

# Read files with progress tracking and memory monitoring
pb <- progress_bar$new(
  format = "Loading files [:bar] :percent :current/:total :eta",
  total = nrow(files_tbl), clear = FALSE
)

cells_list <- vector("list", nrow(files_tbl))
memory_ok <- TRUE

for (i in seq_len(nrow(files_tbl))) {
  pb$tick()
  
  # Check memory before loading each file
  if (!check_memory()) {
    message("Memory pressure detected. Skipping remaining files...")
    memory_ok <- FALSE
    break
  }
  
  tryCatch({
    cells_list[[i]] <- read_one_optimized(files_tbl$file[i], files_tbl$file_name[i])
  }, error = function(e) {
    message("Error loading file ", basename(files_tbl$file[i]), ": ", e$message)
    cells_list[[i]] <- NULL
  })
  
  # Force garbage collection every 5 files
  if (i %% 5 == 0) {
    gc()
  }
}

# Remove NULL entries
cells_list <- cells_list[!sapply(cells_list, is.null)]

if (length(cells_list) == 0) {
  die("No files could be loaded successfully")
}

# Combine all data
message("Combining data from ", length(cells_list), " files...")

# Use rbindlist for faster binding with data.table
cells_dt <- rbindlist(cells_list, fill = TRUE)
cells <- as.data.frame(cells_dt)  # Convert back to data.frame for compatibility

rm(cells_list, cells_dt)  # Free memory immediately
gc()

# check columns exist
required_cols <- c(CONFIG$umap_x, CONFIG$umap_y)
missing_umap <- setdiff(required_cols, names(cells))
if (length(missing_umap) > 0) {
  message("Warning: Selected UMAP columns not found: ", paste(missing_umap, collapse = ", "))
  message("Available columns: ", paste(names(cells)[1:min(20, length(names(cells)))], collapse = ", "))
}

# markers present?
if (length(CONFIG$markers_of_interest) > 0) {
  present_markers <- intersect(CONFIG$markers_of_interest, names(cells))
  if (length(present_markers) == 0) {
    message("Warning: None of the selected markers were found in the data.")
    message("Available columns: ", paste(names(cells)[1:min(20, length(names(cells)))], collapse = ", "))
  } else if (length(present_markers) < length(CONFIG$markers_of_interest)) {
    message("Note: Some selected markers are missing. Using present subset: ",
            paste(present_markers, collapse = ", "))
  }
  CONFIG$markers_of_interest <- present_markers
}

# join conditions
cells <- cells %>%
  left_join(conds, by = "file_name")

if (any(is.na(cells$condition))) {
  die("Some rows have missing condition after join; check your mapping file_name values.")
}

message("Combined dataset: ", nrow(cells), " cells from ", length(unique(cells$file_name)), " files")
message("Conditions: ", paste(table(cells$condition), names(table(cells$condition)), sep = " ", collapse = ", "))

# =========================
# Identify and validate cluster columns
# =========================
cluster_cols <- names(cells)[str_detect(names(cells), CONFIG$cluster_regex)]
if (length(cluster_cols) == 0) {
  message("No cluster columns matched regex: ", CONFIG$cluster_regex)
  message("Available columns: ", paste(names(cells)[1:min(20, length(names(cells)))], collapse = ", "))
  die("Please check your cluster_regex pattern.")
}

message("Found ", length(cluster_cols), " cluster columns matching pattern '", CONFIG$cluster_regex, "'")

# Ensure binary 0/1 (coerce logical / numeric)
coerce01 <- function(x) {
  if (is.logical(x)) as.integer(x)
  else as.integer(round(as.numeric(x)))
}
cells[cluster_cols] <- lapply(cells[cluster_cols], coerce01)

# Validate cluster columns
for (col in cluster_cols) {
  unique_vals <- unique(cells[[col]])
  unique_vals <- unique_vals[!is.na(unique_vals)]
  if (!all(unique_vals %in% c(0, 1))) {
    message("Warning: Cluster column '", col, "' contains non-binary values: ", 
            paste(unique_vals, collapse = ", "))
  }
}

# Optionally drop clusters with near-zero membership
cluster_totals <- sapply(cluster_cols, function(col) sum(cells[[col]], na.rm = TRUE))
total_cells <- nrow(cells)
keep_clusters <- names(cluster_totals)[(cluster_totals / total_cells) >= CONFIG$min_cluster_frac]
dropped <- setdiff(cluster_cols, keep_clusters)
if (length(dropped) > 0) {
  message("Dropping ", length(dropped), " rare clusters (below min_cluster_frac): ",
          paste(head(dropped, 10), collapse = ", "), 
          if (length(dropped) > 10) " ..." else "")
}
cluster_cols <- keep_clusters
if (length(cluster_cols) == 0) die("All clusters were filtered out by min_cluster_frac.")

message("Analyzing ", length(cluster_cols), " clusters")

# =========================
# Check for Existing Results (Resume Functionality)
# =========================
results_file <- file.path(CONFIG$out_dir, "cluster_enrichment_results.csv")
if (file.exists(results_file)) {
  message("Found existing cluster enrichment results, loading...")
  results <- read_csv(results_file, show_col_types = FALSE)
  message("Loaded results for ", nrow(results), " clusters")
  message("Significant clusters: ", sum(results$significant, na.rm = TRUE))
  
  # Check if by_file_counts exists (optional, for additional analysis)
  by_file_counts_file <- file.path(CONFIG$out_dir, "by_file_cluster_counts.csv")
  if (file.exists(by_file_counts_file)) {
    message("Loading by-file counts...")
    by_file_counts <- read_csv(by_file_counts_file, show_col_types = FALSE)
  }
  
  message("Skipping statistical analysis, proceeding to visualizations...")
} else {
  message("No existing results found, performing fresh statistical analysis...")

# =========================
# Statistical Analysis
# =========================
message("Performing statistical analysis using method: ", CONFIG$stat_method)

make_outdir(CONFIG$out_dir)

if (CONFIG$stat_method == "pooled") {
  # Pooled Fisher's exact tests (faster, treats all cells independently)
  message("Running pooled Fisher's exact tests...")
  
  pb <- progress_bar$new(
    format = "Testing clusters [:bar] :percent :current/:total :eta",
    total = length(cluster_cols), clear = FALSE
  )
  
  results <- map_dfr(cluster_cols, function(cc) {
    pb$tick()
    stats <- do_fisher(cells[[cc]], cells$condition, levels(conds$condition))
    stats$cluster <- cc
    stats
  })
  
  # Calculate additional metrics
  results <- results %>%
    mutate(
      total_A = in_A + not_in_A,
      total_B = in_B + not_in_B,
      frac_A  = in_A / pmax(total_A, 1),
      frac_B  = in_B / pmax(total_B, 1),
      log2_enrichment_A_vs_B = log2((frac_A + 1e-9) / (frac_B + 1e-9)),
      p_adj = p.adjust(p_value, method = CONFIG$p_adjust_method),
      direction = case_when(
        log2_enrichment_A_vs_B > 0 ~ paste0(levels(conds$condition)[1], "↑"),
        log2_enrichment_A_vs_B < 0 ~ paste0(levels(conds$condition)[2], "↑"),
        TRUE ~ "tie"
      ),
      significant = p_adj < 0.05 & !is.na(p_adj)
    ) %>%
    select(cluster, starts_with("in_"), starts_with("not_in_"),
           starts_with("total_"), starts_with("frac_"),
           odds_ratio, ci_lower, ci_upper, log2_enrichment_A_vs_B, 
           p_value, p_adj, direction, significant, method) %>%
    arrange(p_adj, desc(abs(log2_enrichment_A_vs_B)))

} else if (CONFIG$stat_method == "mixed_effects") {
  # Mixed-effects models (accounts for file-level variation)
  message("Running mixed-effects models...")
  
  pb <- progress_bar$new(
    format = "Testing clusters [:bar] :percent :current/:total :eta",
    total = length(cluster_cols), clear = FALSE
  )
  
  me_results <- map_dfr(cluster_cols, function(cc) {
    pb$tick()
    stats <- do_mixed_effects(cc, cells)
    stats$cluster <- cc
    stats
  })
  
  # Calculate proportions for log2 enrichment
  cluster_props <- map_dfr(cluster_cols, function(cc) {
    cells %>%
      group_by(condition) %>%
      summarise(
        total = n(),
        in_cluster = sum(.data[[cc]], na.rm = TRUE),
        frac = in_cluster / total,
        .groups = "drop"
      ) %>%
      mutate(cluster = cc)
  })
  
  prop_wide <- cluster_props %>%
    select(cluster, condition, frac) %>%
    pivot_wider(names_from = condition, values_from = frac, 
                names_prefix = "frac_")
  
  # Combine results
  results <- me_results %>%
    left_join(prop_wide, by = "cluster") %>%
    mutate(
      log2_enrichment_A_vs_B = log2((get(paste0("frac_", CONFIG$condition_levels[1])) + 1e-9) / 
                                   (get(paste0("frac_", CONFIG$condition_levels[2])) + 1e-9)),
      p_adj = p.adjust(p_value, method = CONFIG$p_adjust_method),
      direction = case_when(
        log2_enrichment_A_vs_B > 0 ~ paste0(CONFIG$condition_levels[1], "↑"),
        log2_enrichment_A_vs_B < 0 ~ paste0(CONFIG$condition_levels[2], "↑"),
        TRUE ~ "tie"
      ),
      significant = p_adj < 0.05 & !is.na(p_adj)
    ) %>%
    arrange(p_adj, desc(abs(log2_enrichment_A_vs_B)))
}

# =========================
# Export Results
# =========================
write_csv(results, file.path(CONFIG$out_dir, "cluster_enrichment_results.csv"))

# Export cell-level data with cluster assignments for downstream analysis
message("Exporting cell-level cluster assignments...")
cell_data_export <- cells %>%
  mutate(cell_index = row_number()) %>%  # Add unique cell index
  select(
    cell_index,
    file_name, 
    condition,
    all_of(if(!any(missing_umap)) c(CONFIG$umap_x, CONFIG$umap_y) else character(0)),  # UMAP coordinates if available
    all_of(cluster_cols)  # All cluster assignments
  )

write_csv(cell_data_export, file.path(CONFIG$out_dir, "cell_cluster_assignments.csv"))
message("Cell-level data exported: ", nrow(cell_data_export), " cells with ", length(cluster_cols), " cluster assignments")

# Export per-file cluster counts for reference
by_file_counts <- cells %>%
  select(file_name, condition, all_of(cluster_cols)) %>%
  pivot_longer(all_of(cluster_cols), names_to = "cluster", values_to = "in_cluster") %>%
  group_by(file_name, condition, cluster) %>%
  summarise(n_in_cluster = sum(in_cluster, na.rm = TRUE),
            n_total = n(), .groups = "drop") %>%
  mutate(frac = n_in_cluster / n_total)

write_csv(by_file_counts, file.path(CONFIG$out_dir, "by_file_cluster_counts.csv"))

# Summary statistics
summary_stats <- results %>%
  summarise(
    total_clusters = n(),
    significant_clusters = sum(significant, na.rm = TRUE),
    mean_log2_enrichment = mean(abs(log2_enrichment_A_vs_B), na.rm = TRUE),
    median_p_value = median(p_value, na.rm = TRUE)
  )

write_csv(summary_stats, file.path(CONFIG$out_dir, "analysis_summary.csv"))

} # End of statistical analysis and export section

# =========================
# Visualizations
# =========================
message("Generating visualizations...")

# 1) Volcano plot: log2 enrichment vs -log10(p_adj)
volcano <- results %>%
  filter(!is.na(p_adj) & !is.na(log2_enrichment_A_vs_B)) %>%
  mutate(sig = significant) %>%
  ggplot(aes(x = log2_enrichment_A_vs_B, y = -log10(p_adj))) +
  geom_point(aes(color = sig, size = sig), alpha = 0.7) +
  geom_text_repel(
    data = . %>% slice_min(p_adj, n = 15),
    aes(label = cluster), 
    size = 3, max.overlaps = 20, box.padding = 0.3
  ) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = paste("Cluster Enrichment Analysis (", CONFIG$stat_method, ")", sep = ""),
    subtitle = paste0(CONFIG$condition_levels[1], " vs ", CONFIG$condition_levels[2]),
    x = paste0("log2 enrichment (", CONFIG$condition_levels[1], " / ", CONFIG$condition_levels[2], ")"),
    y = "-log10(adjusted p-value)",
    color = "Significant",
    size = "Significant"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(CONFIG$out_dir, "volcano_enrichment.png"), volcano, 
       width = 10, height = 8, dpi = 300)

# 2) Barplot of top significant clusters
n_top <- min(25, sum(results$significant, na.rm = TRUE))
if (n_top > 0) {
  bar_top <- results %>%
    filter(significant) %>%
    slice_min(p_adj, n = n_top) %>%
    mutate(cluster = fct_reorder(cluster, log2_enrichment_A_vs_B)) %>%
    ggplot(aes(x = cluster, y = log2_enrichment_A_vs_B, fill = direction)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = paste0("Top ", n_top, " significant clusters"),
      y = paste0("log2 enrichment (", CONFIG$condition_levels[1], " / ", CONFIG$condition_levels[2], ")"),
      x = "Cluster",
      fill = "Direction"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(CONFIG$out_dir, "top_clusters_bar.png"), bar_top, 
         width = 8, height = max(6, n_top * 0.3), dpi = 300)
}

# 3) UMAP visualization for top significant cluster
if (!any(missing_umap) && length(cluster_cols) > 0) {
  sig_clusters <- results %>% filter(significant) %>% arrange(p_adj)
  if (nrow(sig_clusters) > 0) {
    best_cluster <- sig_clusters$cluster[1]
    
    # Sample data for UMAP if too large
    plot_data <- cells
    if (nrow(plot_data) > 50000) {
      plot_data <- plot_data %>% slice_sample(n = 50000)
      message("Sampling 50,000 cells for UMAP visualization")
    }
    
    umap_df <- plot_data %>%
      select(all_of(c(CONFIG$umap_x, CONFIG$umap_y, "condition", best_cluster))) %>%
      rename(in_cluster = all_of(best_cluster)) %>%
      mutate(in_cluster = factor(in_cluster, levels = c(0,1), labels = c("No","Yes")))
    
    umap_plot <- ggplot(umap_df, aes_string(x = CONFIG$umap_x, y = CONFIG$umap_y)) +
      geom_point(aes(color = in_cluster), alpha = 0.6, size = 0.3) +
      scale_color_manual(values = c("No" = "lightgray", "Yes" = "red")) +
      facet_grid(in_cluster ~ condition) +
      coord_equal() +
      labs(
        title = paste0("UMAP for most significant cluster: ", best_cluster),
        subtitle = paste0("p_adj = ", scientific(sig_clusters$p_adj[1], digits = 3)),
        x = CONFIG$umap_x, y = CONFIG$umap_y,
        color = "In cluster"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
    
    ggsave(file.path(CONFIG$out_dir, paste0("umap_", make.names(best_cluster), ".png")),
           umap_plot, width = 10, height = 8, dpi = 300)
  }
}

# 4) Marker expression in best cluster (if markers available)
if (length(CONFIG$markers_of_interest) > 0 && length(cluster_cols) > 0) {
  sig_clusters <- results %>% filter(significant) %>% arrange(p_adj)
  if (nrow(sig_clusters) > 0) {
    best_cluster <- sig_clusters$cluster[1]
    
    # Sample data if too large
    plot_data <- cells
    if (nrow(plot_data) > 10000) {
      plot_data <- plot_data %>% slice_sample(n = 10000)
    }
    
    marker_long <- plot_data %>%
      select(all_of(c("condition", best_cluster, CONFIG$markers_of_interest))) %>%
      rename(in_cluster = all_of(best_cluster)) %>%
      mutate(in_cluster = factor(in_cluster, levels = c(0,1), labels = c("No","Yes"))) %>%
      pivot_longer(all_of(CONFIG$markers_of_interest), names_to = "marker", values_to = "value")
    
    marker_violin <- ggplot(marker_long, aes(x = in_cluster, y = value)) +
      geom_violin(aes(fill = condition), scale = "width", trim = TRUE, alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) +
      facet_wrap(~marker, scales = "free_y", ncol = 3) +
      scale_fill_brewer(type = "qual", palette = "Set2") +
      labs(
        title = paste0("Marker expression in/out of cluster: ", best_cluster),
        x = "In cluster?", y = "Marker intensity",
        fill = "Condition"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
    
    ggsave(file.path(CONFIG$out_dir, paste0("markers_", make.names(best_cluster), ".png")),
           marker_violin, width = 12, height = 8, dpi = 300)
  }
}

# 5) Summary heatmap of enrichment patterns
if (nrow(results) > 1) {
  # Select top clusters for heatmap
  n_heatmap <- min(50, nrow(results))
  heatmap_data <- results %>%
    slice_head(n = n_heatmap) %>%
    mutate(
      cluster = fct_reorder(cluster, log2_enrichment_A_vs_B),
      significance = case_when(
        significant ~ "Significant",
        TRUE ~ "Not significant"
      )
    )
  
  enrichment_heatmap <- ggplot(heatmap_data, aes(x = 1, y = cluster)) +
    geom_tile(aes(fill = log2_enrichment_A_vs_B), color = "white", size = 0.1) +
    geom_point(data = . %>% filter(significant), aes(x = 1.3), 
               shape = 8, size = 2, color = "black") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "log2\nenrichment"
    ) +
    labs(
      title = "Cluster Enrichment Summary",
      subtitle = paste0("★ indicates significance (adj. p < 0.05)"),
      y = "Cluster",
      x = ""
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave(file.path(CONFIG$out_dir, "enrichment_heatmap.png"), enrichment_heatmap,
         width = 6, height = max(6, n_heatmap * 0.2), dpi = 300)
}

# =========================
# Main Analysis with Error Handling
# =========================

# Wrap the entire analysis in error handling
tryCatch({
  
  message("=== Analysis Complete ===")
  message("Results saved to: ", CONFIG$out_dir)
  message("Files generated:")
  message("- cluster_enrichment_results.csv: Main statistical results")
  message("- by_file_cluster_counts.csv: Per-file cluster counts")
  message("- analysis_summary.csv: Summary statistics")
  message("- volcano_enrichment.png: Volcano plot")
  if (n_top > 0) message("- top_clusters_bar.png: Top significant clusters")
  if (!any(missing_umap)) message("- umap_*.png: UMAP visualization")
  if (length(CONFIG$markers_of_interest) > 0) message("- markers_*.png: Marker distributions")
  message("- enrichment_heatmap.png: Summary heatmap")

  # Print key findings
  sig_results <- results %>% filter(significant) %>% arrange(p_adj)
  message("\nKey findings (", nrow(sig_results), " significant clusters at adj. p < 0.05):")

  if (nrow(sig_results) > 0) {
    for (i in seq_len(min(10, nrow(sig_results)))) {
      row <- sig_results[i, ]
      direction_text <- if (row$log2_enrichment_A_vs_B > 0) {
        paste("enriched in", CONFIG$condition_levels[1])
      } else {
        paste("enriched in", CONFIG$condition_levels[2])
      }
      message(sprintf("  %d. %s: %.2f-fold %s (adj. p = %.2e)",
                     i, row$cluster, 2^abs(row$log2_enrichment_A_vs_B), 
                     direction_text, row$p_adj))
    }
    if (nrow(sig_results) > 10) {
      message("  ... and ", nrow(sig_results) - 10, " more significant clusters")
    }
  } else {
    message("  No statistically significant cluster enrichments found.")
  }

  message("\nAnalysis method: ", CONFIG$stat_method)
  message("Multiple testing correction: ", CONFIG$p_adjust_method)
  message("Total clusters analyzed: ", length(cluster_cols))
  message("Total cells: ", nrow(cells))
  message("Files: ", length(unique(cells$file_name)))
  
}, error = function(e) {
  message("\n=== ANALYSIS FAILED ===")
  message("Error: ", e$message)
  message("\nPossible solutions:")
  message("1. Reduce CONFIG$chunk_size (currently ", CONFIG$chunk_size, ")")
  message("2. Reduce CONFIG$max_memory_gb (currently ", CONFIG$max_memory_gb, ")")
  message("3. Set CONFIG$sample_large_files = TRUE")
  message("4. Reduce CONFIG$max_cells_per_file (currently ", CONFIG$max_cells_per_file, ")")
  message("5. Run analysis on a subset of files first")
  
  # Save partial results if they exist
  if (exists("results") && nrow(results) > 0) {
    make_outdir(CONFIG$out_dir)
    write_csv(results, file.path(CONFIG$out_dir, "partial_results.csv"))
    message("\nPartial results saved to: ", file.path(CONFIG$out_dir, "partial_results.csv"))
  }
  
  # Print memory info
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -h", intern = TRUE)
    message("\nCurrent memory status:")
    message(paste(mem_info, collapse = "\n"))
  }
  
  quit(status = 1)
})
