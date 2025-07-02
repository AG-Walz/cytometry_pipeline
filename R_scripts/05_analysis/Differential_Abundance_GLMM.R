#!/usr/bin/env Rscript

# =========================
# Differential Abundance Analysis with Negative Binomial GLMMs
# Following diffcyt methodology (Weber et al.) for cytometry data
# =========================

# =========================
# CONFIG — EDIT THIS BLOCK
# =========================
CONFIG <- list(
  # 1) Folder with per-cell CSVs
    data_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR",
  
  # 2) File condition mapping CSV
  conditions_csv = "/data/Cytometry/ICK/ICK_T1_P+_poR/Output/file_conditions.csv",  # Path to conditions file
  
  # 3) Regex to pick cluster columns (binary 0/1 membership).
  #    Example: "^cluster_" or "(?i)meta.*cluster" or "^(?:C[0-9]+)$"
  cluster_regex = "Clustering annotated",
  
  # 4) Expression markers for characterization
  markers_of_interest = NULL,  # Will be selected interactively
  
  # 5) Two condition labels (reference level first)
  condition_levels = c("pCR", "npCR"),
  
  # 6) Output directory
  out_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/Output/GLMM_Results/",
  
  # 7) Minimum cells per cluster to retain
  min_cluster_size = 10,
  
  # 8) Multiple testing correction
  p_adjust_method = "BH",
  
  # 9) Random seed
  seed = 123,
  
  # 10) Memory management
  chunk_size = 25000,      # Smaller chunks for large datasets
  max_memory_gb = 2,       # Conservative memory limit
  max_cells_per_file = 100000,  # Much smaller per-file limit
  sample_large_files = TRUE,
  max_total_cells = 2000000,    # Maximum total cells to keep in memory
  progressive_processing = TRUE,  # Process files progressively to create abundance matrix
  
  # 11) Processing strategy
  process_in_batches = TRUE,     # Process files in small batches
  batch_size = 5,               # Number of files per batch
  save_intermediate = TRUE,      # Save intermediate results
  
  # 11) GLMM-specific parameters
  include_batch_effect = TRUE,  # Model file as random effect
  use_observation_level_re = TRUE,  # Include observation-level random effects for overdispersion
  count_offset = "log_total_cells",  # Offset for total cell count normalization
  
  # 12) Model selection and diagnostics
  model_family = "nb",  # "nb" (negative binomial), "poisson", or "auto"
  check_overdispersion = TRUE,  # Check for overdispersion and adjust model accordingly
  convergence_max_iter = 20000,
  
  # 13) Filtering parameters
  min_samples_per_cluster = 3,  # Minimum samples that must contain the cluster
  min_proportion_threshold = 0.0001  # 0.01% minimum proportion
)

# =========================
# Libraries
# =========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(stringr)
  library(lme4)          # For glmer.nb
  library(MASS)          # For glm.nb
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(viridis)
  library(progress)
  library(tcltk)
  
  # Optional packages - load if available
  tryCatch(library(broom.mixed), error = function(e) message("broom.mixed not available - using manual extraction"))
  tryCatch(library(emmeans), error = function(e) message("emmeans not available"))
  tryCatch(library(DHARMa), error = function(e) message("DHARMa not available - skipping advanced diagnostics"))
  tryCatch(library(performance), error = function(e) message("performance not available - using manual overdispersion check"))
})

set.seed(CONFIG$seed)

# =========================
# Helper Functions (same as original)
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

# Clean cluster names (from 01_split_clusters.R)
clean_cluster_name <- function(cluster_name) {
  if (is.na(cluster_name) || cluster_name == "") {
    return("Unknown")
  }
  
  # Remove "AI_Cluster_" prefix
  cleaned <- gsub("^AI_Cluster_", "", cluster_name)
  
  # Remove " (Clustering annotated)" suffix
  cleaned <- gsub(" \\(Clustering annotated\\)$", "", cleaned)
  
  return(cleaned)
}

# Extract cluster group from cluster name (from 01_split_clusters.R)
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
      if (grepl(pattern, cleaned_name, perl = TRUE)) {
        # Extract number if present
        if (grepl("_[0-9]+$", cleaned_name)) {
          number_part <- gsub(".*_([0-9]+)$", "\\1", cleaned_name)
          return(paste0(group_name, "_", number_part))
        } else {
          return(group_name)
        }
      }
    }
  }
  
  # If no pattern matches, return cleaned name
  return(cleaned_name)
}

# Convert binary cluster columns to single cluster assignment column
convert_binary_to_single_cluster <- function(cells_data, cluster_cols) {
  message("Converting ", length(cluster_cols), " binary cluster columns to single cluster assignment...")
  
  # Create a new column for cluster assignment
  cells_data$cluster_assignment <- NA_character_
  
  # For each row, find which cluster(s) the cell belongs to
  for (i in seq_len(nrow(cells_data))) {
    cluster_memberships <- character(0)
    
    for (col in cluster_cols) {
      if (!is.na(cells_data[[col]][i]) && cells_data[[col]][i] == 1) {
        # Clean and shorten the cluster name
        clean_name <- extract_cluster_group(col)
        cluster_memberships <- c(cluster_memberships, clean_name)
      }
    }
    
    if (length(cluster_memberships) == 0) {
      cells_data$cluster_assignment[i] <- "Unassigned"
    } else if (length(cluster_memberships) == 1) {
      cells_data$cluster_assignment[i] <- cluster_memberships[1]
    } else {
      # Handle multiple cluster assignments - take the first one or create composite name
      cells_data$cluster_assignment[i] <- cluster_memberships[1]  # Take first for simplicity
    }
  }
  
  return(cells_data)
}

# Memory-efficient CSV reading with chunking for large files
read_csv_chunked <- function(file_path, chunk_size = CONFIG$chunk_size) {
  file_size_mb <- file.info(file_path)$size / (1024^2)
  
  if (file_size_mb > CONFIG$max_memory_gb * 1000) {
    message("Large file detected (", round(file_size_mb, 1), " MB). Using chunked reading...")
    
    header <- fread(file_path, nrows = 0)
    
    total_rows <- tryCatch({
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
      tryCatch({
        length(readLines(file_path))
      }, error = function(e2) {
        message("Line counting failed, estimating from file size...")
        round(file.size(file_path) / 1000)
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

check_memory <- function() {
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -m", intern = TRUE)
    mem_line <- mem_info[2]
    mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
    used_pct <- mem_values[2] / mem_values[1]
    
    if (used_pct > 0.85) {
      message("Warning: High memory usage (", round(used_pct * 100, 1), "%). Running garbage collection...")
      gc()
      return(FALSE)
    }
  }
  return(TRUE)
}

# GUI for column selection
select_columns_gui <- function(column_names, title = "Select Columns", multiple = TRUE) {
  require(tcltk)
  
  if (length(column_names) == 0) {
    message("No columns available for selection")
    return(character(0))
  }
  
  tt <- tktoplevel()
  tkwm.title(tt, title)
  tkwm.geometry(tt, "500x400")
  
  main_frame <- tkframe(tt)
  tkpack(main_frame, fill = "both", expand = TRUE, padx = 10, pady = 10)
  
  instr_text <- if (multiple) "Select multiple columns (Ctrl+click):" else "Select one column:"
  instr_label <- tklabel(main_frame, text = instr_text, 
                         font = tkfont.create(size = 10, weight = "bold"))
  tkpack(instr_label, pady = 5)
  
  list_frame <- tkframe(main_frame)
  tkpack(list_frame, fill = "both", expand = TRUE)
  
  scrollbar <- tkscrollbar(list_frame, repeatinterval = 5,
                          command = function(...) tkyview(lb, ...))
  select_mode <- if (multiple) "extended" else "single"
  lb <- tklistbox(list_frame, yscrollcommand = function(...) tkset(scrollbar, ...),
                  selectmode = select_mode, height = 15)
  
  tkpack(scrollbar, side = "right", fill = "y")
  tkpack(lb, side = "left", fill = "both", expand = TRUE)
  
  for (col in column_names) {
    tkinsert(lb, "end", col)
  }
  
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
  
  tkfocus(lb)
  tkwait.window(tt)
  
  return(selected_items)
}

get_file_info <- function(file_path) {
  tryCatch({
    header <- fread(file_path, nrows = 0)
    
    row_count <- tryCatch({
      if (file.exists(file_path)) {
        wc_result <- system2("wc", args = c("-l", shQuote(file_path)), stdout = TRUE, stderr = FALSE)
        if (length(wc_result) > 0 && !is.na(wc_result)) {
          as.numeric(str_extract(wc_result, "^\\d+")) - 1
        } else {
          stop("wc command failed")
        }
      } else {
        stop("File does not exist")
      }
    }, error = function(e) {
      message("wc command failed for ", basename(file_path), ", using R fallback...")
      tryCatch({
        nrow(fread(file_path, select = 1))
      }, error = function(e2) {
        message("Row counting failed, estimating...")
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
# Progressive Data Processing Functions
# =========================

# Process single file and return cluster counts directly
process_file_to_counts <- function(path, file_name, cluster_regex, condition) {
  tryCatch({
    file_info <- get_file_info(path)
    
    # Determine sampling strategy
    if (file_info$too_large && CONFIG$sample_large_files) {
      message("Large file detected (", file_info$n_rows, " rows, ", 
              round(file_info$size_mb, 1), " MB). Sampling ", CONFIG$max_cells_per_file, " cells...")
      
      sample_size <- min(CONFIG$max_cells_per_file, file_info$n_rows)
      skip_rows <- sort(sample(2:(file_info$n_rows + 1), sample_size))
      
      header <- fread(path, nrows = 0)
      dt <- fread(path, skip = min(skip_rows) - 1, nrows = length(skip_rows))
      names(dt) <- names(header)
      
    } else if (file_info$size_mb > CONFIG$max_memory_gb * 200) {
      dt <- read_csv_chunked(path)
    } else {
      dt <- fread(path)
    }
    
    if (is.data.table(dt)) dt <- as.data.frame(dt)
    
    # Find cluster columns
    cluster_cols <- names(dt)[str_detect(names(dt), cluster_regex)]
    if (length(cluster_cols) == 0) {
      message("Warning: No cluster columns found in ", file_name)
      return(NULL)
    }
    
    # Coerce to binary 0/1
    coerce01 <- function(x) {
      if (is.logical(x)) as.integer(x)
      else as.integer(round(as.numeric(x)))
    }
    dt[cluster_cols] <- lapply(dt[cluster_cols], coerce01)
    
    # Create cluster assignment column
    cluster_assignment <- apply(dt[cluster_cols], 1, function(row) {
      which_clusters <- which(row == 1)
      if (length(which_clusters) == 0) {
        return("Unassigned")
      } else if (length(which_clusters) == 1) {
        return(cluster_cols[which_clusters])
      } else {
        # Multiple assignments - take first one or create combined name
        return(cluster_cols[which_clusters[1]])
      }
    })
    
    # Clean cluster names
    cluster_assignment <- clean_cluster_names(cluster_assignment)
    
    # Count clusters
    cluster_counts <- table(cluster_assignment)
    total_cells <- nrow(dt)
    
    # Return as tibble
    tibble(
      file_name = file_name,
      condition = condition,
      cluster_id = names(cluster_counts),
      count = as.numeric(cluster_counts),
      total_cells = total_cells,
      proportion = count / total_cells
    )
    
  }, error = function(e) {
    message("Error processing file ", file_name, ": ", e$message)
    return(NULL)
  })
}

# Process files in batches to avoid memory overload
process_files_progressive <- function(files_tbl, conds, cluster_regex) {
  message("Processing ", nrow(files_tbl), " files progressively...")
  
  # Join conditions to file table
  files_with_conds <- files_tbl %>%
    left_join(conds, by = "file_name")
  
  # Initialize results list
  abundance_parts <- list()
  
  # Process in batches
  n_batches <- ceiling(nrow(files_with_conds) / CONFIG$batch_size)
  
  pb <- progress_bar$new(
    format = "Processing batches [:bar] :percent :current/:total :eta",
    total = n_batches, clear = FALSE
  )
  
  for (batch_i in seq_len(n_batches)) {
    pb$tick()
    
    # Get batch indices
    start_idx <- (batch_i - 1) * CONFIG$batch_size + 1
    end_idx <- min(batch_i * CONFIG$batch_size, nrow(files_with_conds))
    batch_files <- files_with_conds[start_idx:end_idx, ]
    
    message("Processing batch ", batch_i, "/", n_batches, " (", nrow(batch_files), " files)")
    
    # Process each file in batch
    batch_results <- map2(batch_files$file, batch_files$file_name, function(path, fname) {
      condition <- batch_files$condition[batch_files$file_name == fname]
      process_file_to_counts(path, fname, cluster_regex, condition)
    })
    
    # Remove NULL results
    batch_results <- batch_results[!sapply(batch_results, is.null)]
    
    if (length(batch_results) > 0) {
      # Combine batch results
      batch_abundance <- bind_rows(batch_results)
      abundance_parts[[batch_i]] <- batch_abundance
      
      # Save intermediate results if configured
      if (CONFIG$save_intermediate) {
        write_csv(batch_abundance, file.path(CONFIG$out_dir, paste0("batch_", batch_i, "_abundance.csv")))
      }
    }
    
    # Force garbage collection after each batch
    gc()
    
    # Check memory
    if (!check_memory()) {
      message("Memory pressure detected after batch ", batch_i)
      gc()
    }
  }
  
  # Combine all batches
  message("Combining results from ", length(abundance_parts), " batches...")
  final_abundance <- bind_rows(abundance_parts)
  
  # Clean up intermediate files if they exist
  if (CONFIG$save_intermediate) {
    intermediate_files <- list.files(CONFIG$out_dir, pattern = "^batch_.*_abundance\\.csv$", full.names = TRUE)
    file.remove(intermediate_files)
  }
  
  return(final_abundance)
}

# =========================
# Memory-Efficient Data Processing
# =========================

# Note: The create_abundance_matrix is now handled by process_files_progressive

# Check for overdispersion
check_overdispersion <- function(count_data) {
  mean_count <- mean(count_data, na.rm = TRUE)
  var_count <- var(count_data, na.rm = TRUE)
  
  overdispersion_ratio <- var_count / mean_count
  
  # Rule of thumb: ratio > 1.5 suggests overdispersion
  is_overdispersed <- overdispersion_ratio > 1.5
  
  list(
    mean = mean_count,
    variance = var_count,
    ratio = overdispersion_ratio,
    overdispersed = is_overdispersed
  )
}

# Fit negative binomial GLMM for differential abundance
fit_nb_glmm <- function(cluster_data, cluster_id) {
  # Prepare data for this cluster
  model_data <- cluster_data %>%
    mutate(
      condition = factor(condition, levels = CONFIG$condition_levels),
      file_name = factor(file_name),
      obs_id = factor(1:n())  # Observation-level random effect
    )
  
  # Check if cluster has sufficient data
  total_count <- sum(model_data$count)
  if (total_count < CONFIG$min_cluster_size) {
    return(tibble(
      cluster_id = cluster_id,
      estimate = NA, std_error = NA, statistic = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      mean_count = 0, overdispersion_ratio = NA,
      n_samples = nrow(model_data),
      method = "insufficient_data"
    ))
  }
  
  # Check for variation in conditions
  if (length(unique(model_data$condition)) < 2) {
    return(tibble(
      cluster_id = cluster_id,
      estimate = NA, std_error = NA, statistic = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      mean_count = mean(model_data$count),
      overdispersion_ratio = NA,
      n_samples = nrow(model_data),
      method = "insufficient_conditions"
    ))
  }
  
  # Check if cluster appears in minimum number of samples
  samples_with_cluster <- sum(model_data$count > 0)
  if (samples_with_cluster < CONFIG$min_samples_per_cluster) {
    return(tibble(
      cluster_id = cluster_id,
      estimate = NA, std_error = NA, statistic = NA, p_value = NA,
      odds_ratio = NA, ci_lower = NA, ci_upper = NA,
      mean_count = mean(model_data$count),
      overdispersion_ratio = NA,
      n_samples = nrow(model_data),
      method = "too_sparse"
    ))
  }
  
  # Check for overdispersion
  overdispersion_check <- check_overdispersion(model_data$count)
  
  # Determine model family
  use_nb <- TRUE
  if (CONFIG$model_family == "poisson") {
    use_nb <- FALSE
  } else if (CONFIG$model_family == "auto") {
    use_nb <- overdispersion_check$overdispersed
  }
  
  # Build model formula
  if (CONFIG$include_batch_effect && CONFIG$use_observation_level_re) {
    # Full model: condition + random(file) + random(observation)
    formula_str <- "count ~ condition + offset(log_total_cells) + (1|file_name) + (1|obs_id)"
  } else if (CONFIG$include_batch_effect) {
    # Batch effect only
    formula_str <- "count ~ condition + offset(log_total_cells) + (1|file_name)"
  } else if (CONFIG$use_observation_level_re) {
    # Observation-level RE only
    formula_str <- "count ~ condition + offset(log_total_cells) + (1|obs_id)"
  } else {
    # Simple GLM (no random effects)
    formula_str <- "count ~ condition + offset(log_total_cells)"
  }
  
  tryCatch({
    # Fit appropriate model
    if (str_detect(formula_str, "\\|")) {
      # Mixed model
      if (use_nb) {
        model <- glmer.nb(as.formula(formula_str), 
                          data = model_data,
                          control = glmerControl(optimizer = "bobyqa", 
                                               optCtrl = list(maxfun = CONFIG$convergence_max_iter)))
      } else {
        model <- glmer(as.formula(formula_str), 
                       data = model_data,
                       family = poisson,
                       control = glmerControl(optimizer = "bobyqa", 
                                            optCtrl = list(maxfun = CONFIG$convergence_max_iter)))
      }
    } else {
      # Fixed effects only
      if (use_nb) {
        model <- glm.nb(as.formula(str_replace(formula_str, " \\+ offset\\(log_total_cells\\)", "")), 
                        data = model_data,
                        offset = model_data$log_total_cells)
      } else {
        model <- glm(as.formula(str_replace(formula_str, " \\+ offset\\(log_total_cells\\)", "")), 
                     data = model_data,
                     family = poisson,
                     offset = model_data$log_total_cells)
      }
    }
    
    # Extract results
    if (inherits(model, c("glmerMod", "glmerModLmerTest"))) {
      # Mixed model results
      fixed_effects <- fixef(model)
      se <- sqrt(diag(vcov(model)))
      
      # Get condition effect
      cond_effect_idx <- grep("condition", names(fixed_effects))
      if (length(cond_effect_idx) == 0) {
        return(tibble(
          cluster_id = cluster_id,
          estimate = NA, std_error = NA, statistic = NA, p_value = NA,
          odds_ratio = NA, ci_lower = NA, ci_upper = NA,
          mean_count = mean(model_data$count),
          overdispersion_ratio = overdispersion_check$ratio,
          n_samples = nrow(model_data),
          method = "no_condition_effect"
        ))
      }
      
      estimate <- fixed_effects[cond_effect_idx]
      std_error <- se[cond_effect_idx]
      statistic <- estimate / std_error
      p_value <- 2 * (1 - pnorm(abs(statistic)))
      
    } else {
      # GLM results
      coef_summary <- summary(model)$coefficients
      cond_effect_idx <- grep("condition", rownames(coef_summary))
      
      if (length(cond_effect_idx) == 0) {
        return(tibble(
          cluster_id = cluster_id,
          estimate = NA, std_error = NA, statistic = NA, p_value = NA,
          odds_ratio = NA, ci_lower = NA, ci_upper = NA,
          mean_count = mean(model_data$count),
          overdispersion_ratio = overdispersion_check$ratio,
          n_samples = nrow(model_data),
          method = "no_condition_effect"
        ))
      }
      
      estimate <- coef_summary[cond_effect_idx, "Estimate"]
      std_error <- coef_summary[cond_effect_idx, "Std. Error"]
      statistic <- coef_summary[cond_effect_idx, if(use_nb) "z value" else "z value"]
      p_value <- coef_summary[cond_effect_idx, if(use_nb) "Pr(>|z|)" else "Pr(>|z|)"]
    }
    
    # Calculate fold change and confidence interval
    fold_change <- exp(estimate)
    ci_lower <- exp(estimate - 1.96 * std_error)
    ci_upper <- exp(estimate + 1.96 * std_error)
    
    # Model diagnostics (optional)
    model_method <- paste0(
      if (str_detect(formula_str, "\\|")) "glmm_" else "glm_",
      if (use_nb) "nb" else "poisson"
    )
    
    tibble(
      cluster_id = cluster_id,
      estimate = estimate,
      std_error = std_error,
      statistic = statistic,
      p_value = p_value,
      odds_ratio = fold_change,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      mean_count = mean(model_data$count),
      overdispersion_ratio = overdispersion_check$ratio,
      n_samples = nrow(model_data),
      method = model_method
    )
    
  }, error = function(e) {
    message("Model failed for cluster ", cluster_id, ": ", e$message)
    
    # Fallback to Wilcoxon test on proportions
    tryCatch({
      prop_data <- model_data %>%
        mutate(proportion = count / total_cells) %>%
        select(condition, proportion)
      
      if (length(unique(prop_data$condition)) == 2) {
        wilcox_test <- wilcox.test(proportion ~ condition, data = prop_data)
        
        # Calculate median difference as estimate
        medians <- prop_data %>%
          group_by(condition) %>%
          summarise(median_prop = median(proportion), .groups = "drop")
        
        median_diff <- diff(medians$median_prop)
        
        tibble(
          cluster_id = cluster_id,
          estimate = median_diff,
          std_error = NA,
          statistic = as.numeric(wilcox_test$statistic),
          p_value = wilcox_test$p.value,
          odds_ratio = NA,
          ci_lower = NA,
          ci_upper = NA,
          mean_count = mean(model_data$count),
          overdispersion_ratio = overdispersion_check$ratio,
          n_samples = nrow(model_data),
          method = "wilcoxon_fallback"
        )
      } else {
        tibble(
          cluster_id = cluster_id,
          estimate = NA, std_error = NA, statistic = NA, p_value = NA,
          odds_ratio = NA, ci_lower = NA, ci_upper = NA,
          mean_count = mean(model_data$count),
          overdispersion_ratio = overdispersion_check$ratio,
          n_samples = nrow(model_data),
          method = "model_failed"
        )
      }
    }, error = function(e2) {
      tibble(
        cluster_id = cluster_id,
        estimate = NA, std_error = NA, statistic = NA, p_value = NA,
        odds_ratio = NA, ci_lower = NA, ci_upper = NA,
        mean_count = mean(model_data$count),
        overdispersion_ratio = overdispersion_check$ratio,
        n_samples = nrow(model_data),
        method = "all_failed"
      )
    })
  })
}

# =========================
# Condition File Handling (same as original)
# =========================
check_conditions_file_exists <- function() {
  message("Checking for existing conditions file...")
  
  if (!is.null(CONFIG$conditions_csv)) {
    message("✓ Conditions CSV already specified: ", CONFIG$conditions_csv)
    return(CONFIG$conditions_csv)
  }
  
  potential_paths <- c(
    file.path(CONFIG$out_dir, "file_conditions.csv"),
    file.path(getwd(), "file_conditions.csv"),
    file.path("~", "file_conditions.csv"),
    file.path(dirname(CONFIG$data_dir), "file_conditions.csv")
  )
  
  for (path in potential_paths) {
    if (file.exists(path)) {
      message("Found file: ", path)
      tryCatch({
        test_data <- read.csv(path, stringsAsFactors = FALSE, nrows = 5)
        if ("condition" %in% names(test_data) && "file_name" %in% names(test_data)) {
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

safe_scaffold_check <- function() {
  existing_file <- check_conditions_file_exists()
  
  if (!is.null(existing_file)) {
    message("Using existing conditions file: ", existing_file)
    return(TRUE)
  }
  
  message("\n", strrep("=", 50))
  message("SETUP REQUIRED: No conditions file found")
  message(strrep("=", 50))
  message("Please create a file called 'file_conditions.csv' with this structure:")
  message("")
  message("file_name,condition")
  message("your_file1,pCR")
  message("your_file2,npCR")
  message("...")
  message("")
  message("RECOMMENDED LOCATION:")
  message("Save to: ", file.path(CONFIG$out_dir, "file_conditions.csv"))
  message("")
  message("Then set: CONFIG$conditions_csv <- '", file.path(CONFIG$out_dir, "file_conditions.csv"), "'")
  message("And re-run this script.")
  message(strrep("=", 50))
  
  return(FALSE)
}

# =========================
# Data Loading and Processing
# =========================

# Load files with error handling
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

# Check if we can proceed
message("Step: Checking conditions file setup...")
can_proceed <- safe_scaffold_check()

if (!can_proceed) {
  message("Script stopping safely - please create conditions file as instructed above.")
  stop("Please create conditions file and re-run script. See instructions above.")
}

# Read and validate conditions file
conds <- tryCatch({
  temp_conds <- read_csv(CONFIG$conditions_csv, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(across(everything(), as.character))
  
  if (all(is.na(temp_conds$condition))) {
    stop("All conditions are NA. Please fill in the condition column with valid values.")
  }
  
  valid_conditions <- temp_conds$condition[!is.na(temp_conds$condition)]
  if (length(valid_conditions) == 0) {
    stop("No valid conditions found. Please fill in the condition column.")
  }
  
  na_count <- sum(is.na(temp_conds$condition))
  if (na_count > 0) {
    message("Warning: ", na_count, " files have NA conditions and will be excluded:")
    na_files <- temp_conds$file_name[is.na(temp_conds$condition)]
    message(paste(head(na_files, 5), collapse = ", "))
    if (length(na_files) > 5) message("... and ", length(na_files) - 5, " more")
    
    temp_conds <- temp_conds[!is.na(temp_conds$condition), ]
  }
  
  temp_conds
}, error = function(e) {
  message("Error reading conditions CSV file: ", e$message)
  message("File: ", CONFIG$conditions_csv)
  quit(save = "no", status = 1)
})

if (!all(c("file_name","condition") %in% names(conds))) {
  die("`conditions_csv` must have columns: file_name, condition")
}

missing <- setdiff(files_tbl$file_name, conds$file_name)
if (length(missing) > 0) {
  die("These file_names are missing from the condition map: ",
      paste(missing, collapse = ", "))
}

unique_conditions <- sort(unique(conds$condition))
if (length(unique_conditions) != 2) {
  die("Condition mapping must contain exactly two unique labels; found: ",
      paste(unique_conditions, collapse = ", "))
}

conds <- conds %>%
  mutate(condition = factor(condition, levels = CONFIG$condition_levels)) %>%
  arrange(file_name)

# =========================
# Interactive Column Selection
# =========================
message("Analyzing available columns...")

first_file_info <- get_file_info(files_tbl$file[1])
all_columns <- first_file_info$columns

if (length(all_columns) == 0) {
  die("Could not read columns from first file: ", files_tbl$file[1])
}

message("Found ", length(all_columns), " columns in data files")

# Interactive selection of markers
message("Please select markers of interest for characterization...")
marker_candidates <- grep("CD[0-9]|HLA|Ki67|IL|TNF|IFN|Foxp3", all_columns, value = TRUE, ignore.case = TRUE)
if (length(marker_candidates) > 0) {
  message("Suggested marker candidates: ", paste(head(marker_candidates, 10), collapse = ", "))
}

CONFIG$markers_of_interest <- select_columns_gui(all_columns, "Select Markers of Interest", multiple = TRUE)
if (length(CONFIG$markers_of_interest) == 0) {
  message("No markers selected - proceeding without marker characterization")
}

message("Selected markers: ", paste(CONFIG$markers_of_interest, collapse = ", "))

# =========================
# Progressive Data Processing (Memory-Efficient)
# =========================
message("Processing files progressively to avoid memory issues...")

make_outdir(CONFIG$out_dir)

# Check if we should use progressive processing
use_progressive <- CONFIG$progressive_processing || nrow(files_tbl) > 10

if (use_progressive) {
  message("Using progressive processing for ", nrow(files_tbl), " files...")
  
  # Get first file to determine cluster regex pattern
  first_file_info <- get_file_info(files_tbl$file[1])
  all_columns <- first_file_info$columns
  
  # Find cluster columns using regex
  cluster_cols_sample <- names(all_columns)[str_detect(names(all_columns), CONFIG$cluster_regex)]
  if (length(cluster_cols_sample) == 0) {
    die("No cluster columns found matching regex: ", CONFIG$cluster_regex)
  }
  
  message("Found ", length(cluster_cols_sample), " cluster columns matching pattern")
  message("Sample cluster columns: ", paste(head(cluster_cols_sample, 5), collapse = ", "))
  
  # Process files progressively to create abundance matrix directly
  abundance_data <- process_files_progressive(files_tbl, conds, CONFIG$cluster_regex)
  
  if (nrow(abundance_data) == 0) {
    die("No abundance data could be generated")
  }
  
  # Add log offset for GLMM
  abundance_data <- abundance_data %>%
    mutate(log_total_cells = log(total_cells))
  
  message("Generated abundance matrix: ", nrow(abundance_data), " cluster-sample combinations")
  
} else {
  # Original approach for smaller datasets
  message("Using original approach for small dataset...")
  
  # [Keep original code here if needed for small datasets]
  stop("Original approach not implemented in this version - use progressive processing")
}

# =========================
# Filter and Validate Abundance Data
# =========================

# Filter clusters with insufficient data
cluster_totals <- abundance_data %>%
  group_by(cluster_id) %>%
  summarise(
    total_count = sum(count), 
    n_samples_with_cluster = sum(count > 0),
    max_proportion = max(proportion),
    .groups = "drop"
  )

keep_clusters <- cluster_totals %>%
  filter(
    total_count >= CONFIG$min_cluster_size,
    n_samples_with_cluster >= CONFIG$min_samples_per_cluster,
    max_proportion >= CONFIG$min_proportion_threshold
  ) %>%
  pull(cluster_id)

dropped_clusters <- setdiff(cluster_totals$cluster_id, keep_clusters)
if (length(dropped_clusters) > 0) {
  message("Dropping ", length(dropped_clusters), " clusters due to insufficient data:")
  message("- Total count < ", CONFIG$min_cluster_size, ": ", 
          sum(cluster_totals$total_count < CONFIG$min_cluster_size))
  message("- Samples with cluster < ", CONFIG$min_samples_per_cluster, ": ",
          sum(cluster_totals$n_samples_with_cluster < CONFIG$min_samples_per_cluster))
  message("- Max proportion < ", CONFIG$min_proportion_threshold, ": ",
          sum(cluster_totals$max_proportion < CONFIG$min_proportion_threshold))
}

abundance_data <- abundance_data %>%
  filter(cluster_id %in% keep_clusters)

message("Analyzing ", length(keep_clusters), " clusters with sufficient data")

# Run differential abundance testing
message("Running negative binomial GLMM analysis...")
message("Model configuration:")
message("- Family: ", CONFIG$model_family)
message("- Batch effects: ", CONFIG$include_batch_effect)
message("- Observation-level RE: ", CONFIG$use_observation_level_re)
message("- Check overdispersion: ", CONFIG$check_overdispersion)

pb <- progress_bar$new(
  format = "Testing clusters [:bar] :percent :current/:total :eta",
  total = length(keep_clusters), clear = FALSE
)

results_list <- map(keep_clusters, function(cluster) {
  pb$tick()
  
  cluster_data <- abundance_data %>%
    filter(cluster_id == cluster)
  
  fit_nb_glmm(cluster_data, cluster)
})

# Combine results
results <- bind_rows(results_list)

# Calculate additional metrics
results <- results %>%
  mutate(
    p_adj = p.adjust(p_value, method = CONFIG$p_adjust_method),
    log2_fold_change = log2(odds_ratio),
    direction = case_when(
      is.na(estimate) ~ "No test",
      estimate > 0 ~ paste0(CONFIG$condition_levels[2], "↑"),
      estimate < 0 ~ paste0(CONFIG$condition_levels[1], "↑"),
      TRUE ~ "tie"
    ),
    significant = p_adj < 0.05 & !is.na(p_adj),
    effect_size = case_when(
      abs(log2_fold_change) >= 1 ~ "Large",
      abs(log2_fold_change) >= 0.5 ~ "Medium", 
      abs(log2_fold_change) >= 0.2 ~ "Small",
      TRUE ~ "Negligible"
    )
  ) %>%
  arrange(p_adj, desc(abs(log2_fold_change)))

# =========================
# Export Results
# =========================
write_csv(results, file.path(CONFIG$out_dir, "differential_abundance_results.csv"))

# Export abundance matrix
write_csv(abundance_data, file.path(CONFIG$out_dir, "cluster_abundance_matrix.csv"))

# Export cluster summary
cluster_summary <- abundance_data %>%
  group_by(cluster_id, condition) %>%
  summarise(
    n_samples = n(),
    total_count = sum(count),
    mean_count = mean(count),
    median_count = median(count),
    mean_proportion = mean(proportion),
    median_proportion = median(proportion),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(n_samples, total_count, mean_count, median_count, mean_proportion, median_proportion),
    names_sep = "_"
  )

write_csv(cluster_summary, file.path(CONFIG$out_dir, "cluster_summary_by_condition.csv"))

# Export model diagnostics
model_diagnostics <- results %>%
  count(method, name = "n_clusters") %>%
  mutate(proportion = n_clusters / sum(n_clusters))

write_csv(model_diagnostics, file.path(CONFIG$out_dir, "model_diagnostics.csv"))

# Summary statistics
summary_stats <- results %>%
  summarise(
    total_clusters_tested = n(),
    successful_tests = sum(!is.na(p_value)),
    significant_clusters = sum(significant, na.rm = TRUE),
    mean_overdispersion_ratio = mean(overdispersion_ratio, na.rm = TRUE),
    median_p_value = median(p_value, na.rm = TRUE),
    mean_abs_log2fc = mean(abs(log2_fold_change), na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_stats, file.path(CONFIG$out_dir, "analysis_summary.csv"))

# =========================
# Visualizations
# =========================
message("Generating visualizations...")

# 1) Volcano plot
volcano_data <- results %>%
  filter(!is.na(p_adj) & !is.na(log2_fold_change))

if (nrow(volcano_data) > 0) {
  volcano <- volcano_data %>%
    mutate(sig = significant) %>%
    ggplot(aes(x = log2_fold_change, y = -log10(p_adj))) +
    geom_point(aes(color = sig, size = sig), alpha = 0.7) +
    geom_text_repel(
      data = . %>% filter(significant) %>% slice_min(p_adj, n = 15),
      aes(label = cluster_id), 
      size = 3, max.overlaps = 20, box.padding = 0.3
    ) +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = "Differential Abundance Analysis (Negative Binomial GLMM)",
      subtitle = paste0(CONFIG$condition_levels[1], " vs ", CONFIG$condition_levels[2]),
      x = paste0("log2 fold-change (", CONFIG$condition_levels[2], " / ", CONFIG$condition_levels[1], ")"),
      y = "-log10(adjusted p-value)",
      color = "Significant",
      size = "Significant"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(CONFIG$out_dir, "volcano_plot.png"), volcano, 
         width = 10, height = 8, dpi = 300)
}

# 2) Bar plot of significant results
sig_results <- results %>% filter(significant)
n_sig <- nrow(sig_results)

if (n_sig > 0) {
  n_top <- min(25, n_sig)
  
  bar_plot <- sig_results %>%
    slice_min(p_adj, n = n_top) %>%
    mutate(cluster_id = fct_reorder(cluster_id, log2_fold_change)) %>%
    ggplot(aes(x = cluster_id, y = log2_fold_change, fill = direction)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = paste0("Top ", n_top, " significant clusters"),
      y = paste0("log2 fold-change (", CONFIG$condition_levels[2], " / ", CONFIG$condition_levels[1], ")"),
      x = "Cluster ID",
      fill = "Direction"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(CONFIG$out_dir, "significant_clusters_bar.png"), bar_plot, 
         width = 8, height = max(6, n_top * 0.3), dpi = 300)
}

# 3) Abundance heatmap for significant clusters
if (n_sig > 0) {
  n_heatmap <- min(30, n_sig)
  
  heatmap_clusters <- sig_results %>%
    slice_min(p_adj, n = n_heatmap) %>%
    pull(cluster_id)
  
  heatmap_data <- abundance_data %>%
    filter(cluster_id %in% heatmap_clusters) %>%
    group_by(cluster_id, condition) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop")
  
  abundance_heatmap <- heatmap_data %>%
    ggplot(aes(x = condition, y = cluster_id, fill = mean_proportion)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_viridis_c(name = "Mean\nProportion") +
    labs(
      title = "Cluster Abundance by Condition",
      subtitle = paste0("Top ", n_heatmap, " significant clusters"),
      x = "Condition",
      y = "Cluster ID"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  ggsave(file.path(CONFIG$out_dir, "abundance_heatmap.png"), abundance_heatmap,
         width = 6, height = max(6, n_heatmap * 0.2), dpi = 300)
}

# 4) Model diagnostics plots
if (nrow(results) > 0) {
  # Overdispersion vs count relationship
  overdispersion_plot <- results %>%
    filter(!is.na(overdispersion_ratio) & !is.na(mean_count)) %>%
    ggplot(aes(x = mean_count, y = overdispersion_ratio)) +
    geom_point(aes(color = significant), alpha = 0.6) +
    geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
    labs(
      title = "Overdispersion vs Mean Count",
      subtitle = "Red line shows overdispersion threshold (1.5)",
      x = "Mean Count (log10)",
      y = "Overdispersion Ratio (log10)",
      color = "Significant"
    ) +
    theme_minimal()
  
  ggsave(file.path(CONFIG$out_dir, "overdispersion_diagnostic.png"), overdispersion_plot,
         width = 8, height = 6, dpi = 300)
  
  # Method usage summary
  method_plot <- model_diagnostics %>%
    ggplot(aes(x = reorder(method, n_clusters), y = n_clusters)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste0(n_clusters, "\n(", round(proportion * 100, 1), "%)")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Model Methods Used",
      x = "Method",
      y = "Number of Clusters"
    ) +
    theme_minimal()
  
  ggsave(file.path(CONFIG$out_dir, "method_usage.png"), method_plot,
         width = 8, height = 6, dpi = 300)
}

# =========================
# Final Summary
# =========================
message("=== Differential Abundance Analysis Complete ===")
message("Results saved to: ", CONFIG$out_dir)
message("Files generated:")
message("- differential_abundance_results.csv: Main statistical results")
message("- cluster_abundance_matrix.csv: Cluster-sample count matrix")
message("- cluster_summary_by_condition.csv: Summary statistics by condition")
message("- model_diagnostics.csv: Model fitting diagnostics")
message("- analysis_summary.csv: Overall summary statistics")
message("- volcano_plot.png: Volcano plot of results")
if (n_sig > 0) {
  message("- significant_clusters_bar.png: Bar plot of significant clusters")
  message("- abundance_heatmap.png: Abundance heatmap")
}
message("- overdispersion_diagnostic.png: Overdispersion diagnostic plot")
message("- method_usage.png: Model method usage summary")

message("\nKey findings (", n_sig, " significant clusters at adj. p < 0.05):")

if (n_sig > 0) {
  for (i in seq_len(min(10, n_sig))) {
    row <- sig_results[i, ]
    direction_text <- if (row$log2_fold_change > 0) {
      paste("enriched in", CONFIG$condition_levels[2])
    } else {
      paste("enriched in", CONFIG$condition_levels[1])
    }
    message(sprintf("  %d. Cluster %s: %.2f-fold %s (adj. p = %.2e, effect: %s)",
                   i, row$cluster_id, 2^abs(row$log2_fold_change), 
                   direction_text, row$p_adj, row$effect_size))
  }
  if (n_sig > 10) {
    message("  ... and ", n_sig - 10, " more significant clusters")
  }
} else {
  message("  No statistically significant differential abundance found.")
}

message("\nModel details:")
message("- Method: Negative binomial GLMM")
message("- Random effects: ", if (CONFIG$include_batch_effect) "File (batch)" else "None", 
        if (CONFIG$use_observation_level_re) " + Observation-level" else "")
message("- Offset: Total cells per sample (log scale)")
message("- Multiple testing correction: ", CONFIG$p_adjust_method)
message("- Total clusters tested: ", length(keep_clusters))
message("- Successful model fits: ", sum(!is.na(results$p_value)))
message("- Mean overdispersion ratio: ", round(summary_stats$mean_overdispersion_ratio, 2))

message("\nThis analysis follows the diffcyt methodology for differential abundance testing in cytometry data.")
message("For more details, see Weber et al. (2019) 'diffcyt: Differential discovery in high-dimensional cytometry'")