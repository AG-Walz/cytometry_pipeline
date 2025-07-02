# =========================== Packages ======================================

pkgs <- c("data.table", "lme4", "lmerTest", "emmeans",
          "future.apply", "progress", "tcltk", "stringr")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing))
  install.packages(missing, repos = "https://cloud.r-project.org")

invisible(lapply(pkgs, library, character.only = TRUE))

# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("Functions.R")

# Parallelisation ---------------------------------------------------------------
future::plan(multisession, workers = max(1, parallel::detectCores() - 1))

# =========================== Parameters ======================================

input_folder <- paste0("/data/Cytometry/ICK/ICK_T1_P+_poR/")
  # folder with CSVs
output_folder <- paste0(input_folder, "UMAP/")
dir_create(output_folder)

# Regex that identifies gating columns (comes from the original script)
gating_pat    <- "\\(Manual gating\\)$"

# Ensure folders exist ------------------------------------------------------
stopifnot(dir.exists(input_folder), dir.exists(output_folder))
files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0L)
message(sprintf("Found %d CSV files", length(files)))

# ---------------------------------------------------------------------------
# 1) Discover columns & interactive selection
# ---------------------------------------------------------------------------
colnames_all <- names(data.table::fread(files[1], nrows = 0))

gate_cols    <- grep(gating_pat, colnames_all, value = TRUE)
candidate_marker_cols <- setdiff(colnames_all, gate_cols)

# Debug: Check for NAs in candidate columns
message("Debugging column selection...")
message("Total candidate marker columns: ", length(candidate_marker_cols))
message("Any NAs in candidate_marker_cols: ", any(is.na(candidate_marker_cols)))
message("Any empty strings: ", any(candidate_marker_cols == ""))

# Remove any NA or empty values
candidate_marker_cols <- candidate_marker_cols[!is.na(candidate_marker_cols) & candidate_marker_cols != ""]
message("After cleaning - candidate marker columns: ", length(candidate_marker_cols))

# Interactive chooser for markers -------------------------------------------------
sel_markers <- tk_multiselect(candidate_marker_cols, "Select markers to analyse")
if (!length(sel_markers)) stop("No markers selected – aborting.")

message("Selected markers: ", length(sel_markers))
print(sel_markers)

# Debug gate columns
message("Total gate columns: ", length(gate_cols))
gate_display <- gsub(gating_pat, "", gate_cols)

# Interactive chooser for gates ---------------------------------------------------
sel_gates_display <- tk_multiselect(gate_display, "Select cell populations to analyse")
if (!length(sel_gates_display)) stop("No gates selected – aborting.")

message("Selected gate displays: ", length(sel_gates_display))
print(sel_gates_display)

# Map back to original gate names
sel_gates <- gate_cols[gate_display %in% sel_gates_display]
message("Selected gates (original names): ", length(sel_gates))
print(sel_gates)

# ---------------------------------------------------------------------------
# 2) Stream‑read only selected columns   +   append provenance column
# ---------------------------------------------------------------------------
need_cols <- c(sel_markers, sel_gates)
message(sprintf("\nReading %d files – keeping %d columns …", length(files), length(need_cols)))

# Debug: Check for issues with column names
message("Debugging column names...")
message("sel_markers length: ", length(sel_markers))
message("sel_gates length: ", length(sel_gates))
message("need_cols length: ", length(need_cols))

# Clean need_cols
need_cols <- need_cols[!is.na(need_cols) & need_cols != ""]
message("After cleaning need_cols length: ", length(need_cols))

# Read first file to check column availability
message("Checking first file columns...")
first_file_cols <- names(data.table::fread(files[1], nrows = 0))
missing_cols <- setdiff(need_cols, first_file_cols)
if (length(missing_cols) > 0) {
  message("WARNING: These columns are missing from the data:")
  print(missing_cols)
  need_cols <- intersect(need_cols, first_file_cols)
  message("Proceeding with available columns: ", length(need_cols))
}

expr_data <- data.table::rbindlist(
  lapply(files, function(f) {
    tryCatch({
      dt <- data.table::fread(f, select = need_cols, showProgress = FALSE)
      dt[, file := basename(f)]          # **provenance column**
      dt
    }, error = function(e) {
      message("Error reading file: ", f)
      message("Error: ", e$message)
      return(NULL)
    })
  }),
  use.names = TRUE, fill = FALSE
)

gc()  # free memory early

# ---------------------------------------------------------------------------
# 3) Convert gating columns → long format (population label per row)
# ---------------------------------------------------------------------------
expr_long_list <- list()

for (gate_idx in seq_along(sel_gates)) {
  gate_col <- sel_gates[gate_idx]
  gate_name <- sel_gates_display[gate_idx]
  
  # Filter only rows where this gate == 1
  gate_data <- expr_data[get(gate_col) == 1, 
                         c("file", sel_markers), 
                         with = FALSE]
  gate_data[, population := gate_name]
  
  expr_long_list[[gate_idx]] <- gate_data
}

expr_long <- data.table::rbindlist(expr_long_list)
rm(expr_long_list)
gc()
expr_df     <- expr_long            # pipeline downstream expects this name
marker_cols <- sel_markers

# ---------------------------------------------------------------------------
# 4) Aggregate to per‑file × population medians
# ---------------------------------------------------------------------------
message("Aggregating medians …")
agg_dt <- expr_df[, lapply(.SD, median, na.rm = TRUE),
                  by = .(file, population),
                  .SDcols = marker_cols]

# ---------------------------------------------------------------------------
# 5) Mixed‑effects model for every marker + log2FC + AUC
# ---------------------------------------------------------------------------
analyse_marker <- function(marker) {
  # Debug: Check marker input
  if (is.na(marker) || marker == "") {
    message("ERROR: Invalid marker input: ", marker)
    return(data.frame())
  }
  
  message("Analyzing marker: ", marker)
  
  # Create results list for each population vs total lymphocytes
  results_list <- list()
  
  # Get all populations
  all_populations <- unique(agg_dt$population)
  message("Populations found: ", paste(all_populations, collapse = ", "))
  
  for (pop in all_populations) {
    if (is.na(pop) || pop == "") {
      message("WARNING: Skipping invalid population: ", pop)
      next
    }
    
    tryCatch({
      # Get data for this specific population
      pop_data <- agg_dt[agg_dt$population == pop, ][[marker]]
      
      # Get data for ALL populations combined (total lymphocytes)
      all_lymph_data <- agg_dt[[marker]]
      
      if (length(pop_data) > 0 && length(all_lymph_data) > 0) {
        # Calculate statistics
        median_pop <- median(pop_data, na.rm = TRUE)
        median_lymph <- median(all_lymph_data, na.rm = TRUE)
        
        # Delta median
        delta_median <- median_pop - median_lymph
        
        # Log2 fold change
        log2FC <- log2((median_pop + 1e-9) / (median_lymph + 1e-9))
        
        # AUC using Wilcoxon test
        AUC_val <- NA
        p_value_val <- NA
        
        if (length(pop_data) > 1 && length(all_lymph_data) > 1) {
          tryCatch({
            w_test <- wilcox.test(pop_data, all_lymph_data, exact = FALSE)
            U <- as.numeric(w_test$statistic)
            n1 <- length(pop_data)
            n2 <- length(all_lymph_data)
            AUC_val <- U / (n1 * n2)
            p_value_val <- w_test$p.value
          }, error = function(e) {
            message("Wilcoxon test failed for ", pop, " - ", marker, ": ", e$message)
            AUC_val <- NA
            p_value_val <- NA
          })
        }
        
        # Create result row
        result_row <- data.frame(
          contrast = paste(pop, "- Total Lymphocytes"),
          estimate = delta_median,
          SE = NA,  # Not available from Wilcoxon test
          df = NA,  # Not available from Wilcoxon test
          t.ratio = NA,  # Not available from Wilcoxon test
          p.value = p_value_val,
          log2FC = log2FC,
          AUC = AUC_val,
          delta_median = delta_median,
          marker = marker,
          population = pop,
          stringsAsFactors = FALSE
        )
        
        results_list[[length(results_list) + 1]] <- result_row
      }
    }, error = function(e) {
      message("ERROR analyzing population ", pop, " for marker ", marker, ": ", e$message)
    })
  }
  
  # Combine all results for this marker
  if (length(results_list) > 0) {
    return(do.call(rbind, results_list))
  } else {
    message("No results generated for marker: ", marker)
    return(data.frame())
  }
}

message("Fitting mixed‑effects models …")
pb <- progress_bar$new(total = length(marker_cols),
                       format = "  :current/:total [:bar] :percent")

# Debug marker_cols before processing
message("Total markers to analyze: ", length(marker_cols))
message("Any NAs in marker_cols: ", any(is.na(marker_cols)))
message("Any empty strings in marker_cols: ", any(marker_cols == ""))

# Clean marker_cols
marker_cols <- marker_cols[!is.na(marker_cols) & marker_cols != ""]
message("After cleaning marker_cols: ", length(marker_cols))

if (length(marker_cols) == 0) {
  stop("No valid markers to analyze after cleaning!")
}

# Sequential processing with progress bar
results_list <- lapply(marker_cols, function(mk) {
  pb$tick()
  
  # Additional check per marker
  if (is.na(mk) || mk == "") {
    message("Skipping invalid marker: ", mk)
    return(data.frame())
  }
  
  # Check if marker exists in agg_dt
  if (!mk %in% names(agg_dt)) {
    message("Marker ", mk, " not found in aggregated data!")
    message("Available columns: ", paste(names(agg_dt), collapse = ", "))
    return(data.frame())
  }
  
  analyse_marker(mk)
})

results_dt <- data.table::rbindlist(results_list)

# ---------------------------------------------------------------------------
# 6) Additional analysis: Between-file comparisons
# ---------------------------------------------------------------------------
message("Analyzing between-file differences ...")

# Debug: Check data structure first
message("Debugging file variability analysis...")
message("Unique populations: ", paste(unique(agg_dt$population), collapse = ", "))
message("Unique files: ", paste(unique(agg_dt$file), collapse = ", "))
message("Total rows in agg_dt: ", nrow(agg_dt))

# Check sample sizes per population
pop_counts <- table(agg_dt$population)
message("Samples per population:")
print(pop_counts)

# For file-to-file variability, we need to analyze the structure differently
# Since agg_dt contains one median per file-population, we should analyze variation across files
file_variability_list <- list()

for (pop in unique(agg_dt$population)) {
  message("Processing population: ", pop)
  
  for (marker in marker_cols) {
    pop_data <- agg_dt[agg_dt$population == pop, ]
    
    message("  Marker: ", marker, " - Population: ", pop, " - Files: ", length(unique(pop_data$file)))
    
    if (length(unique(pop_data$file)) >= 3) {  # Need at least 3 different files
      marker_values <- pop_data[[marker]]
      file_names <- pop_data$file
      
      if (length(marker_values) > 2 && length(unique(marker_values)) > 1) {
        # Calculate file-to-file statistics
        file_means <- tapply(marker_values, file_names, mean, na.rm = TRUE)
        
        # Calculate Coefficient of Variation across files
        cv <- sd(file_means, na.rm = TRUE) / mean(file_means, na.rm = TRUE)
        
        # Calculate Intraclass Correlation using one-way ANOVA approach
        # This is appropriate when we have one observation per file-population
        tryCatch({
          # Create a simple ANOVA to estimate between-file variance
          aov_result <- aov(marker_values ~ factor(file_names))
          mse_between <- summary(aov_result)[[1]][1, "Mean Sq"]  # Between files
          mse_within <- summary(aov_result)[[1]][2, "Mean Sq"]   # Within files (error)
          
          # For balanced design (which we have), ICC calculation
          n_per_file <- 1  # We have 1 observation per file after aggregation
          icc <- (mse_between - mse_within) / (mse_between + (n_per_file - 1) * mse_within)
          icc <- max(0, icc)  # ICC cannot be negative
          
          file_variability_list[[length(file_variability_list) + 1]] <- data.frame(
            population = pop,
            marker = marker,
            file_variance = mse_between,
            residual_variance = mse_within,
            ICC = icc,
            CV = cv,
            n_files = length(unique(pop_data$file)),
            stringsAsFactors = FALSE
          )
          
          message("    SUCCESS: ICC = ", round(icc, 3), ", CV = ", round(cv, 3))
          
        }, error = function(e) {
          message("    ANOVA ERROR: ", e$message)
          
          # Fallback: Just CV calculation
          file_variability_list[[length(file_variability_list) + 1]] <<- data.frame(
            population = pop,
            marker = marker,
            file_variance = NA,
            residual_variance = NA,
            ICC = NA,
            CV = cv,
            n_files = length(unique(pop_data$file)),
            stringsAsFactors = FALSE
          )
          
          message("    FALLBACK: CV = ", round(cv, 3))
        })
      } else {
        message("    SKIPPED: No variation in marker values or insufficient data")
      }
    } else {
      message("    SKIPPED: Too few files (n=", length(unique(pop_data$file)), ")")
    }
  }
}

message("File variability analysis completed. Results: ", length(file_variability_list))

if (length(file_variability_list) > 0) {
  file_variability_dt <- data.table::rbindlist(file_variability_list, fill = TRUE)
} else {
  # Create empty data.table with correct structure
  file_variability_dt <- data.table::data.table(
    population = character(0),
    marker = character(0),
    file_variance = numeric(0),
    residual_variance = numeric(0),
    ICC = numeric(0),
    CV = numeric(0),
    n_files = integer(0)
  )
  message("WARNING: No file variability results generated!")
}

# Add effect size categories and interpretations
results_dt[, effect_size := fifelse(abs(log2FC) < 0.5, "small",
                           fifelse(abs(log2FC) < 1.0, "medium", "large"))]

results_dt[, AUC_interpretation := fifelse(is.na(AUC), "undefined",
                                  fifelse(AUC < 0.3 | AUC > 0.7, "strong",
                                  fifelse(AUC < 0.4 | AUC > 0.6, "moderate", "weak")))]

# Global BH FDR --------------------------------------------------------------
results_dt[, q_value := p.adjust(p.value, method = "BH")]

# ---------------------------------------------------------------------------
# 7) Save outputs
# ---------------------------------------------------------------------------
med_out  <- file.path(output_folder, "file_population_medians.csv")
stat_out <- file.path(output_folder, "mixedlm_results.csv")
file_var_out <- file.path(output_folder, "file_variability_analysis.csv")

data.table::fwrite(agg_dt, med_out)
data.table::fwrite(results_dt, stat_out)

# Only write file variability if we have results
if (nrow(file_variability_dt) > 0) {
  data.table::fwrite(file_variability_dt, file_var_out)
  message("✓ Saved file variability    →  ", file_var_out)
} else {
  message("No file variability results to save (insufficient data or model failures)")
}

message("✓ Saved per‑file medians     →  ", med_out)
message("✓ Saved mixed‑model stats   →  ", stat_out)

# Quick glance ----------------------------------------------------------------
message("\nTop hits (q ≤ 0.05):")
top_hits <- results_dt[q_value <= 0.05][order(q_value)]
print(top_hits[, .(contrast, marker, estimate, p.value, q_value, log2FC, AUC, delta_median, effect_size)])

# Summary statistics
message("\nEffect size summary (significant results):")
if (nrow(top_hits) > 0) {
  print(table(top_hits$effect_size))
  
  message("\nAUC interpretation summary (significant results):")
  print(table(top_hits$AUC_interpretation))
  
  message("\nFile variability summary:")
  if (nrow(file_variability_dt) > 0) {
    # Show ICC results if available
    if ("ICC" %in% names(file_variability_dt)) {
      high_icc <- file_variability_dt[!is.na(ICC) & ICC > 0.3][order(-ICC)]
      if (nrow(high_icc) > 0) {
        message("High ICC (> 0.3) results:")
        print(head(high_icc[, .(population, marker, ICC, CV, n_files)], 10))
      } else {
        message("No markers show high file-to-file variability (ICC > 0.3)")
      }
    }
    
    # Show CV results (always available now)
    if ("CV" %in% names(file_variability_dt)) {
      high_cv <- file_variability_dt[!is.na(CV) & CV > 0.2][order(-CV)]
      if (nrow(high_cv) > 0) {
        message("High coefficient of variation (> 0.2) results:")
        print(head(high_cv[, .(population, marker, CV, n_files)], 10))
      } else {
        message("All markers show low file-to-file variability (CV ≤ 0.2)")
      }
      
      # Overall CV statistics
      message("CV statistics across all markers:")
      cv_stats <- file_variability_dt[!is.na(CV), .(
        mean_CV = mean(CV, na.rm = TRUE),
        median_CV = median(CV, na.rm = TRUE),
        max_CV = max(CV, na.rm = TRUE),
        min_CV = min(CV, na.rm = TRUE)
      )]
      print(cv_stats)
    }
  } else {
    message("No file variability analysis results available")
  }
} else {
  message("No significant results found (q > 0.05)")
}