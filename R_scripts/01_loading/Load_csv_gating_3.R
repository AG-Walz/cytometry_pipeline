# =============================================================================
# Load CSV Gating Data Script
# =============================================================================
# Description: Loads and processes CSV files exported from flow cytometry 
#              analysis software with gating information. Supports interactive
#              marker and gate selection.
#
# Dependencies: data.table, lme4, lmerTest, emmeans, future.apply, progress, 
#               tcltk, stringr, Functions.R
#
# TODOs:
#   - [ ] Add support for multiple gating column patterns (FlowJo, Cytobank, etc.)
#   - [ ] Implement automatic detection of gating format
#   - [ ] Add validation for column data types (numeric markers, boolean gates)
#   - [ ] Support for hierarchical gate structures
#   - [ ] Add progress bar for large file loading
#   - [ ] Implement memory-efficient chunked reading for very large files
#   - [ ] Add option to cache loaded data as RDS files
# =============================================================================

# =========================== Packages ======================================

pkgs <- c("data.table", "lme4", "lmerTest", "emmeans",
          "future.apply", "progress", "tcltk", "stringr")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing))
  install.packages(missing, repos = "https://cloud.r-project.org")

invisible(lapply(pkgs, library, character.only = TRUE))

# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Source from parent directory
source("../Functions.R")

# Parallelisation ---------------------------------------------------------------
future::plan(multisession, workers = max(1, parallel::detectCores() - 1))

# =========================== Parameters ======================================

input_folder <- paste0("/data/Cytometry/ICK/files_ICK18_80_pre-gating_11_18_2025/")
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

