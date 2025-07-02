#!/usr/bin/env Rscript
# ================================================================
# Analyse single-cell expression markers vs. gating groups
# ================================================================

# Rename the library path if needed
#.libPaths("Q:/IM/AG_Märklin/AG Märklin_Software/library/")

# Required packages
pkgs <- c("data.table", "lme4", "lmerTest", "emmeans", "stringr",
          "future.apply", "progress", "dplyr", "tcltk", "readr", "purrr",
          "pbapply", "matrixTests", "matrixStats", "R.utils")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install,
                                         repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Parallelism (set to available cores – 1 to keep UI responsive)
plan(multisession, workers = max(1, parallel::detectCores() - 1))

# ========================= Parameters =====================================
input_folder <- "/data/Cytometry/ICK/ICK_T1_PR/files_ICK_T1_P+_poR_pre-gating_9_15_2025/"
output_folder <- paste0(input_folder, "Results/")
gating_pat  <- "\\(Manual gating\\)$"

# ========================== Source packages ==================================
source("Packages.R")
source("Functions.R")
dir_create(output_folder)


# ========================== Loading data ======================================

files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
message(sprintf("Reading %d files ...", length(files)))
colnames_all <- names(data.table::fread(files[1], nrows = 0)) 

# ----------------------------------------------------------------
# 2. Detect marker & gating columns
# ----------------------------------------------------------------

gate_cols   <- grep(gating_pat, colnames_all, value = TRUE) 
candidate_marker_cols <- setdiff(colnames_all, gate_cols) 

# ----------------------------------------------------------------
# 3. Interactive column selection
# ----------------------------------------------------------------
# MARKER SELECTION
sel_markers <- tk_multiselect(candidate_marker_cols, "Select markers to analyze")
if (!length(sel_markers)) stop("No markers selected – aborting.")

# GATE SELECTION 
gate_display <- gsub(gating_pat, "", gate_cols)
sel_gates_display <- tk_multiselect(gate_display, "Select cell populations to analyze")
if (!length(sel_gates_display)) stop("No gates selected – aborting.")

# Map back to original gate names
sel_gates <- gate_cols[gate_display %in% sel_gates_display]

# ----------------------------------------------------------------
# 2.  Stream-read only the columns the user wants  (⭑ key change ⭑)
# ----------------------------------------------------------------
need_cols <- c(sel_markers, sel_gates)        # nothing else!
message(sprintf("\nReading %d files – keeping %d columns …", length(files), length(need_cols)))

expr_data <- data.table::rbindlist(
  lapply(files, function(f) data.table::fread(f, select = need_cols, showProgress = FALSE)),
  use.names = TRUE, fill = FALSE
)

gc()  # prompt garbage collection

