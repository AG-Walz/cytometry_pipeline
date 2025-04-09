# ===================== CYTOMETRY ANALYSIS PIPELINE MAIN =============================
#
# Title: Cytometry Analysis Pipeline
# Author: Ignacio Garcia Ribelles
# Description: This script orchestrates the flow of the cytometry analysis
#              pipeline. It sets up the environment, loads necessary packages,
#              and executes various stages of data processing and analysis.
#
# Date Created: 2022-05-25
# Date Modified: 2024-11-11
# Version: 0.0
# Dependencies: Listed in the Packages.R file
# Usage: Edit the parameter_template.R, renaming it properly and setting the
# parameters before running the sections of this script as needed or commenting
# unnecessary steps.
#
# =============================== SETUP ========================================

# Rename the library path if needed
# .libPaths("Q:/IM/AG_Märklin/FSFC_Pipeline/library/")

# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Install and load packages
source("Packages.R", chdir = FALSE)

# Increase memory and use cache for larger files
knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)

# ================================= INPUT ======================================

parameter_file <- "./../test_data/Parameters_template.R"

source(parameter_file, chdir = TRUE)

# Load all files with the specified Loading parameters
source("Load_raw.R", chdir = TRUE)

# Check loaded files (including panel comparison)
source("Check.R", chdir = TRUE)

# Homogenize or downsample
source("Unify.R", chdir = TRUE)

# Downsample
source("Downsample.R", chdir = TRUE)

# plot density plots
source("Plot_raw.R", chdir = TRUE)

# ============================ PREPROCESSING ===================================

# Preprocess the files
## TODO: To improve, it is not running through
source("Preprocess.R", chdir = TRUE)

# Manual quality control to specify parameters more closely
# flow_iQC()

# ============================ ANALYSIS ========================================

# modality <- c(2,2,2,3,2,1,3,1,3)

# choose_antimode <- c(1,1,1,2,1,1,1,1)

# ============================ PLOTTING ========================================

# Plot the preprocessed files
source("UMAP.R", chdir = TRUE)
