# Helper script to install packages not available in standard Conda channels
# This is run AFTER the Conda environment is activated

options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- Install into the ACTIVE conda environment's library, never a global lib ---
# Installing into a shared ~/R/<version>-library would let packages built for a
# different R version leak into (and break) this env. Force the conda env lib.
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_prefix)) {
  env_lib <- file.path(conda_prefix, "lib", "R", "library")
  if (dir.exists(env_lib)) .libPaths(env_lib)
} else {
  warning("CONDA_PREFIX is not set — are you inside the 'cytometry_env' ",
          "environment? Packages may install to the wrong library.")
}
if (!startsWith(normalizePath(.libPaths()[1]), normalizePath(file.path(conda_prefix, "lib")))) {
  warning("Install target ", .libPaths()[1], " is outside the conda env. ",
          "Check for a global ~/.Rprofile / ~/.Renviron pinning .libPaths().")
}
cat("Installing packages to:", .libPaths()[1], "\n")

# 1. Install packages missing from Conda/CRAN
# (Certain packages like sBIC might be on CRAN but tricky in Conda)
pkg_cran <- c("sBIC")
new_cran <- pkg_cran[!(pkg_cran %in% installed.packages()[,"Package"])]
if(length(new_cran)) install.packages(new_cran)

# 1b. Install packages missing from Conda (Bioconductor)
if (!requireNamespace("harmony", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("harmony", update = FALSE, ask = FALSE)
}

# 2. Install GitHub Packages
# These are explicitly requested in your Packages.R and not in Conda

github_repos <- c(
  "DillonHammill/DataEditR",
  "JinmiaoChenLab/Rphenograph",
  "DillonHammill/HeatmapR", 
  "maxsonBraunLab/CITE-Viz",
  "DillonHammill/CytoExploreRData",
  "DillonHammill/CytoExploreR",
  "ImmuneDynamics/Spectre",
  "biosurf/cyCombine", 
  "tidyverse/tidyverse"
)

# Function to install if missing
install_if_missing_gh <- function(repo) {
  pkg_name <- basename(repo)
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat(sprintf("Installing %s from GitHub...\n", repo))
    tryCatch({
      devtools::install_github(repo, upgrade = "never", dependencies = TRUE)
      cat(sprintf("✓ %s installed.\n", pkg_name))
    }, error = function(e) {
      cat(sprintf("✗ Failed to install %s: %s\n", repo, e$message))
    })
  } else {
    cat(sprintf("• %s is already installed.\n", pkg_name))
  }
}

cat("=== Checking GitHub Packages ===\n")
invisible(lapply(github_repos, install_if_missing_gh))
cat("=== installation complete ===\n")
