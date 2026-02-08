# =============================================================================
# config/Packages.R — Package Installation & Loading
# =============================================================================
# Installs and loads all packages required by the cytometry pipeline.
# Run once on a new environment, then source from pipeline scripts as needed.
# =============================================================================

## General CRAN packages
packages <- c(
  "devtools", "remotes", "BiocManager", "Cairo",
  "mclust", "scales", "tibble", "ggplot2", "openxlsx",
  "viridis", "plotly", "data.table", "ggrepel", "janitor",
  "rmarkdown", "RColorBrewer", "dplyr", "reshape2", "Rtsne",
  "diptest", "LaplacesDemon", "ggpubr", "units",
  "mousetrap", "readxl", "magrittr", "ggridges", "emdist", "plyr",
  "Matrix", "irlba", "RSpectra", "bigstatsr", "multimode", "sBIC",
  "tidyverse", "patchwork", "progress", "stringr", "lme4", "MASS",
  "broom", "purrr"
)

req_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(req_packages)) {
  install.packages(req_packages, dependencies = TRUE,
                   repos = "https://cloud.r-project.org")
}
invisible(sapply(packages, library, character.only = TRUE))

## Bioconductor packages
biopackages <- c(
  "openCyto", "flowCore", "flowAI", "PeacoQC", "flowWorkspace",
  "flowWorkspaceData", "uwot", "pheatmap", "tidyr", "ggcyto",
  "flowVS", "preprocessCore", "flowTrans",
  "harmony", "Biobase", "BiocGenerics", "cytolib"
)

req_biopackages <- biopackages[!(biopackages %in% installed.packages()[, "Package"])]
if (length(req_biopackages)) BiocManager::install(req_biopackages)
invisible(sapply(biopackages, library, character.only = TRUE))

## GitHub packages
githubrepo <- c(
  "Bioconductor/BiocParallel", "DillonHammill/openCyto",
  "DillonHammill/DataEditR", "JinmiaoChenLab/Rphenograph",
  "DillonHammill/HeatmapR", "maxsonBraunLab/CITE-Viz",
  "DillonHammill/CytoExploreRData", "DillonHammill/CytoExploreR",
  "ImmuneDynamics/Spectre", "biosurf/cyCombine",
  "tidyverse/tidyverse", "saeyslab/FlowSOM"
)
gitpackages <- c(
  "BiocParallel", "openCyto", "DataEditR", "Rphenograph",
  "HeatmapR", "CITEViz", "CytoExploreRData", "CytoExploreR",
  "Spectre", "cyCombine", "tidyverse", "FlowSOM"
)

req_github <- githubrepo[!(gitpackages %in% installed.packages()[, "Package"])]
if (length(req_github)) {
  devtools::install_github(req_github, dependencies = TRUE, build_vignettes = FALSE)
}
invisible(sapply(gitpackages, library, character.only = TRUE))
