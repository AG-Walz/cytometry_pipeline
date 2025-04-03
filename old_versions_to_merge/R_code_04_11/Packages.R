### ======================= INSTALL & LOAD PACKAGES ======================== ###
# This section of the script is responsible for installing and loading all the necessary packages
# required for the analysis. It is divided into three main parts:
# 1. General packages: These are common R packages that can be installed from CRAN.
# 2. Bioconductor packages: These are specialized packages for bioinformatics that are available from Bioconductor.
# 3. Github repositories: These are packages that are hosted on GitHub and need to be installed using devtools.

# The script first checks if each package is already installed. If a package is not installed,
# it will be installed automatically. After ensuring all packages are installed, the script
# loads each package into the R session.

## =========================== General packages ============================= ##

# List of required packages that can be installed with install.packages()
packages = c("devtools","remotes", "BiocManager",
             "mclust","scales","tibble", "ggplot2",
             "viridis","plotly","data.table",   "ggrepel", 
             "rmarkdown", "RColorBrewer", "dplyr", "reshape2", "Rtsne", 
             "diptest", "LaplacesDemon",  "diptest",  
             "mousetrap", "readxl", "magrittr","ggridges", 'emdist', "plyr",
             "Matrix", "irlba", "RSpectra", "bigstatsr", "multimode","sBIC", 
             "tidyverse")
#

# List of general required packages that are not already installed
req_packages = packages[!(packages %in% installed.packages()[,"Package"])]
if(length(req_packages)){install.packages(req_packages, dependencies = TRUE, repos='http://cran.rstudio.com/')} # install the rest
invisible(sapply(packages, library, character.only = TRUE)) # Load packages

## ========================== Bioconductor packages ========================= ##

# Repeat the same procedure for Bioconductor packages
biopackages = c("openCyto","flowCore", "flowAI", "PeacoQC","flowWorkspace", 
                "uwot", "pheatmap","tidyr", "ggcyto", "flowVS",
                "preprocessCore", "flowTrans", "preprocessCore","harmony", 
                "Biobase","BiocGenerics")
# "FlowSOM", "ggpubr", "phyloseq",
req_biopackages = biopackages[!(biopackages %in% installed.packages()[,"Package"])]
if(length(req_biopackages)){BiocManager::install(req_biopackages)}
invisible(sapply(biopackages, library, character.only = TRUE))

## =========================== Github repositories ========================== ##
# setRepositories(ind = c(1:6, 8))
# Repeat the same procedure for Bioconductor packages
githubrepo = c("DillonHammill/CytoExploreRData","DillonHammill/CytoExploreR",
               "ImmuneDynamics/Spectre","biosurf/cyCombine", 
               "tidyverse/tidyverse","saeyslab/FlowSOM")
gitpackages = c("CytoExploreRData", "CytoExploreR", "Spectre", "cyCombine",
                "tidyverse", "FlowSOM")
req_github = githubrepo[!(githubrepo %in% installed.packages()[,"Package"])]
#if(length(req_github)){devtools::install_github(req_github, dependencies= TRUE, build_vignettes = TRUE)}
#invisible(sapply(gitpackages, library, character.only = TRUE))

# options(repos = BiocManager::repositories()) # required for publishing  


