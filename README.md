# Cytometry pipeline

## Overview
This pipeline is supposed to become a comprehensive tool designed to facilitate the analysis of cytometry data by preprocessing, analyzing, and visualizing cytometry data efficiently using open-source R-packages. 

## Getting started

### Prerequisites
- R(version 4.0 or higher)
- Rstudio (recommended, in some cases comment out rstudioapi and similar steps necessary)
- Proper rights and internet to download packages automatically using `Bioconductor::BiocManager::install()`, `remotes::install_github()` and `devtools::install_github()`

### Setup 

1. Edit the Parameter File:
    - Open parameter_template.R in your preferred text editor or RStudio.
    - Modify the parameters to suit your analysis needs. This includes setting paths to your data, specifying markers of interest, and other analysis options.
    - Save the modified file with a new name, e.g., `parameters_ICK_study.R`.
2. Edit the Main Script:
    - Open Main.R in your preferred text editor or RStudio.
    - Update the source line to point to your customized parameter file:
`source("path/to/my_parameters.R")`

### Running the pipeline
You have two options to run the pipeline:

1. Run Line by Line:
    - Open Main.R in RStudio or your preferred R environment.
    - Execute the script line by line, choosing which steps you wish to run.
2. Run the Entire Script:
    - Ensure that Main.R is properly configured with your parameter file.
    - Comment out any steps you wish to skip.
    - Run the entire script:
`source("Main.R")`
### Running standalone scripts
Some scripts are not yet included in the main workflow and can be run on their own, sometimes needing data produced by other scripts...
