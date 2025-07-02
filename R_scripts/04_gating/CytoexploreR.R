if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("DillonHammill/cytoexplorer")

library(cytoexplorer)

# Load your flow cytometry data (e.g., FCS files)
cyto_data <- cyto_load("path/to/your/fcs/files")

# Launch the interactive gating interface
cyto_gate(cyto_data)