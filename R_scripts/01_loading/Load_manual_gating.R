# Required Libraries
library(flowWorkspace) # for loading FCS files
library(dplyr)         # for data manipulation
library(openxlsx)      # for writing Excel files

# Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 

input_file_path <-"T:/User/Nacho/Data/Cytometry/ICK_data/CD3_gated_tf/"

dir.create(paste0(input_file_path, "/CD3_gated/"))
# Get list of FCS and csv files
fcs_files <- list.files(input_file_path, pattern = "\\.fcs$", full.names = TRUE)
csv_files <- list.files(input_file_path, pattern = "\\.csv$", full.names = TRUE)

# Ensure files are in pairs by name
fcs_files_base <- gsub("\\.fcs$", "", basename(fcs_files))
csv_files_base <- gsub("\\.csv$", "", basename(csv_files))

# Ensure files are in pairs by name, ignoring extensions
if (!setequal(fcs_files_base, csv_files_base)) {
  stop("Mismatch between FCS and CSV file names.")
}

# Initialize a list to store results
results_list <- list()

# Iterate through the FCS and corresponding csv files
for (i in seq_along(fcs_files)) {
  fcs_file <- fcs_files[i]
  fcs_name <- gsub("\\.fcs$", "", basename(fcs_file))
  
  print(paste0("Processing file: ", fcs_file))
  
  # Find the corresponding csv file
  csv_file <- csv_files[which(csv_files_base == fcs_name)]
  
  # Load the FCS file using flowWorkspace
  cytoframe <- load_cytoframe_from_fcs(fcs_file, transformation = FALSE)
  
  # Load the corresponding csv file
  gating_data <- read.csv(csv_file, sep = ";", header = TRUE)
  
  # Check dimensions of the cytoframe
  cytoframe_dim <- dim(cytoframe)
  csv_row_count <- nrow(gating_data)
  
  print(paste0("Cytoframe dimensions: ", paste(cytoframe_dim, collapse = " x ")))
  print(paste0("CSV row count: ", csv_row_count))
  
  # Ensure that the number of rows in CSV matches the number of events (rows) in the cytoframe
  if (csv_row_count != cytoframe_dim[1]) {
    stop("Mismatch between the number of rows in CSV and the number of events 
         in the FCS file.")
  }
  
  # Filter based on the 'CD3 (Manual gating)' column (using row numbers)
  cd3_gated_rows <- gating_data%>%
    dplyr::filter(gating_data$CD3..Manual.gating. == 1) %>%
    row_number()
  
  # Subset the cytoframe based on the row numbers from the CSV
  cytoframe_filtered <- cytoframe[cd3_gated_rows, ]
  
  new_cytoframe_dim <- dim(cytoframe_filtered)
  
  print(paste0("Gated_cytoframe dimensions: ", paste(new_cytoframe_dim, collapse = " x ")))
  
  # Save the results to a list
  results_list[[i]] <- data.frame(
    fcs_file = basename(fcs_file),
    num_cells_before = cytoframe_dim,
    num_cells_after = new_cytoframe_dim
  )
  
  # Save or process the subsetted cytoframe as needed
  # For example, save it to a new FCS file:
  flowCore::write.FCS(cytoframe_filtered, filename  = paste0(input_file_path, "/CD3_gated/", fcs_name, "_CD3_gated.fcs"))
}
# Combine all results into a single data frame
results_df <- bind_rows(results_list)

# Write the results to an Excel file
output_excel <- file.path(data_dir, "gating_results.xlsx")
write.xlsx(results_df, file = output_excel)

# Print message indicating completion
print(paste0("Results saved to: ", output_excel))