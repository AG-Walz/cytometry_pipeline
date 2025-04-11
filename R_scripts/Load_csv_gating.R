# Required Libraries: flowWorkspace, dplyr, openxlsx

input_csv_path <-"/data/Maerklin/ICK/cytolution_gated/files_ICK18_58_pre-gating_2_17_2025/"
input_fcs_path <- "/data/Maerklin/ICK/cytolution_gated/files_ICK18_58_pre-gating_2_17_2025/"
output_path <- paste0(input_fcs_path, "/Gated_fo/")
dir_create(output_path)

selected_gates <- c(
  "CD3...Manual.gating..1",
  "Monocytes..Manual.gating."
)

# Get list of FCS and csv files
fcs_files <- list.files(input_fcs_path, pattern = "\\.fcs$", full.names = TRUE)
csv_files <- list.files(input_csv_path, pattern = "\\.csv$", full.names = TRUE)

# Ensure files are in pairs by name, ignoring extensions
fcs_files_base <- gsub("\\.fcs$", "", basename(fcs_files))
csv_files_base <- gsub("\\.csv$", "", basename(csv_files))
if (!setequal(fcs_files_base, csv_files_base)) {
  stop("Mismatch between FCS and CSV file names.")
}

# Initialize a list to store results
results_list <- list()

# Iterate through the FCS and corresponding csv files
for (i in seq_along(fcs_files)) {

  fcs_file <- fcs_files[i]

  print(paste0("Processing file: ", basename(fcs_file)))

  # Find the corresponding csv file
  # TODO: add check if seems necessary
  csv_file <- csv_files[i]

  # Load the FCS file using flowWorkspace
  cytoframe <- load_cytoframe_from_fcs(fcs_file, transformation = FALSE)

  # Load the corresponding csv file
  gating_data <- read.csv(csv_file, sep = ";", header = TRUE)

  # Check dimensions of the cytoframe
  cytoframe_dim <- dim(cytoframe)
  csv_row_count <- nrow(gating_data)

  print(paste0("Cytoframe dimensions: ", paste(cytoframe_dim, collapse = " x ")))
  print(paste0("CSV row count:        ", csv_row_count))

  # Ensure that the number of rows in CSV matches the number of
  # events (rows) in the cytoframe
  if (csv_row_count != cytoframe_dim[1]) {
    stop("Mismatch between the number of rows in CSV and the number of events 
         in the FCS file.")
  }

  # Dynamically filter based on specified gating criteria
  for (j in seq_along(selected_gates)) {
    gate <- selected_gates[j]
    filtered_rows <- which(gating_data[[gate]] == 1)

    # Subset the cytoframe based on the row numbers from the CSV
    cytoframe_filtered <- cytoframe[filtered_rows, ]

    new_cytoframe_dim <- dim(cytoframe_filtered)

    print(paste0(gate, ": ",
                 paste(new_cytoframe_dim, collapse = " x "),
                " (",round(new_cytoframe_dim[1]/cytoframe_dim[1],3),"%)"))

    # Create a directory for the specific gate
    gate_output_path <- file.path(output_path, gate)
    dir_create(gate_output_path)

    # Save the filtered cytoframe to the respective gate directory
    flowCore::write.FCS(
      cytoframe_filtered,
      filename = file.path(gate_output_path,
                           paste0(fcs_name, "_", gate, "_gated.fcs"))
    )

    # Save the results to a list
    results_list[[paste0(fcs_name, "_", gate)]] <- data.frame(
      fcs_file = basename(fcs_file),
      gate = gate,
      num_cells_before = cytoframe_dim[1],
      num_cells_after = new_cytoframe_dim[1]
    )
  }
}

