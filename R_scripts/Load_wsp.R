# Required Libraries
library(CytoML)
library(flowWorkspace)
library(flowCore)

# -------------------- Parameter Definition --------------------
# Path to the FlowJo workspace file (.wsp)
wsp_file <- paste0("/home/ngr24/Desktop/Maerklin/cytometry_pipeline/",
                   "FlowJo tests/ICK-049-051.wsp")  # Replace with actual path

# Path to the directory containing associated FCS files
fcs_dir <- paste0("/home/ngr24/Desktop/Maerklin/cytometry_pipeline/",
                  "FlowJo tests")  # Replace with the actual path

out_dir <- "/home/ngr24/Desktop/Maerklin/cytometry_pipeline/FlowJo tests/"

# -------------------- Load Workspace --------------------
# Open the FlowJo workspace
flowjo_workspace <- open_flowjo_xml(wsp_file)

# List available sample groups in the FlowJo workspace
sample_groups <- fj_ws_get_sample_groups(flowjo_workspace)
cat("Available sample groups:\n")
print(sample_groups)

# Select a specific sample group (replace "YourGroupName" with the desired
# group name) selected_group <- "All Samples"  # Replace with the name of
# the group you want to use

# Check if the selected group exists
#if (!(selected_group %in% sample_groups$groupName)) {
#  stop(paste("The selected group", selected_group, "does not exist in the
# workspace."))
# }

# Parse the workspace into a GatingSet for the selected group
#gating_set <- flowjo_to_gatingset(flowjo_workspace, name = selected_group,
#                                  path = fcs_dir)

# Parse the workspace into a GatingSet
gating_set <- flowjo_to_gatingset(flowjo_workspace, path = fcs_dir)

# Get the sample names from the GatingSet
wsp_sample_names <- sampleNames(gating_set)

#Get the actual FCS file names in the directory
fcs_files <- list.files(path = fcs_dir, pattern = "\\.fcs$", full.names = TRUE)
fcs_file_names <- basename(fcs_files)

# Map workspace sample names to actual FCS file names
# Remove suffixes from workspace sample names to match FCS file names
mapped_sample_names <- sapply(wsp_sample_names, function(ws_name) {
  # Remove everything after the first underscore
  base_name <- sub("_.*$", "", ws_name)
  # Find the matching FCS file
  match <- grep(base_name, fcs_file_names, value = TRUE)
  if (length(match) == 1) {
    match
  } else {
    warning(paste("No matching FCS file found for workspace sample:", ws_name))
    NA
  }
})

# Remove samples with no matching FCS file
valid_samples <- !is.na(mapped_sample_names)
workspace_sample_names <- wsp_sample_names[valid_samples]
mapped_sample_names <- mapped_sample_names[valid_samples]


# Initialize a list to store classes for each FCS file
classes_list <- list()

# -------------------- Process Each Sample --------------------
for (i in seq_along(wsp_sample_names)) {
  # Get the gating hierarchy for the sample
  gating_hierarchy <- gating_set[[i]]

  # Extract the population names (classes) for the sample
  population_names <- gs_get_pop_paths(gating_set, path = "full")

  # Access the cytoframe for th ecurrent sample
  cytoframe <- load_cytoframe_from_fcs(fcs_files[[i]])

  # Convert the cytoframe to a data frame
  cytoframe_data <- exprs(cytoframe)

  # Get the number of events in the cytoframe
  num_events <- nrow(cytoframe_data)
  classes <- rep(NA, num_events)

  # Assign classes based on population data
  for (pop_name in population_names) {
    # Get the indices of events in the current population
    indices <- gh_pop_get_indices(gating_hierarchy, pop_name)
    classes[indices] <- pop_name
  }

  # Store the classes in the list
  classes_list[[workspace_sample_names[i]]] <- classes
}

# -------------------- Output --------------------
# Print the first few entries of the classes for each sample
for (sample in names(classes_list)) {
  cat("Sample:", sample, "\n")
  # Count occurrences of each class
  class_counts <- table(classes_list[[sample]])
  cat("Class counts:\n")
  print(class_counts)
}

# Save the classes list to an RDS file for later use
saveRDS(classes_list, file = paste0(out_dir, "classes_list.rds"))

# Print completion message
cat(paste0("Classes list saved to '", out_dir, "classes_list.rds'\n"))

# -------------------- Create DataFrames for Each Sample --------------------
# Initialize a list to store dataframes for each sample
dataframes_list <- list()

# Iterate through each sample
for (i in seq_along(workspace_sample_names)) {
  # Get the gating hierarchy for the sample
  gating_hierarchy <- gating_set[[i]]

  # Extract the population names (classes) for the sample
  population_names <- gs_get_pop_paths(gating_set, path = "full")

  # Access the cytoframe for the current sample
  cytoframe <- load_cytoframe_from_fcs(fcs_files[[i]])

  # Convert the cytoframe to a data frame
  cytoframe_data <- exprs(cytoframe)

  # Get the number of events (cells) in the cytoframe
  num_events <- nrow(cytoframe_data)

  # Initialize a matrix to store binary classification for each population
  population_matrix <- matrix(0, nrow = num_events,
                              ncol = length(population_names))
  colnames(population_matrix) <- population_names

  # Assign 1s for each population based on gating indices
  for (pop_name in population_names) {
    # Get the indices of events in the current population
    indices <- gh_pop_get_indices(gating_hierarchy, pop_name)
    population_matrix[indices, pop_name] <- 1
  }

  # Convert the matrix to a data frame
  population_df <- as.data.frame(population_matrix)

  # Add the data frame to the list
  dataframes_list[[workspace_sample_names[i]]] <- population_df
}

# -------------------- Output --------------------
# Save each dataframe to a CSV file
for (wsp_name in names(dataframes_list)) {
  output_file <- paste0(out_dir,
                        sub("\\.fcs$", "", wsp_name), "_populations.csv")
  write.csv(dataframes_list[[wsp_name]], file = output_file, row.names = FALSE)
  cat(paste0("Population dataframe saved to '", output_file, "'\n"))
}
