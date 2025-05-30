# -------------------- Required Libraries --------------------
library(CytoML)
library(flowWorkspace)
library(flowCore)
library(ggplot2)
library(dplyr)

# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# -------------------- Parameter Definition --------------------
# Path to the FlowJo workspace file (.wsp)
# wsp_files <- file.choose()
# List workspace files
wsp_dir <- "./../Gating_comparisons/FlowJo_gatings/"

# Custom name for the workspace (keep identifier as short as possible)
custom_name <- c("Lea", "Jonas", "Samu2")

# Path to the directory containing associated FCS files
fcs_dir <- "./../Gating_comparisons/raw_data/" # Replace with actual path

# Turn plotting barplots on or off
# TODO: change T/F for choose which type of plotting you want to do inkl None
plot_barplot <- FALSE

# Output directory
out_dir <- "./../Output/"

# ------------------------------- Functions -------------------------------
source("Functions.R")
# ------------------------------ Load Workspace --------------------------------
wsp_files <- list.files(path = wsp_dir, pattern = "\\.wsp$", full.names = TRUE)

# Initialize all required lists in one line
list_names <- c("fj_wsp_list", "gset_list", "sample_groups_list",
                "pop_path_list", "samp_name_list", "pop_strings_list",
                "total_counts", "percentages_of_parent", "remainder_counts",
                "cyto_list")
list2env(setNames(replicate(length(list_names), list(),
                            simplify = FALSE), list_names), envir = .GlobalEnv)

# Iterate through each workspace file
for (wsp_file in wsp_files) {
  # Use the workspace file name as a unique identifier
  wsp_id <- tools::file_path_sans_ext(basename(wsp_file))
  wsp_out_dir <- file.path(out_dir, wsp_id)
  dir.create(wsp_out_dir, showWarnings = FALSE, recursive = TRUE)

  # Open the FlowJo workspace
  fj_wsp_list[[wsp_id]] <- open_flowjo_xml(wsp_file)

  # Parse the workspace into a GatingSet
  gset_list[[wsp_id]] <- flowjo_to_gatingset(fj_wsp_list[[wsp_id]],
                                             path = fcs_dir)
  # TODO: If there is only one population sample group to choose from 
  # (be it All or whatever), then automatically choose that option, as I did a 
  # while back - it is possible.
  sample_groups_list[[wsp_id]] <- fj_ws_get_sample_groups(fj_wsp_list[[wsp_id]])

  pop_path_list[[wsp_id]] <- gs_get_pop_paths(gset_list[[wsp_id]],
                                              path = "full")
  samp_name_list[[wsp_id]] <- sampleNames(gset_list[[wsp_id]])

  # Population strings
  pop_strings_list[[wsp_id]] <- get_population_strings_list(gset_list[[wsp_id]],
                                                         samp_name_list[[wsp_id]])
  save_population_strings(pop_strings_list, wsp_out_dir)

  if(plot_barplot){
    # Total counts
    total_counts[[wsp_id]] <- get_population_stats(gating_set = gset_list[[wsp_id]])

    colnames(total_counts[[wsp_id]])[3] <- "Count"
    population_barplot(total_counts[[wsp_id]],
                       paste0(wsp_out_dir, "total_population_counts_barplot.png"),
                       "Total Population Counts by Sample", "Total Count")
    # Percentages of parent
    percentages_of_parent[[wsp_id]] <- get_population_stats(gset_list[[wsp_id]], 
                                                            type = "percent")
    colnames(percentages_of_parent[[wsp_id]])[3] <- "Count"
    population_barplot(percentages_of_parent[[wsp_id]],
                       paste0(wsp_out_dir, "percentages_of_parent_barplot.png"),
                       "Percentages of Parent Population by Sample",
                       "Percentage of Parent (%)")

    # Remainder counts
    remainder_counts[[wsp_id]] <- summarize_remainder_counts(pop_strings_list[[wsp_id]])
    population_barplot(remainder_counts[[wsp_id]], 
                       paste0(wsp_out_dir, "grouped_remainder_counts_barplot.png"),
                       "Remaining Population Counts by Sample", "Remainder Counts")
  }
}

# =============== Plot Venn diagram ================================


# ------------------ Map Workspace Sample Names to FCS Files -------------------
# Get the sample names from the GatingSet
wsp_sample_names <- sampleNames(gating_set)

# Get the actual FCS file names in the directory
fcs_files <- list.files(path = fcs_dir, pattern = "\\.fcs$", full.names = TRUE)
fcs_file_names <- basename(fcs_files)

# Map workspace sample names to actual FCS file names by removing suffixes
# from workspace sample names to match original FCS
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

# -------------------- Create DataFrames for Each Sample --------------------
# Initialize a list to store dataframes per sample
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

# -------------------- Save DataFrames to CSV --------------------
# Save each dataframe to a CSV file
for (wsp_name in names(dataframes_list)) {
  output_file <- paste0(out_dir, sub("\\.fcs", "", wsp_name),
                        "_", custom_name,"_populations.csv")
  print(output_file)
  write.csv(dataframes_list[[wsp_name]], file = output_file, row.names = FALSE)
  cat(paste0("Population dataframe saved to '", output_file, "'\n"))
}

# -------------------- Verify DataFrame Rows --------------------
nrow(dataframes_list[[1]])

