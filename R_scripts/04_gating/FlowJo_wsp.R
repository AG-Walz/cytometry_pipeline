# Set the R_LIBS_USER environment variable
.libPaths("Q:/IM/AG_Märklin/Software/library/")

# Load required libraries
library(flowWorkspace)
library(openCyto)
library(CytoML)

# Define the path to your FlowJo workspaces
workspace_paths <- list.files("F:/ICK/Workspaces", pattern = "\\.wsp$", full.names = TRUE)

# Initialize a list to store gate correspondences
gate_correspondences <- list()

# Loop through each FlowJo workspace
for (workspace_path in workspace_paths) {
  # Load the FlowJo workspace
  ws <- CytoML::open_flowjo_xml(workspace_path)
  
  # Parse the workspace into a GatingSet
  gs <- parseWorkspace(ws, name = "All Samples") # Adjust the group name if needed
  
  # Extract sample names
  sample_names <- sampleNames(gs)
  
  # Loop through each sample in the GatingSet
  for (sample in sample_names) {
    # Get the gating hierarchy for the sample
    gating_hierarchy <- gs[[sample]]
    
    # Extract gate information for each population
    gates <- getPopStats(gating_hierarchy)
    
    # Add the gate correspondences to the list
    gate_correspondences[[paste0(workspace_path, "_", sample)]] <- gates
  }
}

# Combine all gate correspondences into a single data frame
gate_correspondences_df <- do.call(rbind, lapply(names(gate_correspondences), function(name) {
  data <- gate_correspondences[[name]]
  data$Sample <- name
  return(data)
}))

# Save the gate correspondences to a CSV file
write.csv(gate_correspondences_df, "gate_correspondences.csv", row.names = FALSE)

# Print a message
cat("Gate correspondences exported to 'gate_correspondences.csv'.\n")