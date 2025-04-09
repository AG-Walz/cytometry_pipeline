if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("flowWorkspace", "openCyto", "ggcyto"))

library(flowWorkspace)
library(ggcyto)

# Load the FlowJo workspace
wsp_file <- "path/to/your/workspace.wsp"
ws <- open_flowjo_xml(wsp_file)

# Parse the workspace to extract gating information
gs <- parseWorkspace(ws, name = "Sample Group")  # Replace "Sample Group" with your group name

# Plot a specific gate
ggcyto(gs[[1]], aes(x = FSC-A, y = SSC-A)) + geom_gate("GateName")  # Replace "GateName" with your gate name

# Get the gate
gate <- getGate(gs[[1]], "GateName")

# Modify the gate (example: change the boundary)
gate@boundaries <- matrix(c(100, 200, 300, 400), ncol = 2)

# Set the modified gate back
setGate(gs[[1]], "GateName", gate)
recompute(gs[[1]])  # Recompute the gating hierarchy

save_gs(gs, path = "path/to/save/updated_workspace")

# Extract population statistics
stats <- getPopStats(gs)

# Plot percentage shifts
library(ggplot2)
ggplot(stats, aes(x = Population, y = Percent)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Population Percentage Shifts", x = "Population", y = "Percentage")