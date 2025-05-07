# Required Libraries
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Example Input: Two lists of the same size
# Generate two lists of 35,000 entries each
set.seed(123)  # Set seed for reproducibility

# Define the possible strings
categories <- c("CD4TCs", "CD8TCs", "gdTCs", "BCs", "NCs", "NKTCs")

# Generate random entries for each list
list1 <- sample(categories, size = 35000, replace = TRUE)
list2 <- sample(categories, size = 35000, replace = TRUE)

# Print the first few entries of each list to verify
print(head(list1))
print(head(list2))

# Verify the lengths of the lists
print(length(list1))  # Should be 35000
print(length(list2))  # Should be 35000

# Ensure the lists are of the same size
if (length(list1) != length(list2)) {
  stop("The two lists must have the same size.")
}

# Get all unique strings across both lists
unique_strings <- unique(c(list1, list2))

# Create a directory to save the Venn diagrams
output_dir <- "venn_diagrams"
dir.create(output_dir, showWarnings = FALSE)

# Generate a Venn diagram for each unique string
for (string in unique_strings) {
  # Identify rows where the string appears in each list
  in_list1 <- which(list1 == string)
  in_list2 <- which(list2 == string)

  # Create the Venn diagram
  venn <- draw.pairwise.venn(
    area1 = length(in_list1),
    area2 = length(in_list2),
    cross.area = length(intersect(in_list1, in_list2)),
    category = c("Manual", "Automatic"),
    fill = c("blue", "red"),
    alpha = 0.5,
    cat.pos = c(-20, 20),
    cat.dist = c(0.05, 0.05),
    cex = 2,  # Increase text size in the diagram
    cat.cex = 2  # Increase category label size
  )

  # Save the Venn diagram to a separate file
  output_file <- file.path(output_dir, paste0("venn_", string, ".png"))
  png(filename = output_file, width = 600, height = 400)  # Adjusted size
  grid.draw(venn)
  grid.text(
    label = paste(string, "population comparison"),
    x = 0.5,
    y = 0.95,
    gp = gpar(fontsize = 24),  # Match font size to VennDiagram
    just = "center"
  )
  dev.off()
}

# Print completion message
print(paste("Venn diagrams saved to:", output_dir))

