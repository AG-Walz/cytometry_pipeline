library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(pheatmap)

# Compact heatmap function
create_heatmap <- function(x, y, values, title = "Heatmap") {
  if (length(x) != length(y) || length(x) != length(values)) stop("All vectors must have same length")
  
  data.frame(x = x, y = y, values = values, stringsAsFactors = FALSE) %>%
    mutate(across(1:2, ~ str_trim(as.character(.x)))) %>%
    group_by(across(1:2)) %>%
    summarise(value = mean(values, na.rm = TRUE), .groups = 'drop') %>%
    setNames(c("x_var", "y_var", "value")) %>%
    pivot_wider(names_from = "x_var", values_from = "value", values_fill = NA) %>%
    column_to_rownames("y_var") %>%
    as.matrix()
}

# Load and process data
df <- read.csv("/home/ngr24/Desktop/Maerklin/Cytometry/gate_stats.csv") %>%
  mutate(
    gate = unlist(strsplit(gate, split=' (Manual gating)', fixed = TRUE)),
    marker = sub("\\|.*$", "", marker),
    marker = str_trim(marker),
    gate = str_trim(gate)
  ) %>%
  select(marker, gate, log2FC)

# Create matrix and determine colors
mat <- create_heatmap(df$marker, df$gate, df$log2FC)
data_range <- range(mat, na.rm = TRUE)

if (min(data_range) < 0 && max(data_range) > 0) {
  cols <- colorRampPalette(c("navy", "white", "firebrick"))(50)
} else if (max(data_range) <= 0) {
  cols <- colorRampPalette(c("navy", "lightblue"))(50)
} else {
  cols <- colorRampPalette(c("white", "firebrick"))(50)
}

# Create and save heatmap
p <- pheatmap(mat, color = cols, cluster_rows = TRUE, cluster_cols = TRUE,
              border_color = "white", cellwidth = 20, cellheight = 20,
              main = "Expression Heatmap: log2FC by Marker and Gate",
              show_rownames = TRUE, show_colnames = TRUE, na_col = "grey90",
              display_numbers = TRUE, number_format = "%.1f", fontsize_number = 8,
              silent = TRUE)

# Save outputs
output_dir <- dirname("/home/ngr24/Desktop/Maerklin/Cytometry/gate_stats.csv")
base_name <- paste0("heatmap_gate_stats_", Sys.Date())
png_file <- file.path(output_dir, paste0(base_name, ".png"))

png(png_file, width = 12, height = 8, units = "in", res = 300)
grid::grid.draw(p$gtable)
dev.off()

message("✓ Heatmap saved as: ", png_file)
message("Data range: ", round(min(data_range), 3), " to ", round(max(data_range), 3))
message("Matrix: ", nrow(mat), " rows × ", ncol(mat), " columns")

