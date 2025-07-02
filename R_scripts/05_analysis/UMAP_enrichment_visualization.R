# UMAP Enrichment Visualization
# Colors points by cluster enrichment: pCR-enriched (red), npCR-enriched (blue), non-significant (grey)

# =========================
# Configuration - MODIFY THESE
# =========================
CONFIG <- list(
  # 1) Path to the saved results CSV 
  results_csv = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/cluster_enrichment_results.csv",
  
  # 2) Path to the cell-level cluster assignments CSV (NEW - exported by Fine_clustering script)
  cell_data_csv = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/cell_cluster_assignments.csv",
  
  # 3) Output directory for plots
  out_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/",
  
  # 4) UMAP column names - will be selected interactively via GUI
  # Set to NULL to trigger interactive selection, or specify column names directly
  umap_x = NULL,  # e.g., "UMAP1" or NULL for GUI selection
  umap_y = NULL,  # e.g., "UMAP2" or NULL for GUI selection
  
  # 5) Condition levels (must match your analysis)
  condition_levels = c("pCR", "npCR"),
  
  # 6) Memory settings
  max_total_cells = 100000,    # Maximum total cells to plot
  
  # 7) Plot settings
  plot_width = 12,
  plot_height = 10,
  plot_dpi = 300,
  point_size = 0.3,
  point_alpha = 0.6
)

# =========================
# Load Required Libraries
# =========================
message("Loading libraries...")
suppressMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(scales)
  library(RColorBrewer)
  library(tcltk)        # For GUI column selection
})

# =========================
# Helper Functions
# =========================
make_outdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created output directory: ", path)
  }
}

safe_sample_data <- function(data, max_rows) {
  if (nrow(data) > max_rows) {
    message("Sampling ", max_rows, " rows from ", nrow(data), " total")
    return(data %>% slice_sample(n = max_rows))
  }
  return(data)
}

# Interactive column selection GUI
select_columns_gui <- function(column_names, title = "Select Columns", multiple = TRUE) {
  require(tcltk)
  
  if (length(column_names) == 0) {
    message("No columns available for selection")
    return(character(0))
  }
  
  # Create main window
  tt <- tktoplevel()
  tkwm.title(tt, title)
  tkwm.geometry(tt, "500x400")
  
  # Main frame
  main_frame <- tkframe(tt)
  tkpack(main_frame, fill = "both", expand = TRUE, padx = 10, pady = 10)
  
  # Instructions
  instr_text <- if (multiple) "Select multiple columns (Ctrl+click):" else "Select one column:"
  instr_label <- tklabel(main_frame, text = instr_text, 
                         font = tkfont.create(size = 10, weight = "bold"))
  tkpack(instr_label, pady = 5)
  
  # Listbox with scrollbar
  list_frame <- tkframe(main_frame)
  tkpack(list_frame, fill = "both", expand = TRUE)
  
  scrollbar <- tkscrollbar(list_frame, repeatinterval = 5,
                          command = function(...) tkyview(lb, ...))
  select_mode <- if (multiple) "extended" else "single"
  lb <- tklistbox(list_frame, yscrollcommand = function(...) tkset(scrollbar, ...),
                  selectmode = select_mode, height = 15)
  
  tkpack(scrollbar, side = "right", fill = "y")
  tkpack(lb, side = "left", fill = "both", expand = TRUE)
  
  # Add columns to listbox
  for (col in column_names) {
    tkinsert(lb, "end", col)
  }
  
  # Buttons
  button_frame <- tkframe(main_frame)
  tkpack(button_frame, pady = 10)
  
  selected_items <- character(0)
  
  ok_function <- function() {
    selections <- as.integer(tkcurselection(lb))
    if (length(selections) > 0) {
      selected_items <<- column_names[selections + 1]
    }
    tkdestroy(tt)
  }
  
  cancel_function <- function() {
    selected_items <<- character(0)
    tkdestroy(tt)
  }
  
  ok_button <- tkbutton(button_frame, text = "OK", command = ok_function)
  cancel_button <- tkbutton(button_frame, text = "Cancel", command = cancel_function)
  
  tkpack(ok_button, side = "left", padx = 5)
  tkpack(cancel_button, side = "left", padx = 5)
  
  # Focus and wait
  tkfocus(lb)
  tkwait.window(tt)
  
  return(selected_items)
}

# =========================
# Load and Process Results
# =========================
message("Loading enrichment results...")

if (!file.exists(CONFIG$results_csv)) {
  stop("Results file not found: ", CONFIG$results_csv)
}

results <- read_csv(CONFIG$results_csv, show_col_types = FALSE)
message("Loaded results for ", nrow(results), " clusters")

# Classify clusters by enrichment direction
cluster_colors <- results %>%
  mutate(
    enrichment_type = case_when(
      significant & log2_enrichment_A_vs_B > 0 ~ paste0(CONFIG$condition_levels[1], "_enriched"),
      significant & log2_enrichment_A_vs_B < 0 ~ paste0(CONFIG$condition_levels[2], "_enriched"), 
      TRUE ~ "non_significant"
    ),
    color_group = case_when(
      enrichment_type == paste0(CONFIG$condition_levels[1], "_enriched") ~ CONFIG$condition_levels[1],
      enrichment_type == paste0(CONFIG$condition_levels[2], "_enriched") ~ CONFIG$condition_levels[2],
      TRUE ~ "Non-significant"
    )
  ) %>%
  select(cluster, enrichment_type, color_group, log2_enrichment_A_vs_B, p_adj)

message("Cluster classification:")
message("- ", CONFIG$condition_levels[1], " enriched: ", sum(cluster_colors$color_group == CONFIG$condition_levels[1]))
message("- ", CONFIG$condition_levels[2], " enriched: ", sum(cluster_colors$color_group == CONFIG$condition_levels[2]))
message("- Non-significant: ", sum(cluster_colors$color_group == "Non-significant"))

# =========================
# Load Saved Cell Data
# =========================
message("Loading cell-level cluster assignments...")

# Check if cell data file exists
if (!file.exists(CONFIG$cell_data_csv)) {
  stop("Cell data file not found: ", CONFIG$cell_data_csv,
       "\nMake sure you've run the updated Fine_clustering_cytolution.R script that exports cell_cluster_assignments.csv")
}

# Load cell data (this contains cell_index, file_name, condition, UMAP coordinates, and cluster assignments)
combined_cells <- read_csv(CONFIG$cell_data_csv, show_col_types = FALSE)
message("Loaded cell data: ", nrow(combined_cells), " cells")

# =========================
# Interactive UMAP Column Selection
# =========================
if (is.null(CONFIG$umap_x) || is.null(CONFIG$umap_y)) {
  message("\nInteractive column selection required...")
  
  all_columns <- names(combined_cells)
  message("Available columns: ", paste(all_columns, collapse = ", "))
  
  # Select UMAP X coordinate
  if (is.null(CONFIG$umap_x)) {
    message("\nPlease select UMAP X coordinate...")
    CONFIG$umap_x <- select_columns_gui(all_columns, "Select UMAP X coordinate", multiple = FALSE)
    if (length(CONFIG$umap_x) == 0) {
      stop("No UMAP X column selected")
    }
    message("Selected UMAP X: ", CONFIG$umap_x)
  }
  
  # Select UMAP Y coordinate
  if (is.null(CONFIG$umap_y)) {
    message("\nPlease select UMAP Y coordinate...")
    # Remove already selected X coordinate from options
    remaining_columns <- all_columns[all_columns != CONFIG$umap_x]
    CONFIG$umap_y <- select_columns_gui(remaining_columns, "Select UMAP Y coordinate", multiple = FALSE)
    if (length(CONFIG$umap_y) == 0) {
      stop("No UMAP Y column selected")
    }
    message("Selected UMAP Y: ", CONFIG$umap_y)
  }
  
  message("\nUMAP coordinates selected:")
  message("- X axis: ", CONFIG$umap_x)
  message("- Y axis: ", CONFIG$umap_y)
}

# Check if UMAP columns exist in the data
if (!CONFIG$umap_x %in% names(combined_cells) || !CONFIG$umap_y %in% names(combined_cells)) {
  stop("Selected UMAP columns not found in cell data. Available columns: ", 
       paste(names(combined_cells), collapse = ", "))
}

# Filter to cells with UMAP coordinates (remove NAs)
cells_with_umap <- combined_cells %>%
  filter(!is.na(.data[[CONFIG$umap_x]]) & !is.na(.data[[CONFIG$umap_y]]))

message("Cells with UMAP coordinates: ", nrow(cells_with_umap), " / ", nrow(combined_cells))

if (nrow(cells_with_umap) == 0) {
  stop("No cells have valid UMAP coordinates!")
}

# Sample if too large for visualization
if (nrow(cells_with_umap) > CONFIG$max_total_cells) {
  message("Sampling ", CONFIG$max_total_cells, " cells from ", nrow(cells_with_umap), " for visualization")
  cells_with_umap <- cells_with_umap %>% slice_sample(n = CONFIG$max_total_cells)
}

# =========================
# Prepare Data for Plotting
# =========================
message("Preparing data for visualization...")

# Find cluster columns in the cell data
cluster_cols <- names(cells_with_umap)[str_detect(names(cells_with_umap), "^cluster_|^c_")]
if (length(cluster_cols) == 0) {
  stop("No cluster columns found in cell data! Expected columns starting with 'cluster_' or 'c_'")
}

message("Found ", length(cluster_cols), " cluster columns")
message("Using ", nrow(cells_with_umap), " cells for visualization")

# Create a mapping for each cell to its enrichment category
# Each cell gets colored by its most significant enriched cluster
cell_enrichment_mapping <- cells_with_umap %>%
  select(all_of(c("cell_index", CONFIG$umap_x, CONFIG$umap_y, "file_name", cluster_cols))) %>%
  # Convert to long format to find which clusters each cell belongs to
  pivot_longer(all_of(cluster_cols), names_to = "cluster", values_to = "in_cluster") %>%
  filter(in_cluster == 1) %>%  # Only keep clusters the cell belongs to
  # Join with enrichment information
  left_join(cluster_colors, by = "cluster") %>%
  # For cells in multiple clusters, prioritize by significance
  group_by(cell_index) %>%
  arrange(p_adj) %>%  # Most significant first
  slice_head(n = 1) %>%  # Take most significant cluster
  ungroup()

# For cells not in any significant cluster or not found, assign "Non-significant"
all_cell_coords <- cells_with_umap %>%
  select(cell_index, all_of(c(CONFIG$umap_x, CONFIG$umap_y, "file_name", "condition")))

final_plot_data <- all_cell_coords %>%
  left_join(cell_enrichment_mapping %>% select(cell_index, color_group), by = "cell_index") %>%
  mutate(
    color_group = if_else(is.na(color_group), "Non-significant", color_group),
    color_group = factor(color_group, levels = c(CONFIG$condition_levels[1], CONFIG$condition_levels[2], "Non-significant"))
  )

# Summary of cell assignments
message("\nCell color assignments:")
color_summary <- final_plot_data %>% count(color_group, sort = TRUE)
for (i in 1:nrow(color_summary)) {
  message("- ", color_summary$color_group[i], ": ", format(color_summary$n[i], big.mark = ","), " cells")
}

# =========================
# Create UMAP Plot
# =========================
message("\nGenerating UMAP enrichment plot...")

make_outdir(CONFIG$out_dir)

# Define colors
plot_colors <- c(
  setNames("red", CONFIG$condition_levels[1]),
  setNames("blue", CONFIG$condition_levels[2]), 
  "Non-significant" = "lightgray"
)

# Create the plot
umap_enrichment_plot <- ggplot(final_plot_data, aes_string(x = CONFIG$umap_x, y = CONFIG$umap_y)) +
  geom_point(aes(color = color_group), size = CONFIG$point_size, alpha = CONFIG$point_alpha) +
  scale_color_manual(
    values = plot_colors,
    name = "Cluster Enrichment",
    breaks = c(CONFIG$condition_levels[1], CONFIG$condition_levels[2], "Non-significant")
  ) +
  coord_equal() +
  labs(
    title = "UMAP Colored by Cluster Enrichment",
    subtitle = paste0("Red = ", CONFIG$condition_levels[1], " enriched clusters, Blue = ", CONFIG$condition_levels[2], " enriched clusters, Gray = Non-significant"),
    x = CONFIG$umap_x,
    y = CONFIG$umap_y,
    caption = paste0("Total cells: ", format(nrow(final_plot_data), big.mark = ","))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

# Save plot
plot_filename <- file.path(CONFIG$out_dir, "UMAP_cluster_enrichment.png")
ggsave(plot_filename, umap_enrichment_plot, 
       width = CONFIG$plot_width, height = CONFIG$plot_height, dpi = CONFIG$plot_dpi)

message("✓ UMAP enrichment plot saved: ", plot_filename)

# =========================
# Create Faceted Plot by Condition
# =========================
message("Creating faceted plot by condition...")

faceted_plot <- ggplot(final_plot_data, aes_string(x = CONFIG$umap_x, y = CONFIG$umap_y)) +
  geom_point(aes(color = color_group), size = CONFIG$point_size, alpha = CONFIG$point_alpha) +
  scale_color_manual(
    values = plot_colors,
    name = "Cluster Enrichment"
  ) +
  facet_wrap(~condition, scales = "free") +
  coord_equal() +
  labs(
    title = "UMAP Colored by Cluster Enrichment (by Condition)",
    x = CONFIG$umap_x,
    y = CONFIG$umap_y
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

faceted_filename <- file.path(CONFIG$out_dir, "UMAP_cluster_enrichment_by_condition.png")
ggsave(faceted_filename, faceted_plot, 
       width = CONFIG$plot_width, height = CONFIG$plot_height, dpi = CONFIG$plot_dpi)

message("✓ Faceted UMAP plot saved: ", faceted_filename)

# =========================
# Summary Statistics
# =========================
message("\n=== UMAP ENRICHMENT VISUALIZATION COMPLETE ===")
message("Output directory: ", CONFIG$out_dir)
message("Generated files:")
message("- UMAP_cluster_enrichment.png: Main enrichment visualization")
message("- UMAP_cluster_enrichment_by_condition.png: Faceted by condition")

message("\nVisualization summary:")
message("- Total cells plotted: ", format(nrow(final_plot_data), big.mark = ","))
message("- Files represented: ", length(unique(final_plot_data$file_name)))
message("- Significant clusters: ", sum(results$significant))
message("- Color scheme: ", CONFIG$condition_levels[1], " = red, ", CONFIG$condition_levels[2], " = blue, Non-significant = gray")
message("- Data source: cell_cluster_assignments.csv (from Fine_clustering analysis)")

# Print most significant clusters
sig_clusters <- results %>% 
  filter(significant) %>% 
  arrange(p_adj) %>%
  slice_head(n = 10)

if (nrow(sig_clusters) > 0) {
  message("\nTop significant clusters visualized:")
  for (i in 1:nrow(sig_clusters)) {
    direction <- if (sig_clusters$log2_enrichment_A_vs_B[i] > 0) CONFIG$condition_levels[1] else CONFIG$condition_levels[2]
    message(sprintf("  %d. %s (enriched in %s, p_adj = %.2e)", 
                   i, sig_clusters$cluster[i], direction, sig_clusters$p_adj[i]))
  }
}

message("\n✓ UMAP enrichment visualization completed successfully!")