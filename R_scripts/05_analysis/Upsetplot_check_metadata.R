# Check_metadata.R
# Script to analyze metadata CSV file and generate summary statistics and visualizations
# Reads all_metadata.csv and provides summaries for each column (excluding identifiers)

# Set library path
.libPaths("Q:/IM/AG_Märklin/AG Märklin_Software/library/")

# Load required libraries
library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Install and load UpSetR if not available
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  cat("Installing UpSetR package...\n")
  install.packages("UpSetR", lib = "Q:/IM/AG_Märklin/AG Märklin_Software/library/")
}
library(UpSetR)

# ============================== Parameters ====================================

# Path to the metadata file
metadata_file <- "Q:/IM/AG_Märklin/AG Märklin_Software/FSFC_Pipeline/01_data_loading/all_metadata.csv"
output_dir <- "Q:/IM/AG_Märklin/AG Märklin_Software/FSFC_Pipeline/Results/"

# Define columns to skip
skip_columns <- c("*file_name", "*processing", "Patients", "Samples", "*batch","*batch_passive")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================== Functions ====================================

# Function to clean and label empty values
clean_column <- function(x) {
  # Convert to character and handle various empty representations
  x <- as.character(x)
  x[is.na(x) | x == "" | x == " " | tolower(x) == "na"] <- "?"
  return(x)
}

# Function to create summary statistics for a column
summarize_column <- function(column_data, column_name) {
  # Clean the data
  clean_data <- clean_column(column_data)
  
  # Count frequencies
  freq_table <- table(clean_data, useNA = "no")
  freq_df <- data.frame(
    Value = names(freq_table),
    Count = as.numeric(freq_table),
    Percentage = round(as.numeric(freq_table) / length(clean_data) * 100, 2)
  )
  
  # Sort by count (descending)
  freq_df <- freq_df[order(freq_df$Count, decreasing = TRUE), ]
  
  return(list(
    column_name = column_name,
    total_entries = length(clean_data),
    unique_values = nrow(freq_df),
    frequency_table = freq_df
  ))
}

# Function to create barplot for a column
create_barplot <- function(column_data, column_name, output_dir) {
  # Clean the data
  clean_data <- clean_column(column_data)
  
  # Create frequency table
  freq_table <- table(clean_data)
  freq_df <- data.frame(
    Value = names(freq_table),
    Count = as.numeric(freq_table)
  )
  
  # Sort by count for better visualization
  freq_df <- freq_df[order(freq_df$Count, decreasing = TRUE), ]
  freq_df$Value <- factor(freq_df$Value, levels = freq_df$Value)
  
  # Create color palette for different bars
  n_colors <- nrow(freq_df)
  if (n_colors <= 8) {
    # Use distinct colors for small number of categories
    base_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#17becf")
  } else {
    # Use color palette for larger number of categories
    base_colors <- rainbow(n_colors, alpha = 0.7)
  }
  
  # Assign colors, making "Unassigned" always grey
  colors <- character(n_colors)
  color_index <- 1
  for (i in 1:n_colors) {
    if (freq_df$Value[i] == "?") {
      colors[i] <- "#808080"  # Grey color for Unassigned
    } else {
      colors[i] <- base_colors[color_index]
      color_index <- color_index + 1
      if (color_index > length(base_colors)) {
        color_index <- 1  # Wrap around if we run out of colors
      }
    }
  }
  
  # Create the plot
  p <- ggplot(freq_df, aes(x = Value, y = Count, fill = Value)) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.8) +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = Count), vjust = -0.3, size = 3, color = "black") +
    labs(
      title = paste("Distribution of", column_name),
      x = column_name,
      y = "Count",
      subtitle = paste("Total entries:", sum(freq_df$Count), "| Unique values:", nrow(freq_df))
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),  # Horizontal labels
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "none",  # Remove legend since colors are just for visual distinction
      panel.grid.major.x = element_blank(),  # Remove vertical grid lines
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.ticks.y = element_line(color = "gray50")  # Keep y-axis ticks
    )
  
  # Save the plot
  filename <- paste0(output_dir, "barplot_", gsub("[^A-Za-z0-9]", "_", column_name), ".png")
  ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)
  
  return(p)
}

# Function to create upset plot for combinatorial analysis 
if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
  install.packages("ComplexUpset")
}
library(ComplexUpset)
library(ggplot2)

create_upset_plot <- function(metadata, output_dir) {
  message("Preparing data for ComplexUpset plot...")

  # Select relevant columns
  required_cols <- c("Timepoint","Pembro","Menopause","T-stadium","N-stadium",
                     "irAE","Response postOP","Response radio")

  md <- metadata[, ..required_cols]
  
  # Memory optimization: Limit to maximum 1000 rows to prevent memory overflow
  if (nrow(md) > 1000) {
    message("Sampling 1000 rows to prevent memory issues...")
    set.seed(123)  # For reproducibility
    md <- md[sample(nrow(md), 1000), ]
  }

  # Keep all rows, treat "?" as valid level
  md[] <- lapply(md, function(x) {
    x <- as.character(x)
    x[x == "" | is.na(x)] <- "?"
    factor(x)
  })

  # Rename for safety in formula
  setnames(md,
           old = c("T-stadium","N-stadium","Response postOP","Response radio"),
           new = c("T_stadium","N_stadium","Response_postOP","Response_radio"))

  message("Rows used: ", nrow(md))

  # Memory-efficient one-hot encoding: Process in smaller chunks
  tryCatch({
    # Create smaller model matrix
    mm <- model.matrix(~ Timepoint + Pembro + Menopause + T_stadium +
                         N_stadium + irAE + Response_postOP + Response_radio - 1,
                       data = md)
  }, error = function(e) {
    # If still too large, reduce further
    message("Model matrix too large, reducing to 500 rows...")
    md_small <- md[sample(nrow(md), min(500, nrow(md))), ]
    mm <- model.matrix(~ Timepoint + Pembro + Menopause + T_stadium +
                         N_stadium + irAE + Response_postOP + Response_radio - 1,
                       data = md_small)
  })
  
  # Clean column names consistently
  colnames(mm) <- gsub("[^A-Za-z0-9]+", "_", colnames(mm))

  # Convert to upset format
  upset_df <- as.data.frame(mm)
  upset_df[] <- lapply(upset_df, function(x) as.integer(x > 0))

  # Order sets alphabetically so categories group visually
  upset_df <- upset_df[, sort(names(upset_df)), drop = FALSE]

  # Drop any empty sets to avoid downstream errors
  non_empty <- colSums(upset_df, na.rm = TRUE) > 0
  upset_df  <- upset_df[, non_empty, drop = FALSE]
  if (ncol(upset_df) == 0) stop("All sets are empty after filtering.")

  # Build intersections using the cleaned column names
  rows_to_sets <- apply(upset_df == 1L, 1, function(row) {
    # Get the cleaned column names where value is 1
    colnames(upset_df)[row]
  })
  
  # Create combination strings using cleaned names
  comb_strings <- vapply(rows_to_sets, function(s) paste(sort(s), collapse = "|"), character(1))
  comb_tab     <- sort(table(comb_strings), decreasing = TRUE)
  
  # Limit to top 20 intersections to reduce complexity and avoid errors
  comb_tab <- head(comb_tab[comb_tab >= 1], 20)

  # Convert back to a list of character vectors using cleaned names
  observed_intersections <- lapply(names(comb_tab), function(cs) {
    sets <- strsplit(cs, "\\|", fixed = TRUE)[[1]]
    # Verify all sets exist in our data
    valid_sets <- sets[sets %in% colnames(upset_df)]
    if (length(valid_sets) > 0) {
      return(valid_sets)
    } else {
      return(NULL)
    }
  })
  
  # Remove NULL intersections
  observed_intersections <- observed_intersections[!sapply(observed_intersections, is.null)]
  
  message("Found ", length(observed_intersections), " valid intersections")
  if (length(observed_intersections) == 0) {
    message("No valid intersections found, using automatic intersection detection")
    observed_intersections <- NULL  # Let ComplexUpset auto-detect
  }

  # Create plot with reduced complexity
  tryCatch({
    p <- ComplexUpset::upset(
      data               = upset_df,
      intersect          = colnames(upset_df),
      intersections      = observed_intersections,
      sort_sets          = "ascending",
      sort_intersections = "descending",
      width_ratio        = 0.30,

      # Simplified settings for better performance
      set_sizes = upset_set_size(position = "left"),

      base_annotations = list(
        "Intersection size" = intersection_size(
          text = list(size = 3)
        ) + labs(y = "Intersection size")
      ),

      themes = upset_modify_themes(list(
        "intersections" = theme(
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 10)
        ),
        "intersections_matrix" = theme(
          panel.grid  = element_blank(),
          axis.text.y = element_text(size = 8)
        )
      ))
    ) + labs(x = "Combination of features")
  }, error = function(e) {
    message("Error with custom intersections, using auto-detection: ", e$message)
    # Fallback: let ComplexUpset auto-detect intersections
    p <- ComplexUpset::upset(
      data               = upset_df,
      intersect          = colnames(upset_df),
      # No intersections parameter - let it auto-detect
      sort_sets          = "ascending",
      sort_intersections = "descending",
      width_ratio        = 0.30,
      min_size           = 1,  # Show intersections with at least 1 member
      n_intersections    = 20, # Limit to top 20

      set_sizes = upset_set_size(position = "left"),

      base_annotations = list(
        "Intersection size" = intersection_size(
          text = list(size = 3)
        ) + labs(y = "Intersection size")
      ),

      themes = upset_modify_themes(list(
        "intersections" = theme(
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 10)
        ),
        "intersections_matrix" = theme(
          panel.grid  = element_blank(),
          axis.text.y = element_text(size = 8)
        )
      ))
    ) + labs(x = "Combination of features")
  })

  # Save with reduced DPI for faster rendering
  ggsave(file.path(output_dir, "upset_plot_metadata_combinations.png"),
         p, width = 12, height = 7, dpi = 200)

  message("UpSet plot saved successfully! Total rows: ", nrow(upset_df))
  invisible(p)
}

# ============================== Main Analysis ====================================

cat("=== Metadata Analysis Report ===\n")
cat("Date:", Sys.Date(), "\n")
cat("File:", metadata_file, "\n\n")

# Read the metadata file
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

# Read CSV with semicolon separator
metadata <- fread(metadata_file, sep = ";", header = TRUE)

cat("Dataset Overview:\n")
cat("- Total rows:", nrow(metadata), "\n")
cat("- Total columns:", ncol(metadata), "\n")
cat("- Column names:", paste(colnames(metadata), collapse = ", "), "\n\n")

# Get columns to analyze
analyze_columns <- setdiff(colnames(metadata), skip_columns)

cat("Columns to analyze:", paste(analyze_columns, collapse = ", "), "\n")
cat("Skipped columns:", paste(skip_columns, collapse = ", "), "\n\n")

# Initialize lists to store results
summary_results <- list()
plots <- list()

# Analyze each column
for (col in analyze_columns) {
  cat("Analyzing column:", col, "\n")
  
  # Get summary statistics
  col_summary <- summarize_column(metadata[[col]], col)
  summary_results[[col]] <- col_summary
  
  # Create barplot
  plot <- create_barplot(metadata[[col]], col, output_dir)
  plots[[col]] <- plot
  
  # Print summary to console
  cat("  - Total entries:", col_summary$total_entries, "\n")
  cat("  - Unique values:", col_summary$unique_values, "\n")
  cat("  - Top 3 values:\n")
  
  top3 <- head(col_summary$frequency_table, 3)
  for (i in 1:nrow(top3)) {
    cat("    ", top3$Value[i], ":", top3$Count[i], "(", top3$Percentage[i], "%)\n")
  }
  cat("\n")
}


# ============================== Generate UpSet Plot ====================================

cat("Creating UpSet plot for combinatorial analysis...\n")

# Check if UpSetR is available
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  cat("Warning: UpSetR package not available. Skipping UpSet plot.\n")
  cat("To install UpSetR, run: install.packages('UpSetR')\n")
} else {
  tryCatch({
    upset_plot <- create_upset_plot(metadata, output_dir)
    cat("UpSet plot saved successfully!\n")
  }, error = function(e) {
    cat("Error creating UpSet plot:", e$message, "\n")
    cat("Continuing with other analyses...\n")
  })
}

# ============================== Save Detailed Results ====================================

# Save detailed frequency tables
for (col in names(summary_results)) {
  result <- summary_results[[col]]
  filename <- paste0(output_dir, "frequency_table_", gsub("[^A-Za-z0-9]", "_", col), ".csv")
  write.csv(result$frequency_table, file = filename, row.names = FALSE)
}

# Create a combined plot with all barplots
if (length(plots) > 0) {
  # Arrange plots in a grid
  n_plots <- length(plots)
  n_cols <- min(3, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  combined_plot <- do.call(grid.arrange, c(plots, ncol = n_cols))
  
  # Save combined plot
  ggsave(paste0(output_dir, "combined_metadata_plots.png"), 
         plot = combined_plot, 
         width = n_cols * 5, 
         height = n_rows * 4, 
         dpi = 300)
}

cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", output_dir, "\n")
cat("- Summary table: metadata_summary_table.csv\n")
cat("- Individual frequency tables: frequency_table_*.csv\n")
cat("- Individual barplots: barplot_*.png\n")
cat("- Combined plot: combined_metadata_plots.png\n")
cat("- UpSet plot: upset_plot_metadata_combinations.png\n")
