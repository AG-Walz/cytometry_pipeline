# Simple UpSet Plot for Metadata Analysis
# Optimized for handling 236 files with 7 features efficiently

# Set library path
.libPaths("Q:/IM/AG_Märklin/AG Märklin_Software/library/")

# Load required libraries
library(data.table)
library(UpSetR)

# ============================== Parameters ====================================

# Path to the metadata file
metadata_file <- "Q:/IM/AG_Märklin/AG Märklin_Software/FSFC_Pipeline/01_data_loading/all_metadata.csv"
output_dir <- "Q:/IM/AG_Märklin/AG Märklin_Software/FSFC_Pipeline/Results/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================== Simple UpSet Plot Function ====================================

create_complex_upset <- function(metadata_file, output_dir) {
  cat("Loading metadata for complex UpSet analysis...\n")
  
  # Read metadata
  metadata <- fread(metadata_file, sep = ";", header = TRUE)
  
  # Select only the metadata columns (exclude file identifiers)
  meta_cols <- c("Timepoint", "Pembro", "Menopause", "T-stadium", "N-stadium", "irAE", "Response postOP", "Response radio")
  
  # Include ALL rows, handle missing values explicitly
  clean_data <- metadata[, ..meta_cols]
  
  # Clean missing values but keep them as a category
  clean_data[] <- lapply(clean_data, function(x) {
    x <- as.character(x)
    x[is.na(x) | x == "" | x == " "] <- "Unknown"
    return(x)
  })
  
  cat("Total rows included:", nrow(clean_data), "\n")
  
  # Increase sample size for more complexity while keeping it manageable
  if (nrow(clean_data) > 200) {
    cat("Sampling 200 rows for detailed analysis...\n")
    set.seed(123)
    clean_data <- clean_data[sample(nrow(clean_data), 200)]
  }
  
  # Create comprehensive binary matrix with individual categories
  # Include all actual values, not just simplified groups
  upset_data <- data.frame(
    # All Timepoint values
    T1 = as.numeric(clean_data$Timepoint == "T1"),
    T2 = as.numeric(clean_data$Timepoint == "T2"),
    T3 = as.numeric(clean_data$Timepoint == "T3"),
    T4 = as.numeric(clean_data$Timepoint == "T4"),
    T5 = as.numeric(clean_data$Timepoint == "T5"),
    Timepoint_Unknown = as.numeric(clean_data$Timepoint == "Unknown"),
    
    # Pembro values
    Pembro_With = as.numeric(clean_data$Pembro == "With"),
    Pembro_Without = as.numeric(clean_data$Pembro == "Without"),
    Pembro_Unknown = as.numeric(clean_data$Pembro == "Unknown"),
    
    # Menopause values
    Menopause_Post = as.numeric(clean_data$Menopause == "Post"),
    Menopause_Pre = as.numeric(clean_data$Menopause == "Pre"),
    Menopause_Peri = as.numeric(clean_data$Menopause == "Peri"),
    Menopause_Unknown = as.numeric(clean_data$Menopause == "Unknown"),
    
    # T-stadium values
    T_Stage_1 = as.numeric(clean_data$`T-stadium` == "1"),
    T_Stage_2 = as.numeric(clean_data$`T-stadium` == "2"),
    T_Stage_3 = as.numeric(clean_data$`T-stadium` == "3"),
    T_Stage_4 = as.numeric(clean_data$`T-stadium` == "4"),
    T_Stage_Unknown = as.numeric(clean_data$`T-stadium` == "Unknown"),
    
    # N-stadium values
    N_Stage_0 = as.numeric(clean_data$`N-stadium` == "0"),
    N_Stage_1 = as.numeric(clean_data$`N-stadium` == "1"),
    N_Stage_2 = as.numeric(clean_data$`N-stadium` == "2"),
    N_Stage_Plus = as.numeric(clean_data$`N-stadium` == "+"),
    N_Stage_X = as.numeric(clean_data$`N-stadium` == "x"),
    N_Stage_Unknown = as.numeric(clean_data$`N-stadium` == "Unknown"),
    
    # irAE values
    irAE_Yes = as.numeric(clean_data$irAE == "Yes"),
    irAE_None = as.numeric(clean_data$irAE == "None"),
    irAE_Unknown = as.numeric(clean_data$irAE == "Unknown"),
    
    # Response postOP values
    PostOP_pCR = as.numeric(clean_data$`Response postOP` == "pCR"),
    PostOP_nonpCR = as.numeric(clean_data$`Response postOP` == "nonpCR"),
    PostOP_Unknown = as.numeric(clean_data$`Response postOP` == "Unknown"),
    
    # Response radio values
    Radio_CR = as.numeric(clean_data$`Response radio` == "CR"),
    Radio_PR = as.numeric(clean_data$`Response radio` == "PR"),
    Radio_Unknown = as.numeric(clean_data$`Response radio` == "Unknown")
  )
  
  # Remove columns that are all zeros (categories that don't exist in the data)
  non_empty_cols <- colSums(upset_data) > 0
  upset_data <- upset_data[, non_empty_cols, drop = FALSE]
  
  cat("Creating complex UpSet plot with", nrow(upset_data), "samples and", ncol(upset_data), "categories...\n")
  cat("Categories included:", paste(colnames(upset_data), collapse = ", "), "\n")
  
  # Create UpSet plot with enhanced settings for complexity
  png(paste0(output_dir, "complex_upset_plot.png"), width = 1600, height = 1000, res = 200)
  
  upset_plot <- upset(upset_data, 
                     sets = colnames(upset_data),
                     keep.order = FALSE,
                     order.by = "freq",
                     nsets = ncol(upset_data),  # Show all sets
                     nintersects = 40,          # Show top 40 intersections for more detail
                     mb.ratio = c(0.7, 0.3),
                     text.scale = c(1.4, 1.4, 1, 1, 1.6, 1.3),  # Adjusted for readability
                     point.size = 2.5,
                     line.size = 1.2,
                     mainbar.y.label = "Intersection Size",
                     sets.x.label = "Set Size",
                     # Color code by category groups
                     sets.bar.color = c(
                       rep("#1f77b4", 6),  # Timepoint (blues)
                       rep("#ff7f0e", 3),  # Pembro (oranges)
                       rep("#2ca02c", 4),  # Menopause (greens)
                       rep("#d62728", 5),  # T-stadium (reds)
                       rep("#9467bd", 6),  # N-stadium (purples)
                       rep("#8c564b", 3),  # irAE (browns)
                       rep("#e377c2", 3),  # Response postOP (pinks)
                       rep("#7f7f7f", 3)   # Response radio (grays)
                     )[1:ncol(upset_data)])  # Adjust to actual number of columns
  
  print(upset_plot)
  dev.off()
  
  # Also create a summary of the data
  summary_data <- data.frame(
    Category = colnames(upset_data),
    Count = colSums(upset_data),
    Percentage = round(colSums(upset_data) / nrow(upset_data) * 100, 1)
  )
  
  write.csv(summary_data, paste0(output_dir, "upset_category_summary.csv"), row.names = FALSE)
  
  cat("Complex UpSet plot saved successfully!\n")
  cat("Main plot:", paste0(output_dir, "complex_upset_plot.png"), "\n")
  cat("Summary data:", paste0(output_dir, "upset_category_summary.csv"), "\n")
  
  return(list(plot = upset_plot, summary = summary_data))
}

# ============================== Main Execution ====================================

cat("=== Complex UpSet Plot Generation ===\n")
cat("Input file:", metadata_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Check if file exists
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

# Create the complex plot
upset_result <- create_complex_upset(metadata_file, output_dir)

cat("\n=== Analysis Summary ===\n")
print(head(upset_result$summary, 10))

cat("\n=== Complete ===\n")


