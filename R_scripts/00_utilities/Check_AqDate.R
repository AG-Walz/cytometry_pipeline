# dieses R Skript soll in einem Ordner, alle fcs files auch in unterordnern laden 
# und das datum der quisition ausgeben, sowie identische fcs files mit unterschiedlichem Namen eine gleiche nummer zuornen

# Rename the library path if needed
.libPaths("Q:/IM/AG_Märklin/AG Märklin_Software/library/")

library(flowCore)
library(data.table)

# =============================== Functions ====================================

# Function to get acquisition date from single FCS file
get_fcs_keyword <- function(fcs_path, keyword_name) {
  tryCatch({
    ff <- read.FCS(fcs_path, transformation = FALSE, which.lines = 1, truncate_max_range = FALSE)
    keyword_value <- keyword(ff, keyword_name)
    if (is.null(keyword_value) || is.na(keyword_value) || keyword_value == "") {
      return("Unknown")
    }
    return(keyword_value)
  }, error = function(e) {
    cat("Error reading file:", fcs_path, "- Error:", e$message, "\n")
    return("Error")
  })
}

# Function to explore available keywords in an FCS file (for debugging)
explore_keywords <- function(fcs_path) {
  tryCatch({
    ff <- read.FCS(fcs_path, transformation = FALSE, which.lines = 1, truncate_max_range = FALSE)
    keywords <- keyword(ff)
    
    cat("File:", basename(fcs_path), "\n")
    
    # Look for data type info
    if ("$DATATYPE" %in% names(keywords)) {
      cat("  $DATATYPE:", keywords[["$DATATYPE"]], "\n")
    }
    
    # Look for bit-related keywords
    bit_keywords <- names(keywords)[grepl("P[0-9]+B", names(keywords))]
    if (length(bit_keywords) > 0) {
      cat("  Parameter bit depths:\n")
      for (kw in bit_keywords) {
        cat("    ", kw, ":", keywords[[kw]], "\n")
      }
    }
    
    # Look for range keywords
    range_keywords <- names(keywords)[grepl("P[0-9]+R", names(keywords))]
    if (length(range_keywords) > 0) {
      cat("  Parameter ranges:\n")
      for (kw in range_keywords[1:min(3, length(range_keywords))]) {  # Show first 3
        range_val <- as.numeric(keywords[[kw]])
        if (!is.na(range_val)) {
          calculated_bits <- ceiling(log2(range_val + 1))
          cat("    ", kw, ":", keywords[[kw]], "(~", calculated_bits, "bits)\n")
        } else {
          cat("    ", kw, ":", keywords[[kw]], "\n")
        }
      }
    }
    
    return(keywords)
  }, error = function(e) {
    cat("Error exploring file:", fcs_path, "- Error:", e$message, "\n")
    return(NULL)
  })
}

# Function to get acquisition date and file size from single FCS file
get_fcs_metadata <- function(fcs_path) {
  tryCatch({
    ff <- read.FCS(fcs_path, transformation = FALSE, which.lines = 1, truncate_max_range = FALSE)
    
    # Get acquisition date
    acq_date <- keyword(ff, "$DATE")
    if (is.null(acq_date) || is.na(acq_date) || acq_date == "") {
      acq_date <- "Unknown"
    }
    
    # Get file size in bytes
    file_size <- file.size(fcs_path)
    if (is.na(file_size)) {
      file_size <- "Unknown"
    } else {
      file_size <- as.character(file_size)
    }
    
    return(c(acq_date, file_size))
  }, error = function(e) {
    cat("Error reading file:", fcs_path, "- Error:", e$message, "\n")
    return(c("Error", "Error"))
  })
}

# ============================== Parameters ====================================

# Parent folder that contains all subfolders with FCS files
parent_folder <- "Q:/IM/AG_Märklin/AG Märklin_ICK Studie Gyn/01 Panel Design/Referenzsammlung/"
output_filename_1 <- "Referenzen_aq_dates.csv"
output_filename_2 <- "Referenzen_aq_dates_filesize.csv"

#===============================  Loading ====================================== 

# Check if parent folder exists
if (!dir.exists(parent_folder)) {
  stop("Parent folder does not exist: ", parent_folder)
}

# Get all FCS files recursively (including from subfolders)
all_fcs_files <- list.files(parent_folder, pattern = "\\.fcs$", full.names = TRUE, recursive = TRUE)

# Check if any FCS files were found
if (length(all_fcs_files) == 0) {
  stop("No FCS files found in: ", parent_folder)
}

cat("Found", length(all_fcs_files), "FCS files\n")

# Explore keywords in first few files to understand what's available
cat("\n=== Exploring keywords in first 3 files ===\n")
for (i in 1:min(3, length(all_fcs_files))) {
  explore_keywords(all_fcs_files[i])
  cat("\n")
}

# Get relative file paths relative to parent folder
relative_paths <- normalizePath(all_fcs_files)
parent_folder_norm <- normalizePath(parent_folder)
# Make sure parent folder ends with a slash, then remove it from the full paths
if (!endsWith(parent_folder_norm, "/") && !endsWith(parent_folder_norm, "\\")) {
  parent_folder_norm <- paste0(parent_folder_norm, "/")
}
relative_paths <- substring(relative_paths, nchar(parent_folder_norm) + 1)

#================================ Main code ====================================

# Extract metadata for all files
cat("Extracting metadata (date and file size)...\n")
metadata_list <- t(sapply(all_fcs_files, function(file) {
  cat("Processing metadata:", basename(file), "\n")
  get_fcs_metadata(file)
}))

# Debug: Check dimensions
cat("Metadata list dimensions:", dim(metadata_list), "\n")
cat("Number of files:", length(all_fcs_files), "\n")
cat("Number of relative paths:", length(relative_paths), "\n")

# Set column names
colnames(metadata_list) <- c("acquisition_date", "file_size")

# Make final dataframe
acquisition_info <- data.frame(
  file_relative_path = relative_paths,
  acquisition_date = as.character(metadata_list[, "acquisition_date"]),
  file_size = as.character(metadata_list[, "file_size"]),
  stringsAsFactors = FALSE
)

# Assign unique IDs for files with the same acquisition date & file size
# Handle NA values properly before creating factors
paste_values <- paste(acquisition_info$acquisition_date, acquisition_info$file_size, sep = "_")
# Remove any NA values that might cause issues
paste_values[is.na(paste_values)] <- "Unknown_Unknown"

acquisition_info$unique_file_id <- as.integer(as.factor(paste_values))

# Sort by unique_file_id and optionally by acquisition_date
acquisition_info <- acquisition_info[order(acquisition_info$unique_file_id, acquisition_info$acquisition_date), ]

# Write to CSV
fwrite(acquisition_info, file = file.path(parent_folder, output_filename_2))
print(acquisition_info)

