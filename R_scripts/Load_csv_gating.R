# Required Libraries: flowWorkspace, dplyr, openxlsx
# !For RStudio only: Set working directory to the same path as the current file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# =============================== Parameters ===================================

# Taget FCS filename
target_fcs <- "ICK-049-02_Unmixed.fcs"

input_fcs_path <- paste0("/home/ngr24/Desktop/Maerklin/cytometry_pipeline/",
                         "Gating_comparisons/raw_data")

input_csv_path <- paste0("/home/ngr24/Desktop/Maerklin/cytometry_pipeline/",
                         "Gating_comparisons/Cytolution_broad_vs_narrow_gating/")

csv_label_dict <- c(
  "CD3n..Manual.gating." = "ncCD3-",
  "NKCs..Manual.gating." = "NK cells",
  "gdTCRnTCs..Manual.gating." = "ncT cells",
  "CD4TCs..Manual.gating." = "CD4T cells",
  "CD8TCs..Manual.gating." = "CD8T cells",
  "BCs..Manual.gating." = "B cells",
  "gdTCRnCD56nTCs..Manual.gating." = "ncT cells",
  "gdTCs..Manual.gating." = "gdT cells",
  "NKTCs..Manual.gating." = "NKT cells",
  "root" = "Removed"
)

gating_tree <- c(
  "Other" = 0,
  "T cells" = 1, "CD3-" = 1, "gdT cells" = 1,
  "NK cells" = 2, "NKT cells" = 2, "B cells" = 2,
  "CD8T cells" = 3, "CD4T cells" = 3
)

output_path <- paste0(input_fcs_path, "Gated_fo/")

# ================================ Functions ===================================

source("Packages.R")
source("Functions.R")
dir_create(output_path)

#================================= Load files ==================================
# Get list of FCS and csv files
fcs_files <- list.files(input_fcs_path, pattern = "\\.fcs$", full.names = TRUE)
csv_files <- list.files(input_csv_path, pattern = "\\.csv$", full.names = TRUE)

# Ensure files are in pairs by name, ignoring extensions
fcs_files_base <- gsub("\\.fcs$", "", basename(fcs_files))
csv_files_base <- gsub("\\.csv$", "", basename(csv_files))

#========================= Create CSV to FCS mapping ==========================

# TODO: Instead of extracting the FCS from the base names, extract from within

# Create a mapping from CSV to FCS
csv_to_fcs <- character(length(csv_files))
names(csv_to_fcs) <- csv_files

for (i in seq_along(csv_files)) {
  csv_file <- csv_files[i]
  # Normalize CSV and FCS base names (remove hyphens and underscores)
  csv_base <- gsub("\\.csv$", "", basename(csv_file))
  csv_base_norm <- gsub("[-_]", "", csv_base)
  fcs_bases_norm <- gsub("[-_]", "", gsub("\\.fcs$", "", basename(fcs_files)))
  # Find FCS base that is a substring of the CSV base (ignoring prefixes)
  match_idx <- which(sapply(fcs_bases_norm, function(fcs_norm) {
    grepl(fcs_norm, csv_base_norm, fixed = TRUE)}))
  if (length(match_idx) != 1) {
    stop(sprintf("Could not uniquely match FCS file for CSV '%s' (found %d matches)",
                 csv_base, length(match_idx)))
  }
  csv_to_fcs[csv_file] <- fcs_files[match_idx]
}

# Find the FCS of each csv using csv_to_fcs[csv_file]
# Use this mapping on downstream code

# Group CSVs by FCS
fcs_to_csvs <- split(names(csv_to_fcs), csv_to_fcs)

# Print mapping
cat("CSV to FCS mapping:\n")
for (fcs in names(fcs_to_csvs)) {
  cat(basename(fcs), ":\n")
  for (csv in fcs_to_csvs[[fcs]]) {
    cat("  -", basename(csv), "\n")
  }
}
# When processing, you can use csv_to_fcs to know which FCS each CSV belongs to.
# ================================ Main body =================================
# Read the first CSV file to get the gating columns
gating_data_example <- read.csv(list.files(input_csv_path, pattern = "\\.csv$",
                                           full.names = TRUE)[1], sep = ",",
                                header = TRUE)
check_gates <- colnames(gating_data_example)

# Check that all CSV files have the same column names
csv_data <- list()
label_lists <- list()
flat_list <- list()
idx_list <- list()
idx_labels <- list()
for (csv_file in csv_files) {
  separator <- detect_csv_sep(csv_file)
  csv_data[[csv_file]] <- read.csv(csv_file, sep = separator, header = TRUE)
  if (!identical(colnames(csv_data[[csv_file]]), check_gates)) {
    diff1 <- setdiff(check_gates, colnames(csv_data[[csv_file]] ))
    diff2 <- setdiff(colnames(csv_data[[csv_file]] ), check_gates)
    cat("Difference in columns for", basename(csv_file), ":\n")
    if (length(diff1) > 0) cat("  In selected_gates but not in file:",
                               paste(diff1, collapse = ", "), "\n")
    if (length(diff2) > 0) cat("  In file but not in selected_gates:",
                               paste(diff2, collapse = ", "), "\n")
    stop(paste("Column names in", basename(csv_file),
               "do not match the others."))
  }
  idx_list[[csv_file]] <-  csv_data[[csv_file]][, 1]
  colnames(csv_data[[csv_file]]) <- translate_labels(
    colnames(csv_data[[csv_file]]), csv_label_dict)
  label_lists[[csv_file]] <- labels_from_binary_df(csv_data[[csv_file]])
  flat_list[[csv_file]] <- flatten_label_list(label_lists[[csv_file]],
                                              gating_tree)
  # Extract the correct CSV-associated FCS file from the CSV mapping
  fcs_file <- csv_to_fcs[[csv_file]]
  num_events <- nrow(flowCore::read.FCS(fcs_file, truncate_max_range = FALSE))
  idx_labels_vec <- rep("", num_events)
  idx_labels_vec[idx_list[[csv_file]]+1] <- translate_labels(flat_list[[csv_file]],
                                                           csv_label_dict)
  # Initialize or add column
  idx_labels[[csv_file]] <- data.frame(idx_labels_vec, stringsAsFactors = FALSE)
}

