# Required Libraries: flowWorkspace, dplyr, openxlsx, data.table
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# =============================== Parameters ===================================
input_fcs_path <- paste0("/data/Cytometry/ICK/ICK_T1_PR/files_ICK_T1_P+_poR_pre-gating_9_15_2025/")
input_csv_path <- input_fcs_path

col_pattern <- "Manual gating"

csv_label_dict <- c(
  "CD3n..Manual.gating." = "ncCD3-", "NKCs..Manual.gating." = "NK cells",
  "gdTCRnTCs..Manual.gating." = "ncT cells", "CD4TCs..Manual.gating." = "CD4T cells",
  "CD8TCs..Manual.gating." = "CD8T cells", "BCs..Manual.gating." = "B cells",
  "gdTCRnCD56nTCs..Manual.gating." = "ncT cells", "gdTCs..Manual.gating." = "gdT cells",
  "NKTCs..Manual.gating." = "NKT cells", "root" = "Removed"
)

gating_tree <- c(
  "Other" = 0, "T cells" = 1, "CD3-" = 1, "gdT cells" = 1,
  "NK cells" = 2, "NKT cells" = 2, "B cells" = 2,
  "CD8T cells" = 3, "CD4T cells" = 3
)

output_path <- paste0(input_fcs_path, "UMAP")

# ================================ Load Functions ==============================
# Source from parent directory
source("../Packages.R")
source("../Functions.R")
dir_create(output_path)

# ================================ Main Processing ============================
# Get file lists
fcs_files <- list.files(input_fcs_path, pattern = "\\.fcs$", full.names = TRUE)
csv_files <- list.files(input_csv_path, pattern = "\\.csv$", full.names = TRUE)

# Filter to valid CSV files with manual gating columns
valid_csvs <- filter_valid_csv_files(csv_files, col_pattern)

# Create mapping
csv_to_fcs <- create_csv_fcs_mapping(valid_csvs, fcs_files)
cat("CSV to FCS mapping created for", length(csv_to_fcs), "files\n")

# Process all CSV files
results <- vector("list", length(valid_csvs))
names(results) <- basename(valid_csvs)

for (i in seq_along(valid_csvs)) {
csv_file <- valid_csvs[i]
cat(sprintf("\n[%d/%d] ", i, length(valid_csvs)))

tryCatch({
    result_file <- process_single_csv(csv_files = valid_csvs,
                                      pattern = col_pattern,
                                      csv_to_fcs = csv_to_fcs,
                                      csv_label_dict = csv_label_dict,
                                      gating_tree = gating_tree,
                                      output_path = output_path,
                                      chunk_size = 100000 )
    results[[basename(csv_file)]] <- result_file
}, error = function(e) {
    cat("ERROR:", e$message, "\n")
    results[[basename(csv_file)]] <- NULL
})

gc()  # Force cleanup between files
}

cat("\nProcessing complete. Results saved to:", output_path, "\n")

# Summary of results
successful_results <- results[!sapply(results, is.null)]
validate_csv_results(successful_results)

cat("Processing complete. Results:", length(results), "files processed\n")
