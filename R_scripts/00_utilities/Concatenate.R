# ========================= PARAMETERS =========================

fallback_dir <- "/home/ngr24/Desktop/Maerklin/Cytometry/R_scripts"
input_file_path <- "/data/Maerklin/ICK/ICK_T1_PR/"
results_path<- "/data/Maerklin/ICK/ICK_T1_PR/Results/"
project_name <- "ICK_T1_PR"

# Marker selection options:
# - "popup": Open GUI window for interactive marker selection
# - "terminal": Choose from all filtered markers in terminal
# - "all": Use all filtered markers
# - vector of marker names: Use specific markers (e.g. c("CD4", "CD8", "TNF"))
# - vector of column numbers: Use specific columns (e.g. c(5, 7, 9))
MARKER_SELECTION_MODE <- "popup"  # Change this to control marker selection

# Patterns to use to filter out markers
EXCLUDE_PATTERNS <- c("^$", "ALD", "Column", "index", "Time", "FSC", "SSC")

# Output format options: Choose which formats to generate
# Set to TRUE for formats you want, FALSE to skip
# You can choose any combination: none, one, two, or all three formats
OUTPUT_FORMATS <- list(
  pdf = TRUE,    # Vector format for publications (scalable, small file size)
  png = TRUE,    # Raster format for pwp (high quality, larger file size)
  svg = TRUE     # Vector format for web (scalable, editable)
)

# Plot dimensions (inches)
PLOT_WIDTH <- 8
PLOT_HEIGHT <- 6
PNG_DPI <- 300  # Resolution for PNG files

# ========================= SETUP =========================
# Set working directory to script location
setwd(fallback_dir)

source_files <- c("Packages.R", "Functions.R", "UMAP_functions.R",
                  "Downsampling_functions.R") # UMAP.R

lapply(source_files, source, chdir = TRUE)

# Setup output directory
umap_dir <- paste0(results_path, "/umap_plots/")
dir_create(c(results_path, umap_dir))

# ========================= DIRECT FCS LOADING =========================

# Load FCS files into flow frames and extract sample names
flow_frames <- load_fcs_folder(input_file_path)
sample_names <- sapply(list.files(input_file_path,
                                  pattern = "\\.fcs$",
                                  full.names = FALSE),
                       function(x) sub("\\.fcs$", "", x))

# ========================= MARKER SELECTION =========================
# Get marker information and filter markers
marker_info <- explore_data(flow_frames[[1]])

# Filter out empty/unnamed markers and uninteresting ones
all_markers <- marker_info$markers[marker_info$mcols]
all_indices <- marker_info$mcols

# Apply filtering using base R to avoid namespace conflicts
exclude_regex <- paste(EXCLUDE_PATTERNS, collapse = "|")
keep_mask <- !grepl(exclude_regex, all_markers, ignore.case = TRUE)
filtered_markers <- all_markers[keep_mask]
filtered_indices <- all_indices[keep_mask]

# Simplified marker selection
marker_indices_of_interest <- select_markers_by_mode(MARKER_SELECTION_MODE,
                                                     filtered_markers,
                                                     filtered_indices)

cat("Final markers to use for UMAP:", length(marker_indices_of_interest),
    "markers | Indices:", paste(marker_indices_of_interest, collapse = ", "),
    "\n")

# ========================= CONCATENATE ============================

concatenate_to_file(flow_frames,
                    backingfile = file.path(results_path, "concatenated_data"),
                    col_indices = marker_indices_of_interest,
                    add_origin_col = TRUE,
                    overwrite = TRUE,
                    set_names = TRUE)

fbm <- bigstatsr::big_attach(paste0(file.path(results_path, "concatenated_data"), ".rds"))
cat("Concatenated data dimensions:", nrow(fbm), "cells x", ncol(fbm), "markers\n")

# ========================= EXPORT AS CSV ============================
feat_out <- file.path(results_path, "concatenated_data.csv")
lab_out  <- file.path(results_path, "labels.csv")

obj <- readRDS(file.path(results_path, "concatenated_data.rds"))

n <- nrow(obj); p <- ncol(obj)
if (p < 2) stop("Need at least 2 columns (last is label)")

feat_cols <- seq_len(p - 1)
label_col <- p

bigstatsr::big_apply(
  obj,
  a.FUN = function(X, ind) {
    chunk <- X[ind, , drop = FALSE]            # matrix (rows = ind, all cols)
    features <- chunk[, feat_cols, drop = FALSE]
    labels   <- chunk[, label_col, drop = FALSE]

    # Append without headers/quotes
    utils::write.table(features, feat_out, sep = ",", row.names = FALSE,
                       col.names = FALSE, append = TRUE, quote = FALSE)
    utils::write.table(labels, lab_out, sep = ",", row.names = FALSE,
                       col.names = FALSE, append = TRUE, quote = FALSE)
    NULL
  },
  a.combine = 'c',
  ind = rows_along(obj),          # all rows
  block.size = 10000              # adjust if you want bigger/smaller chunks
)

cat("Wrote:\n", feat_out, "\n", lab_out, "\n")