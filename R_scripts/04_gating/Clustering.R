# ========================= PARAMETERS =========================

fallback_dir <- "/home/ngr24/Desktop/Maerklin/Cytometry/R_scripts"
input_file_path <- "/data/Maerklin/ICK/ICK_T1_PR/"
results_path <- "./Results"
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