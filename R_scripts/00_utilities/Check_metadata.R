# compact_metadata_barplots_fixed.R
# One compact faceted plot; per-facet x scales so columns are never combined.

# Set library path
.libPaths("Q:/IM/AG_Märklin/AG Märklin_Software/library/")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

# ============== EDIT HERE ======================================================
metadata_file <- "Q:/IM/AG_Märklin/AG Märklin_Software/FSFC_Pipeline/01_data_loading/all_metadata_datum.csv"
output_dir    <- "Q:/IM/AG_Märklin/AG Märklin_Software/FSFC_Pipeline/Results/"
sep_char      <- ","

# Skip set (now includes "Datum" and "File index")
skip_columns <- c("Datum", "File index", "*file_name", "*processing","T-stadium", "N-stadium",
                  "Patients", "Samples", "*batch", "*batch_passive", "Menopause")

# Only analyze these (after skipping); NULL => all non-skipped
only_cols <- NULL

single_plot_name <- "consensus_barplots_compact.png"
summary_csv      <- "metadata_summary_table.csv"

# Facet density (columns per row)
facets_per_row <- NULL

# Sorting before splitting
# Set the name of the column to use for sorting rows prior to creating splits.
# Change this value to match your metadata column (e.g. "Datum"). If the column can be parsed 
# as Date/POSIX it will be sorted accordingly; otherwise it's sorted as character. 
sort_column <- "Datum"

# Set sort_desc = TRUE to sort in descending order. Otherwise ascending will be the default
sort_desc   <- FALSE
# ==============================================================================

# "None" is NOT Unassigned; keep only true empties / NA-like placeholders
placeholders <- c("?", "na", "n/a", "null", "-")  # lower-case compare

to_unassigned <- function(x) {
  x_chr  <- as.character(x)
  x_trim <- str_trim(x_chr)
  blank  <- x_trim == "" | is.na(x_trim)
  low    <- tolower(x_trim)
  place  <- low %in% placeholders
  out <- x_trim
  out[blank | place] <- "Unassigned"
  out
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
stopifnot(file.exists(metadata_file))

md <- suppressWarnings(fread(metadata_file, sep = sep_char, header = TRUE, data.table = FALSE))

# --- Optional sorting by a chosen column before splitting ---
if (!is.null(sort_column) && nzchar(sort_column)) {
  if (!(sort_column %in% colnames(md))) {
    warning(sprintf("Sort column '%s' not found in metadata; skipping sort.", sort_column))
  } else {
    col_vals <- md[[sort_column]]
    sorted <- FALSE

    # Try common date formats explicitly (DD/MM/YYYY, DD-MM-YYYY, etc.)
    date_formats <- c("%d/%m/%Y", "%d-%m-%Y", "%d.%m.%Y", "%Y-%m-%d", "%Y/%m/%d", 
                      "%m/%d/%Y", "%m-%d-%Y")
    
    for (fmt in date_formats) {
      parsed_date <- suppressWarnings(as.Date(col_vals, format = fmt))
      # Check if at least some dates were parsed successfully (not all NA)
      if (sum(!is.na(parsed_date)) > 0) {
        ord <- order(parsed_date, decreasing = sort_desc, na.last = TRUE)
        md  <- md[ord, , drop = FALSE]
        message(sprintf("Sorted metadata by column '%s' as Date (format: %s).", sort_column, fmt))
        sorted <- TRUE
        break
      }
    }
    
    # If date parsing failed, try POSIXct
    if (!sorted) {
      parsed_posix <- suppressWarnings(as.POSIXct(col_vals))
      if (!all(is.na(parsed_posix))) {
        ord <- order(parsed_posix, decreasing = sort_desc, na.last = TRUE)
        md  <- md[ord, , drop = FALSE]
        message(sprintf("Sorted metadata by column '%s' as POSIXct.", sort_column))
        sorted <- TRUE
      }
    }
    
    # Fallback to character sort
    if (!sorted) {
      ord <- order(as.character(col_vals), decreasing = sort_desc, na.last = TRUE)
      md  <- md[ord, , drop = FALSE]
      message(sprintf("Sorted metadata by column '%s' as character.", sort_column))
    }
  }
}
all_cols    <- colnames(md)
analyzable  <- setdiff(all_cols, skip_columns)

if (is.null(only_cols)) {
  analyze_cols <- analyzable
} else {
  missing <- setdiff(only_cols, analyzable)
  if (length(missing)) stop("Requested only_cols not analyzable (missing or in skip list): ",
                            paste(missing, collapse = ", "))
  analyze_cols <- only_cols
}
if (!length(analyze_cols)) stop("No columns selected to analyze.")

# Long format, per-column unassigned handling
long <- md |>
  dplyr::select(all_of(analyze_cols)) |>
  dplyr::mutate(across(everything(), ~ as.character(.))) |>
  tidyr::pivot_longer(
    everything(),
    names_to  = "column",
    values_to = "value_raw"
  ) |>
  dplyr::mutate(
    value         = to_unassigned(value_raw),
    is_unassigned = value == "Unassigned"
  )

suppressPackageStartupMessages({ library(scales) })  # for hue_pal()

# Order within each facet: by frequency desc, Unassigned last
counts <- long |>
  dplyr::count(column, value, name = "Count") |>
  dplyr::group_by(column) |>
  dplyr::mutate(Total = sum(Count)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    is_unassigned = value == "Unassigned",
    value_tmp     = dplyr::if_else(is_unassigned, "zzz__Unassigned__", value)
  ) |>
  dplyr::group_by(column) |>
  dplyr::mutate(value_tmp = forcats::fct_reorder(value_tmp, Count, .desc = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(value_tmp = forcats::fct_recode(value_tmp, Unassigned = "zzz__Unassigned__"))

# Build a fill key so "Unassigned" = grey, others each get a distinct color
counts <- counts |>
  dplyr::mutate(fill_key = dplyr::if_else(value_tmp == "Unassigned", "Unassigned", as.character(value_tmp)))

# Create a manual palette: grey for Unassigned + hue palette for the rest
fill_levels <- unique(counts$fill_key)
other_levels <- setdiff(fill_levels, "Unassigned")
pal <- c("Unassigned" = "grey70")
if (length(other_levels) > 0) {
  pal <- c(pal, setNames(hue_pal()(length(other_levels)), other_levels))
}

# Facet into exactly 2 rows
k_cols <- length(unique(counts$column))
nrow_facets <- 2
ncol_facets <- ceiling(k_cols / nrow_facets)

p <- ggplot(counts, aes(x = value_tmp, y = Count, fill = fill_key)) +
  geom_col() +
  # counts atop each bar
  geom_text(aes(label = Count), vjust = -0.25, size = 3) +
  facet_wrap(~ column, scales = "free_x", nrow = nrow_facets) +
  scale_fill_manual(values = pal, name = "Category") +
  labs(x = NULL, y = "Count") +
  # leave a little headroom for labels
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Save with width based on number of facet columns; height for exactly 2 rows
plot_width  <- max(8, 2.1 * ncol_facets)
plot_height <- max(6, 3.0 * nrow_facets)

ggsave(filename = file.path(output_dir, single_plot_name),
       plot = p, width = plot_width, height = plot_height, dpi = 300)

cat("Saved compact 2-row plot: ", file.path(output_dir, single_plot_name), "\n", sep = "")

# ====================== N PLOT: N splits, fixed batch colours (+ grey for Unassigned) ======================
# Set number of splits (equal-sized by row index)
n_splits <- 3   # <-- change to 2, 3, 4, 5, ...

fixed_batch_plot <- sprintf("consensus_barplots_compact_bySplit_%d.png", n_splits)

# --- Helpers for colours ---
# Fixed batch colours (extend/cycle as needed)
batch_base <- c(
  "Split 1" = "#D62728", # red
  "Split 2" = "#2CA02C", # green
  "Split 3" = "#1F77B4", # blue
  "Split 4" = "#FF7F0E", # orange
  "Split 5" = "#9467BD", # purple
  "Split 6" = "#FFFF00", # yellow
  "Split 7" = "#8C564B", # brown
  "Split 8" = "#17BECF", # cyan
  "Split 9" = "#E377C2", # pink
  "Split 10" = "#90EE91" # lightgreen
)
# ensure we have enough colours for any n_splits
batch_levels <- paste0("Split ", seq_len(n_splits))
batch_cols   <- setNames(rep(batch_base, length.out = n_splits), batch_levels)

# --- Build splits by index ---
md_split_b <- md |>
  dplyr::mutate(
    .row_id = dplyr::row_number(),
    split   = paste0("Split ", dplyr::ntile(.row_id, n_splits))
  )

# Long format with split (force character before pivot)
long_bfix <- md_split_b |>
  dplyr::select(all_of(analyze_cols), split) |>
  dplyr::mutate(across(all_of(analyze_cols), ~ as.character(.))) |>
  tidyr::pivot_longer(
    cols      = all_of(analyze_cols),
    names_to  = "column",
    values_to = "value_raw"
  ) |>
  dplyr::mutate(
    value         = to_unassigned(value_raw),
    is_unassigned = value == "Unassigned",
    split         = factor(split, levels = batch_levels)
  )

# Counts per (column, value, split)
counts_bfix <- long_bfix |>
  dplyr::count(column, value, split, name = "Count")

# Order categories per column by TOTAL freq across splits; keep Unassigned last
totals_bfix <- counts_bfix |>
  dplyr::group_by(column, value) |>
  dplyr::summarise(Total = sum(Count), .groups = "drop")

counts_bfix <- counts_bfix |>
  dplyr::left_join(totals_bfix, by = c("column", "value")) |>
  dplyr::mutate(
    is_unassigned = value == "Unassigned",
    value_tmp     = dplyr::if_else(is_unassigned, "zzz__Unassigned__", value)
  ) |>
  dplyr::group_by(column) |>
  dplyr::mutate(value_tmp = forcats::fct_reorder(value_tmp, Total, .desc = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(value_tmp = forcats::fct_recode(value_tmp, Unassigned = "zzz__Unassigned__"))

# Label totals for the bar tops
label_bfix <- totals_bfix |>
  dplyr::mutate(
    value_tmp = dplyr::if_else(value == "Unassigned", "Unassigned", value),
    value_tmp = dplyr::if_else(value_tmp == "Unassigned", "zzz__Unassigned__", value_tmp)
  ) |>
  dplyr::group_by(column) |>
  dplyr::mutate(value_tmp = forcats::fct_reorder(value_tmp, Total, .desc = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(value_tmp = forcats::fct_recode(value_tmp, Unassigned = "zzz__Unassigned__"))

# --- One-layer stacking with fixed batch colours for all values (including Unassigned) ---
# Use the same batch colors for both assigned and unassigned values
# Add transparency to unassigned values to make assigned values pop out more
counts_bfix <- counts_bfix |>
  dplyr::mutate(
    fill_key = as.character(split),
    alpha_key = dplyr::if_else(value_tmp == "Unassigned", "Unassigned", "Assigned")
  )

# Build the palette: use batch_cols for all splits
pal_bfix <- batch_cols

# Ensure stack order: Batch 1 on TOP, then Batch 2, ...
# (reverse=TRUE stacks in reverse of the factor order; with levels Split 1..N,
#  reverse=TRUE puts Split 1 on top.)
stack_pos <- position_stack(reverse = TRUE)

# Facet in exactly 2 rows
k_cols <- length(unique(counts_bfix$column))
nrow_facets <- 2
ncol_facets <- ceiling(k_cols / nrow_facets)

p_bfix <- ggplot(counts_bfix, aes(x = value_tmp, y = Count, fill = fill_key, alpha = alpha_key)) +
  geom_col(position = stack_pos) +
  # total counts atop stacks
  geom_text(
    data = label_bfix,
    aes(x = value_tmp, y = Total, label = Total),
    vjust = -0.25, size = 3, inherit.aes = FALSE
  ) +
  facet_wrap(~ column, scales = "free_x", nrow = nrow_facets) +
  scale_fill_manual(
    values = pal_bfix,
    name = "Batches",
    breaks = names(batch_cols),  # Show only the main batch colors in legend
    labels = names(batch_cols)   # Use the batch names as labels
  ) +
  scale_alpha_manual(
    values = c("Assigned" = 1.0, "Unassigned" = 0.7),
    guide = "none"  # Don't show alpha in legend
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "Count") +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

plot_width_fb  <- max(8, 2.1 * ncol_facets)
plot_height_fb <- max(6, 3.0 * nrow_facets)

ggsave(
  filename = file.path(output_dir, fixed_batch_plot),
  plot = p_bfix, width = plot_width_fb, height = plot_height_fb, dpi = 300
)

cat("Saved fixed-batch-colour plot: ", file.path(output_dir, fixed_batch_plot), "\n", sep = "")


# ====================== UPSET PLOT: Combinations across full dataset ======================

# Load UpSetR package
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  cat("Installing UpSetR package...\n")
  install.packages("UpSetR", lib = "Q:/IM/AG_Märklin/AG Märklin_Software/library/", repos = "https://cran.r-project.org")
}
suppressPackageStartupMessages(library(UpSetR))

cat("\n=== Creating UpSet plot for metadata combinations ===\n")

# SIMPLIFIED APPROACH: Create one set per column_value combination
# Select a subset of most important columns to avoid explosion
upset_cols_to_use <- analyze_cols[1:min(5, length(analyze_cols))]

cat("Using columns for UpSet:", paste(upset_cols_to_use, collapse = ", "), "\n")

# Prepare data - select these columns and convert to character
upset_df <- md |>
  dplyr::select(all_of(upset_cols_to_use)) |>
  dplyr::mutate(across(everything(), ~ as.character(.)))

# For each column, get ALL non-empty values (or limit to reasonable number)
upset_sets <- list()
set_names <- c()

for (col_name in upset_cols_to_use) {
  col_values <- upset_df[[col_name]]
  
  # Remove empty/NA values
  col_values_clean <- col_values[!is.na(col_values) & col_values != "" & 
                                   col_values != "?" & col_values != "Unassigned"]
  
  if (length(col_values_clean) > 0) {
    # Get all unique values (or top 10 if too many)
    freq_table <- sort(table(col_values_clean), decreasing = TRUE)
    
    # Use all values, but limit to top 10 per column to prevent explosion
    max_values_per_col <- 10
    top_values <- names(freq_table)[1:min(max_values_per_col, length(freq_table))]
    
    cat("  Column '", col_name, "': using ", length(top_values), " values: ", 
        paste(top_values, collapse = ", "), "\n", sep = "")
    
    for (val in top_values) {
      set_name <- paste0(col_name, ":", val)
      set_names <- c(set_names, set_name)
      upset_sets[[set_name]] <- which(upset_df[[col_name]] == val)
    }
  }
}

cat("\nTotal sets created:", length(set_names), "\n")
cat("Set names:", paste(head(set_names, 10), collapse = ", "), "...\n")

# Print set sizes for debugging
cat("\nSet sizes:\n")
for (i in 1:min(10, length(set_names))) {
  cat("  ", set_names[i], ":", length(upset_sets[[set_names[i]]]), "\n")
}

if (length(upset_sets) < 2) {
  cat("Warning: Not enough sets for UpSet plot (need at least 2). Skipping.\n")
} else {
  # Create UpSet plot using fromList
  upset_plot_file <- file.path(output_dir, "upset_plot_metadata_combinations.png")
  
  cat("\nGenerating UpSet plot...\n")
  
  # Convert to binary matrix format (this is what UpSetR actually needs)
  cat("Converting to binary matrix format...\n")
  upset_binary <- fromList(upset_sets)
  
  cat("Binary matrix dimensions:", nrow(upset_binary), "rows x", ncol(upset_binary), "columns\n")
  cat("Column names:", paste(head(colnames(upset_binary), 10), collapse = ", "), "\n")
  cat("First few row sums:", paste(head(rowSums(upset_binary)), collapse = ", "), "\n")
  
  # Check if matrix is valid
  if (nrow(upset_binary) == 0 || ncol(upset_binary) == 0) {
    cat("ERROR: Binary matrix is empty!\n")
  } else if (sum(upset_binary) == 0) {
    cat("ERROR: Binary matrix has no 1s!\n")
  } else {
    cat("Binary matrix looks valid. Total 1s:", sum(upset_binary), "\n")
    
    # Smart sorting: group by column name, then sort values alphabetically within each column
    # Extract column and value from set names (format: "Column:Value")
    col_names <- colnames(upset_binary)
    set_info <- data.frame(
      full_name = col_names,
      column = sapply(strsplit(col_names, ":"), function(x) x[1]),
      value = sapply(strsplit(col_names, ":"), function(x) ifelse(length(x) > 1, x[2], "")),
      stringsAsFactors = FALSE
    )
    
    # Sort by column first, then by value (in descending order for reverse display)
    set_info <- set_info[order(set_info$column, set_info$value, decreasing = TRUE), ]
    sorted_names <- set_info$full_name
    
    # Reorder binary matrix columns
    upset_binary <- upset_binary[, sorted_names, drop = FALSE]
    
    cat("Columns after smart sort (by column then value, reversed):\n")
    cat(paste("  ", head(colnames(upset_binary), 15), collapse = "\n"), "\n")
    
    png(filename = upset_plot_file, width = 1800, height = 1000, res = 120)
    
    tryCatch({
      # Create the plot directly with the binary matrix
      # Filter out intersections with size < 2
      # Sets will appear grouped by column, with values sorted alphabetically
      print(upset(
        upset_binary,
        nsets = ncol(upset_binary),
        nintersects = 50,  # Show more intersections now that we have more sets
        order.by = "freq",
        decreasing = TRUE,
        mb.ratio = c(0.6, 0.4),
        text.scale = c(1.5, 1.4, 1.2, 1.3, 1.6, 1.3),
        point.size = 4,
        line.size = 1.5,
        mainbar.y.label = "Intersection Size",
        sets.x.label = "Set Size",
        show.numbers = TRUE,
        cutoff = 2,  # Only show intersections with 2 or more samples
        sets = colnames(upset_binary),  # Use our custom sort order
        keep.order = TRUE  # Keep the order we specified
      ))
    }, error = function(e) {
      cat("ERROR creating UpSet plot:", e$message, "\n")
      cat("Traceback:", paste(capture.output(traceback()), collapse = "\n"), "\n")
      
      # Create diagnostic plot
      par(mar = c(2, 2, 2, 2))
      plot(1, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
           xlab = "", ylab = "", axes = FALSE)
      text(5, 7, "UpSet Plot Generation Failed", cex = 2, font = 2, col = "red")
      text(5, 5, paste("Error:", e$message), cex = 1.2)
      text(5, 3, paste("Sets:", ncol(upset_binary)), cex = 1.2)
      text(5, 2, paste("Samples:", nrow(upset_binary)), cex = 1.2)
    })
    
    dev.off()
    
    cat("Saved UpSet plot to:", upset_plot_file, "\n")
  }
}

cat("\n=== UpSet plot section complete ===\n\n")


# ====================== TEMPORAL HEATMAP: Feature values over time ======================

cat("\n=== Creating temporal heatmap of feature values ===\n")

# Use the date-sorted metadata (md is already sorted by Datum at this point)
# Extract Datum column for x-axis
if ("Datum" %in% colnames(md)) {
  datum_col <- md$Datum
  
  # Parse dates (try multiple formats as before)
  date_formats <- c("%d/%m/%Y", "%d-%m-%Y", "%d.%m.%Y", "%Y-%m-%d", "%Y/%m/%d", 
                    "%m/%d/%Y", "%m-%d-%Y")
  
  parsed_dates <- NULL
  for (fmt in date_formats) {
    parsed_dates <- suppressWarnings(as.Date(datum_col, format = fmt))
    if (sum(!is.na(parsed_dates)) > 0) {
      break
    }
  }
  
  # If date parsing failed, use row index
  if (is.null(parsed_dates) || all(is.na(parsed_dates))) {
    cat("Warning: Could not parse Datum column as dates. Using row indices.\n")
    timepoint_labels <- seq_len(nrow(md))
    timepoint_ids <- timepoint_labels
  } else {
    # Format dates as "MMM YYYY" for x-axis labels (compact)
    timepoint_labels <- format(parsed_dates, "%b %Y")
    timepoint_ids <- seq_len(nrow(md))
  }
  
  # Create long format with timepoint information
  heatmap_data <- md |>
    dplyr::mutate(
      timepoint_id = timepoint_ids,
      timepoint_label = timepoint_labels
    ) |>
    dplyr::select(timepoint_id, timepoint_label, all_of(analyze_cols)) |>
    dplyr::mutate(across(all_of(analyze_cols), ~ as.character(.))) |>
    tidyr::pivot_longer(
      cols = all_of(analyze_cols),
      names_to = "column",
      values_to = "value"
    ) |>
    dplyr::mutate(
      value = to_unassigned(value),
      feature = paste0(column, ":", value)
    )
  
  # Count unique features
  n_features <- length(unique(heatmap_data$feature))
  n_timepoints <- length(unique(heatmap_data$timepoint_id))
  
  cat("Total features (column:value pairs):", n_features, "\n")
  cat("Total timepoints:", n_timepoints, "\n")
  
  # Create binary presence matrix (1 if feature exists at that timepoint, 0 otherwise)
  heatmap_matrix <- heatmap_data |>
    dplyr::distinct(timepoint_id, timepoint_label, feature) |>
    dplyr::mutate(present = 1)
  
  # For better visualization, we'll create a data frame instead of full matrix expansion
  # This is more efficient for ggplot2
  
  # Get feature order: group by column, then by frequency within column
  feature_order <- heatmap_data |>
    dplyr::count(column, feature, name = "freq") |>
    dplyr::arrange(column, desc(freq)) |>
    dplyr::pull(feature)
  
  heatmap_matrix <- heatmap_matrix |>
    dplyr::mutate(feature = factor(feature, levels = rev(feature_order)))
  
  # Add column information for color coding
  heatmap_matrix <- heatmap_matrix |>
    dplyr::mutate(
      column = sapply(strsplit(as.character(feature), ":"), function(x) x[1]),
      value = sapply(strsplit(as.character(feature), ":"), function(x) ifelse(length(x) > 1, x[2], "")),
      is_unassigned = value == "Unassigned"
    )
  
  # Create color palette: grey for Unassigned, distinct colors per column for others
  unique_cols <- unique(heatmap_matrix$column)
  n_cols <- length(unique_cols)
  
  # Generate distinct colors for each column
  col_colors <- setNames(scales::hue_pal()(n_cols), unique_cols)
  
  # Assign colors: grey for Unassigned, column color for assigned
  heatmap_matrix <- heatmap_matrix |>
    dplyr::mutate(
      fill_color = ifelse(is_unassigned, "grey70", col_colors[column])
    )
  
  # Determine x-axis breaks: show every Nth label to avoid overcrowding
  # For 245 timepoints, show ~12-24 labels (every 10-20 timepoints)
  label_interval <- max(1, floor(n_timepoints / 20))
  x_breaks <- seq(1, n_timepoints, by = label_interval)
  
  # Get labels for those breaks
  label_lookup <- heatmap_matrix |>
    dplyr::distinct(timepoint_id, timepoint_label)
  
  x_labels <- label_lookup$timepoint_label[x_breaks]

  # Create heatmap
  # Make the plotting area smaller and allocate more space to the left for labels
  # Height: smaller overall but keep some spacing per feature
  plot_height_heat <- min(20, max(8, n_features * 0.25))
  # Width: compact to reduce plot area
  plot_width_heat <- 10

  cat("Creating heatmap with dimensions:", plot_width_heat, "x", plot_height_heat, "inches\n")

  p_heat <- ggplot(heatmap_matrix, aes(x = timepoint_id, y = feature, fill = fill_color)) +
    geom_tile(color = NA, linewidth = 0) +  # Very thin tiles, no borders
    scale_fill_identity() +  # Use the colors we assigned directly
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "Timepoint (Date)",
      y = "Feature (Column:Value)",
      title = "Temporal Distribution of Metadata Features"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      # keep x labels readable but small since width is reduced
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      # larger, more readable y labels (take more visual space)
      axis.text.y = element_text(size = 10, hjust = 1, margin = margin(r = 6)),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      # expand left margin so labels aren't clipped and take more space visually
      plot.margin = margin(10, 10, 10, 60)
    )
  
  heatmap_file <- file.path(output_dir, "temporal_feature_heatmap.png")
  
  ggsave(
    filename = heatmap_file,
    plot = p_heat,
    width = plot_width_heat,
    height = plot_height_heat,
    dpi = 300,
    limitsize = FALSE  # Allow large plots
  )
  
  cat("Saved temporal heatmap to:", heatmap_file, "\n")
  
} else {
  cat("Warning: 'Datum' column not found in metadata. Skipping temporal heatmap.\n")
}

cat("\n=== Temporal heatmap section complete ===\n\n")
