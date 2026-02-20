[< Back to Documentation Index](README.md)

# R_shared API Reference

Complete reference for the 8 shared utility modules in `src/R_shared/`. These modules provide the foundational functions used across all cytometry pipeline subprojects.

---

## Dependency Graph

```
io_utils.R (no dependencies)
  └── used by fcs_utils.R (requires io_utils for detect_csv_sep)
gui_utils.R (independent, tcltk optional)
memory_utils.R (independent, Linux-specific)
cluster_utils.R (independent)
plot_utils.R (independent, ggplot2)
flowframe_utils.R (independent, flowCore/Biobase)
stats_utils.R (independent, lme4/dplyr)
```

All modules are sourced via:

```r
shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
```

---

## io_utils.R

I/O utility functions for file management, CSV handling, package loading, and column classification.

**Required packages:** data.table, tibble, dplyr, purrr

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `is_directory` | `is_directory(path)` | Check whether a path points to a directory. | `Logical` |
| `dir_create` | `dir_create(dir_path_list)` | Create multiple directories from a list. Silently succeeds if they already exist. | `NULL` |
| `detect_csv_sep` | `detect_csv_sep(file)` | Auto-detect CSV separator by inspecting the first line of the file. | `";"` or `","` |
| `get_file_size_mb` | `get_file_size_mb(file_path)` | Return the file size in megabytes. | `Numeric` |
| `save_and_cleanup` | `save_and_cleanup(result, output_file, cleanup_vars = NULL)` | Save a data.table to disk, optionally remove variables from the calling environment, and run garbage collection. | `Invisible path` |
| `time_formatter` | `time_formatter(x)` | Format a duration in seconds into a human-readable string. | `Character` (e.g. `"45s"`, `"2min"`, `"1.5h"`) |
| `list_csvs` | `list_csvs(dir_path)` | List all CSV files in a directory with metadata. | `tibble(file_path, file_name, file_size_mb)` |
| `ensure_and_load_packages` | `ensure_and_load_packages(req_packages)` | Install missing packages and load all of them. | `NULL` |
| `remove_columns_from_csvs` | `remove_columns_from_csvs(folder_path, colname_pattern, output_folder = NULL)` | Remove columns matching a regex pattern from every CSV in a folder. Writes results in-place or to `output_folder`. | `List of results` |
| `classify_csv_columns` | `classify_csv_columns(col_names, gating_pattern = "\\(Manual gating\\)$\|\\(Clustering annotated\\)$\|^Clustering raw$")` | Automatically classify Cytolution CSV column names into semantic groups. | `Named list(index, scatter, markers, gating, dimreduction, autofluorescence)` |
| `die` | `die(...)` | Print a fatal error message and terminate the R session with exit code 1. | Never returns |
| `make_outdir` | `make_outdir(path)` | Create a directory if it does not exist. Safe wrapper around `dir.create`. | `NULL` |

### Usage Examples

```r
source(file.path(shared_dir, "io_utils.R"))

# Detect separator and read a CSV
sep <- detect_csv_sep("data/panel_A.csv")
dt <- data.table::fread("data/panel_A.csv", sep = sep)

# List all CSVs in a folder with sizes
csv_info <- list_csvs("data/exported/")
large_files <- csv_info[csv_info$file_size_mb > 100, ]

# Classify columns from a Cytolution export
cols <- classify_csv_columns(names(dt))
marker_data <- dt[, cols$markers, with = FALSE]

# Create output directory and save results
make_outdir("results/cluster_output")
save_and_cleanup(dt, "results/cluster_output/final.csv",
                 cleanup_vars = c("dt", "large_files"))

# Format elapsed time
elapsed <- 3725
message("Completed in ", time_formatter(elapsed))
# "Completed in 1h"

# Fatal error if required file is missing
if (!file.exists("config.R")) die("config.R not found in working directory")
```

---

## gui_utils.R

GUI dialog utilities using tcltk with automatic CLI fallback for headless environments.

**Required packages:** tcltk (optional), purrr

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `is_gui_available` | `is_gui_available()` | Check whether the tcltk graphical toolkit is available in the current session. | `Logical` |
| `tk_multiselect` | `tk_multiselect(choices, title = "Select items", multiple = TRUE, cli_fallback = TRUE)` | Present a multi-select dialog. Falls back to numbered CLI menu when tcltk is unavailable and `cli_fallback = TRUE`. | `Character vector` or `NULL` |
| `select_markers` | `select_markers(marker_info, title = "Select Markers")` | Display a scrollable checkbox dialog for marker selection. `marker_info` must be a list with `indices` and `names` elements. | `Integer vector of indices` or `NULL` |
| `select_markers_by_mode` | `select_markers_by_mode(mode = c("all", "popup", "terminal"), filtered_markers, filtered_indices)` | Dispatch marker selection by mode: return all, show GUI popup, or use terminal input. | `Integer vector` |
| `select_columns_gui` | `select_columns_gui(column_names, title = "Select Columns", multiple = TRUE)` | Convenience wrapper around `tk_multiselect` for column name selection. | `Character vector` |

### Usage Examples

```r
source(file.path(shared_dir, "gui_utils.R"))

# Check if GUI is available before attempting popup
if (is_gui_available()) {
  selected <- tk_multiselect(c("CD3", "CD4", "CD8", "CD19"),
                             title = "Choose markers")
} else {
  message("No GUI available, using all markers")
}

# Marker selection with index tracking
marker_info <- list(indices = 1:5, names = c("CD3", "CD4", "CD8", "CD19", "CD56"))
chosen_idx <- select_markers(marker_info, title = "Select T-cell markers")

# Mode-based dispatch (useful in pipeline configs)
idx <- select_markers_by_mode(mode = "all",
                              filtered_markers = c("CD3", "CD4"),
                              filtered_indices = c(1, 2))
```

---

## memory_utils.R

Memory monitoring and management utilities. Designed for Linux systems using `/proc/meminfo`.

**Required packages:** data.table, stringr, progress (optional)

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `check_memory` | `check_memory(warn_threshold = 0.85)` | Check current memory usage against a threshold. Issues a warning and triggers garbage collection if usage exceeds the threshold. | `Logical` |
| `check_memory_usage` | `check_memory_usage(warn_threshold = 0.8)` | Print a summary of current memory usage (total, used, available) and warn if above threshold. | `Invisible NULL` |
| `read_csv_chunked` | `read_csv_chunked(file_path, chunk_size = 50000, max_memory_gb = 4)` | Read a large CSV file in chunks to avoid exhausting available RAM. Checks memory between chunks and pauses for GC as needed. | `data.table` |
| `get_available_memory_gb` | `get_available_memory_gb()` | Query available system memory from `/proc/meminfo`. Linux-specific. | `Numeric` or `NA_real_` (non-Linux) |

### Usage Examples

```r
source(file.path(shared_dir, "memory_utils.R"))

# Check memory before a large operation
if (check_memory(warn_threshold = 0.9)) {
  message("Memory is getting tight, consider downsampling")
}

# Print a memory status report
check_memory_usage()
# Total: 64.0 GB | Used: 42.3 GB | Available: 21.7 GB

# Load a very large CSV without running out of memory
big_data <- read_csv_chunked("data/all_samples_merged.csv",
                             chunk_size = 100000,
                             max_memory_gb = 8)

# Query available memory programmatically
avail <- get_available_memory_gb()
if (!is.na(avail) && avail < 2) {
  die("Less than 2 GB of memory available. Aborting.")
}
```

---

## cluster_utils.R

Cluster name parsing, label manipulation, and binary-to-categorical conversion utilities.

**Required packages:** none (base R only)

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `clean_cluster_name` | `clean_cluster_name(cluster_name)` | Strip the `"AI_Cluster_"` prefix and `"(Clustering annotated)"` suffix from a cluster name string. | `Character` |
| `extract_cluster_group` | `extract_cluster_group(cluster_name)` | Map a cluster name to a standardized immune cell group (e.g. `"BCs"`, `"CD4TCs"`, `"CD8TCs"`, `"NKCs"`). | `Character` |
| `create_cluster_id` | `create_cluster_id(cluster_name, max_length = 30)` | Generate a filesystem-safe identifier from a cluster name, truncated to `max_length`. | `Character` |
| `flatten_label_list` | `flatten_label_list(label_lists, gating_tree)` | Resolve multi-label assignments to a single label per cell using gating tree depth. Deeper nodes take priority. | `Character vector` |
| `labels_from_binary_df` | `labels_from_binary_df(df)` | Convert a binary data.frame (0/1 columns) to a list of character vectors containing the column names where value equals 1 for each row. | `List of character vectors` |
| `labels_from_binary_dt` | `labels_from_binary_dt(dt)` | Same as `labels_from_binary_df` but for data.table input. Assumes first column is an index column and skips it. | `List of character vectors` |
| `convert_binary_to_single_cluster` | `convert_binary_to_single_cluster(cells_data, cluster_cols)` | Convert a binary cluster assignment matrix into a single `cluster_assignment` column. Handles multi-positive cells by concatenation. | `data.frame` with `cluster_assignment` column |
| `translate_labels` | `translate_labels(string_list, dictionary)` | Apply dictionary-based label translation to a character vector. Unmatched labels are returned unchanged. | `Character vector` |

### Usage Examples

```r
source(file.path(shared_dir, "cluster_utils.R"))

# Clean a raw cluster name from Cytolution
raw <- "AI_Cluster_CD4+ T cells (Clustering annotated)"
clean <- clean_cluster_name(raw)
# "CD4+ T cells"

# Map to standardized group
group <- extract_cluster_group("CD4+ central memory T cells")
# "CD4TCs"

# Create a safe filename for a cluster
id <- create_cluster_id("CD8+ effector memory T cells (re-expressing CD45RA)")
# "CD8_effector_memory_T_cells"

# Convert binary gating matrix to labels
binary_df <- data.frame(CD3 = c(1, 0, 1), CD4 = c(1, 1, 0), CD8 = c(0, 0, 1))
labels <- labels_from_binary_df(binary_df)
# [[1]] c("CD3", "CD4")
# [[2]] c("CD4")
# [[3]] c("CD3", "CD8")

# Translate labels with a dictionary
dict <- c("BCs" = "B cells", "CD4TCs" = "CD4+ T cells", "NKCs" = "NK cells")
translated <- translate_labels(c("BCs", "NKCs", "Unknown"), dict)
# c("B cells", "NK cells", "Unknown")
```

---

## plot_utils.R

Plotting utilities for multi-format saving, barplots, and filename sanitization.

**Required packages:** ggplot2, rlang

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `save_plot_formats` | `save_plot_formats(plot, filename_base, width = 10, height = 8, formats = list(pdf = TRUE, png = TRUE, svg = FALSE), png_dpi = 300)` | Save a ggplot object in multiple formats (PDF, PNG, SVG) based on the `formats` toggle list. | `Character vector` of formats saved |
| `population_barplot` | `population_barplot(population_counts, out_file, title, y_label)` | Create and save a grouped barplot showing population counts across samples or conditions. | `Invisible NULL` |
| `plot_per_sample` | `plot_per_sample(df, value_col, out_dir, out_prefix, y_label, title)` | Generate faceted barplots, one panel per sample, and save to `out_dir`. | `Invisible NULL` |
| `create_safe_filename` | `create_safe_filename(name, max_length = 100)` | Sanitize a string for use as a filename by removing or replacing special characters. | `Character` |
| `create_compact_names` | `create_compact_names(names, max_length = 25)` | Truncate long names and append an ellipsis if they exceed `max_length`. Useful for plot axis labels. | `Character vector` |

### Usage Examples

```r
source(file.path(shared_dir, "plot_utils.R"))

# Save a ggplot in PDF and PNG
library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
saved <- save_plot_formats(p, "results/scatter_wt_mpg",
                           width = 8, height = 6,
                           formats = list(pdf = TRUE, png = TRUE, svg = TRUE))
# saved: c("pdf", "png", "svg")

# Create population barplot
pop_counts <- data.frame(
  sample = rep(c("S1", "S2"), each = 3),
  population = rep(c("CD4", "CD8", "NK"), 2),
  count = c(1200, 800, 400, 1500, 600, 350)
)
population_barplot(pop_counts, "results/pop_barplot", "Cell Counts", "N cells")

# Sanitize filenames from cluster names
safe <- create_safe_filename("CD4+ T cells / memory (activated)")
# "CD4_T_cells_memory_activated"

# Compact long names for axis labels
long_names <- c("CD4+ central memory T cells", "CD8+ effector memory RA+", "NK")
compact <- create_compact_names(long_names, max_length = 20)
# c("CD4+ central memor...", "CD8+ effector memor...", "NK")
```

---

## flowframe_utils.R

FlowFrame and FlowSet manipulation utilities for working with flow cytometry data objects.

**Required packages:** flowCore, Biobase

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `explore_data` | `explore_data(flowframe)` | Extract marker metadata (dye names, marker names, marker-channel mapping) from a flowFrame object. | `list(dyes, markers, mcols)` |
| `compare_panels` | `compare_panels(flowframes_list)` | Compare marker panels across a list of flowFrames to identify unique panel compositions. | `List of unique marker sets` |
| `extract_exprmtx_list` | `extract_exprmtx_list(flow_frames_list, sample_names, markers_of_interest)` | Extract expression matrices from a list of flowFrames, subsetting to specified markers. | `Named list of expression matrices` |
| `build_combined_expr_matrix` | `build_combined_expr_matrix(flow_frames_list, sample_names, markers_of_interest)` | Extract expression matrices from flowFrames, row-bind into a single combined matrix, and generate matching sample labels. | `list(combined_matrix, sample_labels)` |
| `generate_cellID` | `generate_cellID(flowframe, sample_name)` | Generate unique cell identifiers by combining a sample name with sequential cell indices. | `Character vector` (e.g. `"sample_name_cell_1"`) |
| `from_mtx_to_ffr` | `from_mtx_to_ffr(mtx, marker_names = NULL)` | Convert a numeric matrix into a flowFrame object. Optionally set marker names in the parameter description. | `flowFrame` |
| `from_ffl_to_fSet` | `from_ffl_to_fSet(flow_frame_list, sample_names = NULL)` | Convert a list of flowFrame objects into a single flowSet. | `flowSet` |

### Usage Examples

```r
source(file.path(shared_dir, "flowframe_utils.R"))

library(flowCore)

# Explore a flowFrame
ff <- read.FCS("data/sample_001.fcs")
info <- explore_data(ff)
cat("Markers:", info$markers, "\n")
cat("Dyes:", info$dyes, "\n")

# Compare panels across multiple FCS files
fcs_files <- list.files("data/", pattern = "\\.fcs$", full.names = TRUE)
ff_list <- lapply(fcs_files, read.FCS)
panels <- compare_panels(ff_list)
# Returns distinct panel configurations found

# Extract expression matrices for selected markers
expr_list <- extract_exprmtx_list(ff_list,
                                  sample_names = c("S1", "S2", "S3"),
                                  markers_of_interest = c("CD3", "CD4", "CD8"))

# Generate cell IDs for tracking
cell_ids <- generate_cellID(ff, "patient_01")
# c("patient_01_cell_1", "patient_01_cell_2", ...)

# Convert matrix back to flowFrame
mtx <- exprs(ff)
ff_new <- from_mtx_to_ffr(mtx, marker_names = c("CD3", "CD4", "CD8", "CD19"))

# Create a flowSet from individual flowFrames
fset <- from_ffl_to_fSet(ff_list, sample_names = c("S1", "S2", "S3"))
```

---

## fcs_utils.R

FCS file handling utilities including low-level FCS parsing, CSV-FCS mapping, and batch loading.

**Required packages:** flowCore, flowWorkspace, data.table, purrr

> **Note:** This module depends on `io_utils.R`. Always source `io_utils.R` before sourcing `fcs_utils.R`.

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `readFCStext` | `readFCStext(con, offsets, emptyValue = TRUE, cpp = TRUE, ...)` | Parse the TEXT segment of an FCS file from a connection and byte offsets. Low-level function. | `Named character vector` |
| `readFCSgetPar` | `readFCSgetPar(x, pnam, strict = TRUE)` | Retrieve specific FCS parameters from a parsed TEXT segment by parameter name pattern. | `Character vector` |
| `check_tot` | `check_tot(flowframe, index = "")` | Validate that the `$TOT` keyword in an FCS file header matches the actual number of events in the flowFrame. | `Invisible NULL` |
| `create_csv_fcs_mapping` | `create_csv_fcs_mapping(csv_files, fcs_files)` | Create a mapping between CSV filenames and FCS filenames based on string similarity. | `Named character vector` (CSV -> FCS) |
| `get_fcs_event_count` | `get_fcs_event_count(fcs_file)` | Read the event count from an FCS file header without loading the full expression data. | `Integer` |
| `has_manual_gating` | `has_manual_gating(csv_file, pattern = "Manual gating")` | Check whether a CSV file contains columns matching the manual gating pattern. | `Logical` |
| `filter_valid_csv_files` | `filter_valid_csv_files(csv_files, pattern = "Manual gating")` | Filter a vector of CSV file paths to only those containing manual gating columns. | `Character vector` |
| `read_csv_filtered` | `read_csv_filtered(csv_file, pattern = "Manual gating")` | Read a CSV file and return only the index column and columns matching the gating pattern. | `data.table` or `NULL` |
| `load_fcs_list` | `load_fcs_list(fcs_file_list)` | Load a list of FCS file paths into a list of flowFrame objects. | `List of flowFrames` |
| `load_fcs_folder` | `load_fcs_folder(input_file_path)` | Load all FCS files from a directory into a list of flowFrames. | `List of flowFrames` |

### Usage Examples

```r
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "fcs_utils.R"))

# Validate $TOT for a loaded flowFrame
ff <- flowCore::read.FCS("data/sample_001.fcs")
check_tot(ff, index = "sample_001")

# Quick event count without full load
n_events <- get_fcs_event_count("data/sample_001.fcs")
message("Events: ", n_events)

# Map CSVs to their corresponding FCS files
csv_files <- list.files("data/csv/", pattern = "\\.csv$", full.names = TRUE)
fcs_files <- list.files("data/fcs/", pattern = "\\.fcs$", full.names = TRUE)
mapping <- create_csv_fcs_mapping(csv_files, fcs_files)

# Filter CSVs to only those with manual gating
gated_csvs <- filter_valid_csv_files(csv_files, pattern = "Manual gating")

# Read only gating columns from a CSV
gating_dt <- read_csv_filtered("data/csv/sample_001.csv")
# Returns index + manual gating columns only

# Batch load all FCS files from a folder
all_ff <- load_fcs_folder("data/fcs/")
message("Loaded ", length(all_ff), " FCS files")
```

---

## stats_utils.R

Statistical analysis functions including enrichment testing, mixed-effects models, dimensionality reduction, and summary statistics.

**Required packages:** lme4, dplyr, tibble, tidyr, rlang, uwot, performance (optional)

| Name | Signature | Description | Return |
|------|-----------|-------------|--------|
| `calculate_enrichment` | `calculate_enrichment(cluster_col, cells_data, cond_levels)` | Perform Fisher's exact test to determine if a cluster is enriched in one condition versus another. | `tibble(in_A, not_in_A, in_B, not_in_B, odds_ratio, p_value, ci_lower, ci_upper, method)` |
| `do_mixed_effects` | `do_mixed_effects(cluster_col, cells_data, condition_levels)` | Fit a mixed-effects logistic regression for cluster membership as a function of condition, with sample as random effect. | `tibble(estimate, std_error, p_value, odds_ratio, ci_lower, ci_upper, method)` |
| `check_overdispersion` | `check_overdispersion(model)` | Test a fitted model for overdispersion using the ratio of residual deviance to residual degrees of freedom. | `list(ratio, p_value, is_overdispersed)` |
| `fit_nb_glmm` | `fit_nb_glmm(count_data, formula = count ~ condition + (1\|sample_id), max_iter = 20000)` | Fit a negative binomial GLMM for count data. Falls back to Poisson GLMM if negative binomial fails to converge. | Model object or `NULL` |
| `run_pca` | `run_pca(expr_matrix, n_components = 2, scale = TRUE)` | Run Principal Component Analysis on an expression matrix. | `list(coordinates, variance_explained, loadings, pca_object)` |
| `run_umap` | `run_umap(expr_matrix, n_neighbors = 15, min_dist = 0.1, n_components = 2, metric = "euclidean", seed = 42)` | Run UMAP dimensionality reduction on an expression matrix using the uwot package. | `Matrix` of UMAP coordinates (n_cells x n_components) |
| `calculate_statistics` | `calculate_statistics(expr_data, cluster_col, marker_cols)` | Calculate per-cluster, per-marker summary statistics. | `tibble(cluster, marker, mean, median, sd, q25, q75, n)` |
| `get_top_features` | `get_top_features(expr_matrix, n_top = 20, method = "variance")` | Identify the top N features (columns) by variance or mean expression. | `Character vector` |

### Usage Examples

```r
source(file.path(shared_dir, "stats_utils.R"))

# Fisher's exact enrichment test
enrichment <- calculate_enrichment(
  cluster_col = "CD4_Tcells",
  cells_data = cells_df,
  cond_levels = c("Healthy", "Disease")
)
cat("Odds ratio:", enrichment$odds_ratio, "p =", enrichment$p_value, "\n")

# Mixed-effects logistic regression
me_result <- do_mixed_effects(
  cluster_col = "CD4_Tcells",
  cells_data = cells_df,
  condition_levels = c("Healthy", "Disease")
)

# Negative binomial GLMM for count data
count_df <- data.frame(
  count = c(120, 85, 200, 150, 90, 180),
  condition = rep(c("Ctrl", "Treated"), each = 3),
  sample_id = paste0("S", 1:6)
)
model <- fit_nb_glmm(count_df, formula = count ~ condition + (1 | sample_id))

# Run PCA
pca_res <- run_pca(expr_matrix, n_components = 3, scale = TRUE)
cat("Variance explained:", pca_res$variance_explained, "\n")

# Run UMAP
umap_coords <- run_umap(expr_matrix, n_neighbors = 30, min_dist = 0.2, seed = 123)
plot(umap_coords[, 1], umap_coords[, 2], pch = 16, cex = 0.3)

# Per-cluster statistics
stats <- calculate_statistics(expr_data, cluster_col = "cluster", marker_cols = c("CD3", "CD4"))
# Returns: cluster | marker | mean | median | sd | q25 | q75 | n

# Select top variable markers
top_markers <- get_top_features(expr_matrix, n_top = 15, method = "variance")
```

---

## Quick Reference: Module Selection Guide

| Task | Module | Key Function |
|------|--------|-------------|
| Read/write CSV files | `io_utils.R` | `detect_csv_sep`, `save_and_cleanup` |
| Classify CSV columns | `io_utils.R` | `classify_csv_columns` |
| Interactive marker selection | `gui_utils.R` | `tk_multiselect`, `select_markers` |
| Monitor RAM usage | `memory_utils.R` | `check_memory`, `get_available_memory_gb` |
| Load large files safely | `memory_utils.R` | `read_csv_chunked` |
| Parse cluster names | `cluster_utils.R` | `clean_cluster_name`, `extract_cluster_group` |
| Convert binary gating | `cluster_utils.R` | `labels_from_binary_df`, `convert_binary_to_single_cluster` |
| Save plots multi-format | `plot_utils.R` | `save_plot_formats` |
| Inspect FCS marker panels | `flowframe_utils.R` | `explore_data`, `compare_panels` |
| Convert matrices to flowFrames | `flowframe_utils.R` | `from_mtx_to_ffr` |
| Validate FCS headers | `fcs_utils.R` | `check_tot`, `get_fcs_event_count` |
| Map CSV to FCS files | `fcs_utils.R` | `create_csv_fcs_mapping` |
| Cluster enrichment testing | `stats_utils.R` | `calculate_enrichment`, `do_mixed_effects` |
| Dimensionality reduction | `stats_utils.R` | `run_pca`, `run_umap` |
| Marker summary statistics | `stats_utils.R` | `calculate_statistics`, `get_top_features` |
