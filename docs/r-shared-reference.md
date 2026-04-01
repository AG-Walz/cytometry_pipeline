[< Back to Documentation Index](README.md)

# R_shared API Reference

Complete reference for the 12 shared utility modules in `src/R_shared/`. These modules provide the foundational functions used across all cytometry pipeline subprojects.

---

## Table of Contents

- [Dependency Graph](#dependency-graph)
- [Sourcing Pattern](#sourcing-pattern)
- [Module Reference](#module-reference)
  - [io_utils.R](#io_utilsr) -- File I/O, CSV handling, directory management
  - [gui_utils.R](#gui_utilsr) -- Interactive selection dialogs (tcltk / CLI)
  - [memory_utils.R](#memory_utilsr) -- Memory monitoring and chunked CSV reading
  - [cluster_utils.R](#cluster_utilsr) -- Cluster name parsing and label manipulation
  - [plot_utils.R](#plot_utilsr) -- Multi-format plot saving, barplots, filename sanitization
  - [flowframe_utils.R](#flowframe_utilsr) -- FlowFrame/FlowSet creation and conversion
  - [fcs_utils.R](#fcs_utilsr) -- FCS file parsing, validation, and batch loading
  - [stats_utils.R](#stats_utilsr) -- Statistical tests, GLMM, PCA, UMAP, summary statistics
  - [clustering_dispatch.R](#clustering_dispatchr) -- Clustering method dispatcher
  - [dimred_dispatch.R](#dimred_dispatchr) -- Dimensionality reduction dispatcher
  - [batch_correction_dispatch.R](#batch_correction_dispatchr) -- Batch correction dispatcher
  - [pipeline_utils.R](#pipeline_utilsr) -- Pipeline initialization and parameter builders
- [Quick Reference: Module Selection Guide](#quick-reference-module-selection-guide)

---

## Dependency Graph

```
io_utils.R (no dependencies)
  +-- fcs_utils.R (requires io_utils.R for detect_csv_sep)
gui_utils.R (independent; tcltk optional)
memory_utils.R (independent; cross-platform Linux/macOS/Windows)
cluster_utils.R (independent)
plot_utils.R (independent; requires ggplot2)
flowframe_utils.R (independent; requires flowCore, Biobase)
stats_utils.R (independent; requires lme4, dplyr, uwot)
clustering_dispatch.R (sources clustering/*.R at runtime)
dimred_dispatch.R (sources visualization/*.R at runtime)
batch_correction_dispatch.R (sources harmonization/*.R and flowframe_utils.R at runtime)
pipeline_utils.R (sources io_utils.R and memory_utils.R at runtime)
```

---

## Sourcing Pattern

All modules are sourced via the standard pattern:

```r
shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
# ... source only what you need
```

---

## Module Reference

---

### io_utils.R

I/O utility functions for file management, CSV handling, package loading, and column classification.

**Required packages:** data.table, tibble, dplyr, purrr

**Used by:** All subprojects (loading, preprocessing, harmonization, clustering, analysis, visualization, metadata_qc)

#### `is_directory`

Check whether a path points to an existing directory.

```r
is_directory(path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | character | File system path to check |

**Returns:** `logical` -- `TRUE` if `path` exists and is a directory, `FALSE` otherwise.

```r
if (is_directory("data/output")) cat("Output dir exists\n")
```

---

#### `dir_create`

Create multiple directories from a character vector of paths. Silently succeeds if directories already exist.

```r
dir_create(dir_path_list)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `dir_path_list` | character vector | Directory paths to create |

**Returns:** `NULL` (invisible). Side effect: directories created on disk.

```r
dir_create(c("results/plots", "results/tables", "results/models"))
```

---

#### `detect_csv_sep`

Auto-detect CSV separator by inspecting the first line of the file. Returns `";"` if a semicolon is found, otherwise `","`.

```r
detect_csv_sep(file)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `file` | character | Path to the CSV file |

**Returns:** `character` -- `";"` or `","`.

```r
sep <- detect_csv_sep("data/panel_A.csv")
dt <- data.table::fread("data/panel_A.csv", sep = sep)
```

---

#### `get_file_size_mb`

Return the file size in megabytes.

```r
get_file_size_mb(file_path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `file_path` | character | Path to the file |

**Returns:** `numeric` -- file size in MB.

```r
size <- get_file_size_mb("data/large_export.csv")
if (size > 500) message("Large file: ", round(size), " MB")
```

---

#### `save_and_cleanup`

Save a data.table to a semicolon-delimited CSV, optionally remove variables from the calling environment, and trigger garbage collection.

```r
save_and_cleanup(result, output_file, cleanup_vars = NULL)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `result` | data.table | The data to write |
| `output_file` | character | Output file path |
| `cleanup_vars` | character vector or `NULL` | Variable names to `rm()` from the parent frame |

**Returns:** `character` (invisible) -- the `output_file` path.

```r
save_and_cleanup(dt, "results/final.csv", cleanup_vars = c("dt", "raw_data"))
```

---

#### `time_formatter`

Format a numeric vector of durations in seconds into human-readable strings.

```r
time_formatter(x)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `x` | numeric vector | Duration values in seconds |

**Returns:** `character vector` -- e.g. `"45s"`, `"2min"`, `"1.5h"`.

```r
elapsed <- 3725
message("Completed in ", time_formatter(elapsed))
# "Completed in 1.0h"
```

---

#### `list_csvs`

List all `.csv` files in a directory with file metadata. Stops with an error if no CSV files are found.

```r
list_csvs(dir_path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `dir_path` | character | Directory to scan |

**Returns:** `tibble` with columns `file_path`, `file_name`, `file_size_mb`, sorted by `file_name`.

```r
csv_info <- list_csvs("data/exported/")
large_files <- csv_info[csv_info$file_size_mb > 100, ]
```

---

#### `ensure_and_load_packages`

Install any missing packages from CRAN, then load all packages.

```r
ensure_and_load_packages(req_packages)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `req_packages` | character vector | Package names to install and load |

**Returns:** `NULL` (invisible). Side effect: packages installed and attached.

```r
ensure_and_load_packages(c("ggplot2", "data.table", "uwot"))
```

---

#### `remove_columns_from_csvs`

Remove columns matching a regex pattern from every CSV file in a folder. Writes filtered results to `output_folder` (defaults to a `filtered/` subdirectory).

```r
remove_columns_from_csvs(folder_path, colname_pattern, output_folder = NULL)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `folder_path` | character | Path to folder containing CSV files |
| `colname_pattern` | character | Regex pattern matching column names to remove |
| `output_folder` | character or `NULL` | Output directory (default: `folder_path/filtered`) |

**Returns:** `list` (invisible) -- named list (one entry per file) with `original_cols`, `filtered_cols`, `removed_cols`.

```r
remove_columns_from_csvs("data/exports/", colname_pattern = "^cyto_2d_")
```

---

#### `classify_csv_columns`

Automatically classify Cytolution CSV column names into semantic groups based on known naming patterns.

```r
classify_csv_columns(col_names, gating_pattern = "\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `col_names` | character vector | Column names from a Cytolution CSV export |
| `gating_pattern` | character | Regex pattern identifying gating/clustering columns |

**Returns:** Named `list` with character vectors: `index`, `scatter`, `markers`, `gating`, `dimreduction`, `autofluorescence`.

```r
cols <- classify_csv_columns(names(dt))
marker_data <- dt[, cols$markers, with = FALSE]
```

---

#### `die`

Print a fatal error message and terminate the R session with exit code 1. Used in scripts for unrecoverable errors.

```r
die(...)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `...` | any | Message parts passed to `message()` |

**Returns:** Never returns (calls `quit(status = 1)`).

```r
if (!file.exists("config.R")) die("config.R not found in working directory")
```

---

#### `make_outdir`

Create a directory if it does not exist. Safe wrapper around `dir.create(recursive = TRUE)`.

```r
make_outdir(path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | character | Directory path to create |

**Returns:** `NULL` (invisible). Side effect: directory created on disk.

```r
make_outdir("results/cluster_output")
```

---

### gui_utils.R

GUI dialog utilities using tcltk with automatic CLI fallback for headless environments.

**Required packages:** tcltk (optional), purrr

**Used by:** loading, analysis, clustering

#### `is_gui_available`

Check whether the tcltk graphical toolkit is available in the current session. Returns `TRUE` only if the session is interactive and tcltk can be loaded.

```r
is_gui_available()
```

**Returns:** `logical`.

```r
if (is_gui_available()) {
  selected <- tk_multiselect(markers, title = "Choose markers")
}
```

---

#### `tk_multiselect`

Present a multi-select dialog. Shows a tcltk listbox when GUI is available; falls back to a numbered CLI menu when `cli_fallback = TRUE` and no GUI is present.

```r
tk_multiselect(choices, title = "Select items", multiple = TRUE, cli_fallback = TRUE)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `choices` | character vector | Items to choose from |
| `title` | character | Dialog window title |
| `multiple` | logical | Allow multiple selections (`TRUE`) or single selection (`FALSE`) |
| `cli_fallback` | logical | Use CLI prompts when GUI is unavailable |

**Returns:** `character vector` of selected items, or `NULL` if cancelled.

```r
selected <- tk_multiselect(c("CD3", "CD4", "CD8", "CD19"), title = "Choose markers")
```

---

#### `select_markers`

Display a scrollable checkbox dialog for marker selection. The `marker_info` argument must be a list with `indices` (integer vector) and `names` (character vector) elements.

```r
select_markers(marker_info, title = "Select Markers")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `marker_info` | list | List with `indices` (integer vector) and `names` (character vector) |
| `title` | character | Dialog title |

**Returns:** `integer vector` of selected marker indices, or `NULL` if cancelled.

```r
marker_info <- list(indices = 1:5, names = c("CD3", "CD4", "CD8", "CD19", "CD56"))
chosen_idx <- select_markers(marker_info, title = "Select T-cell markers")
```

---

#### `select_markers_by_mode`

Dispatch marker selection by mode: return all markers immediately, show a GUI popup, or prompt in the terminal.

```r
select_markers_by_mode(mode = c("all", "popup", "terminal"), filtered_markers, filtered_indices)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `mode` | character | One of `"all"`, `"popup"`, `"terminal"` |
| `filtered_markers` | character vector | Marker names available for selection |
| `filtered_indices` | integer vector | Column indices corresponding to `filtered_markers` |

**Returns:** `integer vector` of selected marker column indices.

```r
idx <- select_markers_by_mode(mode = "all",
                              filtered_markers = c("CD3", "CD4"),
                              filtered_indices = c(1, 2))
```

---

#### `select_columns_gui`

Convenience wrapper around `tk_multiselect` for column name selection.

```r
select_columns_gui(column_names, title = "Select Columns", multiple = TRUE)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `column_names` | character vector | Available column names |
| `title` | character | Dialog title |
| `multiple` | logical | Allow multiple selections |

**Returns:** `character vector` of selected column names.

```r
cols <- select_columns_gui(names(dt), title = "Select analysis columns")
```

---

### memory_utils.R

Memory monitoring and management utilities. Cross-platform: supports Linux, macOS, and Windows.

**Required packages:** data.table, stringr, progress (optional)

**Used by:** analysis, loading

#### `check_memory`

Check current memory usage against a threshold. Issues a warning and triggers garbage collection if usage exceeds the threshold.

```r
check_memory(warn_threshold = 0.85)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `warn_threshold` | numeric | Fraction (0--1) at which to warn and trigger GC |

**Returns:** `logical` -- `TRUE` if memory usage is OK, `FALSE` if under pressure.

```r
if (!check_memory(warn_threshold = 0.9)) {
  message("Memory pressure detected, consider downsampling")
}
```

---

#### `check_memory_usage`

Print a formatted summary of current system memory usage (used / total MB, percentage) and warn if above threshold.

```r
check_memory_usage(warn_threshold = 0.8)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `warn_threshold` | numeric | Fraction (0--1) at which to trigger GC |

**Returns:** `invisible(NULL)`. Side effect: prints memory summary to console.

```r
check_memory_usage()
# Memory: 43264 / 65536 MB (66.0%)
```

---

#### `read_csv_chunked`

Read a large CSV file in chunks to avoid exhausting available RAM. Only engages chunked reading when the file size exceeds `max_memory_gb * 1000` MB; otherwise falls back to a simple `data.table::fread()`.

```r
read_csv_chunked(file_path, chunk_size = 50000, max_memory_gb = 4)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `file_path` | character | Path to the CSV file |
| `chunk_size` | integer | Number of rows per chunk |
| `max_memory_gb` | numeric | File size threshold (in GB) above which chunked reading is used |

**Returns:** `data.table` with the full file contents.

```r
big_data <- read_csv_chunked("data/all_samples_merged.csv",
                             chunk_size = 100000, max_memory_gb = 8)
```

---

#### `get_available_memory_gb`

Query available system memory. Cross-platform (Linux via `free`, macOS via `vm_stat`, Windows via `wmic`).

```r
get_available_memory_gb()
```

**Returns:** `numeric` -- available memory in GB, or `NA_real_` if it cannot be determined.

```r
avail <- get_available_memory_gb()
if (!is.na(avail) && avail < 2) die("Less than 2 GB available. Aborting.")
```

---

### cluster_utils.R

Cluster name parsing, label manipulation, and binary-to-categorical conversion utilities.

**Required packages:** none (base R only)

**Used by:** analysis, clustering, loading

#### `clean_cluster_name`

Strip the `"AI_Cluster_"` prefix and `"(Clustering annotated)"` suffix from a cluster name string.

```r
clean_cluster_name(cluster_name)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cluster_name` | character | Raw cluster name string |

**Returns:** `character` -- cleaned cluster name, or `"Unknown"` if input is `NA` or empty.

```r
clean_cluster_name("AI_Cluster_CD4+ T cells (Clustering annotated)")
# "CD4+ T cells"
```

---

#### `extract_cluster_group`

Map a cluster name to a standardized immune cell group label. Recognizes patterns for B cells, CD3-, CD4 T cells, CD8 T cells, NK cells, NKT cells, gamma-delta T cells, T cells, and lymphoid cells.

```r
extract_cluster_group(cluster_name)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cluster_name` | character | Cluster name to classify |

**Returns:** `character` -- standardized group label (e.g. `"BCs"`, `"CD4TCs"`, `"NKCs"`), or the cleaned name if no pattern matches.

```r
extract_cluster_group("CD4TCs_memory_activated")
# "CD4TCs"
```

---

#### `create_cluster_id`

Generate a filesystem-safe identifier from a cluster name. Replaces non-alphanumeric characters with underscores and truncates to `max_length`.

```r
create_cluster_id(cluster_name, max_length = 30)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cluster_name` | character | Full cluster name |
| `max_length` | integer | Maximum length for the ID |

**Returns:** `character` -- safe identifier string.

```r
create_cluster_id("CD8+ effector memory T cells (re-expressing CD45RA)")
# "CD8__effector_memory_T_cells_"
```

---

#### `flatten_label_list`

Resolve multi-label assignments to a single label per cell using gating tree depth. When a cell has multiple labels, the label with the deepest position in the gating tree wins.

```r
flatten_label_list(label_lists, gating_tree)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `label_lists` | list | List where each element is a character vector of labels |
| `gating_tree` | named numeric vector | Label depths (names = labels, values = depth) |

**Returns:** `character vector` -- single label per cell.

```r
labels <- list(c("Lymphocytes", "CD3+"), c("Lymphocytes"), character(0))
tree <- c("Lymphocytes" = 1, "CD3+" = 2)
flatten_label_list(labels, tree)
# c("CD3+", "Lymphocytes", "root")
```

---

#### `labels_from_binary_df`

Convert a binary (0/1) data.frame to a list of character vectors containing the column names where value equals 1 for each row.

```r
labels_from_binary_df(df)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `df` | data.frame | Data frame with binary (0/1) columns |

**Returns:** `list` of character vectors.

```r
binary_df <- data.frame(CD3 = c(1, 0, 1), CD4 = c(1, 1, 0))
labels_from_binary_df(binary_df)
# [[1]] c("CD3", "CD4")  [[2]] c("CD4")  [[3]] c("CD3")
```

---

#### `labels_from_binary_dt`

Same as `labels_from_binary_df` but for `data.table` input. Assumes the first column is an index column and skips it.

```r
labels_from_binary_dt(dt)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `dt` | data.table | data.table with first column as index, remaining as binary |

**Returns:** `list` of character vectors.

```r
dt <- data.table::data.table(index = 1:3, CD3 = c(1, 0, 1), CD4 = c(1, 1, 0))
labels_from_binary_dt(dt)
```

---

#### `convert_binary_to_single_cluster`

Convert binary cluster columns into a single `cluster_assignment` column. For multi-positive cells, the first matching group is used.

```r
convert_binary_to_single_cluster(cells_data, cluster_cols)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cells_data` | data.frame | Data frame with binary cluster columns |
| `cluster_cols` | character vector | Names of binary cluster columns |

**Returns:** `data.frame` -- input with an added `cluster_assignment` column. Unassigned cells are labeled `"Unassigned"`.

```r
cells_data <- convert_binary_to_single_cluster(cells_df, c("CD4TCs", "CD8TCs", "NKCs"))
table(cells_data$cluster_assignment)
```

---

#### `translate_labels`

Apply dictionary-based label translation to a character vector. Labels not found in the dictionary are returned unchanged.

```r
translate_labels(string_list, dictionary)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `string_list` | character vector | Labels to translate |
| `dictionary` | named character vector | Mapping (names = original, values = standardized) |

**Returns:** `character vector` -- translated labels.

```r
dict <- c("BCs" = "B cells", "CD4TCs" = "CD4+ T cells", "NKCs" = "NK cells")
translate_labels(c("BCs", "NKCs", "Unknown"), dict)
# c("B cells", "NK cells", "Unknown")
```

---

### plot_utils.R

Plotting utilities for multi-format saving, grouped barplots, and filename sanitization.

**Required packages:** ggplot2, rlang, grDevices

**Used by:** visualization, analysis, metadata_qc

#### `save_plot_formats`

Save a ggplot object in multiple formats (PDF, PNG, SVG) based on the `formats` toggle list.

```r
save_plot_formats(plot, filename_base, width = 14, height = 10,
                  formats = list(pdf = TRUE, png = TRUE, svg = TRUE),
                  png_dpi = 200)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `plot` | ggplot | The ggplot object to save |
| `filename_base` | character | Base filename without extension |
| `width` | numeric | Plot width in inches |
| `height` | numeric | Plot height in inches |
| `formats` | named logical list | Which formats to produce (keys: `pdf`, `png`, `svg`) |
| `png_dpi` | integer | DPI for PNG output |

**Returns:** `character vector` -- uppercase format names that were saved (e.g. `c("PDF", "PNG")`).

```r
p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
save_plot_formats(p, "results/scatter", formats = list(pdf = TRUE, png = TRUE, svg = FALSE))
```

---

#### `save_base_plot_formats`

Save a base R plot (passed as a zero-argument function) in PNG and SVG formats. Use this instead of `save_plot_formats` when the plot is not a ggplot object.

```r
save_base_plot_formats(filename_base, plot_fn, width = 14, height = 10, png_dpi = 200)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `filename_base` | character | Base filename without extension |
| `plot_fn` | function | Zero-argument function that draws the plot |
| `width` | numeric | Plot width in inches |
| `height` | numeric | Plot height in inches |
| `png_dpi` | integer | DPI for PNG output |

**Returns:** `NULL` (invisible). Side effect: writes `{filename_base}.png` and `{filename_base}.svg`.

```r
save_base_plot_formats("results/heatmap", function() {
  heatmap(as.matrix(expr_data), scale = "column")
})
```

---

#### `population_barplot`

Create and save a grouped barplot showing population counts across samples or conditions.

```r
population_barplot(population_counts, out_file, title, y_label)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `population_counts` | data.frame | Must have `Population`, `Count`, and `Sample` columns |
| `out_file` | character | Output file path (extension stripped for multi-format saving) |
| `title` | character | Plot title |
| `y_label` | character | Y-axis label |

**Returns:** `NULL` (invisible). Side effect: plot saved via `save_plot_formats()`.

```r
pop_df <- data.frame(Population = c("CD4", "CD8"), Count = c(1200, 800), Sample = c("S1", "S1"))
population_barplot(pop_df, "results/pop_barplot.pdf", "Cell Counts", "N cells")
```

---

#### `plot_per_sample`

Generate one barplot per sample (faceted by `Workspace`) and save each to disk.

```r
plot_per_sample(df, value_col, out_dir, out_prefix, y_label, title)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `df` | data.frame | Must have `Population`, `Sample`, and `Workspace` columns |
| `value_col` | character | Name of the column to plot on the Y-axis |
| `out_dir` | character | Output directory |
| `out_prefix` | character | Filename prefix |
| `y_label` | character | Y-axis label |
| `title` | character | Plot title prefix (sample name appended) |

**Returns:** `NULL` (invisible). Side effect: one multi-format plot per unique sample.

```r
plot_per_sample(comparison_df, "count", "results/per_sample/", "comparison", "N", "Counts -")
```

---

#### `create_safe_filename`

Sanitize a string for use as a filename. Replaces special characters with underscores and collapses consecutive underscores.

```r
create_safe_filename(name, max_length = 100)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `name` | character | String to sanitize |
| `max_length` | integer | Maximum filename length |

**Returns:** `character` -- sanitized filename.

```r
create_safe_filename("CD4+ T cells / memory (activated)")
# "CD4_T_cells_memory_activated"
```

---

#### `create_compact_names`

Truncate long names and append `"..."` if they exceed `max_length`. Useful for plot axis labels.

```r
create_compact_names(names, max_length = 25)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `names` | character vector | Names to compact |
| `max_length` | integer | Maximum display length |

**Returns:** `character vector` -- compacted names.

```r
create_compact_names(c("CD4+ central memory T cells", "NK"), max_length = 20)
# c("CD4+ central memor...", "NK")
```

---

### flowframe_utils.R

FlowFrame and FlowSet manipulation utilities for working with flow cytometry data objects.

**Required packages:** flowCore, Biobase

**Used by:** loading, preprocessing, harmonization, batch_correction_dispatch

#### `explore_data`

Extract marker metadata from a flowFrame: dye (channel) names, marker (antigen) names, and the column indices of channels that have marker annotations.

```r
explore_data(flowframe)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `flowframe` | flowFrame | A flowCore flowFrame object |

**Returns:** `list` with `dyes` (character vector), `markers` (character vector), `mcols` (integer vector of marker column indices).

```r
ff <- flowCore::read.FCS("data/sample_001.fcs")
info <- explore_data(ff)
cat("Markers:", info$markers[info$mcols], "\n")
```

---

#### `compare_panels`

Compare marker panels across a list of flowFrames to identify unique panel compositions.

```r
compare_panels(flowframes_list)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `flowframes_list` | list | List of flowFrame objects |

**Returns:** `list` -- each element is a character vector representing a unique marker panel.

```r
panels <- compare_panels(ff_list)
cat("Found", length(panels), "distinct panels\n")
```

---

#### `extract_exprmtx_list`

Extract expression matrices from a list of flowFrames, subsetting to specified markers of interest.

```r
extract_exprmtx_list(flow_frames_list, sample_names, markers_of_interest)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `flow_frames_list` | list | List of flowFrame objects |
| `sample_names` | character vector | Sample names (one per flowFrame) |
| `markers_of_interest` | character vector | Marker names to extract |

**Returns:** Named `list` of numeric matrices (one per sample, rows = cells, columns = markers).

```r
expr_list <- extract_exprmtx_list(ff_list, c("S1", "S2"), c("CD3", "CD4", "CD8"))
```

---

#### `build_combined_expr_matrix`

Extract expression matrices from flowFrames, row-bind into a single combined matrix, and generate a matching sample-label vector. This is a key bridge function connecting flowFrame-based modules to matrix-based modules.

```r
build_combined_expr_matrix(flow_frames_list, sample_names, markers_of_interest)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `flow_frames_list` | list | List of flowFrame objects |
| `sample_names` | character vector | Sample names (one per flowFrame) |
| `markers_of_interest` | character vector | Marker names to extract |

**Returns:** `list` with `combined_matrix` (numeric matrix, cells x markers) and `sample_labels` (character vector, one per row).

```r
result <- build_combined_expr_matrix(ff_list, sample_names, markers)
dim(result$combined_matrix)  # [n_cells, n_markers]
table(result$sample_labels)  # cells per sample
```

---

#### `generate_cellID`

Generate unique cell identifiers by combining a sample name with sequential cell indices.

```r
generate_cellID(flowframe, sample_name)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `flowframe` | flowFrame | A flowFrame object |
| `sample_name` | character | Sample name prefix |

**Returns:** `character vector` -- e.g. `c("S1_cell_1", "S1_cell_2", ...)`.

```r
cell_ids <- generate_cellID(ff, "patient_01")
```

---

#### `from_mtx_to_ffr`

Convert a numeric matrix into a flowFrame object. Creates minimal FCS-compliant annotation data with range parameters.

```r
from_mtx_to_ffr(mtx, marker_names = NULL)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `mtx` | numeric matrix | Expression matrix (rows = cells, columns = markers) |
| `marker_names` | character vector or `NULL` | Optional marker names to set as column names |

**Returns:** `flowFrame` object.

```r
ff_new <- from_mtx_to_ffr(expr_matrix, marker_names = c("CD3", "CD4", "CD8"))
```

---

#### `from_ffl_to_fSet`

Convert a list of flowFrame objects into a single flowSet.

```r
from_ffl_to_fSet(flow_frame_list, sample_names = NULL)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `flow_frame_list` | list | List of flowFrame objects |
| `sample_names` | character vector or `NULL` | Optional sample names |

**Returns:** `flowSet` object.

```r
fset <- from_ffl_to_fSet(ff_list, sample_names = c("S1", "S2", "S3"))
```

---

### fcs_utils.R

FCS file handling utilities including low-level FCS parsing, CSV-FCS mapping, and batch loading.

**Required packages:** flowCore, flowWorkspace, data.table, purrr

**Dependency:** This module requires `io_utils.R` to be sourced first (for `detect_csv_sep()`).

**Used by:** loading

#### `readFCStext`

Parse the TEXT segment of an FCS file from a raw connection and byte offsets. Low-level function; most users should use `flowCore::read.FCS()` instead.

```r
readFCStext(con, offsets, emptyValue = TRUE, cpp = TRUE, ...)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `con` | connection | Open connection to the FCS file |
| `offsets` | named numeric vector | Must include `textstart`, `textend`, `FCSversion` |
| `emptyValue` | logical | Include empty values in output |
| `cpp` | logical | Use the C++ parser if available |

**Returns:** Named `character vector` of text segment key-value pairs (FCS keywords).

```r
con <- file("data/sample.fcs", open = "rb")
offsets <- c(textstart = 256, textend = 4096, FCSversion = 3)
txt <- readFCStext(con, offsets)
close(con)
```

---

#### `readFCSgetPar`

Retrieve specific parameters from a parsed FCS TEXT segment by parameter name.

```r
readFCSgetPar(x, pnam, strict = TRUE)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `x` | named character vector | Parsed FCS text segment (from `readFCStext()`) |
| `pnam` | character vector | Parameter names to retrieve |
| `strict` | logical | If `TRUE`, stop on missing parameters; if `FALSE`, return `NA` for missing |

**Returns:** `character vector` of parameter values.

```r
tot <- readFCSgetPar(txt, "$TOT")
```

---

#### `check_tot`

Validate that the `$TOT` keyword in an FCS file header matches the actual number of events in the flowFrame.

```r
check_tot(flowframe, index = "")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `flowframe` | flowFrame | A loaded flowFrame object |
| `index` | character | Optional sample identifier for log messages |

**Returns:** `NULL` (invisible). Side effect: prints confirmation message or issues a warning on mismatch.

```r
ff <- flowCore::read.FCS("data/sample_001.fcs")
check_tot(ff, index = "sample_001")
```

---

#### `create_csv_fcs_mapping`

Create a mapping between CSV filenames and FCS filenames based on substring matching. Normalizes filenames by removing hyphens, underscores, and extensions before matching.

```r
create_csv_fcs_mapping(csv_files, fcs_files)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `csv_files` | character vector | Full paths to CSV files |
| `fcs_files` | character vector | Full paths to FCS files |

**Returns:** Named `character vector` (names = CSV paths, values = matched FCS paths). Stops with an error if any CSV cannot be uniquely matched.

```r
mapping <- create_csv_fcs_mapping(csv_files, fcs_files)
fcs_for_csv <- mapping[csv_files[1]]
```

---

#### `get_fcs_event_count`

Read the event count from an FCS file header without loading the full expression data. Falls back to loading the full file if the `$TOT` keyword is unavailable.

```r
get_fcs_event_count(fcs_file)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `fcs_file` | character | Path to FCS file |

**Returns:** `integer` -- event (cell) count.

```r
n_events <- get_fcs_event_count("data/sample_001.fcs")
message("Events: ", n_events)
```

---

#### `has_manual_gating`

Check whether a CSV file contains columns matching a gating pattern.

```r
has_manual_gating(csv_file, pattern = "Manual gating")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `csv_file` | character | Path to CSV file |
| `pattern` | character | Regex pattern for gating columns |

**Returns:** `logical`.

```r
if (has_manual_gating("data/sample.csv")) cat("Gating columns found\n")
```

---

#### `filter_valid_csv_files`

Filter a vector of CSV file paths to only those containing columns matching a gating pattern.

```r
filter_valid_csv_files(csv_files, pattern = "Manual gating")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `csv_files` | character vector | CSV file paths to filter |
| `pattern` | character | Regex pattern for gating columns |

**Returns:** `character vector` -- filtered file paths.

```r
gated_csvs <- filter_valid_csv_files(csv_files, pattern = "Manual gating")
```

---

#### `read_csv_filtered`

Read a CSV file and return only the index column and columns matching the gating pattern.

```r
read_csv_filtered(csv_file, pattern = "Manual gating")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `csv_file` | character | Path to CSV file |
| `pattern` | character | Regex pattern for columns to keep |

**Returns:** `data.table` with filtered columns, or `NULL` if no matching columns found.

```r
gating_dt <- read_csv_filtered("data/sample_001.csv")
# Returns index + manual gating columns only
```

---

#### `load_fcs_list`

Load a vector of FCS file paths into a list of flowFrame objects. Uses `flowWorkspace::load_cytoframe_from_fcs()` with `transformation = FALSE`.

```r
load_fcs_list(fcs_file_list)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `fcs_file_list` | character vector | Full paths to FCS files |

**Returns:** `list` of flowFrame objects.

```r
fcs_paths <- list.files("data/fcs/", pattern = "\\.fcs$", full.names = TRUE)
ff_list <- load_fcs_list(fcs_paths)
```

---

#### `load_fcs_folder`

Load all FCS files from a directory into a list of flowFrames. Derives sample names from the first two underscore-delimited parts of each filename.

```r
load_fcs_folder(input_file_path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `input_file_path` | character | Path to folder containing FCS files |

**Returns:** `list` of flowFrame objects.

```r
all_ff <- load_fcs_folder("data/fcs/")
message("Loaded ", length(all_ff), " FCS files")
```

---

### stats_utils.R

Statistical analysis functions including enrichment testing, mixed-effects models, dimensionality reduction, and summary statistics.

**Required packages:** lme4, dplyr, tibble, tidyr, rlang, uwot, performance (optional)

**Used by:** analysis, visualization

#### `calculate_enrichment`

Perform Fisher's exact test to determine if a cluster is enriched in one condition versus another. Constructs a 2x2 contingency table from binary cluster membership and condition labels.

```r
calculate_enrichment(cluster_col, cells_data, cond_levels)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cluster_col` | character | Column name with binary (0/1) cluster membership |
| `cells_data` | data.frame | Must contain `cluster_col` and a `condition` column |
| `cond_levels` | character vector (length 2) | The two condition levels to compare |

**Returns:** `tibble` with columns `in_A`, `not_in_A`, `in_B`, `not_in_B`, `odds_ratio`, `p_value`, `ci_lower`, `ci_upper`, `method`.

```r
enrichment <- calculate_enrichment("CD4_Tcells", cells_df, c("Healthy", "Disease"))
cat("OR:", enrichment$odds_ratio, "p =", enrichment$p_value, "\n")
```

---

#### `do_mixed_effects`

Fit a mixed-effects logistic regression (`cluster ~ condition + (1 | file_name)`) for cluster membership as a function of condition, with sample as random intercept. Uses bobyqa optimizer.

```r
do_mixed_effects(cluster_col, cells_data, condition_levels)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cluster_col` | character | Column name with binary cluster membership |
| `cells_data` | data.frame | Must contain `cluster_col`, `condition`, and `file_name` columns |
| `condition_levels` | character vector (length 2) | The two condition levels |

**Returns:** `tibble` with columns `estimate`, `std_error`, `p_value`, `odds_ratio`, `ci_lower`, `ci_upper`, `method`.

```r
me_result <- do_mixed_effects("CD4_Tcells", cells_df, c("Healthy", "Disease"))
```

---

#### `check_overdispersion`

Test a fitted GLM or GLMM model for overdispersion. Uses the `performance` package if available; otherwise falls back to a manual deviance/df ratio check.

```r
check_overdispersion(model)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `model` | glm/glmer | A fitted model object |

**Returns:** `list` with `ratio` (numeric), `p_value` (numeric), `is_overdispersed` (logical).

```r
od <- check_overdispersion(poisson_model)
if (od$is_overdispersed) message("Consider negative binomial model")
```

---

#### `fit_nb_glmm`

Fit a negative binomial GLMM for differential abundance testing. Falls back to Poisson GLMM if the NB fit fails to converge.

```r
fit_nb_glmm(count_data, formula = count ~ condition + (1 | sample_id), max_iter = 20000)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `count_data` | data.frame | Must contain columns referenced in `formula` |
| `formula` | formula | Model formula (default: `count ~ condition + (1 \| sample_id)`) |
| `max_iter` | integer | Maximum optimizer iterations |

**Returns:** Fitted model object (`glmerMod`), or `NULL` if both NB and Poisson fits fail.

```r
model <- fit_nb_glmm(count_df, formula = count ~ condition + (1 | sample_id))
if (!is.null(model)) summary(model)
```

---

#### `run_pca`

Run Principal Component Analysis on an expression matrix using `prcomp()`.

```r
run_pca(expr_matrix, n_components = 2, scale = TRUE)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_matrix` | numeric matrix | Rows = cells/samples, columns = markers/features |
| `n_components` | integer | Number of principal components to compute |
| `scale` | logical | Scale variables before PCA |

**Returns:** `list` with `coordinates` (matrix), `variance_explained` (numeric vector), `loadings` (matrix), `pca_object` (prcomp result).

```r
pca_res <- run_pca(expr_matrix, n_components = 3)
cat("Variance explained:", round(pca_res$variance_explained, 3), "\n")
```

---

#### `run_umap`

Run UMAP dimensionality reduction on an expression matrix using the `uwot` package.

```r
run_umap(expr_matrix, n_neighbors = 15, min_dist = 0.1,
         n_components = 2, metric = "euclidean", seed = 42)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_matrix` | numeric matrix | Rows = cells, columns = markers |
| `n_neighbors` | integer | Number of nearest neighbors |
| `min_dist` | numeric | Minimum distance parameter |
| `n_components` | integer | Number of output dimensions |
| `metric` | character | Distance metric (e.g. `"euclidean"`, `"cosine"`) |
| `seed` | integer | Random seed for reproducibility |

**Returns:** `numeric matrix` -- UMAP coordinates (n_cells x n_components).

```r
umap_coords <- run_umap(expr_matrix, n_neighbors = 30, min_dist = 0.2, seed = 123)
plot(umap_coords[, 1], umap_coords[, 2], pch = 16, cex = 0.3)
```

---

#### `calculate_statistics`

Calculate per-cluster, per-marker summary statistics (mean, median, SD, quartiles, count).

```r
calculate_statistics(expr_data, cluster_col, marker_cols)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_data` | data.frame | Must contain `cluster_col` and all `marker_cols` |
| `cluster_col` | character | Name of the cluster/group column |
| `marker_cols` | character vector | Names of marker columns to summarize |

**Returns:** `tibble` with columns: cluster, marker, mean, median, sd, q25, q75, n.

```r
stats <- calculate_statistics(expr_data, "cluster", c("CD3", "CD4", "CD8"))
```

---

#### `get_top_features`

Identify the top N features (columns) by variance or mean expression.

```r
get_top_features(expr_matrix, n_top = 20, method = "variance")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_matrix` | numeric matrix | Expression matrix |
| `n_top` | integer | Number of top features to return |
| `method` | character | Ranking criterion: `"variance"` or `"mean"` |

**Returns:** `character vector` -- names of the top features, in descending order.

```r
top_markers <- get_top_features(expr_matrix, n_top = 15, method = "variance")
```

---

### clustering_dispatch.R

Method dispatcher that routes clustering to Leiden, FlowSOM, PhenoGraph, or DBSCAN via a single function call. All methods are normalized to return an integer vector of cluster IDs with length equal to the number of input cells.

**Required packages:** Depends on chosen method (loaded at runtime from `src/clustering/`)

**Used by:** pipelines (A, E), analysis

#### `run_clustering`

Dispatch clustering to the specified method. Sources the corresponding module from `src/clustering/` at runtime. DBSCAN metadata (noise indices, eps used) is attached as the `"dbscan_metadata"` attribute on the returned vector. Leiden output is coerced from numeric to integer.

```r
run_clustering(expr_matrix, method = "leiden", params = list())
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_matrix` | numeric matrix | Expression matrix (cells x markers) |
| `method` | character | One of `"leiden"`, `"flowsom"`, `"phenograph"`, `"dbscan"` |
| `params` | named list | Method-specific parameters (see table below) |

**Method-specific parameters in `params`:**

| Key | Methods | Default | Description |
|-----|---------|---------|-------------|
| `n_neighbors` | leiden, phenograph | 30 | Number of nearest neighbors for graph construction |
| `resolution` | leiden | 1.0 | Leiden resolution parameter |
| `n_clusters` | flowsom | 20 | Number of meta-clusters for FlowSOM |
| `grid_xdim` | flowsom | 10 | SOM grid X dimension |
| `grid_ydim` | flowsom | 10 | SOM grid Y dimension |
| `k` | phenograph | 30 | K nearest neighbors (alias for `n_neighbors`) |
| `eps` / `dbscan_eps` | dbscan | auto | DBSCAN epsilon; `NULL` triggers automatic estimation |
| `min_pts` / `dbscan_min_pts` | dbscan | 5 | DBSCAN minimum points |
| `seed` | all | 42 | Random seed |

**Returns:** `integer vector` of cluster IDs (length == `nrow(expr_matrix)`). For DBSCAN, noise points are labeled `0`; metadata is in `attr(result, "dbscan_metadata")`.

```r
source(file.path(shared_dir, "clustering_dispatch.R"))
ids <- run_clustering(expr_matrix, method = "flowsom",
                      params = list(n_clusters = 20, seed = 42))
table(ids)
```

---

### dimred_dispatch.R

Method dispatcher that routes dimensionality reduction to UMAP, t-SNE, or diffusion map via a single function call. All methods return a standardized list with coordinates, method name, optional plots, and extra metadata.

**Required packages:** Depends on chosen method (loaded at runtime from `src/visualization/`)

**Used by:** pipelines (A, E), visualization

#### `run_dimred`

Dispatch dimensionality reduction to the specified method. Sources the corresponding visualization module at runtime. If `output_dir` is provided, plots are generated; otherwise only coordinates are computed.

```r
run_dimred(expr_matrix, method = "umap",
           sample_labels = NULL, cluster_labels = NULL,
           condition_labels = NULL, output_dir = NULL,
           params = list())
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_matrix` | numeric matrix | Expression matrix (cells x markers) |
| `method` | character | One of `"umap"`, `"tsne"`, `"diffusion_map"` |
| `sample_labels` | character vector or `NULL` | Sample label per cell (length == nrow) |
| `cluster_labels` | integer/character vector or `NULL` | Cluster label per cell (length == nrow) |
| `condition_labels` | character vector or `NULL` | Condition label per cell (length == nrow) |
| `output_dir` | character or `NULL` | Directory for plot output; `NULL` = no plots |
| `params` | named list | Method-specific parameters (see table below) |

**Method-specific parameters in `params`:**

| Key | Methods | Default | Description |
|-----|---------|---------|-------------|
| `n_neighbors` / `umap_n_neighbors` | umap | 15 | Number of UMAP neighbors |
| `umap_min_dist` / `min_dist` | umap | 0.1 | UMAP minimum distance |
| `umap_coords` | umap | `NULL` | Pre-computed UMAP coordinates (skip recomputation) |
| `tsne_perplexity` / `perplexity` | tsne | 30 | t-SNE perplexity |
| `tsne_max_iter` / `max_iter` | tsne | 1000 | t-SNE maximum iterations |
| `tsne_coords` | tsne | `NULL` | Pre-computed t-SNE coordinates |
| `dm_n_components` / `n_components` | diffusion_map | 3 | Number of diffusion components |
| `dm_sigma` / `sigma` | diffusion_map | `"local"` | Diffusion map kernel width |
| `output_formats` | all | `list(pdf=FALSE, png=TRUE, svg=FALSE)` | Plot format toggles |
| `plot_width` | all | 10 | Plot width in inches |
| `plot_height` | all | 8 | Plot height in inches |
| `seed` | all | 42 | Random seed |

**Returns:** `list` with:
- `coords` -- numeric matrix (n_cells x n_dimensions)
- `method` -- character string (`"umap"`, `"tsne"`, or `"diffusion_map"`)
- `plots` -- list of ggplot objects (or `NULL` if `output_dir` is `NULL`)
- `extra` -- list with additional metadata (`eigenvalues` for diffusion_map; `NULL` otherwise)

```r
source(file.path(shared_dir, "dimred_dispatch.R"))
result <- run_dimred(expr_matrix, method = "tsne",
                     sample_labels = labels, output_dir = "results/",
                     params = list(perplexity = 30))
dim(result$coords)  # [n_cells, 2]
```

---

### batch_correction_dispatch.R

Method dispatcher that routes batch correction to Harmony, CytoNorm, ComBat, or a pass-through. Accepts either a list of flowFrame objects or a numeric matrix, converting between formats as needed for each backend method.

**Required packages:** Depends on chosen method (loaded at runtime from `src/harmonization/`)

**Used by:** pipelines (A, E), harmonization

#### `run_batch_correction_dispatch`

Dispatch batch correction to the specified method. Bridges input types: if the backend expects flowFrames but receives a matrix (or vice versa), conversion is handled automatically.

```r
run_batch_correction_dispatch(data, method = "harmony",
                               batch_labels,
                               markers_of_interest = NULL,
                               sample_names = NULL,
                               params = list())
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `data` | list of flowFrames OR numeric matrix | Input data; if matrix, `sample_names` must provide per-row sample labels |
| `method` | character | One of `"harmony"`, `"cytonorm"`, `"combat"`, `"none"` |
| `batch_labels` | character vector | Batch assignments -- one per flowFrame (if list input) or one per row (if matrix input) |
| `markers_of_interest` | character vector or `NULL` | Marker names to correct (required for CytoNorm) |
| `sample_names` | character vector or `NULL` | Sample identifiers -- one per flowFrame or one per row |
| `params` | named list | Method-specific parameters (see table below) |

**Method-specific parameters in `params`:**

| Key | Methods | Default | Description |
|-----|---------|---------|-------------|
| `batch_do_pca` | harmony | `FALSE` | Run PCA before Harmony |
| `combat_method` | combat | `"parametric"` | ComBat estimation method |
| `combat_covariates` | combat | `NULL` | Covariates to preserve during ComBat |
| `cytonorm_train_indices` | cytonorm | all frames | Integer vector: which frames to use for training |
| `cytonorm_n_clusters` | cytonorm | 10 | Number of FlowSOM clusters for CytoNorm |
| `cytonorm_goal` | cytonorm | `"mean"` | CytoNorm normalization goal |
| `seed` | cytonorm | 42 | Random seed |

**Returns:** `list` with:
- `corrected_matrix` -- numeric matrix (cells x markers)
- `sample_labels` -- character vector (length == nrow)
- `batch_labels` -- character vector (length == nrow)
- `markers` -- character vector of marker names
- `corrected_frames` -- list of corrected flowFrames (CytoNorm only; `NULL` otherwise)

```r
source(file.path(shared_dir, "batch_correction_dispatch.R"))
result <- run_batch_correction_dispatch(
  ff_list, method = "harmony",
  batch_labels = c("batch1", "batch1", "batch2", "batch2"),
  markers_of_interest = c("CD3", "CD4", "CD8"),
  sample_names = c("S1", "S2", "S3", "S4")
)
dim(result$corrected_matrix)
```

---

### pipeline_utils.R

Pipeline boilerplate utilities for config loading, module sourcing, and parameter list construction. Eliminates repeated setup code from pipeline scripts.

**Required packages:** none (sources io_utils.R and memory_utils.R internally)

**Used by:** pipelines (A, B, C, D, E)

#### `init_pipeline`

Initialize a pipeline: parse `--config` from command-line arguments, source the config file, source core shared modules (io_utils.R and memory_utils.R), and record the start time.

```r
init_pipeline(pipeline_name, default_config)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `pipeline_name` | character | Short name for logging (e.g. `"Pipeline A"`) |
| `default_config` | character | Relative path from project root to the default config file (e.g. `"config/pipeline_A_config.R"`) |

**Returns:** Named `list` with `base_dir` (character), `shared_dir` (character), `start_time` (POSIXct).

```r
source(file.path(shared_dir, "pipeline_utils.R"))
env <- init_pipeline("Pipeline A", "config/pipeline_A_config.R")
# Prints: === Pipeline A ===
# Sources config and core modules
```

---

#### `finish_pipeline`

Print a pipeline completion message with elapsed time in minutes.

```r
finish_pipeline(pipeline_name, start_time)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `pipeline_name` | character | Short name (e.g. `"Pipeline A"`) |
| `start_time` | POSIXct | The `start_time` value from `init_pipeline()` |

**Returns:** `NULL` (invisible). Side effect: prints completion message.

```r
finish_pipeline("Pipeline A", env$start_time)
# === Pipeline A complete (3.2 minutes) ===
```

---

#### `clustering_params`

Build a params list for `run_clustering()` from a user config. Extracts and applies defaults for all clustering method parameters.

```r
clustering_params(params)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `params` | named list | The user's config `params` list |

**Returns:** Named `list` suitable for passing to `run_clustering(..., params = result)`.

**Extracted keys and defaults:**

| Key | Default |
|-----|---------|
| `n_neighbors` | 15 |
| `resolution` | 0.5 |
| `n_clusters` | 20 |
| `grid_xdim` | 10 |
| `grid_ydim` | 10 |
| `k` | `params$n_neighbors` or 30 |
| `dbscan_eps` | `NULL` |
| `dbscan_min_pts` | 5 |
| `seed` | 42 |

```r
cluster_p <- clustering_params(params)
ids <- run_clustering(expr_matrix, method = params$clustering_method, params = cluster_p)
```

---

#### `dimred_params`

Build a params list for `run_dimred()` from a user config. Extracts and applies defaults for all dimensionality reduction parameters.

```r
dimred_params(params)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `params` | named list | The user's config `params` list |

**Returns:** Named `list` suitable for passing to `run_dimred(..., params = result)`.

**Extracted keys and defaults:**

| Key | Default |
|-----|---------|
| `n_neighbors` | 15 |
| `umap_min_dist` | 0.1 |
| `tsne_perplexity` | 30 |
| `tsne_max_iter` | 1000 |
| `dm_n_components` | 3 |
| `dm_sigma` | `"local"` |
| `output_formats` | `list(pdf=FALSE, png=TRUE, svg=FALSE)` |
| `plot_width` | 10 |
| `plot_height` | 8 |
| `seed` | 42 |

```r
dr_p <- dimred_params(params)
result <- run_dimred(expr_matrix, method = params$dimred_method, params = dr_p)
```

---

#### `batch_params`

Build a params list for `run_batch_correction_dispatch()` from a user config.

```r
batch_params(params)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `params` | named list | The user's config `params` list |

**Returns:** Named `list` suitable for passing to `run_batch_correction_dispatch(..., params = result)`.

**Extracted keys and defaults:**

| Key | Default |
|-----|---------|
| `batch_do_pca` | `FALSE` |
| `combat_method` | `"parametric"` |
| `combat_covariates` | `NULL` |
| `cytonorm_train_indices` | `NULL` |
| `cytonorm_n_clusters` | 10 |
| `cytonorm_goal` | `"mean"` |
| `seed` | 42 |

```r
bp <- batch_params(params)
result <- run_batch_correction_dispatch(data, method = params$batch_method, params = bp)
```

---

## Quick Reference: Module Selection Guide

| Task | Module | Key Function |
|------|--------|-------------|
| Read/write CSV files | `io_utils.R` | `detect_csv_sep`, `save_and_cleanup` |
| Classify Cytolution CSV columns | `io_utils.R` | `classify_csv_columns` |
| Create output directories | `io_utils.R` | `make_outdir`, `dir_create` |
| Fatal error and exit | `io_utils.R` | `die` |
| Interactive marker selection | `gui_utils.R` | `tk_multiselect`, `select_markers` |
| Mode-based marker dispatch | `gui_utils.R` | `select_markers_by_mode` |
| Monitor RAM usage | `memory_utils.R` | `check_memory`, `get_available_memory_gb` |
| Load large files safely | `memory_utils.R` | `read_csv_chunked` |
| Parse cluster names | `cluster_utils.R` | `clean_cluster_name`, `extract_cluster_group` |
| Convert binary gating | `cluster_utils.R` | `labels_from_binary_df`, `convert_binary_to_single_cluster` |
| Translate labels | `cluster_utils.R` | `translate_labels` |
| Save plots in multiple formats | `plot_utils.R` | `save_plot_formats`, `save_base_plot_formats` |
| Create grouped barplots | `plot_utils.R` | `population_barplot`, `plot_per_sample` |
| Sanitize filenames | `plot_utils.R` | `create_safe_filename` |
| Inspect FCS marker panels | `flowframe_utils.R` | `explore_data`, `compare_panels` |
| Build combined expression matrix | `flowframe_utils.R` | `build_combined_expr_matrix` |
| Convert matrices to flowFrames | `flowframe_utils.R` | `from_mtx_to_ffr`, `from_ffl_to_fSet` |
| Validate FCS headers | `fcs_utils.R` | `check_tot`, `get_fcs_event_count` |
| Map CSV to FCS files | `fcs_utils.R` | `create_csv_fcs_mapping` |
| Batch load FCS files | `fcs_utils.R` | `load_fcs_list`, `load_fcs_folder` |
| Cluster enrichment testing | `stats_utils.R` | `calculate_enrichment`, `do_mixed_effects` |
| Differential abundance (NB GLMM) | `stats_utils.R` | `fit_nb_glmm`, `check_overdispersion` |
| Dimensionality reduction | `stats_utils.R` | `run_pca`, `run_umap` |
| Marker summary statistics | `stats_utils.R` | `calculate_statistics`, `get_top_features` |
| Swappable clustering (dispatch) | `clustering_dispatch.R` | `run_clustering` |
| Swappable dim. reduction (dispatch) | `dimred_dispatch.R` | `run_dimred` |
| Swappable batch correction (dispatch) | `batch_correction_dispatch.R` | `run_batch_correction_dispatch` |
| Pipeline initialization | `pipeline_utils.R` | `init_pipeline`, `finish_pipeline` |
| Build dispatcher param lists | `pipeline_utils.R` | `clustering_params`, `dimred_params`, `batch_params` |
