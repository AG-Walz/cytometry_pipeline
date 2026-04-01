[< Back to Documentation Index](../README.md)

# src/metadata_qc/ -- Metadata & Quality Assurance

Four scripts for inspecting study metadata (clinical variables, acquisition dates), comparing gating strategies across FlowJo workspaces, and computing inter-rater agreement statistics. These scripts are used primarily by Pipeline D and can also be used standalone for quality checks.

## Summary Table

| Script | Functions | Description |
|--------|-----------|-------------|
| `check_metadata.R` | `check_metadata()`, `plot_temporal_heatmap()` | Metadata barplots, batch-split plots, temporal heatmaps |
| `check_acquisition.R` | `get_fcs_keyword()`, `explore_fcs_keywords()`, `get_fcs_metadata()`, `check_acquisition_dates()` | FCS acquisition dates and duplicate detection |
| `compare_gatings.R` | `compare_gatings()` | Compare FlowJo gating strategies via Venn/Euler diagrams |
| `overlap_statistics.R` | `compute_fleiss_kappa()` | Fleiss' kappa inter-rater agreement |

---

## check_metadata.R

Creates compact metadata barplots for all clinical variables and temporal heatmaps. Useful for spotting data quality issues and batch effects before running analysis.

### `check_metadata()`

```r
check_metadata(metadata_file, output_dir, sep = ",",
               skip_columns = NULL, only_cols = NULL,
               sort_column = NULL, sort_desc = FALSE,
               n_splits = 3,
               placeholders = c("?", "na", "n/a", "null", "-"))
```

Reads a metadata CSV, classifies values (converting blanks and placeholders to "Unassigned"), and generates compact faceted barplots showing the distribution of each clinical variable. Optionally generates batch-split barplots to visualize potential confounding between acquisition order and metadata.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `metadata_file` | character | *required* | Path to metadata CSV |
| `output_dir` | character | *required* | Output directory for plots and CSV |
| `sep` | character | `","` | CSV separator |
| `skip_columns` | character vector | `NULL` | Columns to exclude from analysis |
| `only_cols` | character vector | `NULL` | Columns to include; `NULL` = all non-skipped |
| `sort_column` | character | `NULL` | Column name to sort rows by before splitting (e.g., `"Datum"` for dates). Auto-detects date formats. |
| `sort_desc` | logical | `FALSE` | Sort in descending order |
| `n_splits` | integer | `3` | Number of splits for batch-colored plot |
| `placeholders` | character vector | `c("?", "na", "n/a", "null", "-")` | Strings to treat as "Unassigned" |

**Returns:** Invisible `NULL`.

**Output Files:**

| File | Description |
|------|-------------|
| `consensus_barplots_compact.png` | Faceted barplot of all metadata variables with count labels |
| `consensus_barplots_compact_bySplit_{n}.png` | Batch-split stacked barplot (when `n_splits > 1`) |
| `metadata_summary_table.csv` | Long-format counts table with `column`, `value`, `Count` |

**Input from:** Clinical metadata CSV on disk.
**Output to:** Standalone QC reporting.

**Example:**

```r
source("check_metadata.R")
check_metadata(
  metadata_file = "/data/study_metadata.csv",
  output_dir    = "/output/metadata_qc",
  sep           = ";",
  skip_columns  = c("ID", "internal_code"),
  sort_column   = "Datum",
  n_splits      = 4
)
```

### `plot_temporal_heatmap()`

```r
plot_temporal_heatmap(metadata_file, output_dir, sep = ",",
                       skip_columns = NULL, date_column = "Datum",
                       placeholders = c("?", "na", "n/a", "null", "-"))
```

Creates a temporal heatmap showing the distribution of all metadata feature values over time (acquisition order). Each row is a `column:value` combination, each column is a timepoint. Color-coded by metadata column. Useful for detecting batch effects and confounding.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `metadata_file` | character | *required* | Path to metadata CSV |
| `output_dir` | character | *required* | Output directory |
| `sep` | character | `","` | CSV separator |
| `skip_columns` | character vector | `NULL` | Columns to exclude |
| `date_column` | character | `"Datum"` | Name of the date column for temporal ordering. Multiple date formats are auto-detected. |
| `placeholders` | character vector | `c("?", "na", "n/a", "null", "-")` | Strings to treat as "Unassigned" |

**Returns:** Invisible `NULL`.

**Output File:** `temporal_feature_heatmap.png`.

**Example:**

```r
source("check_metadata.R")
plot_temporal_heatmap(
  metadata_file = "/data/clinical_metadata.csv",
  output_dir    = "/output/metadata_qc",
  date_column   = "Datum"
)
```

**Dependencies:** `io_utils`, `plot_utils`. Requires `data.table`, `ggplot2`, `dplyr`, `tidyr`, `forcats`, `stringr`, `scales`.

---

## check_acquisition.R

Extracts FCS file metadata (acquisition dates, file sizes) for quality assurance and duplicate detection.

### `get_fcs_keyword()`

```r
get_fcs_keyword(fcs_path, keyword_name)
```

Extracts a single FCS keyword from a file. Reads only the header (1 event) for efficiency.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fcs_path` | character | *required* | Path to FCS file |
| `keyword_name` | character | *required* | FCS keyword to extract (e.g., `"$DATE"`, `"$CYT"`, `"$SRC"`) |

**Returns:** The keyword value as a character string, `"Unknown"` if not found, or `"Error"` on read failure.

**Input from:** FCS files on disk.
**Output to:** Standalone keyword lookup.

### `explore_fcs_keywords()`

```r
explore_fcs_keywords(fcs_path)
```

Reads all keywords from an FCS file and prints a summary including data type (`$DATATYPE`) and parameter bit depths (`$PnB`).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fcs_path` | character | *required* | Path to FCS file |

**Returns:** Named list of all keywords (invisible).

**Example:**

```r
source("check_acquisition.R")
keywords <- explore_fcs_keywords("/data/sample_001.fcs")
# File: sample_001.fcs
#   $DATATYPE: F
#   Parameter bit depths:
#     $P1B : 32
#     $P2B : 32
```

### `get_fcs_metadata()`

```r
get_fcs_metadata(fcs_path)
```

Extracts acquisition date (`$DATE` keyword) and file size from an FCS file.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fcs_path` | character | *required* | Path to FCS file |

**Returns:** Named character vector: `c(acquisition_date = "...", file_size = "...")`.

### `check_acquisition_dates()`

```r
check_acquisition_dates(parent_folder, output_dir)
```

Recursively scans for `.fcs` files, extracts acquisition dates and file sizes, and assigns unique IDs for files with matching date+size (potential duplicates).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `parent_folder` | character | *required* | Root folder to scan recursively |
| `output_dir` | character | *required* | Output directory for CSV results |

**Returns:** `data.frame` with columns:

| Column | Type | Description |
|--------|------|-------------|
| `file_relative_path` | character | Relative path from `parent_folder` |
| `acquisition_date` | character | `$DATE` keyword value |
| `file_size` | character | File size in bytes |
| `unique_file_id` | integer | ID grouping files with same date+size (likely duplicates) |

**Output File:** `acquisition_dates.csv`.

**Input from:** FCS files on disk.
**Output to:** Standalone QC reporting.

**Example:**

```r
source("check_acquisition.R")

# Quick keyword lookup
date <- get_fcs_keyword("/data/sample_001.fcs", "$DATE")
cat("Acquired:", date, "\n")
# Acquired: 15-MAR-2024

# Full acquisition check
dates <- check_acquisition_dates(
  parent_folder = "/data/all_fcs/",
  output_dir    = "/output/acquisition"
)
cat("Found", nrow(dates), "files\n")
cat("Unique file groups:", max(dates$unique_file_id), "\n")
```

**Dependencies:** `io_utils`. Requires `flowCore`.

---

## compare_gatings.R

Compares gating strategies across FlowJo workspaces by generating Venn/Euler diagrams for each population.

### `compare_gatings()`

```r
compare_gatings(pop_strings_list, target_fcs, label_dict,
                 custom_name = NULL, output_dir,
                 idx_labels = NULL, scale_size = TRUE)
```

Loads label lists from multiple FlowJo workspaces (pre-extracted via `load_multiple_workspaces()`), translates labels using a dictionary, and generates one Venn or Euler diagram per population showing overlap between analysts/methods.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pop_strings_list` | named list | *required* | Population string lists from `load_multiple_workspaces()`. Keys = workspace names, values = named lists of per-sample label vectors. |
| `target_fcs` | character | *required* | FCS filename to compare across workspaces |
| `label_dict` | named character vector | *required* | Maps gate paths to simplified labels (e.g., `c("/Live/CD3+" = "CD3+")`) |
| `custom_name` | character vector | `NULL` | Custom names for each workspace in diagrams |
| `output_dir` | character | *required* | Output directory for plots and CSVs |
| `idx_labels` | named list | `NULL` | Optional CSV-derived label vectors for additional comparison (e.g., from Pipeline C) |
| `scale_size` | logical | `TRUE` | `TRUE` = area-proportional Euler, `FALSE` = non-proportional Venn |

**Returns:** Invisible `data.frame` of combined label lists.

**Process:**

1. Extracts labels for `target_fcs` from each workspace
2. Translates labels via `label_dict` to a common vocabulary
3. Optionally merges CSV-derived labels for cross-method comparison
4. Generates one Venn/Euler diagram per unique population

**Output Files:**

| File | Description |
|------|-------------|
| `fj_label_lists_{target_fcs}.csv` | FlowJo-derived label lists per workspace |
| `label_lists_{target_fcs}.csv` | Combined label lists (FlowJo + optional CSV) |
| `subsetted {population}_{target_fcs}.png` | One Venn/Euler diagram per population |

**Input from:** `load_multiple_workspaces()` (loading/load_flowjo_wsp.R).
**Output to:** `compute_fleiss_kappa()` (overlap_statistics.R), standalone visualization.

**Example:**

```r
source("compare_gatings.R")

# Pre-load workspaces
pop_strings_list <- load_multiple_workspaces(
  wsp_paths = c("/data/analyst_A.wsp", "/data/analyst_B.wsp"),
  names     = c("AnalystA", "AnalystB")
)

# Define label dictionary
label_dict <- c(
  "/Live/Single Cells/CD3+"              = "CD3+",
  "/Live/Single Cells/CD3+/CD4+"         = "CD4+ T cells",
  "/Live/Single Cells/CD3+/CD8+"         = "CD8+ T cells",
  "/Live/Single Cells/CD19+"             = "B cells"
)

compare_gatings(
  pop_strings_list = pop_strings_list,
  target_fcs       = "sample_001.fcs",
  label_dict       = label_dict,
  custom_name      = c("Analyst A", "Analyst B"),
  output_dir       = "/output/gating_comparison",
  scale_size       = TRUE
)
```

**Dependencies:** `io_utils`, `cluster_utils`. Sources `visualization/venn_diagrams.R` at runtime.

---

## overlap_statistics.R

Computes Fleiss' kappa for inter-rater agreement assessment.

### `compute_fleiss_kappa()`

```r
compute_fleiss_kappa(label_lists, output_dir,
                      subset_specs = NULL,
                      subset_labels = NULL,
                      row_range = NULL)
```

Given a data frame where each column represents a rater's labels, computes Fleiss' kappa for the full dataset and optionally for user-specified subsets (e.g., pairwise rater comparisons). Uses `irr::kappam.fleiss()` with detail output.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `label_lists` | data.frame | *required* | Label assignments (columns = raters, rows = cells). At least 2 columns required. |
| `output_dir` | character | *required* | Output directory for log file |
| `subset_specs` | named list | `NULL` | List of column index vectors for subsets. If `NULL`, only the full dataset is analyzed. |
| `subset_labels` | character vector | `NULL` | Labels for each subset (same length as `subset_specs`) |
| `row_range` | integer vector | `NULL` | Row indices to subset before analysis |

**Returns:** Named list of `irr::kappam.fleiss` result objects (one per subset). Each result contains:

- `$value` -- the kappa statistic
- `$statistic` -- the z-statistic
- `$p.value` -- the p-value
- `$detail` -- per-category kappa values

**Output File:** `fleiss_kappa_log.txt` -- detailed log with kappa values, z-statistics, elapsed times, and per-category detail.

**Input from:** `compare_gatings()` (compare_gatings.R), `load_multiple_workspaces()` (loading/load_flowjo_wsp.R).
**Output to:** Standalone inter-rater agreement report.

**Kappa Interpretation:**

| Kappa | Agreement Level |
|-------|----------------|
| < 0.00 | Less than chance |
| 0.01 -- 0.20 | Slight |
| 0.21 -- 0.40 | Fair |
| 0.41 -- 0.60 | Moderate |
| 0.61 -- 0.80 | Substantial |
| 0.81 -- 1.00 | Almost perfect |

**Example:**

```r
source("overlap_statistics.R")

# label_lists: data.frame from compare_gatings() output
label_lists <- read.csv("/output/gating_comparison/label_lists_sample_001.csv")

# Full kappa (all raters)
results <- compute_fleiss_kappa(
  label_lists = label_lists,
  output_dir  = "/output/overlap_stats"
)
cat("Kappa:", results[["Full data"]]$value, "\n")

# Pairwise subsets
results <- compute_fleiss_kappa(
  label_lists  = label_lists,
  output_dir   = "/output/overlap_stats",
  subset_specs = list(c(1, 2), c(1, 3), c(2, 3)),
  subset_labels = c("A vs B", "A vs C", "B vs C")
)
for (nm in names(results)) {
  cat(nm, ": kappa =", round(results[[nm]]$value, 4), "\n")
}
```

**Dependencies:** `io_utils`. Requires `irr` package.
