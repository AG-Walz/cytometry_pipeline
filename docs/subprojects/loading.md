[< Back to Documentation Index](../README.md)

# src/loading/ -- Data Loading

6 scripts for ingesting flow cytometry data from FCS files, CSV expression matrices, FlowJo workspaces, and manual gating labels. All scripts source shared utilities via:

```r
shared_dir <- file.path(dirname(getwd()), "R_shared")
```

---

## load_fcs_raw.R

Loads raw FCS files from a directory into flowFrame objects.

### Main Function

```r
load_fcs_raw(input_dir, sample_names = NULL, truncation = FALSE,
             min_limit = NULL, extract_csv = FALSE, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory containing `.fcs` files |
| `sample_names` | character vector | `NULL` | Subset of sample names to load; `NULL` loads all |
| `truncation` | logical | `FALSE` | Whether to truncate extreme values |
| `min_limit` | numeric | `NULL` | Minimum cell count threshold; files below this are skipped |
| `extract_csv` | logical | `FALSE` | If `TRUE`, also exports expression matrices as CSVs |
| `output_dir` | character | `NULL` | Directory for extracted CSVs (required when `extract_csv = TRUE`) |

### Output

```r
list(
  flow_frames  = list(),    # Named list of flowFrame objects
  sample_names = character(), # Vector of sample identifiers
  input_files  = character()  # Full paths to loaded .fcs files
)
```

### CLI Usage

```bash
Rscript load_fcs_raw.R <input_dir> [output_dir]
```

### src/R_shared Dependencies

`io_utils`, `fcs_utils`, `flowframe_utils`

---

## load_csv_gating.R

Loads CSV files containing gating information alongside FCS files. Auto-classifies columns, supports interactive marker selection, and validates CSV vs FCS row counts.

### Main Function

```r
load_csv_gating(input_fcs_path, input_csv_path = input_fcs_path,
                output_path = ..., gating_pattern = "\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$",
                csv_label_dict = NULL, gating_tree = NULL,
                processing_mode = "automatic", output_format = "wide",
                chunk_size = 100000)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_fcs_path` | character | *required* | Path to directory with FCS files |
| `input_csv_path` | character | `input_fcs_path` | Path to directory with CSV files |
| `output_path` | character | auto-generated | Output directory for processed files |
| `gating_pattern` | character | see signature | Regex to identify gating columns |
| `csv_label_dict` | named list | `NULL` | Dictionary for label translation |
| `gating_tree` | list | `NULL` | Hierarchical gating structure |
| `processing_mode` | character | `"automatic"` | `"automatic"` or `"interactive"` |
| `output_format` | character | `"wide"` | `"wide"` (per-population CSVs) or `"long"` (single CSV with population column) |
| `chunk_size` | integer | `100000` | Rows per chunk for memory-efficient processing |

### Helper Functions

```r
# Process a single CSV file with gating data
process_single_csv_gating(csv_file, gating_pattern, expr_cols = NULL, ...)

# Sanitize a name for use as a filename
create_safe_filename(name)
```

### Output

```r
list(
  output_dir  = character(), # Path to output directory
  results     = list(),      # Per-file processing results
  n_processed = integer(),   # Number of files successfully processed
  n_total     = integer()    # Total number of files attempted
)
```

- **Wide format**: One CSV per population in the output directory
- **Long format**: Single CSV with an added `population` column

### Features

- Auto-classifies columns via `classify_csv_columns()`
- Interactive marker selection when `processing_mode = "interactive"`
- Validates that CSV and FCS row counts match

### CLI Usage

```bash
Rscript load_csv_gating.R <input_fcs_path> [output_path]
```

### src/R_shared Dependencies

`io_utils`, `gui_utils`, `fcs_utils`, `cluster_utils`, `memory_utils`

---

## load_csv_expression.R

Loads a directory of CSV expression matrices into data.table objects.

### Main Function

```r
load_csv_expression(input_dir, sep = NULL, max_memory_gb = 4)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory containing expression CSV files |
| `sep` | character | `NULL` | Column separator; `NULL` for auto-detection |
| `max_memory_gb` | numeric | `4` | Maximum memory budget in GB |

### Output

Named list of `data.table` objects, one per CSV file.

### CLI Usage

```bash
Rscript load_csv_expression.R <input_dir>
```

### src/R_shared Dependencies

`io_utils`, `memory_utils`

---

## load_flowjo_wsp.R

Loads FlowJo workspace (`.wsp`) files and extracts gating set information.

### Main Function

```r
load_flowjo_wsp(wsp_path, fcs_dir = dirname(wsp_path), name = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wsp_path` | character | *required* | Path to the `.wsp` file |
| `fcs_dir` | character | `dirname(wsp_path)` | Directory containing the referenced FCS files |
| `name` | character | `NULL` | Optional label for this workspace |

### Output

```r
list(
  gating_set       = GatingSet,    # flowWorkspace GatingSet object
  sample_names     = character(),   # Vector of sample identifiers
  population_names = character(),   # All population/gate names
  workspace_name   = character()    # Label for this workspace
)
```

### Helper Functions

```r
# Get population statistics from a GatingSet
get_population_stats(gs, type = NULL)

# Extract population label strings per sample
get_population_strings_list(gs, sample_names)

# Summarize population counts across samples
summarize_population_counts(pop_strings_list)

# Save population strings to disk
save_population_strings(pop_strings_list, out_dir)
```

### Notes

- No CLI interface. Used as a library function.
- Requires `flowWorkspace` and `CytoML` packages.

### src/R_shared Dependencies

`io_utils`, `flowframe_utils`

---

## load_manual_gating.R

Loads manual gating labels from CSV files, optionally matched to FCS data.

### Main Function

```r
load_manual_gating(csv_dir, fcs_dir = csv_dir,
                   gating_pattern = "Manual gating",
                   label_dict = NULL, gating_tree = NULL,
                   output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_dir` | character | *required* | Directory containing CSV files with gating labels |
| `fcs_dir` | character | `csv_dir` | Directory containing matching FCS files |
| `gating_pattern` | character | `"Manual gating"` | Pattern to identify gating columns |
| `label_dict` | named list | `NULL` | Dictionary for translating label names |
| `gating_tree` | list | `NULL` | Hierarchical gating structure |
| `output_dir` | character | `NULL` | If set, saves CSVs with `cell_id` + `label` columns |

### Output

Named list of label vectors per sample. When `output_dir` is specified, also writes CSV files with `cell_id` and `label` columns.

### Notes

- No CLI interface. Used as a library function.

### src/R_shared Dependencies

`io_utils`, `fcs_utils`, `cluster_utils`

---

## create_flowframes.R

Creates flowFrame objects from expression matrices, and optionally merges them into a flowSet.

### Main Function

```r
create_flowframes_from_matrices(expr_list, markers_of_interest = NULL,
                                 save_fcs = FALSE, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_list` | named list | *required* | Named list of expression matrices or data.tables |
| `markers_of_interest` | character vector | `NULL` | Subset of markers to include; `NULL` keeps all |
| `save_fcs` | logical | `FALSE` | Whether to write FCS files to disk |
| `output_dir` | character | `NULL` | Directory for saved FCS files |

### Helper Function

```r
# Merge multiple flowFrames into a single flowSet
create_merged_flowset(flow_frames)  # Returns: flowSet
```

### Notes

- No CLI interface. Used as a library function.

### src/R_shared Dependencies

`io_utils`, `flowframe_utils`
