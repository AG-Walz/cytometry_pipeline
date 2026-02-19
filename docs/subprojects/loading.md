[< Back to Documentation Index](../README.md)

# src/loading/ — Data Loading

Six scripts for ingesting flow cytometry data from FCS files, CSV expression matrices, FlowJo workspaces, and manual gating labels. Every script defines one or more functions that you call after sourcing. All scripts load shared utilities from `src/R_shared/`.

---

## load_fcs_raw.R

Loads raw FCS files from a directory into flowFrame objects. This is the primary FCS loader used by Pipelines A and E.

### Function

```r
load_fcs_raw(input_dir, sample_names = NULL, truncation = FALSE,
             min_limit = NULL, extract_csv = FALSE, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory containing `.fcs` files |
| `sample_names` | character vector | `NULL` | Specific samples to load; `NULL` loads all |
| `truncation` | logical | `FALSE` | Truncate values exceeding `$PnR` keyword range |
| `min_limit` | numeric | `NULL` | Skip files with fewer than this many cells |
| `extract_csv` | logical | `FALSE` | Also export expression matrices as CSVs |
| `output_dir` | character | `NULL` | Directory for extracted CSVs |

### Returns

```r
list(
  flow_frames  = list(sample_001 = <flowFrame>, sample_002 = <flowFrame>, ...),
  sample_names = c("sample_001", "sample_002", "sample_003"),
  input_files  = c("/path/to/sample_001.fcs", "/path/to/sample_002.fcs", ...)
)
```

### Example

```r
setwd("/path/to/cytometry/src/loading")
source("load_fcs_raw.R")

result <- load_fcs_raw(
  input_dir  = "/data/fcs_files",
  truncation = FALSE
)

cat("Loaded", length(result$flow_frames), "samples\n")
cat("Markers:", colnames(flowCore::exprs(result$flow_frames[[1]])), "\n")
```

**Console output:**

```
  Loading: sample_001.fcs ... 50000 cells x 8 markers
  Loading: sample_002.fcs ... 48000 cells x 8 markers
  Loading: sample_003.fcs ... 52000 cells x 8 markers
Loaded 3 samples
```

### CLI Usage

```bash
Rscript load_fcs_raw.R <input_dir> [output_dir]
```

### Dependencies

`io_utils`, `fcs_utils`, `flowframe_utils`

---

## load_csv_expression.R

Loads a directory of CSV expression matrices into data.table objects. Used by Pipeline B.

### Function

```r
load_csv_expression(input_dir, sep = NULL, max_memory_gb = 4)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with CSV files |
| `sep` | character | `NULL` | Separator; `NULL` auto-detects |
| `max_memory_gb` | numeric | `4` | Memory budget for large files |

### Returns

Named list of `data.table` objects, one per CSV file.

### Example Input

```
csv_expression/
├── sample_001_CD4_T_cells.csv
├── sample_001_CD8_T_cells.csv
└── sample_002_CD4_T_cells.csv
```

Each CSV:

```
CD3;CD4;CD8;CD19;CD56;HLA-DR
3.94;4.31;5.33;3.23;2.35;2.78
2.75;4.80;3.52;3.42;2.00;4.06
```

### Example

```r
source("load_csv_expression.R")
expr_data <- load_csv_expression(input_dir = "/data/csv_expression")

names(expr_data)
# [1] "sample_001_CD4_T_cells" "sample_001_CD8_T_cells" "sample_002_CD4_T_cells"

head(expr_data[[1]])
#      CD3   CD4   CD8  CD19  CD56 HLA-DR
# 1: 3.94  4.31  5.33  3.23  2.35   2.78
# 2: 2.75  4.80  3.52  3.42  2.00   4.06
```

### CLI Usage

```bash
Rscript load_csv_expression.R <input_dir>
```

### Dependencies

`io_utils`, `memory_utils`

---

## load_csv_gating.R

Loads CSV files containing gating information alongside FCS files. Automatically classifies columns, validates CSV vs FCS row counts, and splits cells by population. Used by Pipeline C.

### Function

```r
load_csv_gating(input_fcs_path, input_csv_path = input_fcs_path,
                output_path, gating_pattern,
                csv_label_dict = NULL, gating_tree = NULL,
                processing_mode = "automatic", output_format = "wide",
                chunk_size = 100000)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_fcs_path` | character | *required* | Directory with FCS files |
| `input_csv_path` | character | `input_fcs_path` | Directory with CSV files |
| `output_path` | character | auto-generated | Output directory for split files |
| `gating_pattern` | character | see below | Regex to identify gating columns |
| `processing_mode` | character | `"automatic"` | `"automatic"` (GUI), `"manual"`, `"expression_only"` |
| `output_format` | character | `"wide"` | `"wide"` (per-population CSVs) or `"long"` (single CSV) |

Default gating pattern:

```r
"\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$"
```

### Returns

```r
list(
  output_dir  = "/path/to/output",
  results     = list(...),     # Per-file processing results
  n_processed = 3,
  n_total     = 3
)
```

### Example Input

A directory with paired CSV + FCS files:

```
csv_gating/
├── sample_001.csv    # Expression + gating columns
├── sample_001.fcs    # Matching FCS for validation
├── sample_002.csv
└── sample_002.fcs
```

Each CSV:

```
CD3;CD4;CD8;CD19;CD56;HLA-DR;CD4+ (Clustering annotated);CD8+ (Clustering annotated)
2.46;4.01;2.15;1.29;2.37;1.22;1;0
3.89;2.12;4.56;1.87;3.21;3.45;0;1
```

### Example Output

The function splits each CSV by population:

```
output/01_loaded/
├── sample_001_CD4_.csv     # Expression data for CD4+ cells only
├── sample_001_CD8_.csv     # Expression data for CD8+ cells only
├── sample_002_CD4_.csv
└── sample_002_CD8_.csv
```

### Example

```r
source("load_csv_gating.R")
result <- load_csv_gating(
  input_fcs_path  = "/data/csv_gating",
  input_csv_path  = "/data/csv_gating",
  output_path     = "/output/loaded",
  gating_pattern  = "\\(Clustering annotated\\)$",
  processing_mode = "manual"
)
cat("Processed", result$n_processed, "of", result$n_total, "files\n")
```

### CLI Usage

```bash
Rscript load_csv_gating.R <input_fcs_path> [output_path]
```

### Dependencies

`io_utils`, `gui_utils`, `fcs_utils`, `cluster_utils`, `memory_utils`

---

## load_flowjo_wsp.R

Loads FlowJo workspace (`.wsp`) files and extracts gating set information. Used by Pipeline D.

### Function

```r
load_flowjo_wsp(wsp_path, fcs_dir = dirname(wsp_path), name = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wsp_path` | character | *required* | Path to `.wsp` file |
| `fcs_dir` | character | `dirname(wsp_path)` | Directory with referenced FCS files |
| `name` | character | `NULL` | Label for this workspace |

### Returns

```r
list(
  gating_set       = <GatingSet>,
  sample_names     = c("sample_001", "sample_002"),
  population_names = c("/Live", "/Live/Single Cells", "/Live/Single Cells/CD3+", ...),
  workspace_name   = "my_workspace"
)
```

### Helper Functions

```r
get_population_stats(gs, type = NULL)           # Get gate statistics
get_population_strings_list(gs, sample_names)   # Extract per-event labels
summarize_population_counts(pop_strings_list)   # Count cells per population
save_population_strings(pop_strings_list, out_dir)  # Save to disk
```

### Example

```r
source("load_flowjo_wsp.R")
ws <- load_flowjo_wsp("/data/analyst_A.wsp", fcs_dir = "/data/fcs/")
cat("Populations:", paste(ws$population_names, collapse = ", "), "\n")
```

### Dependencies

`io_utils`, `flowframe_utils`. Requires `flowWorkspace` and `CytoML` packages.

---

## load_manual_gating.R

Loads manual gating labels from CSV files, optionally matched to FCS data.

### Function

```r
load_manual_gating(csv_dir, fcs_dir = csv_dir,
                   gating_pattern = "Manual gating",
                   label_dict = NULL, gating_tree = NULL,
                   output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_dir` | character | *required* | Directory with CSV files |
| `fcs_dir` | character | `csv_dir` | Directory with matching FCS files |
| `gating_pattern` | character | `"Manual gating"` | Column pattern for gating labels |
| `label_dict` | named list | `NULL` | Label translation dictionary |
| `output_dir` | character | `NULL` | If set, writes per-sample label CSVs |

### Returns

Named list of label vectors per sample. When `output_dir` is set, also writes CSV files with `cell_id` and `label` columns.

### Dependencies

`io_utils`, `fcs_utils`, `cluster_utils`

---

## create_flowframes.R

Creates flowFrame objects from expression matrices, and optionally merges them into a flowSet.

### Functions

```r
create_flowframes_from_matrices(expr_list, markers_of_interest = NULL,
                                 save_fcs = FALSE, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_list` | named list | *required* | List of expression matrices or data.tables |
| `markers_of_interest` | character vector | `NULL` | Markers to include; `NULL` keeps all |
| `save_fcs` | logical | `FALSE` | Write FCS files to disk |
| `output_dir` | character | `NULL` | Directory for saved FCS files |

```r
create_merged_flowset(flow_frames)   # Merge flowFrames into a flowSet
```

### Example

```r
source("create_flowframes.R")

# Convert matrices to flowFrames
expr_list <- list(
  sample_001 = matrix(rnorm(300), ncol = 6, dimnames = list(NULL, c("CD3","CD4","CD8","CD19","CD56","HLA-DR"))),
  sample_002 = matrix(rnorm(300), ncol = 6, dimnames = list(NULL, c("CD3","CD4","CD8","CD19","CD56","HLA-DR")))
)

ff_list <- create_flowframes_from_matrices(expr_list, save_fcs = TRUE, output_dir = "/output/fcs")
# Creates sample_001.fcs and sample_002.fcs
```

### Dependencies

`io_utils`, `flowframe_utils`
