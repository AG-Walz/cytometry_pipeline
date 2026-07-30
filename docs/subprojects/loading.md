[< Back to Documentation Index](../README.md)

# src/loading/ -- Data Loading

Seven scripts for ingesting flow cytometry data from FCS files, CSV expression matrices, FlowJo workspaces, and manual gating labels. Every script defines one or more functions that you call after sourcing. All scripts load shared utilities from `src/R_shared/`.

## Summary Table

| Script | Functions | Description |
|--------|-----------|-------------|
| `load_fcs_raw.R` | `load_fcs_raw()` | Load raw FCS files into flowFrame objects |
| `load_csv_expression.R` | `load_csv_expression()` | Load CSV expression matrices into data.tables |
| `load_csv_gating.R` | `load_csv_gating()`, `process_single_csv_gating()`, `create_safe_filename()` | Load CSV+FCS gating data, split by population |
| `load_flowjo_wsp.R` | `load_flowjo_wsp()`, `get_population_stats()`, `get_population_strings_list()`, `summarize_population_counts()`, `save_population_strings()`, `load_multiple_workspaces()` | Load FlowJo workspaces and extract gating data |
| `flowjo_workspace.R` | `find_workspace_files()`, `workspace_labels()`, `index_fcs_dir()`, `resolve_fcs_names()`, `wsp_group_table()`, `wsp_sample_files()`, `resolve_group()`, `resolve_workspace_specs()`, `format_resolution_report()`, `build_workspace_gatingset()`, `disambiguate_sample_labels()` | Find many `.wsp` under a folder and pair each with its FCS files by name |
| `load_manual_gating.R` | `load_manual_gating()` | Load manual gating labels from CSVs |
| `create_flowframes.R` | `create_flowframes_from_matrices()`, `create_merged_flowset()` | Create flowFrame/flowSet objects from matrices |

---

## load_fcs_raw.R

Loads raw FCS files from a directory into flowFrame objects. This is the primary FCS loader used by Pipelines A and E.

### `load_fcs_raw()`

```r
load_fcs_raw(input_dir, sample_names = NULL, truncation = FALSE,
             min_limit = NULL, extract_csv = FALSE, output_dir = NULL)
```

Reads all `.fcs` files from a directory using `flowCore::read.FCS()`, returning named flowFrame objects. Optionally exports expression matrices as CSVs.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory containing `.fcs` files |
| `sample_names` | character vector | `NULL` | Specific sample names; `NULL` derives from filenames |
| `truncation` | logical | `FALSE` | Truncate values exceeding `$PnR` keyword range |
| `min_limit` | numeric | `NULL` | Lower truncation limit passed to `read.FCS()` |
| `extract_csv` | logical | `FALSE` | Also export expression matrices as CSVs |
| `output_dir` | character | `NULL` | Directory for extracted CSVs (auto-generated if `extract_csv = TRUE` and `NULL`) |

**Returns:**

```r
list(
  flow_frames  = list(sample_001 = <flowFrame>, sample_002 = <flowFrame>, ...),
  sample_names = c("sample_001", "sample_002", "sample_003"),
  input_files  = c("/path/to/sample_001.fcs", "/path/to/sample_002.fcs", ...)
)
```

**Input from:** `.fcs` files on disk.
**Output to:** `apply_compensation()` (preprocessing/compensation.R), `run_preprocessing()` (preprocessing/preprocess_pipeline.R), `downsample_flowframes()` (preprocessing/downsampling.R), `plot_umap()` (visualization/umap.R).

**Example:**

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

**CLI Usage:**

```bash
Rscript load_fcs_raw.R <input_dir> [output_dir]
```

**Dependencies:** `io_utils`, `fcs_utils`, `flowframe_utils`

---

## load_csv_expression.R

Loads a directory of CSV expression matrices into data.table objects. Used by Pipeline B.

### `load_csv_expression()`

```r
load_csv_expression(input_dir, sep = NULL, max_memory_gb = 4)
```

Reads all CSV files from a directory using `data.table::fread()`. Auto-detects the separator when `sep = NULL`. Uses chunked reading for files exceeding the memory budget.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with CSV files |
| `sep` | character | `NULL` | Separator; `NULL` auto-detects via `detect_csv_sep()` |
| `max_memory_gb` | numeric | `4` | Memory budget; files larger than this use chunked reading |

**Returns:** Named list of `data.table` objects, one per CSV file (keyed by filename without extension).

**Input from:** CSV expression files on disk.
**Output to:** `create_flowframes_from_matrices()` (loading/create_flowframes.R), `combine_expression_matrices()` (analysis/expression_matrices.R).

**Example:**

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

**CLI Usage:**

```bash
Rscript load_csv_expression.R <input_dir>
```

**Dependencies:** `io_utils`, `memory_utils`

---

## load_csv_gating.R

Loads CSV files containing gating information alongside FCS files. Automatically classifies columns, validates CSV vs FCS row counts, and splits cells by population. Used by Pipeline C.

### `load_csv_gating()`

```r
load_csv_gating(input_fcs_path, input_csv_path = input_fcs_path,
                output_path, gating_pattern,
                csv_label_dict = NULL, gating_tree = NULL,
                processing_mode = "automatic", output_format = "wide",
                chunk_size = 100000)
```

Discovers paired CSV+FCS files, classifies columns into expression/gating/scatter categories, and splits expression data by population using binary gating columns.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_fcs_path` | character | *required* | Directory with FCS files |
| `input_csv_path` | character | `input_fcs_path` | Directory with CSV files |
| `output_path` | character | auto-generated | Output directory for split files |
| `gating_pattern` | character | `"\\(Manual gating\\)$\|\\(Clustering annotated\\)$\|^Clustering raw$"` | Regex to identify gating columns |
| `csv_label_dict` | named vector | `NULL` | Label translation dictionary |
| `gating_tree` | named numeric vector | `NULL` | Gating hierarchy depths for label flattening |
| `processing_mode` | character | `"automatic"` | `"automatic"` (GUI), `"manual"`, `"expression_only"` |
| `output_format` | character | `"wide"` | `"wide"` (per-population CSVs) or `"long"` (single CSV with population column) |
| `chunk_size` | integer | `100000` | Rows per processing chunk |

**Returns:**

```r
list(
  output_dir  = "/path/to/output",
  results     = list(...),     # Per-file processing results
  n_processed = 3,
  n_total     = 3
)
```

**Input from:** CSV + FCS files on disk.
**Output to:** `split_by_cluster()` (analysis/cluster_splitting.R), `run_differential_abundance()` (analysis/differential_abundance.R).

**Example:**

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

**CLI Usage:**

```bash
Rscript load_csv_gating.R <input_fcs_path> [output_path]
```

### `process_single_csv_gating()` (internal)

```r
process_single_csv_gating(csv_file, gating_pattern, expr_cols, csv_label_dict,
                           gating_tree, csv_to_fcs, output_path, output_format, chunk_size)
```

Processes a single CSV file: reads data, identifies gating/expression columns, validates against matched FCS row count, and writes per-population or long-format output files.

**Returns:** `list(file, n_cells, n_gating_cols)`

### `create_safe_filename()`

```r
create_safe_filename(name)
```

Sanitizes a string for use as a filename by replacing non-alphanumeric characters with underscores.

**Returns:** Character string safe for filenames.

**Dependencies:** `io_utils`, `gui_utils`, `fcs_utils`, `cluster_utils`, `memory_utils`

---

## load_flowjo_wsp.R

Loads FlowJo workspace (`.wsp`) files and extracts gating set information. Used by Pipeline D.

### `load_flowjo_wsp()`

```r
load_flowjo_wsp(wsp_path, fcs_dir = dirname(wsp_path), name = NULL)
```

Parses a FlowJo `.wsp` file using CytoML, imports the gating hierarchy as a GatingSet, and extracts sample and population names.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wsp_path` | character | *required* | Path to `.wsp` file |
| `fcs_dir` | character | `dirname(wsp_path)` | Directory with referenced FCS files |
| `name` | character | `NULL` | Label for this workspace (defaults to filename without extension) |

**Returns:**

```r
list(
  gating_set       = <GatingSet>,
  sample_names     = c("sample_001", "sample_002"),
  population_names = c("/Live", "/Live/Single Cells", "/Live/Single Cells/CD3+", ...),
  workspace_name   = "my_workspace"
)
```

**Input from:** `.wsp` + `.fcs` files on disk.
**Output to:** `get_population_stats()`, `get_population_strings_list()`, `load_multiple_workspaces()` (all in this file).

**Example:**

```r
source("load_flowjo_wsp.R")
ws <- load_flowjo_wsp("/data/analyst_A.wsp", fcs_dir = "/data/fcs/")
cat("Populations:", paste(ws$population_names, collapse = ", "), "\n")
```

### `get_population_stats()`

```r
get_population_stats(gating_set, type = NULL)
```

Retrieves per-sample gate statistics (counts or frequencies) from a GatingSet.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gating_set` | GatingSet | *required* | A GatingSet object |
| `type` | character | `NULL` | `NULL` for counts, `"freq"` for frequencies |

**Returns:** `data.frame` with columns `Population`, `Sample`, `Value`.

**Input from:** `load_flowjo_wsp()`.
**Output to:** Standalone reporting.

### `get_population_strings_list()`

```r
get_population_strings_list(gating_set, sample_names)
```

Extracts per-event population path assignments for each cell in each sample. Iterates through all populations and assigns the deepest matching gate path to each cell.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gating_set` | GatingSet | *required* | A GatingSet object |
| `sample_names` | character vector | *required* | Sample names to process |

**Returns:** Named list of character vectors (one per sample), where each element is the population path for a cell.

**Input from:** `load_flowjo_wsp()`.
**Output to:** `load_multiple_workspaces()` (this file), `compare_gatings()` (metadata_qc/compare_gatings.R).

### `summarize_population_counts()`

```r
summarize_population_counts(population_strings_list)
```

Counts cells per population per sample from a population strings list.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `population_strings_list` | named list | *required* | Output from `get_population_strings_list()` |

**Returns:** `data.frame` with columns `Sample`, `Population`, `Count`.

### `save_population_strings()`

```r
save_population_strings(population_strings_list, out_dir)
```

Saves a population strings list to disk as an RDS file.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `population_strings_list` | named list | *required* | Population strings to save |
| `out_dir` | character | *required* | Output directory |

**Returns:** Invisible `NULL`. Writes `population_strings_list.rds`.

### `load_multiple_workspaces()`

```r
load_multiple_workspaces(wsp_paths, names, fcs_dirs = NULL)
```

Batch-loads multiple FlowJo workspaces and extracts population strings from each. Wraps `load_flowjo_wsp()` and `get_population_strings_list()` in a loop, returning a named list ready for `compare_gatings()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wsp_paths` | character vector | *required* | Paths to `.wsp` files |
| `names` | character vector | *required* | Labels for each workspace (same length as `wsp_paths`) |
| `fcs_dirs` | character vector | `NULL` | FCS directories per workspace; `NULL` uses `dirname()` of each `.wsp`. If length 1, recycled for all. |

**Returns:** Named list of population string lists, one entry per workspace (keyed by `names`).

**Input from:** `.wsp` files.
**Output to:** `compare_gatings()` (metadata_qc/compare_gatings.R), `compute_fleiss_kappa()` (metadata_qc/overlap_statistics.R).

**Example:**

```r
source("load_flowjo_wsp.R")

wsp_paths <- c("/data/analyst_A.wsp", "/data/analyst_B.wsp")
names     <- c("AnalystA", "AnalystB")

pop_strings_list <- load_multiple_workspaces(wsp_paths, names = names)
# Loading workspace 1 of 2 : /data/analyst_A.wsp
# Loading workspace 2 of 2 : /data/analyst_B.wsp
# Loaded 2 workspaces: AnalystA, AnalystB
```

**Dependencies:** `io_utils`, `flowframe_utils`. Requires `flowWorkspace` and `CytoML` packages.

---

## flowjo_workspace.R

Turns two **folder roots** — one holding `.wsp` workspaces, one holding `.fcs` files, both searched recursively — into ready-to-load workspace specs. Used by [`umap_from_flowjo()`](../scripts/umap-from-flowjo.md#many-workspaces-at-once) to put many workspaces into one shared embedding.

**Why it exists.** A `.wsp` records a full `uri` per sample, but it is the path on the Windows machine that built the workspace, so it does not exist on the analysis host. What a workspace reliably gives is the **list of FCS file names** — so those are read out of the XML first, the FCS root is indexed recursively by file name, and the two are matched. No directory convention is imposed on the data, and one workspace's files may live in several directories.

Everything is resolved **before any GatingSet is built**, and problems across all workspaces are reported together: a GatingSet load is slow, so discovering one mislaid file per load is not a workflow.

| Function | Signature | Returns |
|----------|-----------|---------|
| `find_workspace_files()` | `(root)` | Character vector of `.wsp` paths. `root` may be a `.wsp` file or a directory searched recursively |
| `workspace_labels()` | `(wsp_paths)` | Short display label per workspace; colliding base names are qualified by parent directory, then by full path |
| `index_fcs_dir()` | `(root)` | Named list: lower-cased file name -> every full path carrying it. Built once, shared by all workspaces |
| `resolve_fcs_names()` | `(names, fcs_index)` | `list(file, missing, ambiguous)`. Matching is case-insensitive, tolerates an omitted `.fcs`, and reduces Windows paths to their base name. Pure — no CytoML, no file system |
| `wsp_group_table()` | `(ws)` | `data.frame(index, groupName, n_samples)`. Isolates the CytoML accessor |
| `wsp_sample_files()` | `(ws, group_index)` | `data.frame(sampleID, name)` for one group |
| `resolve_group()` | `(group, group_table, where = "")` | `list(index, name, n_samples)`. `group` is an integer index **or** a group name |
| `resolve_workspace_specs()` | `(wsp_paths, labels, fcs_index, group, open_ws = NULL, read_samples = NULL)` | `list(specs, problems)` — never stops on failure. `open_ws`/`read_samples` are injectable, so the orchestration is testable without CytoML |
| `format_resolution_report()` | `(problems, fcs_root = "")` | One consolidated error message, or `""` when nothing is wrong |
| `build_workspace_gatingset()` | `(spec)` | A `GatingSet`. Passes a directory when the workspace's files share one, else a per-sample `data.frame(sampleID, file)` path map |
| `disambiguate_sample_labels()` | `(clean_names, wsp_labels)` | Display labels: unchanged when unique across the run, `"<wsp>: <sample>"` only for names that actually collide |

```r
source("src/loading/flowjo_workspace.R")

wsp    <- find_workspace_files("data/workspaces")   # recursive
labels <- workspace_labels(wsp)
index  <- index_fcs_dir("data/fcs")                 # recursive

res <- resolve_workspace_specs(wsp, labels, index, group = "All Samples")
report <- format_resolution_report(res$problems, fcs_root = "data/fcs")
if (nzchar(report)) stop(report)

gs_list <- lapply(res$specs, build_workspace_gatingset)
```

**Dependencies:** none for the discovery, indexing and labelling helpers (base R only, which is what makes them unit-testable in `test_data/test_umap_multi_wsp.R`). `wsp_group_table()`, `wsp_sample_files()` and `build_workspace_gatingset()` require `flowWorkspace` and `CytoML`.

---

## load_manual_gating.R

Loads manual gating labels from CSV files, optionally matched to FCS data.

### `load_manual_gating()`

```r
load_manual_gating(csv_dir, fcs_dir = csv_dir,
                   gating_pattern = "Manual gating",
                   label_dict = NULL, gating_tree = NULL,
                   output_dir = NULL)
```

Reads CSV files with binary gating columns, extracts per-cell labels, optionally flattens hierarchical labels using a gating tree, and translates labels via a dictionary.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_dir` | character | *required* | Directory with CSV files |
| `fcs_dir` | character | `csv_dir` | Directory with matching FCS files for event count validation |
| `gating_pattern` | character | `"Manual gating"` | Regex pattern identifying gating columns |
| `label_dict` | named list | `NULL` | Label translation dictionary |
| `gating_tree` | named numeric vector | `NULL` | Gating hierarchy depths for flattening |
| `output_dir` | character | `NULL` | If set, writes per-sample label CSVs with `cell_id` and `label` columns |

**Returns:** Named list of label vectors per sample. When `output_dir` is set, also writes CSV files.

**Input from:** CSV + FCS files on disk.
**Output to:** Downstream analysis, cluster comparison.

**Example:**

```r
source("load_manual_gating.R")
labels <- load_manual_gating(
  csv_dir        = "/data/gating_csvs",
  gating_pattern = "Manual gating",
  output_dir     = "/output/labels"
)
cat("Loaded labels for", length(labels), "samples\n")
```

**Dependencies:** `io_utils`, `fcs_utils`, `cluster_utils`

---

## create_flowframes.R

Creates flowFrame objects from expression matrices, and optionally merges them into a flowSet. This is a bridge function enabling CSV expression data to flow into flowFrame-based scripts.

### `create_flowframes_from_matrices()`

```r
create_flowframes_from_matrices(expr_list, markers_of_interest = NULL,
                                 save_fcs = FALSE, output_dir = NULL)
```

Converts a named list of expression matrices (or data.frames/data.tables) into flowFrame objects. Optionally subsets to specified markers and writes FCS files.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_list` | named list | *required* | List of expression matrices or data.frames |
| `markers_of_interest` | character vector | `NULL` | Markers to include; `NULL` keeps all |
| `save_fcs` | logical | `FALSE` | Write FCS files to disk |
| `output_dir` | character | `NULL` | Directory for saved FCS files |

**Returns:** Named list of flowFrame objects.

**Input from:** `load_csv_expression()` (loading/load_csv_expression.R), `combine_expression_matrices()` (analysis/expression_matrices.R).
**Output to:** `plot_umap()` (visualization/umap.R), `plot_density()` (visualization/density_plots.R), `homogenize_flowframes()` (harmonization/panel_homogenization.R), `downsample_flowframes()` (preprocessing/downsampling.R).

**Example:**

```r
source("create_flowframes.R")

expr_list <- list(
  sample_001 = matrix(rnorm(300), ncol = 6,
    dimnames = list(NULL, c("CD3","CD4","CD8","CD19","CD56","HLA-DR"))),
  sample_002 = matrix(rnorm(300), ncol = 6,
    dimnames = list(NULL, c("CD3","CD4","CD8","CD19","CD56","HLA-DR")))
)

ff_list <- create_flowframes_from_matrices(expr_list, save_fcs = TRUE,
                                            output_dir = "/output/fcs")
# Creates sample_001.fcs and sample_002.fcs
```

### `create_merged_flowset()`

```r
create_merged_flowset(flow_frames)
```

Merges a named list of flowFrames into a single flowSet. Validates that all frames have the same number of columns.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | named list | *required* | List of flowFrame objects |

**Returns:** A `flowSet` object.

**Input from:** `create_flowframes_from_matrices()`.
**Output to:** FlowSOM and other flowSet-based analyses.

**Dependencies:** `io_utils`, `flowframe_utils`
