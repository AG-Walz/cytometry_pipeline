[< Back to Documentation Index](../README.md)

# src/harmonization/ — Panel Comparison, Homogenization, Batch Correction

Three scripts for comparing marker panels across samples, unifying flowFrames to a common marker set, and correcting batch effects. These are essential when combining data from different instruments or staining panels.

---

## panel_comparison.R

Compares marker panels across multiple flowFrame objects to identify shared and divergent markers. Run this before homogenization to understand what your data looks like.

### Function

```r
compare_panels_detailed(flow_frames)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |

### Returns

```r
list(
  common_markers     = character(),  # Markers present in ALL samples
  all_markers        = character(),  # Union of all markers
  not_common         = character(),  # Markers missing from at least one sample
  presence_matrix    = matrix(),     # Logical: samples x markers
  panel_groups       = list(),       # Groups of samples with identical panels
  per_sample_markers = list()        # Markers per sample
)
```

### Example

```r
source("panel_comparison.R")
panel_info <- compare_panels_detailed(flow_frames)

cat("Common markers:", paste(panel_info$common_markers, collapse = ", "), "\n")
# Common markers: CD3, CD4, CD8, CD19

cat("Missing from some:", paste(panel_info$not_common, collapse = ", "), "\n")
# Missing from some: CD56, HLA-DR

cat("Panel groups:", length(panel_info$panel_groups), "\n")
# Panel groups: 2 (some samples have 6 markers, others have 4)
```

### Example Presence Matrix

```
             CD3  CD4  CD8  CD19  CD56  HLA-DR
sample_001   TRUE TRUE TRUE TRUE  TRUE  TRUE
sample_002   TRUE TRUE TRUE TRUE  TRUE  TRUE
sample_003   TRUE TRUE TRUE TRUE  FALSE TRUE
```

---

## panel_homogenization.R

Unifies flowFrames to a common set of markers in a consistent column order. Missing markers are filled with `NA`. This is the key step that makes it possible to combine data from different instruments or panels.

### Function

```r
homogenize_flowframes(flow_frames, sample_names, markers_of_interest,
                      markername_replacing = list(),
                      save_csv = FALSE, save_fcs = FALSE,
                      output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `markers_of_interest` | character vector | *required* | Target marker set |
| `markername_replacing` | named list | `list()` | Rename map: `list("old" = "new")` |
| `save_csv` | logical | `FALSE` | Save homogenized CSVs (`hg_*.csv`) |
| `save_fcs` | logical | `FALSE` | Save homogenized FCS files (`hg_*.fcs`) |
| `output_dir` | character | `NULL` | Directory for saved files |

### What It Does

1. Applies `markername_replacing` to rename markers across all flowFrames
2. For each flowFrame, subsets and reorders columns to match `markers_of_interest`
3. Fills missing markers with `NA` values
4. Ensures all flowFrames have identical column names in identical order

### Returns

Named list of homogenized `flowFrame` objects.

### Example

```r
source("panel_homogenization.R")

# Rename inconsistent markers, then homogenize to 4 markers
homogenized <- homogenize_flowframes(
  flow_frames         = my_frames,
  sample_names        = c("S1", "S2", "S3"),
  markers_of_interest = c("CD3", "CD4", "CD8", "PD1"),
  markername_replacing = list(`PD1 (CD279)` = "PD1", `CD8alpha` = "CD8"),
  save_csv            = TRUE,
  output_dir          = "/output/harmonized"
)
```

**Console output:**

```
  Homogenizing 3 samples to 4 markers
    Homogenizing S1 ...  50000 cells x 4 markers
    Homogenizing S2 ...  48000 cells x 4 markers  (PD1 was missing, filled with NA)
    Homogenizing S3 ...  52000 cells x 4 markers
```

**Output files (when `save_csv = TRUE`):**

```
output/harmonized/
├── hg_S1.csv    # 50000 rows x 4 columns (CD3;CD4;CD8;PD1)
├── hg_S2.csv    # 48000 rows x 4 columns (PD1 column is all NA)
└── hg_S3.csv
```

### Helper Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `translate_vector` | `(vector, dictionary)` | Apply a name mapping to a character vector |
| `replace_markernames` | `(flowframe, replace_map)` | Rename markers in a flowFrame |
| `reorder_exprmtx` | `(flowframe, marker_name_list)` | Reorder columns; fill missing with `NA` |
| `generate_cellID_list` | `(flow_frames)` | Generate unique cell IDs across all frames |

### Dependencies

`io_utils`, `flowframe_utils` (`explore_data`, `from_mtx_to_ffr`)

---

## batch_correction.R

Corrects batch effects across samples using Harmony. Use this when samples were processed in different batches and you want to remove batch-associated variation before clustering.

### Function

```r
run_batch_correction(flow_frames, sample_names, batch_labels,
                     markers_of_interest = NULL, do_pca = FALSE, ...)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `batch_labels` | character vector | *required* | Batch assignment per sample |
| `markers_of_interest` | character vector | `NULL` | Markers to correct; `NULL` uses all |
| `do_pca` | logical | `FALSE` | Run PCA before Harmony |
| `...` | | | Additional args to `harmony::HarmonyMatrix()` |

### Algorithm

Uses `harmony::HarmonyMatrix(vars_use = "batch")` to iteratively remove batch-associated variation while preserving biological signal.

### Returns

```r
list(
  corrected_matrix = matrix(),    # Batch-corrected expression matrix
  sample_labels    = character(), # Sample label per cell
  batch_labels     = character(), # Batch label per cell
  markers          = character()  # Marker names
)
```

### Example

```r
source("batch_correction.R")

result <- run_batch_correction(
  flow_frames         = my_frames,
  sample_names        = c("S1", "S2", "S3", "S4"),
  batch_labels        = c("batch1", "batch1", "batch2", "batch2"),
  markers_of_interest = c("CD3", "CD4", "CD8"),
  do_pca              = TRUE
)

# result$corrected_matrix is ready for downstream clustering
dim(result$corrected_matrix)
# [1] 200000      3
```

### Dependencies

`io_utils`, `flowframe_utils`. Requires `harmony` package.
