[< Back to Documentation Index](../README.md)

# src/harmonization/ -- Panel Comparison, Homogenization, Batch Correction

3 scripts for comparing marker panels across samples, unifying flowFrames to a common marker set, and correcting batch effects.

---

## panel_comparison.R

Compares marker panels across multiple flowFrame objects to identify shared and divergent markers.

### Main Function

```r
compare_panels_detailed(flow_frames)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |

### Output

```r
list(
  common_markers     = character(),  # Markers present in ALL samples
  all_markers        = character(),  # Union of all markers across samples
  not_common         = character(),  # Markers missing from at least one sample
  presence_matrix    = matrix(),     # Logical matrix: samples x markers
  panel_groups       = list(),       # Groups of samples with identical panels
  per_sample_markers = list()        # Named list of marker vectors per sample
)
```

### Logic

1. Extracts markers from each flowFrame via `explore_data()`
2. Computes the intersection (common markers) and union (all markers)
3. Builds a presence matrix indicating which markers appear in which samples
4. Groups samples that share identical marker panels

### Example

```r
panel_info <- compare_panels_detailed(flow_frames)

cat("Common markers:", paste(panel_info$common_markers, collapse = ", "), "\n")
cat("Missing from some samples:", paste(panel_info$not_common, collapse = ", "), "\n")
cat("Number of panel groups:", length(panel_info$panel_groups), "\n")
```

---

## panel_homogenization.R

Unifies flowFrames to a common set of markers in a consistent column order. Missing markers are filled with `NA`.

### Main Function

```r
homogenize_flowframes(flow_frames, sample_names, markers_of_interest,
                      markername_replacing = list(),
                      save_csv = FALSE, save_fcs = FALSE,
                      output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Corresponding sample identifiers |
| `markers_of_interest` | character vector | *required* | Target marker set (defines output columns) |
| `markername_replacing` | named list | `list()` | Rename map: `list(old_name = "new_name", ...)` |
| `save_csv` | logical | `FALSE` | Save homogenized expression matrices as CSVs (`hg_*.csv`) |
| `save_fcs` | logical | `FALSE` | Save homogenized flowFrames as FCS files (`hg_*.fcs`) |
| `output_dir` | character | `NULL` | Directory for saved files |

### Helper Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `translate_vector` | `(vector, dictionary)` | Applies a name mapping dictionary to a character vector |
| `replace_markernames` | `(flowframe, replace_map)` | Renames markers in a flowFrame using a mapping |
| `reorder_exprmtx` | `(flowframe, marker_name_list)` | Reorders columns to match `marker_name_list`; fills missing markers with `NA` |
| `generate_cellID_list` | `(flow_frames)` | Generates unique cell IDs across all flowFrames |

### Output

Named list of homogenized `flowFrame` objects, all with identical columns in the order specified by `markers_of_interest`.

### Example

```r
# Rename inconsistent marker names, then homogenize
rename_map <- list("CD3e" = "CD3", "CD8alpha" = "CD8a")

hg_frames <- homogenize_flowframes(
  flow_frames, sample_names,
  markers_of_interest = c("CD3", "CD4", "CD8a", "CD19"),
  markername_replacing = rename_map,
  save_csv = TRUE,
  output_dir = "output/harmonized"
)
```

### src/R_shared Dependencies

`io_utils`, `flowframe_utils` (`explore_data`, `from_mtx_to_ffr`)

---

## batch_correction.R

Corrects batch effects across samples using Harmony.

### Main Function

```r
run_batch_correction(flow_frames, sample_names, batch_labels,
                     markers_of_interest = NULL, do_pca = FALSE, ...)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `batch_labels` | character vector | *required* | Batch assignment per sample (same length as `sample_names`) |
| `markers_of_interest` | character vector | `NULL` | Markers to correct; `NULL` uses all common markers |
| `do_pca` | logical | `FALSE` | Whether to run PCA before Harmony |
| `...` | | | Additional arguments passed to `harmony::HarmonyMatrix()` |

### Algorithm

Uses `harmony::HarmonyMatrix(vars_use = "batch")` to iteratively remove batch-associated variation while preserving biological signal.

### Output

```r
list(
  corrected_matrix = matrix(),    # Batch-corrected expression matrix
  sample_labels    = character(), # Sample label per cell
  batch_labels     = character(), # Batch label per cell
  markers          = character()  # Marker names (columns of corrected_matrix)
)
```

### Example

```r
result <- run_batch_correction(
  flow_frames, sample_names,
  batch_labels = c("batch1", "batch1", "batch2", "batch2"),
  markers_of_interest = c("CD3", "CD4", "CD8a"),
  do_pca = TRUE
)

# result$corrected_matrix is ready for downstream clustering
```

### src/R_shared Dependencies

`io_utils`, `flowframe_utils`
