[< Back to Documentation Index](../README.md)

# Pipeline A: Raw FCS -> Full Analysis

## Purpose

Complete end-to-end analysis starting from raw FCS files. This is the most comprehensive pipeline, covering every stage from initial loading through quality control, transformation, panel harmonization, downsampling, unsupervised clustering, differential abundance analysis, and visualization.

## Pipeline Flow

```
Load FCS -> QC -> Transform -> Homogenize -> Downsample -> Cluster -> Analyze -> Visualize
```

```
  .fcs files
      |
      v
  [1] load_fcs_raw()
      |
      v
  [2] run_preprocessing()       # QC + transform
      |
      v
  [3] homogenize_flowframes()   # Unify marker panels
      |
      v
  [4] downsample_flowframes()   # Reduce event count
      |
      v
  [5] run_leiden_with_umap()    # Cluster + embed
      |
      v
  [6] split_by_cluster()        # Per-cluster analysis
      run_differential_abundance()
      |
      v
  [7] plot_umap()               # Final figures
      plot_heatmap()
```

## Input

A directory containing `.fcs` files. All files in the directory will be loaded and processed together. Files may have different marker panels -- the harmonization step resolves mismatches.

## Configuration

```bash
Rscript pipelines/pipeline_A_raw_fcs.R --config config/study_params/my_study.R
```

If `--config` is not provided, the pipeline falls back to `config/Parameters_template.R`.

## Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fcs_folder` | character | -- | Path to directory containing `.fcs` files |
| `truncate_max_range` | logical | `TRUE` | Truncate values exceeding the `$PnR` keyword range |
| `markers` | character vector | -- | Marker names to include in analysis |
| `qc_method` | character | `"PeacoQC"` | Quality control method: `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `transform_method` | character | `"arcsinh"` | Transformation: `"arcsinh"`, `"logicle"`, `"biexponential"`, or `"None"` |
| `cofactor` | numeric | `5` | Cofactor for arcsinh transformation |
| `downsample_method` | character | `"Maximal"` | `"Maximal"` (keep as many as smallest file) or `"Consistent"` (fixed N per file) |
| `downsample_n` | integer | `NULL` | Number of events per file (required when `downsample_method = "Consistent"`) |
| `n_neighbors` | integer | `15` | Number of neighbors for UMAP and Leiden graph construction |
| `resolution` | numeric | `0.8` | Leiden clustering resolution |
| `output_dir` | character | -- | Root directory for all output subdirectories |

## Steps in Detail

### Step 1: Load FCS Files

```r
# src/loading/load_fcs_raw.R
flow_frames <- load_fcs_raw(
  fcs_folder = params$fcs_folder,
  truncate_max_range = params$truncate_max_range
)
```

Reads all `.fcs` files from the specified folder into a list of `flowFrame` objects. When `truncate_max_range = TRUE`, values exceeding the per-channel `$PnR` keyword are capped.

### Step 2: Preprocessing (QC + Transform)

```r
# src/preprocessing/preprocess_pipeline.R
flow_frames <- run_preprocessing(
  flow_frames,
  qc_method    = params$qc_method,
  transform    = params$transform_method,
  cofactor     = params$cofactor
)
```

Applies quality control to remove aberrant events (margin events, signal anomalies), then transforms fluorescence/mass channels. Scatter channels are left untransformed.

### Step 3: Panel Harmonization

```r
# src/harmonization/panel_homogenization.R
flow_frames <- homogenize_flowframes(
  flow_frames,
  markers = params$markers
)
```

Unifies marker panels across files. Files missing specific markers receive `NA` columns; extra markers are dropped. Ensures all downstream steps operate on a consistent column set.

### Step 4: Downsampling

```r
# src/preprocessing/downsampling.R
flow_frames <- downsample_flowframes(
  flow_frames,
  method = params$downsample_method,
  n      = params$downsample_n
)
```

Reduces event counts to balance representation across samples and manage memory. `"Maximal"` equalizes to the smallest file; `"Consistent"` takes a fixed number from each file.

### Step 5: Leiden Clustering + UMAP

```r
# src/clustering/leiden.R
results <- run_leiden_with_umap(
  flow_frames,
  markers      = params$markers,
  n_neighbors  = params$n_neighbors,
  resolution   = params$resolution
)
```

Builds a shared nearest-neighbor graph, applies Leiden community detection, and computes a 2D UMAP embedding for visualization.

### Step 6: Analysis

```r
# src/analysis/cluster_splitting.R
clusters <- split_by_cluster(results$expr_matrix, results$cluster_ids)

# src/analysis/differential_abundance.R
da_results <- run_differential_abundance(
  clusters,
  condition_col = params$condition_col,
  sample_col    = params$sample_col
)
```

Splits the expression matrix by cluster assignment, then runs differential abundance testing (negative binomial GLMM or Fisher's exact test) to identify clusters enriched in specific conditions.

### Step 7: Visualization

```r
# src/visualization/umap.R
plot_umap(results$umap_coords, results$cluster_ids, output_dir = params$output_dir)

# src/visualization/heatmaps.R
plot_heatmap(results$expr_matrix, results$cluster_ids, output_dir = params$output_dir)
```

Generates UMAP scatter plots colored by cluster, condition, and individual markers, plus median-expression heatmaps summarizing cluster phenotypes.

## Output Directory Structure

```
output_dir/
  01_preprocessing/       # QC reports, transformation logs
  02_harmonization/       # Unified panel summary
  03_downsampled/         # Downsampled expression matrices
  04_clustering/          # Cluster assignments, Leiden graph
  05_analysis/            # DA results, enrichment tables
  06_visualization/       # UMAP plots, heatmaps, barplots
```

## Notes

- The pipeline logs elapsed time for each step and prints a total runtime summary at the end.
- For large datasets (>10M events total), consider using `downsample_method = "Consistent"` with a moderate `downsample_n` to keep memory usage manageable.
- If marker panels are identical across all files, the harmonization step is a no-op but still validates consistency.
