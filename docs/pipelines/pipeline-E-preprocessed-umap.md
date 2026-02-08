[< Back to Documentation Index](../README.md)

# Pipeline E: Preprocessed FCS -> UMAP / Clustering

## Purpose

Quick unsupervised clustering and UMAP visualization for FCS files that have already been quality-controlled and transformed. This pipeline skips QC and transformation entirely, going straight from loading to panel harmonization, downsampling, Leiden clustering, and UMAP embedding. Use it when your FCS files come from an upstream preprocessing tool and you want to explore cluster structure without re-running QC.

## Pipeline Flow

```
Load FCS (no QC) -> Homogenize -> Downsample -> Cluster -> UMAP
```

```
  Preprocessed .fcs files
      |
      v
  [1] load_fcs_raw(truncate_max_range=FALSE)
      |
      v
  [2] homogenize_flowframes()
      |
      v
  [3] downsample_flowframes()
      |
      v
  [4] run_leiden_with_umap()
      |
      v
  [5] plot_umap()
```

## Input

Already preprocessed `.fcs` files in a directory. These files should have undergone QC (removal of margin events, doublets, dead cells) and appropriate transformation (arcsinh, logicle, etc.) before being passed to this pipeline.

## Configuration

```bash
Rscript pipelines/pipeline_E_preprocessed_umap.R --config config/study_params/my_study.R
```

If `--config` is not provided, the pipeline falls back to `config/Parameters_template.R`.

## Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fcs_folder` | character | -- | Path to directory containing preprocessed `.fcs` files |
| `markers` | character vector | -- | Marker names to include in clustering and UMAP |
| `downsample_method` | character | `"Maximal"` | `"Maximal"` (equalize to smallest file) or `"Consistent"` (fixed N per file) |
| `downsample_n` | integer | `NULL` | Events per file (required when `downsample_method = "Consistent"`) |
| `n_neighbors` | integer | `15` | Number of neighbors for UMAP and Leiden graph construction |
| `resolution` | numeric | `0.8` | Leiden clustering resolution |
| `output_dir` | character | -- | Root directory for output subdirectories |

## Key Difference from Pipeline A

Pipeline E differs from Pipeline A in several important ways:

| Aspect | Pipeline A | Pipeline E |
|--------|-----------|-----------|
| `truncate_max_range` | `TRUE` | `FALSE` |
| QC step | Yes (`PeacoQC` / `flowAI`) | No |
| Transform step | Yes (`arcsinh` / `logicle`) | No |
| Analysis step | Yes (differential abundance) | No |
| Focus | Full end-to-end analysis | Unsupervised discovery |

Pipeline E sets `truncate_max_range = FALSE` because the input files are already preprocessed -- truncating would potentially discard valid transformed values that exceed the original instrument range keywords.

## Steps in Detail

### Step 1: Load FCS Files (No QC)

```r
# src/loading/load_fcs_raw.R
flow_frames <- load_fcs_raw(
  fcs_folder          = params$fcs_folder,
  truncate_max_range  = FALSE
)
```

Loads all `.fcs` files without applying any quality control or transformation. The `truncate_max_range = FALSE` flag preserves all values as-is from the preprocessed files.

### Step 2: Panel Harmonization

```r
# src/harmonization/panel_homogenization.R
flow_frames <- homogenize_flowframes(
  flow_frames,
  markers = params$markers
)
```

Unifies marker panels across files. Even preprocessed files may originate from different instrument configurations or staining panels. Missing markers are filled with `NA`; extra markers are dropped.

### Step 3: Downsampling

```r
# src/preprocessing/downsampling.R
flow_frames <- downsample_flowframes(
  flow_frames,
  method = params$downsample_method,
  n      = params$downsample_n
)
```

Reduces event counts for computational tractability. Particularly important when preprocessed files retain large event counts.

### Step 4: Leiden Clustering + UMAP

```r
# src/clustering/leiden.R
results <- run_leiden_with_umap(
  flow_frames,
  markers      = params$markers,
  n_neighbors  = params$n_neighbors,
  resolution   = params$resolution
)
```

Constructs a shared nearest-neighbor graph from the marker expression matrix, applies Leiden community detection to identify clusters, and computes a 2D UMAP embedding.

### Step 5: UMAP Visualization

```r
# src/visualization/umap.R
plot_umap(
  results$umap_coords,
  results$cluster_ids,
  output_dir = params$output_dir
)
```

Generates UMAP scatter plots colored by cluster assignment, sample identity, and individual marker expression.

## Output Directory Structure

```
output_dir/
  01_harmonized/      # Unified panel summary, harmonized expression matrices
  02_downsampled/     # Downsampled expression matrices
  03_clustering/      # Cluster assignments, Leiden parameters, UMAP coordinates
  04_visualization/   # UMAP plots (by cluster, by sample, by marker)
```

## Notes

- This pipeline is focused on unsupervised discovery: it does not run differential abundance, enrichment, or any statistical comparison between conditions. Use Pipeline A or Pipeline C for condition-aware analysis.
- Since no QC is applied, the quality of results depends entirely on the preprocessing applied to the input FCS files. Verify that dead cells, doublets, and debris have been removed upstream.
- For rapid iteration on clustering parameters (`resolution`, `n_neighbors`), this pipeline's lean structure makes it faster than Pipeline A.
- The pipeline logs elapsed time for each step.
