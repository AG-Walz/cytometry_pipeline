[< Back to Documentation Index](../README.md)

# Pipeline C: CSV Gating -> Cluster Analysis

## Purpose

Analyze CSV files that already contain binary gating or clustering columns alongside expression data. This is the standard pipeline for Cytolution exports and similar tools that produce per-event classification columns (0/1 values indicating gate membership). The pipeline splits events by cluster, runs differential abundance or enrichment analysis, and generates UMAP and barplot visualizations.

## Pipeline Flow

```
Load gating CSV -> Split clusters -> GLMM / Enrichment -> UMAP / Barplots
```

```
  CSV files (expression + gating columns)
      |
      v
  [1] load_csv_gating()         # Column classification + marker selection
      |
      v
  [2] split_by_cluster()        # Separate events by gate/cluster
      |
      v
  [3] run_differential_abundance()   # or run_cluster_enrichment()
      |
      v
  [4] plot_umap() / barplots()  # Visualization
```

## Input

CSV files containing both expression columns and binary gating columns (0/1 values). Optionally, matching FCS files can be provided for validation. Files are semicolon-separated by default with auto-detection.

Example CSV structure:

```
Index;FSC-A;SSC-A;CD3;CD4;CD8;CD3+ (Manual gating);CD4+ T cells (Clustering annotated);Clustering raw;UMAP_1;UMAP_2
1;65432;23456;2.31;1.87;0.42;1;0;3;-5.2;3.1
2;54321;34567;0.11;0.05;0.23;0;0;7;2.8;-1.4
...
```

## Column Classification

The `classify_csv_columns()` function automatically categorizes every column in the input:

| Category | Detection Logic | Examples |
|----------|----------------|----------|
| **Index** | Column named `Index` or row-number-like integers | `Index` |
| **Scatter** | Names matching `FSC-*`, `SSC-*` | `FSC-A`, `SSC-H` |
| **Markers** | Numeric columns not matching other categories | `CD3`, `CD4`, `CD8`, `CD19` |
| **Gating** | Binary (0/1) columns matching `gating_pattern` | `CD3+ (Manual gating)`, `Clustering raw` |
| **Dim. reduction** | Names matching `UMAP_*`, `tSNE_*`, `PC_*` | `UMAP_1`, `UMAP_2` |
| **Autofluorescence** | Names matching `AF*`, `Autofluorescence` | `AF488`, `AF700` |

The default `gating_pattern` regex:

```r
gating_pattern <- "\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$"
```

## Configuration

```bash
Rscript pipelines/pipeline_C_csv_gating.R --config config/study_params/my_study.R
```

If `--config` is not provided, the pipeline falls back to `config/Parameters_template.R`.

## Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_folder` | character | -- | Path to directory containing gating CSV files |
| `output_dir` | character | -- | Root directory for output subdirectories |
| `condition_col` | character | `"condition"` | Column identifying the experimental condition/group |
| `sample_col` | character | `"sample_id"` | Column identifying the sample of origin |
| `analysis_method` | character | `"enrichment"` | `"enrichment"` (Fisher's exact) or `"glmm"` (mixed-effects logistic regression) |
| `gating_pattern` | character | *(see above)* | Regex pattern to identify gating columns |
| `markers` | character vector | *(auto-detected)* | Override auto-detected marker columns |

## Steps in Detail

### Step 1: Load CSV Gating Files

```r
# loading/load_csv_gating.R
gating_data <- load_csv_gating(
  csv_folder     = params$csv_folder,
  gating_pattern = params$gating_pattern,
  markers        = params$markers
)
```

Reads all CSV files, runs `classify_csv_columns()` to separate expression data from gating assignments, and validates binary column integrity (warns if gating columns contain values other than 0 and 1). Users can override auto-detected markers with an explicit list.

### Step 2: Split by Cluster

```r
# analysis/cluster_splitting.R
clusters <- split_by_cluster(
  gating_data$expr_matrix,
  gating_data$gating_columns
)
```

Separates events into per-cluster subsets based on the binary gating columns. Each gate column produces a positive and negative subset. For the `Clustering raw` column, events are split by their integer cluster label.

### Step 3: Differential Abundance or Enrichment

```r
# analysis/differential_abundance.R
da_results <- run_differential_abundance(
  clusters,
  condition_col = params$condition_col,
  sample_col    = params$sample_col
)

# -- or --

# analysis/cluster_enrichment.R
enrich_results <- run_cluster_enrichment(
  clusters,
  condition_col = params$condition_col
)
```

When `analysis_method = "enrichment"`, Fisher's exact test is applied per cluster to test for over-/under-representation in each condition. When `analysis_method = "glmm"`, a mixed-effects logistic regression model accounts for sample-level random effects.

### Step 4: Visualization

```r
# visualization/umap.R
plot_umap(gating_data$umap_coords, gating_data$gating_columns, output_dir = params$output_dir)

# visualization/barplots.R
barplots(clusters, condition_col = params$condition_col, output_dir = params$output_dir)
```

UMAP plots are colored by gating assignment if dimension reduction columns are present in the input. Barplots show cluster proportions per condition and per sample.

## Output Directory Structure

```
output_dir/
  01_loaded/      # Column classification report, merged expression + gating table
  02_clusters/    # Per-cluster expression subsets
  03_analysis/    # DA or enrichment results, p-value tables
```

## Notes

- The column classification step is critical: review the printed summary to verify that gating columns were detected correctly. Adjust `gating_pattern` if your column naming convention differs.
- If UMAP coordinates are not present in the CSV, the visualization step skips UMAP plots and produces only barplots and heatmaps.
- FCS files are optional; they are used only for metadata validation (event count cross-check) when provided alongside CSVs.
