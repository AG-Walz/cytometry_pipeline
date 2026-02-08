[< Back to Documentation Index](../README.md)

# Pipeline B: CSV Expression -> Statistics

## Purpose

Lightweight analysis pipeline for expression CSV files. This pipeline bypasses all FCS-specific processing (QC, transformation, panel harmonization) and goes directly from tabular expression data to marker statistics and visualization. Ideal when data has already been exported from another tool or when working with non-FCS data sources.

## Pipeline Flow

```
Load CSV -> Marker Statistics -> Density / Heatmap
```

```
  CSV files (expression matrices)
      |
      v
  [1] load_csv_expression()
      |
      v
  [2] run_marker_statistics()
      |
      v
  [3] density_plots() / heatmaps()
```

## Input

A directory of CSV files containing expression columns. Files are expected to be semicolon-separated by default, but the loader auto-detects the delimiter (semicolon, comma, or tab). Each file typically represents one sample.

Example CSV structure:

```
sample_id;CD3;CD4;CD8;CD19;CD56;condition
S001;2.31;1.87;0.42;0.11;0.05;Treatment
S001;1.98;0.23;2.54;0.08;0.12;Treatment
...
```

## Configuration

```bash
Rscript pipelines/pipeline_B_csv_expression.R --config config/study_params/my_study.R
```

If `--config` is not provided, the pipeline falls back to `config/Parameters_template.R`.

## Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_folder` | character | -- | Path to directory containing CSV expression files |
| `markers` | character vector | -- | Column names to treat as markers for statistical analysis |
| `condition_col` | character | `"condition"` | Column name identifying the experimental condition/group |
| `sample_col` | character | `"sample_id"` | Column name identifying the sample of origin |
| `output_dir` | character | -- | Root directory for output subdirectories |

## Steps in Detail

### Step 1: Load CSV Expression Files

```r
# src/loading/load_csv_expression.R
expr_data <- load_csv_expression(
  csv_folder = params$csv_folder,
  markers    = params$markers
)
```

Reads all CSV files from the specified directory using `data.table::fread()` with automatic delimiter detection. Validates that the specified marker columns exist and coerces them to numeric. Non-marker columns (metadata, sample identifiers) are preserved as annotations.

### Step 2: Marker Statistics

```r
# src/analysis/marker_statistics.R
stats <- run_marker_statistics(
  expr_data,
  markers       = params$markers,
  condition_col = params$condition_col,
  sample_col    = params$sample_col
)
```

Computes per-marker summary statistics (mean, median, SD, IQR) grouped by condition and sample. Runs statistical comparisons between conditions (Wilcoxon rank-sum or t-test depending on group count and normality).

### Step 3: Visualization

```r
# src/visualization/density_plots.R
density_plots(expr_data, markers = params$markers, output_dir = params$output_dir)

# src/visualization/heatmaps.R
plot_heatmap(expr_data, markers = params$markers, output_dir = params$output_dir)
```

Generates per-marker density plots (overlaid by condition) and summary heatmaps showing median expression across markers and samples.

## Output Directory Structure

```
output_dir/
  01_loaded/        # Merged expression table, column summary
  02_statistics/    # Per-marker stats tables, comparison results, density plots, heatmaps
```

## Notes

- This is the simplest pipeline in the framework -- no FCS processing overhead, no clustering, no dimensional reduction.
- Useful as a quick sanity check on exported expression data before committing to a full Pipeline A run.
- The auto-detection of delimiters handles mixed-format inputs gracefully; a warning is issued if files within the same folder use different separators.
