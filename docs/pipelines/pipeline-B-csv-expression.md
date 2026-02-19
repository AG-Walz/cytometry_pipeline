[< Back to Documentation Index](../README.md)

# Pipeline B: CSV Expression Files — Marker Statistics

Pipeline B is the simplest pipeline in the framework. It takes pre-exported expression CSV files (e.g. from FlowJo, Cytolution, or any tool that produces per-cell expression matrices), computes per-marker statistics using linear mixed-effects models, and generates a clustered expression heatmap. There is no clustering, no UMAP, no differential testing — just straightforward marker-level analysis.

This pipeline is useful as a quick sanity check on exported expression data, or when you only need per-marker variance decomposition (how much of the variance is between-sample vs. within-sample).

---

## When to Use Pipeline B

Use this pipeline when:

- Your data is already exported as **expression CSVs** (one per sample-population combination)
- You want **per-marker statistics**: mean, between-sample variance, within-sample variance, ICC
- You don't need clustering, UMAP, or differential abundance testing
- You want a quick **heatmap** overview of expression patterns

If your CSVs also contain binary gating columns (0/1), use [Pipeline C](pipeline-C-csv-gating.md) instead.

---

## Pipeline Flow

```
  CSV files (one per sample × population)
       │
       ▼
  ┌────────────────────────────┐
  │ 1. load_csv_expression()   │  Read CSVs, auto-detect separator
  └─────────┬──────────────────┘
            ▼
  ┌────────────────────────────┐
  │ 2. run_marker_statistics() │  Per-marker lmer: mean, variance, ICC
  └─────────┬──────────────────┘
            ▼
  ┌────────────────────────────┐
  │ 3. plot_heatmap()          │  Clustered expression heatmap (pheatmap)
  └────────────────────────────┘
```

---

## Input

A directory containing semicolon-separated CSV files. Each file represents one sample-population combination, with filenames encoding sample and population: `{sample}_{population}.csv`. Every column is an expression marker — no metadata columns, no gating columns.

### Directory Layout

```
csv_expression/
├── sample_001_CD4_T_cells.csv
├── sample_001_CD8_T_cells.csv
├── sample_002_CD4_T_cells.csv
├── sample_002_CD8_T_cells.csv
├── sample_003_CD4_T_cells.csv
└── sample_003_CD8_T_cells.csv
```

### Example CSV Content

Each CSV is semicolon-separated with marker names as column headers and one row per cell:

```
CD3;CD4;CD8;CD19;CD56;HLA-DR
3.94;4.31;5.33;3.23;2.35;2.78
2.75;4.80;3.52;3.42;2.00;4.06
1.89;3.12;2.67;4.01;3.55;1.92
4.12;2.98;1.45;2.87;4.23;3.11
0.56;3.67;4.89;1.23;2.78;2.45
...
```

The separator is auto-detected (semicolon, comma, or tab). A warning is issued if files within the same folder use different separators.

### Filename Convention

The default regex for extracting sample and population from filenames is:

```r
filename_regex = "^(.+)_(.+)\\.csv$"
```

This means `sample_001_CD4_T_cells.csv` is parsed as:
- **Sample:** `sample_001`
- **Population:** `CD4_T_cells`

You can customize this regex in the config.

---

## Configuration

### Running the Pipeline

```bash
cd pipelines
Rscript pipeline_B_csv_expression.R --config ../config/pipeline_B_config.R
```

### Example Config

```r
# config/pipeline_B_config.R

params <- list(
  csv_folder = "/path/to/csv/expression/files",
  output_dir = "/path/to/output/pipeline_B_results",
  markers    = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),
  test_type  = "lmer",    # "lmer" (mixed-effects model)
  seed       = 42
)
```

### Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_folder` | character | *required* | Directory containing expression CSV files |
| `output_dir` | character | *required* | Root output directory |
| `markers` | character vector | *required* | Column names to treat as markers |
| `test_type` | character | `"lmer"` | Statistical method for marker analysis |
| `seed` | integer | `42` | Random seed for reproducibility |

---

## Steps in Detail

### Step 1: Load CSV Expression Files

```r
source(file.path(base_dir, "src", "loading", "load_csv_expression.R"))
expr_data <- load_csv_expression(
  input_dir = params$csv_folder
)
```

**What it does:** Reads every CSV file in the directory using `data.table::fread()` with automatic separator detection. Returns a named list of data.tables — one per file.

**Function signature:**

```r
load_csv_expression(input_dir, sep = NULL, max_memory_gb = 4)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with CSV files |
| `sep` | character | `NULL` | Separator; `NULL` for auto-detection |
| `max_memory_gb` | numeric | `4` | Memory budget (for very large files) |

**Returns:** Named list of `data.table` objects.

```r
# Example:
# expr_data$sample_001_CD4_T_cells  -> data.table with 50 rows x 6 cols
# expr_data$sample_002_CD8_T_cells  -> data.table with 50 rows x 6 cols
```

**Console output:**

```
--- Step 1: Loading CSV expression data ---
  Scanning: /path/to/csv_expression/
  Found 6 CSV files
  Loading sample_001_CD4_T_cells.csv (0.01 MB)...
  Loading sample_001_CD8_T_cells.csv (0.01 MB)...
  Loading sample_002_CD4_T_cells.csv (0.01 MB)...
  Loading sample_002_CD8_T_cells.csv (0.01 MB)...
  Loading sample_003_CD4_T_cells.csv (0.01 MB)...
  Loading sample_003_CD8_T_cells.csv (0.01 MB)...
Loaded 6 expression matrices
```

---

### Step 2: Marker Statistics

```r
source(file.path(base_dir, "src", "analysis", "marker_statistics.R"))
stats_result <- run_marker_statistics(
  input_folder  = params$csv_folder,
  output_folder = file.path(params$output_dir, "01_statistics"),
  markers       = params$markers,
  test_type     = params$test_type
)
```

**What it does:** For each population (parsed from filenames) and each marker, fits a linear mixed-effects model:

```
marker_expression ~ 1 + (1 | sample_id)
```

This decomposes variance into:
- **Between-sample variance**: How much does the marker's mean differ across samples?
- **Within-sample variance**: How much does the marker vary within a single sample?
- **ICC** (intraclass correlation coefficient): The proportion of total variance that is between-sample. An ICC near 0 means most variation is cell-to-cell; near 1 means most variation is sample-to-sample.

**Function signature:**

```r
run_marker_statistics(input_folder, output_folder,
                      filename_regex = "^(.+)_(.+)\\.csv$",
                      markers = NULL, test_type = "lmer",
                      n_workers = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_folder` | character | *required* | Directory with CSVs |
| `output_folder` | character | *required* | Output directory for results |
| `filename_regex` | character | `"^(.+)_(.+)\\.csv$"` | Regex to extract sample + population |
| `markers` | character vector | `NULL` | Markers to analyze; `NULL` uses all |
| `test_type` | character | `"lmer"` | `"lmer"` for mixed-effects model |

**Output** — `marker_statistics.csv`:

```
population;marker;mean;between_sample_var;within_sample_var;icc;n_cells;n_samples;method
CD4_T_cells;CD3;3.02;0.00;0.97;0.00;150;3;lmer
CD4_T_cells;CD4;3.08;0.003;1.10;0.003;150;3;lmer
CD4_T_cells;CD8;2.95;0.024;1.04;0.022;150;3;lmer
CD4_T_cells;CD19;3.02;0.008;1.02;0.008;150;3;lmer
CD4_T_cells;CD56;3.03;0.000;1.01;0.000;150;3;lmer
CD4_T_cells;HLA-DR;2.98;0.006;1.05;0.006;150;3;lmer
CD8_T_cells;CD3;2.92;0.001;0.98;0.001;150;3;lmer
...
```

**Reading the results:**

- **High ICC (>0.1):** This marker varies meaningfully between samples — worth investigating for biological signal
- **Near-zero ICC:** All variation is cell-to-cell within each sample — the marker behaves consistently across samples
- **method = "simple_mean"**: The mixed model failed to converge (too few samples or too little variance); a simple average was used as fallback

---

### Step 3: Visualization — Expression Heatmap

```r
source(file.path(base_dir, "src", "visualization", "heatmaps.R"))
plot_heatmap(
  expr_matrix = t(combined_matrix[sample(nrow(combined_matrix), min(500, nrow(combined_matrix))), ]),
  output_path = file.path(viz_dir, "expression_heatmap.png"),
  title       = "Expression Heatmap (sampled cells)"
)
```

**What it does:** Generates a clustered heatmap (using `pheatmap`) showing expression patterns across markers and cells. If the combined matrix has more than 500 cells, a random subset of 500 is drawn for readability. Rows (markers) and columns (cells) are hierarchically clustered.

**Function signature:**

```r
plot_heatmap(expr_matrix, output_path,
             title = "Expression Heatmap",
             cluster_rows = TRUE, cluster_cols = TRUE,
             scale = "row", color_palette = NULL,
             width = 10, height = 8)
```

**Output file:** `expression_heatmap.png` — a heatmap where each row is a marker and each column is a cell, colored by z-scored expression. Row and column dendrograms show hierarchical clustering.

---

## Output Directory Structure

```
output/pipeline_B_results/
├── 01_statistics/
│   └── marker_statistics.csv      # Per-population, per-marker statistics
└── 02_visualization/
    ├── expression_heatmap.png     # Clustered heatmap
    └── expression_heatmap.svg     # Vector format (if svg enabled)
```

---

## Tips

- **Filename convention matters.** The regex `^(.+)_(.+)\\.csv$` parses `sample_001_CD4_T_cells.csv` into sample=`sample_001` and population=`CD4_T_cells`. If your filenames don't follow this pattern, set `filename_regex` in the config.
- **Auto-detection of separators.** You can mix semicolon- and comma-separated files in the same directory (though a warning is issued). The loader handles each file independently.
- **Memory.** For very large CSVs (>1 GB per file), the `max_memory_gb` parameter triggers chunked reading.
- **Quick sanity check.** This pipeline is ideal for a fast look at your data before committing to a full Pipeline A run.
