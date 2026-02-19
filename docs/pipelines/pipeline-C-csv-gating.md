[< Back to Documentation Index](../README.md)

# Pipeline C: CSV Gating Data — Cluster Analysis

Pipeline C processes CSV files that contain both expression data and binary gating columns. This is the standard pipeline for Cytolution exports and similar tools that produce per-event classification columns (0/1 values indicating gate membership). The pipeline loads the gating data, splits cells by cluster, tests for enrichment between experimental conditions, and generates composition barplots and abundance visualizations.

---

## When to Use Pipeline C

Use this pipeline when:

- Your data is **CSV files with gating columns** (binary 0/1 membership)
- The CSVs come from **Cytolution**, manual gating exports, or any tool that produces `(Clustering annotated)` / `(Manual gating)` columns
- You want **Fisher's exact test** or **mixed-effects logistic regression** for enrichment
- You want **abundance barplots** showing cluster proportions across samples

If your CSVs contain only expression columns (no gating), use [Pipeline B](pipeline-B-csv-expression.md) instead. If you have raw FCS files, use [Pipeline A](pipeline-A-raw-fcs.md).

---

## Pipeline Flow

```
  CSV + FCS file pairs (from Cytolution / manual gating)
       │
       ▼
  ┌──────────────────────────────┐
  │ 1. load_csv_gating()         │  Classify columns, split per population
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 2. split_by_cluster()        │  Count cells per cluster per sample
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 3. run_cluster_enrichment()  │  Fisher's exact / mixed-effects between conditions
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 4. plot_composition_barplot()│  Cluster composition barplot
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 5. build_tracking_matrix()   │  Abundance stacked / proportional / heatmap
  │    plot_abundance()          │
  └──────────────────────────────┘
```

---

## Input

### CSV + FCS File Pairs

A directory containing paired CSV and FCS files. The CSVs have expression columns plus binary gating columns (0/1 values). Each CSV should have a matching FCS file with the same base name — the FCS is used only for event count validation.

```
csv_gating/
├── sample_001.csv        # Expression + gating columns
├── sample_001.fcs        # Matching FCS (for validation)
├── sample_002.csv
├── sample_002.fcs
├── sample_003.csv
└── sample_003.fcs
```

### Example CSV Content

Each CSV is semicolon-separated with expression markers and binary gating columns:

```
CD3;CD4;CD8;CD19;CD56;HLA-DR;CD4+ (Clustering annotated);CD8+ (Clustering annotated);Clustering raw
2.46;4.01;2.15;1.29;2.37;1.22;0;0;1
2.57;4.40;3.12;4.38;4.07;2.33;0;0;0
3.89;2.12;4.56;1.87;3.21;3.45;1;0;1
0.98;1.56;2.78;3.45;1.23;2.67;0;1;0
1.34;3.78;1.99;2.56;4.12;3.89;1;0;1
...
```

### Column Auto-Classification

The `classify_csv_columns()` function automatically categorizes every column:

| Category | Detection Pattern | Examples |
|----------|-------------------|----------|
| **Index** | Column named `Index` or `index` | `Index` |
| **Scatter/Time** | `FSC-*`, `SSC-*`, `TIME`, `Time` | `FSC-A`, `SSC-H`, `TIME \| Time` |
| **Expression** | Numeric columns not matching other patterns | `CD3`, `CD4 \| BV421.A` |
| **Gating** | Ends with `(Clustering annotated)`, `(Manual gating)`, or is `Clustering raw` | `CD4+ (Clustering annotated)` |
| **Dim. reduction** | `cyto_2d_*`, `UMAP_*`, `tSNE_*` | `cyto_2d_Samples_All_X` |
| **Autofluorescence** | Marker name equals channel name | `AF.P1 \| AF.P1.A` |

The default gating pattern regex:

```r
gating_pattern = "\\(Manual gating\\)$|\\(Clustering annotated\\)$|^Clustering raw$"
```

### Conditions CSV

A separate file mapping samples to experimental groups (stored **outside** the data directory to avoid being picked up as a data file):

```
file_name;condition
sample_001;Control
sample_002;Control
sample_003;Treatment
```

---

## Configuration

### Running the Pipeline

```bash
cd pipelines
Rscript pipeline_C_csv_gating.R --config ../config/pipeline_C_config.R
```

### Example Config

```r
# config/pipeline_C_config.R

params <- list(
  # --- Paths ---
  csv_folder             = "/path/to/csv/gating/files",
  output_dir             = "/path/to/output/pipeline_C_results",
  conditions_csv         = "/path/to/conditions.csv",    # outside csv_folder!

  # --- Column detection ---
  gating_pattern         = "\\(Clustering annotated\\)$|^Clustering raw$",
  cluster_regex          = "Clustering annotated",
  processing_mode        = "manual",    # "automatic" (GUI), "manual", "expression_only"

  # --- Analysis ---
  analysis_method        = "enrichment",    # "enrichment" (Fisher) or "glmm" (NB-GLMM)
  condition_levels       = c("Control", "Treatment"),
  stat_method            = "pooled",        # "pooled" or "mixed_effects"
  p_adjust_method        = "BH",

  # --- Abundance plotting ---
  exclude_patterns       = c(),             # patterns to exclude from tracking
  cytolution_metadata_file = NULL,          # optional metadata CSV for grouping

  # --- Reproducibility ---
  seed                   = 42
)
```

### Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_folder` | character | *required* | Directory with CSV + FCS files |
| `conditions_csv` | character | *required* | Path to conditions file (keep outside `csv_folder`) |
| `gating_pattern` | character | see above | Regex to identify gating columns |
| `cluster_regex` | character | `"Clustering annotated"` | Regex for cluster label columns |
| `processing_mode` | character | `"automatic"` | `"automatic"` (GUI), `"manual"`, or `"expression_only"` |
| `analysis_method` | character | `"enrichment"` | `"enrichment"` (Fisher) or `"glmm"` (NB-GLMM) |
| `stat_method` | character | `"mixed_effects"` | `"pooled"` or `"mixed_effects"` for enrichment |
| `condition_levels` | character vector | -- | Two-element vector: `c("Control", "Treatment")` |
| `exclude_patterns` | character vector | `c()` | Patterns to exclude from tracking matrix |
| `cytolution_metadata_file` | character | `NULL` | Optional metadata CSV for grouped abundance plots |

**Important:** Keep the `conditions_csv` file **outside** the `csv_folder` directory. The loader scans all CSVs in the data directory, and having a conditions CSV there will cause errors (it has no matching FCS file).

---

## Steps in Detail

### Step 1: Load CSV Gating Files

```r
source(file.path(base_dir, "src", "loading", "load_csv_gating.R"))
gating_data <- load_csv_gating(
  input_fcs_path  = params$csv_folder,
  input_csv_path  = params$csv_folder,
  output_path     = file.path(params$output_dir, "01_loaded"),
  gating_pattern  = params$gating_pattern,
  processing_mode = params$processing_mode
)
```

**What it does:**

1. Scans the directory for CSV and FCS file pairs
2. For each CSV, runs `classify_csv_columns()` to separate expression data from gating columns
3. Validates that CSV row counts match FCS event counts
4. Splits each CSV by population — one output CSV per population

**Function signature:**

```r
load_csv_gating(input_fcs_path, input_csv_path = input_fcs_path,
                output_path, gating_pattern,
                csv_label_dict = NULL, gating_tree = NULL,
                processing_mode = "automatic", output_format = "wide",
                chunk_size = 100000)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `input_fcs_path` | character | Directory with FCS files |
| `input_csv_path` | character | Directory with CSV files (defaults to same as FCS) |
| `output_path` | character | Output directory for split files |
| `gating_pattern` | character | Regex identifying gating columns |
| `processing_mode` | character | `"automatic"` (GUI), `"manual"`, or `"expression_only"` |

**Returns:**

```r
list(
  output_dir  = "/path/to/output/01_loaded",
  results     = list(...),     # Per-file processing results
  n_processed = 3,
  n_total     = 3
)
```

**Output** — per-population CSVs in `01_loaded/`:

```
01_loaded/
├── sample_001_CD4_.csv          # Only CD4+ cells from sample_001
├── sample_001_CD8_.csv          # Only CD8+ cells from sample_001
├── sample_001_Clustering_raw.csv
├── sample_002_CD4_.csv
├── sample_002_CD8_.csv
├── sample_002_Clustering_raw.csv
├── sample_003_CD4_.csv
├── sample_003_CD8_.csv
└── sample_003_Clustering_raw.csv
```

Each per-population CSV contains only the expression columns for cells belonging to that population.

---

### Step 2: Split by Cluster

```r
source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
split_result <- split_by_cluster(
  input_dir     = params$csv_folder,     # reads from ORIGINAL CSVs (with gating columns)
  output_dir    = file.path(params$output_dir, "02_clusters"),
  cluster_regex = params$cluster_regex
)
```

**What it does:** Reads the **original** CSV files (which still have gating columns), identifies cluster columns matching `cluster_regex`, and produces:

1. Per-cluster CSV files for each sample
2. A tracking matrix counting cells per cluster per sample

**Function signature:**

```r
split_by_cluster(input_dir, output_dir,
                 cluster_regex = "Clustering annotated",
                 markers_of_interest = NULL,
                 min_cells_per_cluster = 10)
```

**Returns:**

```r
list(
  tracking_matrix  = data.frame(),   # Sample-by-cluster cell counts
  cluster_summary  = data.frame(),
  output_files     = character(),
  markers_used     = character()
)
```

**Output** — `cluster_tracking_matrix.csv`:

```
file;CD4+;CD8+
sample_001;20;15
sample_002;20;15
sample_003;20;15
```

This matrix shows how many cells from each sample belong to each cluster. It is the foundation for all abundance analyses and plots.

---

### Step 3: Statistical Analysis (Enrichment or GLMM)

The analysis method is controlled by `params$analysis_method`:

#### Option A: Fisher's Exact Test (enrichment)

```r
source(file.path(base_dir, "src", "analysis", "cluster_enrichment.R"))
enrichment <- run_cluster_enrichment(
  data_dir         = params$csv_folder,
  conditions_csv   = params$conditions_csv,
  output_dir       = file.path(params$output_dir, "03_analysis"),
  cluster_regex    = params$cluster_regex,
  condition_levels = params$condition_levels,
  stat_method      = params$stat_method,
  p_adjust_method  = params$p_adjust_method,
  seed             = params$seed
)
```

**Function signature:**

```r
run_cluster_enrichment(data_dir, conditions_csv, output_dir,
                       cluster_regex = "Clustering annotated",
                       condition_levels = c("pCR", "npCR"),
                       stat_method = "mixed_effects",
                       min_cluster_frac = 0.001,
                       p_adjust_method = "BH",
                       resume = FALSE, seed = 123)
```

**Output** — `enrichment_results.csv`:

```
cluster;in_A;not_in_A;in_B;not_in_B;odds_ratio;p_value;ci_lower;ci_upper;method;p_adjusted;significant;log2_enrichment
CD4+;40;60;20;30;1.00;1.00;0.47;2.13;fisher;1.00;FALSE;0.00
CD8+;30;70;15;35;1.00;1.00;0.45;2.27;fisher;1.00;FALSE;0.00
```

Column meanings:

| Column | Meaning |
|--------|---------|
| `in_A` / `not_in_A` | Cells in/out of this cluster in condition A (first level) |
| `in_B` / `not_in_B` | Cells in/out of this cluster in condition B (second level) |
| `odds_ratio` | Fisher's exact test odds ratio |
| `p_value` | Raw p-value from Fisher's exact test |
| `ci_lower` / `ci_upper` | 95% confidence interval for odds ratio |
| `p_adjusted` | Benjamini-Hochberg adjusted p-value |
| `significant` | `TRUE` if `p_adjusted < 0.05` |
| `log2_enrichment` | log2 of the proportion ratio between conditions |

#### Option B: Negative Binomial GLMM (differential abundance)

```r
source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))
abundance <- run_differential_abundance(
  data_dir         = params$csv_folder,
  conditions_csv   = params$conditions_csv,
  output_dir       = file.path(params$output_dir, "03_analysis"),
  cluster_regex    = params$cluster_regex,
  condition_levels = params$condition_levels,
  min_cluster_size = params$min_cluster_size,
  p_adjust_method  = params$p_adjust_method,
  seed             = params$seed
)
```

Output: `glmm_results.csv` with columns `cluster`, `estimate`, `std_error`, `z_value`, `p_value`, `log2fc`, `odds_ratio`, `method`, `p_adjusted`, `significant`.

---

### Step 4: Composition Barplot

```r
source(file.path(base_dir, "src", "visualization", "barplots.R"))
plot_composition_barplot(
  data       = long_data,
  output_dir = file.path(params$output_dir, "04_visualization"),
  x_col      = "cluster",
  fill_col   = "file",
  title      = "Cell Composition by Cluster"
)
```

**What it does:** Converts the tracking matrix to long format and generates a stacked barplot showing how many cells from each sample belong to each cluster.

**Function signature:**

```r
plot_composition_barplot(data, output_dir, x_col = "cluster",
                         fill_col = "sample", position = "stack",
                         title = "Cell Composition")
```

**Output:** `barplot_cluster_file.png` — a stacked barplot with clusters on the x-axis and colored segments for each sample.

---

### Step 5: Abundance Analysis

```r
source(file.path(base_dir, "src", "analysis", "abundance_plotting.R"))

# Build tracking matrix from original CSVs
tracking_matrix <- build_tracking_matrix(
  input_dir        = params$csv_folder,
  cluster_regex    = params$cluster_regex,
  exclude_patterns = params$exclude_patterns,
  output_dir       = file.path(params$output_dir, "05_abundance")
)

# Generate abundance plots
plot_abundance(
  output_dir      = file.path(params$output_dir, "05_abundance"),
  tracking_matrix = tracking_matrix,
  metadata_file   = params$cytolution_metadata_file
)
```

**What it does:**

1. `build_tracking_matrix()` scans the input directory, counts cells per cluster per sample, and saves the matrix
2. `plot_abundance()` generates multiple visualization styles from the tracking matrix

**Function signatures:**

```r
build_tracking_matrix(input_dir, cluster_regex, exclude_patterns, output_dir)
```

```r
plot_abundance(tracking_matrix_path = NULL, output_dir,
               tracking_matrix = NULL, metadata_file = NULL)
```

**Output** — four plots in `05_abundance/`:

| File | Description |
|------|-------------|
| `abundance_stacked.png` | Stacked barplot — samples sorted by total cells, clusters stacked |
| `abundance_proportional.png` | Proportional barplot with Ward.D2 dendrogram for cluster ordering |
| `total_cells_per_cluster.png` | Simple bar chart showing total cells per cluster across all samples |
| `cluster_proportion_heatmap.png` | Biclustered heatmap (inferno palette) showing cluster proportions |
| `cluster_tracking_matrix.csv` | The tracking matrix itself |

If `metadata_file` is provided (a Cytolution metadata sheet), additional grouped plots are generated based on metadata columns (e.g. grouped by treatment, response, etc.).

---

## Output Directory Structure

```
output/pipeline_C_results/
├── 01_loaded/                           # Per-population expression CSVs
│   ├── sample_001_CD4_.csv
│   ├── sample_001_CD8_.csv
│   ├── sample_001_Clustering_raw.csv
│   └── ...
├── 02_clusters/
│   ├── cluster_tracking_matrix.csv      # Sample × cluster cell counts
│   ├── sample_001_CD4+.csv             # Expression data for CD4+ cells
│   └── ...
├── 03_analysis/
│   └── enrichment_results.csv           # Fisher / mixed-effects results
├── 04_visualization/
│   └── barplot_cluster_file.png         # Composition barplot
└── 05_abundance/
    ├── cluster_tracking_matrix.csv      # Rebuilt tracking matrix
    ├── abundance_stacked.png            # Stacked barplot by total count
    ├── abundance_proportional.png       # Proportional barplot + dendrogram
    ├── total_cells_per_cluster.png      # Total cells per cluster
    └── cluster_proportion_heatmap.png   # Biclustered proportion heatmap
```

---

## Tips

- **Conditions file location:** Store the conditions CSV **outside** the `csv_folder` directory. The loader scans all CSVs in the data directory, and an extra CSV with no matching FCS file will cause an error.
- **Custom gating patterns:** If your gating columns use a different naming convention (e.g. `Gate: CD4+` instead of `CD4+ (Clustering annotated)`), update `gating_pattern` in the config.
- **Fisher vs GLMM:** Use `analysis_method = "enrichment"` with `stat_method = "pooled"` for simple proportional testing (works well with small sample sizes). Switch to `"glmm"` when you have enough samples per condition (>5) to benefit from the random-effects structure.
- **Excluding clusters:** Use `exclude_patterns = c("Clustering raw")` to skip raw numeric cluster columns from the tracking matrix.
- **Metadata grouping:** Set `cytolution_metadata_file` to a Cytolution metadata sheet CSV to get additional grouped abundance plots (by treatment, response, etc.).
