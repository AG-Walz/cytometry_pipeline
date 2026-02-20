[< Back to Documentation Index](../README.md)

# Pipeline A: Raw FCS Files — Full End-to-End Analysis

Pipeline A is the most comprehensive pipeline in the framework. It takes raw, unprocessed `.fcs` files directly from the cytometer and runs them through every stage of analysis: quality control, mathematical transformation, panel harmonization across files, downsampling, unsupervised Leiden clustering, differential abundance testing between experimental conditions, and UMAP visualization. By the end, you have cluster assignments, statistical test results, and publication-ready plots.

---

## When to Use Pipeline A

Use this pipeline when:

- Your data is **raw FCS** straight from the instrument (not yet QC'd or transformed)
- You want **unsupervised clustering** (Leiden) rather than manual gating
- You want **differential abundance** testing between experimental conditions
- Your samples may have **different marker panels** that need harmonization

If your FCS files have already been QC'd and transformed, use [Pipeline E](pipeline-E-preprocessed-umap.md) instead — it skips preprocessing entirely.

---

## Pipeline Steps

1. **`load_fcs_raw()`** — Read FCS files into flowFrame objects
2. **`run_preprocessing()`** — QC (PeacoQC/flowAI) + arcsinh/logicle transform
3. **`homogenize_flowframes()`** — Rename markers, reorder to common panel
4. **`downsample_flowframes()`** — Reduce to N cells per sample
5. **`build_combined_expr_matrix()`** + **`run_leiden_with_umap()`** — Combined expression matrix, KNN graph, Leiden communities, UMAP embedding
6. **`clusters_to_gating_csv()`** + **`split_by_cluster()`** + **`run_differential_abundance()`** — Per-sample gating CSVs, tracking matrix, NB-GLMM between conditions
7. **`plot_umap()`** — UMAP plots per sample + per marker

---

## Input

A directory containing `.fcs` files — one per sample. All files in the directory are loaded and processed together. Files may have different marker panels; the harmonization step in Step 3 resolves any mismatches.

```
fcs_samples/
├── sample_001.fcs     # 50,000 cells, 8 markers
├── sample_002.fcs     # 48,000 cells, 8 markers
└── sample_003.fcs     # 52,000 cells, 7 markers  (missing CD56 — handled by harmonization)
```

You also need a **conditions CSV** that maps sample filenames to experimental groups:

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
Rscript pipeline_A_raw_fcs.R --config ../config/pipeline_A_config.R
```

If `--config` is not provided, the pipeline falls back to `config/pipeline_A_config.R`.

### Example Config

```r
# config/pipeline_A_config.R — key parameters

params <- list(
  # --- Paths ---
  fcs_folder          = "/path/to/fcs/files",
  output_dir          = "/path/to/output/pipeline_A_results",
  conditions_csv      = "/path/to/conditions.csv",

  # --- Markers ---
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),
  markername_replacing = list(),   # e.g. list(`PD1 (CD279)` = "PD1")

  # --- Preprocessing ---
  pp_order            = "QT",           # "QT" (QC first) or "TQ" (Transform first)
  qc_type             = "None",         # "PeacoQC", "flowAI", or "None"
  transform_type      = "arcsinh",      # "arcsinh", "logicle", "biexponential", "None"
  offset_list         = rep(0, 6),      # one value per marker (cofactor offset)
  width_list          = rep(5, 6),      # one value per marker (cofactor width)
  truncation          = FALSE,          # truncate values to $PnR range
  save_pp_fcs         = FALSE,          # save preprocessed FCS files
  save_pp_csv         = FALSE,          # save preprocessed CSVs

  # --- Downsampling ---
  downsample_method   = "Consistent",   # "Consistent" (fixed N per sample) or "Maximal"
  downsample_n        = 10000,          # target cell count

  # --- Clustering ---
  n_neighbors         = 15,             # KNN neighbors for UMAP + Leiden
  resolution          = 0.5,            # Leiden resolution (higher = more clusters)

  # --- Analysis ---
  condition_levels    = c("Control", "Treatment"),
  cluster_regex       = "Clustering annotated",
  min_cluster_size    = 10,
  p_adjust_method     = "BH",

  # --- Reproducibility ---
  seed                = 42,

  # --- Plot formats ---
  output_formats      = list(pdf = FALSE, png = TRUE, svg = TRUE),
  plot_width          = 8,
  plot_height         = 6
)
```

### Key Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `fcs_folder` | character | Path to directory containing `.fcs` files |
| `markers_of_interest` | character vector | Marker names to include in analysis |
| `markername_replacing` | named list | Rename map: `list("old_name" = "new_name")` |
| `pp_order` | character | `"QT"` (QC then Transform) or `"TQ"` (reverse) |
| `qc_type` | character | `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `transform_type` | character | `"arcsinh"`, `"logicle"`, `"biexponential"`, or `"None"` |
| `offset_list` | numeric vector | Per-marker offset for arcsinh (length must equal markers) |
| `width_list` | numeric vector | Per-marker width for arcsinh (length must equal markers) |
| `downsample_method` | character | `"Consistent"` (equal N per sample) or `"Maximal"` (proportional) |
| `downsample_n` | integer | Target cell count per sample |
| `n_neighbors` | integer | Nearest neighbors for KNN graph (UMAP + Leiden) |
| `resolution` | numeric | Leiden resolution (higher = more, smaller clusters) |
| `conditions_csv` | character | Path to conditions file mapping samples to groups |
| `condition_levels` | character vector | Two-element vector: `c("Control", "Treatment")` |

---

## Steps in Detail

### Step 1: Load FCS Files

```r
source(file.path(base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  input_dir  = params$fcs_folder,
  truncation = params$truncation
)
```

**What it does:** Reads every `.fcs` file in the directory into a `flowFrame` object (the standard R representation of a single flow cytometry sample). Returns a named list of flowFrames plus sample names extracted from filenames.

**Function signature:**

```r
load_fcs_raw(input_dir, sample_names = NULL, truncation = FALSE,
             min_limit = NULL, extract_csv = FALSE, output_dir = NULL)
```

**Returns:**

```r
list(
  flow_frames  = list(sample_001 = <flowFrame>, sample_002 = <flowFrame>, ...),
  sample_names = c("sample_001", "sample_002", "sample_003"),
  input_files  = c("/path/to/sample_001.fcs", "/path/to/sample_002.fcs", ...)
)
```

**Console output:**

```
--- Step 1: Loading FCS files ---
  Loading: sample_001.fcs ... 50 cells x 6 markers
  Loading: sample_002.fcs ... 50 cells x 6 markers
  Loading: sample_003.fcs ... 50 cells x 6 markers
Loaded 3 samples
```

---

### Step 2: Preprocessing (QC + Transform)

```r
source(file.path(base_dir, "src", "preprocessing", "preprocess_pipeline.R"))
preprocessed <- run_preprocessing(
  flow_frames = flow_data$flow_frames,
  params      = params,
  output_dir  = file.path(params$output_dir, "01_preprocessing")
)
```

**What it does:** Applies quality control to remove aberrant events (margin events, signal anomalies), then transforms fluorescence/mass channels. The processing order is controlled by `pp_order`: `"QT"` runs QC first then transforms, `"TQ"` transforms first then runs QC.

**Function signature:**

```r
run_preprocessing(flow_frames, params, output_dir = NULL)
```

The function reads these keys from the `params` list:

| Key | Type | Description |
|-----|------|-------------|
| `pp_order` | character | `"QT"` or `"TQ"` — order of QC and transformation |
| `qc_type` | character | `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `markers_of_interest` | character vector | Which markers to process |
| `offset_list` | numeric vector | Per-marker offsets for transformation |
| `width_list` | numeric vector | Per-marker widths for transformation |
| `transform_type` | character | Transformation type |
| `save_pp_fcs` | logical | Save preprocessed FCS files to disk |
| `save_pp_csv` | logical | Save preprocessed expression matrices as CSVs |

**Returns:** Named list of preprocessed flowFrame objects.

**Console output:**

```
--- Step 2: Preprocessing ---
  Processing order: QT (QC then Transform)
  QC: None (skipping quality control)
  Transforming sample_001 with arcsinh (6 markers)...
  Transforming sample_002 with arcsinh (6 markers)...
  Transforming sample_003 with arcsinh (6 markers)...
```

---

### Step 3: Panel Harmonization

```r
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames         = preprocessed,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  markername_replacing = params$markername_replacing,
  output_dir          = file.path(params$output_dir, "02_harmonization")
)
```

**What it does:** Ensures all flowFrames have the same marker columns in the same order. If different FCS files came from different instrument configurations, some markers may be missing or named differently. This step:

1. Applies `markername_replacing` to rename markers (e.g. `PD1 (CD279)` → `PD1`)
2. Subsets each flowFrame to only the `markers_of_interest`
3. Reorders columns to match across all samples
4. Fills missing markers with `NA`

**Function signature:**

```r
homogenize_flowframes(flow_frames, sample_names, markers_of_interest,
                      markername_replacing = list(),
                      save_csv = FALSE, save_fcs = FALSE, output_dir = NULL)
```

**Returns:** Named list of homogenized flowFrame objects.

**Console output:**

```
--- Step 3: Panel Homogenization ---
  Homogenizing 3 samples to 6 markers
    Homogenizing sample_001 ...  50 cells x 6 markers
    Homogenizing sample_002 ...  50 cells x 6 markers
    Homogenizing sample_003 ...  50 cells x 6 markers
```

---

### Step 4: Downsampling

```r
source(file.path(base_dir, "src", "preprocessing", "downsampling.R"))
downsampled <- downsample_flowframes(
  flow_frames  = homogenized,
  sample_names = flow_data$sample_names,
  target_cells = params$downsample_n,
  method       = params$downsample_method,
  seed         = params$seed
)
```

**What it does:** Reduces cell counts to balance representation across samples and manage memory. Two methods are available:

| Method | Behavior |
|--------|----------|
| `"Consistent"` | Takes exactly `target_cells / n_samples` cells from each sample |
| `"Maximal"` | Proportional to original counts; total = `target_cells` |

**Function signature:**

```r
downsample_flowframes(flow_frames, sample_names, target_cells = 10000,
                      method = "Maximal", seed = 42,
                      save_fcs = FALSE, output_dir = NULL)
```

**Returns:** Named list of downsampled flowFrame objects.

**Console output:**

```
--- Step 4: Downsampling ---
  Downsampling 3 samples (total: 150 cells)
  Method: Consistent | Target: 100 cells
  Consistent: taking 33 cells per sample
    sample_001 : 50 -> 33
    sample_002 : 50 -> 33
    sample_003 : 50 -> 33
  Total after downsampling: 99 cells
```

---

### Step 5: Leiden Clustering + UMAP

```r
source(file.path(base_dir, "src", "clustering", "leiden.R"))

# Build combined expression matrix from all samples
expr_data <- build_combined_expr_matrix(
  downsampled, flow_data$sample_names, params$markers_of_interest
)
clusters <- run_leiden_with_umap(
  expr_matrix   = expr_data$combined_matrix,
  sample_labels = expr_data$sample_labels,
  output_dir    = file.path(params$output_dir, "03_clustering"),
  n_neighbors   = params$n_neighbors,
  resolution    = params$resolution,
  seed          = params$seed
)
```

**What it does:**

1. `build_combined_expr_matrix()` extracts expression matrices from all downsampled flowFrames and combines them into a single cells x markers matrix with matching sample labels
2. Builds a K-nearest-neighbor graph (using `uwot`)
3. Converts to a sparse adjacency matrix
4. Runs Leiden community detection to assign clusters
5. Computes a 2D UMAP embedding for visualization

**`build_combined_expr_matrix()` signature:**

```r
build_combined_expr_matrix(flow_frames, sample_names, markers_of_interest)
```

Returns a list with `combined_matrix` (cells x markers matrix) and `sample_labels` (character vector identifying the source sample for each row).

The `resolution` parameter controls cluster granularity: lower values produce fewer, larger clusters; higher values produce more, smaller clusters. Typical range: 0.1–2.0.

**Function signature:**

```r
run_leiden_with_umap(expr_matrix, sample_labels = NULL, output_dir = NULL,
                     n_neighbors = 30, resolution = 1.0, seed = 42)
```

**Returns:** A `data.frame` with columns `UMAP1`, `UMAP2`, `cluster`, `sample`.

**Output file** — `leiden_clusters.csv`:

```
"UMAP1","UMAP2","cluster","sample"
1.84,0.85,1,"sample_001"
-0.85,-1.51,1,"sample_001"
0.96,-0.36,1,"sample_001"
3.37,-0.17,1,"sample_003"
3.12,-0.20,1,"sample_003"
...
```

Each row is one cell. `cluster` is an integer (1, 2, 3...) assigned by Leiden. `sample` identifies which input file the cell came from.

**Console output:**

```
--- Step 5: Clustering ---
  Computing UMAP (n_neighbors=15)...
  Building KNN graph (k=15)...
  Running Leiden clustering (resolution=0.5)...
  Found 1 cluster
  Cluster sizes: 99
  Results saved to: output/pipeline_A_results/03_clustering/leiden_clusters.csv
```

---

### Step 6: Analysis (Cluster Splitting + Differential Abundance)

This step has two parts. First, the Leiden cluster assignments are converted to per-sample CSV files with binary gating-style columns, then those CSVs are split by cluster and tested for differential abundance.

**Part A: Save cluster data as per-sample CSVs**

```r
cluster_data_dir <- file.path(params$output_dir, "04_analysis", "cluster_data")
clusters_to_gating_csv(
  combined_matrix = expr_data$combined_matrix,
  sample_labels   = expr_data$sample_labels,
  cluster_ids     = clusters$cluster,
  output_dir      = cluster_data_dir
)
```

Each sample gets a CSV with expression columns plus binary cluster membership columns (1 = cell belongs to this cluster, 0 = does not).

**`clusters_to_gating_csv()` signature:**

```r
clusters_to_gating_csv(combined_matrix, sample_labels, cluster_ids, output_dir)
```

Splits the combined matrix back into per-sample data.frames, appends binary gating columns for each cluster (e.g. `Cluster_1 (Clustering annotated)`), and writes one semicolon-separated CSV per sample.

**Part B: Split by cluster**

```r
source(file.path(base_dir, "src", "analysis", "cluster_splitting.R"))
split_data <- split_by_cluster(
  input_dir     = cluster_data_dir,
  output_dir    = file.path(params$output_dir, "04_analysis", "split_clusters"),
  cluster_regex = params$cluster_regex
)
```

**Function signature:**

```r
split_by_cluster(input_dir, output_dir,
                 cluster_regex = "Clustering annotated",
                 markers_of_interest = NULL,
                 min_cells_per_cluster = 10)
```

**Output** — `cluster_tracking_matrix.csv`:

```
file;Cluster_1;Cluster_2;Cluster_3
sample_001;33;12;5
sample_002;28;18;4
sample_003;35;8;7
```

This tracking matrix shows how many cells from each sample ended up in each cluster.

**Part C: Differential abundance**

```r
source(file.path(base_dir, "src", "analysis", "differential_abundance.R"))
da_results <- run_differential_abundance(
  data_dir         = cluster_data_dir,
  conditions_csv   = params$conditions_csv,
  output_dir       = file.path(params$output_dir, "04_analysis", "differential_abundance"),
  cluster_regex    = params$cluster_regex,
  condition_levels = params$condition_levels,
  min_cluster_size = params$min_cluster_size,
  p_adjust_method  = params$p_adjust_method,
  seed             = params$seed
)
```

**What it does:** Fits a negative binomial GLMM (via `lme4`) for each cluster, modeling cell count as a function of condition with sample as a random effect:

```
count ~ condition + offset(log_total) + (1 | sample_id)
```

**Function signature:**

```r
run_differential_abundance(data_dir, conditions_csv, output_dir,
                           cluster_regex = "Clustering annotated",
                           condition_levels = c("pCR", "npCR"),
                           markers_of_interest = NULL,
                           min_cluster_size = 10,
                           p_adjust_method = "BH",
                           max_cells_per_file = 100000,
                           seed = 123)
```

**Output** — `glmm_results.csv`:

```
cluster;estimate;std_error;z_value;p_value;log2fc;odds_ratio;method;p_adjusted;significant
Cluster_1;-0.23;0.15;-1.53;0.126;-0.33;0.79;glmerMod;0.126;FALSE
Cluster_2;0.87;0.21;4.14;0.00003;1.26;2.39;glmerMod;0.0001;TRUE
```

Column meanings:

| Column | Meaning |
|--------|---------|
| `estimate` | Model coefficient for condition effect |
| `std_error` | Standard error of the estimate |
| `z_value` | Wald z-statistic |
| `p_value` | Raw p-value |
| `log2fc` | Log2 fold change between conditions |
| `odds_ratio` | Exponentiated estimate |
| `p_adjusted` | Benjamini-Hochberg adjusted p-value |
| `significant` | `TRUE` if `p_adjusted < 0.05` |

Also produces `volcano_plot.png`: a scatter plot of log2FC vs -log10(adjusted p-value), with significant clusters highlighted.

---

### Step 7: Visualization

```r
source(file.path(base_dir, "src", "visualization", "umap.R"))
plot_umap(
  flow_frames         = downsampled,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  output_dir          = file.path(params$output_dir, "05_visualization"),
  n_neighbors         = params$n_neighbors,
  output_formats      = params$output_formats,
  seed                = params$seed
)
```

**What it does:** Computes a UMAP embedding from the downsampled flowFrames and generates scatter plots:

- `umap_sample.png` — colored by sample origin
- `umap_CD3.png`, `umap_CD4.png`, ... — one plot per marker, colored by expression intensity (viridis palette)

**Function signature:**

```r
plot_umap(flow_frames, sample_names, markers_of_interest = NULL,
          output_dir, n_neighbors = 15, min_dist = 0.1,
          marker_selection_mode = "all",
          output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
          plot_width = 8, plot_height = 6, seed = 42)
```

**Returns:**

```r
list(
  umap_coords   = matrix(),    # N x 2 UMAP coordinates
  sample_labels = character(),  # Sample label per cell
  plots         = list()        # Named list of ggplot objects
)
```

---

## Output Directory Structure

After a successful run with 3 samples and 6 markers:

```
output/pipeline_A_results/
├── 01_preprocessing/                    # (created if save_pp_fcs/csv enabled)
├── 02_harmonization/                    # (created if save_csv/fcs enabled)
├── 03_clustering/
│   └── leiden_clusters.csv              # UMAP1, UMAP2, cluster, sample
├── 04_analysis/
│   ├── cluster_data/                    # Per-sample CSVs with binary cluster columns
│   │   ├── sample_001.csv
│   │   ├── sample_002.csv
│   │   └── sample_003.csv
│   ├── split_clusters/
│   │   ├── cluster_tracking_matrix.csv  # Sample-by-cluster cell counts
│   │   ├── sample_001_Cluster_1.csv     # Expression data for cluster 1, sample 1
│   │   └── ...
│   └── differential_abundance/
│       ├── glmm_results.csv             # Per-cluster GLMM test results
│       └── volcano_plot.png             # log2FC vs -log10(p-adjusted)
└── 05_visualization/
    ├── umap_sample.png                  # UMAP colored by sample
    ├── umap_CD3.png                     # UMAP colored by CD3 expression
    ├── umap_CD4.png
    ├── umap_CD8.png
    ├── umap_CD19.png
    ├── umap_CD56.png
    └── umap_HLA-DR.png
```

---

## Pipeline vs Pipeline E: Key Differences

| Aspect | Pipeline A | Pipeline E |
|--------|-----------|-----------|
| `truncation` | Configurable (often `TRUE`) | Always `FALSE` |
| QC step | Yes (`PeacoQC` / `flowAI` / `None`) | No |
| Transform step | Yes (`arcsinh` / `logicle` / etc.) | No |
| Analysis step | Yes (differential abundance) | No |
| Focus | Full end-to-end analysis | Quick unsupervised discovery |

Pipeline E is Pipeline A with Steps 2 and 6 removed. Use Pipeline A when you need the complete workflow; use Pipeline E when your data is already preprocessed.

---

## Tips

- **Large datasets (>10M events):** Use `downsample_method = "Consistent"` with a moderate `downsample_n` (10,000–50,000) to keep memory usage manageable.
- **Finding the right resolution:** Start with `resolution = 0.5` for broad clusters, increase to 1.0–2.0 for finer granularity. Run Pipeline E multiple times to iterate quickly.
- **Panel mismatch:** If marker names differ between FCS files (e.g. `PD1 (CD279)` vs `PD1`), use `markername_replacing` to standardize.
- **Transformation cofactors:** For mass cytometry, use `width_list = rep(5, n_markers)`. For fluorescence, try `logicle` transformation instead.
- The pipeline logs elapsed time per step and prints a total runtime at the end.
