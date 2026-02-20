[< Back to Documentation Index](../README.md)

# Pipeline E: Preprocessed FCS — Quick Clustering & UMAP

Pipeline E is a streamlined version of [Pipeline A](pipeline-A-raw-fcs.md) for FCS files that have already been quality-controlled and transformed. It skips all preprocessing steps and goes directly from loading to panel harmonization, downsampling, Leiden clustering, and UMAP visualization. Use it when your FCS files come from an upstream preprocessing tool and you want to explore cluster structure without re-running QC.

---

## When to Use Pipeline E

Use this pipeline when:

- Your FCS files are **already preprocessed** (QC'd, transformed, cleaned)
- You want to **quickly explore clustering** with different parameters
- You only need **UMAP plots** and **cluster assignments** — no differential abundance
- You want to **iterate on `resolution` and `n_neighbors`** without the overhead of full Pipeline A

If your FCS files are raw and unprocessed, use [Pipeline A](pipeline-A-raw-fcs.md) instead. If you need differential abundance testing, also use Pipeline A.

---

## Pipeline Steps

1. **`load_fcs_raw()`** — Load FCS files (truncation=FALSE, no QC applied)
2. **`homogenize_flowframes()`** — Unify marker panels across files
3. **`downsample_flowframes()`** — Reduce cell count for tractability
4. **`run_leiden_with_umap()`** — KNN graph, Leiden clusters, UMAP embedding
5. **`plot_umap()`** — UMAP per sample + per marker

---

## Input

A directory containing preprocessed `.fcs` files. These files should have already undergone:

- Quality control (removal of margin events, doublets, dead cells)
- Transformation (arcsinh, logicle, or similar)

```
preprocessed_fcs/
├── sample_001.fcs     # Already QC'd and transformed
├── sample_002.fcs
└── sample_003.fcs
```

No conditions CSV is needed — Pipeline E does not perform differential testing.

---

## Configuration

### Running the Pipeline

```bash
cd pipelines
Rscript pipeline_E_preprocessed_umap.R --config ../config/pipeline_E_config.R
```

### Example Config

```r
# config/pipeline_E_config.R

params <- list(
  # --- Paths ---
  fcs_folder            = "/path/to/preprocessed/fcs",
  output_dir            = "/path/to/output/pipeline_E_results",

  # --- Markers ---
  markers_of_interest   = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),
  markername_replacing  = list(),     # e.g. list(`PD1 (CD279)` = "PD1")

  # --- Downsampling ---
  downsample_method     = "Consistent",
  downsample_n          = 10000,

  # --- Clustering ---
  n_neighbors           = 15,
  resolution            = 0.5,

  # --- Reproducibility ---
  seed                  = 42,

  # --- Plot formats ---
  output_formats        = list(pdf = FALSE, png = TRUE, svg = TRUE),
  plot_width            = 8,
  plot_height           = 6
)
```

### Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fcs_folder` | character | *required* | Directory with preprocessed FCS files |
| `markers_of_interest` | character vector | *required* | Markers to include in clustering and UMAP |
| `markername_replacing` | named list | `list()` | Rename map for inconsistent marker names |
| `downsample_method` | character | `"Maximal"` | `"Consistent"` or `"Maximal"` |
| `downsample_n` | integer | `10000` | Target cell count |
| `n_neighbors` | integer | `15` | KNN neighbors for UMAP and Leiden |
| `resolution` | numeric | `0.5` | Leiden resolution (higher = more clusters) |
| `seed` | integer | `42` | Random seed for reproducibility |

---

## Steps in Detail

### Step 1: Load FCS Files (No QC)

```r
source(file.path(base_dir, "src", "loading", "load_fcs_raw.R"))
flow_data <- load_fcs_raw(
  input_dir  = params$fcs_folder,
  truncation = FALSE
)
```

**What it does:** Loads all `.fcs` files without applying any quality control or transformation. The `truncation = FALSE` flag is important: since the data is already preprocessed, truncating values to the original instrument `$PnR` range would discard valid transformed values.

**Function signature:**

```r
load_fcs_raw(input_dir, sample_names = NULL, truncation = FALSE,
             min_limit = NULL, extract_csv = FALSE, output_dir = NULL)
```

**Returns:**

```r
list(
  flow_frames  = list(sample_001 = <flowFrame>, ...),
  sample_names = c("sample_001", "sample_002", "sample_003"),
  input_files  = c("/path/sample_001.fcs", ...)
)
```

**Console output:**

```
--- Step 1: Loading preprocessed FCS files ---
  Loading: sample_001.fcs ... 50 cells x 6 markers
  Loading: sample_002.fcs ... 50 cells x 6 markers
  Loading: sample_003.fcs ... 50 cells x 6 markers
Loaded 3 samples
```

---

### Step 2: Panel Harmonization

```r
source(file.path(base_dir, "src", "harmonization", "panel_homogenization.R"))
homogenized <- homogenize_flowframes(
  flow_frames         = flow_data$flow_frames,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  markername_replacing = params$markername_replacing,
  output_dir          = file.path(params$output_dir, "01_harmonized")
)
```

**What it does:** Even preprocessed files may have different marker panels (different instruments, different staining panels). This step unifies all flowFrames to a common set of markers in a consistent column order. Missing markers are filled with `NA`; extra markers are dropped.

**Function signature:**

```r
homogenize_flowframes(flow_frames, sample_names, markers_of_interest,
                      markername_replacing = list(),
                      save_csv = FALSE, save_fcs = FALSE, output_dir = NULL)
```

**Console output:**

```
--- Step 2: Panel Homogenization ---
  Homogenizing 3 samples to 6 markers
    Homogenizing sample_001 ...  50 cells x 6 markers
    Homogenizing sample_002 ...  50 cells x 6 markers
    Homogenizing sample_003 ...  50 cells x 6 markers
```

---

### Step 3: Downsampling

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

**What it does:** Reduces cell counts for computational tractability. Particularly important when preprocessed files retain large event counts from the original acquisition.

**Function signature:**

```r
downsample_flowframes(flow_frames, sample_names, target_cells = 10000,
                      method = "Maximal", seed = 42,
                      save_fcs = FALSE, output_dir = NULL)
```

| Method | Behavior |
|--------|----------|
| `"Consistent"` | Takes exactly `target_cells / n_samples` cells from each sample |
| `"Maximal"` | Proportional to original counts; total = `target_cells` |

---

### Step 4: Leiden Clustering + UMAP

```r
source(file.path(base_dir, "src", "clustering", "leiden.R"))

# Build combined expression matrix from all samples
expr_list <- extract_exprmtx_list(downsampled, flow_data$sample_names, params$markers_of_interest)
combined_matrix <- do.call(rbind, expr_list)
sample_labels <- rep(names(expr_list), sapply(expr_list, nrow))

clusters <- run_leiden_with_umap(
  expr_matrix   = combined_matrix,
  sample_labels = sample_labels,
  output_dir    = file.path(params$output_dir, "02_clustering"),
  n_neighbors   = params$n_neighbors,
  resolution    = params$resolution,
  seed          = params$seed
)
```

**What it does:**

1. Extracts expression matrices from all downsampled flowFrames
2. Combines into a single cells x markers matrix
3. Builds a KNN graph and runs Leiden community detection
4. Computes UMAP embedding for visualization

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
...
```

---

### Step 5: UMAP Visualization

```r
source(file.path(base_dir, "src", "visualization", "umap.R"))
plot_umap(
  flow_frames         = downsampled,
  sample_names        = flow_data$sample_names,
  markers_of_interest = params$markers_of_interest,
  output_dir          = file.path(params$output_dir, "03_visualization"),
  n_neighbors         = params$n_neighbors,
  output_formats      = params$output_formats,
  seed                = params$seed
)
```

**What it does:** Computes UMAP from the downsampled flowFrames and generates:

- `umap_sample.png` — UMAP colored by sample identity
- `umap_CD3.png`, `umap_CD4.png`, ... — one plot per marker, colored by expression (viridis palette)

**Function signature:**

```r
plot_umap(flow_frames, sample_names, markers_of_interest = NULL,
          output_dir, n_neighbors = 15, min_dist = 0.1,
          output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
          plot_width = 8, plot_height = 6, seed = 42)
```

**Output:**

```
03_visualization/
├── umap_sample.png     # Colored by sample
├── umap_CD3.png        # Colored by CD3 expression
├── umap_CD4.png        # Colored by CD4 expression
├── umap_CD8.png
├── umap_CD19.png
├── umap_CD56.png
└── umap_HLA-DR.png
```

---

## Output Directory Structure

```
output/pipeline_E_results/
├── 01_harmonized/           # (created if save_csv/fcs enabled in harmonization)
├── 02_clustering/
│   └── leiden_clusters.csv  # UMAP1, UMAP2, cluster, sample
└── 03_visualization/
    ├── umap_sample.png      # UMAP colored by sample
    ├── umap_CD3.png         # UMAP colored by CD3 expression
    ├── umap_CD4.png
    ├── umap_CD8.png
    ├── umap_CD19.png
    ├── umap_CD56.png
    └── umap_HLA-DR.png
```

---

## Pipeline A vs Pipeline E

| Aspect | Pipeline A | Pipeline E |
|--------|-----------|-----------|
| Input | Raw, unprocessed FCS | Already QC'd and transformed FCS |
| `truncation` | Configurable (`TRUE`/`FALSE`) | Always `FALSE` |
| QC step | Yes (PeacoQC, flowAI, or None) | **Skipped** |
| Transform step | Yes (arcsinh, logicle, etc.) | **Skipped** |
| Differential abundance | Yes (NB-GLMM) | **No** |
| Output focus | Full analysis + statistics | Cluster discovery + UMAP |
| Typical runtime | Longer (7 steps) | Shorter (5 steps) |

Pipeline E is Pipeline A with Steps 2 (preprocessing) and 6 (statistical analysis) removed. Use Pipeline E for rapid iteration on clustering parameters.

---

## Tips

- **Iterating on parameters:** Since Pipeline E is lean (5 steps vs 7 for Pipeline A), it runs much faster. Use it to try different `resolution` and `n_neighbors` values before committing to a full Pipeline A run.
- **No analysis output:** Pipeline E does not run differential abundance, enrichment, or any statistical comparison between conditions. If you need that, use Pipeline A or Pipeline C.
- **truncation = FALSE is intentional.** Preprocessed FCS files often have transformed values that exceed the original instrument `$PnR` range. Truncating would clip valid data. Pipeline E always loads with `truncation = FALSE`.
- **Quality depends on upstream preprocessing.** Since no QC is applied, verify that dead cells, doublets, and debris were removed before the data was exported.
