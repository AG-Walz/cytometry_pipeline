# Composability Guide

How to mix and match scripts from different pipeline stages to build custom
analysis workflows. Every script in this project is designed as a standalone
function that accepts well-defined inputs and produces well-defined outputs.
You are not limited to the five predefined pipelines (A through E); any
script can connect to any other script as long as the data formats match.

This guide covers:

1. [Data Flow Diagram](#1-data-flow-diagram)
2. [Format Contracts](#2-format-contracts)
3. [Bridge Functions](#3-bridge-functions)
4. [Example Workflows](#4-example-workflows)
5. [Decision Tree](#5-decision-tree)
6. [Quick Reference Table](#6-quick-reference-table)
7. [Method Swapping](#7-method-swapping)

---

## 1. Data Flow Diagram

The diagram below shows every script group, the data formats they produce
and consume, and the bridge functions that convert between formats.

```
INPUT FILES
============

  .fcs files          .csv expression       .csv gating + .fcs      .wsp + .fcs
       |                    |                       |                     |
       v                    v                       v                     v
  load_fcs_raw()     load_csv_expression()   load_csv_gating()    load_flowjo_wsp()
       |                    |                       |               load_multiple_workspaces()
       v                    v                       v                     |
  [flowFrame list]   [data.table list]       [gating CSVs on disk]       v
       |                    |                       |              [GatingSet / pop strings]
       |                    |                       |                     |
       |                    |                       |                     |
  PREPROCESSING             |                       |               METADATA QC
  =============             |                       |               ==========
       |                    |                       |                     |
       v                    |                       |                     v
  run_preprocessing()       |                       |             compare_gatings()
  downsample_flowframes()   |                       |             compute_fleiss_kappa()
  transform_flowframe()     |                       |                     |
  quality_control()         |                       |                     v
       |                    |                       |             [Venn/Euler/UpSet]
       v                    |                       |               plot_euler()
  [flowFrame list]          |                       |               plot_upset()
       |                    |                       |
       |                    |                       |
  HARMONIZATION             |                       |
  =============             |                       |
       |                    |                       |
       v                    |                       |
  homogenize_flowframes()   |                       |
  run_batch_correction()    |                       |
       |                    |                       |
       v                    |                       |
  [flowFrame list]          |                       |
       |                    |                       |
       |                    |                       |
       |   BRIDGE FUNCTIONS |                       |
       |   =====================================    |
       |         |                    |              |
       v         v                    v              |
  build_combined_expr_matrix()   create_flowframes_from_matrices()
  (flowFrames --> matrix)        (matrix/data.table --> flowFrames)
       |                              |              |
       v                              v              |
  [numeric matrix +            [flowFrame list]      |
   sample_labels]                     |              |
       |                              |              |
       |  .............................|              |
       |  :                                          |
       v  v                                          |
  CLUSTERING (via run_clustering() dispatcher)        |
  ===============================================    |
       |                                             |
       v                                             |
  run_leiden()        run_flowsom()                   |
  run_phenograph()    run_dbscan()                   |
       |                                             |
       v                                             |
  [cluster labels]                                   |
       |                                             |
       v                                             |
  clusters_to_gating_csv() ----> [gating CSVs] <-----+
       |                              |
       |                              |
       +------------------------------+
                      |
                      v
  ANALYSIS (all operate on gating CSVs or matrices)
  =================================================
       |
       +---> split_by_cluster()         --> [per-cluster CSVs] + [tracking matrix]
       |          |
       |          +--> tracking_matrix_to_long()  --> [long data.frame]
       |          +--> create_expression_matrices()--> [expression matrices]
       |
       +---> build_abundance_matrix()   --> [abundance data.frame]
       +---> build_tracking_matrix()    --> [tracking matrix]
       |
       +---> run_differential_abundance() --> [GLMM results] + volcano plot
       +---> run_cluster_enrichment()     --> [enrichment results]
       +---> run_marker_statistics()      --> [per-marker stats]
       +---> calculate_cell_ratios()      --> [ratio data.frame]
       |
       v
  VISUALIZATION (accept flowFrames, matrices, or data.frames)
  ===========================================================

  flowFrame input:                Matrix input:
    plot_umap()                     plot_umap_from_matrix()
    plot_density()                  plot_density_from_matrix()
                                    plot_pca()
                                    plot_heatmap()
                                    plot_cluster_heatmap()

  Data.frame input:
    plot_composition_barplot()    (long data.frame)
    plot_abundance()              (tracking matrix)
    plot_volcano()                (results data.frame)
    plot_upset() / plot_euler()   (binary matrix or set lists)
    plot_venn()                   (set lists)
```

### Key connection paths at a glance

| From | To | Via |
|------|----|-----|
| flowFrame list | numeric matrix | `build_combined_expr_matrix()` |
| data.table list | flowFrame list | `create_flowframes_from_matrices()` |
| data.table list | numeric matrix | `combine_expression_matrices()` |
| cluster labels + matrix | gating CSVs on disk | `clusters_to_gating_csv()` |
| expression matrix | supercell matrix | `run_supercell()` |
| expression matrix | cluster labels (FlowSOM) | `run_flowsom()` |
| expression matrix | cluster labels (PhenoGraph) | `run_phenograph()` |
| expression matrix | cluster labels (DBSCAN) | `run_dbscan()` |
| clusters + expression | cell type annotations | `annotate_clusters()` |
| clusters + expression | validation metrics | `validate_clusters()` |
| expression matrix | batch-corrected matrix | `run_combat()` |
| flowFrame list + batch labels | corrected flowFrames | `run_cytonorm()` |
| gating CSVs | tracking matrix | `split_by_cluster()` or `build_tracking_matrix()` |
| tracking matrix (wide) | long data.frame | `tracking_matrix_to_long()` |
| gating CSVs | expression matrices | `split_by_cluster()` then `create_expression_matrices()` |

---

## 2. Format Contracts

Every script in the pipeline reads and writes one of the five core formats
below. When two scripts share a format, they can be chained directly.

### 2.1 flowFrame list

**What it is.** A named R list where each element is a `flowCore::flowFrame`
object. Names are sample identifiers.

**Structure:**

```r
flow_frames <- list(
  "Sample_01" = <flowFrame: 50000 cells x 35 channels>,
  "Sample_02" = <flowFrame: 48000 cells x 35 channels>
)
```

**Required properties:**

- Each flowFrame contains an expression matrix accessible via `flowCore::exprs(ff)`.
- Marker names are embedded in the flowFrame parameter annotations (accessed
  via `flowCore::markernames(ff)` or `explore_data(ff)$markers`).
- All flowFrames in the list should share the same marker set after
  harmonization. Before harmonization, panels may differ.

**Produced by:**
`load_fcs_raw()`, `run_preprocessing()`, `downsample_flowframes()`,
`homogenize_flowframes()`, `create_flowframes_from_matrices()`

**Consumed by:**
`run_preprocessing()`, `downsample_flowframes()`,
`homogenize_flowframes()`, `run_batch_correction()`,
`build_combined_expr_matrix()`, `plot_umap()`, `plot_density()`

### 2.2 Numeric matrix

**What it is.** A standard R numeric matrix with rows as cells and columns as
markers. Often accompanied by a parallel `sample_labels` character vector of
the same length as `nrow(matrix)`.

**Structure:**

```r
#          CD3    CD4    CD8   CD19   CD56
# [1,]   3.21   2.45   0.12   0.05   1.87
# [2,]   3.15   0.08   2.91   0.03   0.22
# ...
# dim: 98000 x 5
# colnames: c("CD3", "CD4", "CD8", "CD19", "CD56")

sample_labels <- c("Sample_01", "Sample_01", ..., "Sample_02", ...)
# length: 98000 (one entry per row)
```

**Required properties:**

- Column names must be marker names (character strings).
- Values must be numeric (no NA unless explicitly handled).
- No metadata columns (sample, cluster, etc.) -- those travel as separate
  vectors or in the return list.

**Produced by:**
`build_combined_expr_matrix()` (returns `$combined_matrix` + `$sample_labels`),
`combine_expression_matrices()`, `run_batch_correction()` (returns `$corrected_matrix`)

**Consumed by:**
`run_leiden()`, `run_leiden_with_umap()`, `clusters_to_gating_csv()`,
`plot_umap_from_matrix()`, `plot_density_from_matrix()`, `plot_pca()`,
`plot_heatmap()`, `harmony_normalize()`

### 2.3 Gating CSV (per-sample CSV with binary cluster columns)

**What it is.** One CSV file per sample on disk. Each row is a cell. Columns
are expression markers plus binary (0/1) cluster indicator columns.

**Structure (semicolon-separated):**

```
CD3;CD4;CD8;CD19;CD56;Cluster_1 (Clustering annotated);Cluster_2 (Clustering annotated)
3.21;2.45;0.12;0.05;1.87;1;0
3.15;0.08;2.91;0.03;0.22;0;1
```

**Required properties:**

- Separator: semicolon (`;`) by default, auto-detected by `detect_csv_sep()`.
- Expression columns: numeric, named by marker.
- Cluster columns: binary integer (`0` or `1`). Column names match the
  regex `"Clustering annotated"` (configurable via `cluster_regex`).
- One file per sample, named `<sample_name>.csv`.

**Produced by:**
`clusters_to_gating_csv()`, `load_csv_gating()`

**Consumed by:**
`split_by_cluster()`, `build_abundance_matrix()`, `build_tracking_matrix()`,
`run_differential_abundance()`, `run_cluster_enrichment()`

### 2.4 Tracking matrix (wide format)

**What it is.** A data.frame with one row per sample and one column per
cluster. Cell values are integer cell counts.

**Structure:**

```r
#   file         Cluster_1  Cluster_2  Cluster_3
#   Sample_01        1250        870        430
#   Sample_02        1100        920        510
```

**Required properties:**

- First column is named `file` (sample identifier).
- Remaining columns are cluster names (character).
- Values are non-negative integers.

**Produced by:**
`split_by_cluster()` (returns `$tracking_matrix`), `build_tracking_matrix()`

**Consumed by:**
`tracking_matrix_to_long()`, `plot_abundance()`, `plot_abundance_stacked()`,
`plot_abundance_proportional()`, `plot_abundance_heatmap()`,
`plot_abundance_totals()`, `calculate_cell_ratios()` (via CSV on disk)

### 2.5 Long data.frame

**What it is.** A tidy data.frame with one row per sample-cluster combination.
Used for barplots and ggplot-based visualizations.

**Structure:**

```r
#   file         cluster      count
#   Sample_01    Cluster_1     1250
#   Sample_01    Cluster_2      870
#   Sample_02    Cluster_1     1100
#   Sample_02    Cluster_2      920
```

**Required properties:**

- Column names are configurable via `id_col`, `names_to`, `values_to` in
  `tracking_matrix_to_long()`, but default to `file`, `cluster`, `count`.

**Produced by:**
`tracking_matrix_to_long()`

**Consumed by:**
`plot_composition_barplot()`

---

## 3. Bridge Functions

Bridge functions are the key to cross-pipeline composability. They convert
data between formats so that scripts designed for one input type can accept
data from a different source.

### 3.1 `build_combined_expr_matrix()` -- flowFrame list to matrix

**Location:** `src/R_shared/flowframe_utils.R`

**Purpose:** Concatenates expression data from multiple flowFrames into a
single numeric matrix with a parallel sample-labels vector. This is the
standard way to move from FCS-originated data into matrix-based scripts
(clustering, PCA, matrix-based UMAP).

```r
result <- build_combined_expr_matrix(flow_frames, sample_names, markers_of_interest)
# result$combined_matrix  -- numeric matrix (cells x markers)
# result$sample_labels    -- character vector (length = nrow)
```

### 3.2 `create_flowframes_from_matrices()` -- matrix/data.table to flowFrame list

**Location:** `src/loading/create_flowframes.R`

**Purpose:** Wraps expression matrices (or data.tables) into flowFrame
objects. This lets CSV-loaded data flow into flowFrame-based scripts like
`plot_umap()`, `plot_density()`, `downsample_flowframes()`, and
`homogenize_flowframes()`.

```r
flow_frames <- create_flowframes_from_matrices(expr_list, markers_of_interest)
# expr_list: named list of matrices or data.frames
# Returns: named list of flowFrame objects
```

### 3.3 `combine_expression_matrices()` -- data.table list to matrix

**Location:** `src/analysis/expression_matrices.R`

**Purpose:** Row-binds multiple data.tables into a single numeric matrix,
subsetting to shared markers. A lighter-weight alternative to
`build_combined_expr_matrix()` when you already have data.tables (not
flowFrames).

```r
matrix <- combine_expression_matrices(expr_data_list, markers, max_cells = 50000)
```

### 3.4 `clusters_to_gating_csv()` -- cluster labels to per-sample gating CSVs

**Location:** `src/clustering/leiden.R`

**Purpose:** Takes a combined expression matrix, sample labels, and cluster
assignments, then writes one CSV per sample with expression columns plus
binary cluster indicator columns. This converts in-memory clustering results
into the on-disk gating CSV format consumed by all analysis scripts.

```r
clusters_to_gating_csv(combined_matrix, sample_labels, cluster_ids, output_dir)
# Writes: output_dir/Sample_01.csv, output_dir/Sample_02.csv, ...
```

### 3.5 `plot_umap_from_matrix()` / `plot_density_from_matrix()` -- direct matrix visualization

**Location:** `src/visualization/umap.R`, `src/visualization/density_plots.R`

**Purpose:** Accept a numeric matrix directly, bypassing any need for
flowFrame conversion. If you have matrix data and only need a plot, use
these instead of converting to flowFrame first.

```r
plot_umap_from_matrix(expr_matrix, sample_labels, output_dir = "out/")
plot_density_from_matrix(expr_matrix, sample_labels, markers, output_dir = "out/")
```

---

## 4. Example Workflows

Each workflow below shows the exact function call chain with the data
format at each step. Copy these as templates and modify parameters for
your study.

### 4.1 FCS to preprocess to cluster to abundance analysis to volcano plot

The full Pipeline A path: raw FCS files through preprocessing,
harmonization, Leiden clustering, and differential abundance with
volcano plot output.

```r
# Step 1: Load raw FCS files
#   Input:  .fcs files on disk
#   Output: flowFrame list
loaded <- load_fcs_raw(input_dir = "data/fcs/")
flow_frames <- loaded$flow_frames
sample_names <- loaded$sample_names

# Step 2: Preprocess (QC + transform)
#   Input:  flowFrame list
#   Output: flowFrame list
pp_frames <- run_preprocessing(flow_frames, params, output_dir = "out/preprocessing/")

# Step 3: Downsample
#   Input:  flowFrame list
#   Output: flowFrame list
ds_frames <- downsample_flowframes(pp_frames, sample_names, target_cells = 50000)

# Step 4: Harmonize panels
#   Input:  flowFrame list
#   Output: flowFrame list (unified markers)
hg_frames <- homogenize_flowframes(ds_frames, sample_names, markers_of_interest)

# Step 5: BRIDGE -- flowFrames to matrix
#   Input:  flowFrame list
#   Output: numeric matrix + sample labels
expr_data <- build_combined_expr_matrix(hg_frames, sample_names, markers_of_interest)

# Step 6: Leiden clustering + UMAP
#   Input:  numeric matrix
#   Output: data.frame with UMAP1, UMAP2, cluster, sample
leiden_result <- run_leiden_with_umap(
  expr_data$combined_matrix, expr_data$sample_labels,
  output_dir = "out/clustering/"
)

# Step 7: BRIDGE -- cluster labels to gating CSVs
#   Input:  matrix + labels + clusters
#   Output: per-sample CSVs on disk
clusters_to_gating_csv(
  expr_data$combined_matrix, expr_data$sample_labels,
  leiden_result$cluster, output_dir = "out/gating_csvs/"
)

# Step 8: Differential abundance (GLMM) + volcano
#   Input:  gating CSVs on disk + conditions mapping
#   Output: results data.frame + volcano plot
results <- run_differential_abundance(
  data_dir = "out/gating_csvs/",
  conditions_csv = "config/conditions.csv",
  output_dir = "out/abundance/"
)
```

### 4.2 CSV expression to UMAP + density + PCA

Lightweight exploration of CSV expression data with three visualizations,
no preprocessing or clustering needed.

```r
# Step 1: Load CSVs
#   Input:  .csv files on disk
#   Output: named list of data.tables
expr_list <- load_csv_expression(input_dir = "data/csv_expr/")

# Step 2: BRIDGE -- data.tables to matrix
markers <- c("CD3", "CD4", "CD8", "CD19", "CD56")
expr_matrix <- combine_expression_matrices(expr_list, markers, max_cells = 100000)
sample_labels <- rep(names(expr_list), sapply(expr_list, function(dt) {
  min(nrow(dt), 100000 / length(expr_list))
}))

# Step 3: UMAP (matrix path, no flowFrame needed)
plot_umap_from_matrix(expr_matrix, sample_labels, output_dir = "out/umap/")

# Step 4: Density (matrix path)
plot_density_from_matrix(expr_matrix, sample_labels, markers, output_dir = "out/density/")

# Step 5: PCA
plot_pca(expr_matrix, sample_labels, output_dir = "out/pca/")
```

### 4.3 CSV expression to clustering to gating CSV to enrichment analysis

Start from CSV expression data, cluster it, convert to gating CSVs, then
run enrichment statistics.

```r
# Step 1: Load CSVs
expr_list <- load_csv_expression(input_dir = "data/csv_expr/")

# Step 2: BRIDGE -- data.tables to matrix
markers <- c("CD3", "CD4", "CD8", "CD19", "CD56")
expr_matrix <- combine_expression_matrices(expr_list, markers)
sample_labels <- rep(names(expr_list), sapply(expr_list, nrow))

# Step 3: Leiden clustering
clusters <- run_leiden(expr_matrix, resolution = 0.8)

# Step 4: BRIDGE -- cluster labels to gating CSVs
clusters_to_gating_csv(expr_matrix, sample_labels, clusters,
                        output_dir = "out/gating_csvs/")

# Step 5: Enrichment analysis
enrichment <- run_cluster_enrichment(
  data_dir = "out/gating_csvs/",
  conditions_csv = "config/conditions.csv",
  output_dir = "out/enrichment/"
)

# Optional: volcano plot from enrichment results
plot_volcano(enrichment, output_dir = "out/enrichment/",
             fc_col = "log2_enrichment", title = "Cluster Enrichment")
```

### 4.4 FCS to matrix to heatmap + PCA (skip clustering)

Quick quality check of FCS data: load, harmonize, extract a matrix, and
visualize directly. No clustering required.

```r
# Step 1: Load
loaded <- load_fcs_raw(input_dir = "data/fcs/")

# Step 2: Harmonize
hg_frames <- homogenize_flowframes(
  loaded$flow_frames, loaded$sample_names, markers_of_interest
)

# Step 3: BRIDGE -- flowFrames to matrix
expr_data <- build_combined_expr_matrix(
  hg_frames, loaded$sample_names, markers_of_interest
)

# Step 4: Heatmap (median expression per sample)
# Build a sample-level summary matrix first
sample_medians <- do.call(rbind, lapply(split(
  as.data.frame(expr_data$combined_matrix),
  expr_data$sample_labels
), function(df) apply(df, 2, median)))

plot_heatmap(sample_medians, output_path = "out/heatmap/median_expr.png",
             title = "Median Expression per Sample")

# Step 5: PCA
plot_pca(expr_data$combined_matrix, expr_data$sample_labels,
         output_dir = "out/pca/")
```

### 4.5 CSV gating to abundance + enrichment to volcano

When you already have gating CSVs (from FlowJo export or a previous
pipeline run), go directly to statistical analysis.

```r
# Step 1: Load gating CSVs (produces processed CSVs on disk)
result <- load_csv_gating(
  input_fcs_path = "data/fcs/",
  input_csv_path = "data/csv_gating/",
  output_path = "out/loaded_gating/",
  processing_mode = "manual"
)

# Step 2a: Differential abundance
da_results <- run_differential_abundance(
  data_dir = "out/loaded_gating/",
  conditions_csv = "config/conditions.csv",
  output_dir = "out/abundance/"
)

# Step 2b: Cluster enrichment (can run in parallel with 2a)
enrich_results <- run_cluster_enrichment(
  data_dir = "out/loaded_gating/",
  conditions_csv = "config/conditions.csv",
  output_dir = "out/enrichment/"
)

# Step 3: Volcano plots for both
plot_volcano(da_results, output_dir = "out/abundance/",
             title = "Differential Abundance (GLMM)")
plot_volcano(enrich_results, output_dir = "out/enrichment/",
             fc_col = "log2_enrichment", title = "Cluster Enrichment (Fisher)")
```

### 4.6 FCS to cluster to gating CSV to split to expression matrices to heatmap

Full pipeline from raw FCS files through to per-cluster expression heatmaps.

```r
# Step 1: Load and preprocess
loaded <- load_fcs_raw(input_dir = "data/fcs/")
pp_frames <- run_preprocessing(loaded$flow_frames, params)
hg_frames <- homogenize_flowframes(pp_frames, loaded$sample_names, markers_of_interest)

# Step 2: BRIDGE -- flowFrames to matrix
expr_data <- build_combined_expr_matrix(hg_frames, loaded$sample_names, markers_of_interest)

# Step 3: Cluster
leiden_result <- run_leiden_with_umap(
  expr_data$combined_matrix, expr_data$sample_labels
)

# Step 4: BRIDGE -- clusters to gating CSVs
clusters_to_gating_csv(
  expr_data$combined_matrix, expr_data$sample_labels,
  leiden_result$cluster, output_dir = "out/gating_csvs/"
)

# Step 5: Split by cluster
split_result <- split_by_cluster(
  input_dir = "out/gating_csvs/",
  output_dir = "out/split_clusters/",
  markers_of_interest = markers_of_interest
)

# Step 6: Create per-group expression matrices
expr_matrices <- create_expression_matrices(
  input_dir = "out/split_clusters/",
  output_dir = "out/expr_matrices/"
)

# Step 7: Heatmap per cluster group
for (group_name in names(expr_matrices)) {
  mtx <- as.matrix(expr_matrices[[group_name]])
  plot_heatmap(mtx, output_path = paste0("out/heatmaps/", group_name, ".png"),
               title = paste("Expression:", group_name))
}

# Bonus: Abundance plots from the tracking matrix
plot_abundance(tracking_matrix = split_result$tracking_matrix,
               output_dir = "out/abundance_plots/")
```

### 4.7 FlowJo workspace comparison to Venn/UpSet

Compare gating across multiple FlowJo workspaces.

```r
# Step 1: Load multiple workspaces
pop_strings <- load_multiple_workspaces(
  wsp_paths = c("data/wsp/analyst1.wsp", "data/wsp/analyst2.wsp"),
  names = c("Analyst_1", "Analyst_2"),
  fcs_dirs = "data/fcs/"
)

# Step 2: Build binary membership matrix for each population
# (typically handled by metadata_qc/compare_gatings.R, but you can
#  construct it manually from the pop_strings lists)
for (sample_name in names(pop_strings[["Analyst_1"]])) {
  pops_1 <- pop_strings[["Analyst_1"]][[sample_name]]
  pops_2 <- pop_strings[["Analyst_2"]][[sample_name]]
  # Build index lists for overlap
  all_pops <- union(unique(pops_1), unique(pops_2))
  for (pop in all_pops[!is.na(all_pops)]) {
    sets <- list(
      Analyst_1 = which(pops_1 == pop),
      Analyst_2 = which(pops_2 == pop)
    )
    if (all(sapply(sets, length) > 0)) {
      plot_euler(sets, output_dir = "out/venn/",
                 filename = paste0(sample_name, "_", basename(pop)))
      plot_upset_from_list(sets, output_path = paste0("out/upset/",
                           sample_name, "_", basename(pop), ".png"))
    }
  }
}
```

### 4.8 CSV expression to flowFrame to UMAP (via bridge)

When you want to use the flowFrame-based `plot_umap()` with CSV data,
convert through the bridge function.

```r
# Step 1: Load CSVs
expr_list <- load_csv_expression(input_dir = "data/csv_expr/")

# Step 2: BRIDGE -- data.tables to flowFrames
markers <- c("CD3", "CD4", "CD8", "CD19", "CD56")
flow_frames <- create_flowframes_from_matrices(expr_list, markers_of_interest = markers)

# Step 3: flowFrame-based UMAP
plot_umap(flow_frames, names(flow_frames), markers_of_interest = markers,
          output_dir = "out/umap/")

# Note: for most cases, plot_umap_from_matrix() is simpler. Use this path
# only if you need to chain into other flowFrame-based scripts (e.g.,
# downsampling, harmonization) before plotting.
```

---

## 5. Decision Tree

Use this tree to determine which scripts to chain for your specific
situation.

```
What is your INPUT data?
|
+-- .fcs files
|   |
|   +-- Do you need quality control / transformation?
|   |   |
|   |   +-- YES --> load_fcs_raw() -> run_preprocessing()
|   |   |           -> [continue below with flowFrame list]
|   |   |
|   |   +-- NO  --> load_fcs_raw()
|   |               -> [continue below with flowFrame list]
|   |
|   +-- [flowFrame list] -- What is your GOAL?
|       |
|       +-- UMAP/density visualization
|       |   |
|       |   +-- Option A (direct): plot_umap() / plot_density()
|       |   +-- Option B (via matrix): build_combined_expr_matrix()
|       |       -> plot_umap_from_matrix() / plot_density_from_matrix()
|       |
|       +-- PCA / heatmap
|       |   -> build_combined_expr_matrix() -> plot_pca() / plot_heatmap()
|       |
|       +-- Clustering
|       |   -> build_combined_expr_matrix() -> run_leiden()
|       |      -> [continue below with cluster labels]
|       |
|       +-- Differential abundance / enrichment
|           -> build_combined_expr_matrix() -> run_leiden()
|              -> clusters_to_gating_csv()
|              -> run_differential_abundance() / run_cluster_enrichment()
|
+-- .csv expression files (no gating columns)
|   |
|   +-- load_csv_expression()
|   |
|   +-- What is your GOAL?
|       |
|       +-- UMAP/density/PCA (quick look)
|       |   -> combine_expression_matrices()
|       |      -> plot_umap_from_matrix() / plot_density_from_matrix() / plot_pca()
|       |
|       +-- Heatmap
|       |   -> combine_expression_matrices() -> plot_heatmap()
|       |
|       +-- Clustering then analysis
|       |   -> combine_expression_matrices() -> run_leiden()
|       |      -> clusters_to_gating_csv()
|       |      -> [see gating CSV goals below]
|       |
|       +-- Feed into flowFrame-based scripts
|           -> create_flowframes_from_matrices()
|              -> [treat as flowFrame list above]
|
+-- .csv files with gating/cluster columns
|   |
|   +-- load_csv_gating() or use CSVs directly
|   |
|   +-- What is your GOAL?
|       |
|       +-- Abundance analysis
|       |   -> build_abundance_matrix() or build_tracking_matrix()
|       |      -> run_differential_abundance() / plot_abundance()
|       |
|       +-- Enrichment analysis
|       |   -> run_cluster_enrichment()
|       |
|       +-- Split by cluster then analyze per-cluster
|       |   -> split_by_cluster()
|       |      -> create_expression_matrices() -> plot_heatmap()
|       |      -> run_marker_statistics()
|       |
|       +-- Composition barplots
|       |   -> split_by_cluster() -> tracking_matrix_to_long()
|       |      -> plot_composition_barplot()
|       |
|       +-- Cell type ratios
|       |   -> build_tracking_matrix() -> calculate_cell_ratios()
|       |
|       +-- Volcano plot
|           -> run_differential_abundance() -> plot_volcano()
|           or run_cluster_enrichment() -> plot_volcano(fc_col="log2_enrichment")
|
+-- .wsp FlowJo workspace files
    |
    +-- Single workspace
    |   -> load_flowjo_wsp() -> get_population_stats()
    |
    +-- Compare multiple workspaces
        -> load_multiple_workspaces()
           -> compare_gatings() / compute_fleiss_kappa()
           -> plot_euler() / plot_upset() / plot_venn()
```

---

## 6. Quick Reference Table

For each script, this table shows: the primary function(s) exported, what
they require as input, what they produce as output, and which downstream
scripts they connect to.

### Loading

| Script | Function | Input | Output | Connects to |
|--------|----------|-------|--------|-------------|
| `load_fcs_raw.R` | `load_fcs_raw(input_dir)` | Directory of `.fcs` files | `list(flow_frames, sample_names, input_files)` | `run_preprocessing()`, `downsample_flowframes()`, `homogenize_flowframes()`, `build_combined_expr_matrix()`, `plot_umap()`, `plot_density()` |
| `load_csv_expression.R` | `load_csv_expression(input_dir)` | Directory of `.csv` expression files | Named list of data.tables | `combine_expression_matrices()`, `create_flowframes_from_matrices()` |
| `load_csv_gating.R` | `load_csv_gating(input_fcs_path, ...)` | Directory of `.csv` + `.fcs` files | Gating CSVs on disk (via `output_path`) | `split_by_cluster()`, `build_abundance_matrix()`, `run_differential_abundance()`, `run_cluster_enrichment()` |
| `load_flowjo_wsp.R` | `load_flowjo_wsp(wsp_path)` | `.wsp` file + FCS directory | `list(gating_set, sample_names, population_names)` | `get_population_stats()`, `get_population_strings_list()`, `import_flowjo_gates()` |
| `load_flowjo_wsp.R` | `load_multiple_workspaces(wsp_paths, names)` | Multiple `.wsp` files | Named list of population string lists | `compare_gatings()`, `compute_fleiss_kappa()`, `plot_euler()`, `plot_venn()`, `plot_upset()` |
| `load_manual_gating.R` | `load_manual_gating(csv_dir)` | Directory of CSV files with gating columns | Named list of label vectors (one per sample) | Custom downstream analysis |
| `create_flowframes.R` | `create_flowframes_from_matrices(expr_list)` | Named list of matrices/data.frames | Named list of flowFrame objects | `plot_umap()`, `plot_density()`, `downsample_flowframes()`, `homogenize_flowframes()` |

### Preprocessing

| Script | Function | Input | Output | Connects to |
|--------|----------|-------|--------|-------------|
| `preprocess_pipeline.R` | `run_preprocessing(flow_frames, params)` | Named list of flowFrames + params list | Named list of preprocessed flowFrames | `downsample_flowframes()`, `homogenize_flowframes()`, `build_combined_expr_matrix()`, `plot_umap()`, `plot_density()` |
| `downsampling.R` | `downsample_flowframes(flow_frames, sample_names, target_cells)` | Named list of flowFrames | Named list of downsampled flowFrames | Same as above |
| `transformation.R` | `transform_flowframe(flowframe, sample_name, markers, offsets, widths)` | Single flowFrame | Transformed flowFrame | Used within `run_preprocessing()` |
| `normalization.R` | `z_score_scaling(flowframe, marker_indices)` | Single flowFrame | Scaled flowFrame | Used within `run_preprocessing()` |
| `normalization.R` | `quantile_norm(flowframe, marker_indices)` | Single flowFrame | Normalized flowFrame | Used within `run_preprocessing()` |
| `normalization.R` | `harmony_normalize(expr_matrix, batch_vector)` | Numeric matrix + batch labels | Corrected matrix | `run_leiden()`, `plot_umap_from_matrix()` |
| `quality_control.R` | `quality_control(flowframe, sample_name)` | Single flowFrame | QC-filtered flowFrame | Used within `run_preprocessing()` |
| `compensation.R` | `apply_compensation(flowframe, comp_matrix)` | Single flowFrame | Compensated flowFrame | `quality_control()`, `transform_flowframe()` |
| `spectral_unmixing.R` | `unmix_spectral(flowframe, unmixing_matrix)` | Single flowFrame + unmixing matrix | Unmixed flowFrame | `quality_control()`, `transform_flowframe()` |
| `bead_normalization.R` | `normalize_beads(flowframe, bead_channels)` | Single flowFrame | `list(normalized_ff, bead_indices, normalization_factors)` | `quality_control()`, `transform_flowframe()` |
| `doublet_detection.R` | `detect_doublets(flowframe)` | Single flowFrame | `list(filtered_ff, doublet_indices, singlet_rate)` | `quality_control()`, `transform_flowframe()` |
| `live_dead_gating.R` | `gate_live_cells(flowframe, viability_channel)` | Single flowFrame | `list(live_ff, dead_indices, viability_rate)` | `quality_control()`, `transform_flowframe()` |
| `supercell.R` | `run_supercell(expr_matrix, sample_labels, markers)` | Numeric matrix + sample labels | `list(supercell_matrix, supercell_labels, cell_map, supercell_sizes)` | `run_leiden()`, `run_flowsom()`, `plot_umap_from_matrix()` |
| `automated_gating.R` | `run_auto_gating(flowframe, method, channels)` | Single flowFrame | `list(gated_ff, gate_coords, population_stats)` | Downstream analysis |
| `gate_utils.R` | `apply_gate_sequence(flowframe, gate_list)` | Single flowFrame + gate list | `list(gated_ff, gate_indices, population_stats)` | Downstream analysis |

### Harmonization

| Script | Function | Input | Output | Connects to |
|--------|----------|-------|--------|-------------|
| `panel_homogenization.R` | `homogenize_flowframes(flow_frames, sample_names, markers)` | Named list of flowFrames | Named list of homogenized flowFrames | `build_combined_expr_matrix()`, `run_batch_correction()`, `plot_umap()`, `plot_density()` |
| `batch_correction.R` | `run_batch_correction(flow_frames, sample_names, batch_labels)` | Named list of homogenized flowFrames | `list(corrected_matrix, sample_labels, batch_labels, markers)` | `run_leiden()`, `plot_umap_from_matrix()`, `plot_pca()` |
| `cytonorm.R` | `run_cytonorm(train_frames, validation_frames, batch_labels, markers)` | flowFrame lists + batch labels | `list(corrected_frames, model)` | `transform_flowframe()`, `run_leiden()` |
| `combat.R` | `run_combat(expr_matrix, batch_vector, covariates)` | Numeric matrix + batch vector | Corrected numeric matrix | `run_leiden()`, `plot_umap_from_matrix()` |

### Clustering

| Script | Function | Input | Output | Connects to |
|--------|----------|-------|--------|-------------|
| `leiden.R` | `run_leiden(expr_matrix)` | Numeric matrix | Integer vector of cluster IDs | `clusters_to_gating_csv()` |
| `leiden.R` | `run_leiden_with_umap(expr_matrix, sample_labels)` | Numeric matrix + sample labels | Data.frame with UMAP1, UMAP2, cluster, sample | `clusters_to_gating_csv()`, `plot_umap_from_matrix()` (reuse UMAP coords) |
| `leiden.R` | `clusters_to_gating_csv(matrix, labels, clusters, output_dir)` | Matrix + sample labels + cluster IDs | Per-sample gating CSVs on disk | `split_by_cluster()`, `build_abundance_matrix()`, `run_differential_abundance()`, `run_cluster_enrichment()` |
| `flowjo_import.R` | `import_flowjo_gates(wsp_path)` | `.wsp` file + FCS directory | `list(gating_set, assignments)` with binary population columns | `get_cluster_assignments()` |
| `cluster_assignment.R` | `get_cluster_assignments(source, ...)` | Depends on source type | Data.frame with cell_id, cluster, cluster_name | Custom downstream analysis |
| `flowsom.R` | `run_flowsom(expr_matrix, n_clusters)` | Numeric matrix | Integer vector of cluster IDs | `clusters_to_gating_csv()`, `validate_clusters()` |
| `flowsom.R` | `run_flowsom_flowframes(flow_frames, markers, n_clusters)` | flowFrame list | `list(cluster_ids, flowsom_obj, sample_labels)` | `clusters_to_gating_csv()` |
| `phenograph.R` | `run_phenograph(expr_matrix, k)` | Numeric matrix | Integer vector of cluster IDs | `clusters_to_gating_csv()`, `validate_clusters()` |
| `dbscan.R` | `run_dbscan(expr_matrix, eps, min_pts)` | Numeric matrix | `list(cluster_ids, noise_indices, n_clusters)` | `clusters_to_gating_csv()`, `validate_clusters()` |
| `meta_clustering.R` | `run_resolution_sweep(expr_matrix, method, resolutions)` | Numeric matrix | `list(best_clusters, best_param, scores_df, all_clusterings)` | `clusters_to_gating_csv()` |
| `meta_clustering.R` | `consensus_clusters(cluster_list)` | Named list of cluster vectors | Integer vector of consensus cluster IDs | `clusters_to_gating_csv()` |
| `cluster_validation.R` | `validate_clusters(expr_matrix, cluster_ids)` | Numeric matrix + cluster IDs | Data.frame with metric values | Standalone report |
| `cell_annotation.R` | `annotate_clusters(expr_matrix, cluster_ids, markers)` | Numeric matrix + cluster IDs | Data.frame with cluster, annotation, confidence | Visualization, analysis |

### Analysis

| Script | Function | Input | Output | Connects to |
|--------|----------|-------|--------|-------------|
| `cluster_splitting.R` | `split_by_cluster(input_dir, output_dir)` | Directory of gating CSVs | `list(tracking_matrix, cluster_summary, output_files)` + per-cluster CSVs on disk | `tracking_matrix_to_long()`, `create_expression_matrices()`, `run_marker_statistics()`, `plot_abundance()` |
| `cluster_splitting.R` | `tracking_matrix_to_long(tracking_matrix)` | Tracking matrix (wide) | Long data.frame | `plot_composition_barplot()` |
| `differential_abundance.R` | `build_abundance_matrix(data_dir)` | Directory of gating CSVs | `list(abundance_df, all_clusters, valid_clusters)` | `run_differential_abundance()`, `run_cluster_enrichment()`, `plot_abundance()` |
| `differential_abundance.R` | `run_differential_abundance(data_dir, conditions_csv, output_dir)` | Directory of gating CSVs + conditions mapping | Results data.frame + volcano plot | `plot_volcano()` |
| `differential_abundance.R` | `plot_volcano(results_df, output_dir)` | Data.frame with log2fc, p_adjusted, cluster | Volcano plot on disk | Standalone output |
| `cluster_enrichment.R` | `run_cluster_enrichment(data_dir, conditions_csv, output_dir)` | Directory of gating CSVs + conditions mapping | Results data.frame with Fisher/ME statistics | `plot_volcano(fc_col="log2_enrichment")` |
| `expression_matrices.R` | `create_expression_matrices(input_dir, output_dir)` | Directory of cluster-split CSVs | Named list of per-group data.tables + CSVs on disk | `plot_heatmap()`, `plot_cluster_heatmap()` |
| `expression_matrices.R` | `combine_expression_matrices(expr_list, markers)` | Named list of data.tables | Numeric matrix | `plot_heatmap()`, `plot_pca()`, `plot_umap_from_matrix()`, `run_leiden()` |
| `marker_statistics.R` | `run_marker_statistics(input_folder, output_folder)` | Directory of per-population CSVs | Results data.frame + CSV on disk | Standalone output |
| `abundance_plotting.R` | `build_tracking_matrix(input_dir, cluster_regex, ...)` | Directory of gating CSVs | Tracking matrix (wide data.frame) | `plot_abundance()`, `tracking_matrix_to_long()`, `calculate_cell_ratios()` |
| `abundance_plotting.R` | `plot_abundance(output_dir, tracking_matrix)` | Tracking matrix or path to CSV | Stacked barplot, proportional barplot, heatmap, totals plot | Standalone output |
| `cell_ratio.R` | `calculate_cell_ratios(abundance_csv, num, den)` | Abundance CSV (tracking matrix on disk) | Data.frame with ratio, log2_ratio per sample | Standalone output |
| `differential_expression.R` | `run_differential_expression(input_dir, conditions_csv, output_dir)` | Gating CSVs + conditions mapping | Data.frame with cluster, marker, log2fc, pvalue, padj | `plot_heatmap()` |
| `effect_sizes.R` | `compute_effect_sizes(expr_data, condition_col, value_cols)` | data.table with markers + condition | Data.frame with marker, effect_size, ci, method | Standalone report |
| `correlation_analysis.R` | `run_correlation_analysis(expr_matrix, method)` | Numeric matrix | `list(cor_matrix, p_matrix, significant_pairs)` | `plot_heatmap()` |
| `diffcyt_analysis.R` | `run_diffcyt(flow_frames, experiment_info, markers)` | flowFrame list + experiment info | `list(results, top_table, cluster_assignments)` | Standalone output |
| `survival_analysis.R` | `run_survival_analysis(abundance_matrix, clinical_csv, output_dir)` | Abundance data.frame + clinical CSV | Data.frame with cluster, hr, ci, pvalue, padj | Standalone output |

### Visualization

| Script | Function | Input | Output | Connects to |
|--------|----------|-------|--------|-------------|
| `umap.R` | `plot_umap(flow_frames, sample_names, ...)` | Named list of flowFrames | UMAP plots on disk + `list(umap_coords, plots)` | Standalone output |
| `umap.R` | `plot_umap_from_matrix(expr_matrix, ...)` | Numeric matrix | UMAP plots on disk + `list(umap_coords, plots)` | Standalone output |
| `density_plots.R` | `plot_density(flow_frames, sample_names, markers, ...)` | Named list of flowFrames | Density/ridge plots on disk | Standalone output |
| `density_plots.R` | `plot_density_from_matrix(expr_matrix, ...)` | Numeric matrix | Density/ridge plots on disk | Standalone output |
| `heatmaps.R` | `plot_heatmap(expr_matrix, output_path)` | Numeric matrix (rows x markers) | Heatmap on disk | Standalone output |
| `heatmaps.R` | `plot_cluster_heatmap(expr_data, cluster_col, marker_cols, ...)` | Data.frame with expression + cluster column | Median expression heatmap on disk | Standalone output |
| `barplots.R` | `plot_composition_barplot(data, output_dir)` | Long data.frame with sample, cluster, count | Barplot on disk | Standalone output |
| `pca.R` | `plot_pca(expr_matrix, sample_labels, output_dir)` | Numeric matrix | PCA scatter + loadings plots on disk | Standalone output |
| `upset_plots.R` | `plot_upset(data, set_columns, output_path)` | Data.frame with binary columns | UpSet plot on disk | Standalone output |
| `upset_plots.R` | `plot_upset_from_list(sets, output_path)` | Named list of vectors | UpSet plot on disk | Standalone output |
| `venn_diagrams.R` | `plot_euler(sets, output_dir, filename)` | Named list of index vectors | Euler diagram on disk | Standalone output |
| `venn_diagrams.R` | `plot_venn(sets, filename, output_dir)` | Named list of index vectors | Venn diagram on disk | Standalone output |
| `tsne.R` | `plot_tsne(flow_frames, sample_names, ...)` | Named list of flowFrames | t-SNE plots on disk + `list(tsne_coords, plots)` | Standalone output |
| `tsne.R` | `plot_tsne_from_matrix(expr_matrix, ...)` | Numeric matrix | t-SNE plots on disk + `list(tsne_coords, plots)` | Standalone output |
| `diffusion_map.R` | `plot_diffusion_map(expr_matrix, ...)` | Numeric matrix | Diffusion map plots + `list(dm_coords, eigenvalues, plots)` | Standalone output |
| `biaxial_plots.R` | `plot_biaxial(flowframe, x_marker, y_marker, ...)` | flowFrame or numeric matrix | Biaxial scatter/pseudocolor plot | Standalone output |
| `violin_plots.R` | `plot_violin(expr_data, markers, group_col, ...)` | data.table with markers + grouping col | Violin/box plots on disk | Standalone output |
| `radar_plots.R` | `plot_radar(cluster_medians, markers, ...)` | data.frame (cluster x marker medians) | Radar/spider charts on disk | Standalone output |
| `contour_plots.R` | `plot_contour(coords, group_labels, ...)` | data.frame with x, y columns | 2D contour/density plots on disk | Standalone output |

### Shared Utilities (bridge functions)

| Script | Function | Input | Output | Purpose |
|--------|----------|-------|--------|---------|
| `R_shared/flowframe_utils.R` | `build_combined_expr_matrix(flow_frames, sample_names, markers)` | flowFrame list | `list(combined_matrix, sample_labels)` | flowFrame list to numeric matrix |
| `R_shared/flowframe_utils.R` | `extract_exprmtx_list(flow_frames, sample_names, markers)` | flowFrame list | Named list of matrices | flowFrame list to per-sample matrices |
| `R_shared/flowframe_utils.R` | `from_mtx_to_ffr(mtx)` | Numeric matrix | Single flowFrame | Matrix to flowFrame (single sample) |
| `R_shared/stats_utils.R` | `run_umap(expr_matrix, ...)` | Numeric matrix | 2-column UMAP coordinate matrix | Reusable UMAP computation |
| `R_shared/stats_utils.R` | `run_pca(expr_matrix, ...)` | Numeric matrix | `list(coordinates, variance_explained, loadings)` | Reusable PCA computation |

---

## 7. Method Swapping

Three dispatcher modules in `src/R_shared/` let you swap any method at any
pipeline step by changing a single config parameter. The dispatchers
normalize inputs and outputs so downstream code needs no changes.

### 7.1 Clustering: `run_clustering()`

**Location:** `src/R_shared/clustering_dispatch.R`

Set `clustering_method` in your config to one of:

| Value | Routes to | Key params |
|-------|-----------|------------|
| `"leiden"` | `run_leiden()` | `n_neighbors`, `resolution` |
| `"flowsom"` | `run_flowsom()` | `n_clusters`, `grid_xdim`, `grid_ydim` |
| `"phenograph"` | `run_phenograph()` | `k` (or `n_neighbors`) |
| `"dbscan"` | `run_dbscan()` | `dbscan_eps`, `dbscan_min_pts` |

**Output:** Always an integer vector of `length == nrow(expr_matrix)`.
DBSCAN's list output is normalized: `cluster_ids` extracted, metadata stored
as `attr(ids, "dbscan_metadata")`.

```r
source(file.path(shared_dir, "clustering_dispatch.R"))
ids <- run_clustering(expr_matrix, method = params$clustering_method,
                       params = params)
```

### 7.2 Dimensionality Reduction: `run_dimred()`

**Location:** `src/R_shared/dimred_dispatch.R`

Set `dimred_method` in your config to one of:

| Value | Routes to | Key params |
|-------|-----------|------------|
| `"umap"` | `plot_umap_from_matrix()` | `n_neighbors`, `umap_min_dist` |
| `"tsne"` | `plot_tsne_from_matrix()` | `tsne_perplexity`, `tsne_max_iter` |
| `"diffusion_map"` | `plot_diffusion_map()` | `dm_n_components`, `dm_sigma` |

**Output:** Always `list(coords, method, plots, extra)` where `coords` is
a numeric matrix (cells x dims).

```r
source(file.path(shared_dir, "dimred_dispatch.R"))
result <- run_dimred(expr_matrix, method = params$dimred_method,
                      sample_labels = labels, cluster_labels = ids,
                      output_dir = "out/", params = params)
```

### 7.3 Batch Correction: `run_batch_correction_dispatch()`

**Location:** `src/R_shared/batch_correction_dispatch.R`

Set `batch_method` in your config to one of:

| Value | Routes to | Accepts | Key params |
|-------|-----------|---------|------------|
| `"harmony"` | `run_batch_correction()` | flowFrames or matrix | `batch_do_pca` |
| `"cytonorm"` | `run_cytonorm()` | flowFrames only | `cytonorm_train_indices`, `cytonorm_n_clusters` |
| `"combat"` | `run_combat()` | flowFrames or matrix | `combat_method` |
| `"none"` | pass-through | flowFrames or matrix | — |

The dispatcher automatically bridges between flowFrame and matrix inputs:
if Harmony needs flowFrames and you pass a matrix, it converts via
`from_mtx_to_ffr()`. If ComBat needs a matrix and you pass flowFrames, it
converts via `build_combined_expr_matrix()`.

**Output:** Always `list(corrected_matrix, sample_labels, batch_labels,
markers, corrected_frames)`.

```r
source(file.path(shared_dir, "batch_correction_dispatch.R"))
result <- run_batch_correction_dispatch(
  data = homogenized, method = params$batch_method,
  batch_labels = params$batch_names,
  markers_of_interest = params$markers_of_interest,
  sample_names = sample_names, params = params)
```

### 7.4 Config Example

To switch Pipeline A from Leiden + UMAP to FlowSOM + t-SNE with ComBat
batch correction, change three lines in your config:

```r
params <- list(
  clustering_method = "flowsom",    # was: "leiden"
  batch_method      = "combat",     # was: "none"
  dimred_method     = "tsne",       # was: "umap"
  # ... rest of config unchanged
)
```
