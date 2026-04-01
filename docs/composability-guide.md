# Composability Guide

This guide explains how to combine scripts from `src/` into custom workflows outside of the five built-in pipelines. It covers the data flow between pipeline stages, the bridge functions that convert between data types, the contracts each stage expects, and complete example workflows you can adapt.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Data Flow Diagram](#2-data-flow-diagram)
3. [Data Type Contracts](#3-data-type-contracts)
4. [Bridge Functions](#4-bridge-functions)
5. [Method Dispatchers](#5-method-dispatchers)
6. [Example Workflows](#6-example-workflows)
7. [Common Pitfalls](#7-common-pitfalls)

---

## 1. Architecture Overview

The pipeline has six processing stages. Each stage is a directory under `src/` containing standalone scripts. Every script documents its upstream ("Input from") and downstream ("Output to") connections in its roxygen headers.

The five built-in pipelines (`pipelines/pipeline_A_raw_fcs.R` through `pipeline_E_preprocessed_umap.R`) are thin orchestrators that chain these scripts in a fixed order. You are not limited to those orderings. Any script can be called independently, and data can flow between scripts from different pipelines, as long as you match the expected data types at each boundary.

**Bridge functions** handle the type conversions between stages. For example, `build_combined_expr_matrix()` converts a list of flowFrame objects into the numeric matrix that clustering and visualization functions expect.

---

## 2. Data Flow Diagram

```
                         INPUT FORMATS
           +-------------+-------------+-------------+
           |             |             |             |
       Raw FCS     Preprocessed   CSV Expression  FlowJo WSP
       files       FCS files      matrices        workspaces
           |             |             |             |
           v             v             v             v
  +--------+------+ +----+----+ +-----+------+ +----+-----+
  | load_fcs_raw  | | load_   | | load_csv_  | | load_    |
  |               | | fcs_raw | | expression | | flowjo_  |
  +-------+-------+ +----+----+ +-----+------+ | wsp      |
          |              |             |        +----+-----+
          |              |             |             |
          v              |             |             v
  +-------+-------+     |             |      +------+------+
  | PREPROCESSING |     |             |      | METADATA_QC |
  | qc, transform,|     |             |      | fleiss_kappa|
  | doublets,     |     |             |      | gating_comp |
  | live/dead,    |     |             |      +------+------+
  | normalization |     |             |             |
  +-------+-------+     |             |             v
          |              |             |      Venn/Euler/UpSet
          v              v             |
  +-------+--------------+---+        |
  |     HARMONIZATION        |        |
  | compare_panels           |        |
  | homogenize_flowframes    |        |
  | batch correction         |        |
  | (Harmony/CytoNorm/ComBat)|       |
  +-------+------------------+        |
          |                           |
          v                           |
  +-------+-------+                   |
  | build_combined| <----- BRIDGE --- + create_flowframes
  | _expr_matrix  |       FUNCTIONS     _from_matrices
  +-------+-------+                       (reverse)
          |
          +------- run_supercell() -------> compressed matrix
          |        (optional)                     |
          v                                       v
  +-------+----------------------------+---------+
  |              CLUSTERING                       |
  | run_leiden / run_flowsom / run_phenograph /    |
  | run_dbscan / run_clustering (dispatcher)      |
  | run_resolution_sweep / consensus_clusters     |
  +---+---+-----------+-----------+-------+-------+
      |   |           |           |       |
      |   |           v           |       |
      |   |   validate_clusters   |       |
      |   |           |           |       |
      |   |           v           |       |
      |   |   annotate_clusters   |       |
      |   |                       |       |
      v   v                       v       v
  +---+---+--------+    +--------+-------+-------+
  | clusters_to_   |    |     VISUALIZATION       |
  | gating_csv     |    | plot_umap_from_matrix   |
  +-------+--------+    | plot_tsne_from_matrix   |
          |              | plot_diffusion_map      |
          v              | plot_density_from_matrix|
  +-------+--------+    | plot_heatmap            |
  |    ANALYSIS    |    | plot_biaxial            |
  | cluster_split  |    | plot_violin / radar     |
  | diff_abundance |    | plot_contour / barplots |
  | enrichment     |    | plot_pca / upset / venn |
  | marker_stats   |    +-------------------------+
  | diff_expression|
  | effect_sizes   |
  | correlations   |
  | diffcyt        |
  | survival       |
  +----------------+
```

### Key data type boundaries

| Boundary | Upstream type | Bridge function | Downstream type |
|----------|--------------|-----------------|-----------------|
| Loading to Clustering | `list(flowFrame)` | `build_combined_expr_matrix()` | `numeric matrix` + `character vector` (labels) |
| CSV to flowFrame path | `list(matrix)` | `create_flowframes_from_matrices()` | `list(flowFrame)` |
| Clustering to Analysis | `integer vector` + `numeric matrix` | `clusters_to_gating_csv()` | Per-sample CSV files |
| Any to Visualization | `numeric matrix` | `plot_umap_from_matrix()`, etc. | Plot files |
| Full-size to Compressed | `numeric matrix` | `run_supercell()` | Smaller `numeric matrix` |
| Clustering to Validation | `integer vector` + `numeric matrix` | `validate_clusters()` | `data.frame` of metrics |
| Clustering to Annotation | `integer vector` + `numeric matrix` | `annotate_clusters()` | `data.frame` with labels |

---

## 3. Data Type Contracts

Every script in the pipeline expects and produces one of a small number of standard data types. Matching these types at each boundary is the key to composability.

### 3.1 flowFrame list

**What it is.** A named R list where each element is a `flowCore::flowFrame` object. Names are sample identifiers.

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
`build_combined_expr_matrix()`, `plot_umap()`, `plot_tsne()`,
`plot_density()`, `run_flowsom_flowframes()`, `run_cytonorm()`,
`run_diffcyt()`

### 3.2 Numeric expression matrix

**What it is.** A standard R numeric matrix with rows as cells and columns as markers. Often accompanied by a parallel `sample_labels` character vector of the same length as `nrow(matrix)`.

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
`build_combined_expr_matrix()$combined_matrix`,
`run_supercell()$supercell_matrix`,
`run_batch_correction_dispatch()$corrected_matrix`,
`run_combat()`, `flowCore::exprs(ff)`

**Consumed by:**
`run_leiden()`, `run_flowsom()`, `run_phenograph()`, `run_dbscan()`,
`run_clustering()` (dispatcher),
`plot_umap_from_matrix()`, `plot_tsne_from_matrix()`,
`plot_density_from_matrix()`, `plot_diffusion_map()`,
`plot_heatmap()`, `plot_pca()`,
`validate_clusters()`, `annotate_clusters()`,
`run_supercell()`, `clusters_to_gating_csv()`

### 3.3 Cluster ID vector

**What it is.** An integer vector of cluster assignments, one per cell.

**Shape:** `integer(n_cells)` where values are 1-based cluster IDs. DBSCAN uses 0 for noise points. Length must equal `nrow(expr_matrix)`.

**Produced by:** `run_leiden()`, `run_flowsom()`, `run_phenograph()`,
`run_dbscan()$cluster_ids`, `run_clustering()` (dispatcher),
`consensus_clusters()`, `run_resolution_sweep()$best_clusters`

**Consumed by:** `clusters_to_gating_csv()`, `validate_clusters()`,
`annotate_clusters()`,
`plot_umap_from_matrix(cluster_labels = ...)`,
`plot_tsne_from_matrix(cluster_labels = ...)`,
`run_supercell_per_cluster()`

**Important:** Leiden returns `numeric`, not `integer`. The `run_clustering()` dispatcher coerces to `integer` automatically. If calling `run_leiden()` directly, wrap with `as.integer()`.

### 3.4 Sample labels vector

**What it is.** A character vector mapping each cell row to its source sample.

**Shape:** `character(n_cells)` where `length(sample_labels) == nrow(expr_matrix)`.

**Produced by:** `build_combined_expr_matrix()$sample_labels`,
`run_batch_correction_dispatch()$sample_labels`,
`run_supercell()$supercell_labels`

**Consumed by:** `clusters_to_gating_csv()`, visualization functions,
`run_supercell()`

### 3.5 Per-sample gating CSVs

**What it is.** Semicolon-separated CSV files, one per sample, with expression columns plus binary indicator columns for each cluster (named `Cluster_X (Clustering annotated)`).

**Structure (semicolon-separated):**

```
CD3;CD4;CD8;CD19;CD56;Cluster_1 (Clustering annotated);Cluster_2 (Clustering annotated)
3.21;2.45;0.12;0.05;1.87;1;0
3.15;0.08;2.91;0.03;0.22;0;1
```

**Produced by:** `clusters_to_gating_csv()`, `load_csv_gating()`

**Consumed by:** `split_by_cluster()`, `run_differential_abundance()`,
`run_cluster_enrichment()`, `build_abundance_matrix()`,
`build_tracking_matrix()`

### 3.6 Validation metrics data frame

**What it is.** Data frame of cluster quality metrics.

**Shape:** `data.frame(metric, value, description)` with one row per computed metric.

**Produced by:** `validate_clusters()`

**Consumed by:** User inspection, reporting scripts. Can be saved with `data.table::fwrite()`.

### 3.7 Annotation data frame

**What it is.** Data frame mapping clusters to cell type labels.

**Shape:** `data.frame(cluster, annotation, confidence, median_marker1, median_marker2, ...)` with one row per cluster.

**Produced by:** `annotate_clusters()`

**Consumed by:** Can be joined to cluster ID vectors for labeled plots, or used to rename clusters before `clusters_to_gating_csv()`.

---

## 4. Bridge Functions

### 4.1 `build_combined_expr_matrix()`

**Location:** `src/R_shared/flowframe_utils.R`

**Purpose:** Converts a list of flowFrame objects into a single numeric matrix with parallel sample labels. This is the primary bridge from the flowFrame world (loading, preprocessing, harmonization) into the matrix world (clustering, visualization, analysis).

**Signature:**

```r
build_combined_expr_matrix(flow_frames_list, sample_names, markers_of_interest)
# Returns: list(combined_matrix, sample_labels)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `flow_frames_list` | `list(flowFrame)` | Named list of flowFrame objects |
| `sample_names` | `character` | One name per flowFrame |
| `markers_of_interest` | `character` | Marker names to extract (must match `markernames()` in the flowFrames) |

**Returns:** A list with two elements:
- `combined_matrix`: Numeric matrix (cells x markers), all samples row-bound
- `sample_labels`: Character vector (length = `nrow(combined_matrix)`) tracking which sample each row came from

**Example:**

```r
source("src/R_shared/flowframe_utils.R")
result <- build_combined_expr_matrix(flow_frames, sample_names, markers)
# result$combined_matrix: 50000 x 12 numeric matrix
# result$sample_labels:   character vector of length 50000
```

---

### 4.2 `create_flowframes_from_matrices()`

**Location:** `src/loading/create_flowframes.R`

**Purpose:** The reverse of `build_combined_expr_matrix()`. Converts expression matrices (or data.frames) into flowFrame objects, enabling CSV-loaded data to enter the flowFrame-based path (preprocessing, harmonization, flowFrame-based visualization).

**Signature:**

```r
create_flowframes_from_matrices(expr_list, markers_of_interest = NULL,
                                save_fcs = FALSE, output_dir = NULL)
# Returns: list(flowFrame)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_list` | `named list(matrix)` | Named list of expression matrices or data.frames |
| `markers_of_interest` | `character` or `NULL` | Subset columns to these markers (NULL = keep all) |
| `save_fcs` | `logical` | If TRUE, also write FCS files to disk |
| `output_dir` | `character` or `NULL` | Directory for FCS output (required if `save_fcs = TRUE`) |

**Returns:** Named list of flowFrame objects.

**Example:**

```r
source("src/loading/create_flowframes.R")
# Load CSVs, then convert to flowFrames
expr_list <- list(
  sample_A = as.matrix(data.table::fread("sample_A.csv")),
  sample_B = as.matrix(data.table::fread("sample_B.csv"))
)
flow_frames <- create_flowframes_from_matrices(expr_list, markers_of_interest = markers)
# Now these can go into homogenize_flowframes(), plot_umap(), etc.
```

---

### 4.3 `clusters_to_gating_csv()`

**Location:** `src/clustering/leiden.R`

**Purpose:** Takes a combined expression matrix, sample labels, and cluster IDs, and writes one CSV per sample with binary indicator columns for each cluster (format: `Cluster_X (Clustering annotated)`). This converts clustering output into the CSV gating format that Pipeline C analysis scripts expect.

**Signature:**

```r
clusters_to_gating_csv(combined_matrix, sample_labels, cluster_ids,
                       output_dir, sep = ";")
# Returns: character vector of file paths (invisibly)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `combined_matrix` | `numeric matrix` | Cells x markers |
| `sample_labels` | `character` | One per row, identifies source sample |
| `cluster_ids` | `integer` or `character` | One per row, cluster assignment |
| `output_dir` | `character` | Directory for output CSVs |
| `sep` | `character` | CSV separator (default `";"`) |

**Returns:** Character vector of written file paths (invisibly).

**Output format:** Each CSV contains all original expression columns plus one binary column per cluster. For example, if Leiden found 8 clusters, each CSV will have 8 additional columns named `Cluster_1 (Clustering annotated)` through `Cluster_8 (Clustering annotated)`, where each cell has a `1` in exactly one cluster column and `0` in all others.

**Example:**

```r
source("src/clustering/leiden.R")
clusters_to_gating_csv(
  combined_matrix = expr_data$combined_matrix,
  sample_labels   = expr_data$sample_labels,
  cluster_ids     = cluster_ids,
  output_dir      = "output/gating_csvs/"
)
# Writes: output/gating_csvs/sample_A.csv, output/gating_csvs/sample_B.csv, ...
```

---

### 4.4 `run_supercell()`

**Location:** `src/preprocessing/supercell.R`

**Purpose:** Compresses single-cell data into supercells (metacells) using SuperCellCyto. Reduces the number of rows by a factor of `gam` (default 20x) while preserving population structure. Insert this between `build_combined_expr_matrix()` and clustering/visualization to handle large datasets.

**Signature:**

```r
run_supercell(expr_matrix, sample_labels, markers,
              gam = 20, k_knn = 30, n_pc = 10,
              do_pca = TRUE, seed = 42)
# Returns: list(supercell_matrix, supercell_labels, cell_map, supercell_sizes)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_matrix` | `numeric matrix` | Cells x markers |
| `sample_labels` | `character` | One per row |
| `markers` | `character` | Marker column names to use for compression |
| `gam` | `numeric` | Graining level: higher = fewer supercells (20 = ~1/20 compression) |
| `k_knn` | `integer` | Nearest neighbors for graph construction |
| `n_pc` | `integer` | Number of PCs for initial reduction |
| `do_pca` | `logical` | Run PCA before supercell construction |
| `seed` | `integer` | Random seed |

**Returns:** A list with:
- `supercell_matrix`: Compressed numeric matrix (supercells x markers)
- `supercell_labels`: Sample labels for each supercell
- `cell_map`: `data.table` mapping original cell IDs to supercell IDs
- `supercell_sizes`: `data.table` with the number of cells in each supercell

**Also available:** `run_supercell_per_cluster()` applies compression independently within each cluster, preserving cluster boundaries. Useful when clusters have very different densities or when you need cluster-aware metacells.

**Example:**

```r
source("src/preprocessing/supercell.R")
sc <- run_supercell(expr_matrix, sample_labels, markers, gam = 20)
# sc$supercell_matrix: 2500 x 12 (down from 50000 x 12)
# Now cluster the compressed matrix:
cluster_ids <- run_leiden(sc$supercell_matrix)
```

---

### 4.5 `plot_umap_from_matrix()` / `plot_tsne_from_matrix()`

**Location:** `src/visualization/umap.R` / `src/visualization/tsne.R`

**Purpose:** Generate dimensionality reduction visualizations directly from a numeric matrix, without requiring flowFrame objects. These are the matrix-based counterparts to `plot_umap()` and `plot_tsne()`.

**Signature (`plot_umap_from_matrix`):**

```r
plot_umap_from_matrix(expr_matrix, sample_labels = NULL,
                      markers_of_interest = NULL,
                      cluster_labels = NULL,
                      condition_labels = NULL,
                      output_dir, n_neighbors = 15,
                      min_dist = 0.1, umap_coords = NULL,
                      output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                      plot_width = 10, plot_height = 8, seed = 42)
# Returns: list(umap_coords, sample_labels, plots)
```

**Signature (`plot_tsne_from_matrix`):**

```r
plot_tsne_from_matrix(expr_matrix, sample_labels = NULL,
                      cluster_labels = NULL,
                      condition_labels = NULL,
                      output_dir, perplexity = 30,
                      max_iter = 1000, tsne_coords = NULL,
                      output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                      plot_width = 10, plot_height = 8, seed = 42)
# Returns: list(tsne_coords, sample_labels, plots)
```

**Key features:**
- Accept optional `cluster_labels` to overlay cluster assignments on the embedding
- Accept optional `condition_labels` (UMAP) for faceted condition plots
- Accept optional pre-computed coordinates (`umap_coords` / `tsne_coords`) to skip recomputation
- Automatically exclude Time, FSC, SSC, and index columns from computation
- Generate per-marker expression overlay plots

**Example:**

```r
source("src/visualization/umap.R")
umap_result <- plot_umap_from_matrix(
  expr_matrix    = expr_data$combined_matrix,
  sample_labels  = expr_data$sample_labels,
  cluster_labels = cluster_ids,
  output_dir     = "output/umap/"
)
# Saves: umap_sample.png, umap_cluster.png, umap_<marker>.png for each marker
# Returns umap_coords for reuse in other visualizations
```

---

### 4.6 Clustering Functions

**Location:** `src/clustering/` (individual files) or `src/R_shared/clustering_dispatch.R` (dispatcher)

All four clustering methods accept a numeric expression matrix and return cluster assignments:

| Function | File | Key Parameters | Return type | Notes |
|----------|------|----------------|-------------|-------|
| `run_leiden()` | `clustering/leiden.R` | `n_neighbors` (30), `resolution` (1.0), `seed` (42) | `numeric` vector | Coerce with `as.integer()` when calling directly |
| `run_flowsom()` | `clustering/flowsom.R` | `n_clusters` (20), `grid_xdim` (10), `grid_ydim` (10), `seed` (42) | `integer` vector | Deterministic cluster count via ConsensusClusterPlus |
| `run_phenograph()` | `clustering/phenograph.R` | `k` (30), `seed` (42) | `integer` vector | KNN + SNN graph + Louvain community detection |
| `run_dbscan()` | `clustering/dbscan.R` | `eps` (auto), `min_pts` (5), `auto_eps` (TRUE), `seed` (42) | **list** with `$cluster_ids` | Cluster 0 = noise; use `$cluster_ids` to extract vector |

**Contract:** All methods accept `matrix(cells x markers)` and produce cluster assignments of `length == nrow(expr_matrix)`. The `run_clustering()` dispatcher enforces a uniform `integer` vector output for all methods.

**FlowSOM also provides** `run_flowsom_flowframes()` which accepts a flowFrame list directly and handles the matrix conversion internally, returning `list(cluster_ids, flowsom_obj, sample_labels)`.

---

### 4.7 `validate_clusters()`

**Location:** `src/clustering/cluster_validation.R`

**Purpose:** Compute internal validation metrics (silhouette width, Calinski-Harabasz index, Davies-Bouldin index) and, when ground-truth labels are available, external metrics (adjusted Rand index). Works with any clustering method's output.

**Signature:**

```r
validate_clusters(expr_matrix, cluster_ids,
                  true_labels = NULL,
                  metrics = c("silhouette", "calinski_harabasz"),
                  sample_size = 5000, seed = 42)
# Returns: data.frame(metric, value, description)
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `expr_matrix` | `numeric matrix` | Cells x markers |
| `cluster_ids` | `integer` | Cluster assignments (length = nrow) |
| `true_labels` | `integer`/`character` or `NULL` | Ground-truth labels for ARI computation |
| `metrics` | `character` | Which metrics: `"silhouette"`, `"calinski_harabasz"`, `"davies_bouldin"` |
| `sample_size` | `integer` | Max cells for metric computation (subsamples if exceeded) |

**Interpretation of metrics:**

| Metric | Range | Better is | Measures |
|--------|-------|-----------|----------|
| Silhouette | -1 to 1 | Higher | How well each cell fits its cluster vs. nearest neighbor cluster |
| Calinski-Harabasz | 0 to Inf | Higher | Ratio of between-cluster to within-cluster dispersion |
| Davies-Bouldin | 0 to Inf | Lower | Average similarity between each cluster and its most similar cluster |
| Adjusted Rand Index | -1 to 1 | Higher (1 = perfect) | Agreement between predicted and true labels |

**Example:**

```r
source("src/clustering/cluster_validation.R")
quality <- validate_clusters(expr_matrix, cluster_ids,
                             metrics = c("silhouette", "calinski_harabasz", "davies_bouldin"))
print(quality)
#              metric     value                                          description
# 1        n_clusters        12                                  Number of clusters
# 2           n_cells     50000                              Total number of cells
# 3        silhouette    0.3421  Mean silhouette width (-1 to 1; higher = better)
# 4 calinski_harabasz   245.120  Calinski-Harabasz index (higher = better defined)
# 5    davies_bouldin    1.2340  Davies-Bouldin index (lower = better separation)
```

---

### 4.8 `annotate_clusters()`

**Location:** `src/clustering/cell_annotation.R`

**Purpose:** Assign cell type labels to clusters based on their marker expression profiles. Three annotation methods are available.

**Signature:**

```r
annotate_clusters(expr_matrix, cluster_ids,
                  markers = NULL, reference = "auto",
                  method = "threshold")
# Returns: data.frame(cluster, annotation, confidence, median_marker1, ...)
```

**Annotation methods:**

| Method | How it works | When to use |
|--------|-------------|-------------|
| `"threshold"` | Markers above the global median are called positive (+), below negative (-). Builds phenotype strings like `CD3+CD4+CD8-`. Confidence = fraction of markers with clear (abs(z) > 1) calls. | Quick exploratory annotation, well-separated populations |
| `"correlation"` | Cosine similarity between cluster median profiles and a user-provided reference matrix. Assigns the best-matching reference cell type. | When you have a known reference panel (e.g., from the literature or a previous experiment) |
| `"manual"` | Returns median expression per cluster without automatic labels. User inspects the medians and assigns labels manually. | Complex panels, ambiguous populations |

**Example:**

```r
source("src/clustering/cell_annotation.R")
annotations <- annotate_clusters(expr_matrix, cluster_ids, method = "threshold")
# annotations$cluster:    1, 2, 3, ...
# annotations$annotation: "CD3+CD4+", "CD3+CD8+", "CD19+", ...
# annotations$confidence: 0.85, 0.92, 0.78, ...

# With a reference matrix for correlation-based annotation:
ref_profiles <- matrix(c(1, 1, 0, 0, 0,    # T helper
                          1, 0, 1, 0, 0,    # cytotoxic T
                          0, 0, 0, 1, 0),   # B cell
                        nrow = 3, byrow = TRUE,
                        dimnames = list(c("T_helper", "Cytotoxic_T", "B_cell"),
                                        c("CD3", "CD4", "CD8", "CD19", "CD56")))
annotations <- annotate_clusters(expr_matrix, cluster_ids,
                                 reference = ref_profiles, method = "correlation")
```

---

## 5. Method Dispatchers

Three dispatcher functions in `src/R_shared/` allow swapping methods via a single config parameter. Each dispatcher validates input, sources the appropriate module, normalizes output format, and handles fallbacks.

### 5.1 Clustering Dispatcher: `run_clustering()`

**Location:** `src/R_shared/clustering_dispatch.R`

**Config parameter:** `clustering_method` (or pass `method` directly)

```r
source(file.path(shared_dir, "clustering_dispatch.R"))
cluster_ids <- run_clustering(expr_matrix, method = "flowsom",
                              params = list(n_clusters = 20, seed = 42))
```

| Method | Config value | Key params | Package |
|--------|-------------|------------|---------|
| Leiden | `"leiden"` | `n_neighbors` (30), `resolution` (1.0) | leiden, igraph, uwot |
| FlowSOM | `"flowsom"` | `n_clusters` (20), `grid_xdim` (10), `grid_ydim` (10) | FlowSOM |
| PhenoGraph | `"phenograph"` | `k` / `n_neighbors` (30) | igraph, uwot |
| DBSCAN | `"dbscan"` | `dbscan_eps` / `eps` (auto), `dbscan_min_pts` / `min_pts` (5) | dbscan |

**Output normalization:**
- All methods return `integer(n_cells)` (Leiden's numeric output is coerced).
- DBSCAN's list output is unwrapped; metadata (noise indices, eps used) is stored as `attr(cluster_ids, "dbscan_metadata")`.
- Post-conditions enforced: `is.integer(cluster_ids)`, `length(cluster_ids) == nrow(expr_matrix)`.

### 5.2 Dimensionality Reduction Dispatcher: `run_dimred()`

**Location:** `src/R_shared/dimred_dispatch.R`

**Config parameter:** `dimred_method` (or pass `method` directly)

```r
source(file.path(shared_dir, "dimred_dispatch.R"))
result <- run_dimred(expr_matrix, method = "tsne",
                     sample_labels = labels, output_dir = "output/",
                     params = list(perplexity = 30))
# result$coords: matrix of coordinates (cells x dims)
# result$method: "tsne"
# result$plots:  list of ggplot objects (or NULL if output_dir is NULL)
# result$extra:  NULL (or list(eigenvalues = ...) for diffusion_map)
```

| Method | Config value | Key params | Package |
|--------|-------------|------------|---------|
| UMAP | `"umap"` | `n_neighbors` / `umap_n_neighbors` (15), `umap_min_dist` / `min_dist` (0.1) | uwot |
| t-SNE | `"tsne"` | `tsne_perplexity` / `perplexity` (30), `tsne_max_iter` / `max_iter` (1000) | Rtsne |
| Diffusion Map | `"diffusion_map"` | `dm_n_components` / `n_components` (3), `dm_sigma` / `sigma` ("local") | destiny |

**Output normalization:** All methods return `list(coords, method, plots, extra)` where `nrow(coords) == nrow(expr_matrix)`.

### 5.3 Batch Correction Dispatcher: `run_batch_correction_dispatch()`

**Location:** `src/R_shared/batch_correction_dispatch.R`

**Config parameter:** `batch_method` (or pass `method` directly)

```r
source(file.path(shared_dir, "batch_correction_dispatch.R"))
result <- run_batch_correction_dispatch(
  data                = flow_frames,   # or expr_matrix
  method              = "combat",
  batch_labels        = batch_labels,
  markers_of_interest = markers,
  sample_names        = sample_names
)
# result$corrected_matrix: numeric matrix (cells x markers)
# result$sample_labels:    character vector
# result$batch_labels:     character vector
# result$markers:          character vector
# result$corrected_frames: list(flowFrame) (CytoNorm only; NULL otherwise)
```

| Method | Config value | Accepts | Key params | Package |
|--------|-------------|---------|------------|---------|
| Harmony | `"harmony"` | flowFrame list or matrix | `batch_do_pca` (FALSE) | harmony |
| CytoNorm | `"cytonorm"` | flowFrame list **only** | `cytonorm_train_indices`, `cytonorm_n_clusters` (10), `cytonorm_goal` ("mean") | CytoNorm |
| ComBat | `"combat"` | flowFrame list or matrix | `combat_covariates` (NULL), `combat_method` ("parametric") | sva |
| None | `"none"` | flowFrame list or matrix | (pass-through) | -- |

**Input flexibility:** The dispatcher automatically converts between flowFrame list and matrix as needed for each backend. If you pass a matrix to Harmony, it creates temporary flowFrames via `from_mtx_to_ffr()`. If you pass flowFrames to ComBat, it extracts the expression matrix via `build_combined_expr_matrix()`.

**CytoNorm caveat:** Requires flowFrame input and only produces corrected frames for the validation subset (not training frames). Specify `cytonorm_train_indices` in your config.

### 5.4 Config Example

To switch Pipeline A from Leiden + UMAP + no batch correction to FlowSOM + t-SNE + ComBat, change three lines in your config:

```r
params <- list(
  clustering_method = "flowsom",    # was: "leiden"
  batch_method      = "combat",     # was: "none"
  dimred_method     = "tsne",       # was: "umap"
  n_clusters        = 15,           # FlowSOM-specific
  tsne_perplexity   = 40,           # t-SNE-specific
  # ... rest of config unchanged
)
```

---

## 6. Example Workflows

### 6.1 Load CSV Data, Cluster, and Visualize

A minimal workflow for expression CSV files that skips the flowFrame path entirely.

```r
# --- Setup ---
shared_dir <- "/path/to/cytometry/src/R_shared"
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(shared_dir, "clustering_dispatch.R"))

# --- Step 1: Load CSV expression data ---
source("src/loading/load_csv_expression.R")
expr_data <- load_csv_expression(input_dir = "data/csv_expression/")
# expr_data: named list of data.frames, one per sample

# --- Step 2: Combine into a single matrix ---
all_matrices <- lapply(expr_data, as.matrix)
combined <- do.call(rbind, all_matrices)
sample_labels <- rep(names(all_matrices), sapply(all_matrices, nrow))

# --- Step 3: Select markers of interest ---
markers <- c("CD3", "CD4", "CD8", "CD19", "CD14", "CD56")
expr_matrix <- combined[, markers, drop = FALSE]

# --- Step 4: Cluster via dispatcher ---
cluster_ids <- run_clustering(expr_matrix, method = "leiden",
                              params = list(resolution = 0.8))

# --- Step 5: Visualize with UMAP ---
source("src/visualization/umap.R")
plot_umap_from_matrix(expr_matrix, sample_labels = sample_labels,
                      cluster_labels = cluster_ids,
                      output_dir = "output/umap/")

# --- Step 6: Export gating CSVs for downstream analysis ---
source("src/clustering/leiden.R")
clusters_to_gating_csv(expr_matrix, sample_labels, cluster_ids,
                       output_dir = "output/gating_csvs/")
```

---

### 6.2 Load FCS, Preprocess, Compress with SuperCell, Cluster, Annotate, Analyze

A full workflow for large datasets that uses SuperCellCyto to reduce computational cost before clustering and analysis.

```r
# --- Setup ---
shared_dir <- "/path/to/cytometry/src/R_shared"
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(shared_dir, "clustering_dispatch.R"))

markers <- c("CD3", "CD4", "CD8", "CD19", "CD14", "CD56",
             "CD16", "CD25", "HLA-DR", "CD45RA", "CD127", "CD38")

# --- Step 1: Load raw FCS ---
source("src/loading/load_fcs_raw.R")
flow_data <- load_fcs_raw(input_dir = "data/fcs_files/")

# --- Step 2: Preprocess (QC + transformation) ---
source("src/preprocessing/preprocess_pipeline.R")
preprocessed <- run_preprocessing(
  flow_frames = flow_data$flow_frames,
  params      = list(qc_type = "flowAI", transform_type = "arcsinh",
                     cofactor = 5, pp_order = "QT"),
  output_dir  = "output/preprocessing/"
)

# --- Step 3: Homogenize panels across samples ---
source("src/harmonization/panel_homogenization.R")
homogenized <- homogenize_flowframes(
  flow_frames         = preprocessed,
  sample_names        = flow_data$sample_names,
  markers_of_interest = markers,
  output_dir          = "output/harmonized/"
)

# --- Step 4: Build combined expression matrix (BRIDGE) ---
expr_data <- build_combined_expr_matrix(homogenized, flow_data$sample_names, markers)

# --- Step 5: Compress with SuperCellCyto ---
source("src/preprocessing/supercell.R")
sc <- run_supercell(expr_data$combined_matrix, expr_data$sample_labels,
                    markers, gam = 20)
# 500,000 cells -> ~25,000 supercells

# --- Step 6: Cluster the compressed data ---
cluster_ids <- run_clustering(sc$supercell_matrix, method = "flowsom",
                              params = list(n_clusters = 15))

# --- Step 7: Annotate clusters with cell type labels ---
source("src/clustering/cell_annotation.R")
annotations <- annotate_clusters(sc$supercell_matrix, cluster_ids,
                                 markers = markers, method = "threshold")
cat("Cluster annotations:\n")
print(annotations[, c("cluster", "annotation", "confidence")])

# --- Step 8: Validate clustering quality ---
source("src/clustering/cluster_validation.R")
quality <- validate_clusters(sc$supercell_matrix, cluster_ids,
                             metrics = c("silhouette", "calinski_harabasz"))
print(quality)

# --- Step 9: Visualize ---
source("src/visualization/umap.R")
plot_umap_from_matrix(sc$supercell_matrix, sample_labels = sc$supercell_labels,
                      cluster_labels = cluster_ids,
                      output_dir = "output/umap/")

# --- Step 10: Export gating CSVs and run differential abundance ---
source("src/clustering/leiden.R")
clusters_to_gating_csv(sc$supercell_matrix, sc$supercell_labels, cluster_ids,
                       output_dir = "output/gating_csvs/")

source("src/analysis/differential_abundance.R")
da_results <- run_differential_abundance(
  data_dir       = "output/gating_csvs/",
  conditions_csv = "config/conditions.csv",
  output_dir     = "output/abundance/"
)
```

---

### 6.3 Compare Multiple Clustering Methods with Validation

Run Leiden, FlowSOM, PhenoGraph, and DBSCAN on the same data, compare their quality metrics, and build a consensus clustering.

```r
# --- Setup ---
shared_dir <- "/path/to/cytometry/src/R_shared"
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
source(file.path(shared_dir, "clustering_dispatch.R"))

# --- Assume expr_matrix and sample_labels are already prepared ---
# (from build_combined_expr_matrix() or CSV loading)

# --- Run all four methods via the dispatcher ---
leiden_ids <- run_clustering(expr_matrix, method = "leiden",
                             params = list(resolution = 1.0, n_neighbors = 30))

flowsom_ids <- run_clustering(expr_matrix, method = "flowsom",
                              params = list(n_clusters = 15))

phenograph_ids <- run_clustering(expr_matrix, method = "phenograph",
                                 params = list(k = 30))

dbscan_ids <- run_clustering(expr_matrix, method = "dbscan",
                              params = list(min_pts = 5))

# --- Validate each method ---
source("src/clustering/cluster_validation.R")
methods <- list(leiden = leiden_ids, flowsom = flowsom_ids,
                phenograph = phenograph_ids, dbscan = dbscan_ids)
all_metrics <- list()

for (name in names(methods)) {
  cat("\n=== Validating:", name, "===\n")
  metrics <- validate_clusters(expr_matrix, methods[[name]],
                               metrics = c("silhouette", "calinski_harabasz",
                                           "davies_bouldin"))
  metrics$method <- name
  all_metrics[[name]] <- metrics
}

# Compare all methods side by side
comparison <- do.call(rbind, all_metrics)
print(comparison[, c("method", "metric", "value")])

# --- Build consensus from three graph-based methods ---
source("src/clustering/meta_clustering.R")
consensus_ids <- consensus_clusters(
  cluster_list = list(leiden = leiden_ids, flowsom = flowsom_ids,
                      phenograph = phenograph_ids),
  method = "majority"
)

# --- Visualize each clustering on the same UMAP coordinates ---
source("src/visualization/umap.R")

# Compute UMAP once
umap_result <- plot_umap_from_matrix(expr_matrix, sample_labels = sample_labels,
                                     cluster_labels = leiden_ids,
                                     output_dir = "output/umap_leiden/")

# Reuse coordinates for other methods (avoids recomputing UMAP)
for (name in c("flowsom", "phenograph", "dbscan", "consensus")) {
  ids <- if (name == "consensus") consensus_ids else methods[[name]]
  plot_umap_from_matrix(expr_matrix, sample_labels = sample_labels,
                        cluster_labels = ids,
                        umap_coords = umap_result$umap_coords,
                        output_dir = paste0("output/umap_", name, "/"))
}

# --- Optional: compare to ground truth if available ---
# quality_vs_truth <- validate_clusters(
#   expr_matrix, leiden_ids, true_labels = manual_labels,
#   metrics = c("silhouette", "calinski_harabasz")
# )
```

---

### 6.4 Load FlowJo Workspaces, Compare Gatings, Compute Agreement

Compare gating strategies from multiple FlowJo analysts on the same FCS files.

```r
# --- Setup ---
shared_dir <- "/path/to/cytometry/src/R_shared"
source(file.path(shared_dir, "io_utils.R"))

# --- Step 1: Load multiple FlowJo workspaces ---
source("src/loading/load_flowjo_wsp.R")
workspaces <- load_multiple_workspaces(
  wsp_paths = c("data/analyst_A.wsp", "data/analyst_B.wsp", "data/analyst_C.wsp"),
  names     = c("Analyst_A", "Analyst_B", "Analyst_C"),
  fcs_dirs  = "data/fcs_files/"
)

# --- Step 2: Compare gating across analysts ---
source("src/metadata_qc/gating_comparison.R")
comparison <- compare_gating_strategies(
  gating_sets            = workspaces$gating_sets,
  analyst_names          = c("Analyst_A", "Analyst_B", "Analyst_C"),
  populations_of_interest = c("CD3+", "CD4+", "CD8+", "CD19+"),
  output_dir             = "output/gating_comparison/"
)

# --- Step 3: Compute inter-rater agreement (Fleiss' kappa) ---
source("src/metadata_qc/fleiss_kappa.R")
kappa_results <- compute_fleiss_kappa(
  gating_matrices = comparison$gating_matrices,
  output_dir      = "output/agreement/"
)
cat("Fleiss' kappa:", kappa_results$kappa, "\n")
# Interpretation: < 0.20 poor, 0.21-0.40 fair, 0.41-0.60 moderate,
#                 0.61-0.80 substantial, > 0.80 almost perfect

# --- Step 4: Visualize overlap with Euler/Venn diagrams ---
source("src/visualization/venn_diagrams.R")
for (pop in c("CD3+", "CD4+", "CD8+")) {
  plot_euler(
    sets      = comparison$population_sets[[pop]],
    set_names = c("Analyst_A", "Analyst_B", "Analyst_C"),
    title     = paste("Gating agreement:", pop),
    output_dir = "output/venn/"
  )
}

# --- Step 5: UpSet plot for multi-population intersection ---
source("src/visualization/upset_plots.R")
plot_upset_from_list(
  sets        = comparison$population_sets,
  output_path = "output/upset/all_populations.png"
)
```

---

## 7. Common Pitfalls

### 7.1 Type mismatches at boundaries

**Problem:** Passing a flowFrame list to a function that expects a matrix (or vice versa).

**Symptoms:** Errors like `"is.matrix(expr_matrix) is not TRUE"` or `"object of type 'S4' is not subsettable"`.

**Fix:** Use the appropriate bridge function. If you have flowFrames and need a matrix, call `build_combined_expr_matrix()`. If you have a matrix and need flowFrames, call `create_flowframes_from_matrices()`.

### 7.2 Leiden returns numeric, not integer

**Problem:** `run_leiden()` returns a numeric vector. Downstream functions that check `is.integer()` will fail.

**Fix:** Use the `run_clustering()` dispatcher, which coerces Leiden's output to integer. If calling `run_leiden()` directly, wrap the result: `cluster_ids <- as.integer(run_leiden(...))`.

### 7.3 DBSCAN returns a list, not a vector

**Problem:** `run_dbscan()` returns a list with `$cluster_ids`, `$noise_indices`, `$n_clusters`, `$eps_used`, and `$noise_fraction`. Passing the list directly to downstream functions fails.

**Fix:** Extract the integer vector: `result <- run_dbscan(expr_matrix); ids <- result$cluster_ids`. The `run_clustering()` dispatcher handles this automatically and stores DBSCAN metadata as `attr(ids, "dbscan_metadata")`. Also note that DBSCAN uses cluster ID 0 for noise points, which some downstream functions may not expect.

### 7.4 Missing panel homogenization

**Problem:** Different FCS files have different marker sets. Combining them with `build_combined_expr_matrix()` may produce a matrix with fewer markers than expected, or emit `"No matching markers found"` warnings.

**Fix:** Always run `homogenize_flowframes()` before `build_combined_expr_matrix()` when working with multiple FCS files from different instruments or panels. The homogenization step renames markers to a common set and drops channels that cannot be harmonized.

### 7.5 Forgetting sample_labels after compression

**Problem:** After supercell compression, the original `sample_labels` vector no longer matches the rows of the compressed matrix.

**Fix:** Always use the labels that travel with the output. After `run_supercell()`, use `sc$supercell_labels` (not the original `sample_labels`). After `run_batch_correction_dispatch()`, use `result$sample_labels`. The labels vector and the matrix must always have the same number of rows.

### 7.6 CytoNorm requires flowFrame input

**Problem:** `run_batch_correction_dispatch()` with `method = "cytonorm"` fails if you pass a matrix.

**Symptoms:** Error: `"CytoNorm requires flowFrame input. Provide a list of flowFrame objects."`.

**Fix:** Pass a list of flowFrame objects. If your data is in matrix form, use `create_flowframes_from_matrices()` first. Alternatively, use `method = "harmony"` or `method = "combat"`, which accept both input types.

### 7.7 Large datasets exhaust memory

**Problem:** `build_combined_expr_matrix()` concatenates all samples into one matrix. With millions of cells, this can exceed available RAM.

**Fix:** Two options:
1. **Downsample first:** Use `downsample_flowframes()` before building the combined matrix.
2. **Compress with SuperCellCyto:** Call `run_supercell()` on the combined matrix to reduce rows by 10-20x before clustering.

For reading large CSV files, use `read_csv_chunked()` from `src/R_shared/memory_utils.R`.

### 7.8 Shared directory path resolution

**Problem:** Scripts use `shared_dir <- file.path(dirname(getwd()), "src", "R_shared")` to locate shared utilities. If your working directory is not a direct subdirectory of the project root, this path will be wrong.

**Fix:** Set `shared_dir` explicitly before sourcing any scripts:

```r
shared_dir <- "/absolute/path/to/cytometry/src/R_shared"
```

### 7.9 `%||%` operator not defined in R < 4.0

**Problem:** The null-coalescing operator `%||%` is used throughout dispatchers and pipelines. In R versions before 4.0, it is not defined.

**Fix:** The dispatcher files define a polyfill automatically: `if (!exists("%||%")) \`%||%\` <- function(a, b) if (is.null(a)) b else a`. If you write custom scripts that use `%||%`, source any dispatcher file first, or add the polyfill yourself.

### 7.10 Passing the wrong markers to downstream steps

**Problem:** After batch correction or supercell compression, the set of marker columns may differ from what you specified initially (e.g., CytoNorm may reorder columns, or missing markers get dropped).

**Fix:** Always check `colnames()` on the output matrix and use `result$markers` from dispatcher output. Do not assume column order is preserved across bridge functions.
