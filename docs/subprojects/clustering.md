[< Back to Documentation Index](../README.md)

# src/clustering/ -- Clustering, Validation & Annotation

Nine scripts providing multiple clustering algorithms (Leiden, FlowSOM, PhenoGraph, DBSCAN), resolution sweep and consensus clustering, cluster validation metrics, cell type annotation, FlowJo gate import, and a unified cluster assignment interface.

## Summary Table

| Script | Functions | Description |
|--------|-----------|-------------|
| `leiden.R` | `run_leiden()`, `run_leiden_with_umap()` | Leiden community detection + UMAP (CSV export via `R_shared/cluster_utils.R::clusters_to_gating_csv()`) |
| `flowsom.R` | `run_flowsom()`, `run_flowsom_flowframes()` | FlowSOM self-organizing map meta-clustering |
| `phenograph.R` | `run_phenograph()` | PhenoGraph KNN + Louvain community detection |
| `dbscan.R` | `run_dbscan()`, `estimate_eps()` | DBSCAN density-based clustering with auto-eps |
| `meta_clustering.R` | `run_resolution_sweep()`, `consensus_clusters()`, `consensus_via_sampling()`, `compute_silhouette()`, `compute_calinski_harabasz()` | Resolution optimization and consensus clustering |
| `cluster_validation.R` | `validate_clusters()`, `compute_davies_bouldin()`, `compute_ari_manual()` | Internal and external validation metrics |
| `cell_annotation.R` | `annotate_clusters()`, `compute_cluster_medians()`, `cosine_similarity()` | Threshold / correlation / manual cell type annotation |
| `flowjo_import.R` | `import_flowjo_gates()` | Import gates from FlowJo workspaces |
| `cluster_assignment.R` | `get_cluster_assignments()`, `get_clusters_from_leiden()`, `get_clusters_from_flowjo()`, `get_clusters_from_csv()`, `get_clusters_from_flowsom()` | Unified cluster label interface |

---

## leiden.R

Runs Leiden community detection on expression matrices, with optional UMAP embedding. This is the primary clustering method used by Pipelines A and E.

### `run_leiden()`

```r
run_leiden(expr_matrix, n_neighbors = 30, resolution = 1.0, seed = 42)
```

Builds a KNN graph via `uwot`, constructs a sparse adjacency matrix, and runs Leiden community detection via `leiden::leiden()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `n_neighbors` | integer | `30` | KNN graph neighbors |
| `resolution` | numeric | `1.0` | Leiden resolution (higher = more clusters) |
| `seed` | integer | `42` | Random seed |

**Returns:** Integer vector of cluster assignments (one per cell).

**Input from:** `build_combined_expr_matrix()` (R_shared/flowframe_utils.R).
**Output to:** `clusters_to_gating_csv()` (this file), `validate_clusters()` (cluster_validation.R).

### `run_leiden_with_umap()`

```r
run_leiden_with_umap(expr_matrix, sample_labels = NULL, output_dir = NULL,
                     n_neighbors = 30, resolution = 1.0, seed = 42)
```

Combines UMAP embedding and Leiden clustering in one call. Saves results to `leiden_clusters.csv`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell (for output) |
| `output_dir` | character | `NULL` | Directory to save `leiden_clusters.csv` |
| `n_neighbors` | integer | `30` | Neighbors for both KNN graph and UMAP |
| `resolution` | numeric | `1.0` | Leiden resolution |
| `seed` | integer | `42` | Random seed |

**Returns:** `data.frame` with columns `UMAP1`, `UMAP2`, `cluster`, `sample`.

**Output File -- `leiden_clusters.csv`:**

| Column | Type | Description |
|--------|------|-------------|
| `UMAP1` | numeric | First UMAP coordinate |
| `UMAP2` | numeric | Second UMAP coordinate |
| `cluster` | integer | Leiden cluster assignment (1, 2, 3, ...) |
| `sample` | character | Sample of origin |

**Resolution Guide:**

| Resolution | Effect |
|------------|--------|
| 0.1 -- 0.3 | Very few, broad clusters |
| 0.5 -- 0.8 | Moderate number (good starting point) |
| 1.0 -- 1.5 | Many clusters |
| 2.0+ | Fine-grained sub-populations |

**Example:**

```r
source("leiden.R")

result <- run_leiden_with_umap(
  expr_matrix   = combined,
  sample_labels = labels,
  output_dir    = "/output/clustering",
  n_neighbors   = 15,
  resolution    = 0.8,
  seed          = 42
)
```

> `clusters_to_gating_csv()` lives in `src/R_shared/cluster_utils.R` — see [r-shared-reference.md](../r-shared-reference.md#clusters_to_gating_csv) for full docs.

**Packages:** `leiden`, `igraph`, `uwot`, `Matrix`

---

## flowsom.R

FlowSOM self-organizing map with consensus meta-clustering.

### `run_flowsom()`

```r
run_flowsom(expr_matrix, n_clusters = 20, grid_xdim = 10,
            grid_ydim = 10, seed = 42)
```

Builds a SOM grid on the expression data, then applies ConsensusClusterPlus to determine meta-cluster assignments. Falls back to manual SOM + consensus clustering if the primary method fails.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `n_clusters` | integer | `20` | Target number of meta-clusters |
| `grid_xdim` | integer | `10` | SOM grid x dimension |
| `grid_ydim` | integer | `10` | SOM grid y dimension |
| `seed` | integer | `42` | Random seed |

**Returns:** Integer vector of cluster assignments (length = nrow(expr_matrix)).

**Input from:** `build_combined_expr_matrix()`.
**Output to:** `clusters_to_gating_csv()`, `validate_clusters()`, `annotate_clusters()`.

### `run_flowsom_flowframes()`

```r
run_flowsom_flowframes(flow_frames, markers_of_interest,
                       n_clusters = 20, seed = 42)
```

Extracts expression matrices from flowFrames, concatenates, runs FlowSOM, and returns cluster assignments with sample-level tracking.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | List of flowFrame objects |
| `markers_of_interest` | character vector | *required* | Markers to cluster on |
| `n_clusters` | integer | `20` | Target meta-clusters |
| `seed` | integer | `42` | Random seed |

**Returns:**

```r
list(
  cluster_ids   = integer vector,
  flowsom_obj   = <FlowSOM object>,
  sample_labels = character vector
)
```

**Input from:** `load_fcs_raw()`, `homogenize_flowframes()`.
**Output to:** `clusters_to_gating_csv()`, `run_differential_abundance()`.

**Packages:** `FlowSOM`, `flowCore`

---

## phenograph.R

PhenoGraph-style clustering: KNN graph construction with shared nearest-neighbor (SNN) weighting and Louvain community detection.

### `run_phenograph()`

```r
run_phenograph(expr_matrix, k = 30, seed = 42)
```

Builds a shared nearest-neighbor graph from a KNN search using `uwot`, computes Jaccard-weighted edges, and runs Louvain community detection via `igraph::cluster_louvain()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `k` | integer | `30` | Number of nearest neighbors |
| `seed` | integer | `42` | Random seed |

**Returns:** Integer vector of cluster assignments (length = nrow(expr_matrix)).

**Input from:** `build_combined_expr_matrix()`.
**Output to:** `clusters_to_gating_csv()`, `validate_clusters()`, `annotate_clusters()`.

**Packages:** `igraph`, `uwot`

---

## dbscan.R

DBSCAN density-based clustering for rare population detection. Noise points receive cluster ID 0.

### `run_dbscan()`

```r
run_dbscan(expr_matrix, eps = NULL, min_pts = 5,
           auto_eps = TRUE, seed = 42)
```

Performs DBSCAN clustering. When `eps = NULL` and `auto_eps = TRUE`, estimates epsilon from the k-distance plot knee point. Falls back to scaled data if the initial run fails.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `eps` | numeric | `NULL` | Epsilon radius; `NULL` = auto-estimate |
| `min_pts` | integer | `5` | Minimum points to form a dense region |
| `auto_eps` | logical | `TRUE` | Estimate eps from k-distance knee point |
| `seed` | integer | `42` | Random seed for subsampling during eps estimation |

**Returns:**

```r
list(
  cluster_ids    = integer vector,   # 0 = noise
  noise_indices  = integer vector,
  n_clusters     = integer,
  eps_used       = numeric,
  noise_fraction = numeric
)
```

**Input from:** `build_combined_expr_matrix()`.
**Output to:** `clusters_to_gating_csv()`, `validate_clusters()`, `annotate_clusters()`.

### `estimate_eps()`

```r
estimate_eps(expr_matrix, k = 5, max_sample = 10000, seed = 42)
```

Computes k-nearest-neighbor distances, sorts them, and finds the knee point using the maximum curvature method.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers matrix |
| `k` | integer | `5` | Number of nearest neighbors |
| `max_sample` | integer | `10000` | Maximum cells to sample for estimation |
| `seed` | integer | `42` | Random seed |

**Returns:** Numeric scalar: estimated epsilon value.

**Packages:** `dbscan`

---

## meta_clustering.R

Resolution parameter optimization via sweep and consensus clustering from multiple independent runs.

### `run_resolution_sweep()`

```r
run_resolution_sweep(expr_matrix, method = "leiden",
                     resolutions = seq(0.2, 2.0, by = 0.2),
                     metric = "silhouette",
                     n_neighbors = 30, seed = 42)
```

Tries multiple resolution values, evaluates each with an internal validation metric, and returns the best clustering.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers matrix |
| `method` | character | `"leiden"` | Clustering method |
| `resolutions` | numeric vector | `seq(0.2, 2.0, by = 0.2)` | Resolution values to try |
| `metric` | character | `"silhouette"` | `"silhouette"` or `"calinski_harabasz"` |
| `n_neighbors` | integer | `30` | Neighbors for graph-based methods |
| `seed` | integer | `42` | Random seed |

**Returns:**

```r
list(
  best_clusters   = integer vector,
  best_param      = numeric,
  scores_df       = data.frame(resolution, n_clusters, score, metric),
  all_clusterings = named list of cluster vectors
)
```

**Input from:** `build_combined_expr_matrix()`.
**Output to:** `clusters_to_gating_csv()`, `validate_clusters()`.

### `consensus_clusters()`

```r
consensus_clusters(cluster_list, method = "majority",
                   n_consensus_clusters = NULL)
```

Builds consensus clusters from multiple clustering runs. The "majority" method builds a co-association matrix and applies hierarchical clustering. The "voting" method uses cell-level majority vote.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_list` | named list | *required* | List of integer vectors (one per run, all same length) |
| `method` | character | `"majority"` | `"majority"` (co-association + hclust) or `"voting"` (cell-level majority) |
| `n_consensus_clusters` | integer | `NULL` | Target clusters; `NULL` = median across runs |

**Returns:** Integer vector of consensus cluster assignments.

**Input from:** `run_resolution_sweep()`, `run_leiden()`, `run_flowsom()`, `run_phenograph()`.
**Output to:** `clusters_to_gating_csv()`, `validate_clusters()`.

### `consensus_via_sampling()` (internal)

Used internally for large datasets (>10,000 cells). Computes co-association on landmark cells and assigns remaining cells by co-occurrence.

### `compute_silhouette()`

```r
compute_silhouette(dist_obj, cluster_ids)
```

**Returns:** Mean silhouette width (numeric scalar).

### `compute_calinski_harabasz()`

```r
compute_calinski_harabasz(expr_matrix, cluster_ids)
```

**Returns:** Calinski-Harabasz index (numeric scalar, higher = better).

---

## cluster_validation.R

Internal and external validation metrics for evaluating cluster quality.

### `validate_clusters()`

```r
validate_clusters(expr_matrix, cluster_ids, true_labels = NULL,
                  metrics = c("silhouette", "calinski_harabasz"),
                  sample_size = 5000, seed = 42)
```

Computes validation metrics, subsampling large datasets for efficiency. When `true_labels` are provided, also computes the adjusted Rand index.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers matrix |
| `cluster_ids` | integer vector | *required* | Cluster assignments |
| `true_labels` | vector | `NULL` | Ground-truth labels for ARI |
| `metrics` | character vector | `c("silhouette", "calinski_harabasz")` | Metrics to compute. Options: `"silhouette"`, `"calinski_harabasz"`, `"davies_bouldin"` |
| `sample_size` | integer | `5000` | Maximum cells for metric computation |
| `seed` | integer | `42` | Random seed |

**Returns:** `data.frame` with columns `metric`, `value`, `description`.

**Input from:** Any clustering function.
**Output to:** `run_resolution_sweep()`, reporting.

**Example:**

```r
source("cluster_validation.R")
metrics <- validate_clusters(
  expr_matrix = expr,
  cluster_ids = clusters,
  metrics     = c("silhouette", "calinski_harabasz", "davies_bouldin")
)
print(metrics)
#              metric    value                                    description
# 1       n_clusters    5.000                           Number of clusters
# 2          n_cells 5000.000                         Total number of cells
# 3       silhouette    0.342    Mean silhouette width (-1 to 1; higher = ...)
# 4 calinski_harabasz 234.567    Calinski-Harabasz index (higher = better ...)
# 5    davies_bouldin    1.234    Davies-Bouldin index (lower = better ...)
```

### `compute_davies_bouldin()` (internal)

```r
compute_davies_bouldin(expr_matrix, cluster_ids)
```

**Returns:** Davies-Bouldin index (numeric scalar, lower = better).

### `compute_ari_manual()` (internal)

```r
compute_ari_manual(labels1, labels2)
```

Fallback ARI computation when `mclust` is unavailable.

**Returns:** Adjusted Rand Index (numeric scalar, -1 to 1).

---

## cell_annotation.R

Annotates clusters with cell type labels based on marker expression profiles.

### `annotate_clusters()`

```r
annotate_clusters(expr_matrix, cluster_ids, markers = NULL,
                  reference = "auto", method = "threshold")
```

Computes per-cluster median expression and assigns cell type annotations using one of three methods.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers matrix |
| `cluster_ids` | integer vector | *required* | Cluster assignments |
| `markers` | character vector | `NULL` | Marker names; `NULL` uses `colnames(expr_matrix)` |
| `reference` | varies | `"auto"` | `"auto"` for threshold phenotyping, or a reference matrix for correlation |
| `method` | character | `"threshold"` | `"threshold"`, `"correlation"`, or `"manual"` |

**Methods:**

| Method | Description |
|--------|-------------|
| `"threshold"` | Markers above global median are positive (+), others negative (-). Builds phenotype string. |
| `"correlation"` | Cosine similarity to reference cell type profiles. Requires reference matrix. |
| `"manual"` | Returns median expression per cluster for user labeling. |

**Returns:** `data.frame` with columns `cluster`, `annotation`, `confidence`, plus `median_<marker>` columns.

**Input from:** Any clustering function.
**Output to:** `clusters_to_gating_csv()`, abundance analysis.

**Example:**

```r
source("cell_annotation.R")

annotations <- annotate_clusters(
  expr_matrix = combined_expr,
  cluster_ids = leiden_clusters,
  method      = "threshold"
)
print(annotations[, c("cluster", "annotation", "confidence")])
#   cluster       annotation confidence
# 1       1 CD3+CD4+CD8-CD19-     0.75
# 2       2 CD3+CD4-CD8+CD19-     0.80
```

### `compute_cluster_medians()`

```r
compute_cluster_medians(expr_matrix, cluster_ids, markers, unique_clusters)
```

**Returns:** Matrix (clusters x markers) of median expression values.

### `cosine_similarity()`

```r
cosine_similarity(a, b)
```

**Returns:** Cosine similarity between two numeric vectors (-1 to 1).

---

## flowjo_import.R

Imports gating hierarchies from FlowJo workspace files.

### `import_flowjo_gates()`

```r
import_flowjo_gates(wsp_path, fcs_dir = dirname(wsp_path),
                    target_populations = NULL)
```

Parses a `.wsp` file, extracts per-cell population assignments as binary columns.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wsp_path` | character | *required* | Path to `.wsp` file |
| `fcs_dir` | character | `dirname(wsp_path)` | Directory with FCS files |
| `target_populations` | character vector | `NULL` | Populations to extract; `NULL` = all |

**Returns:**

```r
list(
  gating_set   = <GatingSet>,
  assignments  = data.frame(cell_id, sample, + binary population columns),
  populations  = character vector,
  sample_names = character vector
)
```

**Input from:** `.wsp` + `.fcs` files.
**Output to:** `get_cluster_assignments()` (cluster_assignment.R).

**Packages:** `CytoML`, `flowWorkspace`

---

## cluster_assignment.R

Unified interface for obtaining cluster assignments from multiple sources.

### `get_cluster_assignments()`

```r
get_cluster_assignments(source = c("leiden", "flowjo", "csv", "flowsom"), ...)
```

Dispatches to the appropriate loader function based on `source`.

| Source | Internal Function | Key Arguments |
|--------|-------------------|---------------|
| `"leiden"` | `get_clusters_from_leiden(csv_path)` | Path to `leiden_clusters.csv` |
| `"flowjo"` | `get_clusters_from_flowjo(wsp_path, fcs_dir)` | Workspace + FCS dir |
| `"csv"` | `get_clusters_from_csv(csv_path, cluster_col)` | CSV with cluster labels |
| `"flowsom"` | `get_clusters_from_flowsom(fsom_object)` | FlowSOM result object |

**Returns:** All dispatch functions return a standardized `data.frame`:

| Column | Type | Description |
|--------|------|-------------|
| `cell_id` | integer | Sequential cell identifier |
| `cluster` | integer/character | Cluster ID |
| `cluster_name` | character | Human-readable label |

### `get_clusters_from_leiden()`

```r
get_clusters_from_leiden(csv_path)
```

Reads `leiden_clusters.csv` and returns standardized cluster assignments.

### `get_clusters_from_flowjo()`

```r
get_clusters_from_flowjo(wsp_path, fcs_dir = dirname(wsp_path))
```

Loads FlowJo gates and converts binary population columns to single cluster assignments.

### `get_clusters_from_csv()`

```r
get_clusters_from_csv(csv_path, cluster_col = "cluster")
```

Reads a generic CSV with a cluster column. Auto-detects cluster-like columns if the specified column is missing.

### `get_clusters_from_flowsom()`

```r
get_clusters_from_flowsom(fsom_object)
```

Extracts meta-cluster assignments from a FlowSOM result object.

**Example:**

```r
source("cluster_assignment.R")

# From Leiden output
assignments <- get_cluster_assignments("leiden",
  csv_path = "/output/leiden_clusters.csv"
)

# From FlowJo workspace
assignments <- get_cluster_assignments("flowjo",
  wsp_path = "/data/experiment.wsp",
  fcs_dir  = "/data/fcs/"
)

# From a generic CSV
assignments <- get_cluster_assignments("csv",
  csv_path    = "/data/labels.csv",
  cluster_col = "population"
)
```

**Dependencies:** `io_utils`, `cluster_utils`
