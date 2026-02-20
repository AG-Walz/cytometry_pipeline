[< Back to Documentation Index](../README.md)

# src/clustering/ — Leiden Clustering, FlowJo Import, Unified Assignment

Three scripts for assigning cell clusters: Leiden community detection on expression data, importing gate memberships from FlowJo workspaces, and a unified interface that dispatches to the appropriate source.

---

## leiden.R

Runs Leiden community detection on expression matrices, with optional UMAP embedding. This is the primary clustering method used by Pipelines A and E.

### Functions

#### run_leiden

```r
run_leiden(expr_matrix, n_neighbors = 30, resolution = 1.0, seed = 42)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `n_neighbors` | integer | `30` | KNN graph neighbors |
| `resolution` | numeric | `1.0` | Leiden resolution (higher = more clusters) |
| `seed` | integer | `42` | Random seed |

**Returns:** Integer vector of cluster assignments (one per cell).

#### run_leiden_with_umap

```r
run_leiden_with_umap(expr_matrix, sample_labels = NULL, output_dir = NULL,
                     n_neighbors = 30, resolution = 1.0, seed = 42)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell (for output) |
| `output_dir` | character | `NULL` | Directory to save `leiden_clusters.csv` |
| `n_neighbors` | integer | `30` | Neighbors for both KNN graph and UMAP |
| `resolution` | numeric | `1.0` | Leiden resolution |
| `seed` | integer | `42` | Random seed |

**Returns:** `data.frame` with columns `UMAP1`, `UMAP2`, `cluster`, `sample`.

### Algorithm

1. Build KNN graph using `uwot` (nearest neighbor search)
2. Construct sparse adjacency matrix (`Matrix` package)
3. Run Leiden community detection via `leiden::leiden()`
4. (If `run_leiden_with_umap`) Compute 2D UMAP embedding via `uwot::umap()`

### Example

```r
source("leiden.R")

# Build combined matrix from downsampled flowFrames
expr_list <- extract_exprmtx_list(downsampled, sample_names, markers)
combined <- do.call(rbind, expr_list)
labels <- rep(names(expr_list), sapply(expr_list, nrow))

result <- run_leiden_with_umap(
  expr_matrix   = combined,
  sample_labels = labels,
  output_dir    = "/output/clustering",
  n_neighbors   = 15,
  resolution    = 0.8,
  seed          = 42
)
```

**Console output:**

```
  Computing UMAP (n_neighbors=15)...
  Building KNN graph (k=15)...
  Running Leiden clustering (resolution=0.8)...
  Found 5 clusters
  Cluster sizes: 1234, 987, 876, 765, 543
  Results saved to: /output/clustering/leiden_clusters.csv
```

### Output File — `leiden_clusters.csv`

```
"UMAP1","UMAP2","cluster","sample"
1.84,0.85,1,"sample_001"
-0.85,-1.51,1,"sample_001"
0.96,-0.36,2,"sample_001"
3.37,-0.17,3,"sample_003"
...
```

| Column | Type | Description |
|--------|------|-------------|
| `UMAP1` | numeric | First UMAP coordinate |
| `UMAP2` | numeric | Second UMAP coordinate |
| `cluster` | integer | Leiden cluster assignment (1, 2, 3, ...) |
| `sample` | character | Sample of origin |

### Resolution Guide

| Resolution | Effect |
|------------|--------|
| 0.1 – 0.3 | Very few, broad clusters |
| 0.5 – 0.8 | Moderate number (good starting point) |
| 1.0 – 1.5 | Many clusters |
| 2.0+ | Fine-grained sub-populations |

### Python Dependencies

The R `leiden` package calls Python's `leidenalg` via `reticulate`. Required:

```bash
pip install leidenalg igraph numpy
```

If reticulate uses the wrong Python: `export RETICULATE_PYTHON=/path/to/python`

### Packages

`leiden`, `igraph`, `uwot`, `Matrix`

### Dependencies

`io_utils`, `stats_utils` (`run_umap`)

### clusters_to_gating_csv

Converts cluster IDs to per-sample CSVs with binary gating columns. Each output CSV contains the original expression columns plus one `Cluster_N (Clustering annotated)` column per cluster, matching the format expected by downstream analysis scripts.

```r
clusters_to_gating_csv(combined_matrix, sample_labels, cluster_ids,
                       output_dir, sep = ";")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `combined_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `sample_labels` | character vector | *required* | Sample label per cell (same length as rows in `combined_matrix`) |
| `cluster_ids` | integer vector | *required* | Cluster assignment per cell |
| `output_dir` | character | *required* | Directory for output CSVs |
| `sep` | character | `";"` | CSV separator |

**Returns:** Character vector of output file paths (invisible).

#### Example

```r
source("leiden.R")

# Full flow: build matrix, cluster, then export gating CSVs
result <- build_combined_expr_matrix(flow_frames, sample_names, markers)
leiden_result <- run_leiden_with_umap(
  expr_matrix   = result$combined_matrix,
  sample_labels = result$sample_labels,
  output_dir    = "/output/clustering"
)

clusters_to_gating_csv(
  combined_matrix = result$combined_matrix,
  sample_labels   = result$sample_labels,
  cluster_ids     = leiden_result$cluster,
  output_dir      = "/output/gating_csvs"
)
# Wrote 3 cluster CSV files to: /output/gating_csvs
```

**Input from:** `run_leiden_with_umap()`. **Output to:** `split_by_cluster()`, `run_differential_abundance()`

---

## flowjo_import.R

Imports gating hierarchies from FlowJo workspace files and extracts per-cell population assignments as binary columns.

### Function

```r
import_flowjo_gates(wsp_path, fcs_dir = dirname(wsp_path),
                    target_populations = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wsp_path` | character | *required* | Path to `.wsp` file |
| `fcs_dir` | character | `dirname(wsp_path)` | Directory with FCS files |
| `target_populations` | character vector | `NULL` | Populations to extract; `NULL` = all |

### Returns

```r
list(
  gating_set   = <GatingSet>,
  assignments  = data.frame(),   # cell_id, sample, + binary population columns
  populations  = character(),
  sample_names = character()
)
```

### Example Assignments Output

```
cell_id;sample;CD3+;CD4+ T cells;CD8+ T cells;B cells
sample_001_cell_1;sample_001;1;1;0;0
sample_001_cell_2;sample_001;1;0;1;0
sample_001_cell_3;sample_001;0;0;0;1
```

### Packages

`CytoML`, `flowWorkspace`

### Dependencies

`io_utils`, `cluster_utils`

---

## cluster_assignment.R

Unified interface for obtaining cluster assignments from multiple sources (Leiden, FlowJo, CSV, FlowSOM).

### Function

```r
get_cluster_assignments(source = c("leiden", "flowjo", "csv", "flowsom"), ...)
```

### Dispatch Table

| Source | Internal Function | Key Arguments |
|--------|-------------------|---------------|
| `"leiden"` | `get_clusters_from_leiden(csv_path)` | Path to `leiden_clusters.csv` |
| `"flowjo"` | `get_clusters_from_flowjo(wsp_path, fcs_dir)` | Workspace + FCS dir |
| `"csv"` | `get_clusters_from_csv(csv_path, cluster_col)` | CSV with cluster labels |
| `"flowsom"` | `get_clusters_from_flowsom(fsom_object)` | FlowSOM result object |

### Returns

All dispatch functions return a standardized `data.frame`:

| Column | Type | Description |
|--------|------|-------------|
| `cell_id` | character | Unique cell identifier |
| `cluster` | integer | Numeric cluster ID |
| `cluster_name` | character | Human-readable label |

### Example

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

### Dependencies

`io_utils`, `cluster_utils`
