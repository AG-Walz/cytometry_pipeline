[< Back to Documentation Index](../README.md)

# src/clustering/ -- Leiden Clustering, FlowJo Import, Unified Assignment

3 scripts for assigning cell clusters via Leiden community detection, importing FlowJo gates, and providing a unified interface to multiple clustering sources.

---

## leiden.R

Runs Leiden community detection on expression matrices with optional UMAP embedding.

### Functions

#### run_leiden

```r
run_leiden(expr_matrix, n_neighbors = 30, resolution = 1.0, seed = 42)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `n_neighbors` | integer | `30` | Number of nearest neighbors for KNN graph |
| `resolution` | numeric | `1.0` | Leiden resolution parameter (higher = more clusters) |
| `seed` | integer | `42` | Random seed for reproducibility |

**Returns**: Integer vector of cluster assignments (one per cell).

#### run_leiden_with_umap

```r
run_leiden_with_umap(expr_matrix, sample_labels = NULL, output_dir = NULL,
                     n_neighbors = 30, resolution = 1.0, seed = 42)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell (for coloring) |
| `output_dir` | character | `NULL` | Directory to save output CSV |
| `n_neighbors` | integer | `30` | Nearest neighbors for both KNN and UMAP |
| `resolution` | numeric | `1.0` | Leiden resolution parameter |
| `seed` | integer | `42` | Random seed |

**Returns**: `data.frame` with columns `UMAP1`, `UMAP2`, `cluster`, `sample`.

**Output file**: `leiden_clusters.csv` (when `output_dir` is provided).

### Algorithm

1. Build KNN graph using `uwot` (nearest neighbor search)
2. Construct sparse adjacency matrix (`Matrix` package)
3. Run Leiden community detection via `leiden::leiden()`

### Packages

`leiden`, `igraph`, `uwot`, `Matrix`

### src/R_shared Dependencies

`io_utils`, `stats_utils` (`run_umap`)

---

## flowjo_import.R

Imports gating hierarchies from FlowJo workspace files and extracts per-cell population assignments.

### Main Function

```r
import_flowjo_gates(wsp_path, fcs_dir = dirname(wsp_path),
                    target_populations = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wsp_path` | character | *required* | Path to the FlowJo `.wsp` file |
| `fcs_dir` | character | `dirname(wsp_path)` | Directory containing referenced FCS files |
| `target_populations` | character vector | `NULL` | Specific populations to extract; `NULL` extracts all |

### Output

```r
list(
  gating_set   = GatingSet,    # flowWorkspace GatingSet object
  assignments  = data.frame(), # Columns: cell_id, sample, + binary population columns
  populations  = character(),  # Names of extracted populations
  sample_names = character()   # Sample identifiers
)
```

The `assignments` data.frame contains one row per cell. Each population is represented as a binary column (1 = cell belongs to population, 0 = does not).

### Internal Workflow

1. `CytoML::open_flowjo_xml()` -- parse workspace XML
2. `flowWorkspace::flowjo_to_gatingset()` -- build GatingSet with linked FCS data
3. `flowWorkspace::gh_pop_get_indices()` -- extract boolean gate membership per cell

### Packages

`CytoML`, `flowWorkspace`

### src/R_shared Dependencies

`io_utils`, `cluster_utils` (`clean_cluster_name`)

---

## cluster_assignment.R

Unified interface for obtaining cluster assignments from multiple sources (Leiden, FlowJo, CSV, FlowSOM).

### Main Function

```r
get_cluster_assignments(source = c("leiden", "flowjo", "csv", "flowsom"), ...)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `source` | character | *required* | Clustering source to dispatch to |
| `...` | | | Source-specific arguments (see dispatch table) |

### Dispatch Table

| Source | Internal Function | Key Arguments |
|--------|-------------------|---------------|
| `"leiden"` | `get_clusters_from_leiden(csv_path)` | Path to `leiden_clusters.csv` |
| `"flowjo"` | `get_clusters_from_flowjo(wsp_path, fcs_dir)` | Workspace path + FCS directory |
| `"csv"` | `get_clusters_from_csv(csv_path, cluster_col = "cluster")` | Path to CSV with cluster labels |
| `"flowsom"` | `get_clusters_from_flowsom(fsom_object)` | FlowSOM result object |

### Output

All dispatch functions return a standardized `data.frame`:

| Column | Type | Description |
|--------|------|-------------|
| `cell_id` | character | Unique cell identifier |
| `cluster` | integer | Numeric cluster ID |
| `cluster_name` | character | Human-readable cluster label |

### CSV Auto-Detection

The CSV loader (`get_clusters_from_csv`) auto-detects the cluster column by searching for column names matching the pattern:

```
cluster|population|gate
```

If multiple columns match, the first match is used. The `cluster_col` parameter overrides auto-detection.

### Example

```r
# From Leiden clustering output
assignments <- get_cluster_assignments("leiden", csv_path = "output/leiden_clusters.csv")

# From FlowJo workspace
assignments <- get_cluster_assignments("flowjo",
  wsp_path = "data/experiment.wsp",
  fcs_dir  = "data/fcs/"
)

# From a generic CSV file
assignments <- get_cluster_assignments("csv",
  csv_path    = "data/labels.csv",
  cluster_col = "population"
)
```

### src/R_shared Dependencies

`io_utils`, `cluster_utils`
