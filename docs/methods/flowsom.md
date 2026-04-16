# FlowSOM: Self-Organizing Map Clustering for Cytometry

Fast, scalable clustering that combines self-organizing maps with consensus meta-clustering to identify cell populations in high-dimensional flow and mass cytometry data.

---

## When to Use

- You need fast clustering (seconds, not minutes) on large cytometry datasets.
- You have a reasonable estimate of the expected number of cell populations (meta-cluster count).
- You want a deterministic, reproducible partition with fixed seed initialization.
- Your panel has 10--40 markers and 10,000--1,000,000+ cells.

FlowSOM is the best choice when speed matters and the approximate number of populations is known. If you prefer automatic cluster number determination via a resolution parameter, use [Leiden clustering](leiden-clustering.md) instead.

---

## How It Works

FlowSOM operates in two stages: self-organizing map (SOM) construction followed by hierarchical consensus meta-clustering.

In the first stage, a two-dimensional grid of nodes (typically 10x10 = 100 nodes) is trained on the expression data using the SOM algorithm. Each node learns a prototype vector representing a region of the marker expression space. Cells are assigned to their closest grid node based on Euclidean distance. The SOM preserves topological relationships: similar cell phenotypes map to neighboring nodes on the grid. This initial over-clustering produces many small, homogeneous groups that capture fine-grained structure in the data.

In the second stage, the SOM node prototypes are aggregated into a user-specified number of meta-clusters using consensus clustering (ConsensusClusterPlus). This hierarchical approach clusters the 100 SOM nodes rather than the millions of individual cells, making the meta-clustering step computationally trivial. The consensus procedure runs multiple rounds of hierarchical clustering with random subsampling to identify a stable partition. Each cell inherits its meta-cluster label from its assigned SOM node.

The two-stage architecture is what makes FlowSOM exceptionally fast. By reducing millions of cells to ~100 SOM nodes before meta-clustering, it avoids the O(n^2) scaling of direct hierarchical methods. [Van Gassen et al. (2015)](https://doi.org/10.1002/cyto.a.22625) demonstrated that FlowSOM achieved the best speed-accuracy trade-off in cytometry clustering benchmarks, completing in seconds where other methods required minutes to hours.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expr_matrix` | (required) | Numeric matrix, rows = cells, columns = markers |
| `n_clusters` | 20 | Target number of meta-clusters (final populations) |
| `grid_xdim` | 10 | SOM grid x dimension |
| `grid_ydim` | 10 | SOM grid y dimension |
| `seed` | 42 | Random seed for reproducibility |

**Choosing `n_clusters`**: This is the most critical parameter. Set it to the expected number of cell populations. [Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) showed FlowSOM performs best when `n_clusters` is close to the true number of populations and degrades when misspecified. When uncertain, overestimate slightly and merge redundant clusters post hoc.

**Choosing grid dimensions**: The product `grid_xdim * grid_ydim` determines the number of SOM nodes. A 10x10 grid (100 nodes) works well for most panels. Increase to 15x15 or 20x20 for deeply phenotyped panels or when seeking fine-grained subpopulations.

---

## Pipeline Integration

**Script**: `src/clustering/flowsom.R`

**Function**: `run_flowsom(expr_matrix, n_clusters, grid_xdim, grid_ydim, seed)`

**Input from**:
- `build_combined_expr_matrix()` in `src/R_shared/flowframe_utils.R`
- `normalize_expression()` in `src/preprocessing/normalization.R`

**Output to**:
- `clusters_to_gating_csv()` in `src/R_shared/cluster_utils.R` -- exports per-sample cluster CSVs
- `validate_clusters()` in `src/clustering/cluster_validation.R`
- `annotate_clusters()` in `src/clustering/cell_annotation.R`

**Returns**: Integer vector of cluster assignments (one per cell).

**Used in**: Pipelines A and E as an alternative to Leiden clustering. Also used internally by `run_diffcyt()` in `src/analysis/diffcyt_analysis.R` for stratified differential discovery.

---

## Comparison with Leiden

| Aspect | FlowSOM | Leiden |
|--------|---------|--------|
| Speed | Very fast (seconds) | Fast (seconds to minutes) |
| Cluster number | User-specified `n_clusters` | Controlled by resolution parameter |
| Cluster shape | Tends toward compact shapes | Arbitrary (graph-based) |
| Rare populations | Can miss very rare (<0.1%) populations | Better sensitivity to rare populations |
| Reproducibility | Highly stable with fixed seed | Seed-dependent, slightly variable |
| Best for | Known population counts, speed-critical workflows | Exploratory analysis, unknown populations |

---

## References

1. Van Gassen, S., Callebaut, B., Van Helden, M. J. et al. (2015). FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data. *Cytometry Part A*, 87(7), 636-645. [doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625)

2. Weber, L. M. & Robinson, M. D. (2016). Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data. *Cytometry Part A*, 89(12), 1084-1096. [doi:10.1002/cyto.a.22764](https://doi.org/10.1002/cyto.a.22764)

3. Wilkerson, M. D. & Hayes, D. N. (2010). ConsensusClusterPlus: A class discovery tool with confidence assessments and item tracking. *Bioinformatics*, 26(12), 1572-1573. [doi:10.1093/bioinformatics/btq170](https://doi.org/10.1093/bioinformatics/btq170)
