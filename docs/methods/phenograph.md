# PhenoGraph: KNN Graph Clustering for Single-Cell Phenotyping

Graph-based clustering that builds a Jaccard-weighted shared nearest-neighbor graph and applies Louvain community detection to identify cell populations in high-dimensional cytometry data.

---

## When to Use

- You are performing exploratory phenotyping and do not know the number of cell populations in advance.
- You want cluster definitions driven by local neighborhood structure rather than global geometry.
- You need sensitivity to rare cell populations that comprise <1% of the total.
- Your data has been appropriately transformed (arcsinh or logicle) and optionally batch-corrected.

PhenoGraph is the natural choice when you want the algorithm to determine the number of clusters from the data. If you need guaranteed cluster connectivity (no disconnected communities), use [Leiden clustering](leiden-clustering.md), which replaces PhenoGraph's Louvain step with the Leiden algorithm.

---

## How It Works

PhenoGraph operates in three stages: k-nearest-neighbor (KNN) graph construction, shared nearest-neighbor (SNN) weighting, and Louvain community detection.

In the first stage, each cell is connected to its k nearest neighbors in marker expression space using Euclidean distance. The pipeline uses approximate nearest-neighbor search via uwot (Annoy algorithm) for O(n log n) construction. This produces a directed graph where each cell points to its k most similar cells.

In the second stage, the directed KNN graph is converted into an undirected shared nearest-neighbor (SNN) graph with Jaccard-weighted edges. For each pair of cells that share at least one neighbor, the edge weight is computed as the Jaccard index: the number of shared neighbors divided by the size of the union of their neighbor sets. This weighting emphasizes connections between cells embedded in the same local neighborhood and downweights connections between cells at the boundary of different populations. The Jaccard weighting is a key innovation of PhenoGraph that improves community detection in datasets with variable local density.

In the third stage, the Louvain algorithm partitions the SNN graph into communities by iteratively optimizing modularity. Nodes are moved between communities to maximize the density of within-community edges relative to a random null model. The algorithm runs until convergence, producing a flat partition where each cell belongs to exactly one cluster. The number of clusters emerges automatically from the graph structure and is not specified by the user.

[Levine et al. (2015)](https://doi.org/10.1016/j.cell.2015.05.047) demonstrated that PhenoGraph identified 25+ immunophenotypically distinct cell populations in healthy human bone marrow, including rare progenitor subsets that manual gating missed entirely.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expr_matrix` | (required) | Numeric matrix, rows = cells, columns = markers |
| `k` | 30 | Number of nearest neighbors for graph construction |
| `seed` | 42 | Random seed for reproducibility |

**Choosing `k`**: The default of 30 is well-validated for cytometry data. Values of 10--15 produce more granular communities (more clusters). Values of 50--100 produce broader communities (fewer clusters). Levine et al. found PhenoGraph results are robust across k = 10 to 50.

---

## Pipeline Integration

**Script**: `src/clustering/phenograph.R`

**Function**: `run_phenograph(expr_matrix, k, seed)`

**Input from**:
- `build_combined_expr_matrix()` in `src/R_shared/flowframe_utils.R`
- `normalize_expression()` in `src/preprocessing/normalization.R`

**Output to**:
- `clusters_to_gating_csv()` in `src/R_shared/cluster_utils.R` -- exports per-sample cluster CSVs
- `validate_clusters()` in `src/clustering/cluster_validation.R`
- `annotate_clusters()` in `src/clustering/cell_annotation.R`

**Returns**: Integer vector of cluster assignments (one per cell).

**Used in**: Pipelines A and E as an alternative to Leiden or FlowSOM clustering. The implementation builds a KNN graph using uwot, constructs a Jaccard-weighted SNN adjacency matrix, converts to an igraph object, and runs `igraph::cluster_louvain()`.

---

## PhenoGraph vs Leiden

The pipeline's Leiden implementation (`src/clustering/leiden.R`) follows the same KNN-graph architecture as PhenoGraph but replaces Louvain with the Leiden algorithm. The key differences are:

| Aspect | PhenoGraph (Louvain) | Leiden |
|--------|---------------------|--------|
| Connectivity guarantee | No -- can produce disconnected communities | Yes -- all communities guaranteed connected |
| Partition quality | Good | Higher (refinement phase) |
| Convergence speed | Good | Faster on large graphs |
| Edge weighting | Jaccard-weighted SNN | Unweighted (extendable to Jaccard) |

For new analyses, Leiden is generally preferred. PhenoGraph remains useful for comparison with published results that used the original Louvain-based method.

---

## References

1. Levine, J. H., Simonds, E. F., Bendall, S. C. et al. (2015). Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. *Cell*, 162(1), 184-197. [doi:10.1016/j.cell.2015.05.047](https://doi.org/10.1016/j.cell.2015.05.047)

2. Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics*, P10008. [doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)

3. Traag, V. A., Waltman, L. & van Eck, N. J. (2019). From Louvain to Leiden: Guaranteeing well-connected communities. *Scientific Reports*, 9, 5233. [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
