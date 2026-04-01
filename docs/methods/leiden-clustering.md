# Leiden Community Detection Clustering

## Summary

Leiden is a graph-based community detection algorithm that partitions cells into phenotypically similar clusters by optimizing a quality function on a k-nearest-neighbor (KNN) graph built from marker expression data. It [improves upon the widely used Louvain algorithm by guaranteeing that all discovered communities are well-connected](https://doi.org/10.1038/s41598-019-41695-z), resolving a fundamental flaw where Louvain can produce arbitrarily disconnected clusters. In this pipeline, Leiden clustering is the primary unsupervised method for identifying cell populations from transformed expression matrices, implemented in `src/clustering/leiden.R`.

---

## The Method

Leiden clustering operates in three stages: KNN graph construction, modularity optimization, and a refinement phase that distinguishes it from earlier approaches.

### Stage 1: KNN Graph Construction

The algorithm begins by representing cells as nodes in a graph. For each cell, the [k nearest neighbors are identified using Euclidean distance in the marker expression space](https://doi.org/10.1016/j.cell.2015.05.047). In this pipeline, nearest-neighbor search is performed by [uwot](https://doi.org/10.48550/arXiv.1802.03426), which uses the Annoy (Approximate Nearest Neighbors Oh Yeah) algorithm for efficient approximate neighbor lookup. The resulting directed KNN graph is then symmetrized into an undirected adjacency matrix: if cell A is a neighbor of cell B or cell B is a neighbor of cell A, an edge exists between them. This symmetrization step ensures that the neighborhood relationship is mutual, producing a graph where edge weights are binary (present or absent).

### Stage 2: Modularity Optimization

The algorithm [optimizes a quality function --- typically the Constant Potts Model (CPM) or modularity --- by iteratively moving nodes between communities](https://doi.org/10.1038/s41598-019-41695-z). The modularity function, originally formalized by [Newman and Girvan (2004)](https://doi.org/10.1103/PhysRevE.69.026113), quantifies the density of edges within communities compared to what would be expected under a random null model. A resolution parameter gamma scales this comparison, controlling the granularity of the detected communities. Like [Louvain (Blondel et al. 2008)](https://doi.org/10.1088/1742-5468/2008/10/P10008), Leiden performs fast local moves of individual nodes to neighboring communities, accepting moves that increase the quality function.

### Stage 3: Refinement Phase

After the initial partition is obtained, [Leiden introduces a refinement phase that Louvain lacks](https://doi.org/10.1038/s41598-019-41695-z). During refinement, the algorithm considers merging nodes within each community into refined subcommunities, starting from a singleton partition within each community. This ensures that every community in the final partition is a subset of the community found during the initial move phase, and critically, it [guarantees that all communities remain internally connected](https://doi.org/10.1038/s41598-019-41695-z). After refinement, the algorithm aggregates the network based on the refined partition and repeats the process until convergence.

The combination of fast local moves (from Louvain) with the refinement step (unique to Leiden) produces partitions that are both [higher quality and provably well-connected](https://doi.org/10.1038/s41598-019-41695-z).

---

## Why Leiden Over Louvain

The [Louvain algorithm (Blondel et al. 2008)](https://doi.org/10.1088/1742-5468/2008/10/P10008) was the dominant community detection method for over a decade, but it has a well-documented structural flaw.

### The Disconnected Communities Problem

[Traag, Waltman, and van Eck (2019)](https://doi.org/10.1038/s41598-019-41695-z) demonstrated that Louvain can yield communities that are internally disconnected --- meaning a community may contain subsets of nodes with no path between them within the community. This occurs because Louvain's greedy node-moving phase can strand nodes in communities where they are no longer connected to the rest of the community after subsequent moves. In practice, [this means Louvain can report a "cluster" that actually consists of two or more separate cell groups with no phenotypic similarity](https://doi.org/10.1038/s41598-019-41695-z), which is biologically misleading.

### Leiden's Guarantees

Leiden addresses this by:

1. **Guaranteed connectivity**: [Every community in a Leiden partition is guaranteed to be connected as a subgraph](https://doi.org/10.1038/s41598-019-41695-z). There is no possibility of disconnected clusters.
2. **Higher quality partitions**: The refinement phase [explores a broader solution space, consistently finding partitions with higher quality function values than Louvain](https://doi.org/10.1038/s41598-019-41695-z).
3. **Faster convergence**: Despite the additional refinement step, [Leiden converges faster than Louvain on large networks because it visits fewer nodes per iteration](https://doi.org/10.1038/s41598-019-41695-z), using a fast local move strategy that only considers nodes in a queue rather than scanning all nodes.

For cytometry data, where a "cluster" is interpreted as a biologically meaningful cell population, the connectivity guarantee is essential. A disconnected community could conflate unrelated cell types, corrupting downstream differential abundance and enrichment statistics.

---

## Why Graph-Based Clustering for Cytometry

Flow and mass cytometry data present specific challenges that make graph-based clustering methods particularly well-suited compared to classical approaches.

### Limitations of Traditional Methods

- **k-means** assumes spherical, equally-sized clusters and requires pre-specifying the number of clusters. [Cytometry data typically contains populations of highly variable size and irregular shape](https://doi.org/10.1016/j.cell.2015.05.047), making k-means a poor fit. Rare populations (< 1% of cells) are routinely absorbed into neighboring large clusters.
- **Hierarchical clustering** scales poorly with cell count (O(n^2) memory for distance matrices) and requires choosing a cut height, which implicitly assumes a single resolution is appropriate for all populations.
- **Model-based clustering** (e.g., Gaussian mixture models) assumes parametric density shapes that [do not hold for the complex, multimodal distributions observed in high-dimensional cytometry data](https://doi.org/10.1002/cyto.a.22625).

### Graph-Based Advantages

[PhenoGraph (Levine et al. 2015)](https://doi.org/10.1016/j.cell.2015.05.047) demonstrated that graph-based community detection is particularly effective for cytometry because:

1. **No assumption about cluster shape**: Communities are defined by local connectivity, not geometric shape. Elongated, crescent-shaped, or irregular populations are captured naturally.
2. **Automatic cluster number determination**: The resolution parameter controls granularity without requiring the user to specify an exact number of clusters in advance.
3. **Scalability**: KNN graph construction is O(n log n) with approximate methods, and community detection scales near-linearly with graph size. This makes it feasible for datasets of [hundreds of thousands to millions of cells](https://doi.org/10.1016/j.cell.2015.05.047).
4. **Sensitivity to rare populations**: By operating on local neighborhoods, graph-based methods can detect small but distinct clusters that global methods miss. [Levine et al. showed PhenoGraph identified rare progenitor populations that manual gating overlooked](https://doi.org/10.1016/j.cell.2015.05.047).

[FlowSOM (Van Gassen et al. 2015)](https://doi.org/10.1002/cyto.a.22625) takes a different graph-based approach using self-organizing maps, achieving extremely fast runtimes while still leveraging graph topology for metaclustering. Both PhenoGraph and FlowSOM [outperformed traditional methods in systematic benchmarks](https://doi.org/10.1002/cyto.a.22764), confirming the advantage of graph-aware approaches for cytometry.

---

## The Resolution Parameter

The resolution parameter (gamma) is the single most important user-facing control for Leiden clustering. It determines the trade-off between splitting and merging communities.

### How Resolution Works

[In the modularity quality function, the resolution parameter scales the null-model term](https://doi.org/10.1038/s41598-019-41695-z). Higher resolution values penalize large communities more heavily, favoring many small communities. Lower values favor fewer, larger communities. The relationship is monotonic: [increasing resolution never decreases (and usually increases) the number of detected communities](https://doi.org/10.1038/s41598-019-41695-z).

### Guidelines for Cytometry Data

| Resolution Range | Typical Use | Expected Outcome |
|------------------|-------------|------------------|
| 0.1 -- 0.4 | Very coarse grouping | 3--8 broad compartments (e.g., T cells, B cells, myeloid) |
| **0.5 -- 1.0** | **Major populations** | **8--20 canonical immune populations** |
| 1.5 -- 3.0 | Subpopulations | 20--50+ subsets (e.g., Th1, Th2, Th17 within CD4+ T cells) |
| 5.0+ | Fine-grained | Overclustering for consensus approaches or subclustering workflows |

These ranges assume a typical immunophenotyping panel with 10--30 markers and 10,000--500,000 cells. Panels with fewer markers or fewer cells may need lower resolution; deeply phenotyped panels (e.g., CyTOF with 40+ markers) may tolerate higher resolution.

### Validating Your Resolution Choice

1. **Marker heatmaps**: Plot median marker expression per cluster as a heatmap. If two clusters have nearly identical profiles, the resolution may be too high. If one cluster shows a bimodal distribution for a key marker, it may be too low.
2. **Cluster size distribution**: Check for extremely small clusters (< 20 cells). While rare populations are real, a proliferation of tiny clusters often indicates over-resolution.
3. **Biological interpretability**: Clusters should correspond to known or plausible cell types. Domain expertise is the final arbiter.
4. **Stability analysis**: Run clustering at the same resolution with different random seeds. Clusters that appear consistently are robust; clusters that fragment or merge across runs are unstable and may indicate a resolution boundary.
5. **Multi-resolution comparison**: [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) recommended running clustering at multiple resolutions and comparing results, allowing the analyst to select the resolution that best matches the biological question.

---

## KNN Graph Construction

The k-nearest-neighbor graph is the foundation on which Leiden clustering operates. Its construction parameters directly affect the resulting partition.

### The n_neighbors Parameter (k)

The `n_neighbors` parameter specifies how many nearest neighbors each cell is connected to in the initial directed graph. After symmetrization, the effective degree of each node is between k and 2k.

| n_neighbors | Effect | When to Use |
|-------------|--------|-------------|
| 10--15 | Tight local neighborhoods, more granular communities | Small datasets (< 10,000 cells), fine-grained phenotyping |
| **20--30** | **Balanced local/global structure** | **Most immunophenotyping experiments** |
| 50--100 | Broad neighborhoods, smoother community boundaries | Very large datasets (> 500,000 cells), noisy data |

The [PhenoGraph paper (Levine et al. 2015)](https://doi.org/10.1016/j.cell.2015.05.047) used k = 30 as a default for cytometry and found results to be [robust across a range of k values from 10 to 50](https://doi.org/10.1016/j.cell.2015.05.047). This pipeline follows the same default.

### Distance Metric

This pipeline uses Euclidean distance in the transformed marker expression space, which is the [standard choice for cytometry data after appropriate transformation](https://doi.org/10.1016/j.cell.2015.05.047). The Euclidean metric assumes that all marker dimensions contribute equally, which is a reasonable assumption after transformation and optional scaling.

Alternative metrics (cosine, correlation) may be appropriate for specific applications such as RNA-velocity-derived features, but these are not standard for protein-level cytometry data.

### uwot Implementation

The pipeline uses the [uwot R package](https://doi.org/10.48550/arXiv.1802.03426) for nearest-neighbor computation. uwot implements the Annoy algorithm for approximate nearest-neighbor search, which [partitions the space using random projection trees and achieves O(n log n) construction time](https://doi.org/10.48550/arXiv.1802.03426). The approximation introduces a small error rate (typically < 5% of neighbors differ from the exact solution), which [has negligible effect on community detection quality](https://doi.org/10.1016/j.cell.2015.05.047) because the graph structure is robust to minor perturbations.

In this pipeline, `run_leiden()` calls `uwot::umap()` with `ret_nn = TRUE` to simultaneously compute UMAP coordinates and retrieve the nearest-neighbor index, avoiding redundant computation when both UMAP and Leiden clustering are needed.

### Adjacency Matrix Construction

After nearest-neighbor search, the pipeline constructs a sparse adjacency matrix using `Matrix::sparseMatrix()`:

```r
adj_matrix <- Matrix::sparseMatrix(
  i = rep(1:n_cells, each = n_neighbors),
  j = as.vector(nn$nn$euclidean$idx),
  x = 1,
  dims = c(n_cells, n_cells)
)
# Symmetrize: if A->B or B->A, then A<->B
adj_matrix <- adj_matrix + Matrix::t(adj_matrix)
adj_matrix@x[adj_matrix@x > 1] <- 1
```

The symmetrization step follows the [shared nearest neighbor (SNN) convention used in PhenoGraph](https://doi.org/10.1016/j.cell.2015.05.047): an edge exists between two cells if either cell considers the other a neighbor. Edge weights are binary (unweighted). Weighted alternatives (e.g., Jaccard-weighted SNN as in [PhenoGraph](https://doi.org/10.1016/j.cell.2015.05.047)) can improve performance for some datasets but are not currently implemented in this pipeline.

---

## When Your Data Calls For It

Leiden clustering is the right choice when:

### Unsupervised Discovery Without Predefined Gates

You have a new panel, a new tissue, or an exploratory study where the expected populations are not fully known. Manual gating requires prior knowledge of marker co-expression patterns, which may not be available. Leiden clustering [discovers populations directly from the data without requiring predefined gate boundaries](https://doi.org/10.1016/j.cell.2015.05.047), enabling the detection of unexpected or novel cell types.

### Data-Driven Population Identification

You want cluster definitions to emerge from the expression data itself rather than from analyst-imposed hierarchical gating schemes. This is particularly important when [manual gating introduces inter-operator variability](https://doi.org/10.1002/cyto.a.22625), a well-documented problem in clinical cytometry. Algorithmic clustering provides reproducible, operator-independent population definitions.

### Transformed and Optionally Batch-Corrected Input

Your data has been appropriately transformed (arcsinh for CyTOF, logicle for fluorescence) and, if multi-batch, batch-corrected. Leiden clustering operates on the expression space directly, so [the quality of the input transformation directly determines cluster quality](https://doi.org/10.12688/f1000research.11622.3). See [Pre-processing Requirements](#pre-processing-requirements) below.

### Scalable to Hundreds of Thousands of Cells

Your dataset contains tens of thousands to hundreds of thousands of cells. Leiden clustering, with approximate nearest-neighbor search, [scales near-linearly with cell count](https://doi.org/10.1038/s41598-019-41695-z), making it practical for large pooled experiments. [Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) demonstrated that graph-based methods maintain accuracy at high cell counts where model-based approaches become computationally prohibitive.

### Multi-Sample Comparative Studies

You are comparing cell populations across conditions, timepoints, or patient groups. After pooled clustering, each sample's cells retain their cluster labels, enabling [differential abundance testing with appropriate statistical models](https://doi.org/10.12688/f1000research.11622.3) (e.g., GLMM, negative binomial regression). This pipeline supports exactly this workflow via `clusters_to_gating_csv()`, which exports per-sample cluster assignments for downstream analysis.

---

## When to Choose Something Else

Leiden clustering is not always the best tool. Consider alternatives when:

### You Have Predefined Gates

If your lab has established, validated gating strategies and you need to apply them consistently, use **FlowJo import** (`src/clustering/flowjo_import.R`). Leiden will find populations, but they will not correspond one-to-one with manually defined gates. For regulatory or clinical reporting where specific gate definitions are required, manual or template-based gating is more appropriate.

### You Need Deterministic Results

Leiden clustering involves stochastic elements: the Annoy approximate nearest-neighbor search and the random order of node moves during optimization. Results will vary slightly between runs. To mitigate this, **always set a random seed** (the `seed` parameter in `run_leiden()`). With a fixed seed, results are reproducible within the same software environment. If even seed-dependent variation is unacceptable, consider [FlowSOM](https://doi.org/10.1002/cyto.a.22625), which produces more deterministic results when initialized with a fixed seed due to its SOM-based architecture.

### Very Few Cells (< 100)

Graph-based clustering requires sufficient cells to form meaningful neighborhoods. With fewer than ~100 cells, the KNN graph is too sparse for reliable community detection. For very small samples, consider manual gating or simple thresholding.

### High-Dimensional Noise Reduction Needed

If your panel has many markers (> 30) and you suspect that some are adding noise rather than signal, consider [dimensionality reduction with PCA](pca.md) before clustering. Reducing to the top 10--20 principal components can improve cluster quality by removing noisy dimensions. This is particularly relevant for CyTOF panels where some channels may have low signal-to-noise ratios.

### You Need Hierarchical Relationships

Leiden produces a flat partition. If you need to know that Th1 cells are a subset of CD4+ T cells, which are a subset of T cells, you need either hierarchical gating or a multi-resolution clustering approach. One strategy is to run Leiden at a low resolution first, then subcluster specific populations at higher resolution.

---

## Pre-processing Requirements

Leiden clustering operates on numeric expression matrices. The quality of the input data directly determines the quality of the clustering. [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) emphasized that appropriate preprocessing is essential before any unsupervised cytometry analysis.

### Transformation (Required)

Raw fluorescence and mass cytometry intensities span several orders of magnitude and are heavily right-skewed. Clustering on untransformed data will be dominated by the most intense markers and will fail to resolve populations at lower expression levels.

- **CyTOF / Mass Cytometry**: Apply [arcsinh transformation](arcsinh-transform.md) with cofactor 5 (standard) or cofactor 150 for count-based data. The arcsinh function compresses high values while preserving differences at low expression, producing [approximately normally distributed marker intensities](https://doi.org/10.1016/j.cell.2015.05.047).
- **Fluorescence Cytometry**: Apply [logicle (biexponential) transformation](logicle-transform.md), which [handles the negative values common in compensated fluorescence data while providing linear behavior near zero and logarithmic compression at high values](https://doi.org/10.1002/cyto.a.22625).

### Compensation (Fluorescence Only)

Fluorescence data must be compensated for spectral overlap before transformation. Uncompensated data will produce spurious correlations between markers, leading to artificial clusters. Use the compensation matrix embedded in the FCS file or provided by the instrument.

### Batch Correction (Multi-Batch Only)

When combining data from multiple acquisition batches, systematic technical variation can [create batch-specific clusters that do not reflect biological differences](https://doi.org/10.12688/f1000research.11622.3). Apply [batch correction](batch-correction.md) (e.g., via `src/harmonization/batch_correction.R`) before clustering. The pipeline's harmonization step (`src/harmonization/panel_homogenization.R`) addresses panel differences; batch correction addresses intensity shifts within a shared panel.

### Marker Selection

Cluster on phenotypic markers that define cell identity, not functional markers (e.g., cytokines, activation markers) that vary with stimulation state. Including functional markers will [split biologically identical cell types into activation-state-dependent subclusters](https://doi.org/10.12688/f1000research.11622.3), obscuring population structure.

### Downsampling (Optional)

For very large datasets (> 500,000 cells), downsampling before clustering can reduce runtime without substantially affecting population identification. The pipeline's `downsample_flowframes()` function (`src/preprocessing/downsampling.R`) supports random and density-dependent downsampling. [Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) showed that downsampling to 10,000--50,000 cells per sample preserves major population structure.

---

## How to Use in the Pipeline

### From Raw FCS Files (Pipeline A)

After transformation, homogenization, and downsampling:

```r
# Source required modules
source("src/clustering/leiden.R")
source("src/R_shared/flowframe_utils.R")

# Build combined expression matrix from downsampled flowFrames
expr_data <- build_combined_expr_matrix(
  downsampled, flow_data$sample_names, params$markers_of_interest
)

# Run Leiden with UMAP
clusters <- run_leiden_with_umap(
  expr_matrix   = expr_data$combined_matrix,
  sample_labels = expr_data$sample_labels,
  output_dir    = file.path(params$output_dir, "03_clustering"),
  n_neighbors   = 30,
  resolution    = 1.0,
  seed          = 42
)

# Export per-sample cluster CSVs for downstream analysis
clusters_to_gating_csv(
  combined_matrix = expr_data$combined_matrix,
  sample_labels   = expr_data$sample_labels,
  cluster_ids     = clusters$cluster,
  output_dir      = file.path(params$output_dir, "cluster_data")
)
```

### From CSV Expression Matrix (Cross-Pipeline)

When starting from a pre-existing expression matrix (e.g., exported from another tool):

```r
source("src/clustering/leiden.R")

# Load expression data
expr_dt <- data.table::fread("path/to/expression_matrix.csv")
marker_cols <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14")

# Run Leiden clustering
clusters <- run_leiden_with_umap(
  expr_matrix = as.matrix(expr_dt[, marker_cols, with = FALSE]),
  resolution  = 1.0,
  n_neighbors = 30,
  seed        = 42
)
```

### Standalone Leiden (Without UMAP)

If you only need cluster labels and already have UMAP coordinates:

```r
source("src/clustering/leiden.R")

cluster_ids <- run_leiden(
  expr_matrix = expr_matrix,
  n_neighbors = 30,
  resolution  = 1.0,
  seed        = 42
)
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expr_matrix` | (required) | Numeric matrix, rows = cells, columns = markers |
| `n_neighbors` | 30 | Number of nearest neighbors for KNN graph |
| `resolution` | 1.0 | Leiden resolution parameter (higher = more clusters) |
| `seed` | 42 | Random seed for reproducibility |
| `sample_labels` | NULL | Character vector mapping each cell to its sample |
| `output_dir` | NULL | If provided, saves `leiden_clusters.csv` to this directory |

### Output Format

`run_leiden_with_umap()` returns a data frame with columns:

| Column | Type | Description |
|--------|------|-------------|
| UMAP1 | numeric | First UMAP coordinate |
| UMAP2 | numeric | Second UMAP coordinate |
| cluster | integer | Leiden cluster assignment (1-indexed) |
| sample | character | Sample label (if `sample_labels` was provided) |

---

## Validating Clusters

Unsupervised clustering requires careful validation to ensure that the discovered populations are biologically meaningful rather than artifacts of parameter choices or data quality issues.

### Marker Expression Heatmaps

Plot the median expression of each marker across all clusters as a heatmap. Each row should represent a cluster and each column a marker. Well-defined clusters will show [distinct expression signatures that correspond to known or plausible cell phenotypes](https://doi.org/10.12688/f1000research.11622.3). The pipeline's visualization module (`src/visualization/heatmaps.R`) generates these heatmaps automatically.

Look for:
- **Distinct profiles**: Each cluster should have a characteristic combination of high/low markers.
- **Redundant clusters**: Two clusters with nearly identical profiles may indicate over-resolution.
- **Mixed clusters**: A cluster expressing markers from unrelated lineages (e.g., both CD3 and CD19) may indicate under-resolution or a data quality issue.

### UMAP Overlays

Overlay cluster labels on [UMAP](umap.md) coordinates to visualize the spatial organization of clusters. Well-separated clusters should occupy distinct regions of the UMAP embedding. However, note that [UMAP distorts distances and densities](https://doi.org/10.48550/arXiv.1802.03426), so clusters that appear close in UMAP space are not necessarily phenotypically similar. Always confirm UMAP observations with marker expression data.

The pipeline computes UMAP coordinates alongside Leiden clustering in `run_leiden_with_umap()`, ensuring that the same neighbor graph underlies both the embedding and the clustering.

### Per-Cluster Marker Distributions

Beyond median heatmaps, examine the full distribution of key markers within each cluster using density plots or violin plots. [Bimodal distributions within a single cluster suggest that the resolution is too low](https://doi.org/10.12688/f1000research.11622.3) and the cluster should be split.

### Cross-Sample Consistency

For multi-sample experiments, check that each cluster is represented across multiple samples. A cluster appearing in only one sample may reflect a sample-specific artifact (e.g., a doublet population, a staining artifact) rather than a genuine biological population.

### Comparison with Manual Gating

When manual gating results are available, compare Leiden clusters against manually defined populations. [This concordance analysis helps calibrate resolution settings and identify populations that the algorithm may split or merge differently than a human analyst would](https://doi.org/10.1002/cyto.a.22764). The pipeline's metadata QC module (`src/metadata_qc/gating_comparison.R`) supports this comparison.

---

## Comparison with Other Clustering Methods

Several clustering algorithms are commonly used in cytometry. The choice depends on dataset size, desired granularity, and computational constraints.

### Leiden vs. FlowSOM

[FlowSOM (Van Gassen et al. 2015)](https://doi.org/10.1002/cyto.a.22625) uses a self-organizing map (SOM) to create an initial high-resolution grid of nodes, then applies hierarchical consensus metaclustering to group nodes into populations.

| Aspect | Leiden | FlowSOM |
|--------|--------|---------|
| Speed | Fast (seconds to minutes) | Very fast (seconds) |
| Scalability | Excellent (> 1M cells) | Excellent (> 1M cells) |
| Cluster number | Determined by resolution | User specifies metacluster count |
| Cluster shape | Arbitrary | Tends toward compact shapes |
| Rare populations | Good sensitivity | [Can miss very rare populations (< 0.1%) due to SOM grid resolution](https://doi.org/10.1002/cyto.a.22764) |
| Reproducibility | Stochastic (seed-dependent) | More stable with fixed seed |

[Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) found that FlowSOM achieved the best overall performance in their benchmark when the correct number of metaclusters was specified, but [performance degraded when the number was misspecified](https://doi.org/10.1002/cyto.a.22764). Leiden avoids this issue by using a continuous resolution parameter rather than a discrete cluster count.

### Leiden vs. PhenoGraph

[PhenoGraph (Levine et al. 2015)](https://doi.org/10.1016/j.cell.2015.05.047) is architecturally similar to the Leiden pipeline: it builds a KNN graph, computes Jaccard-weighted shared nearest neighbors, and applies Louvain community detection. The key difference is that PhenoGraph uses Louvain, while this pipeline uses Leiden. Replacing Louvain with Leiden in a PhenoGraph-like workflow [provides the same benefits: guaranteed connectivity, higher-quality partitions, and faster convergence](https://doi.org/10.1038/s41598-019-41695-z).

PhenoGraph additionally uses Jaccard-weighted edges (where edge weight equals the number of shared neighbors divided by the size of the union of their neighbor sets), which can improve community detection in datasets with variable local density. This pipeline uses unweighted edges for simplicity but could be extended to Jaccard weighting.

### Leiden vs. SPADE

SPADE (Spanning-tree Progression Analysis of Density-normalized Events) uses density-dependent downsampling followed by agglomerative clustering and minimum spanning tree construction. SPADE was [designed specifically for revealing developmental trajectories and rare populations](https://doi.org/10.1016/j.cell.2015.05.047), but its aggressive downsampling discards information, and the spanning tree representation [can distort relationships between populations that are not arranged along a single trajectory](https://doi.org/10.1002/cyto.a.22764). [Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) found that SPADE performed poorly in cluster identification benchmarks compared to FlowSOM and PhenoGraph.

### Leiden vs. k-means

k-means is fast and simple but assumes spherical clusters of similar size and requires specifying the exact number of clusters. In [comprehensive cytometry benchmarks, k-means consistently underperformed graph-based methods](https://doi.org/10.1002/cyto.a.22764) because cytometry populations are typically non-spherical and span a wide range of sizes. k-means also lacks the resolution parameter's ability to smoothly control granularity.

### Summary Comparison Table

| Method | Graph-Based | Cluster Shape | Scalability | Rare Populations | Key Citation |
|--------|-------------|---------------|-------------|------------------|--------------|
| **Leiden** | Yes | Arbitrary | Excellent | Good | [Traag et al. 2019](https://doi.org/10.1038/s41598-019-41695-z) |
| FlowSOM | Yes (SOM) | Compact | Excellent | Moderate | [Van Gassen et al. 2015](https://doi.org/10.1002/cyto.a.22625) |
| PhenoGraph | Yes (KNN + Louvain) | Arbitrary | Good | Good | [Levine et al. 2015](https://doi.org/10.1016/j.cell.2015.05.047) |
| SPADE | Tree-based | Compact | Moderate | Good (by design) | [Qiu et al. 2011](https://doi.org/10.1038/nbt.1991) |
| k-means | No | Spherical | Excellent | Poor | --- |
| Hierarchical | No | Arbitrary | Poor (O(n^2)) | Moderate | --- |

---

## Related Documentation

- [UMAP Visualization](umap.md) --- Dimensionality reduction for cluster visualization
- [PCA](pca.md) --- Optional pre-reduction before clustering
- [Arcsinh Transform](arcsinh-transform.md) --- Required preprocessing for CyTOF data
- [Logicle Transform](logicle-transform.md) --- Required preprocessing for fluorescence data
- [Batch Correction](batch-correction.md) --- Required preprocessing for multi-batch data
- [Pipeline A: Raw FCS](../pipelines/) --- Full pipeline using Leiden clustering
- [Composability Guide](../composability-guide.md) --- How to combine pipeline modules
