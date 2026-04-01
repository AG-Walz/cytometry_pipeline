# UMAP (Uniform Manifold Approximation and Projection)

## Summary

UMAP is a nonlinear dimensionality reduction algorithm that projects high-dimensional single-cell expression data into two dimensions for visualization. It [constructs a topological representation of the data as a fuzzy simplicial set and then optimizes a low-dimensional layout](https://doi.org/10.48550/arXiv.1802.03426) that preserves both local neighborhood relationships and broader global structure. In this pipeline, UMAP is the primary visualization method for exploring cell populations, overlaying cluster labels, and comparing marker expression across samples.

---

## The Method

### Topological Foundations

UMAP rests on three assumptions about the data: (1) the data is uniformly distributed on a Riemannian manifold, (2) the Riemannian metric is locally constant, and (3) the manifold is locally connected. From these assumptions, [McInnes, Healy, and Melville (2018)](https://doi.org/10.48550/arXiv.1802.03426) derived an algorithm grounded in algebraic topology and category theory, specifically the framework of fuzzy simplicial sets.

### Fuzzy Simplicial Set Construction

The algorithm begins by computing approximate k-nearest neighbors for every data point in the high-dimensional space. For each point, a local distance metric is calibrated such that each point's local neighborhood covers approximately the same number of neighbors (controlled by `n_neighbors`). These local metrics are then [converted into a fuzzy simplicial set](https://doi.org/10.48550/arXiv.1802.03426) -- a weighted graph where edge weights represent the probability that two points are connected. Crucially, each point's distances are normalized with respect to its nearest neighbor, making the representation adaptive to local density. The local fuzzy sets from each point are then merged via a fuzzy union operation, producing a single global topological representation.

### Low-Dimensional Optimization

UMAP initializes a low-dimensional embedding (typically 2D) using spectral initialization on the fuzzy graph Laplacian. It then refines this layout by minimizing the cross-entropy between the high-dimensional fuzzy simplicial set and an analogous fuzzy set constructed in the low-dimensional space. The low-dimensional set uses a smooth approximation family controlled by `min_dist`, which determines [how tightly points are allowed to pack together](https://doi.org/10.48550/arXiv.1802.03426). The optimization uses stochastic gradient descent with negative sampling, making it efficient for large datasets typical in flow cytometry.

### Relationship to Other Methods

UMAP's mathematical framework generalizes several earlier approaches. [The algorithm can be shown to subsume Laplacian eigenmaps and has connections to t-SNE](https://doi.org/10.48550/arXiv.1802.03426), but it differs in using a more principled cross-entropy loss function rather than the KL-divergence used by [t-SNE (van der Maaten and Hinton, 2008)](https://doi.org/10.1007/s10994-008-5065-1). This difference is partly responsible for UMAP's improved preservation of global structure.

---

## Why UMAP for Cytometry

### Advantages Over t-SNE

[Becht et al. (2019)](https://doi.org/10.1038/nbt.4314) systematically benchmarked UMAP against t-SNE on mass cytometry and single-cell RNA-seq datasets and demonstrated three key advantages:

1. **Global structure preservation.** UMAP [retains the relative arrangement of distinct cell populations](https://doi.org/10.1038/nbt.4314), so that cell types which are phenotypically similar in marker space tend to appear near each other in the 2D embedding. t-SNE, in contrast, often [scrambles the arrangement of well-separated clusters](https://doi.org/10.1007/s10994-008-5065-1) because its heavy-tailed Student-t kernel in low dimensions does not enforce any global ordering.

2. **Scalability.** UMAP [runs significantly faster than t-SNE on datasets of comparable size](https://doi.org/10.1038/nbt.4314). For flow cytometry panels that routinely produce hundreds of thousands to millions of events, this speed advantage is substantial. The approximate nearest-neighbor search (using random projection trees and nearest-neighbor descent) and the stochastic gradient descent optimizer both [scale favorably with dataset size](https://doi.org/10.48550/arXiv.1802.03426).

3. **Reproducibility across runs.** While both algorithms are stochastic, UMAP embeddings [show greater consistency across repeated runs with different random seeds](https://doi.org/10.1038/nbt.4314) than t-SNE, particularly in the relative positioning of clusters. In this pipeline, a fixed seed (default `42`) eliminates run-to-run variation entirely.

### Suitability for Immune Cell Phenotyping

Flow cytometry panels for immune phenotyping typically measure 10 to 40 markers simultaneously. This dimensionality is well within UMAP's effective range. [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) established a standard CyTOF analysis workflow that includes dimensionality reduction as a core step for exploring cell population heterogeneity, validating clustering results, and comparing experimental conditions. UMAP has since become the [dominant visualization method in high-dimensional cytometry workflows](https://doi.org/10.1038/nbt.4314), replacing t-SNE in most modern pipelines.

---

## Parameter Tuning

The pipeline exposes four parameters in `run_umap()` (`src/R_shared/stats_utils.R`). Each controls a different aspect of the embedding.

### `n_neighbors` (default: 15)

This parameter controls the size of the local neighborhood used to construct the fuzzy simplicial set. It directly [determines the balance between local and global structure](https://doi.org/10.48550/arXiv.1802.03426) in the final embedding.

| Value Range | Effect | When to Use |
|-------------|--------|-------------|
| **5--10** | Emphasizes fine-grained local structure. Small, tight subclusters become visible. May fragment large populations. | Exploring rare subpopulations; fine-clustering QC; data with many small, distinct cell types. |
| **15** (default) | Balanced trade-off. Populations are well resolved without excessive fragmentation. | General-purpose cytometry visualization. Good starting point for most panels. |
| **30--50** | Emphasizes global topology. Populations appear as broad, connected regions. Small subpopulations may merge visually. | Comparing broad cell lineages (e.g., lymphoid vs. myeloid); initial exploratory overview; datasets with continuous phenotypic gradients. |

**Visual intuition:** Think of `n_neighbors` as the "zoom level" of the manifold reconstruction. Low values build the manifold from small local patches, capturing fine detail but potentially losing the big picture. High values stitch together larger patches, yielding a smoother global view at the cost of local nuance.

### `min_dist` (default: 0.1)

This parameter controls the [minimum distance between points in the low-dimensional embedding](https://doi.org/10.48550/arXiv.1802.03426). It affects the visual compactness of clusters.

| Value Range | Effect | When to Use |
|-------------|--------|-------------|
| **0.01--0.05** | Very tight, compact clusters with clear separation. Points within a population are densely packed. | Publication figures where visual cluster separation is important; identifying discrete populations. |
| **0.1** (default) | Moderate packing. Clusters are compact but not overly so. Internal structure within clusters remains visible. | General-purpose visualization. Balances compactness with readability. |
| **0.3--0.5** | Loose, spread-out arrangement. Clusters blend into each other. Continuous gradients become more apparent. | Visualizing developmental trajectories; data with continuous phenotypic transitions; avoiding the false impression of discrete populations. |

**Visual intuition:** `min_dist` acts like a "packing pressure" knob. At low values, points are squeezed together, producing island-like clusters with visible gaps. At high values, points are allowed to spread, producing a more diffuse, connected landscape.

### `metric` (default: "euclidean")

The distance metric determines how similarity between cells is computed in the original high-dimensional marker space.

| Metric | Formula | When to Use |
|--------|---------|-------------|
| **euclidean** (default) | Standard L2 norm. | Transformed cytometry data (arcsinh or logicle). Appropriate when absolute marker intensities are meaningful. This is the standard choice for cytometry. |
| **cosine** | 1 minus the cosine of the angle between vectors. | When relative marker profiles matter more than absolute intensities. Useful for comparing cell "signatures" independent of overall staining intensity. |
| **correlation** | 1 minus Pearson correlation. | Similar rationale to cosine, but also accounts for mean centering. Rarely needed in standard cytometry workflows. |

For most flow cytometry analyses using [arcsinh-transformed](arcsinh-transform.md) or [logicle-transformed](logicle-transform.md) data, Euclidean distance is appropriate and well validated in [published benchmarks](https://doi.org/10.1038/nbt.4314).

### `seed` (default: 42)

UMAP is a stochastic algorithm: the approximate nearest-neighbor search, spectral initialization, and SGD optimization all involve randomness. The `seed` parameter sets `set.seed()` before computation, ensuring [identical embeddings across runs](https://doi.org/10.48550/arXiv.1802.03426) given the same input data and parameters.

**Always set a seed.** Omitting the seed means different runs produce different layouts, making it impossible to visually compare results or reproduce figures. In this pipeline, the seed is set by default in both `run_umap()` and the plotting wrappers.

---

## When Your Data Calls For It

UMAP is the right tool in the following concrete scenarios:

### Visualize Cell Populations in 2D

The primary use case. After preprocessing and (optionally) batch correction, project all cells into a 2D scatter plot to see how they group. This gives an immediate visual overview of the phenotypic landscape across your panel.

### Overlay Cluster Labels

After running [Leiden clustering](leiden-clustering.md) on the original expression matrix, map cluster labels onto the UMAP coordinates. This validates that clusters correspond to visually coherent regions and helps identify clusters that may need splitting or merging. The pipeline's `plot_umap_from_matrix()` accepts a `cluster_labels` parameter for this purpose.

### Overlay Marker Expression

Color-code each cell by its expression of a specific marker using the per-marker UMAP plots generated by both `plot_umap()` and `plot_umap_from_matrix()`. This reveals which markers drive the separation between populations and is essential for manual annotation of clusters.

### Compare Conditions or Samples

Color-code cells by sample identity or experimental condition to see whether populations shift, expand, or contract. The sample-colored UMAP plot (generated automatically as `umap_sample`) is the first diagnostic output in Pipelines A, B, and E.

### Quality Control of Clustering Results

UMAP provides a visual sanity check for algorithmic clustering. If a cluster appears scattered across the UMAP rather than forming a contiguous region, it may indicate a too-low resolution parameter, a batch effect, or a clustering artifact.

---

## When to Choose Something Else

### Quantitative Distance Measurements

UMAP is a visualization tool, not a quantitative analysis method. [The distances between points and between clusters in UMAP space do not reliably reflect distances in the original high-dimensional space](https://doi.org/10.1038/s41587-023-01986-3). If your question requires measuring how different two populations are, operate directly in the original marker space or use well-defined statistical tests (see `src/analysis/`).

### Variance Decomposition and Feature Importance

If you need to know which markers explain the most variance in your data, or you need a linear projection with interpretable loadings, use [PCA](pca.md). PCA preserves global linear structure and produces loadings that quantify each marker's contribution. UMAP's nonlinear mapping does not provide this.

### Clustering

Never cluster on UMAP coordinates. UMAP distorts distances and densities, so clusters found in UMAP space may not reflect real population structure. Instead, run [Leiden clustering](leiden-clustering.md) in the original high-dimensional marker space (or on a PCA-reduced space), and then overlay those labels onto the UMAP for visualization. The pipeline enforces this separation: `run_leiden()` in `src/clustering/leiden.R` builds its KNN graph from the expression matrix, not from UMAP output.

### Very Small Datasets

With fewer than approximately 100 cells, UMAP's neighborhood-based approach has insufficient data to reconstruct meaningful topology. For small datasets, direct inspection of bivariate scatter plots or PCA is more appropriate.

---

## Pre-processing Requirements

UMAP operates on a numeric matrix of marker expression values. The quality of the embedding depends critically on the preprocessing applied upstream.

### Transform First (Required)

Raw fluorescence or mass cytometry intensities span orders of magnitude and contain near-zero and negative values. Applying UMAP to untransformed data will produce embeddings dominated by a few high-intensity markers. Always transform before running UMAP:

- **[Arcsinh transform](arcsinh-transform.md):** Standard for mass cytometry (CyTOF). Uses a cofactor (typically 5 for CyTOF, 150 for fluorescence) to compress the dynamic range while preserving relative differences.
- **[Logicle / biexponential transform](logicle-transform.md):** Preferred for fluorescence cytometry, particularly when data includes negative values from compensation.

In the pipeline, transformation is handled in `src/preprocessing/` before data reaches the visualization stage.

### Batch Correction (Optional but Recommended)

If your experiment includes multiple batches (different run days, instruments, or staining panels), batch effects can dominate the UMAP embedding, causing cells to cluster by batch rather than by biology. Apply [batch correction](batch-correction.md) (e.g., Harmony via `src/harmonization/`) before UMAP to remove these technical artifacts.

### PCA Pre-reduction (Optional)

For datasets with many markers (>30) or very large cell counts, running PCA first and then applying UMAP to the top principal components can [improve both speed and embedding quality](https://doi.org/10.1038/nbt.4314). The PCA step removes noise in low-variance dimensions and reduces the computational burden of the nearest-neighbor search. This is available via `run_pca()` in `src/R_shared/stats_utils.R`.

### Marker Selection

Both `plot_umap()` and `plot_umap_from_matrix()` automatically exclude scatter parameters (FSC, SSC), Time, and index columns from the UMAP computation. Only phenotypic markers are used. You can further restrict the marker set via the `markers_of_interest` parameter.

---

## How to Use in the Pipeline

### From flowFrames (Pipelines A and E)

For data loaded from FCS files, use `plot_umap()` which handles flowFrame extraction internally:

```r
source("src/visualization/umap.R")

result <- plot_umap(
  flow_frames      = downsampled_frames,   # Named list of flowFrame objects
  sample_names     = sample_names,          # Character vector
  markers_of_interest = markers,            # Markers for embedding
  output_dir       = "output/umap/",
  n_neighbors      = 15,
  min_dist         = 0.1,
  seed             = 42
)

# result$umap_coords  -> Nx2 matrix of coordinates
# result$plots        -> list of ggplot objects
```

### From Expression Matrix (Pipeline B, or Cross-pipeline)

For CSV-derived data or any numeric matrix, use `plot_umap_from_matrix()`:

```r
source("src/visualization/umap.R")

result <- plot_umap_from_matrix(
  expr_matrix      = combined_matrix,       # Numeric matrix (cells x markers)
  sample_labels    = sample_labels,          # Character vector (one per row)
  cluster_labels   = clusters$cluster,       # Optional: overlay cluster colors
  output_dir       = "output/umap/",
  n_neighbors      = 15,
  min_dist         = 0.1,
  seed             = 42
)
```

### Reuse Pre-computed Coordinates

If UMAP coordinates were already computed (e.g., during clustering), pass them directly to avoid redundant computation:

```r
source("src/visualization/umap.R")

result <- plot_umap_from_matrix(
  expr_matrix      = combined_matrix,
  umap_coords      = clusters$umap_coords,  # Pre-computed Nx2 matrix
  cluster_labels   = clusters$cluster,
  output_dir       = "output/umap/"
)
```

### Standalone UMAP Computation (Without Plotting)

If you only need coordinates without automatic plot generation, call `run_umap()` directly:

```r
source("src/R_shared/stats_utils.R")

coords <- run_umap(
  expr_matrix  = my_matrix,
  n_neighbors  = 15,
  min_dist     = 0.1,
  metric       = "euclidean",
  seed         = 42
)
# coords is an Nx2 numeric matrix
```

### Output Files

Both `plot_umap()` and `plot_umap_from_matrix()` automatically save plots to `output_dir` in the configured formats (default: PNG and SVG). Generated files include:

| File | Content |
|------|---------|
| `umap_sample.{png,svg}` | Cells colored by sample identity |
| `umap_cluster.{png,svg}` | Cells colored by cluster label (if `cluster_labels` provided) |
| `umap_{marker}.{png,svg}` | One file per marker, colored by expression intensity (viridis palette) |

---

## Common Pitfalls

### Over-interpreting Distances Between Clusters

The most common misinterpretation of UMAP plots is reading meaning into the distances between clusters. [Chari and Pachter (2023)](https://doi.org/10.1038/s41587-023-01986-3) demonstrated that UMAP embeddings can [place unrelated populations close together and related populations far apart](https://doi.org/10.1038/s41587-023-01986-3) depending on parameter settings and random initialization. The inter-cluster distances in a UMAP plot are not a reliable proxy for biological similarity. Always validate apparent relationships using the original high-dimensional data.

### Over-interpreting Cluster Sizes and Densities

The visual area occupied by a cluster in UMAP space [does not reliably reflect the number of cells in that cluster or the local density in marker space](https://doi.org/10.1038/s41587-023-01986-3). UMAP's optimization can expand sparse regions and compress dense regions to satisfy the embedding constraints. Use abundance barplots (`src/visualization/barplots.R`) or summary tables for quantitative cell counts.

### Over-interpreting Gaps and Continuity

A gap between two regions in a UMAP plot does not necessarily mean those populations are biologically distinct, and a continuous bridge does not necessarily mean there is a phenotypic continuum. These visual features are [artifacts of the embedding parameters and random initialization](https://doi.org/10.1038/s41587-023-01986-3). Validate apparent gaps or continuities by examining marker expression distributions directly.

### Running UMAP on Untransformed Data

Applying UMAP to raw fluorescence intensities will produce misleading embeddings where a few highly expressed markers dominate the structure. Always apply arcsinh or logicle transformation first.

### Using Too Few Cells for Stable Embeddings

UMAP needs sufficient data to build a reliable neighborhood graph. With fewer than a few hundred cells, the embedding becomes unstable and sensitive to individual data points. For very small datasets, use PCA or bivariate scatter plots instead.

### Changing Parameters Without Re-evaluating

Different parameter settings can produce dramatically different-looking embeddings from the same data. Always report your parameter choices, and be cautious about comparing UMAP plots generated with different settings.

---

## UMAP vs t-SNE

Both UMAP and [t-SNE (van der Maaten and Hinton, 2008)](https://doi.org/10.1007/s10994-008-5065-1) are nonlinear dimensionality reduction methods widely used in single-cell analysis. The table below summarizes their key differences as characterized in [comparative evaluations](https://doi.org/10.1038/nbt.4314):

| Aspect | UMAP | t-SNE |
|--------|------|-------|
| **Global structure** | [Preserves relative positioning of distinct populations](https://doi.org/10.1038/nbt.4314). Biologically similar cell types tend to appear near each other. | [Optimizes local neighborhoods but scrambles global arrangement](https://doi.org/10.1007/s10994-008-5065-1). Cluster positions are arbitrary. |
| **Speed** | [Significantly faster](https://doi.org/10.1038/nbt.4314), especially at scale. Handles >10^6 cells without subsampling in practical runtimes. | Slower. Often requires subsampling or approximations (e.g., FIt-SNE) for large datasets. |
| **Determinism** | Stochastic, but [more stable across random seeds](https://doi.org/10.1038/nbt.4314). Combined with fixed seed, produces identical output. | Highly sensitive to random initialization and perplexity. Different runs can look very different. |
| **Loss function** | Cross-entropy between fuzzy simplicial sets. [Grounded in topological data analysis](https://doi.org/10.48550/arXiv.1802.03426). | KL-divergence between conditional probability distributions in high and low dimensions. |
| **Key parameter** | `n_neighbors` (neighborhood size). | Perplexity (related concept but implemented differently). |
| **Cluster appearance** | Compact but varied in size and shape. Internal structure often visible. | Tends to produce [uniformly sized, spherical-looking clusters](https://doi.org/10.1038/nbt.4314) regardless of true population geometry. |
| **Continuous structures** | Better at revealing trajectories and gradients. | Tends to fragment continuous manifolds into discrete blobs. |

### When t-SNE May Still Be Preferred

Despite UMAP's advantages, t-SNE remains useful in specific scenarios: (1) when direct visual comparison with legacy analyses performed using t-SNE is required, (2) when the [perplexity-annealing variants of t-SNE](https://doi.org/10.1007/s10994-008-5065-1) are used to explicitly control local vs. global emphasis, or (3) when the research community or journal reviewers expect t-SNE figures. In this pipeline, UMAP is the default and recommended method.

---

## Reproducibility

### Why Reproducibility Matters

UMAP is a stochastic algorithm. Without controlling for randomness, [repeated runs on the same data with the same parameters will produce visually different embeddings](https://doi.org/10.48550/arXiv.1802.03426). This makes it impossible to compare figures across analyses or verify results.

### How the Pipeline Ensures Reproducibility

Three sources of stochasticity in UMAP are controlled by the `seed` parameter:

1. **Approximate nearest-neighbor search.** The `uwot` package uses random projection trees and nearest-neighbor descent, both of which involve random initialization. Setting `set.seed()` before calling `uwot::umap()` makes these deterministic.

2. **Spectral initialization.** The initial low-dimensional layout is derived from the graph Laplacian, which involves eigenvector computation that can be sensitive to numerical precision. The seed stabilizes any random components.

3. **Stochastic gradient descent.** The optimization loop samples negative examples randomly. A fixed seed ensures the same sequence of samples.

In `run_umap()` (`src/R_shared/stats_utils.R`), `set.seed(seed)` is called immediately before `uwot::umap()`:

```r
run_umap <- function(expr_matrix, n_neighbors = 15, min_dist = 0.1,
                     n_components = 2, metric = "euclidean", seed = 42) {
  set.seed(seed)
  uwot::umap(expr_matrix,
             n_neighbors = n_neighbors,
             min_dist = min_dist,
             n_components = n_components,
             metric = metric)
}
```

### Reporting Checklist

When publishing or sharing UMAP figures, always report:

- **Software and version:** `uwot` R package (specify version from `packageVersion("uwot")`)
- **Parameters:** `n_neighbors`, `min_dist`, `metric`, `n_components`
- **Random seed:** The exact seed used
- **Input data:** Number of cells, number of markers, transformation applied, whether batch correction was applied
- **Preprocessing:** Which markers were included/excluded from the embedding

This information allows any collaborator to regenerate the exact same embedding from the same input data.

---

## Implementation Details

### Package

This pipeline uses the [`uwot`](https://cran.r-project.org/package=uwot) R package, which provides a [fast C++ implementation of UMAP](https://doi.org/10.48550/arXiv.1802.03426) with the Annoy library for approximate nearest-neighbor search. The `uwot` package is the most widely used UMAP implementation in the R ecosystem for single-cell analysis.

### Computational Complexity

- **Nearest-neighbor search:** O(N log N) via Annoy random projection trees, where N is the number of cells.
- **Graph construction:** O(N * k) where k = `n_neighbors`.
- **Optimization:** O(N * n_epochs) with stochastic gradient descent. Default epochs scale with dataset size (larger datasets use fewer epochs per point).

For a typical cytometry dataset of 100,000 cells and 20 markers, UMAP computation takes approximately 30--60 seconds on a modern workstation.

### Memory Usage

UMAP's memory footprint is dominated by the nearest-neighbor graph, which requires O(N * k) storage. For 1 million cells with `n_neighbors = 15`, this is approximately 120 MB. The expression matrix itself requires additional memory proportional to N * number of markers. For very large datasets, consider downsampling via `src/preprocessing/downsampling.R` before UMAP computation.
