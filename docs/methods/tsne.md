# t-SNE: t-Distributed Stochastic Neighbor Embedding

Nonlinear dimensionality reduction that projects high-dimensional single-cell data into 2D by preserving local neighborhood structure through probability distributions with heavy-tailed kernels.

---

## When to Use

- You need 2D visualization of cell populations and want to compare directly with legacy analyses that used t-SNE.
- Your dataset is small to moderate (up to ~100,000 cells) where t-SNE runtime is acceptable.
- You want to emphasize local cluster separation over global arrangement.
- Journal reviewers or collaborators specifically expect t-SNE figures.

For most new analyses, [UMAP](umap.md) is preferred due to better global structure preservation, faster runtime, and greater reproducibility. Use t-SNE when backward compatibility with published results is important.

---

## How It Works

t-SNE operates by matching probability distributions in high-dimensional and low-dimensional spaces through KL-divergence minimization.

In the first stage, the algorithm computes pairwise similarities in the original high-dimensional space. For each cell, a Gaussian distribution is centered on that cell, and the conditional probability that it would pick each other cell as its neighbor is computed. The width of each Gaussian is calibrated so that its effective number of neighbors (measured by perplexity) matches a user-specified target. These conditional probabilities are symmetrized into joint probabilities that represent the similarity structure of the data. Perplexity acts as a smooth analog of the number of effective neighbors and typically ranges from 5 to 50 for cytometry data.

In the second stage, a low-dimensional (2D) embedding is initialized randomly and iteratively refined. In the low-dimensional space, pairwise similarities are computed using a Student's t-distribution with one degree of freedom (a Cauchy distribution) rather than a Gaussian. The heavy tails of the t-distribution allow moderate distances in high-dimensional space to be represented as larger distances in 2D, preventing the "crowding problem" where points collapse into the center. The embedding is optimized by gradient descent on the KL-divergence between the high-dimensional and low-dimensional similarity distributions.

The t-distribution kernel is what gives t-SNE its characteristic appearance: tight, well-separated clusters with clear visual boundaries. However, this same property means that cluster positions in the 2D embedding are largely arbitrary. [Van der Maaten and Hinton (2008)](https://doi.org/10.1007/s10994-008-5065-1) showed t-SNE excels at revealing local neighborhood structure but does not preserve global distances or relative cluster positions.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expr_matrix` | (required) | Numeric matrix, rows = cells, columns = markers |
| `perplexity` | 30 | Effective number of neighbors; balances local vs. global focus |
| `max_iter` | 1000 | Maximum optimization iterations |
| `seed` | 42 | Random seed for reproducibility |

**Choosing perplexity**: Values of 5--15 emphasize fine local structure (small, tight clusters). Values of 30--50 provide a broader view. The embedding is sensitive to perplexity, and different values can produce visually different results from the same data. Always report the perplexity used.

**Iteration count**: 1000 iterations is standard. Fewer iterations may produce an unconverged embedding with artifactual structure. Monitor the KL-divergence to confirm convergence.

---

## Pipeline Integration

**Script**: `src/visualization/tsne.R`

**Functions**:
- `plot_tsne(flow_frames, ...)` -- flowFrame input
- `plot_tsne_from_matrix(expr_matrix, ...)` -- matrix input

**Input from**:
- `load_fcs_raw()` in `src/loading/load_fcs_raw.R`
- `downsample_flowframes()` in `src/preprocessing/downsampling.R`
- `build_combined_expr_matrix()` in `src/R_shared/flowframe_utils.R`

**Output to**: Standalone visualization output (PNG/SVG plots with sample, cluster, and per-marker overlays).

**Used in**: Pipelines A, B, and E as an alternative to UMAP visualization. The implementation uses `Rtsne::Rtsne()` for the Barnes-Hut approximation, which reduces the O(n^2) exact computation to O(n log n).

---

## t-SNE vs UMAP

| Aspect | t-SNE | UMAP |
|--------|-------|------|
| Global structure | Cluster positions are arbitrary | Preserves relative positioning of populations |
| Speed | Slower, especially >100k cells | Significantly faster at scale |
| Reproducibility | Highly sensitive to initialization and perplexity | More stable across random seeds |
| Cluster appearance | Uniformly sized spherical blobs | Variable sizes reflecting actual population structure |
| Continuous structures | Tends to fragment trajectories into discrete blobs | Better at revealing gradients and trajectories |
| Loss function | KL-divergence | Cross-entropy on fuzzy simplicial sets |

Both methods share the same critical limitation: distances between clusters in the 2D embedding do not reliably reflect distances in the original marker space. Never cluster on t-SNE or UMAP coordinates; always cluster in the original expression space and overlay labels on the embedding.

---

## References

1. Van der Maaten, L. & Hinton, G. (2008). Visualizing data using t-SNE. *Journal of Machine Learning Research*, 9, 2579-2605. [doi:10.1007/s10994-008-5065-1](https://doi.org/10.1007/s10994-008-5065-1)

2. Becht, E., McInnes, L., Healy, J. et al. (2019). Dimensionality reduction for visualizing single-cell data using UMAP. *Nature Biotechnology*, 37, 38-44. [doi:10.1038/nbt.4314](https://doi.org/10.1038/nbt.4314)

3. Amir, E. D., Davis, K. L., Tadmor, M. D. et al. (2013). viSNE enables visualization of high dimensional single-cell data and reveals phenotypic heterogeneity of leukemia. *Nature Biotechnology*, 31, 545-552. [doi:10.1038/nbt.2594](https://doi.org/10.1038/nbt.2594)
