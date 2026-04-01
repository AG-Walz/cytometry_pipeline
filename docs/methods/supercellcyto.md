# SuperCellCyto: Metacell Compression for Cytometry

Reduces large cytometry datasets by grouping phenotypically similar cells into "supercells" (metacells), preserving biological signal while compressing cell count by 10--50x to accelerate downstream clustering, visualization, and statistical analysis.

---

## When to Use

- Your dataset is very large (>500,000 cells) and downstream analyses (clustering, UMAP, differential testing) are prohibitively slow or memory-intensive.
- You want to reduce computational cost without the information loss of random downsampling.
- You need a compressed representation that preserves rare populations better than uniform subsampling.
- You want to run ensemble analyses (e.g., multiple clustering resolutions) on a tractable dataset.

SuperCellCyto is preferable to random downsampling when you want to retain the full phenotypic diversity of the data. If your dataset is small enough for direct analysis (<100,000 cells), metacell compression is unnecessary.

---

## How It Works

SuperCellCyto constructs a nearest-neighbor graph on the single-cell data, then partitions it into groups of highly similar cells, replacing each group with a single representative supercell.

In the first stage, the expression matrix is optionally reduced by PCA to the top principal components (default: 10) to denoise and speed up neighbor search. A k-nearest-neighbor graph is then constructed on the (optionally PCA-reduced) data, connecting each cell to its most similar neighbors. The graph captures the local phenotypic structure of the data at single-cell resolution.

In the second stage, the KNN graph is coarsened using a walktrap-based algorithm that groups cells into supercells. The graining level parameter (`gam`) controls the compression ratio: a value of 20 means approximately 1 supercell per 20 original cells, yielding a 20x reduction. Within each supercell, cells are phenotypically homogeneous because they form a tightly connected subgraph in the KNN graph. The expression profile of each supercell is computed as the weighted average of its member cells' expression values, preserving the mean phenotype of the group.

The key advantage over random downsampling is that SuperCellCyto guarantees local coverage of the phenotypic landscape. Random downsampling with a 20x reduction has a high probability of losing rare populations that comprise <5% of the data. SuperCellCyto, by contrast, will create at least one supercell for any population large enough to form a connected neighborhood in the KNN graph. The mapping between supercells and their constituent single cells is retained, allowing results to be projected back to the original resolution after analysis. Supercell sizes (number of member cells) are also recorded, enabling size-weighted statistical analyses.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expr_matrix` | (required) | Numeric matrix (cells x markers) |
| `sample_labels` | (required) | Character vector (one per row) mapping cells to samples |
| `markers` | (required) | Character vector of marker names to use for supercell construction |
| `gam` | 20 | Graining level: higher = fewer supercells (approximate compression ratio) |
| `k_knn` | 30 | Nearest neighbors for graph construction |
| `n_pc` | 10 | Number of PCs for initial dimensionality reduction |
| `do_pca` | TRUE | Whether to run PCA before supercell construction |
| `seed` | 42 | Random seed for reproducibility |

**Choosing `gam`**: A graining level of 10 provides moderate compression (10x), preserving fine-grained structure. A level of 20 (default) provides strong compression suitable for most workflows. Levels of 50--100 provide aggressive compression for very large datasets (>5M cells) but may merge small populations.

**Choosing `k_knn`**: Higher values produce smoother supercells (more averaging). Lower values preserve sharper boundaries between populations. The default of 30 matches the standard KNN parameter used throughout the pipeline.

---

## Pipeline Integration

**Script**: `src/preprocessing/supercell.R`

**Function**: `run_supercell(expr_matrix, sample_labels, markers, gam, k_knn, n_pc, do_pca, seed)`

**Input from**:
- `build_combined_expr_matrix()` in `src/R_shared/flowframe_utils.R`
- `combine_expression_matrices()` in `src/analysis/expression_matrices.R`

**Output to**:
- `run_leiden()` / `run_flowsom()` / `run_phenograph()` in `src/clustering/`
- `plot_umap_from_matrix()` / `plot_tsne_from_matrix()` in `src/visualization/`
- `run_differential_abundance()` / `run_differential_expression()` in `src/analysis/`

**Returns**: List containing `supercell_matrix` (compressed expression matrix), `supercell_labels` (sample assignments), `cell_map` (supercell-to-cell mapping), and `supercell_sizes` (number of cells per supercell).

**Placement in pipeline**: Apply after transformation and optional batch correction, before clustering and visualization. The sequence is: Load -> Transform -> (Batch Correct) -> SuperCellCyto -> Cluster -> Analyze.

**Used in**: Any pipeline (A through E) when the input dataset exceeds practical size limits for direct analysis.

---

## Compression Ratios

| Input Size | `gam` = 10 | `gam` = 20 | `gam` = 50 |
|-----------|------------|------------|------------|
| 100,000 cells | ~10,000 supercells | ~5,000 supercells | ~2,000 supercells |
| 500,000 cells | ~50,000 supercells | ~25,000 supercells | ~10,000 supercells |
| 2,000,000 cells | ~200,000 supercells | ~100,000 supercells | ~40,000 supercells |

---

## References

1. Putri, G. H., Howitt, G., Marsh-Wakefield, F. et al. (2024). SuperCellCyto: Enabling efficient analysis of large scale cytometry datasets. *Bioinformatics*, 40(4), btae125. [doi:10.1093/bioinformatics/btae125](https://doi.org/10.1093/bioinformatics/btae125)

2. Bilous, M., Dessimoz, C. & Gfeller, D. (2022). Metacells untangle large and complex single-cell transcriptome networks. *BMC Bioinformatics*, 23, 336. [doi:10.1186/s12859-022-04861-1](https://doi.org/10.1186/s12859-022-04861-1)

3. Weber, L. M. & Robinson, M. D. (2016). Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data. *Cytometry Part A*, 89(12), 1084-1096. [doi:10.1002/cyto.a.22764](https://doi.org/10.1002/cyto.a.22764)
