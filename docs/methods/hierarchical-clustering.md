# Hierarchical Agglomerative Clustering (Ward.D2)

## Summary

Hierarchical agglomerative clustering is a bottom-up algorithm that iteratively merges the two most similar items (samples, clusters, or features) into a single group until all items belong to one tree-shaped hierarchy called a dendrogram. In this pipeline, the Ward.D2 linkage criterion -- which [minimizes the total within-cluster variance at each merge step using squared Euclidean distances](https://doi.org/10.1007/s00357-014-9161-z) -- is used to order samples and clusters in abundance barplots and heatmaps, producing layouts where similar items appear adjacent. This is **not** a cell-level clustering method; for assigning individual cells to populations, see [Leiden clustering](leiden-clustering.md).

---

## The Method

### Agglomerative (Bottom-Up) Strategy

Hierarchical agglomerative clustering starts with every item as its own singleton cluster and proceeds by repeatedly merging the pair of clusters with the smallest inter-cluster distance. [Ward (1963)](https://doi.org/10.1080/01621459.1963.10500845) proposed a merge criterion based on minimizing the increase in total within-cluster sum of squares (the "error sum of squares") at each step. The algorithm produces a nested sequence of partitions that can be represented as a binary tree -- the dendrogram -- where the height of each internal node corresponds to the dissimilarity at which the two child clusters were merged. This dendrogram provides a [complete record of the clustering hierarchy at every possible granularity](https://doi.org/10.1007/b98835), allowing the analyst to choose a cut level after the fact or, as in this pipeline, to use the leaf ordering directly for visualization.

The algorithm proceeds as follows:

1. **Initialize**: Assign each item (sample or cluster) to its own singleton cluster. Compute the pairwise distance matrix between all items.
2. **Find closest pair**: Identify the two clusters with the smallest inter-cluster distance according to the chosen linkage criterion.
3. **Merge**: Combine the two clusters into a new cluster. Record the merge height in the dendrogram.
4. **Update distances**: Recompute distances between the new cluster and all remaining clusters using the [Lance-Williams recurrence formula](https://doi.org/10.1093/comjnl/9.4.373), which expresses the new distances in terms of the distances before the merge. This avoids recomputing pairwise distances from scratch.
5. **Repeat** steps 2--4 until a single cluster remains.

The computational complexity is O(n^3) for the naive implementation and [O(n^2 log n) with optimized priority-queue-based approaches](https://doi.org/10.1007/b98835), where n is the number of items. For the sample-level and cluster-level matrices used in this pipeline (typically tens to hundreds of items), this cost is negligible.

### Ward's Criterion: Minimizing Within-Cluster Variance

[Ward (1963)](https://doi.org/10.1080/01621459.1963.10500845) formulated the clustering objective as minimizing the total within-cluster sum of squares (WCSS), also called the error sum of squares (ESS). At each merge step, the algorithm evaluates every possible pair of clusters and selects the merge that produces the [smallest increase in WCSS](https://doi.org/10.1080/01621459.1963.10500845). For two clusters A and B with centroids c_A and c_B containing n_A and n_B members respectively, the increase in WCSS upon merging is:

```
delta(A, B) = (n_A * n_B) / (n_A + n_B) * ||c_A - c_B||^2
```

This criterion [favors merging compact, equally-sized clusters](https://doi.org/10.1080/01621459.1963.10500845) and tends to produce clusters of roughly equal variance, which is desirable when the goal is to reveal balanced groupings of similar items rather than chaining together outliers.

Ward's method is [equivalent to k-means-like variance minimization applied greedily at each level of the hierarchy](https://doi.org/10.1007/s00357-014-9161-z). This makes it the natural hierarchical counterpart to the k-means objective, and it produces dendrograms whose structure aligns well with variance-based cluster validity indices.

### Ward.D2: The Correct Implementation

R's `hclust()` function provides two Ward variants: `method = "ward.D"` and `method = "ward.D2"`. The distinction is subtle but important.

[Murtagh and Legendre (2014)](https://doi.org/10.1007/s00357-014-9161-z) demonstrated that the original `"ward.D"` implementation in R (prior to R 3.1.0, when it was the only Ward option) does **not** implement Ward's criterion correctly when given a Euclidean distance matrix computed by `dist()`. The problem arises because `dist()` returns Euclidean distances (not squared distances), while the Lance-Williams update formula for Ward's method [requires squared Euclidean distances as input](https://doi.org/10.1007/s00357-014-9161-z). Passing unsquared distances to the `"ward.D"` algorithm applies the Lance-Williams recurrence to the wrong scale, producing a dendrogram that [does not correspond to Ward's minimum-variance criterion](https://doi.org/10.1007/s00357-014-9161-z).

`method = "ward.D2"` solves this by [internally squaring the input distances before applying the Lance-Williams formula](https://doi.org/10.1007/s00357-014-9161-z), then taking the square root of the merge heights for the dendrogram. This yields a dendrogram that correctly reflects the minimum-variance merges regardless of whether the user provides squared or unsquared distances. As Murtagh and Legendre showed, [the `"ward.D2"` implementation produces dendrograms identical to those from other software packages (SAS, SPSS) that correctly implement Ward's method](https://doi.org/10.1007/s00357-014-9161-z).

**This pipeline uses `"ward.D2"` exclusively** to ensure correct minimum-variance clustering.

---

## Ward.D2 vs Ward.D

This distinction warrants its own section because using the wrong variant is a common source of irreproducible results across software platforms.

### The Problem with Ward.D

When `hclust()` receives a distance matrix `d` from `dist()` and is called with `method = "ward.D"`, it applies the Lance-Williams update formula directly to the values in `d`. However, [Ward's criterion is defined in terms of squared Euclidean distances](https://doi.org/10.1007/s00357-014-9161-z). The Lance-Williams formula for Ward's method is:

```
d(AB, C) = sqrt(
  ((n_A + n_C) * d(A,C)^2 + (n_B + n_C) * d(B,C)^2 - n_C * d(A,B)^2) /
  (n_A + n_B + n_C)
)
```

When `d` contains unsquared Euclidean distances, `"ward.D"` effectively [applies this formula to the wrong values](https://doi.org/10.1007/s00357-014-9161-z), producing a dendrogram that minimizes a different (and generally uninterpretable) objective function.

### The Ward.D2 Correction

[Murtagh and Legendre (2014)](https://doi.org/10.1007/s00357-014-9161-z) proposed two solutions:

1. **Pre-square the distances**: Pass `dist(x)^2` to `hclust(..., method = "ward.D")`. This works but is error-prone and non-obvious to users.
2. **Use `method = "ward.D2"`** (introduced in R 3.1.0): This variant [internally handles the squaring and square-root operations](https://doi.org/10.1007/s00357-014-9161-z), so that `hclust(dist(x), method = "ward.D2")` produces the correct Ward dendrogram regardless of whether the input distances are squared or unsquared.

The practical consequences are:

| Aspect | `ward.D` with `dist()` | `ward.D2` with `dist()` |
|--------|------------------------|--------------------------|
| Implements Ward's criterion | No (unless distances are pre-squared) | Yes |
| Dendrogram heights meaningful | No (arbitrary scale) | Yes (square root of WCSS increase) |
| Matches SAS/SPSS output | No | [Yes](https://doi.org/10.1007/s00357-014-9161-z) |
| Cross-platform reproducibility | Poor | Good |

[Murtagh and Legendre explicitly recommended `"ward.D2"` as the default choice for all new analyses](https://doi.org/10.1007/s00357-014-9161-z), and this pipeline follows that recommendation.

---

## Role in This Pipeline

Hierarchical clustering with Ward.D2 serves a **visualization and ordering** role. It is not used for primary cell population identification -- that task belongs to [Leiden clustering](leiden-clustering.md), which operates on single-cell expression data using graph-based community detection. The distinction is important:

| Task | Method | Scale |
|------|--------|-------|
| Assigning cells to populations | [Leiden](leiden-clustering.md) | 10,000s--1,000,000s of cells |
| Ordering samples in plots | Ward.D2 hierarchical clustering | 10s--100s of samples |
| Ordering clusters in heatmaps | Ward.D2 hierarchical clustering | 5--50 clusters |
| Biclustering heatmap rows and columns | Ward.D2 / pheatmap internal | Both dimensions |

### Ordering Samples in Barplots by Composition Similarity

In `plot_abundance_proportional()` (`src/analysis/abundance_plotting.R`), the sample-by-cluster proportion matrix is hierarchically clustered to determine sample order on the x-axis. Samples with similar cluster compositions appear adjacent, and a dendrogram is drawn above the barplot to show the grouping structure. This makes it visually immediate which samples share similar immune profiles, a pattern that [would be invisible with arbitrary or alphabetical sample ordering](https://doi.org/10.1093/bioinformatics/btr593).

The relevant code:

```r
prop_mat <- mat / rowSums(mat)
hc_samples <- hclust(dist(prop_mat), method = "ward.D2")
dendro_order <- hc_samples$labels[hc_samples$order]
long_df$file <- factor(long_df$file, levels = dendro_order)
```

### Ordering Clusters and Samples in Heatmaps

In `plot_abundance_heatmap()` (`src/analysis/abundance_plotting.R`), both rows (samples) and columns (clusters) are independently hierarchically clustered using Ward.D2. The resulting biclustered heatmap places samples with similar cluster distributions on adjacent rows and clusters with similar cross-sample abundance profiles in adjacent columns. This biclustering approach [reveals block structure in the abundance matrix that corresponds to biologically meaningful sample groups and co-regulated cell populations](https://doi.org/10.1093/bioinformatics/btr593).

```r
hc_samples  <- hclust(dist(prop_mat), method = "ward.D2")
hc_clusters <- hclust(dist(t(prop_mat)), method = "ward.D2")
hm_mat <- prop_mat[hc_samples$order, hc_clusters$order]
```

### Expression Heatmaps via pheatmap

In `plot_heatmap()` and `plot_cluster_heatmap()` (`src/visualization/heatmaps.R`), the `pheatmap` package performs its own internal hierarchical clustering when `cluster_rows = TRUE` and/or `cluster_cols = TRUE`. By default, [pheatmap uses `hclust()` with `method = "complete"` and `dist()` with Euclidean distance](https://cran.r-project.org/package=pheatmap) to cluster rows and columns independently, rendering dendrograms along the margins. This provides an at-a-glance view of which clusters have similar expression profiles and which markers co-vary across clusters.

The `pheatmap` package also accepts custom distance functions and clustering methods via its `clustering_distance_rows`, `clustering_distance_cols`, and `clustering_method` parameters, allowing Ward.D2 to be specified if the default complete linkage is not appropriate for a given analysis.

---

## When Your Data Calls For It

Hierarchical clustering with Ward.D2 is the right tool for the following tasks in this pipeline:

### Organizing Heatmap Rows and Columns

When visualizing an expression matrix or abundance matrix as a heatmap, the ordering of rows and columns determines whether patterns are visible or hidden. [Random or alphabetical ordering of samples and markers produces heatmaps where block structure is obscured by interleaving](https://doi.org/10.1093/bioinformatics/btr593). Hierarchical clustering reorders both dimensions so that similar rows are adjacent and similar columns are adjacent, [revealing clusters-of-clusters that correspond to biological groupings](https://doi.org/10.1093/bioinformatics/btr593). This is the standard approach in [virtually all heatmap implementations used in genomics and cytometry](https://doi.org/10.1093/bioinformatics/btw313).

### Sample Ordering in Abundance Barplots

Stacked or proportional barplots showing cluster composition per sample benefit enormously from hierarchical ordering. When samples are ordered by dendrogram, the resulting barplot immediately shows [which groups of samples have similar immune profiles](https://doi.org/10.12688/f1000research.11622.3), which clusters drive the differences between groups, and whether there are outlier samples with atypical compositions. The paired dendrogram above the barplot provides a visual quantification of sample-to-sample similarity.

### Identifying Groups of Samples with Similar Cluster Compositions

Cutting the sample dendrogram at a chosen height partitions the samples into groups that share similar cluster abundance profiles. While this pipeline does not programmatically cut the dendrogram, the visual grouping [serves as a hypothesis-generating tool](https://doi.org/10.12688/f1000research.11622.3) -- for example, revealing that patients with similar treatment outcomes cluster together based on immune composition, or that samples from the same batch separate from others (indicating a batch effect).

### Identifying Groups of Clusters with Similar Expression Profiles

When the cluster dimension of a heatmap is hierarchically clustered, the column dendrogram shows which clusters have similar marker expression profiles. Clusters that merge at low heights in the dendrogram may represent subpopulations of the same cell type [differing only in one or two markers](https://doi.org/10.12688/f1000research.11622.3), suggesting that the upstream Leiden resolution may be too high for the biological question at hand. Conversely, clusters that remain separated until high merge heights represent phenotypically distinct cell types.

---

## When to Choose Something Else

### Cell-Level Clustering

Hierarchical clustering with Ward.D2 requires computing a full pairwise distance matrix, which scales as O(n^2) in memory and O(n^2 log n) in time. For cell-level clustering on datasets with tens or hundreds of thousands of cells, this is [computationally infeasible](https://doi.org/10.1016/j.cell.2015.05.047). Use [Leiden clustering](leiden-clustering.md), which builds a sparse KNN graph instead and scales near-linearly.

### Need a Specific Number of Flat Clusters

If you need to partition data into exactly k groups (e.g., for downstream modeling that requires a fixed number of categories), consider:

- **k-means**: Directly optimizes the same variance-minimization objective as Ward but [without the greedy hierarchical constraint](https://doi.org/10.1007/b98835), often finding better partitions for a given k.
- **Cutting the dendrogram** with `cutree(hc, k = k)`: This works but the hierarchical constraint means the partition at level k is [a refinement of the partition at level k-1](https://doi.org/10.1007/b98835), which may not be optimal for the target k.

### Non-Euclidean or Graph-Based Data

Ward's criterion is [defined for squared Euclidean distances only](https://doi.org/10.1007/s00357-014-9161-z). If your data requires a non-Euclidean distance metric (e.g., correlation distance for expression profiles, earth mover's distance for distributions), Ward.D2 is not appropriate. Use average linkage or complete linkage instead, which are [valid with any distance metric](https://doi.org/10.1007/b98835).

### Interactive or Streaming Data

Hierarchical clustering is a batch algorithm: all data must be present before clustering begins, and adding a new item requires recomputing the entire hierarchy. For data arriving in batches or interactive exploration, consider incremental methods.

---

## How to Use in the Pipeline

### Abundance Heatmap with Hierarchical Clustering

```r
# Source required modules
source("src/analysis/abundance_plotting.R")

# Generate all abundance plots including the Ward.D2-biclustered heatmap
plot_abundance_heatmap(
  tracking_matrix = tracking_df,
  output_dir = "output/abundance/"
)
```

The heatmap function computes proportions from the raw count tracking matrix, then applies `hclust(dist(...), method = "ward.D2")` independently to both rows (samples) and columns (clusters) before rendering the heatmap in the reordered layout.

### Proportional Barplot with Dendrogram

```r
source("src/analysis/abundance_plotting.R")

# Proportional barplot with Ward.D2 dendrogram on top
plot_abundance_proportional(
  tracking_matrix = tracking_df,
  output_dir = "output/abundance/"
)
```

This function computes sample-level cluster proportions, performs hierarchical clustering with Ward.D2 to determine sample order, and renders the dendrogram above the barplot using `ggdendro` and `patchwork`.

### Expression Heatmap (pheatmap Internal Clustering)

```r
# Source visualization module
source("src/visualization/heatmaps.R")

# Expression heatmap with automatic hierarchical clustering on both axes
plot_heatmap(
  expr_matrix = median_expr,
  output_path = "output/heatmap.png",
  cluster_rows = TRUE,   # hierarchical clustering of rows
  cluster_cols = TRUE,   # hierarchical clustering of columns
  scale = "row"          # z-score normalization per row
)
```

When `cluster_rows` or `cluster_cols` is `TRUE`, pheatmap internally runs `hclust()` on the corresponding dimension, draws the dendrogram in the margin, and reorders the matrix accordingly.

### Cluster Expression Heatmap

```r
source("src/visualization/heatmaps.R")

# Median expression per cluster, with hierarchical clustering
plot_cluster_heatmap(
  expr_data = combined_expression_df,
  cluster_col = "cluster",
  marker_cols = c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14"),
  output_path = "output/cluster_expression.png",
  title = "Median Marker Expression per Cluster"
)
```

This convenience wrapper computes median expression per cluster, then delegates to `plot_heatmap()` with clustering enabled on both axes.

### Custom Ward.D2 Clustering

For analyses not covered by the built-in plotting functions, you can use R's `hclust()` directly:

```r
# Compute proportions
prop_mat <- counts_mat / rowSums(counts_mat)

# Hierarchical clustering with Ward.D2
hc <- hclust(dist(prop_mat), method = "ward.D2")

# Get dendrogram order
ordered_labels <- hc$labels[hc$order]

# Cut into k groups
groups <- cutree(hc, k = 5)

# Plot dendrogram
plot(hc, hang = -1, cex = 0.6, main = "Sample Dendrogram (Ward.D2)")
```

---

## Distance Metrics

The choice of distance metric determines what "similarity" means for the hierarchical clustering. Different metrics are appropriate for different data types encountered in cytometry analysis.

### Euclidean Distance (Default)

Euclidean distance measures the [straight-line distance between two points in multivariate space](https://doi.org/10.1007/b98835). It is the default for `dist()` in R and the natural companion to Ward's criterion, which is [defined in terms of squared Euclidean distances](https://doi.org/10.1007/s00357-014-9161-z).

```r
d <- dist(prop_mat)  # Euclidean (default)
```

**When to use in cytometry**: Euclidean distance is appropriate when the variables are on comparable scales and the absolute differences between values are meaningful. For abundance proportion matrices (all values between 0 and 1, same units), Euclidean distance [directly reflects differences in cluster composition](https://doi.org/10.12688/f1000research.11622.3). For expression heatmaps where the matrix has been row-scaled (z-scored), Euclidean distance measures differences in standardized expression profiles.

**Limitation**: Euclidean distance is [sensitive to the magnitude of individual variables](https://doi.org/10.1007/b98835). If one cluster dominates (e.g., 80% of cells are CD4+ T cells), differences in that cluster will drive the sample ordering. This is usually desirable for abundance analysis but may be undesirable if you want to emphasize rare populations.

### Correlation Distance

Correlation distance is defined as `1 - r`, where r is the Pearson correlation coefficient between two row vectors. It [measures the similarity of the shape of two profiles regardless of their magnitude or offset](https://doi.org/10.1093/bioinformatics/btw313).

```r
d <- as.dist(1 - cor(t(expr_mat)))  # Correlation distance
```

**When to use in cytometry**: Correlation distance is appropriate for expression heatmaps where you care about which markers are relatively high or low within each cluster, not the absolute intensity. Two clusters with proportionally similar expression patterns but different overall staining intensities will have small correlation distance. [This is the standard distance metric for gene expression heatmaps](https://doi.org/10.1093/bioinformatics/btw313) and transfers well to protein expression per cluster.

**Limitation**: Correlation distance discards information about magnitude. Two samples with identical cluster composition shapes but very different total cell counts will appear identical. Ward's criterion is [not formally defined for correlation distance](https://doi.org/10.1007/s00357-014-9161-z), so use average or complete linkage when using correlation distance.

### Manhattan Distance

Manhattan distance (L1 norm) sums the absolute differences across all dimensions, rather than squaring them as Euclidean distance does. It is [less sensitive to large differences in single variables](https://doi.org/10.1007/b98835) compared to Euclidean distance.

```r
d <- dist(prop_mat, method = "manhattan")
```

**When to use in cytometry**: Manhattan distance can be useful when the abundance matrix contains outlier values in individual clusters (e.g., one sample has an anomalously high proportion of a rare population). Because it sums absolute rather than squared differences, [a single large deviation has less influence on the total distance](https://doi.org/10.1007/b98835) than it would with Euclidean distance. However, Manhattan distance is less commonly used in cytometry than Euclidean or correlation distance.

**Limitation**: Ward's criterion assumes squared Euclidean distances. Using Manhattan distance with Ward.D2 [does not correspond to minimizing within-cluster variance](https://doi.org/10.1007/s00357-014-9161-z) and should be used only with average or complete linkage.

### Recommendation for This Pipeline

| Data Type | Recommended Distance | Rationale |
|-----------|---------------------|-----------|
| Abundance proportions (samples x clusters) | **Euclidean** | Directly measures composition differences; compatible with Ward.D2 |
| Row-scaled expression (clusters x markers) | **Euclidean** on z-scores | After scaling, Euclidean distance captures profile shape on standardized data |
| Unscaled expression (clusters x markers) | **Correlation** | Removes intensity bias; focuses on relative marker patterns |
| Presence/absence gating data | **Euclidean** or **binary** | Euclidean on 0/1 data equals square root of Hamming distance |

---

## Linkage Methods

The linkage method defines how the distance between two multi-member clusters is computed from the pairwise distances between their members. The choice of linkage [has a profound effect on the shape of the resulting dendrogram](https://doi.org/10.1007/b98835) and which structures it reveals.

### Ward.D2

Ward.D2 minimizes the [total within-cluster variance at each merge](https://doi.org/10.1007/s00357-014-9161-z), as described above. It tends to produce compact, roughly equally-sized clusters and is the [most commonly used linkage method for heatmap visualization in genomics and cytometry](https://doi.org/10.1093/bioinformatics/btw313).

**Strengths**: Produces balanced, compact dendrograms. Merges reflect meaningful increases in heterogeneity. Good default choice when combined with Euclidean distance.

**Weaknesses**: Requires squared Euclidean distance (use `"ward.D2"` to handle this automatically). [Tends to impose spherical cluster shapes](https://doi.org/10.1007/b98835), which may not reflect the true structure when clusters have elongated or irregular geometry.

### Complete Linkage

Complete linkage defines the inter-cluster distance as the [maximum distance between any pair of members, one from each cluster](https://doi.org/10.1007/b98835). It is also called "farthest neighbor" linkage.

```r
hc <- hclust(dist(mat), method = "complete")
```

**Strengths**: Produces tight, compact clusters. Avoids the chaining problem seen with single linkage. [Works with any distance metric](https://doi.org/10.1007/b98835), not just Euclidean. This is the default method in `pheatmap`.

**Weaknesses**: Sensitive to outliers -- a single distant point in a cluster can inflate the inter-cluster distance, [preventing merges that would otherwise be appropriate](https://doi.org/10.1007/b98835).

**When to use**: When you need a valid clustering with non-Euclidean distances (e.g., correlation distance) and want compact clusters. Good alternative to Ward.D2 when the distance metric is not Euclidean.

### Average Linkage (UPGMA)

Average linkage computes the inter-cluster distance as the [arithmetic mean of all pairwise distances between members of the two clusters](https://doi.org/10.1007/b98835). This is the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) approach.

```r
hc <- hclust(dist(mat), method = "average")
```

**Strengths**: [Robust to outliers compared to both complete and single linkage](https://doi.org/10.1007/b98835) because no single pairwise distance dominates. Works with any distance metric. Produces dendrograms where the merge height [reflects the average dissimilarity between merged groups](https://doi.org/10.1007/b98835), which has a natural interpretation.

**Weaknesses**: Can produce clusters of unequal compactness. [May not recover well-separated compact clusters as cleanly as Ward or complete linkage](https://doi.org/10.1007/b98835).

**When to use**: When working with non-Euclidean distances and outlier robustness is important. Commonly used in [phylogenetic analysis and gene expression studies](https://doi.org/10.1093/bioinformatics/btw313) where the ultrametric assumption (constant evolutionary rate) approximately holds.

### Single Linkage

Single linkage defines the inter-cluster distance as the [minimum distance between any pair of members](https://doi.org/10.1007/b98835). It is also called "nearest neighbor" linkage.

```r
hc <- hclust(dist(mat), method = "single")
```

**Strengths**: Can detect elongated, non-convex cluster shapes. [Corresponds to the minimum spanning tree of the distance matrix](https://doi.org/10.1007/b98835), giving it a graph-theoretic interpretation.

**Weaknesses**: Highly susceptible to the ["chaining" effect](https://doi.org/10.1007/b98835), where a sequence of close intermediate points causes two dissimilar clusters to merge prematurely through a bridge of intermediate items. This produces long, straggling dendrograms that are [rarely useful for heatmap ordering](https://doi.org/10.1093/bioinformatics/btw313).

**When to use**: Rarely appropriate for cytometry heatmaps. May be useful when investigating whether two populations are connected by transitional cells (e.g., developmental trajectories), but [Leiden clustering or trajectory-inference methods are generally more appropriate for this task](https://doi.org/10.1016/j.cell.2015.05.047).

### Comparison Table

| Linkage | Cluster Shape | Outlier Sensitivity | Requires Euclidean | Heatmap Suitability |
|---------|--------------|--------------------|--------------------|---------------------|
| **Ward.D2** | Compact, balanced | Moderate | [Yes](https://doi.org/10.1007/s00357-014-9161-z) | Excellent (default in this pipeline) |
| **Complete** | Compact, tight | High | No | Good (pheatmap default) |
| **Average** | Variable | Low | No | Good |
| **Single** | Elongated, chained | Very high | No | [Poor](https://doi.org/10.1093/bioinformatics/btw313) |

For the abundance proportion matrices and z-scored expression matrices used in this pipeline, **Ward.D2 with Euclidean distance** is the recommended combination. It produces [well-balanced dendrograms where merge heights are directly interpretable as increases in within-cluster heterogeneity](https://doi.org/10.1007/s00357-014-9161-z), and it groups samples and clusters into compact, visually coherent blocks in the heatmap.

---

## Dendrogram Interpretation

The dendrogram produced by hierarchical clustering encodes the complete merge history. Understanding how to read it is essential for interpreting the heatmaps and barplots in this pipeline.

### Merge Heights

Each internal node in the dendrogram has a height corresponding to the dissimilarity at which the two child clusters were merged. With Ward.D2, this height equals the [square root of the increase in total within-cluster variance caused by the merge](https://doi.org/10.1007/s00357-014-9161-z). Large jumps in merge height indicate that the clusters being merged were substantially different -- these jumps correspond to [natural cluster boundaries in the data](https://doi.org/10.1007/b98835).

### Leaf Ordering

The dendrogram constrains which items can be adjacent (they must be in the same subtree at some level), but within those constraints, the left-right ordering of branches at each internal node is arbitrary. R's `hclust()` uses a [deterministic ordering algorithm that places the subtree with the smaller first element on the left](https://doi.org/10.1007/b98835), but this is a convention, not a property of the clustering. When interpreting heatmaps, focus on the tree structure (which items group together) rather than the specific left-right arrangement.

### Cophenetic Correlation

The cophenetic correlation coefficient measures [how faithfully the dendrogram preserves the original pairwise distances](https://doi.org/10.1093/sysbio/12.1.65). It is the Pearson correlation between the original distance matrix and the cophenetic distance matrix (where the cophenetic distance between two items is the merge height at which they first join the same cluster). Values above 0.75 indicate a [good representation of the distance structure](https://doi.org/10.1007/b98835); values below 0.6 suggest that the hierarchical constraint is distorting the relationships. Ward's method typically produces [moderate cophenetic correlations because it optimizes variance rather than distance preservation](https://doi.org/10.1007/s00357-014-9161-z), but for visualization purposes this is rarely a concern because the goal is grouping, not exact distance reproduction.

---

## Related Documentation

- [Leiden Clustering](leiden-clustering.md) -- Primary cell-level clustering method in this pipeline
- [UMAP](umap.md) -- Dimensionality reduction for cell-level visualization
- [PCA](pca.md) -- Linear dimensionality reduction and variance decomposition
- [Batch Correction](batch-correction.md) -- Correcting technical variation before clustering
- [Composability Guide](../composability-guide.md) -- How to combine pipeline modules
