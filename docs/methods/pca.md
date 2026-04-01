# Principal Component Analysis (PCA) in Flow Cytometry

## Summary

Principal Component Analysis is a linear dimensionality reduction technique that
re-expresses a high-dimensional marker expression matrix as a set of orthogonal
axes -- principal components -- ranked by the amount of variance they capture.
In this pipeline, `run_pca()` (in `src/R_shared/stats_utils.R`) wraps R's
`prcomp()` with centering and optional scaling, returning scores, variance
explained, and a loadings (rotation) matrix. The companion `plot_pca()` (in
`src/visualization/pca.R`) produces scatter and loadings plots suitable for
batch-effect QC, marker-importance screening, and pre-clustering noise
reduction.

---

## The Method

### Eigendecomposition and Singular Value Decomposition

PCA finds the directions in which a dataset varies most. Given an *n x p* matrix
**X** (cells x markers), PCA operates on the [covariance (or correlation) matrix
of the centered data](https://doi.org/10.1007/b98835) and extracts its
eigenvectors, which become the principal component axes. In practice, R's
`prcomp()` computes PCA via the [singular value decomposition
(SVD)](https://doi.org/10.1007/b98835) of the centered (and optionally scaled)
data matrix rather than forming the covariance matrix explicitly, because SVD is
[numerically more stable](https://doi.org/10.1007/b98835) for high-dimensional
data.

The SVD factorizes the centered data matrix as:

```
X = U * D * V'
```

where **U** contains the left singular vectors (cell scores, up to a scaling),
**D** is a diagonal matrix of singular values, and **V** contains the right
singular vectors (marker loadings). The principal component scores returned by
`prcomp()` are **X * V**, and the proportion of variance explained by each
component is *d_k^2 / sum(d^2)*.

### Variance Maximization

Each successive principal component is defined as the [linear combination of the
original variables that maximizes the variance of the projected
data](https://doi.org/10.1007/b98835), subject to orthogonality with all
preceding components. This means PC1 captures the single direction of greatest
spread in the marker space, PC2 captures the greatest spread orthogonal to PC1,
and so on.

### Orthogonality

All principal components are mutually orthogonal by construction. This
[guarantees zero correlation between any pair of
PCs](https://doi.org/10.1007/b98835), which is useful both for interpretation
(each PC represents an independent axis of variation) and for downstream methods
that assume uncorrelated inputs (e.g., certain clustering algorithms and
regression models).

---

## Why PCA for Cytometry

Flow and mass cytometry panels routinely measure 10 to 50 protein markers per
cell. Analyzing this many dimensions simultaneously is both computationally
expensive and difficult to interpret. PCA addresses this by compressing the
marker space into a small number of components that retain [most of the
biologically relevant variance](https://doi.org/10.1016/j.cell.2015.05.047).

In their landmark CyTOF study, [Levine et al. (2015)](https://doi.org/10.1016/j.cell.2015.05.047)
used PCA as a preprocessing step before phenograph clustering to reduce
dimensionality and filter noise from mass cytometry expression matrices. This
approach -- reducing to the top *k* PCs before graph-based clustering -- has
since become a [standard practice in high-dimensional cytometry
workflows](https://doi.org/10.12688/f1000research.11622.3).

PCA is also the method of choice for [detecting batch effects in multi-site or
multi-day cytometry experiments](https://doi.org/10.1038/s41592-019-0619-0).
When samples cluster by acquisition batch rather than by biological condition on
a PC1-vs-PC2 plot, this is a strong signal that batch correction (e.g., Harmony)
is needed before downstream analysis.

The [Nowicka et al. (2017) Bioconductor workflow for CyTOF
data](https://doi.org/10.12688/f1000research.11622.3) recommends PCA on
sample-level summary statistics (e.g., median marker expression per cluster per
sample) as a tool for exploratory data analysis and for detecting outlier
samples before differential analysis.

Beyond cytometry, PCA is foundational in single-cell transcriptomics. The
[Scanpy](https://doi.org/10.1186/s13059-017-1382-0) and
[Seurat](https://doi.org/10.1016/j.cell.2019.05.031) frameworks both use PCA as
the default dimensionality reduction step before nearest-neighbor graph
construction, and the same logic applies to high-parameter cytometry data.

---

## Centering and Scaling

### Centering (`center = TRUE`)

Centering subtracts the column mean from each marker, shifting the data cloud so
that its centroid sits at the origin. Without centering, the [first principal
component would point toward the mean of the data rather than in the direction
of maximum variance](https://doi.org/10.1007/b98835). Centering is almost
always required and is enabled by default in both `prcomp()` and `run_pca()`.

### Scaling (`scale. = TRUE`)

Scaling divides each marker by its standard deviation after centering, producing
[unit-variance variables](https://doi.org/10.1007/b98835). This is equivalent
to performing PCA on the correlation matrix rather than the covariance matrix.

**When to scale (the pipeline default).** Flow and mass cytometry markers often
span different dynamic ranges. A marker measured on a 0-to-262144 scale (e.g.,
a bright fluorochrome) will dominate the variance over a marker measured on a
0-to-1000 scale (e.g., a dim marker) if the data is not scaled. Scaling ensures
that [each marker contributes equally to the analysis regardless of its
measurement range](https://doi.org/10.1007/b98835). This is the recommended
default for multi-color panels where markers have heterogeneous dynamic ranges.

**When not to scale.** If all markers have already been transformed to
comparable scales (e.g., arcsinh-transformed CyTOF data where all values
typically fall within a similar range), scaling may not be necessary and can
sometimes [mask real differences in biological
variability](https://doi.org/10.1007/b98835) between markers. If the relative
variance of markers is itself informative (e.g., you want PCA to highlight
highly variable markers), use `scale = FALSE`.

In this pipeline, `run_pca()` defaults to `center = TRUE, scale. = TRUE`:

```r
run_pca <- function(expr_matrix, n_components = 2, scale = TRUE) {
  pca_result <- prcomp(expr_matrix, center = TRUE, scale. = scale,
                       rank. = n_components)
  ...
}
```

Override with `scale = FALSE` when working with pre-normalized data where
differential variance is meaningful.

---

## When Your Data Calls for It

### 1. Understanding the Major Axes of Variation

PCA reveals which directions in marker space carry the most information. A
PC1-vs-PC2 scatter plot provides an immediate summary of the [dominant patterns
in the data](https://doi.org/10.1007/b98835). If PC1 separates T cells from
myeloid cells and PC2 separates activated from resting states, this tells you
the two largest sources of variation are lineage and activation status.

### 2. Batch Effect QC (Before and After Correction)

Plot samples on PC1 vs. PC2 colored by batch. If batch is the primary source of
variation, samples will separate by batch along the early PCs. After applying
batch correction (e.g., [Harmony](https://doi.org/10.1038/s41592-019-0619-0)),
re-run PCA to confirm that [samples now intermingle across
batches](https://doi.org/10.12688/f1000research.11622.3) while retaining
biological separation.

This before-and-after PCA comparison is the standard QC step for batch
correction in cytometry, as recommended by the [Nowicka et al. (2017)
workflow](https://doi.org/10.12688/f1000research.11622.3) and the
[Harmony publication](https://doi.org/10.1038/s41592-019-0619-0).

### 3. Pre-Clustering Noise Reduction (Top N PCs)

High-dimensional expression data contains both signal and measurement noise.
Retaining only the top *k* PCs before clustering [filters out low-variance noise
dimensions that can confuse nearest-neighbor graph
construction](https://doi.org/10.1016/j.cell.2015.05.047). This is the
approach used by PhenoGraph and adopted in single-cell RNA-seq workflows by
[Seurat](https://doi.org/10.1016/j.cell.2019.05.031) and
[Scanpy](https://doi.org/10.1186/s13059-017-1382-0). In this pipeline, you can
cluster on PCA scores:

```r
pca_result <- run_pca(expr_matrix, n_components = 20)
# Feed pca_result$coordinates to Leiden clustering
```

### 4. Loadings Plots to Identify Driving Markers

The loadings matrix (rotation matrix from `prcomp()`) tells you [how much each
original marker contributes to each principal
component](https://doi.org/10.1007/b98835). Markers with large absolute
loadings on a given PC are the primary drivers of that axis of variation.
`plot_pca()` produces a loadings plot where each marker is drawn as an arrow in
PC1-vs-PC2 space. Long arrows far from the origin indicate markers with strong
influence. Arrows pointing in similar directions indicate [correlated
markers](https://doi.org/10.1007/b98835) (e.g., CD3 and CD5 may co-vary as
T-cell lineage markers).

### 5. Outlier Sample Detection

Samples that fall far from the main cluster in PCA space may represent
[technical failures, sample degradation, or acquisition
artifacts](https://doi.org/10.12688/f1000research.11622.3). PCA on
sample-level summary statistics (median marker expression per sample) is
particularly useful for this: an outlier sample will appear as an isolated
point on the PCA plot, flagging it for manual review before it contaminates
downstream analysis.

---

## When to Choose Something Else

PCA is a linear method. When your analysis requires capturing nonlinear
relationships or detecting rare populations, other methods in the pipeline are
more appropriate.

- **Nonlinear manifold structure.** If the relationships between cell
  populations are nonlinear -- populations form arcs, branches, or complex
  topologies in marker space -- PCA will fail to resolve them because it can
  only capture linear axes. Use [UMAP](umap.md) instead, which preserves local
  and global nonlinear structure.

- **Rare populations.** PCA [maximizes total variance and is therefore dominated
  by the most abundant populations](https://doi.org/10.1007/b98835). A rare
  subset comprising 0.1% of cells contributes negligibly to the total variance
  and will not appear as a distinct group on a PCA plot. Use [UMAP](umap.md)
  for visualization when rare populations matter.

- **Clustering.** Do not cluster directly on PC1 and PC2 alone. While PCA can
  serve as a preprocessing step (retaining the top *k* PCs as input features),
  the actual clustering should be performed using a graph-based method such as
  [Leiden](leiden-clustering.md), which is designed to partition a
  nearest-neighbor graph into communities. PCA coordinates can be used as
  *input* to Leiden (via a KNN graph on the top PCs), but running k-means or
  hierarchical clustering on just two PCs discards the information in the
  remaining components.

---

## How to Use in the Pipeline

### From an Expression Matrix (Any Pipeline)

The most common entry point is `plot_pca()` in `src/visualization/pca.R`, which
calls `run_pca()` internally, generates scatter and loadings plots, and saves
them to disk.

```r
source("src/visualization/pca.R")

result <- plot_pca(
  expr_matrix    = combined_matrix,   # numeric matrix: rows = cells, cols = markers
  sample_labels  = sample_labels,     # character vector, one per row
  output_dir     = "output/pca/",
  n_components   = 10,                # compute 10 PCs
  color_by       = "sample",          # or pass a character vector of group labels
  title          = "PCA — All Samples"
)

# Access results programmatically
result$pca_result$variance_explained   # numeric vector
result$pca_result$loadings             # rotation matrix
result$scatter_plot                     # ggplot object
result$loadings_plot                    # ggplot object
```

### Pre-Reduce Before Leiden Clustering

Use `run_pca()` from `src/R_shared/stats_utils.R` to reduce noise before
building the nearest-neighbor graph for Leiden clustering.

```r
source("src/R_shared/stats_utils.R")

pca_result <- run_pca(expr_matrix, n_components = 20, scale = TRUE)

# pca_result$coordinates is an n x 20 matrix of PC scores
# Feed this to Leiden clustering instead of the full marker matrix
source("src/clustering/leiden.R")
clusters <- run_leiden(pca_result$coordinates, resolution = 1.0)
```

### Batch-Effect QC Workflow

```r
source("src/visualization/pca.R")

# Before batch correction
plot_pca(expr_matrix_raw, sample_labels = batch_labels,
         output_dir = "output/qc/pca_before/",
         color_by = batch_labels, title = "PCA — Before Correction")

# After Harmony batch correction
source("src/harmonization/batch_correction.R")
corrected <- run_batch_correction(expr_matrix_raw, batch = batch_labels)

plot_pca(corrected, sample_labels = batch_labels,
         output_dir = "output/qc/pca_after/",
         color_by = batch_labels, title = "PCA — After Correction")
```

### Sample-Level Outlier Detection

```r
source("src/R_shared/stats_utils.R")

# Compute median expression per sample
sample_medians <- aggregate(expr_matrix, by = list(sample = sample_ids), FUN = median)
rownames(sample_medians) <- sample_medians$sample
sample_medians$sample <- NULL

pca_samples <- run_pca(as.matrix(sample_medians), n_components = 5)
# Inspect pca_samples$coordinates for isolated points
```

---

## Interpreting Results

### Scree Plot (Variance Explained)

A scree plot shows the [proportion of variance explained by each successive
PC](https://doi.org/10.1007/b98835). Components before the "elbow" (the point
where the curve flattens) capture structured signal; components after the elbow
capture noise. The pipeline returns `variance_explained` as a numeric vector,
which you can plot:

```r
pca <- run_pca(expr_matrix, n_components = 10)
barplot(pca$variance_explained * 100,
        names.arg = paste0("PC", 1:10),
        ylab = "% Variance Explained",
        main = "Scree Plot")
```

### Loadings (Rotation Matrix)

The loadings matrix has one row per marker and one column per PC. Each entry is
the [weight (coefficient) of that marker in the linear combination defining the
PC](https://doi.org/10.1007/b98835). `plot_pca()` visualizes the first two
columns as a biplot-style arrow diagram, where:

- **Arrow direction** indicates whether the marker is positively or negatively
  associated with each PC.
- **Arrow length** indicates the [strength of the marker's contribution to the
  PC subspace](https://doi.org/10.1007/b98835).
- **Angles between arrows** approximate the [correlation structure among
  markers](https://doi.org/10.1007/b98835): small angles mean positive
  correlation, right angles mean no correlation, and opposite directions mean
  negative correlation.

For cytometry data, loadings are particularly useful for identifying which
surface markers or intracellular stains drive the separation between
experimental groups, as demonstrated in the [Nowicka et al. (2017)
workflow](https://doi.org/10.12688/f1000research.11622.3).

### Biplot

A biplot superimposes the scores (cell or sample positions) and loadings
(marker arrows) on the same plot, providing a [simultaneous visualization of
how observations relate to each other and which variables drive those
relationships](https://doi.org/10.1007/b98835). While `plot_pca()` produces
scores and loadings as separate plots by default, you can construct a biplot
manually from the returned objects if needed.

---

## Variance Explained: How to Decide the Number of Components

Choosing *k*, the number of principal components to retain, is a critical
decision when using PCA as a preprocessing step. Several criteria exist:

### Elbow Rule (Scree Plot)

Inspect the scree plot and select the number of components at the [point where
the curve transitions from steep descent to a
plateau](https://doi.org/10.1007/b98835). Components beyond the elbow explain
diminishing amounts of variance and are increasingly dominated by noise. This is
the most commonly used heuristic in cytometry workflows.

### Kaiser Criterion

Retain only components with [eigenvalues greater than
1](https://doi.org/10.1007/b98835) (equivalently, components that explain more
variance than a single original standardized variable). This rule applies when
PCA is computed on the correlation matrix (i.e., with `scale. = TRUE`). While
simple, it can [over- or under-estimate the true
dimensionality](https://doi.org/10.1007/b98835) depending on the correlation
structure.

### Cumulative Variance Threshold

Retain enough components to explain a target percentage of total variance (e.g.,
80% or 90%). Compute the cumulative sum of `variance_explained` and find the
smallest *k* that crosses the threshold:

```r
cumvar <- cumsum(pca$variance_explained)
k <- which(cumvar >= 0.90)[1]
```

### Practical Recommendations for Cytometry

For [CyTOF panels with 30-40 markers, retaining 15-30 PCs before clustering is
typical](https://doi.org/10.1016/j.cell.2015.05.047). For smaller fluorescence
panels (10-15 markers), fewer components (5-10) usually suffice. When in doubt,
err on the side of retaining more components -- downstream methods like Leiden
clustering are robust to a few extra noisy dimensions, but dropping informative
PCs can [merge biologically distinct populations that differ along
low-variance axes](https://doi.org/10.1016/j.cell.2019.05.031).

---

## Limitations

### Linearity

PCA can only capture [linear relationships between
markers](https://doi.org/10.1007/b98835). If two cell populations are separated
by a curved boundary in marker space, PCA may project them onto overlapping
regions. Nonlinear methods such as [UMAP](umap.md) or t-SNE handle this better,
though at the cost of interpretability and distance preservation.

### Dominated by High-Variance Markers (if Unscaled)

When `scale. = FALSE`, markers with the [largest absolute variance dominate the
first principal components](https://doi.org/10.1007/b98835), regardless of
whether that variance is biologically interesting. A bright fluorochrome with
high spread may monopolize PC1 even if a dim marker carries more discriminating
information. Always consider whether scaling is appropriate for your panel.

### Sensitive to Outliers

Because PCA [maximizes variance, it is sensitive to extreme
values](https://doi.org/10.1007/b98835). A small number of outlier cells with
aberrant expression can distort the principal components. Pre-filtering outlier
events (e.g., removing doublets and debris via the preprocessing pipeline)
mitigates this issue.

### Assumes Multivariate Normality (Approximately)

While PCA does not strictly require normality, its [optimality properties
(maximizing explained variance is equivalent to maximizing likelihood) hold
under multivariate Gaussian
assumptions](https://doi.org/10.1007/b98835). Heavily skewed or
multi-modal marker distributions -- common in cytometry before transformation --
can produce PCs that do not correspond to biologically meaningful axes. Applying
an appropriate transform (arcsinh for CyTOF, logicle for fluorescence) before
PCA improves results, as recommended by [Nowicka et al.
(2017)](https://doi.org/10.12688/f1000research.11622.3).

### Does Not Preserve Local Neighborhood Structure

PCA preserves [global variance structure but does not guarantee that nearby
points in the original space remain nearby in the reduced
space](https://doi.org/10.1038/nbt.4314). Methods like UMAP, which explicitly
optimize for neighborhood preservation, are better suited for visualizing local
cluster structure.

### Rare Populations Are Invisible

Because PCA is [variance-maximizing](https://doi.org/10.1007/b98835), rare cell
types that represent a tiny fraction of total cells contribute almost nothing to
the principal components. These populations may be projected into the dense core
of the plot and become invisible. For rare-event detection, prefer [UMAP](umap.md)
visualization combined with targeted gating or density-based methods.
