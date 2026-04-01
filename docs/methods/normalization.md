# Normalization Methods: Z-Score Scaling and Quantile Normalization

## Summary

Normalization adjusts expression values so that markers and samples can be
compared on a common scale, removing technical variation in measurement
magnitude while preserving biological differences. This pipeline provides two
normalization methods in `src/preprocessing/normalization.R`: **z-score
scaling**, which centers each marker to mean zero and unit variance within a
sample, and **quantile normalization**, which forces the empirical distribution
of every marker (or every sample) to be identical. Harmony batch correction is
a distinct operation with its own page at [Batch Correction (Harmony)](batch-correction.md).

---

## Z-Score Scaling

### What It Does

Z-score scaling transforms each marker channel so that the resulting values
have a mean of zero and a standard deviation of one. For a marker *x* measured
across *n* cells, each value is replaced by:

$$
z_i = \frac{x_i - \bar{x}}{s}
$$

where $\bar{x}$ is the sample mean and $s$ is the sample standard deviation.
This is also called **standardization** and is the transformation performed by
R's built-in `scale()` function. The concept is a direct application of the
[standard score from elementary statistics](https://doi.org/10.1007/978-0-387-32833-1),
where any measurement is expressed as the number of standard deviations it
falls above or below the mean.

### Why It Matters for Cytometry

Flow and mass cytometry panels measure proteins on fundamentally different
intensity scales. A surface marker detected by a bright fluorochrome may
produce values in the tens of thousands, while a dim intracellular marker on
the same panel may peak at a few hundred. When these raw or transformed values
are fed directly into distance-based algorithms (PCA, Leiden clustering, UMAP),
the high-intensity markers
[dominate the distance computation and drive all downstream structure](https://doi.org/10.1007/b98835),
regardless of their biological informativeness. Z-score scaling removes this
scale dependence by putting every marker on a common unitless axis.

This is the same principle behind the recommendation to use the correlation
matrix rather than the covariance matrix for PCA when variables have
[heterogeneous scales](https://doi.org/10.1007/b98835). In single-cell
genomics, analogous per-feature scaling is a
[standard preprocessing step before PCA and graph-based clustering](https://doi.org/10.1186/s13059-017-1382-0)
in workflows such as Scanpy and Seurat.

### Properties

- **Per-marker, within-sample.** Each marker is scaled independently using
  statistics computed from the cells in that sample. Cross-sample alignment is
  not performed.
- **Linear transformation.** Z-score scaling is a shift-and-scale operation
  that preserves the shape of the distribution. If a marker's expression is
  bimodal before scaling, it remains bimodal after scaling; only the axis labels
  change.
- **Reversible.** Given the original mean and standard deviation, the raw
  values can be recovered exactly.
- **Sensitive to outliers.** Because the mean and standard deviation are
  influenced by extreme values, a small number of aberrant cells can inflate
  *s* and compress the z-scores of the remaining population.
  [Robust alternatives](https://doi.org/10.1007/978-0-387-32833-1) such as
  median-and-MAD scaling exist but are not currently implemented in this
  pipeline.

### Mathematical Detail

The `z_score_scaling()` implementation uses `flowCore::transformList()` with
R's `scale` function applied to each selected marker column. Internally,
`scale(x)` computes:

```
z = (x - mean(x)) / sd(x)
```

where `sd()` uses the sample standard deviation with Bessel's correction
(dividing by *n - 1*). For a flowFrame with 50,000 cells and 15 markers, this
produces a 50,000 x 15 matrix where each column has mean approximately zero
and standard deviation exactly one. The operation runs in O(*n*) time per
marker and is
[numerically stable for typical cytometry sample sizes](https://doi.org/10.1007/b98835).

---

## Quantile Normalization

### What It Does

Quantile normalization forces all markers (or all samples) to have the same
empirical distribution. After quantile normalization, the histogram of values
for every normalized column is identical. This is a much stronger adjustment
than z-score scaling: rather than matching only the first two moments (mean and
variance), it matches the entire distribution function.

### Origin and Algorithm

Quantile normalization was originally developed for
[cross-array normalization of Affymetrix microarray data by Bolstad et al. (2003)](https://doi.org/10.1093/bioinformatics/19.1.185).
The motivation was that technical variation between microarray chips caused
systematic differences in intensity distributions that obscured biological
signal. By forcing all arrays to share the same distribution, these technical
differences were removed while preserving the relative ranking of genes within
each array.

The algorithm proceeds in three steps:

1. **Rank.** For each column (marker or sample), sort the values and record
   their ranks. Ties are handled by averaging.

2. **Compute reference distribution.** For each rank position *k*, compute the
   mean (or median) of the *k*-th smallest value across all columns. This
   produces a single reference distribution that represents the "average"
   distribution shape.

3. **Replace.** For each column, replace the value at rank *k* with the
   reference value at rank *k*. The original rank order within each column is
   preserved; only the spacing between values changes.

Bolstad et al. (2003) showed that this procedure
[outperformed other normalization methods for microarray data in terms of bias and variance reduction](https://doi.org/10.1093/bioinformatics/19.1.185),
and it subsequently became the de facto standard for microarray preprocessing.
It is implemented in the Bioconductor
[preprocessCore](https://doi.org/10.18129/B9.bioc.preprocessCore) package,
which is what this pipeline calls via `preprocessCore::normalize.quantiles()`.

### Adaptation to Cytometry

In flow and mass cytometry, quantile normalization is applied when samples
acquired on different days or instruments exhibit
[systematic shifts in marker intensity distributions](https://doi.org/10.1002/cyto.a.22271)
that are not attributable to biology. For example, if two batches of samples
were stained with the same antibody panel but exhibit different median
fluorescence intensities due to laser power drift or lot-to-lot antibody
variation, quantile normalization aligns these distributions so that
downstream clustering and statistical analysis are not confounded by the
technical shift.

Finak et al. (2014) applied quantile normalization as part of a
[standardized preprocessing pipeline for flow cytometry data](https://doi.org/10.1186/1471-2105-15-S10-S2)
in the FlowCAP challenges, demonstrating that distribution alignment improved
the consistency of automated gating across sites. In CyTOF specifically,
[Finck et al. (2013)](https://doi.org/10.1002/cyto.a.22271) developed bead-based
normalization to correct for signal decay over acquisition time; quantile
normalization provides a complementary approach when bead standards are
unavailable or when residual distributional differences remain after bead
correction.

### Properties

- **Cross-column alignment.** Unlike z-score scaling (which operates within
  each column independently), quantile normalization uses information from all
  columns simultaneously to construct the reference distribution.
- **Preserves rank order.** Within each column, the relative ordering of
  values is maintained. If cell A had higher CD4 expression than cell B before
  normalization, this remains true afterward.
- **Destroys distributional differences.** By design, all normalized columns
  share the same histogram shape. If two markers genuinely have different
  distribution shapes (e.g., one is bimodal and the other is unimodal), this
  biological difference is
  [removed by quantile normalization](https://doi.org/10.1093/biostatistics/kxj037).
  Hicks and Irizarry (2015) showed that when groups have
  [genuinely different underlying distributions, quantile normalization can introduce false signals](https://doi.org/10.1093/biostatistics/kxv038)
  by forcing distributions to match when they should not.
- **Not reversible.** The original distribution shapes cannot be recovered from
  quantile-normalized data.

### Mathematical Detail

Given a data matrix **X** with *n* rows (cells) and *p* columns (markers), the
`quantile_norm()` implementation:

1. Extracts the expression matrix via `flowCore::exprs(flowframe)`.
2. Passes the selected marker columns to
   `preprocessCore::normalize.quantiles()`, which implements the
   [Bolstad et al. (2003)](https://doi.org/10.1093/bioinformatics/19.1.185)
   algorithm.
3. Writes the normalized values back into the flowFrame.

The computational cost is O(*np* log *n*) due to the sorting step for each of
the *p* columns.

---

## When Your Data Calls For It

### Z-Score Scaling

Use z-score scaling when:

- **Before PCA.** PCA on unscaled data with heterogeneous marker ranges will
  produce principal components dominated by the highest-variance markers.
  Scaling ensures that
  [each marker contributes proportionally to the variance decomposition](https://doi.org/10.1007/b98835),
  which is equivalent to computing PCA on the correlation matrix. The
  [Jolliffe (2002)](https://doi.org/10.1007/b98835) reference text on PCA
  recommends scaling when variables have different units or scales.

- **Before clustering when markers have different scales.** Leiden clustering
  builds a nearest-neighbor graph from Euclidean distances. If one marker spans
  0--10 and another spans 0--10,000 after transformation, the second marker
  will dominate distances and drive all cluster boundaries. Z-score scaling
  equalizes marker contributions. The
  [Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
  applies this scaling before constructing the SNN graph for clustering.

- **Before distance-based methods generally.** Any algorithm that computes
  inter-cell distances (UMAP, k-means, hierarchical clustering) will be
  affected by differences in marker scale. Z-score scaling is the simplest
  correction. Becht et al. (2019) noted that
  [appropriate preprocessing including scaling is essential for producing UMAP embeddings that reflect biology rather than measurement artifacts](https://doi.org/10.1038/nbt.4314).

- **When comparing marker contributions.** After z-score scaling, a value of
  +2 on any marker means "two standard deviations above the mean," enabling
  direct comparison of expression levels across markers with different original
  scales.

### Quantile Normalization

Use quantile normalization when:

- **Cross-sample distribution alignment.** You have multiple samples that
  should have comparable intensity distributions (e.g., technical replicates,
  or samples stained with the same panel) but show systematic distributional
  shifts due to instrument drift, reagent lot variation, or staining batch
  effects. Quantile normalization
  [forces these distributions to match](https://doi.org/10.1093/bioinformatics/19.1.185).

- **Multi-site studies without bead standards.** When combining cytometry data
  from different instruments or sites, and bead-based normalization
  [(Finck et al. 2013)](https://doi.org/10.1002/cyto.a.22271) was not
  performed at acquisition time, quantile normalization can retrospectively
  align intensity distributions.

- **After transformation, before cross-sample comparison.** Quantile
  normalization is most effective when applied to already-transformed data
  (arcsinh or logicle), so that the distributions being aligned are on a
  biologically meaningful scale.

---

## When to Choose Something Else

Not every situation calls for z-score scaling or quantile normalization. The
table below maps common scenarios to better alternatives.

| Situation | Problem with Normalization | Better Choice |
|-----------|---------------------------|---------------|
| Batch effects from different instruments or acquisition days | Z-score scaling only standardizes per-marker statistics; it does not align samples across batches. Quantile normalization is too aggressive -- it removes all distributional differences, including real ones. | [Batch Correction (Harmony)](batch-correction.md), which [removes batch-specific shifts while preserving biological variation](https://doi.org/10.1038/s41592-019-0619-0) |
| Raw intensity data that spans multiple orders of magnitude | Normalization does not change the fundamental shape of raw intensity distributions. The data needs to be put on a biologically meaningful scale first. | Apply [Arcsinh](arcsinh-transform.md) or [Logicle](logicle-transform.md) transform before any normalization |
| Before differential abundance or differential expression testing | Normalization removes distributional differences between groups. If those differences are the biological signal you are trying to detect, normalization will [eliminate or attenuate real effects](https://doi.org/10.1093/biostatistics/kxv038). | Do not normalize. Perform statistical tests on transformed (but un-normalized) data. Account for batch as a covariate in the model if needed ([negative binomial GLMM](negative-binomial-glmm.md)). |
| Correcting for CyTOF signal decay within a run | Within-run signal drift is a specific instrument artifact that requires time-dependent correction, not global distribution alignment. | [Bead-based normalization](https://doi.org/10.1002/cyto.a.22271) at acquisition time, as described by Finck et al. (2013) |
| Data from a single sample with consistent staining | If there is only one sample and no distributional mismatch to correct, normalization adds complexity without benefit. | Z-score scaling may still be useful before PCA/clustering. Quantile normalization is not applicable (requires multiple columns to align). |

---

## Z-Score vs Quantile vs Harmony

The three normalization-related methods in this pipeline operate at different
levels and make different assumptions. The table below summarizes their key
differences.

| Property | Z-Score Scaling | Quantile Normalization | Harmony Batch Correction |
|----------|----------------|----------------------|--------------------------|
| **Scope** | Per-marker, within-sample | Cross-marker or cross-sample | Cross-sample, batch-aware |
| **What it matches** | Mean (0) and SD (1) per marker | Entire empirical distribution across columns | Cell-level embeddings across batches |
| **What it preserves** | Distribution shape, rank order | Rank order only | [Biological variation while removing batch shifts](https://doi.org/10.1038/s41592-019-0619-0) |
| **What it removes** | Differences in marker scale and location | All distributional differences between columns | [Batch-specific systematic offsets](https://doi.org/10.1038/s41592-019-0619-0) |
| **Requires batch labels** | No | No | Yes |
| **Typical use case** | Before PCA, UMAP, or clustering when markers have different scales | Aligning samples when distributions should be comparable | Integrating data from multiple batches before clustering |
| **Risk of removing biology** | Low (only shifts and scales) | Moderate to high (forces identical distributions; [can introduce artifacts when groups have different true distributions](https://doi.org/10.1093/biostatistics/kxv038)) | Low when batches are balanced across conditions; [high when batch and biology are confounded](https://doi.org/10.1038/s41592-019-0619-0) |
| **Reversible** | Yes | No | No |
| **Pipeline function** | `z_score_scaling()` | `quantile_norm()` | `harmony_normalize()` |
| **Implementation** | `flowCore::transformList` + `scale` | `preprocessCore::normalize.quantiles` | `harmony::HarmonyMatrix` |
| **Cross-link** | This document | This document | [Batch Correction (Harmony)](batch-correction.md) |

### When to Combine Them

These methods are not mutually exclusive. A typical multi-batch workflow might:

1. Transform raw data ([arcsinh](arcsinh-transform.md) or
   [logicle](logicle-transform.md)).
2. Apply Harmony to correct batch effects across samples.
3. Apply z-score scaling to equalize marker contributions before clustering.

However, quantile normalization and Harmony should generally not both be applied
to the same data, because both attempt to align distributions across samples.
Applying both risks over-correction and loss of real biological signal.

---

## How to Use in the Pipeline

### Z-Score Scaling

```r
# Z-score scaling
source("src/preprocessing/normalization.R")
scaled_ff <- z_score_scaling(transformed_ff, marker_indices = 1:15)
```

The function takes a flowFrame object and an integer vector of column indices
specifying which markers to scale. Markers not included in `marker_indices`
(e.g., scatter parameters, time channel) are left unchanged.

### Quantile Normalization

```r
# Quantile normalization
source("src/preprocessing/normalization.R")
qn_ff <- quantile_norm(transformed_ff, marker_indices = 1:15)
```

This function requires the `preprocessCore` Bioconductor package. It extracts
the expression matrix, applies `normalize.quantiles()` to the selected columns,
and writes the result back into the flowFrame.

### Within the Preprocessing Pipeline

Both normalization steps are integrated into the `preprocess_pipeline.R`
orchestrator and controlled by config flags:

```r
source("src/preprocessing/preprocess_pipeline.R")

params <- list(
  pp_order            = "QT",
  qc_type             = "None",
  transform_type      = "arcsinh",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD45", "CD19"),
  offset_list         = rep(0, 5),
  width_list          = rep(5, 5),
  quant_norm          = TRUE,    # Apply quantile normalization
  scaling             = TRUE,    # Apply z-score scaling
  save_pp_fcs         = FALSE,
  save_pp_csv         = FALSE
)

preprocessed <- run_preprocessing(flow_frames, params, output_dir = "output/preprocessing")
```

When both `quant_norm` and `scaling` are `TRUE`, quantile normalization is
applied first, followed by z-score scaling.

---

## Order in the Pipeline

Normalization occurs **after** transformation and **optionally after or before**
downsampling, depending on the analysis goals.

```
Raw FCS
  |
  v
Quality Control (src/preprocessing/quality_control.R)
  |
  v
Transformation — arcsinh or logicle (src/preprocessing/transformation.R)
  |
  v
[Optional] Quantile Normalization (src/preprocessing/normalization.R)
  |
  v
[Optional] Z-Score Scaling (src/preprocessing/normalization.R)
  |
  v
[Optional] Downsampling (src/preprocessing/downsampling.R)
  |
  v
[Optional] Batch Correction — Harmony (src/harmonization/batch_correction.R)
  |
  v
Clustering, Visualization, Analysis
```

### Why After Transformation

Normalization should be applied to transformed data, not raw intensities. Raw
cytometry intensities span a wide, skewed dynamic range where the mean and
standard deviation are dominated by bright events. Transformation
([arcsinh](arcsinh-transform.md) or [logicle](logicle-transform.md))
compresses this range and produces approximately symmetric distributions.
Z-score scaling and quantile normalization are
[most effective when applied to data with approximately symmetric, unimodal or bimodal distributions](https://doi.org/10.12688/f1000research.11622.3)
rather than heavily skewed raw intensities.

### Normalization and Downsampling

Quantile normalization uses the full distribution of each marker to compute
the reference quantiles. Applying it **before** downsampling uses more cells to
estimate the reference distribution, which produces a more accurate alignment.
If normalization is applied **after** downsampling, the reference distribution
is estimated from fewer cells and may be noisier, particularly for rare
populations in the tails.

Z-score scaling is less sensitive to sample size because it depends only on the
first two moments. It can be applied before or after downsampling with minimal
difference in outcome.

### Cross-Links

- **Transformation (prerequisite):** [Arcsinh Transform](arcsinh-transform.md),
  [Logicle Transform](logicle-transform.md)
- **Batch correction (alternative/complement):**
  [Batch Correction (Harmony)](batch-correction.md)
- **Downstream clustering:** [Leiden Clustering](leiden-clustering.md)
- **Downstream dimensionality reduction:** [PCA](pca.md), [UMAP](umap.md)
- **Preprocessing orchestrator:**
  [Preprocessing subproject reference](../subprojects/preprocessing.md)

---

## Cautions

### Do Not Over-Normalize

Applying multiple normalization methods in sequence (e.g., quantile
normalization followed by z-score scaling followed by Harmony) risks removing
real biological variation along with technical noise. Each normalization step
makes the data more uniform, and the cumulative effect can be
[a loss of the very differences you are trying to detect](https://doi.org/10.1093/biostatistics/kxv038).
Use the minimum normalization necessary for your analysis.

### Do Not Normalize Away Real Biological Signal

If your experimental groups have genuinely different marker expression
distributions (e.g., treated cells express more of a marker than untreated
cells), quantile normalization will force these distributions to match and
**erase the treatment effect**. Hicks and Irizarry (2015) demonstrated that
[quantile normalization can introduce systematic bias when applied to data where groups have different true underlying distributions](https://doi.org/10.1093/biostatistics/kxv038).
This is the single most important caution for quantile normalization in
cytometry: it assumes that distributional differences are technical, not
biological.

Z-score scaling is safer in this regard because it only adjusts the mean and
scale, not the distribution shape. However, if one group has a genuinely higher
mean expression of a marker, z-score scaling within that group will center it
at zero just like the control group, removing the mean difference. When
comparing groups, apply z-score scaling to the **combined data** (all groups
together), not separately per group.

### Verify Normalization Visually

Always inspect density plots or histograms before and after normalization.
The [Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
recommends visual QC at every preprocessing step. After z-score scaling, each
marker should have zero mean and unit variance. After quantile normalization,
overlaid density plots of different markers (or samples) should be identical.
If they are not, the normalization may have failed or been applied to the wrong
columns.

### Do Not Normalize Scatter Parameters

Forward scatter (FSC) and side scatter (SSC) channels reflect cell size and
granularity, not protein expression. They operate on a different physical scale
and should be excluded from both z-score scaling and quantile normalization.
Similarly, the Time channel should never be normalized. Pass only marker
indices corresponding to fluorescence or mass cytometry channels to the
normalization functions.

### Consider the Effect on Downstream Statistics

If normalized data is used for differential expression testing (e.g.,
comparing median marker expression between conditions), the normalization
changes the null distribution of the test statistic. Z-score-scaled data
produces test statistics in units of standard deviations rather than original
intensity units. Quantile-normalized data has had its distributional
differences removed, so tests for distributional differences will have
[reduced power or inflated false positives depending on the true data-generating process](https://doi.org/10.1093/biostatistics/kxv038).
For formal statistical testing, it is generally safer to normalize only for
visualization and clustering, and to use the un-normalized (but transformed)
data for statistical inference, including batch as a covariate in the model
when necessary.

---

## Implementation Details

### Source Code Location

- **Z-score scaling:** [`z_score_scaling()`](/src/preprocessing/normalization.R)
  in `src/preprocessing/normalization.R`
- **Quantile normalization:** [`quantile_norm()`](/src/preprocessing/normalization.R)
  in `src/preprocessing/normalization.R`
- **Harmony (cross-link):** [`harmony_normalize()`](/src/preprocessing/normalization.R)
  in `src/preprocessing/normalization.R`
- **Preprocessing orchestrator:**
  [`run_preprocessing()`](/src/preprocessing/preprocess_pipeline.R)
  in `src/preprocessing/preprocess_pipeline.R`

### Dependencies

| Function | Required Package | Bioconductor / CRAN |
|----------|-----------------|---------------------|
| `z_score_scaling()` | `flowCore` | Bioconductor |
| `quantile_norm()` | `flowCore`, `preprocessCore` | Bioconductor |
| `harmony_normalize()` | `harmony` | CRAN |

### Error Handling

- `z_score_scaling()` will produce `NaN` values if a marker has zero variance
  (all cells have the same value), because `scale()` divides by zero standard
  deviation. Filter constant markers before scaling.
- `quantile_norm()` calls `require(preprocessCore)` and will fail with an
  informative error if the package is not installed.
- Both functions modify the expression matrix in place within the flowFrame
  and return the modified flowFrame. The original object is not changed (R uses
  copy-on-modify semantics).

### References by Topic

All references are linked inline throughout this document. The key publications
are:

- **Quantile normalization (foundational):** [Bolstad, B. M., Irizarry, R. A., Astrand, M., & Speed, T. P. (2003)](https://doi.org/10.1093/bioinformatics/19.1.185). A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. *Bioinformatics*, 19(1), 185-193.
- **Quantile normalization cautions:** [Hicks, S. C. & Irizarry, R. A. (2015)](https://doi.org/10.1093/biostatistics/kxv038). quantro: A data-driven approach to guide the choice of an appropriate normalization method. *Biostatistics*, 16(2), 413-426.
- **CyTOF bead normalization:** [Finck, R. et al. (2013)](https://doi.org/10.1002/cyto.a.22271). Normalization of mass cytometry data with bead standards. *Cytometry A*, 83(5), 483-494.
- **CyTOF analysis workflow:** [Nowicka, M. et al. (2017)](https://doi.org/10.12688/f1000research.11622.3). CyTOF workflow: Differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*, 6, 748.
- **Harmony batch correction:** [Korsunsky, I. et al. (2019)](https://doi.org/10.1038/s41592-019-0619-0). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 16(12), 1289-1296.
- **PCA and scaling (textbook):** [Jolliffe, I. T. (2002)](https://doi.org/10.1007/b98835). *Principal Component Analysis* (2nd ed.). Springer.
- **Statistics reference (z-score):** [Dodge, Y. (2008)](https://doi.org/10.1007/978-0-387-32833-1). *The Concise Encyclopedia of Statistics*. Springer.
- **Scanpy single-cell workflow:** [Wolf, F. A., Angerer, P., & Theis, F. J. (2018)](https://doi.org/10.1186/s13059-017-1382-0). SCANPY: Large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
- **UMAP for single-cell data:** [Becht, E. et al. (2019)](https://doi.org/10.1038/nbt.4314). Dimensionality reduction for visualizing single-cell data using UMAP. *Nature Biotechnology*, 37(1), 38-44.
- **Flow cytometry normalization in FlowCAP:** [Finak, G. et al. (2014)](https://doi.org/10.1186/1471-2105-15-S10-S2). Standardizing flow cytometry immunophenotyping analysis from the Human ImmunoPhenotyping Consortium. *BMC Bioinformatics*, 15(Suppl 10), S2.
- **When quantile normalization fails:** [Bullard, J. H. et al. (2010)](https://doi.org/10.1093/biostatistics/kxj037). Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. *BMC Bioinformatics*, 11, 94.
