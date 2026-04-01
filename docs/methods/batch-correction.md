# Harmony Batch Correction

## Summary

Harmony is a fast, scalable batch correction algorithm that [iteratively adjusts cell embeddings using soft clustering and centroid-based correction](https://doi.org/10.1038/s41592-019-0619-0) to remove technical variation between acquisition batches while preserving biological signal. In this pipeline, `harmony::HarmonyMatrix()` operates directly on the expression matrix (with `do_pca = FALSE`), producing a corrected matrix of the same dimensions that can be passed to downstream clustering and visualization without further transformation. Batch correction sits between preprocessing/homogenization and clustering in Pipelines A and E, and should only be applied when acquisition batches are a documented source of unwanted variation.

---

## The Method

### Overview

Harmony was introduced by [Korsunsky et al. (2019)](https://doi.org/10.1038/s41592-019-0619-0) as an integration algorithm for single-cell data. It removes batch effects through an iterative process that alternates between soft clustering of cells and linear correction of cluster centroids, converging when the batch composition within each cluster matches the global batch proportions.

### Soft Clustering

Rather than assigning each cell to a single cluster, Harmony computes a [probabilistic assignment across multiple clusters](https://doi.org/10.1038/s41592-019-0619-0) using a modified k-means objective. Each cell receives a membership weight for every cluster, and these weights are updated at each iteration. This soft assignment avoids the hard boundary artifacts that arise when batch effects shift cells between discrete cluster boundaries. The number of clusters is controlled by the `nclust` parameter (Harmony default: automatic based on dataset size), which governs the granularity of the correction.

### Centroid-Based Correction

Within each soft cluster, Harmony computes per-batch centroids and estimates the [linear displacement needed to align them](https://doi.org/10.1038/s41592-019-0619-0). For each cluster-batch combination, a correction vector is calculated as the difference between the batch-specific centroid and the cluster-wide centroid. These correction vectors are applied to the cell embeddings, weighted by the soft cluster memberships. Because corrections are computed locally (per cluster) rather than globally, [different cell populations can receive different correction magnitudes](https://doi.org/10.1038/s41592-019-0619-0), which is critical when batch effects are not uniform across cell types.

### Iterative Convergence

The soft clustering and centroid correction steps alternate until the [objective function converges](https://doi.org/10.1038/s41592-019-0619-0) or a maximum number of iterations is reached (default: 10). At each iteration, the batch entropy within clusters is evaluated: a well-mixed cluster has high entropy (all batches equally represented), while a batch-dominated cluster has low entropy. The algorithm converges when further iterations no longer meaningfully increase the [within-cluster batch mixing](https://doi.org/10.1038/s41592-019-0619-0). The penalty parameter `theta` (not to be confused with the negative binomial dispersion parameter in [NB-GLMM](negative-binomial-glmm.md)) controls how aggressively Harmony penalizes batch-dominated clusters: higher values produce stronger correction, while lower values preserve more of the original structure.

### Direct Operation on Expression Data

The original Harmony publication applied the algorithm to PCA embeddings of scRNA-seq data. In this cytometry pipeline, Harmony operates directly on the [transformed expression matrix](https://doi.org/10.1038/s41592-019-0619-0) with `do_pca = FALSE`. This is appropriate because flow cytometry panels typically contain 10-40 markers -- a dimensionality low enough that [PCA reduction is unnecessary and could discard informative variation](https://doi.org/10.12688/f1000research.11622.3). The corrected matrix retains the same markers and the same number of cells as the input, with only the expression values adjusted.

---

## Why Harmony for Cytometry

### Speed and Scalability

Harmony was designed for [datasets of hundreds of thousands to millions of cells](https://doi.org/10.1038/s41592-019-0619-0) and completes batch correction in seconds to minutes on typical cytometry datasets. This makes it practical for iterative analysis workflows where correction parameters may need to be tuned. In benchmarks by [Korsunsky et al. (2019)](https://doi.org/10.1038/s41592-019-0619-0), Harmony was the fastest among the methods tested, scaling linearly with cell count.

### Works on Expression Matrices Directly

Unlike methods that require a shared gene set, a reference dataset, or a count matrix with specific distributional assumptions, Harmony accepts any [numeric matrix of cell embeddings](https://doi.org/10.1038/s41592-019-0619-0). This flexibility makes it straightforward to apply to flow cytometry expression data after arcsinh or logicle transformation.

### Preserves Biological Structure

The soft clustering mechanism allows Harmony to [correct batch effects without merging biologically distinct populations](https://doi.org/10.1038/s41592-019-0619-0). In the original benchmarks, Harmony maintained separation between cell types that genuinely differ while removing unwanted technical variation, achieving [high scores on both batch mixing and biological conservation metrics](https://doi.org/10.1038/s41592-019-0619-0).

### Comparison with Other Methods

| Method | Mechanism | Strengths | Limitations | Key Citation |
|--------|-----------|-----------|-------------|--------------|
| **Harmony** | Iterative soft clustering + centroid correction | [Fast, scalable, works on any embedding](https://doi.org/10.1038/s41592-019-0619-0) | Assumes batch effects are approximately linear within clusters | [Korsunsky et al. 2019](https://doi.org/10.1038/s41592-019-0619-0) |
| **ComBat** | Empirical Bayes linear model for location and scale adjustment | [Well-established for microarray/bulk data](https://doi.org/10.1093/biostatistics/kxj037), simple to apply | Assumes a [single global batch effect per feature](https://doi.org/10.1093/biostatistics/kxj037); does not account for cell-type-specific effects | [Johnson et al. 2007](https://doi.org/10.1093/biostatistics/kxj037) |
| **MNN (Mutual Nearest Neighbors)** | Identifies mutual nearest neighbor pairs across batches and computes correction vectors | [Robust to differences in population composition between batches](https://doi.org/10.1038/nbt.4091) | Computationally expensive on large datasets; [sensitive to parameter choices (k)](https://doi.org/10.1038/nbt.4091) | [Haghverdi et al. 2018](https://doi.org/10.1038/nbt.4091) |
| **Scanorama** | Extends MNN with panoramic stitching across multiple batches | [Handles many batches efficiently](https://doi.org/10.1038/s41587-019-0113-3); embeds all batches into a common space | Designed for scRNA-seq; requires adaptation for cytometry marker spaces | [Hie et al. 2019](https://doi.org/10.1038/s41587-019-0113-3) |
| **CytoNorm** | Quantile normalization using bead-standardized control samples | [Purpose-built for cytometry](https://doi.org/10.1002/cyto.a.23904); uses internal controls for ground truth alignment | Requires technical replicates or reference samples across batches; [not applicable without shared controls](https://doi.org/10.1002/cyto.a.23904) | [Van Gassen et al. 2020](https://doi.org/10.1002/cyto.a.23904) |

ComBat applies a [single location-and-scale adjustment per feature per batch](https://doi.org/10.1093/biostatistics/kxj037), which assumes the batch effect is the same magnitude for all cell types. This assumption is often violated in cytometry, where [different cell populations can exhibit different degrees of batch-induced fluorescence shift](https://doi.org/10.12688/f1000research.11622.3). MNN and Harmony both allow for cell-type-specific correction, but Harmony's soft clustering is [substantially faster than MNN's nearest-neighbor search](https://doi.org/10.1038/s41592-019-0619-0). CytoNorm is the gold standard when bead-normalized reference samples are available, but many experimental designs [lack the required technical controls](https://doi.org/10.1002/cyto.a.23904). In those situations, Harmony provides a practical alternative that requires only batch labels.

---

## When Batch Correction Is Needed

### Diagnosing Batch Effects

Batch effects manifest as systematic shifts in marker expression that correlate with technical variables (instrument, acquisition date, reagent lot, operator) rather than biological conditions. The primary diagnostic tools are:

- **PCA**: Plot the first two principal components and color cells by batch. If batches form [visually separable clusters in PC space](https://doi.org/10.12688/f1000research.11622.3), a batch effect is present. See [PCA](pca.md) for implementation details.
- **UMAP**: A [UMAP embedding colored by batch](https://doi.org/10.12688/f1000research.11622.3) provides a nonlinear view of the same phenomenon. Batch-driven separation is visible as spatially segregated groups that do not correspond to known biology. See [UMAP](umap.md) for implementation details.
- **Density plots per marker**: Overlay per-batch density plots for each marker. [Systematic shifts in the mode or spread of distributions](https://doi.org/10.1002/cyto.a.23904) across batches indicate technical variation requiring correction.

### Common Sources of Batch Effects in Flow Cytometry

| Source | Mechanism | Typical Magnitude |
|--------|-----------|-------------------|
| **Different acquisition dates** | Laser alignment drift, [photomultiplier tube (PMT) voltage changes over time](https://doi.org/10.1002/cyto.a.23904) | Moderate to large |
| **Different instruments** | Optical path differences, [detector sensitivity variation](https://doi.org/10.1002/cyto.a.23904) | Large |
| **Different antibody lots** | Variation in [fluorochrome-to-protein ratio and binding affinity](https://doi.org/10.12688/f1000research.11622.3) | Variable |
| **Different operators** | Differences in [sample preparation, staining time, cell concentration](https://doi.org/10.12688/f1000research.11622.3) | Moderate |
| **Longitudinal studies** | All of the above, [accumulated over weeks to months](https://doi.org/10.1002/cyto.a.23904) | Often large |

---

## When Your Data Calls For It

Batch correction with Harmony is appropriate when **all** of the following conditions hold:

### 1. Multiple Acquisition Batches

Your samples were acquired across two or more batches -- different days, instruments, or experimental runs. If all samples were acquired in a single session on a single instrument, [there is no batch structure to correct](https://doi.org/10.1038/s41592-019-0619-0).

### 2. Visible Batch Effects in QC Plots

PCA or UMAP plots colored by batch show [spatial separation that does not correspond to known biological groups](https://doi.org/10.12688/f1000research.11622.3). Per-marker density plots reveal systematic shifts between batches. If batches are already well-mixed in the embedding space, correction is unnecessary and could introduce artifacts.

### 3. Conditions Are Balanced (or Nearly So) Across Batches

Conditions of interest must be distributed across batches so that [batch and condition are not confounded](https://doi.org/10.1038/s41592-019-0619-0). For example, if half of the treatment samples and half of the control samples appear in each batch, Harmony can remove batch variation while preserving the treatment effect. If all treatment samples are in batch 1 and all controls are in batch 2, batch correction will remove the treatment effect along with the batch effect (see the confounding warning below).

### 4. Correction Is Applied Before Clustering and Visualization

Batch correction should be applied [before any step that uses inter-cell distances](https://doi.org/10.12688/f1000research.11622.3): Leiden clustering (which builds a k-nearest-neighbor graph), UMAP embedding, and PCA for visualization. Correcting after clustering risks splitting true biological clusters into batch-specific subclusters.

---

## When to Choose Something Else

### Single Batch

If all samples were acquired in one run on one instrument, there is no batch structure. Applying Harmony when [no batch effect exists can introduce spurious corrections](https://doi.org/10.1038/s41592-019-0619-0) that reduce biological signal. Skip correction entirely.

### Batch Confounded with Biology

When batch perfectly correlates with the experimental condition (all cases in batch 1, all controls in batch 2), no algorithm can separate batch from biology. Harmony will remove [both the batch effect and the biological effect](https://doi.org/10.1038/s41592-019-0619-0). In this scenario, do not apply batch correction. Instead, include batch as a covariate in downstream statistical models such as the [NB-GLMM](negative-binomial-glmm.md):

```r
count ~ condition + batch + offset(log_total) + (1 | sample_id)
```

This allows the model to [estimate and adjust for batch effects while still testing the condition effect](https://doi.org/10.1038/nmeth.4295).

### Biological Covariates Driving Separation

If PCA/UMAP separation is driven by a biological variable (e.g., tissue type, disease subtype) rather than a technical batch, Harmony is not the right tool. Removing this variation would [eliminate biologically meaningful differences](https://doi.org/10.1038/s41592-019-0619-0). Confirm the source of variation by coloring embeddings by both batch labels and biological covariates before deciding to correct.

### Downstream Model Already Accounts for Batch

If your statistical model includes batch as a fixed effect (as in the formula above) or as a random effect, applying Harmony before modeling would [double-correct for the same source of variation](https://doi.org/10.12688/f1000research.11622.3), potentially removing real signal and deflating effect sizes. Choose one approach: either correct the data matrix with Harmony and then run a batch-free model, or skip correction and include batch in the model.

---

## Confounding Warning

If batch and condition are perfectly confounded (every sample from condition A was acquired in batch 1 and every sample from condition B was acquired in batch 2), **no batch correction algorithm can recover the biological signal**. This is a fundamental limitation of observational study design, not a limitation of any particular algorithm. [Korsunsky et al. (2019)](https://doi.org/10.1038/s41592-019-0619-0) explicitly noted that Harmony assumes the user wants to remove the variation associated with the batch variable, and it will do so regardless of whether that variation is technical or biological.

[Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) formalized this concern for cytometry differential abundance analysis: when batch and condition are aliased, the condition effect is absorbed into the batch effect estimate, producing [unbiased batch estimates but biased (typically null) condition estimates](https://doi.org/10.1038/nmeth.4295). The [Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3) recommends that experiments be designed with conditions balanced across batches whenever possible, and that when confounding is unavoidable, the batch variable should be reported transparently and included as a covariate rather than corrected away.

In practice, partial confounding (some conditions overrepresented in certain batches) can still be handled by Harmony, but [the correction will attenuate biological differences that are correlated with batch](https://doi.org/10.1038/s41592-019-0619-0). The degree of attenuation depends on the extent of confounding and the value of the `theta` parameter.

---

## How to Use in the Pipeline

### Batch Correction via `run_batch_correction()` (Pipeline A/E)

The primary interface accepts a named list of flowFrame objects, sample names, and batch labels. It concatenates expression matrices, runs Harmony, and returns the corrected matrix with metadata:

```r
source("src/harmonization/batch_correction.R")

corrected <- run_batch_correction(
  flow_frames         = homogenized_frames,
  sample_names        = sample_names,
  batch_labels        = batch_vector,
  markers_of_interest = params$markers_of_interest,
  do_pca              = FALSE       # default; operate on expression directly
)

# corrected$corrected_matrix  — numeric matrix (cells x markers)
# corrected$sample_labels     — per-cell sample assignment
# corrected$batch_labels      — per-cell batch assignment
# corrected$markers           — marker column names
```

This function lives in `src/harmonization/batch_correction.R` and internally calls `harmony::HarmonyMatrix()` with `vars_use = "batch"`.

### Via the Normalization Wrapper

For simpler workflows where you already have a combined expression matrix and batch vector, the lightweight wrapper in `src/preprocessing/normalization.R` provides a direct matrix-in, matrix-out interface:

```r
source("src/preprocessing/normalization.R")

corrected_matrix <- harmony_normalize(
  expr_matrix  = combined_matrix,
  batch_vector = batch_vector
)
```

This calls `harmony::HarmonyMatrix()` with `do_pca = FALSE` and returns the corrected matrix directly.

### Configuration

Batch labels are defined in the config file via `batch_names` and the derived `sample_batch_dict`:

```r
params <- list(
  sample_names = c("Sample1", "Sample2", "Sample3"),
  batch_names  = c("1", "1", "2"),
  # ...
)
params$sample_batch_dict <- stats::setNames(params$batch_names, params$sample_names)
```

### Placement in the Pipeline

Batch correction should be inserted **after** panel homogenization and **before** downsampling and clustering. In a Pipeline A workflow, the sequence is:

```
Load FCS -> QC -> Transform -> Homogenize panels -> [Batch correct] -> Downsample -> Cluster -> Analyze -> Visualize
```

The batch correction step is optional. Enable it when the diagnostics described above indicate batch effects are present.

---

## Effect on Downstream Analysis

### Apply Before Clustering and Visualization

Harmony correction should precede [Leiden clustering](leiden-clustering.md) and [UMAP](umap.md) because both methods rely on inter-cell distances. Uncorrected batch effects create artificial distances between cells from different batches, leading to [batch-specific subclusters that do not reflect true biology](https://doi.org/10.12688/f1000research.11622.3). The corrected matrix produces a nearest-neighbor graph that reflects biological relationships rather than technical ones, improving both cluster quality and UMAP layout.

### Do Not Double-Correct

If you apply Harmony before clustering, **do not** also include batch as a covariate in the [NB-GLMM](negative-binomial-glmm.md) or [mixed-effects logistic regression](mixed-effects-logistic.md). Doing so [adjusts for the same variation twice](https://doi.org/10.12688/f1000research.11622.3), which can:

- Remove real biological signal that is weakly correlated with batch
- Reduce statistical power by absorbing additional degrees of freedom
- Produce overly conservative p-values and attenuated fold changes

Choose one strategy:

| Strategy | When to Use |
|----------|-------------|
| **Harmony correction + batch-free model** | Batch effects are visible in QC, conditions are balanced across batches |
| **No correction + batch as covariate** | Partial confounding, or you want the model to estimate and adjust for batch effects directly |

### Marker Interpretation

After Harmony correction, the expression values in the corrected matrix are no longer on the original transformed scale. They are [adjusted values that have been shifted to remove batch-specific offsets](https://doi.org/10.1038/s41592-019-0619-0). This means:

- Marker expression heatmaps from corrected data reflect batch-normalized relative expression
- Absolute fluorescence intensity thresholds established on raw data may not apply to corrected values
- For gating or thresholding, use the uncorrected data; for clustering and embedding, use the corrected data

---

## Validating Correction

Batch correction should be validated by comparing the data before and after correction using complementary metrics. The goal is to confirm that (a) batch effects have been removed and (b) biological structure has been preserved.

### Visual Validation

**PCA/UMAP before vs. after correction:** Generate embeddings from both the uncorrected and corrected matrices, colored by batch. Before correction, batches should show [visible spatial separation](https://doi.org/10.12688/f1000research.11622.3). After correction, batches should be well-mixed within each biological cluster. Simultaneously color by condition to verify that [condition-driven separation is preserved](https://doi.org/10.1038/s41592-019-0619-0).

**Per-marker density overlays:** Plot per-batch density curves for each marker. After correction, the [distributions should overlap substantially](https://doi.org/10.1002/cyto.a.23904), while retaining shifts that correspond to real biological differences between conditions.

### Quantitative Metrics

| Metric | What It Measures | Interpretation | Citation |
|--------|-----------------|----------------|----------|
| **kBET (k-nearest-neighbor Batch Effect Test)** | Whether the batch label distribution in local neighborhoods matches the global distribution | [Rejection rate near zero indicates good mixing](https://doi.org/10.1038/s41592-018-0254-1); high rejection indicates residual batch effects | [Buttner et al. 2019](https://doi.org/10.1038/s41592-018-0254-1) |
| **Silhouette score by batch** | How well cells cluster by batch | [Silhouette near zero or negative after correction](https://doi.org/10.1038/s41592-019-0619-0) means batches are no longer separable; high positive values indicate residual batch structure | [Korsunsky et al. 2019](https://doi.org/10.1038/s41592-019-0619-0) |
| **Silhouette score by cell type** | How well cells cluster by biological label | Should [remain high or increase after correction](https://doi.org/10.1038/s41592-019-0619-0), confirming that biological structure is preserved | [Korsunsky et al. 2019](https://doi.org/10.1038/s41592-019-0619-0) |
| **LISI (Local Inverse Simpson Index)** | Diversity of batch or cell-type labels in local neighborhoods | [Higher batch LISI = better mixing; higher cell-type LISI would indicate loss of biological structure](https://doi.org/10.1038/s41592-019-0619-0) | [Korsunsky et al. 2019](https://doi.org/10.1038/s41592-019-0619-0) |
| **Adjusted Rand Index (ARI)** | Agreement between clusters found before and after correction | [High ARI indicates that correction preserved cluster identity](https://doi.org/10.12688/f1000research.11622.3); low ARI suggests overcorrection or undercorrection | [Nowicka et al. 2017](https://doi.org/10.12688/f1000research.11622.3) |

### Validation Checklist

1. **Batch mixing improved.** PCA/UMAP colored by batch shows interspersed points after correction.
2. **Biological separation preserved.** PCA/UMAP colored by condition or cell type retains meaningful groupings.
3. **Per-marker distributions aligned.** Density overlays show batch distributions converging while condition-specific shifts persist.
4. **kBET acceptance rate increased.** Quantitative confirmation that [local neighborhoods reflect global batch composition](https://doi.org/10.1038/s41592-018-0254-1).
5. **No new artifacts introduced.** Check for clusters that appear only in the corrected data and have no clear biological identity.

---

## Comparison with Other Methods

### ComBat

ComBat uses an [empirical Bayes framework to estimate and remove location (mean) and scale (variance) batch effects per feature](https://doi.org/10.1093/biostatistics/kxj037). Originally developed for microarray gene expression data by [Johnson et al. (2007)](https://doi.org/10.1093/biostatistics/kxj037), it applies a single global correction per marker per batch. This works well when batch effects are uniform across all cell types but [fails when different populations experience different batch shifts](https://doi.org/10.1038/s41592-019-0619-0) -- a common scenario in cytometry where fluorescence spillover and compensation errors affect markers differently depending on co-expression patterns. ComBat also assumes approximately [normally distributed expression values within each batch](https://doi.org/10.1093/biostatistics/kxj037), which may not hold for cytometry data with multimodal distributions.

### MNN (Mutual Nearest Neighbors)

MNN, introduced by [Haghverdi et al. (2018)](https://doi.org/10.1038/nbt.4091), identifies pairs of cells from different batches that are each other's nearest neighbors and uses these pairs to estimate batch correction vectors. The key advantage is that MNN [naturally adapts to different batch effects in different parts of the expression space](https://doi.org/10.1038/nbt.4091), similar to Harmony. However, MNN is [computationally expensive for large datasets](https://doi.org/10.1038/s41592-019-0619-0) because it requires computing pairwise distances across batches. For cytometry datasets with hundreds of thousands of cells, this can be prohibitively slow. MNN is also [sensitive to the choice of k](https://doi.org/10.1038/nbt.4091) (number of nearest neighbors) and can produce inconsistent results when batches have very different population compositions.

### CytoNorm

CytoNorm, developed by [Van Gassen et al. (2020)](https://doi.org/10.1002/cyto.a.23904), is a [quantile normalization method specifically designed for flow and mass cytometry](https://doi.org/10.1002/cyto.a.23904). It uses reference samples (bead standards or aliquots of the same biological sample) that are measured across all batches to learn the mapping between batch-specific distributions. This approach has the [strongest theoretical justification when reference samples are available](https://doi.org/10.1002/cyto.a.23904), because the reference provides ground truth for what the same biology should look like across batches. However, many experiments [lack the required technical replicates or bead controls](https://doi.org/10.1002/cyto.a.23904), making CytoNorm inapplicable. When controls are available, CytoNorm is the preferred choice; when they are not, Harmony offers a robust alternative that requires only batch labels.

### Practical Recommendation

| Scenario | Recommended Method |
|----------|--------------------|
| Reference samples available across all batches | [CytoNorm](https://doi.org/10.1002/cyto.a.23904) |
| No reference samples, large dataset (>50k cells), multiple batches | [Harmony](https://doi.org/10.1038/s41592-019-0619-0) |
| No reference samples, small dataset, few batches | [MNN](https://doi.org/10.1038/nbt.4091) or [Harmony](https://doi.org/10.1038/s41592-019-0619-0) |
| Bulk-level normalization before single-cell analysis | [ComBat](https://doi.org/10.1093/biostatistics/kxj037) |
| Batch confounded with condition | No correction; include batch as covariate in [NB-GLMM](negative-binomial-glmm.md) |

---

## Key Parameters

The `harmony::HarmonyMatrix()` function accepts several parameters that affect the correction. The most important ones are:

| Parameter | Default | Description | Guidance |
|-----------|---------|-------------|----------|
| `theta` | 2 | [Diversity penalty for batch mixing](https://doi.org/10.1038/s41592-019-0619-0). Higher values enforce stronger correction. | Start with the default. Increase if batches remain poorly mixed; decrease if biological structure is lost. |
| `nclust` | automatic | Number of soft clusters used for correction. | Harmony selects this automatically. Override only if you have strong prior knowledge about the number of populations. |
| `max.iter.harmony` | 10 | [Maximum number of Harmony iterations](https://doi.org/10.1038/s41592-019-0619-0). | Increase if convergence warnings appear; 20 is usually sufficient. |
| `do_pca` | `FALSE` (in this pipeline) | Whether to run PCA before Harmony. | Keep `FALSE` for cytometry expression matrices. Set `TRUE` only for high-dimensional inputs (>100 features). |
| `sigma` | 0.1 | [Width of the soft clustering kernel](https://doi.org/10.1038/s41592-019-0619-0). | Smaller values produce sharper cluster assignments; larger values produce smoother (more global) corrections. |

These can be passed through the `...` argument in both `run_batch_correction()` and `harmony_normalize()`:

```r
corrected <- run_batch_correction(
  flow_frames         = homogenized_frames,
  sample_names        = sample_names,
  batch_labels        = batch_vector,
  markers_of_interest = params$markers_of_interest,
  theta               = 3,              # stronger correction
  max.iter.harmony    = 20              # more iterations
)
```

---

## Implementation Details

### Source Code Locations

- **Primary function:** [`run_batch_correction()`](/src/harmonization/batch_correction.R) in `src/harmonization/batch_correction.R`
- **Lightweight wrapper:** [`harmony_normalize()`](/src/preprocessing/normalization.R) in `src/preprocessing/normalization.R`
- **Configuration template:** `config/Parameters_template.R` (defines `batch_names` and `sample_batch_dict`)

### How `run_batch_correction()` Works

1. **Concatenates** expression matrices from all flowFrame objects into a single combined matrix.
2. **Subsets** to `markers_of_interest` if provided, using `explore_data()` from `src/R_shared/flowframe_utils.R` to map marker names to column indices.
3. **Builds** a metadata data frame with `batch` and `sample` columns, with one row per cell.
4. **Calls** `harmony::HarmonyMatrix()` with `vars_use = "batch"` and `do_pca = FALSE`.
5. **Returns** a list containing the corrected matrix, per-cell sample labels, per-cell batch labels, and marker names.

### How `harmony_normalize()` Works

1. **Wraps** a pre-built expression matrix and batch vector into a metadata data frame.
2. **Calls** `harmony::HarmonyMatrix()` with `vars_use = "batch"` and `do_pca = FALSE`.
3. **Returns** the corrected matrix directly (no additional metadata).

### Dependencies

- `harmony` (installed via `BiocManager::install("harmony")` in `install_extra_packages.R`)
- `flowCore` (for `exprs()` extraction in `run_batch_correction()`)
- `src/R_shared/io_utils.R` and `src/R_shared/flowframe_utils.R` (for `explore_data()`)

---

## Further Reading

The following resources are linked inline at the point of relevance throughout this document:

- [Korsunsky, I., Millard, N., Fan, J. et al. (2019)](https://doi.org/10.1038/s41592-019-0619-0) --- Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 16, 1289-1296.
- [Johnson, W. E., Li, C. & Rabinovic, A. (2007)](https://doi.org/10.1093/biostatistics/kxj037) --- Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118-127.
- [Haghverdi, L., Lun, A. T. L., Morgan, M. D. & Marioni, J. C. (2018)](https://doi.org/10.1038/nbt.4091) --- Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. *Nature Biotechnology*, 36, 421-427.
- [Hie, B., Bryson, B. & Berger, B. (2019)](https://doi.org/10.1038/s41587-019-0113-3) --- Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. *Nature Biotechnology*, 37, 685-691.
- [Van Gassen, S., Gaudilliere, B., Angst, M. S. et al. (2020)](https://doi.org/10.1002/cyto.a.23904) --- CytoNorm: A normalization algorithm for cytometry data. *Cytometry Part A*, 97(3), 268-278.
- [Buttner, M., Miao, Z., Wolf, F. A., Teichmann, S. A. & Theis, F. J. (2019)](https://doi.org/10.1038/s41592-018-0254-1) --- A test metric for assessing single-cell RNA-seq batch correction. *Nature Methods*, 16, 43-49.
- [Nowicka, M., Krieg, C., Crowell, H. L. et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) --- CyTOF workflow: Differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*, 6, 748.
- [Lun, A. T. L., Richard, A. C. & Marioni, J. C. (2017)](https://doi.org/10.1038/nmeth.4295) --- Testing for differential abundance in mass cytometry data. *Nature Methods*, 14(7), 707-709.
