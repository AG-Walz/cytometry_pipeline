# Arcsinh Transform for Cytometry Data

## Summary

The inverse hyperbolic sine (arcsinh) transform maps raw ion counts or
fluorescence intensities onto a scale that behaves linearly near zero and
logarithmically at large values, preserving zero and negative measurements
while compressing dynamic range. It is
[the standard transformation for CyTOF (mass cytometry) data](https://doi.org/10.1126/science.1198704)
and is controlled by a single parameter -- the cofactor -- which sets the
transition point between linear and logarithmic behavior. In this pipeline,
`transform_flowframe()` in `src/preprocessing/transformation.R` applies the
arcsinh transform to flowFrame objects with per-marker cofactor optimization.

---

## The Method

### Mathematical Definition

The arcsinh transform maps a raw measurement *x* to a transformed value *y*
using:

$$
y = \operatorname{arcsinh}\!\left(\frac{x}{c}\right) = \ln\!\left(\frac{x}{c} + \sqrt{\frac{x^2}{c^2} + 1}\right)
$$

where *c* is the cofactor. This function is defined for all real numbers,
including zero and negative values, which distinguishes it from the logarithmic
transform.

### Behavior Near Zero

For values small relative to the cofactor (|*x*| << *c*), the Taylor expansion
gives:

$$
\operatorname{arcsinh}\!\left(\frac{x}{c}\right) \approx \frac{x}{c}
$$

In this regime the transform is approximately linear, preserving the
proportional spacing of measurements near the noise floor. This property is
critical for cytometry because
[ion counting noise in CyTOF produces values that fluctuate around zero](https://doi.org/10.1126/science.1198704),
and a linear response near zero prevents artificial compression or inflation
of these small signals.

### Behavior for Large Values

For values large relative to the cofactor (*x* >> *c*):

$$
\operatorname{arcsinh}\!\left(\frac{x}{c}\right) \approx \ln\!\left(\frac{2x}{c}\right) = \ln(2) + \ln(x) - \ln(c)
$$

The transform behaves logarithmically, compressing the wide dynamic range of
highly-expressed markers into a manageable scale. This is essential because
cytometry measurements can span
[four to five orders of magnitude](https://doi.org/10.1002/cyto.a.22271).

### Comparison with Log and Logicle

Three transforms are commonly used in cytometry, each suited to different
data characteristics:

**Logarithmic transform.** The classic `log10(x)` is undefined for *x* <= 0
and produces extreme negative values as *x* approaches zero. For fluorescence
data with compensation artifacts that push values negative, the log transform
[creates visual artifacts and distorts population structure near zero](https://doi.org/10.1002/cyto.a.20258).
Parks, Roederer, and Moore (2006) demonstrated that log-transformed data with
negative values produces "data pileup" at the axis boundary that obscures
dim populations.

**Logicle (bi-exponential) transform.** The logicle transform was
[developed specifically for compensated fluorescence data](https://doi.org/10.1002/cyto.a.20258)
to handle the negative values produced by spectral spillover subtraction. It
uses a parametric function that is linear around zero and logarithmic at
both positive and negative extremes, with four parameters (T, W, M, A)
controlling its shape. While more flexible than arcsinh, the additional
parameters make it harder to standardize across experiments, and it was
[designed for the specific statistical properties of compensated fluorescence](https://doi.org/10.1002/cyto.a.20258)
rather than ion counting data.

**Arcsinh transform.** The arcsinh has a single tuning parameter (the cofactor)
and handles zero and negative values naturally. Bendall et al. (2011)
[adopted arcsinh as the default transform for mass cytometry](https://doi.org/10.1126/science.1198704)
because it matches the noise characteristics of CyTOF data: symmetric,
approximately Gaussian near zero, with a long positive tail for highly-expressed
markers. The
[Nowicka et al. (2017) CyTOF analysis workflow](https://doi.org/10.12688/f1000research.11622.3)
subsequently codified this as the standard preprocessing step for CyTOF
analysis pipelines.

---

## The Cofactor

### What the Cofactor Controls

The cofactor *c* sets the crossover point between the linear and logarithmic
regimes of the arcsinh function. Concretely, values much smaller than *c* are
treated approximately linearly, while values much larger than *c* are
compressed logarithmically. A smaller cofactor shifts the transition to lower
intensities, applying logarithmic compression earlier; a larger cofactor
preserves linear spacing over a wider range of values.

### Standard Values

**CyTOF data: cofactor = 5.** Bendall et al. (2011) used a cofactor of 5 in
their [foundational CyTOF study](https://doi.org/10.1126/science.1198704),
and this value has become the field standard. It is appropriate because CyTOF
measures individual ion events with counts typically ranging from 0 to several
thousand; a cofactor of 5 places the linear-to-logarithmic transition at a
point that
[separates noise from true signal](https://doi.org/10.1002/cyto.a.22271)
in this range. The
[Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
uses cofactor 5 as the default in all transformation steps.

**Fluorescence data: cofactor = 150.** Conventional flow cytometry produces
fluorescence intensity values on a much larger scale than CyTOF ion counts.
Finck et al. (2013) demonstrated that when arcsinh is applied to fluorescence
data (for example, in mixed panels or comparison studies), a cofactor of
approximately 150 provides
[appropriate scaling for fluorescence intensity ranges](https://doi.org/10.1002/cyto.a.22271).
This larger cofactor extends the linear regime to encompass the wider noise
floor of photon-based detection.

### How `find_optimal_b_cofactor()` Works

Rather than relying on a fixed cofactor, this pipeline provides a data-driven
optimization function that selects the cofactor maximizing bimodality of the
transformed distribution. The rationale is that immune cell markers are
typically expressed in a bimodal pattern (negative and positive populations),
and the optimal transform should
[maximally separate these two populations](https://doi.org/10.1002/cyto.a.22271).

The implementation in `src/preprocessing/transformation.R`:

1. **Defines a bimodality objective.** For a candidate cofactor *b*, the
   function computes `asinh(x * b)` for all values in the channel, then
   calculates the
   [bimodality coefficient](https://doi.org/10.1080/00031305.2013.847865)
   (BC) of the transformed distribution using the `mousetrap` package. The BC
   ranges from 0 (unimodal) to 1 (strongly bimodal), with values above 0.555
   generally indicating bimodality.

2. **Minimizes the distance to a target BC.** The function uses R's
   `optimize()` to find the *b* value in the interval [0.001, 10000] that
   minimizes `|BC(transformed) - target|`, where the default target is 0.9.

3. **Returns the optimal cofactor.** The *b* value returned is the reciprocal
   of the conventional cofactor (*b* = 1/*c*), matching the parameterization
   used by `flowCore::arcsinhTransform(b = ...)`.

```r
# Source code: src/preprocessing/transformation.R
find_optimal_b_cofactor <- function(data, target_bimodality_coefficient = 0.9) {
  calculate_bc <- function(b_cofactor, data, target) {
    transformed <- asinh(data * b_cofactor)
    bc <- mousetrap::bimodality_coefficient(transformed)
    abs(bc - target)
  }
  result <- optimize(
    f = calculate_bc,
    interval = c(0.001, 10000),
    data = data,
    target = target_bimodality_coefficient,
    maximum = FALSE
  )
  result$minimum
}
```

### Cofactor Sensitivity

The choice of cofactor has measurable effects on downstream analysis.
Weber and Robinson (2016) showed that
[clustering performance in CyTOF workflows is sensitive to transformation parameters](https://doi.org/10.1002/cyto.a.22764),
with suboptimal cofactors compressing or stretching the data in ways that
merge distinct populations or split homogeneous ones. The
[Nowicka et al. (2017) workflow](https://doi.org/10.12688/f1000research.11622.3)
recommends visually inspecting density plots of transformed data to verify
that the cofactor produces well-separated, biologically plausible distributions.

---

## Why Arcsinh for CyTOF

Mass cytometry (CyTOF) data has specific characteristics that make the arcsinh
transform the natural choice.

### Zero and Negative Values Are Meaningful

CyTOF measures individual metal ions via time-of-flight mass spectrometry. For
markers that a cell does not express, the detector registers zero counts.
Background noise, spectral overlap correction, and bead normalization can
produce small negative values. Unlike the log transform, the arcsinh is
[defined for all real numbers](https://doi.org/10.1126/science.1198704)
and maps zero to zero, preserving the distinction between absent expression
and low expression without requiring arbitrary floor values.

### Ion Counting Noise Is Approximately Poisson

The noise in CyTOF measurements arises from stochastic ion detection events
that follow
[approximately Poisson statistics](https://doi.org/10.1002/cyto.a.22271).
For Poisson-distributed data, the variance equals the mean, so low-count
measurements have small absolute noise while high-count measurements have
large absolute noise. The arcsinh transform acts as a
[variance-stabilizing transformation](https://doi.org/10.12688/f1000research.11622.3)
in this regime: its linear behavior near zero avoids amplifying noise in
low-count channels, while its logarithmic compression at high counts
normalizes the variance of abundant signals.

### No Compensation Is Needed

Unlike fluorescence cytometry, CyTOF uses metal isotope labels with minimal
spectral overlap, so
[compensation (spillover subtraction) is either unnecessary or minimal](https://doi.org/10.1126/science.1198704).
This means CyTOF data rarely contains the large, systematic negative values
that motivated the development of the logicle transform for compensated
fluorescence data. The arcsinh is simpler than logicle and sufficient for the
CyTOF noise profile.

### Established Community Standard

Since its use in
[the first large-scale CyTOF study](https://doi.org/10.1126/science.1198704)
by Bendall et al. (2011), arcsinh with cofactor 5 has been adopted by
essentially all major CyTOF analysis pipelines, including
[CATALYST](https://doi.org/10.12688/f1000research.11622.3),
[diffcyt](https://doi.org/10.1038/s42003-019-0415-5),
and the [Bioconductor CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
by Nowicka et al. (2017). Using the same transform facilitates comparison
across studies and ensures compatibility with published analysis protocols.

---

## When Your Data Calls For It

The arcsinh transform is appropriate when your data matches one or more of
the following scenarios.

### CyTOF (Mass Cytometry) Data

This is the canonical use case. If your input consists of FCS files from a
CyTOF instrument, arcsinh with cofactor 5 is
[the standard preprocessing step](https://doi.org/10.1126/science.1198704).
Use it unless you have a specific reason to do otherwise.

### Data Containing Zero or Negative Values

Any expression matrix where channels include zero or negative values benefits
from arcsinh over log. This includes CyTOF data after bead normalization,
fluorescence data after compensation (though logicle may be preferred for
fluorescence -- see below), and any dataset where background subtraction
produces negative values. The arcsinh
[handles these values without requiring ad hoc adjustments](https://doi.org/10.12688/f1000research.11622.3)
such as adding a pseudocount or clipping to zero.

### Need for Linear Behavior Near Zero Plus Logarithmic Compression

When you need to preserve the spacing of low-intensity measurements (to avoid
compressing dim populations) while simultaneously compressing high-intensity
measurements (to visualize the full dynamic range), the arcsinh provides
[exactly this dual behavior](https://doi.org/10.1126/science.1198704)
with a single parameter controlling the transition point.

### Preparing Data for Automated Analysis

When data will be fed into clustering algorithms, dimensionality reduction, or
machine learning models, arcsinh transformation
[normalizes the scale across markers](https://doi.org/10.12688/f1000research.11622.3)
so that all channels contribute proportionally to distance-based computations.
Without transformation, markers with high maximum intensities would dominate
algorithms like Leiden clustering and UMAP that depend on inter-cell distances.

---

## When to Choose Something Else

The arcsinh transform is not universally optimal. The following situations call
for different transforms.

| Situation | Why Not Arcsinh | Better Choice |
|-----------|----------------|---------------|
| Compensated fluorescence data with large spillover corrections | Arcsinh does not account for the asymmetric negative tails created by [compensation matrix subtraction](https://doi.org/10.1002/cyto.a.20258) | [Logicle Transform](logicle-transform.md) |
| Sample-level count data for GLMMs (cells per cluster per sample) | Integer counts need a discrete-distribution model, [not a continuous transform](https://doi.org/10.1038/nmeth.4295) | [Log Transform](log-transform.md) with offset, or use counts directly in a negative binomial model |
| Already pre-transformed data | Applying arcsinh to data that has already been transformed [distorts the scale](https://doi.org/10.12688/f1000research.11622.3) and compresses the dynamic range further | Skip transformation (`transform_type = "None"`) |
| Spectral flow cytometry with high-dimensional unmixing | Unmixed data may have [complex noise structures](https://doi.org/10.1002/cyto.a.20258) that benefit from instrument-specific transforms | [Logicle Transform](logicle-transform.md) or instrument vendor defaults |
| Strictly positive data with no zeros (e.g., gene expression RPKM) | The linear behavior near zero is unnecessary and may reduce resolution compared to log | [Log Transform](log-transform.md) |

---

## How to Use in the Pipeline

### Transform FCS Files (Pipeline A)

The primary entry point is `transform_flowframe()` in
`src/preprocessing/transformation.R`. The function takes per-marker offset and
width lists and applies the arcsinh transform via `flowCore::arcsinhTransform()`.

The parameterization in `flowCore` uses `arcsinhTransform(a, b)`, which
computes `arcsinh(a + b * x)`. In this pipeline, `a = -offset` and
`b = 1/width`, so with default offset = 0, the effective formula is
`arcsinh(x / width)`, where `width` is the cofactor.

```r
# Transform FCS files (Pipeline A)
source("src/preprocessing/transformation.R")

# Using cofactor 5 (standard for CyTOF): width = 5
transformed <- transform_flowframe(
  flowframe          = raw_flowframe,
  sample_name        = "sample_01",
  markers_of_interest = markers_to_transform,
  offset_list        = rep(0, length(markers_to_transform)),
  width_list         = rep(5, length(markers_to_transform)),
  transform_type     = "arcsinh"
)
```

### Auto-Optimize Cofactor

When the standard cofactor of 5 may not be appropriate (e.g., fluorescence
data, non-standard instrument settings), use `find_optimal_b_cofactor()` to
select the cofactor that maximizes bimodality:

```r
source("src/preprocessing/transformation.R")

# Extract expression values for one channel
channel_data <- flowCore::exprs(raw_flowframe)[, "CD4"]

# Find the optimal b parameter (reciprocal of cofactor)
optimal_b <- find_optimal_b_cofactor(channel_data, target_bimodality_coefficient = 0.9)

# Convert to cofactor for interpretation
cat("Optimal cofactor:", 1 / optimal_b, "\n")

# Apply per-marker: use optimal_b as the width parameter
# Note: width = 1/b because arcsinhTransform(b = 1/width) and
#        find_optimal_b_cofactor returns b directly
transformed <- transform_flowframe(
  flowframe          = raw_flowframe,
  sample_name        = "sample_01",
  markers_of_interest = markers_to_transform,
  offset_list        = rep(0, length(markers_to_transform)),
  width_list         = rep(1 / optimal_b, length(markers_to_transform)),
  transform_type     = "arcsinh"
)
```

### Per-Marker Cofactor Optimization

For panels with diverse marker types (e.g., surface markers with high
expression versus intracellular markers with low expression), optimize the
cofactor independently for each channel:

```r
source("src/preprocessing/transformation.R")

markers <- c("CD4", "CD8a", "CD45", "FoxP3", "Ki67")
optimal_widths <- numeric(length(markers))

for (i in seq_along(markers)) {
  channel_data <- flowCore::exprs(raw_flowframe)[, markers[i]]
  b <- find_optimal_b_cofactor(channel_data)
  optimal_widths[i] <- 1 / b
  cat(markers[i], "- optimal cofactor:", optimal_widths[i], "\n")
}

transformed <- transform_flowframe(
  flowframe          = raw_flowframe,
  sample_name        = "sample_01",
  markers_of_interest = markers,
  offset_list        = rep(0, length(markers)),
  width_list         = optimal_widths,
  transform_type     = "arcsinh"
)
```

### Pipeline Configuration

In the config-driven pipeline orchestrators (`pipelines/`), set the
transformation parameters in the config file:

```yaml
preprocessing:
  transform_type: "arcsinh"
  cofactor: 5                # Applied as width; use 150 for fluorescence
  optimize_cofactor: false   # Set true to run find_optimal_b_cofactor()
  target_bimodality: 0.9     # Target BC for optimization
```

---

## Effect on Downstream Analysis

The choice of transform and cofactor directly influences every downstream
step in the pipeline.

### Leiden Clustering

[Leiden clustering](leiden-clustering.md) operates on a shared nearest
neighbor (SNN) graph constructed from Euclidean distances between cells in
marker space. The arcsinh transform
[equalizes the contribution of each marker to inter-cell distances](https://doi.org/10.12688/f1000research.11622.3)
by bringing all channels onto a comparable scale. Without transformation,
markers with higher absolute intensities dominate the distance computation
and drive cluster structure.

Weber and Robinson (2016) showed that
[clustering performance varies with transformation parameters](https://doi.org/10.1002/cyto.a.22764),
with cofactor values that are too small over-compressing high-expression
markers (merging distinct populations) and values that are too large
under-compressing them (fragmenting coherent populations into multiple
clusters). The standard cofactor of 5 for CyTOF was selected in part because
it produces
[biologically interpretable clusters that match expert manual gating](https://doi.org/10.12688/f1000research.11622.3).

### UMAP

[UMAP](umap.md) dimensionality reduction computes distances between cells in
high-dimensional marker space and preserves local neighborhood structure in
the 2D embedding. Because UMAP is
[sensitive to the relative scale of input features](https://doi.org/10.12688/f1000research.11622.3),
applying the arcsinh transform before UMAP ensures that all markers
contribute to the embedding rather than being dominated by one or two
high-intensity channels. The
[Nowicka et al. (2017) workflow](https://doi.org/10.12688/f1000research.11622.3)
applies arcsinh transformation as a mandatory step before any dimensionality
reduction.

### PCA

[PCA](pca.md) identifies axes of maximum variance in marker space. On
untransformed CyTOF data, the first principal components would be driven by
markers with the largest absolute variance, which are not necessarily the most
biologically informative. The arcsinh transform acts as a
[variance-stabilizing transformation](https://doi.org/10.12688/f1000research.11622.3)
that reduces this dominance, allowing PCA to capture biologically relevant
sources of variation rather than scale artifacts.

### Differential Abundance Analysis

After clustering, differential abundance analysis (Fisher's exact test,
mixed-effects logistic regression, negative binomial GLMM) operates on cell
counts per cluster and is therefore
[independent of the transformation applied to expression values](https://doi.org/10.1038/s42003-019-0415-5).
However, the transform indirectly affects abundance results because it
influences cluster assignments: different cofactors produce different cluster
boundaries, which change the cell counts per cluster.

### Marker Expression Statistics

When computing median marker expression per cluster (for heatmaps, violin
plots, or statistical tests), the values are reported on the transformed
scale unless explicitly back-transformed. The
[diffcyt framework](https://doi.org/10.1038/s42003-019-0415-5) recommends
performing differential expression analysis on arcsinh-transformed medians
to satisfy the approximate normality assumptions of linear models.

---

## Comparison: Arcsinh vs Logicle vs Log

| Property | Arcsinh | Logicle | Log |
|----------|---------|---------|-----|
| **Formula** | arcsinh(*x*/*c*) | Parametric bi-exponential (T, W, M, A) | log(*x*) or log(*x* + pseudocount) |
| **Parameters** | 1 (cofactor) | 4 (T, W, M, A) | 0-1 (optional pseudocount) |
| **Handles *x* = 0** | Yes (maps to 0) | Yes | No (undefined; requires pseudocount) |
| **Handles *x* < 0** | Yes (odd function, symmetric) | Yes (designed for negative tails) | No (undefined) |
| **Behavior near zero** | [Linear](https://doi.org/10.1126/science.1198704) | Linear (parameterized width) | Undefined or distorted |
| **Behavior at high values** | Logarithmic | Logarithmic | Logarithmic |
| **Variance stabilization** | [Approximate, for Poisson-like data](https://doi.org/10.1002/cyto.a.22271) | Data-adaptive | [Exact for log-normal data](https://doi.org/10.1002/cyto.a.20258) |
| **Primary use case** | [CyTOF / mass cytometry](https://doi.org/10.1126/science.1198704) | [Compensated fluorescence cytometry](https://doi.org/10.1002/cyto.a.20258) | Positive-only data, count models |
| **Standardization** | [Cofactor 5 is universal for CyTOF](https://doi.org/10.12688/f1000research.11622.3) | Parameters are [estimated per channel per experiment](https://doi.org/10.1002/cyto.a.20258) | No standard parameterization |
| **In this pipeline** | `transform_type = "arcsinh"` | `transform_type = "logicle"` | Not directly supported for expression data; used for count models |
| **Cross-link** | This document | [Logicle Transform](logicle-transform.md) | [Log Transform](log-transform.md) |

### Decision Guide

1. **CyTOF data** -- use arcsinh with cofactor 5. This is the
   [established standard](https://doi.org/10.1126/science.1198704) and
   ensures comparability with published analyses.

2. **Compensated fluorescence data** -- use [logicle](logicle-transform.md).
   The logicle was
   [specifically designed for the negative tails produced by compensation](https://doi.org/10.1002/cyto.a.20258)
   and its four parameters adapt to each channel's noise characteristics.

3. **Fluorescence data without compensation** (rare, e.g., single-color
   controls or uncompensated exploratory panels) -- arcsinh with cofactor
   150 is a reasonable choice, as recommended by
   [Finck et al. (2013)](https://doi.org/10.1002/cyto.a.22271).

4. **Mixed panels or uncertain instrument type** -- use
   `find_optimal_b_cofactor()` to let the data determine the appropriate
   cofactor, then visually inspect the transformed distributions.

5. **Count data for statistical models** -- do not use arcsinh on the
   expression matrix. Instead, compute cluster abundances (cell counts) and
   model them with appropriate discrete distributions, such as the
   [negative binomial GLMM](negative-binomial-glmm.md) recommended by
   [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295).

---

## Implementation Details

### Source Code Location

- **Transform function:** [`transform_flowframe()`](/src/preprocessing/transformation.R)
  in `src/preprocessing/transformation.R`
- **Cofactor optimizer:** [`find_optimal_b_cofactor()`](/src/preprocessing/transformation.R)
  in `src/preprocessing/transformation.R`

### Parameterization Mapping

The `flowCore::arcsinhTransform()` function uses the formula:

$$
y = \operatorname{arcsinh}(a + b \cdot x)
$$

This pipeline maps to the standard cytometry formula `arcsinh(x / cofactor)`
by setting:

| Pipeline parameter | `flowCore` parameter | Relationship |
|--------------------|---------------------|--------------|
| `offset` | `a` | `a = -offset` |
| `width` (cofactor) | `b` | `b = 1 / width` |

With default offset = 0, this produces `arcsinh(x / width)`, which is
equivalent to `arcsinh(x / cofactor)` when `width = cofactor`.

### Error Handling

The `transform_flowframe()` function validates inputs before applying the
transform:

- Marker and offset lists must have the same length (raises `stop()` if not).
- Marker and width lists must have the same length (raises `stop()` if not).
- Markers not found in the flowFrame are skipped with a `warning()`.
- Invalid `transform_type` values raise `stop()`.

### Supported Transform Types

| `transform_type` | `flowCore` function | Parameters |
|-------------------|--------------------| -----------|
| `"arcsinh"` | `arcsinhTransform(a = -offset, b = 1/width)` | Per-marker offset and width |
| `"logicle"` | `logicleTransform()` | Default parameters |
| `"biexponential"` | `biexponentialTransform()` | Default parameters |
| `"None"` | No transform applied | -- |

---

## References by Topic

All references are linked inline throughout this document. The key publications
are:

- **CyTOF with arcsinh transform (foundational):** [Bendall, S. C. et al. (2011)](https://doi.org/10.1126/science.1198704). Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum. *Science*, 332(6030), 687-696.
- **CyTOF normalization and cofactor selection:** [Finck, R. et al. (2013)](https://doi.org/10.1002/cyto.a.22271). Normalization of mass cytometry data with bead standards. *Cytometry A*, 83(5), 483-494.
- **CyTOF analysis workflow:** [Nowicka, M. et al. (2017)](https://doi.org/10.12688/f1000research.11622.3). CyTOF workflow: Differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*, 6, 748.
- **Logicle transform (comparison):** [Parks, D. R., Roederer, M., & Moore, W. A. (2006)](https://doi.org/10.1002/cyto.a.20258). A new "logicle" display method avoids deceptive effects of logarithmic scaling for low signals and compensated data. *Cytometry A*, 69(6), 541-551.
- **Clustering benchmark for CyTOF:** [Weber, L. M. & Robinson, M. D. (2016)](https://doi.org/10.1002/cyto.a.22764). Comparison of clustering methods for high-dimensional single-cell flow and mass cytometry data. *Cytometry A*, 89(12), 1084-1096.
- **diffcyt framework:** [Weber, L. M., Nowicka, M., Soneson, C., & Robinson, M. D. (2019)](https://doi.org/10.1038/s42003-019-0415-5). diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. *Commun. Biol.*, 2, 183.
- **Differential abundance in mass cytometry:** [Lun, A. T. L., Richard, A. C., & Marioni, J. C. (2017)](https://doi.org/10.1038/nmeth.4295). Testing for differential abundance in mass cytometry data. *Nat. Methods*, 14(7), 707-709.
- **Bimodality coefficient:** [Pfister, R. et al. (2013)](https://doi.org/10.1080/00031305.2013.847865). Good things peak in pairs: A note on the bimodality coefficient. *Front. Psychol.*, 4, 700.
