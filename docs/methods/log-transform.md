# Logarithmic Transform for Count Data Modeling

## Summary

The logarithmic transform serves two distinct roles in this cytometry pipeline:
it acts as the canonical link function and offset mechanism in negative binomial
generalized linear mixed models (GLMMs) for
[differential abundance analysis](https://doi.org/10.1038/s42003-019-0415-5),
and it provides the base-2 logarithm used to express fold changes on a symmetric,
interpretable scale for volcano plots. Unlike the logicle or arcsinh transforms
used for fluorescence and mass cytometry intensity data, the log transform here
operates on **count data and ratios**, not on single-cell expression values.

---

## The Method

### Logarithms in Generalized Linear Models

The log function is the
[canonical link function for Poisson and negative binomial regression](https://doi.org/10.1007/978-1-4899-3242-6),
meaning it connects the linear predictor to the expected value of the response
through the relationship:

$$
\log(\mu_i) = \mathbf{x}_i^\top \boldsymbol{\beta}
$$

where $\mu_i$ is the expected count for observation $i$, $\mathbf{x}_i$ is
the vector of covariates, and $\boldsymbol{\beta}$ is the coefficient vector.
[McCullagh and Nelder (1989)](https://doi.org/10.1007/978-1-4899-3242-6)
established that the log link is canonical for the Poisson family because it
arises naturally from the exponential family form of the distribution, and the
same link extends to the negative binomial as an overdispersed generalization
of the Poisson.

The log link ensures that predicted counts are always non-negative -- a
critical constraint for count data -- because exponentiating any real-valued
linear predictor yields a positive number. This is why
[Agresti (2015)](https://doi.org/10.1002/0471249688) describes the log link as
the default choice for modeling count outcomes in the GLM framework.

### The Offset for Normalization

In cytometry, different samples yield different total cell counts due to
variation in sample preparation, acquisition time, and instrument throughput.
To compare cluster abundances fairly across samples, the model must normalize
for these differences. The log of total cells enters the model as an
**offset** -- a covariate with a fixed coefficient of 1:

$$
\log(\mu_{ij}) = \log(N_i) + \mathbf{x}_i^\top \boldsymbol{\beta}
$$

where $N_i$ is the total number of cells in sample $i$ and $\mu_{ij}$ is the
expected count for cluster $j$ in sample $i$. Rearranging:

$$
\log\!\left(\frac{\mu_{ij}}{N_i}\right) = \mathbf{x}_i^\top \boldsymbol{\beta}
$$

This shows that the model is effectively regressing on the log of the
**proportion** of cells in each cluster, which is the quantity of biological
interest. The offset technique is
[standard practice in count-based differential analysis](https://doi.org/10.1093/bioinformatics/btp616)
and was formalized by
[McCullagh and Nelder (1989)](https://doi.org/10.1007/978-1-4899-3242-6) as
the correct way to handle known exposure or size variables in log-linear
models. In RNA-seq, [Robinson et al. (2010)](https://doi.org/10.1093/bioinformatics/btp616)
used the same log-offset approach in edgeR to normalize for library size,
and [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) adopted
this convention for cytometry cell counts in the diffcyt framework.

### Log2 Fold Changes

The pipeline reports differential abundance results as log2 fold changes
(log2FC). The GLMM estimates a coefficient $\hat{\beta}$ on the natural log
scale; the conversion to log2 is:

$$
\text{log2FC} = \frac{\hat{\beta}}{\log(2)}
$$

This places effects on a scale where each unit corresponds to a doubling or
halving of relative abundance: a log2FC of +1 means the cluster proportion is
[twice as high in the test condition](https://doi.org/10.1093/bioinformatics/btp616),
a log2FC of -1 means it is halved, and a log2FC of 0 means no change. The
symmetric, interpretable nature of the log2 scale is why it is the standard
for volcano plots in both
[genomics](https://doi.org/10.1093/bioinformatics/btp616) and
[cytometry](https://doi.org/10.1038/s42003-019-0415-5).

---

## Log as GLM Offset

### Why an Offset and Not a Covariate

When `log(total_cells)` is included as a regular covariate, its coefficient is
free to take any value, meaning the model does not enforce proportional scaling
with sample size. Including it as an offset fixes the coefficient to exactly 1,
which is equivalent to stating that if sample size doubles, the expected count
for every cluster also doubles -- the definition of a
[constant-proportion null model](https://doi.org/10.1007/978-1-4899-3242-6).
This is the correct normalization assumption for differential abundance
testing, where the question is whether a cluster's proportion changes between
conditions, not whether its absolute count changes.

### How This Pipeline Implements It

In `src/analysis/differential_abundance.R`, the `build_abundance_matrix()`
function computes the log offset at line 86:

```r
abundance_df$log_total <- log(abundance_df$total_cells)
```

This is then passed to `fit_nb_glmm()` (from `src/R_shared/stats_utils.R`)
via the model formula at line 223:

```r
model <- fit_nb_glmm(
  count_data,
  formula = count ~ condition + offset(log_total) + (1 | sample_id)
)
```

The `offset()` term in the formula tells `lme4::glmer.nb()` to include
`log_total` with a fixed coefficient of 1, implementing the
[exposure normalization described by McCullagh and Nelder (1989)](https://doi.org/10.1007/978-1-4899-3242-6).

### Relationship to edgeR and diffcyt

The log-offset normalization strategy mirrors the approach used by
[Robinson et al. (2010)](https://doi.org/10.1093/bioinformatics/btp616) in
edgeR for RNA-seq, where the log of library size serves the same role. The
[diffcyt framework by Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5)
explicitly adopted the edgeR paradigm for cytometry, using log-transformed
total cell counts as offsets in negative binomial models to test for
differential cluster abundance. [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295)
further validated this approach for mass cytometry data, demonstrating that
sample-level count models with log offsets control false discovery rates
correctly when biological replicates are available.

---

## Log2 Fold Changes

### Interpretation

The log2 fold change is the fundamental effect size metric for differential
abundance. It quantifies the magnitude and direction of change on a scale
centered at zero:

| log2FC | Fold Change | Interpretation |
|--------|-------------|----------------|
| +2 | 4x increase | Cluster proportion quadrupled in test condition |
| +1 | 2x increase | Cluster proportion doubled |
| +0.585 | 1.5x increase | 50% increase in proportion |
| 0 | No change | Equal proportions |
| -0.585 | 1.5x decrease | 33% decrease in proportion |
| -1 | 2x decrease | Cluster proportion halved |
| -2 | 4x decrease | Cluster proportion quartered |

A key property of the log2 scale is **symmetry**: a doubling (+1) and a
halving (-1) are equidistant from zero.
[This symmetry makes it the natural scale for two-sided hypothesis testing](https://doi.org/10.1093/bioinformatics/btp616)
and for visualizations where both enrichment and depletion should be equally
visible.

### Conversion from Model Coefficients

The negative binomial GLMM uses the natural log link, so the raw coefficient
$\hat{\beta}$ represents the log (base $e$) ratio of expected proportions
between conditions. The pipeline converts this to log2 at line 236 of
`src/analysis/differential_abundance.R`:

```r
log2fc = est / log(2)
```

This is algebraically equivalent to:

$$
\text{log2FC} = \log_2(e^{\hat{\beta}}) = \log_2\!\left(\frac{\mu_{\text{test}}}{\mu_{\text{reference}}}\right)
$$

### Use in Volcano Plots

The `plot_volcano()` function in the same file plots log2FC on the x-axis
against $-\log_{10}(\text{adjusted } p)$ on the y-axis, following the
convention established in microarray analysis and now standard in
[high-dimensional cytometry](https://doi.org/10.1038/s42003-019-0415-5).
This layout places the most biologically meaningful clusters -- those with
both large effect sizes and high statistical significance -- in the upper
corners of the plot. The
[Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
recommends volcano plots as the primary summary visualization for differential
abundance results, precisely because they display both dimensions of evidence
simultaneously.

---

## When Your Data Calls For It

The log transform is the correct choice in the following situations within the
pipeline.

### Count Data Normalization in GLMMs

When you have cell counts per cluster per sample and want to test whether
cluster proportions differ between conditions, the log of total cells should
enter the model as an offset. This applies whenever you are running the
[negative binomial GLMM](negative-binomial-glmm.md) via
`run_differential_abundance()`. The log offset is
[the standard normalization strategy for count-based differential analysis](https://doi.org/10.1093/bioinformatics/btp616)
and is handled automatically by the pipeline.

### Fold Change Computation

Any time you need to express a ratio of abundances on a symmetric scale, log2
is the appropriate transform. The pipeline uses it for:

- **GLMM coefficient interpretation**: Converting natural-log coefficients to
  log2 fold changes for biological interpretation.
- **Enrichment ratios**: The `log2_enrichment` column in enrichment results
  (from `src/analysis/cluster_enrichment.R`) expresses the ratio of cluster
  proportions between conditions on the log2 scale.
- **Volcano plots**: Both `plot_volcano()` in this file and the enrichment
  volcano use log2FC as the x-axis.

### GLM with Log Link

Beyond the GLMM, any generalized linear model for count outcomes in the
pipeline uses the log link function. This includes the Poisson GLMM fallback
in `fit_nb_glmm()` (which activates when the negative binomial model fails to
converge) and any custom models a user might fit using the pipeline's count
data. The log link is
[the canonical choice for both Poisson and negative binomial families](https://doi.org/10.1007/978-1-4899-3242-6).

---

## When to Choose Something Else

The log transform is not appropriate for all data types encountered in
cytometry. The table below maps common situations to better alternatives.

| Data Type | Problem with Log | Better Transform | Reference |
|-----------|-----------------|------------------|-----------|
| Fluorescence intensity (flow cytometry) | [Data span negative values and zero due to compensation](https://doi.org/10.1002/cyto.a.20258); log is undefined for non-positive values | [Logicle Transform](logicle-transform.md) | [Parks et al. (2006)](https://doi.org/10.1002/cyto.a.20258) |
| Mass cytometry (CyTOF) intensity | Many true zeros from ion counting; [log compresses near-zero signals into negative infinity](https://doi.org/10.1038/nbt.2317) | [Arcsinh Transform](arcsinh-transform.md) | [Bendall et al. (2011)](https://doi.org/10.1126/science.1198704) |
| Continuous data with zeros | $\log(0) = -\infty$; adding a pseudocount is [arbitrary and distorts the distribution](https://doi.org/10.1038/nmeth.4302) | Arcsinh (handles zeros natively since $\text{arcsinh}(0) = 0$) | [Bates et al. (2014)](https://doi.org/10.1002/cyto.a.22445) |
| Proportion data (0 to 1) | Log of proportions near zero is highly volatile; [logit is the canonical link for proportional outcomes](https://doi.org/10.1007/978-1-4899-3242-6) | Logit transform or beta regression | [McCullagh and Nelder (1989)](https://doi.org/10.1007/978-1-4899-3242-6) |

### Key Principle

The log transform is designed for **positive count data and ratios**. If your
data can be zero, negative, or represent continuous intensity measurements,
a different transform will produce more reliable results.

---

## How to Use in the Pipeline

The log transform is applied internally by the differential abundance
analysis. Users do not need to apply it manually; the pipeline handles both
the log offset and the log2 fold change conversion.

### Running Differential Abundance

```r
# Log transform is used internally by differential abundance
source("src/analysis/differential_abundance.R")

results <- run_differential_abundance(
  data_dir       = "output/cluster_data/",
  conditions_csv = "metadata/conditions.csv",
  output_dir     = "output/diff_abundance/"
)

# results$log2fc contains log2 fold changes
# results$p_adjusted contains BH-corrected p-values
```

### What Happens Under the Hood

1. `build_abundance_matrix()` computes `log_total <- log(total_cells)` for
   each sample (natural log).
2. `fit_nb_glmm()` includes `offset(log_total)` in the model formula, which
   [normalizes for sample size](https://doi.org/10.1007/978-1-4899-3242-6)
   by fixing the offset coefficient to 1.
3. The GLMM coefficient $\hat{\beta}$ for the condition effect is extracted
   from the model summary.
4. The coefficient is converted to log2: `log2fc = est / log(2)`.
5. `plot_volcano()` displays the log2FC values against adjusted p-values.

### Inspecting Results

The output CSV (`glmm_results.csv`) contains both the raw natural-log
estimate and the log2 fold change:

| Column | Scale | Description |
|--------|-------|-------------|
| `estimate` | $\log_e$ | Raw GLMM coefficient (natural log scale) |
| `log2fc` | $\log_2$ | Fold change on the log2 scale |
| `odds_ratio` | Linear | $e^{\hat{\beta}}$ -- the multiplicative change in proportion |
| `p_value` | Probability | Wald test p-value from the GLMM |
| `p_adjusted` | Probability | [Benjamini-Hochberg](benjamini-hochberg.md) corrected p-value |

---

## Why Not Log for Expression Data

### The Zero Problem

Flow and mass cytometry expression data routinely contain zeros and negative
values. In flow cytometry, fluorescence compensation can produce
[negative values when spectral spillover is subtracted](https://doi.org/10.1002/cyto.a.20258),
and in CyTOF, the discrete ion-counting process means that
[many markers register exactly zero for a given cell](https://doi.org/10.1126/science.1198704).
The logarithm is undefined at zero ($\log(0) = -\infty$) and undefined for
negative values, making it unsuitable as a direct transform for expression
data.

### The Pseudocount Workaround and Its Problems

A common workaround is to add a small constant (pseudocount) before taking the
log: $\log(x + c)$. However, the choice of $c$ is arbitrary and
[materially affects downstream analysis](https://doi.org/10.1038/nmeth.4302).
Small values of $c$ amplify noise near zero; large values compress real
biological variation. There is no principled way to select $c$, which is why
[purpose-built transforms have replaced log+pseudocount in modern cytometry
workflows](https://doi.org/10.1002/cyto.a.20258).

### Compression Near Zero

Even when data are strictly positive, the log transform compresses values near
zero far more aggressively than values further from zero. For cytometry data
where the biologically meaningful distinction between dim-positive and negative
populations occurs precisely in the low-intensity range, this compression
[obscures the populations that matter most](https://doi.org/10.1002/cyto.a.20258).
The logicle transform was specifically designed to provide
[linear-like behavior near zero and logarithmic behavior at high intensities](https://doi.org/10.1002/cyto.a.20258),
preserving resolution across the full dynamic range.

### What to Use Instead

- **Fluorescence intensity data (flow cytometry)**: Use the
  [Logicle Transform](logicle-transform.md), a biexponential function that
  handles negative values and provides
  [optimal visualization and clustering of compensated data](https://doi.org/10.1002/cyto.a.20258).
- **Mass cytometry (CyTOF) intensity data**: Use the
  [Arcsinh Transform](arcsinh-transform.md) with a cofactor (typically 5),
  which maps zero to zero and approximates the log at high values while
  [remaining well-defined across the entire real line](https://doi.org/10.1126/science.1198704).

---

## Relationship to Other Methods

The log transform is one component of a larger statistical framework in this
pipeline. It connects to:

- **[Negative Binomial GLMM](negative-binomial-glmm.md)**: The log transform
  provides both the link function and the offset mechanism for the
  [negative binomial mixed model](https://doi.org/10.1038/s42003-019-0415-5)
  used in differential abundance testing.
- **[Benjamini-Hochberg Correction](benjamini-hochberg.md)**: Applied to the
  p-values produced by the GLMM to control the
  [false discovery rate](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
  when testing multiple clusters.
- **[Fisher's Exact Test](fisher-exact-test.md)**: An alternative to the GLMM
  for pooled or single-sample designs. Fisher's test also produces a log2
  enrichment ratio, but computes it directly from proportions rather than
  through a GLM with log link.
- **Volcano plots**: The joint display of log2FC and $-\log_{10}(p)$ is the
  standard output format for both the GLMM and enrichment analyses, following
  the convention described in the
  [diffcyt framework](https://doi.org/10.1038/s42003-019-0415-5) and the
  [Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3).

---

## Implementation Details

### Source Code Locations

- **Log offset computation**: [`build_abundance_matrix()`](/src/analysis/differential_abundance.R),
  line 86 in `src/analysis/differential_abundance.R`
- **GLMM with log offset**: [`run_differential_abundance()`](/src/analysis/differential_abundance.R),
  line 223 in `src/analysis/differential_abundance.R`
- **NB GLMM fitting**: [`fit_nb_glmm()`](/src/R_shared/stats_utils.R) in
  `src/R_shared/stats_utils.R`
- **Log2FC conversion**: line 236 in `src/analysis/differential_abundance.R`
- **Volcano plot**: [`plot_volcano()`](/src/analysis/differential_abundance.R)
  in `src/analysis/differential_abundance.R`

### Natural Log vs. Log2

The pipeline uses two bases in different contexts:

| Context | Base | Reason |
|---------|------|--------|
| GLM link / offset | $e$ (natural log) | [Canonical link for NB family](https://doi.org/10.1007/978-1-4899-3242-6); coefficients are log-rate ratios |
| Fold change reporting | 2 (log2) | [Each unit = doubling](https://doi.org/10.1093/bioinformatics/btp616); standard in genomics and cytometry |
| Volcano plot x-axis | 2 (log2) | Symmetric scale for visualization |
| Volcano plot y-axis | 10 ($-\log_{10} p$) | Convention for p-value display |

The conversion between bases is: $\log_2(x) = \log_e(x) / \log_e(2)$, which
is exactly what the pipeline computes at line 236:

```r
log2fc = est / log(2)
```

### Edge Cases

| Situation | Behavior | Reason |
|-----------|----------|--------|
| `total_cells = 0` | `log(0) = -Inf` | No cells in sample; the sample should be excluded upstream |
| Cluster count = 0 for all samples | Model may fail to converge | No variation to model; cluster is excluded by `min_cluster_size` filter |
| Very small counts (1-2 cells) | Large variance in log2FC estimates | [The negative binomial handles this better than Poisson](https://doi.org/10.1093/bioinformatics/btp616) due to the dispersion parameter, but interpret with caution |
| All samples have identical total cells | Offset is constant; equivalent to no normalization | Rare in practice; the model still runs correctly |

---

## References by Topic

All references are linked inline throughout this document. The key publications
are:

- **GLM theory and log link**: [McCullagh, P. & Nelder, J. A. (1989)](https://doi.org/10.1007/978-1-4899-3242-6). *Generalized Linear Models*, 2nd ed. Chapman & Hall.
- **edgeR and log-offset normalization**: [Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010)](https://doi.org/10.1093/bioinformatics/btp616). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139-140.
- **diffcyt framework for cytometry**: [Weber, L. M., Nowicka, M., Soneson, C., & Robinson, M. D. (2019)](https://doi.org/10.1038/s42003-019-0415-5). diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. *Commun. Biol.*, 2, 183.
- **Differential abundance in mass cytometry**: [Lun, A. T. L., Richard, A. C., & Marioni, J. C. (2017)](https://doi.org/10.1038/nmeth.4295). Testing for differential abundance in mass cytometry data. *Nat. Methods*, 14(7), 707-709.
- **CyTOF analysis workflow**: [Nowicka, M. et al. (2017)](https://doi.org/10.12688/f1000research.11622.3). CyTOF workflow: Differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*, 6, 748.
- **Logicle transform**: [Parks, D. R., Roederer, M., & Moore, W. A. (2006)](https://doi.org/10.1002/cyto.a.20258). A new "logicle" display method avoids deceptive effects of logarithmic scaling for low signals and compensated data. *Cytometry Part A*, 69A(6), 541-551.
- **CyTOF and arcsinh**: [Bendall, S. C. et al. (2011)](https://doi.org/10.1126/science.1198704). Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum. *Science*, 332(6030), 687-696.
- **FDR correction**: [Benjamini, Y. & Hochberg, Y. (1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B*, 57(1), 289-300.
- **Categorical data analysis**: [Agresti, A. (2015)](https://doi.org/10.1002/0471249688). *Categorical Data Analysis*, 3rd ed. Wiley.
- **Pseudocount problems**: [Lun, A. T. L. (2018)](https://doi.org/10.1038/nmeth.4302). Overcoming systematic errors caused by log-transformation of normalized single-cell RNA sequencing data. *Nat. Methods*, 15, 235-236.
