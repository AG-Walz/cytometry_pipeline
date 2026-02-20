# Mixed-Effects Logistic Regression for Cluster Enrichment

## Summary

Mixed-effects logistic regression is a generalized linear mixed model (GLMM)
that tests whether the probability of a cell belonging to a particular cluster
differs between experimental conditions, while accounting for the fact that
cells from the same sample are correlated. The pipeline implements it in
`src/R_shared/stats_utils.R` (`do_mixed_effects()`) using
[lme4::glmer()](https://doi.org/10.18637/jss.v067.i01), and it is called from
`src/analysis/cluster_enrichment.R` with `stat_method = "mixed_effects"`. This
is the recommended test for binary cluster-membership outcomes when multiple
biological replicates are available per condition.

---

## The Method

### Mathematical Formulation

The model treats each cell's cluster membership as a Bernoulli trial and links
it to the experimental condition through a logit function:

```
logit(P(Y_ij = 1)) = beta_0 + beta_condition * X_i + u_i

where:
  Y_ij     = 1 if cell j in sample i belongs to the cluster, 0 otherwise
  X_i      = condition indicator for sample i (0 or 1)
  u_i      ~ N(0, sigma^2_u)   random intercept for sample i
```

This is a two-level hierarchical model: cells (level 1) are nested within
samples (level 2). The fixed effect `beta_condition` captures the systematic
difference in log-odds of cluster membership between conditions. The random
intercept `u_i` captures sample-specific deviations -- some samples will have
inherently higher or lower baseline probabilities of cells falling into the
cluster, independent of condition.

### Components

| Component | Symbol | Role |
|-----------|--------|------|
| Fixed intercept | beta_0 | Baseline log-odds of cluster membership (reference condition) |
| Condition effect | beta_condition | Change in log-odds when moving from reference to treatment condition |
| Random intercept | u_i ~ N(0, sigma^2_u) | Sample-specific deviation from the population mean |
| Link function | logit | Maps probability [0,1] to linear predictor (-inf, +inf) |
| Distribution | Bernoulli | Each cell is independently in or out of the cluster (conditional on random effect) |

### Why a Logit Link?

The logit function is the canonical link for binomial data. It ensures that
predicted probabilities remain in [0, 1] and provides coefficients that are
directly interpretable as [log-odds ratios](https://doi.org/10.18637/jss.v067.i01).
For a binary outcome like cluster membership, the logit link arises naturally
from the [exponential family formulation of the binomial distribution](https://doi.org/10.1007/978-1-4899-3242-6).

### Estimation

The pipeline fits the model via [Laplace approximation](https://doi.org/10.18637/jss.v067.i01)
as implemented in `lme4::glmer()`. The Laplace approximation integrates over
the random effects by approximating the marginal likelihood with a second-order
Taylor expansion around the mode of the random effects. This is
[computationally efficient and accurate for single random-effect terms](https://doi.org/10.18637/jss.v067.i01),
though it can underestimate variance components in some edge cases. The pipeline
uses the BOBYQA optimizer with a maximum of 20,000 function evaluations:

```r
lme4::glmer(
  cluster ~ condition + (1 | file_name),
  data = model_data, family = binomial,
  control = lme4::glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 20000))
)
```

### Inference

After fitting, the pipeline extracts the fixed-effect estimate for condition
and performs a **Wald z-test**:

```
z = beta_condition / SE(beta_condition)
p = 2 * (1 - Phi(|z|))
```

where Phi is the standard normal CDF. The odds ratio is `exp(beta_condition)`,
and 95% Wald confidence intervals are computed as
`exp(beta_condition +/- 1.96 * SE)`. The Wald test is
[asymptotically valid for GLMMs](https://doi.org/10.18637/jss.v067.i01) and is
the standard inference approach in lme4 for logistic models, since likelihood
ratio tests for GLMMs require careful handling of boundary constraints on
variance parameters.

---

## Why Mixed Effects Matter for Cytometry

### The Pseudoreplication Problem

Flow and mass cytometry experiments produce thousands to millions of cells per
sample, but cells from the same sample share biological and technical
variation. When cells are pooled across samples and tested as if independent,
[the effective sample size is massively inflated](https://doi.org/10.1038/nmeth.4295),
leading to artificially small p-values. This is the classical problem of
**pseudoreplication** -- treating non-independent observations as independent
replicates.

[Lun, Richard, and Marioni (2017)](https://doi.org/10.1038/nmeth.4295) demonstrated
this directly in the cytometry context: naive tests that ignore sample structure
produce [dramatically inflated Type I error rates](https://doi.org/10.1038/nmeth.4295),
rejecting the null hypothesis far more often than the nominal rate would suggest.
A test at alpha = 0.05 might reject 30-50% of true nulls when within-sample
correlation is ignored.

This problem is [well-characterized in ecology and biology more broadly](https://doi.org/10.1016/j.tree.2009.03.008).
Bolker et al. (2009) provide a comprehensive treatment of why GLMMs are
necessary whenever observations are grouped within higher-level units --
whether those units are samples, patients, experimental batches, or field sites.

### The Hierarchical Structure of Cytometry Data

Cytometry data has an inherent multilevel structure:

```
Condition (e.g., treatment vs. control)
  |
  +-- Sample 1 (Patient A, treated)
  |     +-- Cell 1, Cell 2, ..., Cell 50000
  |
  +-- Sample 2 (Patient B, treated)
  |     +-- Cell 1, Cell 2, ..., Cell 30000
  |
  +-- Sample 3 (Patient C, control)
        +-- Cell 1, Cell 2, ..., Cell 45000
```

Cells within a sample share:

- **Biological variation**: The same immune state, disease stage, genetic
  background, and environmental exposures.
- **Technical variation**: The same staining batch, acquisition session,
  instrument settings, and compensation matrix.
- **Compositional coupling**: If one cell type expands within a sample, others
  contract proportionally.

[Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) formalized
this in their diffcyt framework, showing that differential testing in cytometry
must account for the [sample as the unit of replication](https://doi.org/10.1038/s42003-019-0415-5),
not the individual cell. The mixed-effects model achieves this by partitioning
variance into within-sample (residual) and between-sample (random effect)
components.

### What the Random Intercept Does

The random intercept `(1 | sample)` allows each sample to have its own
baseline log-odds of cluster membership. Consider two treated samples where
Patient A has 15% of cells in a cluster and Patient B has 8%. Without a random
intercept, all cells from both patients are pooled, and the model assumes they
come from a single population. With the random intercept, the model recognizes
that [Patient A has an inherently higher baseline for this cluster](https://doi.org/10.1016/j.tree.2009.03.008),
and the condition effect is estimated after accounting for this patient-level
variation.

This [correctly attributes variance to its source](https://doi.org/10.1038/nmeth.4295):
biological differences between patients go into sigma^2_u, and the condition
effect beta_condition reflects only the systematic difference between conditions
after removing patient-level noise. The result is [appropriately calibrated
p-values](https://doi.org/10.1038/s42003-019-0415-5) -- neither too liberal
(as with pooled tests) nor too conservative.

### Intraclass Correlation

The intraclass correlation coefficient (ICC) quantifies how much of the total
variance in cluster membership is attributable to differences between samples:

```
ICC = sigma^2_u / (sigma^2_u + pi^2/3)
```

where `pi^2/3` is the level-1 variance on the logistic scale (the variance of
the standard logistic distribution). A high ICC (e.g., > 0.1) indicates
substantial within-sample correlation, making mixed-effects modeling essential.
In cytometry data, ICCs are [typically moderate to high](https://doi.org/10.1038/nmeth.4295)
because cells from the same patient are biologically similar. The `pi^2/3`
denominator arises from the [latent-variable interpretation of logistic regression](https://doi.org/10.18637/jss.v067.i01),
where the binary outcome is modeled as a thresholded continuous variable with
logistic errors.

---

## When Your Data Calls For It

Mixed-effects logistic regression is appropriate when all of the following
conditions hold:

### 1. Multiple Biological Replicates Per Condition

You need at least two samples per condition for the random effect to be
estimable, and [at least 3-5 per condition for reliable inference](https://doi.org/10.1016/j.tree.2009.03.008).
With fewer samples, the variance of the random effect is poorly estimated, and
the Wald test may not be well-calibrated.

**Typical scenario**: Comparing immune cell phenotypes between responders (n=8)
and non-responders (n=10) in an immunotherapy cohort, where each patient
contributes one PBMC sample measured by flow cytometry.

### 2. Binary Outcome (Cluster Membership)

The outcome for each cell is binary: the cell either belongs to the cluster
(Y=1) or does not (Y=0). This is the natural outcome when working with
binary indicator columns from Leiden clustering or FlowJo gating, as produced
by `clusters_to_gating_csv()` in `src/clustering/leiden.R`.

### 3. Need to Account for Within-Sample Correlation

If you have reason to believe -- and in cytometry you always do -- that cells
from the same sample are more similar to each other than to cells from other
samples, the mixed-effects model is necessary. A Fisher's exact test on pooled
data [treats every cell as an independent observation](https://doi.org/10.1038/nmeth.4295),
which is statistically invalid when cells are nested within samples.

### 4. Sufficient Cluster Prevalence

The cluster must contain enough cells for the model to distinguish signal from
noise. The pipeline enforces a minimum cluster fraction (`min_cluster_frac`,
default 0.001) to filter out extremely rare clusters. If a cluster contains
zero cells in all samples of one condition, the model cannot estimate a
condition effect and returns `no_cluster_variation`.

### 5. Concrete Cytometry Use Cases

- **Immunotherapy response prediction**: Testing whether specific T-cell
  subsets (e.g., PD-1+ CD8+ T cells) are enriched in responders vs.
  non-responders, with multiple patients per group.
- **Disease phenotyping**: Comparing the frequency of regulatory T cells
  (Tregs) between healthy donors and autoimmune disease patients.
- **Vaccine studies**: Evaluating whether a particular B-cell population
  expands post-vaccination relative to baseline, across multiple subjects.
- **Tumor microenvironment**: Assessing whether tumor-infiltrating myeloid
  populations differ between cancer subtypes, with multiple tumor samples per
  subtype.

---

## When to Choose Something Else

### One Sample Per Condition --> [Fisher's Exact Test](fisher-exact-test.md)

If you have only a single sample per condition (or are deliberately pooling all
samples), the random intercept cannot be estimated. Use Fisher's exact test,
which compares pooled proportions in a 2x2 contingency table. This is
statistically valid when there truly is only one observation per group, but
[it cannot account for biological variability between subjects](https://doi.org/10.1038/nmeth.4295).

In the pipeline:

```r
results <- run_cluster_enrichment(data_dir, conditions_csv,
                                   stat_method = "pooled")
```

### Count-Based (Sample-Level) Analysis --> [Negative Binomial GLMM](negative-binomial-glmm.md)

If you want to model the *number* of cells per cluster per sample rather than
cell-level binary outcomes, use the negative binomial GLMM. This approach:

- Operates at the **sample level** (one observation per sample per cluster),
  making it computationally cheaper for large datasets.
- [Naturally handles overdispersion](https://doi.org/10.1038/s42003-019-0415-5)
  through the negative binomial variance function.
- Uses a log-link with an offset for total cell count, making it a model for
  **proportions** on the log scale.

[Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) recommend
this approach in the diffcyt framework for differential abundance testing,
as it [avoids the computational burden of fitting cell-level models](https://doi.org/10.1038/s42003-019-0415-5)
while properly accounting for sample-level variability.

In the pipeline:

```r
source("src/analysis/differential_abundance.R")
results <- run_differential_abundance(data_dir, conditions_csv, output_dir)
```

### Continuous Expression Differences --> [Wilcoxon Rank-Sum](wilcoxon-rank-sum.md)

If your question is about whether a marker's *expression level* differs between
conditions (not whether cells belong to a cluster), use the Wilcoxon rank-sum
test or a mixed-effects linear model. Mixed-effects logistic regression is
designed for binary outcomes, not continuous expression intensities.

In the pipeline:

```r
source("src/analysis/marker_statistics.R")
results <- run_marker_statistics(data_dir, conditions_csv,
                                  test_type = "wilcox")
```

### Model Fails to Converge --> Fallback Strategies

See the [Convergence Issues](#convergence-issues) section below.

---

## Assumptions and Diagnostics

### 1. Conditional Independence

**Assumption**: Given the random effect u_i, cells within a sample are
independent Bernoulli trials. This means that after accounting for the sample's
overall propensity to have cells in the cluster, the remaining cell-to-cell
variation is purely random.

**When it breaks**: When there is additional structure within samples -- for
example, cells from the same spatial region of a tissue, or cells measured at
different time points within the same acquisition. In mass cytometry,
[acquisition effects can create temporal autocorrelation within a sample](https://doi.org/10.12688/f1000research.11622.3),
violating conditional independence.

**Diagnostic**: Plot residuals against acquisition order within each sample.
If systematic trends appear, consider adding a within-sample covariate or
using robust standard errors. In practice, the random intercept absorbs much
of this variation, and [the Wald test remains approximately valid for moderate
violations](https://doi.org/10.1016/j.tree.2009.03.008).

### 2. Normality of Random Effects

**Assumption**: The sample-level random intercepts follow a normal distribution:
u_i ~ N(0, sigma^2_u). This assumption is about the [distribution of
sample-level effects across the population of samples](https://doi.org/10.18637/jss.v067.i01),
not about the distribution of individual cells.

**When it breaks**: With very few samples (< 10 total), normality cannot be
reliably assessed. Outlier samples with extreme cluster frequencies can also
violate this assumption. [Bolker et al. (2009)](https://doi.org/10.1016/j.tree.2009.03.008)
note that GLMMs are moderately robust to departures from normality in the
random effects, especially when the number of groups is not too small.

**Diagnostic**: Extract the estimated random effects with `lme4::ranef(model)`
and create a Q-Q plot against the normal distribution. Substantial departures
(heavy tails, skewness, bimodality) suggest the normality assumption is
problematic.

### 3. Correct Specification of the Linear Predictor

**Assumption**: The logit of the probability is a linear function of the
condition variable. There is no interaction, no nonlinear dose-response, and
no omitted confounders.

**When it breaks**: When important covariates (batch, patient age, disease
severity) are omitted, [the condition effect can be confounded](https://doi.org/10.1016/j.tree.2009.03.008).
The model as implemented in the pipeline includes only condition as a fixed
effect. If batch effects are present and not corrected upstream (via Harmony or
similar), they may bias the condition estimate.

**Diagnostic**: Check whether including additional fixed effects (e.g., batch)
meaningfully changes the condition estimate. If it does, the simpler model is
misspecified.

### 4. Sufficient Variance in Both Outcome and Predictor

**Assumption**: There must be variation in cluster membership across samples
and variation in the condition variable. If a cluster is present in all cells
or absent in all cells (`cluster_sum == 0 || cluster_sum == nrow(model_data)`),
the model is [undefined and the pipeline returns early](https://doi.org/10.18637/jss.v067.i01)
with `method = "no_cluster_variation"`.

**When it breaks**: Extremely rare clusters (< 0.1% of cells) or clusters
that are nearly universal (> 99.9% of cells) provide insufficient signal for
the model. Similarly, if all samples come from the same condition (only one
level present), the condition effect is not identifiable.

**Diagnostic**: Before running the model, check the prevalence of the cluster
in each sample. The pipeline's `min_cluster_frac` filter handles the most
extreme cases.

### 5. No Complete Separation

**Assumption**: The outcome should not be perfectly predicted by the fixed
effects. If all cells in condition A are in the cluster and no cells in
condition B are, the maximum likelihood estimate for `beta_condition` is
[infinite, and the model fails to converge](https://doi.org/10.18637/jss.v067.i01).

**When it breaks**: With small sample sizes and rare clusters, complete or
quasi-complete separation can occur. [Penalized likelihood methods (e.g.,
Firth's bias reduction)](https://doi.org/10.1002/sim.1047) can handle this,
but they are not implemented in the current pipeline.

**Diagnostic**: The pipeline catches convergence failures via `tryCatch` and
returns `method = "mixed_effects_error"` with the error message.

---

## How to Use in the Pipeline

### From Raw FCS Files (Pipeline A)

After running the full Pipeline A through clustering, the cluster assignments
are stored as binary columns in per-cell CSV files. Run enrichment testing:

```r
source("src/analysis/cluster_enrichment.R")

results <- run_cluster_enrichment(
  data_dir       = "output/cluster_data/",
  conditions_csv = "metadata/conditions.csv",
  output_dir     = "output/enrichment/",
  stat_method    = "mixed_effects",
  condition_levels = c("responder", "non-responder"),
  min_cluster_frac = 0.001,
  p_adjust_method  = "BH"
)
```

### From CSV Gating Data (Pipeline C)

If cluster membership comes from pre-existing gating (e.g., FlowJo export or
manual gates), the same function applies:

```r
source("src/analysis/cluster_enrichment.R")

results <- run_cluster_enrichment(
  data_dir       = "data/gating_csvs/",
  conditions_csv = "metadata/conditions.csv",
  output_dir     = "output/enrichment/",
  stat_method    = "mixed_effects",
  condition_levels = c("pCR", "npCR")
)
```

### Input File Format

**Per-cell CSV files** (in `data_dir`): Each row is a cell, with binary columns
indicating cluster membership (1 = in cluster, 0 = not). Column names must
match `cluster_regex` (default: `"Clustering annotated"`).

```
file_name,CD3,CD4,...,Clustering annotated_T_helper,Clustering annotated_Treg
sample_01.fcs,1.23,2.45,...,1,0
sample_01.fcs,0.89,1.67,...,0,1
```

**Conditions CSV** (at `conditions_csv`): Maps file names to conditions.

```
file_name;condition
sample_01.fcs;responder
sample_02.fcs;non-responder
```

### Output

The function returns a data frame and saves `enrichment_results.csv` with
these columns:

| Column | Description |
|--------|-------------|
| cluster | Cleaned cluster name |
| in_A, not_in_A, in_B, not_in_B | 2x2 table counts (from Fisher's test) |
| odds_ratio, p_value | Fisher's exact test results |
| me_estimate | Mixed-effects beta_condition (log-odds) |
| me_std_error | Standard error of the estimate |
| me_p_value | Wald z-test p-value |
| me_odds_ratio | exp(beta_condition) |
| me_ci_lower, me_ci_upper | 95% Wald confidence interval for odds ratio |
| me_method | "mixed_effects" or error description |
| p_adjusted | BH-corrected p-value (from Fisher by default) |
| significant | TRUE if p_adjusted < 0.05 |
| log2_enrichment | log2 of the ratio of proportions |

---

## Convergence Issues

### Why GLMMs Fail to Converge in Cytometry

[GLMM convergence failure is one of the most common practical problems](https://doi.org/10.1016/j.tree.2009.03.008)
when applying mixed-effects models to cytometry data. The causes are specific
to the structure of this data type:

1. **Rare clusters**: When a cluster contains < 1% of cells, many samples may
   have zero or near-zero membership. The log-odds become extreme, and the
   optimizer struggles to find a stable solution. [Sparse binary outcomes are a
   known challenge for logistic GLMMs](https://doi.org/10.18637/jss.v067.i01).

2. **Near-complete separation**: If the cluster is present almost exclusively
   in one condition, the condition coefficient approaches infinity. The
   [maximum likelihood estimate does not exist in this case](https://doi.org/10.1002/sim.1047),
   and iterative algorithms fail to converge.

3. **Small number of groups**: With fewer than 5-6 samples total, the
   random-effect variance is estimated from very few data points.
   [The variance component may be estimated as zero](https://doi.org/10.1016/j.tree.2009.03.008)
   (a boundary solution), which can cause optimizer warnings.

4. **Large datasets**: With hundreds of thousands of cells, the Laplace
   approximation can be [slow and memory-intensive](https://doi.org/10.18637/jss.v067.i01),
   and the optimizer may hit iteration limits.

### How the Pipeline Handles Convergence Failure

The `do_mixed_effects()` function wraps the model fit in `tryCatch`:

```r
tryCatch({
  model <- lme4::glmer(...)
  # ... extract results
}, error = function(e) {
  tibble::tibble(
    estimate = NA, std_error = NA, p_value = NA,
    odds_ratio = NA, ci_lower = NA, ci_upper = NA,
    method = paste0("mixed_effects_error: ", e$message)
  )
})
```

When the model fails, the error message is preserved in the `method` column,
and all result columns are set to NA. The enrichment analysis continues to the
next cluster without stopping.

### Strategies for Resolving Convergence Issues

| Strategy | When to Use | Implementation |
|----------|-------------|----------------|
| Increase iterations | Optimizer runs out of steps | Already set to `maxfun = 20000` in the pipeline |
| Switch optimizer | BOBYQA fails | Try `optimizer = "nloptwrap"` or `optimizer = "Nelder_Mead"` in `glmerControl` |
| Simplify the model | Very few samples | Remove random effect and use Fisher's exact test instead |
| Filter rare clusters | Cluster prevalence < 0.1% | Increase `min_cluster_frac` parameter |
| Use penalized estimation | Near-complete separation | Use [blme::bglmer()](https://doi.org/10.18637/jss.v067.i01) which adds a weakly informative prior on the fixed effects |
| Downsample cells | Memory/speed issues | Subsample cells per sample before fitting (retains variance structure) |
| Fall back to Fisher's test | Persistent failure | The pipeline always computes Fisher alongside mixed-effects; use Fisher p-values when ME fails |

[Bolker et al. (2009)](https://doi.org/10.1016/j.tree.2009.03.008) provide
extensive guidance on diagnosing and resolving GLMM convergence problems. Their
recommendation is to start simple (random intercepts only, as this pipeline
does), increase complexity only if justified, and always compare estimates
across optimizers to check robustness.

---

## Interpreting Results

### Odds Ratios

The primary effect size is the **odds ratio** (OR), computed as
`exp(beta_condition)`:

| Odds Ratio | Interpretation |
|------------|----------------|
| OR = 1 | No difference between conditions |
| OR > 1 | Cluster is enriched in the second condition level |
| OR < 1 | Cluster is enriched in the first (reference) condition level |
| OR = 2 | The odds of being in the cluster are 2x higher in the second condition |
| OR = 0.5 | The odds of being in the cluster are halved in the second condition |

The reference condition is the first element of `condition_levels`. In the
pipeline, `condition` is coded as a factor with the first level as reference,
so `beta_condition` represents the effect of moving from the first to the
second level.

**Important**: Odds ratios from logistic regression are
[non-collapsible](https://doi.org/10.1093/aje/kwt245) -- the marginal odds
ratio (from a model without random effects) will generally differ from the
conditional odds ratio (from the mixed-effects model), even when there is no
confounding. The conditional OR from the GLMM is the correct effect size when
the goal is to estimate the [within-sample effect of condition on cluster
membership](https://doi.org/10.1016/j.tree.2009.03.008).

### Wald Z-Test

The Wald test statistic is z = beta / SE(beta). Under the null hypothesis
(beta_condition = 0), z follows a standard normal distribution asymptotically.
The pipeline computes a two-sided p-value:

```r
z_val <- est / se_val
p_val <- 2 * (1 - pnorm(abs(z_val)))
```

The Wald test is [the default inference method in lme4 for logistic GLMMs](https://doi.org/10.18637/jss.v067.i01).
It is computationally fast but can be
[unreliable when the estimate is far from zero](https://doi.org/10.1016/j.tree.2009.03.008)
(i.e., for very large effects or near-separation). For borderline results,
consider computing profile confidence intervals with
`confint(model, method = "profile")`, which are more accurate but slower.

### Confidence Intervals

The 95% Wald confidence interval for the odds ratio is:

```
CI = exp(beta_condition +/- 1.96 * SE(beta_condition))
```

If the interval excludes 1, the result is significant at alpha = 0.05. The
interval on the log-odds scale is symmetric around the estimate; the
exponentiated interval on the odds-ratio scale is asymmetric (wider on the
upper end). [Profile likelihood confidence intervals](https://doi.org/10.18637/jss.v067.i01)
are an alternative that does not assume symmetry and is more accurate for small
samples, at higher computational cost.

### Multiple Testing

The pipeline applies [Benjamini-Hochberg correction](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
across all clusters to control the false discovery rate. The `p_adjusted`
column and `significant` flag in the output are based on the BH-corrected
p-values. When interpreting results, focus on `p_adjusted < 0.05` rather than
raw p-values.

---

## Comparison with Other Approaches

### Mixed-Effects Logistic Regression vs. Fisher's Exact Test

| Aspect | Mixed-Effects Logistic | Fisher's Exact |
|--------|----------------------|----------------|
| Unit of analysis | Cell (binary) | Cell (binary, pooled) |
| Sample structure | Modeled via random intercept | Ignored -- [cells pooled across samples](https://doi.org/10.1038/nmeth.4295) |
| Type I error control | [Correct when samples are the replicates](https://doi.org/10.1038/nmeth.4295) | [Inflated when samples have correlated cells](https://doi.org/10.1038/nmeth.4295) |
| Minimum samples | 3-5 per condition | 1 per condition |
| Computational cost | Higher (iterative GLMM fitting) | Negligible |
| Covariates | Can include additional fixed effects | Cannot |
| Effect size | Conditional odds ratio | Marginal odds ratio |

**Recommendation**: Use Fisher's exact test only when you have a single sample
per condition or as a quick exploratory screen. Use mixed-effects logistic
regression for any confirmatory analysis with biological replicates.

### Mixed-Effects Logistic Regression vs. Negative Binomial GLMM

| Aspect | Mixed-Effects Logistic | Negative Binomial GLMM |
|--------|----------------------|----------------------|
| Unit of analysis | Cell (binary membership) | Sample (count per cluster) |
| Outcome variable | 0/1 per cell | Integer count per sample |
| Number of observations | N_cells (can be millions) | N_samples (typically 10-50) |
| Computational cost | [High for large datasets](https://doi.org/10.18637/jss.v067.i01) | Low (one observation per sample) |
| Overdispersion handling | Not applicable (binary data) | [Explicit via NB variance function](https://doi.org/10.1038/s42003-019-0415-5) |
| Framework precedent | Used in single-cell RNA-seq GLMM approaches | [diffcyt standard](https://doi.org/10.1038/s42003-019-0415-5) |
| Effect size | Odds ratio (conditional) | Rate ratio / log2 fold change |

**Recommendation**: The [NB-GLMM is generally preferred for differential
abundance testing](https://doi.org/10.1038/s42003-019-0415-5) because it
operates at the sample level, is computationally efficient, and naturally
models count data. The mixed-effects logistic model is appropriate when you
specifically want a cell-level model or when additional cell-level covariates
must be included.

### Mixed-Effects Logistic Regression vs. Permutation Tests

| Aspect | Mixed-Effects Logistic | Permutation Test |
|--------|----------------------|------------------|
| Distributional assumptions | Bernoulli outcome, normal random effects | None (distribution-free) |
| Computational cost | Moderate (GLMM fit per cluster) | High (many permutations per cluster) |
| Effect size | Odds ratio with CI | Typically not available |
| Flexibility | Can include covariates and multiple random effects | Requires careful design of permutation scheme |
| Small-sample behavior | [Can be unreliable with < 5 groups](https://doi.org/10.1016/j.tree.2009.03.008) | Valid for any sample size if permutation scheme is correct |

**Recommendation**: Permutation tests are a useful [non-parametric
alternative](https://doi.org/10.12688/f1000research.11622.3) when you distrust
the distributional assumptions of the GLMM or have very few samples. However,
they do not provide effect sizes or confidence intervals and require careful
thought about [what to permute](https://doi.org/10.1038/nmeth.4295) (samples,
not individual cells) to respect the hierarchical structure.

### Summary Decision Table

| Scenario | Recommended Method | Pipeline Function |
|----------|-------------------|-------------------|
| 1 sample/condition, binary cluster | Fisher's exact test | `run_cluster_enrichment(..., stat_method = "pooled")` |
| Multiple samples/condition, binary cluster | **Mixed-effects logistic** | `run_cluster_enrichment(..., stat_method = "mixed_effects")` |
| Multiple samples/condition, count data | Negative binomial GLMM | `run_differential_abundance(...)` |
| Continuous expression, 2 groups | Wilcoxon rank-sum | `run_marker_statistics(..., test_type = "wilcox")` |
| GLMM convergence failure | Fisher's (fallback) | Automatic in `run_cluster_enrichment()` |

---

## Implementation Details

### Source Code

The core function lives in
[`src/R_shared/stats_utils.R`](/src/R_shared/stats_utils.R) (`do_mixed_effects()`).
The orchestrating function is in
[`src/analysis/cluster_enrichment.R`](/src/analysis/cluster_enrichment.R)
(`run_cluster_enrichment()`).

### Pre-checks Before Model Fitting

The `do_mixed_effects()` function performs two validation checks before
attempting to fit the model:

1. **Sufficient conditions**: At least two distinct levels of `condition` must
   be present. If not, the function returns `method = "insufficient_conditions"`.

2. **Non-degenerate outcome**: The cluster column must have both 0 and 1
   values. If all cells are in the cluster or none are, the function returns
   `method = "no_cluster_variation"`.

### Resume Support

`run_cluster_enrichment()` supports resumable computation via the `resume`
parameter. Intermediate results are saved every 10 clusters to an RDS file. If
the analysis is interrupted, re-running with `resume = TRUE` skips
already-completed clusters. This is important for large panels with hundreds
of clusters where the GLMM fit may take hours.

### Dependencies

- [lme4](https://doi.org/10.18637/jss.v067.i01) (>= 1.1-21): Core GLMM
  fitting engine
- dplyr: Data manipulation within `do_mixed_effects()`
- tibble: Structured output
- data.table: High-performance CSV I/O in the orchestrator
