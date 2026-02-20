# Negative Binomial Generalized Linear Mixed Model (NB-GLMM)

## Summary

The negative binomial GLMM is the primary statistical method in this pipeline for testing whether cell population abundances differ between experimental conditions. It models per-cluster cell counts across biological replicates, accounting for [overdispersion](https://doi.org/10.1007/978-1-4899-3242-6) (variance exceeding the mean) and [sample-level correlation structure](https://doi.org/10.18637/jss.v067.i01) via random intercepts. The implementation lives in `src/R_shared/stats_utils.R` (`fit_nb_glmm()`) and is orchestrated by `src/analysis/differential_abundance.R` (`run_differential_abundance()`), producing log2 fold changes and Benjamini-Hochberg-adjusted p-values that feed directly into volcano plots.

---

## The Method

### Mathematical Formulation

For cluster *k* in sample *j* under condition *i*, let $Y_{ijk}$ denote the observed cell count. The model specifies:

```
Y_ijk ~ NegBin(mu_ijk, theta)

log(mu_ijk) = beta_0 + beta_condition * condition_i + log(total_cells_j) + b_j

b_j ~ N(0, sigma^2)
```

where:

| Symbol | Meaning |
|--------|---------|
| `Y_ijk` | Observed cell count for cluster *k* in sample *j* |
| `mu_ijk` | Expected count (conditional mean) |
| `theta` | Negative binomial dispersion parameter (estimated from data) |
| `beta_0` | Intercept: baseline log expected proportion in reference condition |
| `beta_condition` | Fixed effect of condition: the log fold change on the natural scale |
| `log(total_cells_j)` | Offset term normalizing for differences in total cells acquired per sample |
| `b_j` | Random intercept for sample *j*, capturing sample-level variability |
| `sigma^2` | Between-sample variance of the random intercepts |

In R formula syntax, as implemented in `run_differential_abundance()`:

```r
count ~ condition + offset(log_total) + (1 | sample_id)
```

### The Negative Binomial Distribution

The [negative binomial distribution](https://doi.org/10.1007/978-1-4899-3242-6) is a two-parameter discrete distribution for count data. Its probability mass function is:

```
P(Y = y) = C(y + theta - 1, y) * (mu / (mu + theta))^y * (theta / (mu + theta))^theta
```

The mean is `E[Y] = mu` and the variance is `Var(Y) = mu + mu^2 / theta`. When `theta` approaches infinity, the variance converges to `mu` and the distribution reduces to a Poisson. Smaller `theta` values indicate greater overdispersion. The [diffcyt framework](https://doi.org/10.1038/s42003-019-0415-5) established that this parameterization is well-suited to cytometry cluster counts, where biological and technical variability between samples routinely inflates the variance beyond the Poisson assumption.

### The Log Link and Its Interpretation

The log link function connects the linear predictor to the expected count on a multiplicative scale. The coefficient `beta_condition` represents the [log ratio of expected counts between conditions](https://doi.org/10.1007/978-1-4899-3242-6), controlling for total cell numbers and sample-level effects. Exponentiating `beta_condition` gives the fold change: the multiplicative factor by which the expected proportion changes from the reference condition to the treatment condition. The pipeline converts this to a log2 fold change via `log2fc = estimate / log(2)` for compatibility with standard volcano plot conventions used in genomics and [cytometry differential analysis](https://doi.org/10.1038/s42003-019-0415-5).

### Random Intercepts

The term `(1 | sample_id)` introduces a [random intercept per sample](https://doi.org/10.18637/jss.v067.i01), modeled as normally distributed deviations from the population mean. This accounts for the fact that cells from the same sample are not independent observations. Without this term, the model would treat every cell count as an independent draw, dramatically underestimating standard errors and inflating false positives. [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) demonstrated that failing to account for this sample-level correlation is one of the most common sources of spurious findings in cytometry differential analysis. The random intercept captures inter-sample heterogeneity in overall cell capture efficiency, staining intensity, and biological donor effects.

---

## Why Negative Binomial Over Poisson

### The Overdispersion Problem

In a Poisson distribution, the variance equals the mean. Cytometry cell counts [systematically violate this assumption](https://doi.org/10.1038/s42003-019-0415-5). When you count cells in a cluster across replicate samples from the same condition, the observed variance in those counts is almost always larger than the mean. This phenomenon, called overdispersion, arises from multiple sources:

- **Biological heterogeneity**: Donors differ in immune composition due to genetics, prior antigen exposure, circadian rhythms, and countless unmeasured factors. [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) showed that this donor-to-donor variation is the dominant source of extra-Poisson variability in cytometry experiments.
- **Technical variability**: Differences in sample preparation, staining efficiency, acquisition time, and instrument performance introduce [systematic count variation beyond Poisson noise](https://doi.org/10.12688/f1000research.11622.3).
- **Compositional effects**: When one population expands, others must contract as a fraction of total cells. This induces [negative correlations between cluster counts](https://doi.org/10.1038/nmeth.4295) that the Poisson model cannot represent.

### Consequences of Ignoring Overdispersion

Fitting a Poisson model to overdispersed count data produces standard errors that are [too small by a factor proportional to the square root of the dispersion ratio](https://doi.org/10.1007/978-1-4899-3242-6). In practice, this means:

- Confidence intervals are too narrow
- P-values are artificially small
- False positive rates far exceed the nominal alpha level

The `check_overdispersion()` function in `src/R_shared/stats_utils.R` tests for this by comparing the residual deviance to the residual degrees of freedom. A dispersion ratio substantially greater than 1 (the implementation uses a threshold of 1.5 for the manual check, or defers to the `performance` package) signals that a Poisson model is inappropriate.

### How the Negative Binomial Fixes It

The negative binomial adds the [dispersion parameter theta](https://doi.org/10.1007/978-1-4899-3242-6) that absorbs the extra variance: `Var(Y) = mu + mu^2 / theta`. This gives the model freedom to match the observed mean-variance relationship in the data. The `lme4::glmer.nb()` function [estimates theta jointly with the regression coefficients and random effects](https://doi.org/10.18637/jss.v067.i01) via penalized iteratively reweighted least squares. This approach is analogous to how [edgeR](https://doi.org/10.1093/bioinformatics/btp616) handles overdispersion in RNA-seq count data, and [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) explicitly adopted this parallel when designing the diffcyt framework for cytometry.

---

## The Offset Trick

### Why Not Divide by Total Cells?

A natural instinct is to convert counts to proportions (count / total_cells) and analyze the proportions. This discards information: a proportion of 0.10 based on 100 total cells carries far more uncertainty than 0.10 based on 100,000 cells. Moreover, [proportions are bounded between 0 and 1](https://doi.org/10.1007/978-1-4899-3242-6), violating the assumptions of linear models, and dividing by total cells changes the error distribution in ways that standard variance-stabilizing transforms cannot fully repair.

### How the Offset Works

Instead, we keep the raw count as the response and add `log(total_cells)` as a fixed [offset](https://doi.org/10.1007/978-1-4899-3242-6) in the model:

```r
count ~ condition + offset(log_total) + (1 | sample_id)
```

An offset is a predictor whose coefficient is constrained to 1. Algebraically:

```
log(E[count]) = beta_0 + beta_condition * condition + 1 * log(total_cells) + b_j
```

Rearranging:

```
log(E[count] / total_cells) = beta_0 + beta_condition * condition + b_j
```

This means the model is effectively estimating log proportions while maintaining the correct [count-based error structure](https://doi.org/10.1007/978-1-4899-3242-6). The negative binomial variance function applies to the count, not to a transformed proportion, preserving the natural heteroscedasticity of count data. This offset approach is standard in generalized linear models for rate data, as described by [McCullagh and Nelder (1989)](https://doi.org/10.1007/978-1-4899-3242-6), and was adopted for cytometry differential abundance by [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) and [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5).

### Implementation Detail

In `run_differential_abundance()`, the offset is precomputed and stored in the abundance data frame:

```r
abundance_df$log_total <- log(abundance_df$total_cells)
```

Then passed to `fit_nb_glmm()` as part of the formula. This avoids recomputing the logarithm at each GLMM iteration.

---

## When Your Data Calls For It

The NB-GLMM is appropriate when **all four** of the following conditions hold:

### 1. You Are Comparing Cluster Abundances Across Conditions

You have cell counts per cluster per sample, and the scientific question is whether specific populations are [differentially abundant between conditions](https://doi.org/10.1038/s42003-019-0415-5). For example: "Is the CD8+ T cell cluster enriched in treatment responders versus non-responders?"

### 2. You Have Multiple Biological Replicates Per Condition

The random intercept `(1 | sample_id)` requires [at least 3 samples per condition](https://doi.org/10.1038/nmeth.4295) to estimate between-sample variance. With fewer replicates, the random effect variance estimate becomes unreliable and the model may fail to converge. [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) recommended a minimum of 3 samples per group, though power increases substantially with 5 or more.

### 3. Your Count Data Shows Overdispersion

If `check_overdispersion()` returns a dispersion ratio near 1, a simpler Poisson GLMM may suffice. But in practice, cytometry datasets [almost always exhibit overdispersion](https://doi.org/10.1038/s42003-019-0415-5) because of the biological and technical factors discussed above. The NB-GLMM is the safer default.

### 4. You Want Log2 Fold Changes for Downstream Visualization

The model's condition coefficient on the log scale converts directly to a [fold change amenable to volcano plots](https://doi.org/10.1038/s42003-019-0415-5) and standard effect size reporting conventions. The pipeline automates this conversion and feeds it into `plot_volcano()`.

---

## When to Choose Something Else

| Scenario | Recommended Method | Link |
|----------|--------------------|------|
| Cell-level binary membership (in/out of a gate), with random effects for sample | Mixed-effects logistic regression | [Mixed-Effects Logistic](mixed-effects-logistic.md) |
| Only 1 sample per condition (no replicates) | Fisher's exact test on contingency table | [Fisher's Exact Test](fisher-exact-test.md) |
| Continuous marker expression differences between conditions | Wilcoxon rank-sum or linear mixed model | [Wilcoxon Rank-Sum](wilcoxon-rank-sum.md) |
| Very sparse clusters where most samples have zero counts | Consider collapsing clusters, using zero-inflated models, or increasing `min_cluster_size` | See below |

### Sparse Clusters with Excess Zeros

When a cluster is absent from many samples, the NB-GLMM [may fail to converge](https://doi.org/10.18637/jss.v067.i01) because the likelihood surface becomes flat with respect to the dispersion parameter. The pipeline addresses this with the `min_cluster_size` parameter (default: 10 total cells across all samples), which filters out extremely rare clusters before testing. If zero-inflation is biologically meaningful (e.g., a population that is truly absent in one condition), consider:

- **Zero-inflated negative binomial (ZINB) models**: These [separately model the probability of structural zeros versus count variation](https://doi.org/10.1007/978-1-4899-3242-6), though they require more parameters and larger sample sizes.
- **Cluster merging**: Combine related rare clusters using the hierarchical structure from Leiden clustering at different resolutions.
- **Pseudocounts**: Adding a small constant (e.g., 0.5) to all counts can stabilize estimation, though this introduces [a minor bias toward the null](https://doi.org/10.1093/bioinformatics/btp616).

---

## Assumptions and Diagnostics

### 1. Conditional Negative Binomial Distribution

The counts, given the random effects and covariates, follow a [negative binomial distribution](https://doi.org/10.1007/978-1-4899-3242-6). Violations manifest as systematic patterns in residual plots. After fitting, examine Pearson or deviance residuals against fitted values. Residual plots that show a clear funnel or curvature indicate misspecification.

### 2. Log-Linear Relationship

The [log link assumes a multiplicative relationship](https://doi.org/10.1007/978-1-4899-3242-6) between condition and expected count. If the effect is additive rather than multiplicative, the log link produces biased estimates. This assumption is usually reasonable for proportional data in cytometry because biological processes (proliferation, apoptosis, migration) tend to act as [multiplicative factors on cell population sizes](https://doi.org/10.1038/s42003-019-0415-5).

### 3. Normally Distributed Random Effects

The random intercepts `b_j` are assumed to follow a [normal distribution](https://doi.org/10.18637/jss.v067.i01). Gross violations, such as multimodal random effect distributions from batch effects, can inflate type I error rates. Plot the estimated random effects (available via `ranef(model)`) as a histogram or QQ plot to check. If batch structure is evident, consider adding batch as a fixed effect or using [batch correction in the harmonization step](https://doi.org/10.12688/f1000research.11622.3) before analysis.

### 4. Independence Between Samples

The model assumes that samples are [independent conditional on the fixed and random effects](https://doi.org/10.18637/jss.v067.i01). Paired designs (e.g., pre/post treatment from the same patient) violate this unless a patient-level random effect is included. Extend the formula to `count ~ condition + offset(log_total) + (1 | patient_id/sample_id)` for nested designs.

### 5. No Systematic Compositional Bias

The offset normalizes for total cell count, implicitly assuming that [changes in one cluster do not systematically distort others](https://doi.org/10.1038/nmeth.4295). If a single dominant population changes dramatically (e.g., massive neutrophil expansion), all other clusters appear to change in the opposite direction purely as a compositional artifact. [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) discussed this limitation and suggested sensitivity analyses excluding dominant populations.

### Diagnostic Checklist

| Check | How | Concern If Failed |
|-------|-----|-------------------|
| Overdispersion | `check_overdispersion()` on a Poisson fit | Poisson is inappropriate; NB is needed |
| Convergence | Check for `fit_nb_glmm()` warnings | Estimates may be unreliable |
| Random effects normality | QQ plot of `ranef(model)$sample_id` | Batch effects or outlier samples |
| Residual patterns | Plot Pearson residuals vs. fitted | Model misspecification |
| Influential samples | Cook's distance or leave-one-out | Single sample driving the result |

---

## How to Use in the Pipeline

### From Raw FCS Files (Pipeline A)

After Leiden clustering has produced per-cell CSV files with cluster columns:

```r
source("src/analysis/differential_abundance.R")
results <- run_differential_abundance(
  data_dir = "output/04_analysis/cluster_data/",
  conditions_csv = "metadata/conditions.csv",
  output_dir = "output/04_analysis/differential_abundance/",
  condition_levels = c("Responder", "NonResponder"),
  min_cluster_size = 10,
  p_adjust_method = "BH"
)
```

### From CSV Gating Data (Pipeline C)

When cluster membership comes from gating rather than unsupervised clustering:

```r
source("src/analysis/differential_abundance.R")
results <- run_differential_abundance(
  data_dir = "data/gating_csvs/",
  conditions_csv = "metadata/conditions.csv",
  output_dir = "output/differential_abundance/",
  cluster_regex = "Clustering annotated",
  condition_levels = c("pCR", "npCR")
)
```

### Build the Abundance Matrix Separately

For custom downstream analysis or when you need the abundance matrix for multiple tests:

```r
source("src/analysis/differential_abundance.R")
abundance <- build_abundance_matrix(
  data_dir = "data/",
  conditions_csv = "metadata/conditions.csv",
  min_cluster_size = 10
)

# The abundance data frame has one row per sample:
# file_name | total_cells | cluster_1 | cluster_2 | ... | condition | log_total
head(abundance$abundance_df)

# Valid clusters that passed the minimum size filter
abundance$valid_clusters
```

### Fit a Single Cluster Manually

For debugging or custom analysis, call `fit_nb_glmm()` directly:

```r
source("src/R_shared/stats_utils.R")

count_data <- data.frame(
  count = c(120, 95, 140, 30, 45, 25),
  condition = factor(c("A", "A", "A", "B", "B", "B")),
  sample_id = paste0("S", 1:6),
  log_total = log(c(10000, 9500, 11000, 10200, 9800, 10500))
)

model <- fit_nb_glmm(
  count_data,
  formula = count ~ condition + offset(log_total) + (1 | sample_id)
)

if (!is.null(model)) {
  summary(model)
  coefs <- summary(model)$coefficients
  # Log2 fold change
  coefs["conditionB", "Estimate"] / log(2)
}
```

### Command-Line Interface

```bash
Rscript src/analysis/differential_abundance.R \
  output/04_analysis/cluster_data/ \
  metadata/conditions.csv \
  output/04_analysis/differential_abundance/
```

---

## Interpreting Results

### Output Columns

The results data frame produced by `run_differential_abundance()` contains:

| Column | Description |
|--------|-------------|
| `cluster` | Cluster name (cleaned by `clean_cluster_name()`) |
| `estimate` | Raw coefficient for the condition effect on the natural log scale |
| `std_error` | Standard error of the coefficient |
| `z_value` | Wald z-statistic (`estimate / std_error`) |
| `p_value` | Two-sided p-value from the Wald test |
| `log2fc` | Log2 fold change: `estimate / log(2)` |
| `odds_ratio` | Exponentiated estimate: `exp(estimate)`, i.e., the fold change |
| `method` | Model class used (`glmerMod` for NB, `glmerMod` for Poisson fallback) |
| `p_adjusted` | [Benjamini-Hochberg adjusted p-value](benjamini-hochberg.md) |
| `significant` | Boolean: `p_adjusted < 0.05` |

### Reading the Log2 Fold Change

- `log2fc = 0`: No difference between conditions.
- `log2fc = 1`: The cluster is 2-fold more abundant in the second condition level relative to the first (reference).
- `log2fc = -1`: The cluster is 2-fold less abundant in the second condition.
- The sign depends on the [factor level ordering](https://doi.org/10.18637/jss.v067.i01): the first level in `condition_levels` is the reference. In the default configuration `c("pCR", "npCR")`, a positive log2fc means higher abundance in npCR relative to pCR.

### The Volcano Plot

`plot_volcano()` maps `log2fc` to the x-axis and `-log10(p_adjusted)` to the y-axis. Points are colored red if `p_adjusted < 0.05`. Significant clusters are labeled using `ggrepel` when available. The dashed horizontal line marks the significance threshold. This visualization follows the conventions established in [genomic differential expression analysis](https://doi.org/10.1093/bioinformatics/btp616) and adapted for cytometry by [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5).

### Multiple Testing Correction

With many clusters tested simultaneously, the [Benjamini-Hochberg procedure](benjamini-hochberg.md) controls the false discovery rate (FDR) at the nominal level (default 5%). This is less conservative than Bonferroni correction and is the [standard choice in high-dimensional cytometry analysis](https://doi.org/10.1038/s42003-019-0415-5). The correction is applied in `run_differential_abundance()` via:

```r
results_df$p_adjusted <- p.adjust(results_df$p_value, method = "BH")
```

---

## Fallback Behavior

### When the NB-GLMM Fails

The `fit_nb_glmm()` function implements a two-tier fallback strategy:

```
NB GLMM (lme4::glmer.nb) --[fails]--> Poisson GLMM (lme4::glmer, family=poisson) --[fails]--> NULL
```

Common reasons for NB convergence failure include:

- **Very small counts**: When most samples have 0-2 cells in a cluster, the [negative binomial dispersion parameter cannot be reliably estimated](https://doi.org/10.18637/jss.v067.i01).
- **Near-zero dispersion**: If the data is genuinely Poisson (variance approximately equals mean), the NB likelihood is nearly flat in theta, causing [numerical instability in the optimization](https://doi.org/10.18637/jss.v067.i01).
- **Singular random effects**: When between-sample variance is negligible relative to within-sample variance, the random effect variance estimate hits zero, producing a [singular fit warning](https://doi.org/10.18637/jss.v067.i01).

### Implications of the Poisson Fallback

When `fit_nb_glmm()` falls back to a Poisson GLMM:

- The model still accounts for sample-level correlation via the random intercept.
- The model **does not** account for overdispersion. If the data is truly overdispersed, [p-values will be anti-conservative](https://doi.org/10.1007/978-1-4899-3242-6) (too small), increasing false positives.
- The `method` column in the results will reflect the model class, allowing you to identify which clusters used the fallback.
- Consider increasing `max_iter` (default: 20,000) if the NB model is failing due to slow convergence rather than a fundamental identifiability issue.

### When Both Models Fail

If both NB and Poisson GLMMs fail, `fit_nb_glmm()` returns `NULL`. The `run_differential_abundance()` function records this as a row with `NA` values and `method = "failed"`. Common remedies:

- Increase `min_cluster_size` to filter out problematic sparse clusters.
- Check for samples with zero total cells or degenerate conditions (only one condition level present after filtering).
- Inspect the data for extreme outlier samples that may cause [complete or quasi-complete separation](https://doi.org/10.18637/jss.v067.i01).

---

## Comparison with diffcyt

### What diffcyt Does

The [diffcyt Bioconductor package](https://doi.org/10.1038/s42003-019-0415-5) provides a comprehensive framework for differential discovery in cytometry data. For differential abundance, it uses [edgeR](https://doi.org/10.1093/bioinformatics/btp616) or [limma-voom](https://doi.org/10.1038/nmeth.4295) on cluster cell counts, wrapping the established RNA-seq differential expression machinery. The key insight from [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) was that cluster cell counts in cytometry are analogous to gene read counts in RNA-seq: both are discrete, overdispersed, and require normalization for library size (total cells or total reads).

### How This Pipeline Differs

| Aspect | diffcyt | This pipeline |
|--------|---------|---------------|
| **Dispersion estimation** | [edgeR empirical Bayes shrinkage](https://doi.org/10.1093/bioinformatics/btp616) across clusters | Per-cluster maximum likelihood via `glmer.nb()` |
| **Random effects** | Not explicitly modeled (handled by edgeR's design matrix) | Explicit `(1 \| sample_id)` random intercept via [lme4](https://doi.org/10.18637/jss.v067.i01) |
| **Normalization** | edgeR's TMM normalization or library size factors | Log-total-cells offset (equivalent to [library size normalization](https://doi.org/10.1007/978-1-4899-3242-6)) |
| **Clustering** | FlowSOM | Leiden (from `src/clustering/leiden.R`) |
| **Input format** | `SummarizedExperiment` via `prepareData()` | CSV files with binary cluster columns |
| **Fallback** | limma-voom as alternative | Poisson GLMM |

### When to Prefer diffcyt

- **Information borrowing**: edgeR's [empirical Bayes approach](https://doi.org/10.1093/bioinformatics/btp616) shrinks per-cluster dispersion estimates toward a common value. This stabilizes estimates for clusters with few samples, at the cost of assuming a shared dispersion trend. When you have many clusters and limited samples, this borrowing of information across clusters can [improve power and reduce false discovery rates](https://doi.org/10.1038/s42003-019-0415-5).
- **Established benchmarking**: diffcyt was [benchmarked extensively](https://doi.org/10.1038/s42003-019-0415-5) against simulated and real CyTOF data with known ground truth.

### When to Prefer This Pipeline's Approach

- **Explicit random effects**: The `(1 | sample_id)` term directly models sample-level heterogeneity in the generative model, which is [conceptually cleaner for designs with true biological replicates](https://doi.org/10.1038/nmeth.4295).
- **Per-cluster dispersion**: Each cluster gets its own dispersion parameter, which is more flexible when different cell populations have genuinely different levels of biological variability.
- **Integration with existing pipeline**: The implementation reads directly from the Leiden clustering output and conditions CSV, avoiding the data structure conversions required by diffcyt. The [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) CyTOF workflow demonstrated that both approaches yield comparable results on benchmark datasets.

---

## Further Reading

The following resources provide additional context for the methods described above. Each is hyperlinked inline throughout this document at the point of relevance:

- [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) --- diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. *Communications Biology*.
- [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) --- Testing for differential abundance in mass cytometry data. *Nature Methods*.
- [Bates et al. (2015)](https://doi.org/10.18637/jss.v067.i01) --- Fitting linear mixed-effects models using lme4. *Journal of Statistical Software*.
- [McCullagh and Nelder (1989)](https://doi.org/10.1007/978-1-4899-3242-6) --- Generalized Linear Models, 2nd edition. Chapman & Hall.
- [Robinson et al. (2010)](https://doi.org/10.1093/bioinformatics/btp616) --- edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*.
- [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) --- CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*.
