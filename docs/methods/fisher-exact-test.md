# Fisher's Exact Test for Cluster Enrichment

## Summary

Fisher's exact test determines whether a cell cluster is statistically over- or
under-represented in one experimental condition relative to another, using a 2x2
contingency table of cluster membership versus condition. It computes an
[exact p-value from the hypergeometric distribution](https://doi.org/10.2307/2340521)
rather than relying on asymptotic approximations, making it reliable regardless
of cell count. In this pipeline, it is the default enrichment test when data are
pooled across samples or when only a single sample per condition is available.

---

## The Method

### Mathematical Foundation

Fisher's exact test operates on a 2x2 contingency table constructed from
single-cell data:

|                | Condition A | Condition B | Row Total   |
|----------------|-------------|-------------|-------------|
| In cluster     | *a*         | *c*         | *a + c*     |
| Not in cluster | *b*         | *d*         | *b + d*     |
| Column Total   | *a + b*     | *c + d*     | *N*         |

The null hypothesis is that cluster membership is independent of condition --
that is, the proportion of cells in the cluster is the same in both conditions.

Under this null, and [conditioning on the observed marginal totals](https://doi.org/10.1002/0471249688),
the cell count *a* (cells in the cluster belonging to Condition A) follows a
**hypergeometric distribution**:

$$
P(X = a) = \frac{\binom{a+b}{a}\;\binom{c+d}{c}}{\binom{N}{a+c}}
$$

where *N = a + b + c + d* is the total number of cells. This is
[the formula originally described by Fisher (1922)](https://doi.org/10.2307/2340521)
for exact inference on contingency tables without requiring large-sample
approximations.

The two-sided p-value is computed by summing the probabilities of all tables
that are [as extreme as or more extreme than the observed table](https://doi.org/10.1002/0471249688),
given the fixed marginals. R's `fisher.test()` uses this formulation by
default, and additionally returns a conditional maximum likelihood estimate of
the odds ratio with a corresponding confidence interval based on the
[exact conditional distribution](https://doi.org/10.1002/0471249688).

### What the Pipeline Computes

The `calculate_enrichment()` function in `src/R_shared/stats_utils.R`:

1. Constructs the 2x2 table from binary cluster membership (0/1) and the two
   condition labels.
2. Calls `fisher.test()` on the matrix.
3. Returns the four cell counts (*a*, *b*, *c*, *d*), the odds ratio, the
   p-value, and the 95% confidence interval bounds.

The calling function in `src/analysis/cluster_enrichment.R` then applies
[Benjamini-Hochberg correction](benjamini-hochberg.md) across all clusters and
computes a log2 enrichment ratio from the observed proportions.

---

## When Your Data Calls For It

Fisher's exact test is appropriate when your experimental design matches all of
the following conditions.

### Exactly Two Conditions

The test is inherently pairwise. Your metadata file must assign each file (or
each cell) to exactly one of two condition labels -- for example, responders
versus non-responders, pre-treatment versus post-treatment, or healthy donors
versus patients. This two-group design is
[the classical setting for 2x2 contingency analysis](https://doi.org/10.1002/0471249688).

**Input format:** A conditions CSV with columns `file_name` and `condition`,
where `condition` takes exactly two distinct values:

```
file_name;condition
patient_01.csv;pCR
patient_02.csv;pCR
patient_03.csv;npCR
patient_04.csv;npCR
```

### Single Sample per Condition or Pooled Data

Fisher's exact test treats all cells from each condition as a single pool. It
does not model sample-to-sample variability. This makes it appropriate when you
have one sample per condition, or when you are deliberately pooling multiple
samples and accepting the assumption that
[inter-sample variation is negligible](https://doi.org/10.1038/nmeth.4295).
Lun et al. (2017) caution that ignoring sample-level variation can inflate
false-positive rates when multiple biological replicates exist; the test is
most defensible when pooling is scientifically justified.

### Binary Cluster Membership

Each cell is classified as belonging (1) or not belonging (0) to a given
cluster. The pipeline identifies cluster columns using a regex pattern (default:
`"Clustering annotated"`). Clusters below a minimum frequency threshold
(`min_cluster_frac`, default 0.001) are excluded to avoid testing near-empty
clusters where
[the hypergeometric distribution becomes degenerate](https://doi.org/10.1002/0471249688).

### Concrete Cytometry Scenarios

1. **Pilot study with one sample per arm.** You have two FCS files -- one from a
   treated patient and one from an untreated control -- both clustered by
   Leiden. Fisher's test is the natural choice because there are no replicates
   to model as random effects.

2. **Pooled comparison from a published dataset.** You downloaded a public CyTOF
   dataset and want a quick enrichment screen across dozens of clusters. Fisher's
   test provides a fast, assumption-light first pass before investing in more
   complex models.

3. **Gating-based binary classification.** Your CSV gating data (Pipeline C)
   already encodes whether each cell is positive or negative for a gate. The
   gated populations are binary by construction, directly suitable for a 2x2
   table.

---

## When to Choose Something Else

Fisher's exact test is not always the right tool. The table below maps common
situations to better alternatives.

| Situation | Problem with Fisher's | Better Method |
|-----------|----------------------|---------------|
| Multiple biological replicates per condition | Ignores inter-sample variability, [inflating false positives](https://doi.org/10.1038/nmeth.4295) | [Mixed-Effects Logistic Regression](mixed-effects-logistic.md) |
| Sample-level count data (cells per cluster per sample) | Loses sample-level information by pooling | [Negative Binomial GLMM](negative-binomial-glmm.md) |
| Continuous marker expression differences | Cannot test continuous outcomes | [Wilcoxon Rank-Sum](wilcoxon-rank-sum.md) |
| Testing many clusters simultaneously | Raw p-values will contain false positives | Apply [Benjamini-Hochberg](benjamini-hochberg.md) correction (the pipeline does this automatically) |
| More than two conditions | Fisher's test is strictly pairwise | Use pairwise Fisher with multiplicity correction, or switch to a GLMM |
| Confounders (batch, patient covariates) | No mechanism to adjust for covariates | [Mixed-Effects Logistic Regression](mixed-effects-logistic.md) with additional fixed effects |

As a rule of thumb: if you have three or more biological replicates per
condition, the [diffcyt framework](https://doi.org/10.1038/s42003-019-0415-5)
and [Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) both recommend
mixed-effects or generalized linear mixed models over pooled tests like Fisher's.

---

## Assumptions and Limitations

### 1. Independence of Observations

Fisher's exact test assumes that each cell is an [independent observation](https://doi.org/10.1002/0471249688).
In cytometry, cells from the same sample share biological and technical
variation, violating this assumption. When multiple samples are pooled per
condition, this pseudo-replication
[can dramatically underestimate p-values](https://doi.org/10.1038/nmeth.4295),
because the effective sample size is the number of biological replicates, not
the number of cells. This is the most important limitation to understand:
statistical significance from Fisher's test on pooled multi-sample data should
be interpreted as descriptive rather than inferential.

### 2. Fixed Marginals

The exact probability calculation conditions on [both sets of marginal totals being fixed](https://doi.org/10.2307/2340521).
In practice, only the column margins (total cells per condition) are fixed by
the experimental design; the row margins (total cells in/out of cluster) are
random. This is a
[well-known conservative feature](https://doi.org/10.1002/0471249688) of
Fisher's test -- it can be slightly conservative in settings where only one
margin is truly fixed, though the impact is usually small.

### 3. Correct Classification

The test assumes that cluster assignments are accurate. Misclassification of
cells between clusters adds noise and
[biases the odds ratio toward the null](https://doi.org/10.1002/0471249688),
reducing statistical power. In unsupervised clustering workflows,
[cluster quality should be validated](https://doi.org/10.12688/f1000research.11622.3)
using marker expression heatmaps before running enrichment tests.

### 4. No Confounders

Fisher's test has no mechanism to adjust for confounding variables such as
patient age, batch, or instrument settings. If confounders are
[correlated with both condition and cluster membership](https://doi.org/10.1038/s42003-019-0415-5),
the test can produce spurious associations. The diffcyt framework by Weber et
al. (2019) recommends regression-based approaches that can include covariates.

### 5. Sensitivity to Large Cell Counts

With tens of thousands of cells per condition, Fisher's test has very high
statistical power and will detect
[even biologically trivial differences in proportions](https://doi.org/10.1186/s12859-015-0751-5).
A cluster that is 5.1% of Condition A and 4.9% of Condition B might reach
statistical significance at p < 0.05 with sufficient cell counts, despite
the difference being meaningless. Always examine effect sizes (odds ratios and
log2 enrichment) alongside p-values.

---

## How to Use in the Pipeline

### Pipeline A: From Raw FCS Files

After running the full Pipeline A path (loading, preprocessing,
harmonization, clustering), the cluster assignment CSVs will be in your output
directory. Run enrichment analysis using Fisher's test by setting
`stat_method = "pooled"`:

```r
# Source the enrichment analysis script
source("src/analysis/cluster_enrichment.R")

# Run Fisher's exact test on all clusters
results <- run_cluster_enrichment(
  data_dir        = "output/cluster_data/",
  conditions_csv  = "metadata/conditions.csv",
  output_dir      = "output/enrichment/",
  stat_method     = "pooled",         # Fisher only (no mixed effects)
  condition_levels = c("pCR", "npCR"),
  min_cluster_frac = 0.001,           # Exclude clusters below 0.1%
  p_adjust_method  = "BH"             # Benjamini-Hochberg correction
)

# Inspect significant clusters
sig <- results[results$significant == TRUE, ]
print(sig[, c("cluster", "odds_ratio", "p_value", "p_adjusted", "log2_enrichment")])
```

### Pipeline C: From CSV Gating Data

If your input is pre-gated CSV data with binary cluster columns, the same
function works directly:

```r
# Source the enrichment analysis script
source("src/analysis/cluster_enrichment.R")

# Run enrichment on gating data
results <- run_cluster_enrichment(
  data_dir        = "data/gating_csvs/",
  conditions_csv  = "metadata/conditions.csv",
  output_dir      = "output/enrichment_gating/",
  cluster_regex   = "^gate_",         # Adjust regex to match your gate columns
  stat_method     = "pooled",
  condition_levels = c("responder", "non_responder"),
  p_adjust_method  = "BH"
)
```

### Calling the Low-Level Function Directly

For fine-grained control, you can call `calculate_enrichment()` from
`src/R_shared/stats_utils.R` on a single cluster:

```r
source("src/R_shared/stats_utils.R")

# Prepare data: needs a binary cluster column and a 'condition' column
result <- calculate_enrichment(
  cluster_col = "Clustering annotated CD4+ T cells",
  cells_data  = my_cells_df,
  cond_levels = c("pCR", "npCR")
)

# result is a tibble with columns:
#   in_A, not_in_A, in_B, not_in_B, odds_ratio, p_value, ci_lower, ci_upper, method
print(result)
```

### Output Format

The pipeline writes `enrichment_results.csv` (semicolon-delimited) to the
output directory with the following columns:

| Column | Description |
|--------|-------------|
| `cluster` | Cleaned cluster name |
| `in_A` | Cells in cluster, Condition A |
| `not_in_A` | Cells not in cluster, Condition A |
| `in_B` | Cells in cluster, Condition B |
| `not_in_B` | Cells not in cluster, Condition B |
| `odds_ratio` | Conditional maximum likelihood estimate |
| `p_value` | Two-sided exact p-value |
| `ci_lower` | Lower bound of 95% confidence interval for the odds ratio |
| `ci_upper` | Upper bound of 95% confidence interval for the odds ratio |
| `method` | `"fisher"`, `"fisher_failed"`, or `"fisher_error"` |
| `p_adjusted` | BH-corrected p-value |
| `significant` | `TRUE` if `p_adjusted < 0.05` |
| `log2_enrichment` | log2 of the proportion ratio (Condition A / Condition B) |

---

## Interpreting Results

### Odds Ratio

The odds ratio (OR) quantifies the strength of association between cluster
membership and condition. It is defined as:

$$
\text{OR} = \frac{a \cdot d}{b \cdot c}
$$

where *a*, *b*, *c*, *d* are the cells of the 2x2 table. An OR of 1 indicates
no association. An OR > 1 means the cluster is enriched in Condition A; an OR
< 1 means enrichment in Condition B. The pipeline reports the
[conditional maximum likelihood estimate](https://doi.org/10.1002/0471249688),
which is the standard estimator produced by R's `fisher.test()` and differs
slightly from the sample OR above in small-sample situations.

### Confidence Intervals

The 95% confidence interval for the odds ratio is computed from the
[exact conditional distribution](https://doi.org/10.1002/0471249688), not from
a normal approximation. If the interval excludes 1, the enrichment is
statistically significant at the 5% level. Wide intervals indicate low
precision, typically due to small cell counts in one or more cells of the
contingency table. The [Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
recommends reporting confidence intervals alongside p-values to convey the
uncertainty in effect size estimates.

### P-Values

The p-value is the probability of observing a table as extreme as or more
extreme than the one obtained, given the null hypothesis of independence and
the fixed marginals. The pipeline applies
[Benjamini-Hochberg correction](benjamini-hochberg.md) across all tested
clusters, producing an adjusted p-value (`p_adjusted`) that controls the
[false discovery rate](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) at
5%. Clusters with `p_adjusted < 0.05` are flagged as significant.

### Log2 Enrichment

The pipeline computes a log2 enrichment ratio from the observed proportions:

$$
\log_2\!\left(\frac{a / (a + b)}{c / (c + d)}\right)
$$

This quantity is more interpretable than the odds ratio for cytometry
applications: a value of +1 means the proportion in Condition A is double that
in Condition B. Unlike the odds ratio, it is undefined when a proportion is
zero; the pipeline returns `NA` in that case.

### Practical Significance Checklist

When reviewing enrichment results, consider each of the following:

1. **Is the adjusted p-value significant?** Check `p_adjusted`, not `p_value`.
2. **Is the odds ratio meaningfully different from 1?** A statistically
   significant OR of 1.02 is unlikely to be biologically relevant.
3. **Is the confidence interval narrow enough to be useful?** A CI of
   [0.5, 200] provides little practical guidance.
4. **What is the absolute proportion difference?** A cluster that is 0.01% vs
   0.02% of total cells, while nominally doubled, may not be actionable.
5. **Was the cluster validated biologically?** Check marker expression heatmaps
   to confirm the cluster represents a coherent cell population, as recommended
   by the [Nowicka et al. (2017) workflow](https://doi.org/10.12688/f1000research.11622.3).

---

## Comparison with Chi-Square

The Pearson chi-square test is the other classical approach for testing
independence in contingency tables. It approximates the distribution of the
test statistic using a [chi-square distribution with one degree of freedom](https://doi.org/10.1002/0471249688).
Fisher's exact test is preferred in this pipeline for several reasons.

### Small and Unbalanced Tables

The chi-square approximation is
[unreliable when expected cell counts are below 5](https://doi.org/10.1002/sim.2832),
a situation that arises frequently in cytometry when testing rare cell
populations. Campbell (2007) demonstrated through extensive simulations that the
chi-square test can produce substantially incorrect p-values for sparse tables,
while Fisher's exact test remains valid at any sample size. In flow cytometry,
cluster sizes can range from a handful of cells to hundreds of thousands,
meaning some clusters will inevitably fall into the sparse regime where the
chi-square approximation breaks down.

### Exactness vs. Approximation

Fisher's test computes the [exact probability under the null hypothesis](https://doi.org/10.2307/2340521)
by enumerating all possible tables with the same marginals. There is no
approximation involved, so the p-value is valid for any combination of cell
counts. The chi-square test, by contrast, relies on a
[large-sample normal approximation to the multinomial distribution](https://doi.org/10.1002/0471249688)
that can be inaccurate even at moderate sample sizes when the table is
imbalanced.

### Computational Considerations

For 2x2 tables, Fisher's exact test is computationally trivial -- R's
`fisher.test()` runs in microseconds. The computational advantage of the
chi-square approximation is irrelevant for 2x2 tables and only becomes
meaningful for larger tables. Since the enrichment analysis in this pipeline is
strictly 2x2 (one cluster, two conditions), there is no reason to accept the
approximation error of chi-square.

### When Chi-Square Might Be Acceptable

For very large, balanced tables (all expected counts >> 5), the chi-square test
produces p-values [nearly identical to Fisher's](https://doi.org/10.1002/sim.2832).
In such cases the choice is immaterial. However, since Fisher's test is exact
and equally fast for 2x2 tables, there is no practical benefit to switching.

### Summary Table

| Property | Fisher's Exact | Chi-Square |
|----------|---------------|------------|
| P-value validity | [Exact for all sample sizes](https://doi.org/10.2307/2340521) | [Approximate; unreliable when expected counts < 5](https://doi.org/10.1002/sim.2832) |
| Computation (2x2) | Microseconds | Microseconds |
| Assumptions | Fixed marginals | Large-sample approximation |
| Rare clusters | Fully valid | [May produce incorrect p-values](https://doi.org/10.1002/sim.2832) |
| Pipeline default | Yes | No |

---

## Effect Size Considerations

### The Problem of Overpowered Tests

Flow cytometry routinely generates datasets with tens of thousands to millions
of cells. At these scales, Fisher's exact test (and indeed any hypothesis test)
becomes
[extremely sensitive to trivially small differences](https://doi.org/10.1186/s12859-015-0751-5).
Lin et al. (2015) specifically noted that in high-throughput biological assays,
statistical significance and biological significance can diverge dramatically
when cell counts are large.

Consider a concrete example: if Cluster X contains 5,100 of 100,000 cells in
Condition A (5.10%) and 4,900 of 100,000 cells in Condition B (4.90%), Fisher's
exact test yields p = 0.024, OR = 1.04. The 0.2 percentage point difference is
statistically significant but almost certainly not biologically meaningful.

### Recommended Effect Size Metrics

To guard against over-interpreting significant p-values, examine the following
alongside the p-value:

1. **Odds ratio magnitude.** The
   [conditional maximum likelihood odds ratio](https://doi.org/10.1002/0471249688)
   is the natural effect size for a 2x2 table. Rules of thumb from
   epidemiological literature suggest OR < 1.5 (or > 0.67) often reflects a
   weak association.

2. **Log2 enrichment.** The pipeline computes this from raw proportions. In
   cytometry, a |log2 enrichment| of at least 0.5 (approximately 40% change
   in proportion) is a reasonable threshold for further investigation, though
   the appropriate threshold depends on the biological context.

3. **Absolute proportion difference.** Report the actual percentages. A cluster
   shifting from 0.01% to 0.02% (log2 enrichment = 1) is harder to validate
   experimentally than one shifting from 5% to 10%.

4. **Confidence interval width.** Narrow intervals around an OR near 1 confirm
   the absence of a meaningful effect. Wide intervals indicate insufficient
   data to draw conclusions.

### Cytometry-Specific Guidance

The [diffcyt framework](https://doi.org/10.1038/s42003-019-0415-5) addresses
the effect-size problem by using sample-level aggregate statistics rather than
cell-level tests, which naturally reduces the effective sample size and avoids
the overpowered-test issue. When using Fisher's exact test on pooled data, the
[Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
recommends always visualizing enrichment results with volcano-style plots
(log2 enrichment versus -log10 adjusted p-value) so that both statistical and
practical significance are visible at a glance.

[Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) further recommend that
when the goal is formal inference about differential abundance, sample-level
modeling (such as the negative binomial GLMM or edgeR-based approaches) should
replace cell-level pooled tests, precisely because the effective sample size is
then determined by the number of biological replicates rather than the number of
cells.

### When Large Effects Are Expected

In some cytometry applications, large effect sizes are genuinely expected --
for instance, a drug that ablates a specific T cell subset, or a disease state
that dramatically expands a myeloid population. In these scenarios, Fisher's
exact test is well-suited: the effect is real, large, and the test will detect
it with high confidence. The concern about trivial significance applies
primarily to subtle phenotypic shifts across many clusters in exploratory
screens.

---

## Implementation Details

### Source Code Location

- **Core function:** [`calculate_enrichment()`](/src/R_shared/stats_utils.R) in
  `src/R_shared/stats_utils.R` (lines 7-49)
- **Pipeline orchestrator:** [`run_cluster_enrichment()`](/src/analysis/cluster_enrichment.R)
  in `src/analysis/cluster_enrichment.R`

### Error Handling

The implementation wraps `fisher.test()` in a `tryCatch` block. Three outcomes
are possible:

| `method` value | Meaning |
|----------------|---------|
| `"fisher"` | Test ran successfully |
| `"fisher_failed"` | Input validation failed (negative counts or empty condition group) |
| `"fisher_error"` | R's `fisher.test()` threw an error (e.g., non-finite values) |

When the test fails, all numeric output columns (`odds_ratio`, `p_value`,
`ci_lower`, `ci_upper`) are set to `NA`. Downstream functions that apply
BH correction handle `NA` p-values gracefully via `p.adjust()`.

### Resume Functionality

For large datasets with many clusters, the pipeline saves intermediate results
every 10 clusters to `intermediate_results.rds`. If the run is interrupted, set
`resume = TRUE` to continue from where it left off. The intermediate file is
deleted upon successful completion.

---

## References by Topic

All references are linked inline throughout this document. The key publications
are:

- **Fisher's exact test (original):** [Fisher, R. A. (1922)](https://doi.org/10.2307/2340521). On the interpretation of chi-square from contingency tables, and the calculation of P. *J. R. Stat. Soc.*, 85(1), 87-94.
- **Categorical data analysis (comprehensive):** [Agresti, A. (2002)](https://doi.org/10.1002/0471249688). *Categorical Data Analysis*, 2nd ed. Wiley.
- **Chi-square vs. Fisher's for small samples:** [Campbell, I. (2007)](https://doi.org/10.1002/sim.2832). Chi-squared and Fisher-Irwin tests of two-by-two tables with small sample recommendations. *Stat. Med.*, 26(19), 3661-3675.
- **Differential abundance in mass cytometry:** [Lun, A. T. L., Richard, A. C., & Marioni, J. C. (2017)](https://doi.org/10.1038/nmeth.4295). Testing for differential abundance in mass cytometry data. *Nat. Methods*, 14(7), 707-709.
- **diffcyt framework:** [Weber, L. M., Nowicka, M., Soneson, C., & Robinson, M. D. (2019)](https://doi.org/10.1038/s42003-019-0415-5). diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. *Commun. Biol.*, 2, 183.
- **CyTOF analysis workflow:** [Nowicka, M. et al. (2017)](https://doi.org/10.12688/f1000research.11622.3). CyTOF workflow: Differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*, 6, 748.
- **Statistical power in high-throughput assays:** [Lin, Y. et al. (2015)](https://doi.org/10.1186/s12859-015-0751-5). Sample size and power analysis for high-throughput biological experiments. *BMC Bioinformatics*, 16, 304.
- **FDR correction:** [Benjamini, Y. & Hochberg, Y. (1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). Controlling the false discovery rate. *J. R. Stat. Soc. B*, 57(1), 289-300.
