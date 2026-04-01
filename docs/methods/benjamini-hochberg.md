# Benjamini-Hochberg False Discovery Rate Correction

## Summary

The Benjamini-Hochberg (BH) procedure is a step-up method that [controls the expected proportion of falsely rejected hypotheses](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) --- the false discovery rate (FDR) --- rather than the probability of making any single false rejection. In this pipeline, it is the default multiple testing correction (`p_adjust_method = "BH"`) applied after simultaneously testing 20--50 immune cell clusters for enrichment or differential abundance. By tolerating a small, controlled fraction of false positives instead of guarding against every possible error, [FDR-based correction preserves substantially more statistical power](https://doi.org/10.1111/1467-9868.00346) than classical family-wise error rate methods, making it the standard approach in high-dimensional biological experiments.

---

## The Method

### Mathematical Definition of the False Discovery Rate

Consider *m* hypothesis tests conducted simultaneously. Let *V* denote the number of false positives (null hypotheses that were incorrectly rejected) and *R* denote the total number of rejections. The false discovery rate is [defined as the expected proportion of false discoveries among all rejections](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):

```
FDR = E[V / R]
```

where V/R is defined as 0 when R = 0 (i.e., when no hypotheses are rejected, the false discovery proportion is zero by convention).

This quantity differs fundamentally from the family-wise error rate (FWER), which is P(V >= 1) --- the probability of making *any* false rejection at all.

### The Step-Up Algorithm

The BH procedure, [introduced by Benjamini and Hochberg in 1995](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x), works as follows:

1. **Order the p-values.** Given *m* raw p-values from *m* simultaneous tests, sort them in ascending order: p_(1) <= p_(2) <= ... <= p_(m).

2. **Find the largest significant index.** Choose a target FDR level *q* (commonly 0.05). Find the largest integer *k* such that:

   ```
   p_(k) <= (k / m) * q
   ```

3. **Reject hypotheses.** Reject all hypotheses H_(1), H_(2), ..., H_(k) corresponding to p_(1), p_(2), ..., p_(k).

This is called a **step-up** procedure because it starts from the smallest p-value and steps upward, in contrast to step-down methods like [Holm's procedure](https://doi.org/10.2307/4615733) that start from the largest.

### Adjusted P-Values

In practice, rather than applying the algorithm to find a cutoff *k*, R's `p.adjust()` function computes adjusted p-values directly. [The BH-adjusted p-value for the *i*-th ordered test](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) is:

```
p_adjusted_(i) = min over j >= i of { m * p_(j) / j }
```

enforced to be monotonically non-decreasing. A hypothesis is rejected at FDR level *q* if and only if its adjusted p-value is less than or equal to *q*. This formulation is exactly what `p.adjust(p, method = "BH")` returns in R.

---

## FDR vs FWER

### Family-Wise Error Rate (FWER)

Classical methods such as [Bonferroni](https://doi.org/10.2307/2283989) and [Holm](https://doi.org/10.2307/4615733) control the probability that *any* null hypothesis is incorrectly rejected: FWER = P(V >= 1). This is appropriate when even a single false positive would be catastrophic (e.g., clinical trial endpoints), but it comes at a steep cost in power. With Bonferroni, the significance threshold for each test becomes alpha/m; testing 30 clusters at alpha = 0.05 gives a per-test threshold of 0.0017.

### False Discovery Rate (FDR)

The FDR instead [controls the expected *fraction* of false discoveries](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). At FDR = 0.05, you accept that on average 5% of the discoveries you declare significant may be false positives. The remaining 95% are expected to be true findings.

### Why FDR Became the Standard in High-Dimensional Biology

In genomics, proteomics, and cytometry, researchers routinely test hundreds to thousands of features simultaneously. [Storey (2002) demonstrated](https://doi.org/10.1111/1467-9868.00346) that FDR is a more natural and powerful error rate for large-scale testing: as the number of tests grows, the FDR framework scales gracefully, while FWER control becomes increasingly conservative, potentially masking genuine biological signals. [Storey and Tibshirani (2003) formalized the q-value](https://doi.org/10.1073/pnas.1530509100), the FDR analogue of the p-value, and showed that exploiting the estimated proportion of true nulls can yield even more discoveries while maintaining the same FDR guarantee. These advances established FDR as the dominant multiple testing framework across computational biology, including [CyTOF and flow cytometry differential analysis pipelines](https://doi.org/10.12688/f1000research.11622.3).

| Property | FWER (e.g., Bonferroni) | FDR (Benjamini-Hochberg) |
|----------|------------------------|--------------------------|
| Controls | P(any false positive) | E(fraction of false positives) |
| Power | Low when *m* is large | Substantially higher |
| Typical use | Few critical tests | Many exploratory tests |
| Cytometry fit | Rarely appropriate | Default choice |

---

## Why It Matters for Cytometry

### The Multiple Testing Burden

A typical flow or mass cytometry experiment in this pipeline identifies 20--50 immune cell clusters via Leiden clustering. Each cluster is tested for differential abundance between conditions (e.g., responders vs non-responders). Without any correction, the expected number of false positives at alpha = 0.05 with 30 tests is 30 * 0.05 = 1.5. With 50 clusters, you expect 2.5 spurious discoveries by chance alone.

In [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3), a benchmark study of cytometry differential analysis workflows, the authors applied BH correction after testing multiple cell populations and demonstrated that uncorrected p-values led to inflated false positive rates in simulation studies, while BH adjustment maintained FDR control without excessive loss of power.

### Concrete Impact

Consider a realistic scenario: you test 40 clusters for enrichment in a treatment response study. Without correction, 2 of your 5 "significant" findings may be artifacts. After BH correction at FDR = 0.05, the 3 clusters that remain significant have a collective false discovery rate of at most 5%. This means you can confidently pursue validation experiments on those clusters.

---

## When Your Data Calls For It

BH correction should be applied whenever you perform **more than one** hypothesis test on the same dataset. In this pipeline, that includes:

### Multiple Cluster Tests

The primary use case. When `run_cluster_enrichment()` or `run_differential_abundance()` loops over all valid clusters and produces a vector of p-values, [BH correction prevents the accumulation of false positives across clusters](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):

- **Cluster enrichment**: Fisher's exact test or mixed-effects logistic regression per cluster, testing whether each cluster is over-represented in one condition versus another.
- **Differential abundance**: Negative binomial GLMM per cluster, testing whether the cell count in each cluster differs between conditions after accounting for sample-level variability.

### Multiple Marker Tests

When `run_marker_statistics()` in `src/analysis/marker_statistics.R` tests multiple surface markers for differential expression, each marker constitutes a separate hypothesis test. BH correction across all markers prevents inflated false positive rates.

### Any Multi-Hypothesis Scenario

More broadly, any analysis that produces a table of p-values from related tests on the same data warrants multiple testing correction:

- Testing multiple gating strategies against a reference
- Comparing cluster frequencies across more than two conditions (pairwise tests)
- Screening marker panels for discriminatory markers

---

## When to Choose Something Else

BH is the default for good reason, but certain situations call for alternatives:

### Strict FWER Control: Bonferroni or Holm

When a single false positive is unacceptable --- for instance, confirming a clinical biomarker before a regulatory submission --- use [Bonferroni correction](https://doi.org/10.2307/2283989), which divides alpha by the number of tests, or [Holm's step-down procedure](https://doi.org/10.2307/4615733), which is uniformly more powerful than Bonferroni while still controlling the FWER:

```r
results <- run_cluster_enrichment(..., p_adjust_method = "bonferroni")
results <- run_cluster_enrichment(..., p_adjust_method = "holm")
```

### Arbitrary Dependence: Benjamini-Yekutieli

The standard BH procedure [assumes independence or positive regression dependency among the test statistics](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). When the tests exhibit arbitrary dependence structures --- for example, when clusters share many of the same cells or markers are highly correlated --- [Benjamini and Yekutieli (2001)](https://doi.org/10.1214/aos/1013699998) proposed a modified procedure that controls FDR under any dependency by replacing the denominator *m* with *m * sum(1/i for i in 1..m)*. This is more conservative but provides a guaranteed FDR bound regardless of the correlation structure:

```r
results <- run_cluster_enrichment(..., p_adjust_method = "BY")
```

### Adaptive Procedures: Storey's Q-Value

[Storey (2002)](https://doi.org/10.1111/1467-9868.00346) and [Storey and Tibshirani (2003)](https://doi.org/10.1073/pnas.1530509100) introduced an adaptive procedure that first estimates the proportion of true null hypotheses (pi_0) and then uses this estimate to produce tighter FDR bounds. When many hypotheses are truly non-null (common in large-scale studies), this yields more rejections at the same FDR level. The R package `qvalue` implements this method. Use it when you have a large number of tests (hundreds or more) and suspect a substantial fraction are truly non-null.

### Single Test: No Correction Needed

If you are testing exactly one hypothesis --- for example, a single pre-specified cluster of interest --- no multiple testing correction is needed. The raw p-value is valid as-is.

---

## Assumptions

The BH procedure's FDR guarantee rests on specific assumptions about the input p-values and their dependency structure.

### Valid P-Values Under the Null

Each p-value must be [uniformly distributed on [0, 1] under its respective null hypothesis](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). In practice, this means the underlying test (Fisher's exact, GLMM likelihood ratio, Wilcoxon rank-sum) must be correctly specified. Violations occur when:

- The distributional assumptions of the test are wrong (e.g., using a parametric test on heavily skewed data without transformation).
- The test is applied to data that violate independence assumptions (e.g., pseudoreplication where cells from the same sample are treated as independent observations). This pipeline's use of [mixed-effects models with sample-level random effects](https://doi.org/10.12688/f1000research.11622.3) addresses this by properly accounting for within-sample correlation.

### Positive Regression Dependency on a Subset (PRDS)

[Benjamini and Yekutieli (2001)](https://doi.org/10.1214/aos/1013699998) proved that the BH procedure controls FDR at exactly level *q* under independence, and at a level at most *q* under a condition they called **positive regression dependency on a subset (PRDS)**. Informally, PRDS holds when the test statistics are positively correlated or independent --- seeing a large value for one statistic makes it more (or equally) likely to see a large value for another. This condition is [commonly satisfied in biological experiments](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) where features tend to be positively correlated (e.g., co-expressed markers, overlapping cell populations).

In flow cytometry, PRDS is a reasonable assumption because:

- Cell clusters identified by Leiden clustering may overlap in marker space, inducing positive correlation among test statistics.
- Samples that drive enrichment in one cluster often drive enrichment in related clusters.

When you suspect strong negative dependencies (rare in cytometry), use the [Benjamini-Yekutieli procedure](https://doi.org/10.1214/aos/1013699998) instead (`p_adjust_method = "BY"`).

### Number of Tests

The BH procedure is [valid for any number of tests *m >= 1*](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x), though its power advantage over FWER methods grows with *m*. For very small *m* (e.g., 2--3 tests), the practical difference between BH and Bonferroni is minimal.

---

## How to Use in the Pipeline

BH correction is applied automatically in both main analysis functions. You do not need to call `p.adjust()` yourself --- the pipeline handles it.

### Default Usage (BH Correction)

```r
# BH is the default in cluster enrichment
results <- run_cluster_enrichment(data_dir = "output/cluster_data/",
                                   conditions_csv = "metadata/conditions.csv",
                                   p_adjust_method = "BH")  # default

# Also in differential abundance
results <- run_differential_abundance(data_dir = "output/cluster_data/",
                                       conditions_csv = "metadata/conditions.csv",
                                       p_adjust_method = "BH")  # default
```

### Switching to a Different Correction

```r
# Strict FWER control
results <- run_cluster_enrichment(..., p_adjust_method = "bonferroni")
results <- run_cluster_enrichment(..., p_adjust_method = "holm")

# FDR under arbitrary dependence
results <- run_cluster_enrichment(..., p_adjust_method = "BY")
```

### Where It Happens in the Code

In `src/analysis/cluster_enrichment.R` (line 133):

```r
results_df$p_adjusted <- p.adjust(results_df$p_value, method = p_adjust_method)
results_df$significant <- results_df$p_adjusted < 0.05
```

In `src/analysis/differential_abundance.R` (line 258):

```r
results_df$p_adjusted <- p.adjust(results_df$p_value, method = p_adjust_method)
results_df$significant <- results_df$p_adjusted < 0.05
```

Both functions collect all per-cluster raw p-values into a single vector, apply `p.adjust()` once across all tests, and flag results below 0.05 as significant.

### Available Methods in R's `p.adjust()`

The `p_adjust_method` parameter accepts [any method supported by R's `p.adjust()` function](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x):

| Method string | Full name | Controls |
|---------------|-----------|----------|
| `"BH"` or `"fdr"` | Benjamini-Hochberg | FDR |
| `"BY"` | Benjamini-Yekutieli | FDR (arbitrary dependence) |
| `"bonferroni"` | Bonferroni | FWER |
| `"holm"` | Holm | FWER |
| `"hochberg"` | Hochberg | FWER |
| `"hommel"` | Hommel | FWER |
| `"none"` | No correction | Per-test alpha |

---

## Interpreting Adjusted P-Values

### What an Adjusted P-Value Means

A BH-adjusted p-value (sometimes called a q-value, though [Storey's q-value](https://doi.org/10.1073/pnas.1530509100) is a distinct, related quantity) represents the minimum FDR at which that hypothesis would be rejected. Concretely:

- If a cluster has an adjusted p-value of 0.03, then including this cluster among the rejected set would keep the overall FDR at or below 3%.
- If a cluster has an adjusted p-value of 0.12, you would need to accept a 12% false discovery rate to call it significant.

### Common Thresholds

| Threshold | Interpretation | Typical use |
|-----------|---------------|-------------|
| q < 0.05 | At most 5% of discoveries are expected to be false | Standard threshold for publication-quality results in [cytometry differential analysis](https://doi.org/10.12688/f1000research.11622.3) |
| q < 0.10 | At most 10% of discoveries are expected to be false | Used in exploratory analyses or when sample sizes are small; [common in genomics discovery studies](https://doi.org/10.1073/pnas.1530509100) |
| q < 0.20 | At most 20% of discoveries are expected to be false | Aggressive threshold for hypothesis generation only |

This pipeline uses **q < 0.05** as the default significance cutoff (the `results_df$significant <- results_df$p_adjusted < 0.05` line in both analysis functions).

### Adjusted P-Values Are Not Per-Test Error Rates

A common misinterpretation: an adjusted p-value of 0.03 does **not** mean there is a 3% chance that this specific cluster is a false positive. Rather, [it means that among all clusters with adjusted p-values at or below 0.03, the expected fraction of false positives is at most 3%](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). The guarantee is collective, not individual.

### Reading the Output

The pipeline produces CSV output (semicolon-delimited) with these key columns:

| Column | Meaning |
|--------|---------|
| `p_value` | Raw p-value from the per-cluster test (Fisher's exact or GLMM) |
| `p_adjusted` | BH-adjusted p-value (FDR-controlled) |
| `significant` | Boolean: `TRUE` if `p_adjusted < 0.05` |
| `log2fc` or `log2_enrichment` | Effect size (log2 fold-change or enrichment ratio) |

Always report both the adjusted p-value *and* the effect size. A statistically significant cluster with a tiny fold-change may not be biologically meaningful.

---

## Relationship to Other Methods in This Pipeline

The BH procedure is not a statistical test itself --- it is a **post-hoc correction** applied to the p-values produced by the actual tests. In this pipeline, BH adjustment is applied after:

- **[Fisher's Exact Test](fisher-exact-test.md)** --- Used in `run_cluster_enrichment()` to test whether each cluster is enriched in one condition versus another based on a 2x2 contingency table of cell counts. BH correction is applied across all cluster-level Fisher p-values.

- **[Mixed-Effects Logistic Regression](mixed-effects-logistic.md)** --- Also used in `run_cluster_enrichment()` when `stat_method = "mixed_effects"`. The logistic GLMM accounts for within-sample correlation by including a random intercept per sample. BH correction is applied to the condition-effect p-values across all clusters.

- **[Negative Binomial GLMM](negative-binomial-glmm.md)** --- Used in `run_differential_abundance()` to model cluster cell counts with overdispersion. Each cluster produces a p-value for the condition effect; BH correction is applied across all clusters.

- **[Wilcoxon Rank-Sum Test](wilcoxon-rank-sum.md)** --- Used in `run_marker_statistics()` for non-parametric per-marker comparisons. When testing multiple markers, BH correction prevents the accumulation of false positives across the marker panel.

In every case, the workflow is the same: run the individual test per cluster (or per marker), collect all raw p-values into a vector, and apply `p.adjust(p_values, method = "BH")` once to obtain adjusted p-values.
