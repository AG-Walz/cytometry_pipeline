# Wilcoxon Rank-Sum Test (Mann-Whitney U) for Marker Expression Analysis

## Summary

The Wilcoxon rank-sum test is a non-parametric method for comparing continuous marker expression intensities between two independent groups, applied per-marker and per-population in flow cytometry analysis. It [ranks all observations from both groups jointly and evaluates whether one group's ranks are systematically higher](https://doi.org/10.2307/3001968), making it robust to the skewed, heavy-tailed distributions that are characteristic of cytometry data. In this pipeline it is invoked through `run_marker_statistics()` in `src/analysis/marker_statistics.R` with `test_type = "wilcox"`.

---

## The Method

### Origins

The test was [introduced by Frank Wilcoxon in 1945](https://doi.org/10.2307/3001968) as a rank-based alternative to the two-sample t-test. [Mann and Whitney (1947) independently derived the same test in a more general framework](https://doi.org/10.1214/aoms/1177730491), proving its theoretical properties and establishing the U statistic formulation. The two formulations are [mathematically equivalent and interchangeable](https://doi.org/10.1214/aoms/1177730491); the name "Wilcoxon rank-sum" and "Mann-Whitney U" refer to the same procedure.

### Ranking procedure

Given two independent samples of sizes n1 and n2:

1. **Pool and rank.** Combine all N = n1 + n2 observations into a single sequence. Assign ranks 1 through N from smallest to largest. Tied observations receive the [average of the ranks they span](https://doi.org/10.2307/3001968).

2. **Sum ranks.** Compute R1, the sum of ranks assigned to group 1. Under the null hypothesis that both groups are drawn from identical distributions, the [expected value of R1 is n1(N + 1) / 2](https://doi.org/10.1214/aoms/1177730491).

3. **Compute the U statistic.** The Mann-Whitney U is defined as:

   ```
   U = R1 - n1(n1 + 1) / 2
   ```

   Equivalently, U counts the number of pairwise comparisons in which a value from group 1 exceeds a value from group 2. [Under the null hypothesis, U has expected value n1 * n2 / 2](https://doi.org/10.1214/aoms/1177730491).

4. **Determine significance.** For small samples, exact p-values are computed by enumerating all possible rank assignments. R's `wilcox.test()` uses this [exact permutation distribution by default when sample sizes are small and there are no ties](https://doi.org/10.2307/3001968).

### Normal approximation for large samples

When both n1 and n2 are sufficiently large (conventionally n1, n2 > 20), the [distribution of U is approximately normal](https://doi.org/10.1214/aoms/1177730491) with:

```
mean(U) = n1 * n2 / 2
var(U)  = n1 * n2 * (N + 1) / 12
```

A [continuity-corrected z-statistic](https://doi.org/10.1214/aoms/1177730491) is then:

```
z = (U - mean(U)) / sqrt(var(U))
```

In flow cytometry, where per-population cell counts routinely reach thousands or tens of thousands, the normal approximation is effectively exact. R's `wilcox.test()` automatically switches to this approximation for large samples.

### Tie correction

When tied values exist (common in cytometry after transformation and rounding), the variance formula is [adjusted with a tie-correction factor](https://doi.org/10.2307/3001968):

```
var(U) = (n1 * n2 / 12) * [ (N + 1) - sum(t_j^3 - t_j) / (N(N - 1)) ]
```

where t_j is the number of observations in the j-th group of ties. R handles this correction automatically.

---

## Why Non-Parametric for Cytometry

Flow cytometry expression data exhibits distributional properties that make parametric assumptions unreliable:

- **Skewed distributions.** Even after [arcsinh or biexponential (logicle) transformation](https://doi.org/10.1002/cyto.a.20258), marker expression distributions frequently remain asymmetric, with long right tails corresponding to high-expressing subpopulations. [Parks et al. (2006)](https://doi.org/10.1002/cyto.a.20258) demonstrated that no single transformation renders all markers simultaneously normal.

- **Multimodality.** Many markers exhibit [bimodal or multimodal expression patterns](https://doi.org/10.1038/nri2206) corresponding to negative and positive populations. A t-test's normality assumption is fundamentally violated when the distribution has multiple modes.

- **Zero inflation.** Markers with low expression levels often produce a [spike near zero with a separate positive tail](https://doi.org/10.1038/nmeth.2365), violating the continuous symmetry assumed by parametric tests. [Aghaeepour et al. (2013)](https://doi.org/10.1038/nmeth.2365) highlighted this as a common feature across cytometry benchmarking datasets.

- **Heterogeneous variance.** Different experimental conditions and different markers can have [substantially different variances](https://doi.org/10.1002/cyto.a.22625), violating the equal-variance assumption required by the standard t-test. [Finak et al. (2014)](https://doi.org/10.1002/cyto.a.22625) showed that variance heterogeneity is pervasive across cytometry experiments.

The Wilcoxon rank-sum test [makes no assumptions about distributional shape](https://doi.org/10.2307/3001968). By operating on ranks rather than raw values, it is insensitive to monotonic transformations, outliers, and skewness, making it a natural default for exploratory marker screening in cytometry.

---

## When Your Data Calls For It

### Comparing marker expression between two groups

The primary use case: you have two experimental conditions (e.g., responder vs. non-responder, pre-treatment vs. post-treatment, healthy vs. disease) and want to know whether a given marker is expressed at different levels. The Wilcoxon test [directly addresses whether one group's values tend to be larger than the other's](https://doi.org/10.1214/aoms/1177730491).

### Non-normal distributions

When you examine your expression data with histograms or Q-Q plots and see skewness, heavy tails, bimodality, or zero-inflation, the Wilcoxon test is the appropriate default. As discussed above, [cytometry data rarely meets normality assumptions](https://doi.org/10.1002/cyto.a.20258) even after best-practice transformations.

### Exploratory marker screening

When screening dozens of markers across multiple populations, the Wilcoxon test provides a computationally efficient, assumption-light first pass. It can be combined with [Benjamini-Hochberg correction](benjamini-hochberg.md) to control the false discovery rate across all markers tested.

### Per-population, per-marker comparisons

In this pipeline, marker statistics are computed within each identified cell population (cluster). Because each population represents a phenotypically coherent subset, within-population expression distributions tend to be [more unimodal than the full dataset](https://doi.org/10.1038/s42003-019-0415-5). The Wilcoxon test is well-suited to these focused comparisons, providing [robust inference even when within-population distributions remain non-normal](https://doi.org/10.2307/3001968).

---

## When to Choose Something Else

| Scenario | Recommended method | Rationale |
|---|---|---|
| Paired samples (same patient before/after) | **Wilcoxon signed-rank test** | Accounts for within-subject correlation; [analogous to the rank-sum test but for matched pairs](https://doi.org/10.2307/3001968) |
| More than two groups | **Kruskal-Wallis test** | The [non-parametric extension of one-way ANOVA to k groups](https://doi.org/10.1080/01621459.1952.10483441); follow with Dunn's post-hoc test for pairwise comparisons |
| Hierarchical sample structure (cells nested in samples) | **Mixed-effects linear model** (`test_type = "lmer"` in this pipeline) | [Accounts for non-independence of cells within samples](https://doi.org/10.1038/nmeth.4295) via random intercepts; see the pseudoreplication warning below |
| Binary cluster membership outcomes | [**Fisher's Exact Test**](fisher-exact-test.md) | Tests whether cluster proportions differ between conditions using a [2x2 contingency table](https://doi.org/10.1111/j.1469-1809.1934.tb02105.x) |
| Count data (cells per cluster per sample) | [**Negative Binomial GLMM**](negative-binomial-glmm.md) | Models [overdispersed count data with sample-level random effects](https://doi.org/10.1038/s42003-019-0415-5) |
| Multiple testing across many markers | [**Benjamini-Hochberg correction**](benjamini-hochberg.md) | Controls the [false discovery rate when testing many hypotheses simultaneously](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) |

---

## Assumptions and Limitations

### Independence of observations

The Wilcoxon test assumes that [all observations are mutually independent](https://doi.org/10.1214/aoms/1177730491). In cytometry, cells from the same biological sample share sample-level characteristics (donor effects, batch effects, processing variation), which means they are not truly independent. This assumption is routinely violated when cell-level data is tested directly. See the dedicated section on pseudoreplication below.

### Ordinal measurement scale

The test requires that the data be [at least ordinal -- that is, values can be meaningfully ranked](https://doi.org/10.2307/3001968). Fluorescence intensity values, whether raw or transformed, satisfy this requirement.

### Equal shape under H0

The standard interpretation of the Wilcoxon test as a test of location shift (median difference) requires that [the two distributions have the same shape and differ only in location](https://doi.org/10.1214/aoms/1177730491). When this assumption holds, rejection implies that one group's median is shifted relative to the other. When the shapes differ (e.g., different variances or different skewness), the test becomes a more general test of [stochastic dominance](https://doi.org/10.1080/00031305.2017.1305291) -- it detects any tendency for one group's values to exceed the other's, which may not correspond to a simple median shift. In cytometry, where different conditions can alter both the level and spread of expression, [results should be interpreted as stochastic ordering rather than pure location shift](https://doi.org/10.1080/00031305.2017.1305291) whenever distributional shapes visibly differ.

### Sensitivity to ties

When many observations share the same value (extensive ties), the test loses [discriminating power because rank differences are compressed](https://doi.org/10.2307/3001968). In cytometry, ties typically arise from discretization at low signal intensities or from ceiling effects. The tie-corrected variance formula mitigates this, but severely tied data may warrant alternative approaches.

### No covariate adjustment

The basic Wilcoxon test cannot adjust for covariates (batch, age, sex). When confounders are present, [stratified analyses or regression-based approaches](https://doi.org/10.1038/nmeth.4295) are preferable.

---

## The Pseudoreplication Problem

**WARNING: Applying the Wilcoxon rank-sum test to individual cells across conditions will produce misleadingly small p-values when you have multiple samples per condition.**

### Why it matters

In a typical cytometry experiment, each biological sample contributes thousands to millions of cells. When you test condition A vs. condition B at the cell level, you treat each cell as an independent observation. But [cells from the same sample share donor-specific biology, processing variation, and staining intensity](https://doi.org/10.1038/nmeth.4295), so they are not independent. This is the classic problem of **pseudoreplication**: using the number of cells as if it were the number of independent replicates.

[Lun et al. (2017)](https://doi.org/10.1038/nmeth.4295) demonstrated that cell-level tests in cytometry can produce p-values orders of magnitude smaller than justified, leading to false discoveries at alarming rates. [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) showed that methods accounting for sample structure (such as mixed models or aggregation to sample-level summaries) dramatically reduce false positive rates while retaining power to detect genuine differences.

### How it inflates significance

Consider 5 samples per condition, each contributing 10,000 cells. A cell-level Wilcoxon test treats this as n1 = n2 = 50,000, when the effective number of independent observations is only 5 per group. The [resulting test is massively overpowered for detecting inter-sample variation that has nothing to do with the condition of interest](https://doi.org/10.1038/nmeth.4295).

### Recommended solutions

1. **Aggregate to sample-level medians first (preferred for simplicity).** Compute the median expression of each marker within each sample, then run the Wilcoxon test on these sample-level medians. This correctly treats each sample as the unit of replication. [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) recommended this aggregation approach as a practical default for differential expression analysis in cytometry.

   ```r
   # Aggregate to sample-level medians
   sample_medians <- pop_data[, lapply(.SD, median, na.rm = TRUE),
                              by = .(sample_id, condition),
                              .SDcols = markers]

   # Test on sample medians (correct unit of replication)
   wilcox.test(CD4 ~ condition, data = sample_medians)
   ```

2. **Use a mixed-effects model.** Set `test_type = "lmer"` in the pipeline to fit a linear mixed-effects model with sample as a random effect. This [explicitly models the within-sample correlation](https://doi.org/10.1038/nmeth.4295) and produces valid inference at the sample level while using all cell-level data.

   ```r
   # Mixed-effects approach (accounts for sample structure)
   results <- run_marker_statistics(
     input_folder = "output/cluster_data/",
     output_folder = "output/marker_stats/",
     markers = c("CD3", "CD4", "CD8", "CD19"),
     test_type = "lmer"
   )
   ```

3. **Use a bootstrap or permutation test at the sample level.** Permute sample labels (not cell labels) to build a null distribution. This approach is [valid under minimal assumptions and naturally respects the sample-level structure](https://doi.org/10.12688/f1000research.11622.3).

**Bottom line:** The cell-level Wilcoxon test is useful for exploratory screening and visualization, but any findings should be validated with a sample-aware method before being reported as statistically significant. [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) provide a complete workflow demonstrating this two-stage approach.

---

## How to Use in the Pipeline

### Basic invocation

```r
# Source the analysis module
source("src/analysis/marker_statistics.R")

# Compare marker expression between conditions
results <- run_marker_statistics(
  input_folder  = "output/cluster_data/",
  output_folder = "output/marker_stats/",
  filename_regex = "^(.+)_(.+)\\.csv$",
  markers       = c("CD3", "CD4", "CD8", "CD19"),
  test_type     = "wilcox",
  n_workers     = 4
)
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `input_folder` | (required) | Directory containing per-population expression CSVs |
| `output_folder` | (required) | Directory for results output |
| `filename_regex` | `"^(.+)_(.+)\\.csv$"` | Regex with two capture groups: (1) sample ID, (2) population ID |
| `markers` | `NULL` (interactive) | Character vector of marker columns to test; `NULL` triggers GUI selection |
| `test_type` | `"lmer"` | `"wilcox"` for Wilcoxon rank-sum, `"lmer"` for mixed-effects, or `"both"` |
| `n_workers` | `detectCores() - 1` | Number of parallel workers |

### Output

Results are written to `output_folder/marker_statistics.csv` (semicolon-delimited) with one row per marker per population, containing:

- `population`: the cell population (cluster) tested
- `marker`: the marker name
- `mean`: group mean expression
- `n_cells`: total number of cells tested
- `n_samples`: number of contributing samples
- `method`: the test method used (`"wilcox"`, `"lmer"`, or `"simple_mean"` on fallback)

### Recommended workflow

```r
# Step 1: Exploratory cell-level Wilcoxon (fast screening)
results_wilcox <- run_marker_statistics(
  input_folder  = "output/cluster_data/",
  output_folder = "output/marker_stats_wilcox/",
  markers       = c("CD3", "CD4", "CD8", "CD19", "CD25", "FoxP3"),
  test_type     = "wilcox"
)

# Step 2: Confirm hits with mixed-effects model (sample-aware)
results_lmer <- run_marker_statistics(
  input_folder  = "output/cluster_data/",
  output_folder = "output/marker_stats_lmer/",
  markers       = c("CD3", "CD4", "CD8", "CD19", "CD25", "FoxP3"),
  test_type     = "lmer"
)

# Step 3: Apply Benjamini-Hochberg correction
results_lmer$p_adj <- p.adjust(results_lmer$p_value, method = "BH")
significant <- results_lmer[results_lmer$p_adj < 0.05, ]
```

---

## Effect Size

Statistical significance alone is insufficient, especially with the large sample sizes typical of cytometry. Effect size measures quantify the magnitude of the difference.

### Rank-biserial correlation (r)

The rank-biserial correlation is [the natural effect size for the Wilcoxon rank-sum test](https://doi.org/10.1037/h0077767). It is computed directly from the U statistic:

```
r = 1 - (2U) / (n1 * n2)
```

The value ranges from -1 to +1:
- r = +1: every value in group 1 exceeds every value in group 2
- r = 0: no systematic difference in ranks
- r = -1: every value in group 2 exceeds every value in group 1

[Kerby (2014)](https://doi.org/10.2466/11.IT.3.1) provides a tutorial on computing and interpreting the rank-biserial correlation as a complement to the Wilcoxon test.

Conventional thresholds following [Cohen's guidelines adapted for r](https://doi.org/10.4324/9780203771587):

| |r| | Interpretation |
|---|---|
| < 0.1 | Negligible |
| 0.1 -- 0.3 | Small |
| 0.3 -- 0.5 | Medium |
| > 0.5 | Large |

### Cliff's delta

[Cliff's delta](https://doi.org/10.1177/014662168500900302) is an alternative effect size measure that estimates the probability that a randomly selected value from group 1 exceeds a randomly selected value from group 2, minus the reverse probability:

```
delta = (2U) / (n1 * n2) - 1
```

Note that Cliff's delta = -r (the rank-biserial correlation with sign flipped). It ranges from -1 to +1 with the same interpretation but opposite sign convention.

In R:

```r
# Compute rank-biserial correlation from wilcox.test
wt <- wilcox.test(expression ~ condition, data = pop_data)
U  <- wt$statistic
n1 <- sum(pop_data$condition == levels(pop_data$condition)[1])
n2 <- sum(pop_data$condition == levels(pop_data$condition)[2])

r_rank_biserial <- 1 - (2 * U) / (n1 * n2)
cliffs_delta    <- (2 * U) / (n1 * n2) - 1
```

Alternatively, the `effectsize` R package provides [dedicated functions for computing these measures with confidence intervals](https://doi.org/10.21105/joss.02815):

```r
library(effectsize)
rank_biserial(expression ~ condition, data = pop_data)
```

---

## Large Sample Considerations

Flow cytometry datasets routinely contain tens of thousands to millions of cells. This scale creates a specific statistical problem that must be understood.

### Statistical significance is not biological significance

With large N, the Wilcoxon test has [extremely high power to detect arbitrarily small differences](https://doi.org/10.1214/aoms/1177730491). A shift of 0.01 in median fluorescence intensity -- far below the measurement noise floor of most cytometers -- can yield p < 10^-50 when tested across 100,000 cells. This p-value is statistically valid but biologically meaningless.

[Lin et al. (2013)](https://doi.org/10.1038/nature12873) and other large-scale single-cell studies have noted that with sufficient cell counts, essentially every marker will test as significantly different between conditions, rendering unadjusted p-values uninformative for prioritizing biological findings.

### Practical recommendations

1. **Always report effect sizes alongside p-values.** Use the rank-biserial correlation or Cliff's delta (see above). A finding with p < 0.001 but |r| < 0.1 is unlikely to represent a biologically meaningful shift.

2. **Set a minimum effect size threshold.** For exploratory screening, consider requiring both p_adj < 0.05 (after [Benjamini-Hochberg correction](benjamini-hochberg.md)) and |r| > 0.2 to flag a marker as differentially expressed.

3. **Aggregate to sample medians.** As described in the pseudoreplication section, testing on sample-level summaries reduces N to the number of biological replicates, which produces p-values that more accurately reflect the evidence for a condition effect. [Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) demonstrated that sample-aggregated tests maintain good power for genuine effects while avoiding inflation from large cell counts.

4. **Visualize distributions.** Pair statistical tests with density plots or ridge plots of marker expression per condition. Visual inspection often reveals whether a significant p-value reflects a meaningful shift or a trivial one amplified by sample size.

   ```r
   # Ridge plot for visual validation
   library(ggplot2)
   library(ggridges)
   ggplot(pop_data, aes(x = CD4, y = condition, fill = condition)) +
     geom_density_ridges(alpha = 0.6) +
     theme_minimal() +
     labs(title = "CD4 expression by condition",
          x = "CD4 (transformed)", y = "")
   ```

5. **Subsample for computational efficiency.** When cell counts exceed 100,000 per group, subsampling to a fixed number of cells per sample (e.g., 5,000) before running the Wilcoxon test can reduce computation time without meaningfully affecting the detection of large effects. The `src/preprocessing/downsampling.R` module in this pipeline provides utilities for this.

---

## Summary of Recommendations

| Step | Action |
|---|---|
| 1. Screen | Run `test_type = "wilcox"` for a fast, nonparametric scan across all markers and populations |
| 2. Correct | Apply [Benjamini-Hochberg](benjamini-hochberg.md) correction to adjust for multiple testing |
| 3. Filter | Require both p_adj < 0.05 and |effect size| > 0.2 |
| 4. Validate | Confirm significant markers with `test_type = "lmer"` to account for sample structure |
| 5. Visualize | Generate density or ridge plots for all significant markers |
| 6. Report | Present effect sizes, adjusted p-values, and visualizations; avoid reporting cell-level p-values as definitive |
