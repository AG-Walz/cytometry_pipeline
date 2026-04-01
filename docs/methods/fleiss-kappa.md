# Fleiss' Kappa for Inter-Rater Agreement

## Summary

Fleiss' kappa is a statistical measure of inter-rater agreement for categorical
data when two or more raters independently assign categories to the same set of
items. It [extends Cohen's kappa to any number of raters](https://doi.org/10.1037/h0031619),
quantifying how much the observed agreement exceeds what would be expected by
chance alone. In this pipeline, `compute_fleiss_kappa()` (in
`src/metadata_qc/overlap_statistics.R`) wraps `irr::kappam.fleiss()` to assess
whether multiple analysts produce consistent gating assignments when manually
gating the same FCS data in FlowJo.

---

## The Method

### The Agreement Problem

When multiple analysts gate the same flow cytometry sample, each analyst assigns
every event (cell) to a population label -- for example, "CD4+ T cells",
"CD8+ T cells", "B cells", or "ungated". The question Fleiss' kappa answers
is: do the analysts agree on these assignments more than chance would predict?

Simple percent agreement (the fraction of events where all raters assign the
same label) is [misleading because it ignores the agreement expected by
chance](https://doi.org/10.1037/h0031619). If 90% of events are "ungated" and
all raters assign "ungated" to everything, percent agreement is 90% -- but
this tells you nothing about whether they agree on the interesting populations.
Kappa corrects for this by subtracting the chance-expected agreement.

### Mathematical Formulation

[Fleiss (1971)](https://doi.org/10.1037/h0031619) defined kappa for *N* subjects
(events), *n* raters (analysts), and *k* categories (population labels) as
follows.

**Step 1: Build the assignment matrix.** Construct an *N x k* matrix where
entry *n_ij* is the number of raters who assigned event *i* to category *j*.
Each row sums to *n* (the number of raters).

**Step 2: Compute per-event agreement.** For each event *i*, the proportion of
rater pairs that agree is:

```
P_i = (1 / (n * (n - 1))) * (sum_j(n_ij^2) - n)
```

This counts all concordant rater pairs for event *i* and divides by the total
number of rater pairs.

**Step 3: Compute observed agreement.** The overall observed agreement is the
mean of the per-event agreements:

```
P_bar = (1 / N) * sum_i(P_i)
```

**Step 4: Compute expected agreement.** The expected agreement under the null
hypothesis (raters assign categories independently at random, with category
proportions matching the overall marginals) is:

```
P_e = sum_j(p_j^2)
```

where *p_j* is the overall proportion of assignments to category *j* across
all raters and events:

```
p_j = (1 / (N * n)) * sum_i(n_ij)
```

**Step 5: Compute kappa.** Fleiss' kappa is the
[chance-corrected agreement coefficient](https://doi.org/10.1037/h0031619):

```
kappa = (P_bar - P_e) / (1 - P_e)
```

When the raters agree perfectly, *P_bar = 1* and *kappa = 1*. When the raters
agree at exactly the rate expected by chance, *P_bar = P_e* and *kappa = 0*.
When agreement is worse than chance, kappa is negative.

### Standard Error and Significance Testing

[Fleiss (1971)](https://doi.org/10.1037/h0031619) also derived the
[large-sample standard error of kappa under the null hypothesis of no agreement
beyond chance](https://doi.org/10.1037/h0031619). The `irr::kappam.fleiss()`
function uses this standard error to compute a z-statistic and associated
p-value, testing whether the observed kappa is significantly different from
zero. In this pipeline, `compute_fleiss_kappa()` calls `kappam.fleiss()` with
`detail = TRUE`, which additionally reports per-category kappa values --
telling you which specific populations have high or low agreement.

---

## Interpreting Kappa

### The Landis and Koch Scale

[Landis and Koch (1977)](https://doi.org/10.2307/2529310) proposed the following
benchmark scale for interpreting kappa values, which has become the
[most widely cited interpretation framework in biomedical
research](https://doi.org/10.2307/2529310):

| Kappa Range | Agreement Level |
|-------------|-----------------|
| < 0.00      | Less than chance agreement |
| 0.01 -- 0.20 | Slight agreement |
| 0.21 -- 0.40 | Fair agreement |
| 0.41 -- 0.60 | Moderate agreement |
| 0.61 -- 0.80 | Substantial agreement |
| 0.81 -- 1.00 | Almost perfect agreement |

These benchmarks are guidelines, not absolute thresholds.
[Landis and Koch themselves described them as useful but
"clearly arbitrary"](https://doi.org/10.2307/2529310). The acceptable level of
agreement depends on the application. For clinical diagnostic gating where
patient treatment decisions depend on population frequencies, kappa values
above 0.80 are generally expected. For exploratory research gating, values
above 0.60 may be acceptable.

### Per-Category Kappa

When `detail = TRUE` (the pipeline default), `irr::kappam.fleiss()` reports a
kappa value for each category individually. This is essential in cytometry
because [agreement often varies dramatically across
populations](https://doi.org/10.1002/cyto.b.21239). Analysts may agree almost
perfectly on well-separated populations like CD3+ T cells (kappa > 0.90) while
showing only moderate agreement on populations that require subjective gate
placement, such as transitional cell states or dim-positive populations.

### What Kappa Does Not Tell You

Kappa measures whether raters agree with each other, not whether they are
correct. If all analysts consistently misidentify a population (e.g., all gate
too broadly on CD4+ T cells), kappa will be high even though the gating is
wrong. Kappa should be combined with biological validation -- for example,
checking median marker expression within gated populations to confirm that the
assigned labels match expected phenotypes.

---

## Why It Matters for Cytometry

### The Gating Variability Problem

Manual gating in flow cytometry is inherently subjective.
[Maecker et al. (2005)](https://doi.org/10.1186/1471-2172-6-13) demonstrated
in a multi-site study that inter-operator variability in manual gating was one
of the largest sources of variation in flow cytometry assays, often exceeding
biological variability between samples. Even experienced analysts
[place gate boundaries differently](https://doi.org/10.1002/cyto.b.21239),
particularly for populations without sharp bimodal separation.

[Finak et al. (2014)](https://doi.org/10.1038/nmeth.2702), in the
OpenCyto framework paper, noted that manual gating variability undermines the
reproducibility of flow cytometry studies and that standardized, automated
approaches are needed to address this problem. Before implementing automated
alternatives, however, it is essential to quantify how much analysts actually
disagree -- and that is exactly what Fleiss' kappa provides.

### The FlowCAP Challenge

The [FlowCAP (Flow Cytometry: Critical Assessment of Population Identification
Methods) consortium](https://doi.org/10.1038/nmeth.2365) conducted a
systematic evaluation of automated gating algorithms against manual gating by
multiple experts. A key finding was that [the variability among human experts
was comparable to the variability between automated methods and human
experts](https://doi.org/10.1038/nmeth.2365), highlighting that the "gold
standard" of manual gating is itself inconsistent. Fleiss' kappa is the
natural statistic for quantifying this inconsistency in a rigorous,
chance-corrected manner.

### Specific Sources of Gating Disagreement

Several aspects of flow cytometry gating are particularly prone to
inter-operator variability:

1. **Gate boundary placement on continuous distributions.**
   [Maecker et al. (2005)](https://doi.org/10.1186/1471-2172-6-13) showed
   that when populations lack clear bimodal separation, analysts may place
   boundaries hundreds of fluorescence units apart, substantially affecting
   the reported population frequency.

2. **Hierarchical gating order.** The order in which parent gates are applied
   affects downstream populations. An analyst who draws a tight lymphocyte gate
   will report different CD4+ T cell frequencies than one who draws a liberal
   gate, even if their CD4 vs CD8 boundary is identical.

3. **Dim-positive populations.** Populations defined by low-level expression of
   a marker (e.g., PD-1dim, FoxP3low) are [notoriously difficult to gate
   consistently](https://doi.org/10.1002/cyto.b.21239) because the boundary
   between positive and negative is ambiguous.

4. **Rare events.** When a population represents less than 1% of total events,
   [small shifts in gate placement cause large relative changes in the reported
   frequency](https://doi.org/10.1186/1471-2172-6-13), amplifying
   disagreement.

---

## When Your Data Calls for It

### Comparing Gating Strategies Across FlowJo Workspaces

The primary use case in this pipeline. You have two or more FlowJo `.wsp` files
from different analysts (or different gating approaches) applied to the same FCS
data. Pipeline D loads the workspaces, extracts per-event population labels,
and feeds them to `compute_fleiss_kappa()` to produce a single number
summarizing agreement. This answers the question: do different analysts assign
events to the same populations?

### Multiple Analysts Gating the Same Samples

In multi-center clinical trials, [standardizing gating across
sites is critical](https://doi.org/10.1186/1471-2172-6-13). Fleiss' kappa
provides an objective metric for evaluating whether training efforts and
standard operating procedures have achieved acceptable inter-site agreement.
If kappa is below an acceptable threshold (e.g., 0.80 for a clinical endpoint),
the gating protocol may need refinement before data can be pooled across sites.

### QC of Manual Gating Consistency

Even within a single lab, gating consistency can drift over time. Running
Fleiss' kappa periodically on duplicate-gated samples serves as a quality
control check. A declining kappa over successive time points may indicate
[analyst fatigue, instrument drift affecting gate
placement](https://doi.org/10.1186/1471-2172-6-13), or changes in reagent
lot that subtly shift population boundaries.

### Validating Automated Gating Against Manual

When evaluating an automated gating algorithm (e.g., OpenCyto, FlowSOM-based
gating), the automated assignments can be treated as one "rater" alongside
one or more human analysts. Fleiss' kappa then quantifies whether the algorithm
agrees with the human consensus. [Aghaeepour et al. (2013)](https://doi.org/10.1038/nmeth.2365)
used a related approach in the FlowCAP challenge, comparing automated methods
against a panel of expert manual gaters.

---

## When to Choose Something Else

Fleiss' kappa is not universally appropriate. The following table maps common
situations to better alternatives.

| Situation | Problem with Fleiss' Kappa | Better Method |
|-----------|---------------------------|---------------|
| Exactly 2 raters | Fleiss' kappa works but is unnecessarily general | [Cohen's kappa](https://doi.org/10.1177/001316446002000104) provides an equivalent measure specifically designed for the two-rater case, with better-studied properties |
| Ordinal categories (e.g., -, +, ++, +++) | Treats all disagreements equally regardless of degree | [Weighted kappa](https://doi.org/10.1037/h0026256) assigns partial credit for near-agreements (e.g., + vs ++ is a smaller disagreement than - vs +++) |
| Continuous measurements (e.g., MFI values) | Requires categorical data | [Intraclass correlation coefficient (ICC)](https://doi.org/10.1037/1082-989X.1.1.30) measures agreement for continuous outcomes and is the standard choice for comparing marker expression intensities or population frequencies across analysts |
| Unequal number of raters per event | Assumes every event is rated by exactly *n* raters | Krippendorff's alpha handles [missing data and variable numbers of raters per item](https://doi.org/10.1177/0049124104268644) |
| Large class imbalance (prevalence problem) | Kappa can be misleadingly low even with high percent agreement (see [Limitations](#limitations)) | Supplement kappa with the [prevalence-adjusted bias-adjusted kappa (PABAK)](https://doi.org/10.1016/0895-4356(93)90018-V) or report kappa alongside percent agreement and prevalence index |

---

## How to Use in the Pipeline

### Compare FlowJo Workspaces (Pipeline D)

The standard workflow runs `compare_gatings()` first to generate visual
comparisons, then `compute_fleiss_kappa()` for the statistical summary:

```r
# Load workspaces (Step 1 of Pipeline D)
source("src/loading/load_flowjo_wsp.R")
all_labels <- load_multiple_workspaces(
  wsp_paths = c("analyst_A.wsp", "analyst_B.wsp", "analyst_C.wsp"),
  names     = c("Analyst_A", "Analyst_B", "Analyst_C")
)

# Compare gatings visually (Step 2)
source("src/metadata_qc/compare_gatings.R")
combined_labels <- compare_gatings(
  pop_strings_list = all_labels,
  target_fcs       = "sample_001.fcs",
  label_dict       = label_dict,
  custom_name      = c("Analyst_A", "Analyst_B", "Analyst_C"),
  output_dir       = "output/qc/gating_comparison/"
)

# Compute Fleiss' kappa (Step 3)
source("src/metadata_qc/overlap_statistics.R")
results <- compute_fleiss_kappa(
  label_lists   = combined_labels,  # data.frame, cols = raters
  output_dir    = "output/qc/kappa/",
  subset_specs  = list(c(1,2), c(1,3), c(2,3)),
  subset_labels = c("Analyst_A_vs_B", "Analyst_A_vs_C", "Analyst_B_vs_C")
)
```

### Understanding the Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `label_lists` | data.frame | Each column is one rater's categorical labels. All columns must have the same number of rows (events). |
| `output_dir` | character | Directory where `fleiss_kappa_log.txt` is written |
| `subset_specs` | list of integer vectors | Column indices for pairwise or group comparisons. `list(c(1,2))` compares raters 1 and 2 only. |
| `subset_labels` | character vector | Human-readable labels for each subset (used in log output) |
| `row_range` | integer vector | Optional row indices to restrict analysis to a subset of events |

### Reading the Output

The function returns a named list of `irr::kappam.fleiss` result objects and
writes a detailed log to `fleiss_kappa_log.txt`. Each result object contains:

| Field | Description |
|-------|-------------|
| `$value` | The kappa statistic |
| `$statistic` | z-statistic for the test against H0: kappa = 0 |
| `$p.value` | p-value for the significance test |
| `$detail` | Per-category kappa values (when `detail = TRUE`) |

Example of accessing results programmatically:

```r
# Overall kappa for all three raters
results[["Full data"]]$value
# [1] 0.823

# Per-category breakdown
results[["Full data"]]$detail
#                      Kappa         z   p.value
# B cells          0.912  45.2  < 2e-16
# CD4+ T cells     0.847  38.1  < 2e-16
# CD8+ T cells     0.794  33.6  < 2e-16
# NK cells         0.738  28.9  < 2e-16

# Pairwise kappa: Analyst A vs Analyst B
results[["Analyst_A_vs_B"]]$value
# [1] 0.856
```

---

## Relationship to Gating Comparison

Fleiss' kappa provides the statistical backbone to the visual gating comparison
tools in Pipeline D. It is designed to be used alongside the outputs of
`compare_gatings()` (in `src/metadata_qc/compare_gatings.R`), which produces
Venn/Euler diagrams and UpSet plots showing the per-population overlap between
analysts.

### Complementary Information

- **[Set Visualization](set-visualization.md)** (Venn/Euler diagrams, UpSet
  plots) shows *which specific events* are assigned differently and *which
  populations* have the most disagreement. This is qualitative and visual.
- **Fleiss' kappa** gives a single *quantitative summary* of overall agreement,
  corrected for chance. It tells you *how much* the analysts agree but not
  *where* they disagree.

Together, they provide a complete picture: kappa tells you whether agreement
is acceptable, and the set diagrams tell you where to focus harmonization
efforts if it is not.

### Typical Workflow

1. Run `compare_gatings()` to generate Venn/Euler diagrams for each population.
2. Run `compute_fleiss_kappa()` to get the overall and per-category kappa.
3. If kappa is below an acceptable threshold for a population, inspect the
   corresponding Venn/Euler diagram to identify which raters disagree and how.
4. Refine gating SOPs or retrain analysts on the problematic populations.
5. Re-gate and re-run the comparison to confirm improvement.

---

## Limitations

### The Prevalence Problem (Kappa Paradox)

The most widely discussed limitation of kappa is its
[sensitivity to the prevalence of the categories being
rated](https://doi.org/10.1006/cldy.1996.0044). When one category dominates
(e.g., 95% of events are "ungated"), even high percent agreement produces
a low kappa. This occurs because the expected chance agreement *P_e* is
itself high when category prevalence is skewed, leaving little room for kappa
to rise above zero. [Feinstein and Cicchetti (1990)](https://doi.org/10.1006/cldy.1996.0044)
called this the "paradox" of kappa: two raters can agree on 95% of cases
yet produce a kappa of only 0.40.

In cytometry, this paradox is directly relevant.
[Most events in a typical sample are assigned to a small number of dominant
populations](https://doi.org/10.1186/1471-2172-6-13) (e.g., CD4+ T cells
at 40%, ungated at 30%), while the populations of interest may each represent
only 1--5% of events. The overall kappa will be driven by agreement on the
dominant categories. To mitigate this:

- **Always examine per-category kappa** (the `detail` output), which gives
  a kappa for each population individually.
- **Consider prevalence-adjusted bias-adjusted kappa (PABAK)**
  [(Byrt, Bishop & Carlin, 1993)](https://doi.org/10.1016/0895-4356(93)90018-V)
  as a supplementary metric when prevalence imbalance is extreme.
- **Use `row_range` to restrict analysis** to events within a parent gate (e.g.,
  only lymphocytes), which reduces the prevalence of the "ungated" category.

### Sensitivity to Number of Categories

Kappa values tend to [decrease as the number of categories
increases](https://doi.org/10.1037/h0031619), because more categories create
more opportunities for disagreement. A cytometry panel with 20 fine-grained
population labels will typically produce lower kappa values than one with
5 broad lineage labels, even if the underlying agreement is similar. When
comparing kappa across studies, ensure that the number and granularity of
categories are comparable.

### Equal Weighting of All Disagreements

Fleiss' kappa treats all disagreements as equally severe. Assigning a cell to
"CD4+ T cells" instead of "CD8+ T cells" is penalized the same as assigning
it to "CD4+ T cells" instead of "ungated". In cytometry, these are very
different errors -- the first confuses two biologically related populations,
while the second is a gross misclassification.
[Weighted kappa](https://doi.org/10.1037/h0026256) addresses this for ordinal
categories, but for nominal categories (the typical case in gating), there is
no standard weighting scheme. Analysts should supplement kappa with
population-specific confusion matrices to understand the nature of
disagreements, not just their frequency.

### Assumes Fixed Raters

Fleiss' kappa assumes that all *N* events are rated by the same set of *n*
raters. In flow cytometry, this assumption is typically met because each
analyst gates the entire FCS file. However, if some events are excluded by one
analyst (e.g., due to different debris-removal gates) but not others, the
assumption is violated. The pipeline handles this by requiring all columns of
`label_lists` to have the same number of rows; events removed by one analyst
must be explicitly labeled (e.g., as "Removed") rather than dropped.

### Does Not Quantify Bias Direction

Kappa tells you whether raters agree, but not who is the outlier or in which
direction they deviate. If Analyst A consistently gates CD4+ T cells more
broadly than Analysts B and C, kappa will detect the disagreement but will
not identify Analyst A as the source. Pairwise kappa comparisons (via
`subset_specs`) can help isolate which analyst pairs disagree, and the
Venn/Euler diagrams from `compare_gatings()` provide the visual detail
needed to identify the direction of disagreement.

### Large-Sample Approximation for Significance Testing

The p-value reported by `irr::kappam.fleiss()` relies on a
[large-sample normal approximation to the distribution of
kappa](https://doi.org/10.1037/h0031619). For very small datasets (fewer
than ~30 events), this approximation may be unreliable. In cytometry, where
event counts typically range from thousands to hundreds of thousands, this is
rarely a practical concern.

---

## Related Documentation

- [Set Visualization](set-visualization.md) -- Venn/Euler diagrams and UpSet plots from `compare_gatings.R`
- [Pipeline D: FlowJo Comparison](../pipelines/pipeline-D-flowjo-comparison.md) -- Full pipeline using Fleiss' kappa
- [Composability Guide](../composability-guide.md) -- How to combine pipeline modules
