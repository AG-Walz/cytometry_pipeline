# Downsampling Methods for Flow Cytometry Data

## Summary

Downsampling reduces the number of cells in a dataset by random subsampling, making computationally intensive analyses (UMAP, clustering) tractable while preserving population structure. This pipeline implements two strategies in `src/preprocessing/downsampling.R`: **Consistent** downsampling, which draws an equal number of cells from every sample for balanced cross-sample comparison, and **Maximal** (proportional) downsampling, which preserves each sample's relative contribution to the total dataset.

---

## Why Downsampling

Modern flow and mass cytometry instruments routinely acquire [millions of events per experiment across dozens of samples](https://doi.org/10.1016/j.cell.2015.05.047), producing datasets that can exceed the memory and computational limits of downstream algorithms. Downsampling addresses three practical constraints:

### Memory

High-dimensional cytometry datasets stored as expression matrices consume memory proportional to cells x markers. A 20-sample experiment with 500,000 cells per sample and 30 markers produces a matrix of 300 million entries. Algorithms that require pairwise distance computation or dense matrix operations can [exhaust available RAM at these scales](https://doi.org/10.1016/j.cell.2015.05.047). Downsampling to a manageable total cell count (e.g., 50,000--100,000 cells) reduces the matrix size by orders of magnitude.

### Visualization Clarity

Dimensionality reduction methods like [UMAP (McInnes et al. 2018)](https://doi.org/10.48550/arXiv.1802.03426) and [t-SNE (van der Maaten and Hinton 2008)](https://doi.org/10.1007/s10994-008-5065-1) produce 2D scatter plots where individual points represent cells. When millions of cells are plotted, [overplotting obscures population structure and rare subsets become invisible beneath dense clouds of abundant populations](https://doi.org/10.1038/nbt.4314). [Becht et al. (2019)](https://doi.org/10.1038/nbt.4314) demonstrated that UMAP embeddings of 50,000--100,000 cells capture the same population structure visible at full scale, while producing substantially clearer visualizations.

### Computational Tractability

Graph-based algorithms central to modern cytometry workflows scale superlinearly with cell count. The KNN graph construction underlying both [UMAP](https://doi.org/10.48550/arXiv.1802.03426) and [Leiden clustering](https://doi.org/10.1038/s41598-019-41695-z) requires O(N log N) time for approximate nearest-neighbor search. The subsequent optimization steps scale as O(N x k) for graph construction and O(N x epochs) for embedding optimization. At one million cells, [UMAP computation can take tens of minutes even with approximate methods](https://doi.org/10.1038/nbt.4314), while Leiden community detection on the resulting graph adds further computational cost. [Levine et al. (2015)](https://doi.org/10.1016/j.cell.2015.05.047) noted that PhenoGraph runtime increases substantially with cell count and recommended downsampling for initial exploratory analyses.

---

## Consistent Downsampling

Consistent downsampling draws an **equal number of cells from every sample**, regardless of the original sample sizes. The per-sample cell count is calculated as:

```
per_sample = min(target_cells / n_samples, min(cell_counts))
```

This ensures no sample contributes more cells than any other, and no sample is asked to contribute more cells than it contains.

### How It Works

1. The target total cell count is divided equally among all samples.
2. If any sample has fewer cells than its equal share, the per-sample count is capped at the smallest sample's cell count to maintain strict equality.
3. Each sample is randomly subsampled (without replacement) to exactly `per_sample` cells.

### Example

Given 5 samples with cell counts [80,000, 120,000, 45,000, 200,000, 95,000] and a target of 50,000 total cells:

- Equal share: 50,000 / 5 = 10,000 per sample
- Smallest sample: 45,000 (larger than 10,000, so no cap needed)
- Result: 10,000 cells from each sample, 50,000 total

### Advantages

- **Balanced representation.** Every sample has identical weight in downstream analyses. This prevents large samples from dominating UMAP layouts, clustering results, or marker expression summaries. [Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) noted that unbalanced sample representation can bias clustering algorithms toward populations abundant in the largest samples.
- **Fair cross-sample comparison.** When overlaying sample identity on UMAP or comparing cluster compositions across samples, equal representation ensures visual and statistical comparability.
- **Prevents batch-size confounding.** In experiments where acquisition volume varies across conditions (e.g., tumor samples yield fewer cells than blood), consistent downsampling removes the confound between sample size and biological condition.

### Disadvantages

- **Constrained by the smallest sample.** If one sample has very few cells, the per-sample count for all samples is capped at that minimum. A single low-yield sample can drastically reduce the total number of cells available for analysis.
- **Loss of rare populations.** Aggressive downsampling of large samples may remove cells from rare populations (< 1% frequency). A population present at 0.5% in a sample of 200,000 cells (1,000 cells) is reduced to approximately 5 cells when downsampled to 10,000, making it [statistically undetectable and potentially absent entirely due to random sampling](https://doi.org/10.1002/cyto.a.22764).
- **Discards information.** Equal representation means that large, information-rich samples are truncated to match the smallest sample, wasting acquired data.

---

## Maximal (Proportional) Downsampling

Maximal downsampling draws cells from each sample **in proportion to its original size**, preserving the natural distribution of cells across samples. The per-sample cell count is calculated as:

```
fraction = cell_counts / total_cells
per_sample = round(fraction * target_cells)
```

Each sample's allocation is capped at its actual cell count (no upsampling).

### How It Works

1. Each sample's fraction of the total dataset is computed.
2. The target total is distributed according to these fractions.
3. If a sample's allocation exceeds its actual cell count (possible when one sample is very small relative to the target), it contributes all its cells.
4. Each sample is randomly subsampled (without replacement) to its allocated count.

### Example

Given 5 samples with cell counts [80,000, 120,000, 45,000, 200,000, 95,000] (total: 540,000) and a target of 50,000:

| Sample | Original | Fraction | Allocation |
|--------|----------|----------|------------|
| S1 | 80,000 | 0.148 | 7,407 |
| S2 | 120,000 | 0.222 | 11,111 |
| S3 | 45,000 | 0.083 | 4,167 |
| S4 | 200,000 | 0.370 | 18,519 |
| S5 | 95,000 | 0.176 | 8,796 |
| **Total** | 540,000 | 1.000 | **50,000** |

### Advantages

- **Preserves natural sample sizes.** The relative contribution of each sample to the pooled dataset reflects the original experiment. This is important when [sample size itself carries biological meaning](https://doi.org/10.1016/j.cell.2015.05.047) --- for instance, a tumor biopsy yielding more cells may genuinely contain more immune infiltrate.
- **Better rare population retention in large samples.** Large samples retain more cells, preserving low-frequency populations that would be lost under consistent downsampling. [Levine et al. (2015)](https://doi.org/10.1016/j.cell.2015.05.047) found that proportional representation helped retain rare progenitor populations during exploratory PhenoGraph analysis.
- **No cap from small samples.** A single low-yield sample does not penalize the entire dataset. Small samples contribute proportionally fewer cells while large samples retain their information content.

### Disadvantages

- **Imbalanced representation.** Large samples dominate the pooled data. In UMAP visualizations, [populations from large samples will be overrepresented, visually obscuring populations from smaller samples](https://doi.org/10.1038/nbt.4314). In Leiden clustering, the algorithm will preferentially discover substructure in large samples because they contribute more cells to the KNN graph.
- **Confounds sample size with biology.** Downstream abundance analyses must account for unequal sample sizes. If a cluster appears to be more abundant in one condition simply because those samples contributed more cells, the result is an artifact.
- **Small samples may be underrepresented.** A sample contributing only 4% of the total may have too few cells after downsampling to reliably represent its population structure.

---

## Choosing Between Methods

The choice between Consistent and Maximal downsampling depends on the analytical question.

### Use Consistent Downsampling When

| Scenario | Rationale |
|----------|-----------|
| **Comparing populations across samples** | Equal representation ensures that differences in cluster composition reflect biology, not sample size |
| **Differential abundance analysis** | Statistical tests comparing population frequencies assume comparable denominators across samples. [Mixed-effects logistic regression and negative binomial GLMMs](https://doi.org/10.12688/f1000research.11622.3) benefit from balanced group sizes for stable parameter estimation |
| **UMAP colored by sample** | Equal cell counts per sample produce visually balanced overlays where no single sample dominates the plot |
| **Multi-condition experiments** | Prevents the condition with more cells from dominating the pooled analysis |

### Use Maximal (Proportional) Downsampling When

| Scenario | Rationale |
|----------|-----------|
| **Preserving natural distribution** | When sample size is biologically meaningful (e.g., tumor immune infiltrate) |
| **Initial clustering / population discovery** | Retaining more cells from large samples increases sensitivity to rare populations. [PhenoGraph](https://doi.org/10.1016/j.cell.2015.05.047) benefits from larger input sizes for stable graph construction |
| **Abundance estimation** | When downstream analysis requires that relative sample contributions reflect reality |
| **Exploratory analysis** | When the priority is maximum information retention and visualization balance is secondary |

### Decision Flowchart

1. **Are you comparing populations across samples or conditions?**
   - Yes -> Consistent
   - No -> Continue
2. **Does sample size carry biological meaning?**
   - Yes -> Maximal
   - No -> Continue
3. **Is one sample much smaller than the others?**
   - Yes, and you need balanced analysis -> Consistent (but consider excluding the small sample)
   - Yes, but you want maximum data retention -> Maximal
4. **Default:** Maximal (the pipeline default) for general-purpose exploratory analysis.

---

## Target Cell Count

The target cell count determines the trade-off between computational efficiency and analytical resolution. There is no single correct value; it depends on the analysis stage and the biological question.

### Guidelines

| Target | Use Case | Rationale |
|--------|----------|-----------|
| **5,000--10,000** | Quick exploration, parameter tuning | Fast UMAP and clustering iteration. Sufficient to see major populations. Good for testing pipeline parameters before committing to a full run. |
| **30,000--50,000** | Standard analysis | Balances computational cost with resolution. [Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) showed that clustering algorithms [reliably recover major and most minor populations at these cell counts](https://doi.org/10.1002/cyto.a.22764). |
| **50,000--100,000** | Publication-quality UMAP | Sufficient density for visually compelling embeddings without overplotting. [Becht et al. (2019)](https://doi.org/10.1038/nbt.4314) used datasets in this range for UMAP benchmarking. |
| **100,000--500,000** | Rare population detection | Needed when populations of interest are < 1% frequency. At 0.5% frequency, 100,000 cells yields approximately 500 cells from the rare population --- [enough for reliable cluster detection](https://doi.org/10.1016/j.cell.2015.05.047). |
| **Full dataset** | Statistical testing, final abundance quantification | Statistical tests (Fisher's exact, GLMM) should use all available cells for maximum statistical power. See [When NOT to Downsample](#when-not-to-downsample). |

### Practical Considerations

- **Memory ceiling.** A 100,000 x 30 expression matrix (float64) consumes approximately 24 MB --- trivial for modern workstations. At 1 million cells, the matrix itself is 240 MB, but KNN graph construction and UMAP optimization may require several GB of working memory.
- **Runtime scaling.** UMAP runtime scales approximately as O(N^1.2) for the nearest-neighbor search phase. Doubling the cell count roughly doubles the wall-clock time. [Becht et al. (2019)](https://doi.org/10.1038/nbt.4314) reported [UMAP runtimes of under one minute for 100,000 cells and several minutes for 1 million cells](https://doi.org/10.1038/nbt.4314).
- **Per-sample minimum.** Regardless of the total target, ensure each sample retains at least a few hundred cells. Below this threshold, per-sample population structure cannot be reliably estimated. The Consistent method enforces this automatically by capping at the smallest sample.

---

## When Your Data Calls For It

Downsampling is appropriate in the following concrete scenarios:

### Before UMAP or t-SNE

Dimensionality reduction is the most common reason to downsample. [UMAP's approximate nearest-neighbor search scales as O(N log N)](https://doi.org/10.48550/arXiv.1802.03426), and the stochastic gradient descent optimization scales linearly with cell count. Beyond computational cost, [visual overplotting in UMAP embeddings of millions of cells obscures population boundaries and rare subsets](https://doi.org/10.1038/nbt.4314). Downsampling to 50,000--100,000 cells produces clearer, more interpretable embeddings without loss of major population structure.

### Before Leiden Clustering (Large Datasets)

While [Leiden clustering scales well to large datasets](https://doi.org/10.1038/s41598-019-41695-z), the KNN graph construction step requires computing distances between all cell pairs within each neighborhood. For datasets exceeding 500,000 cells, [downsampling followed by label propagation to the full dataset](https://doi.org/10.1016/j.cell.2015.05.047) is a practical alternative: cluster on the downsampled subset, then assign each remaining cell to the cluster of its nearest neighbor in the downsampled data.

### For Exploratory Parameter Tuning

When iterating on UMAP parameters (`n_neighbors`, `min_dist`) or Leiden resolution, running the full dataset at each iteration is wasteful. Downsample to 10,000 cells for rapid parameter sweeps, then run the final analysis at full scale or a higher target.

### For Multi-Panel Experiments

When combining data from multiple panels with different marker sets, [panel homogenization](../subprojects/harmonization.md) may produce a shared marker matrix. Downsampling ensures that panels contributing different numbers of cells are represented appropriately in the merged analysis.

---

## When NOT to Downsample

Downsampling is counterproductive in specific situations. Using all available cells is preferable when:

### Statistical Tests

[Fisher's exact test](fisher-exact-test.md), [mixed-effects logistic regression](mixed-effects-logistic.md), and [negative binomial GLMMs](negative-binomial-glmm.md) derive their power from sample size. Downsampling reduces the number of cells available for testing, [inflating p-values and reducing the ability to detect true differences in population frequencies](https://doi.org/10.12688/f1000research.11622.3). Statistical tests should always be performed on the full dataset. In this pipeline, the analysis modules in `src/analysis/` operate on complete cluster-assignment tables produced by `clusters_to_gating_csv()`, not on downsampled data.

### Rare Population Detection

A population present at 0.1% frequency contains 1,000 cells in a sample of 1,000,000 but only 10 cells after downsampling to 10,000. At 10 cells, the population is [below the detection threshold for most clustering algorithms](https://doi.org/10.1002/cyto.a.22764) and may be entirely absent due to random sampling. If rare populations are the primary analytical target, either skip downsampling or use density-dependent downsampling (see [Alternative Approaches](#alternative-approaches)) that [preferentially retains cells in low-density regions](https://doi.org/10.1038/nbt.1991).

### Small Datasets

Datasets with fewer than approximately 50,000 total cells across all samples are already within the computational budget of UMAP and Leiden clustering. Downsampling provides no performance benefit and only removes information. The pipeline's downsampling function will operate correctly on small datasets (returning all cells when the target exceeds the actual count), but calling it is unnecessary.

### Final Abundance Quantification

When reporting the frequency of each cluster per sample (e.g., "cluster 5 comprises 12.3% of sample A"), use all cells for accurate denominators. Downsampled frequencies have [higher variance due to the reduced denominator and may not match the true proportions](https://doi.org/10.1002/cyto.a.22764), especially for rare populations.

---

## Reproducibility

Both downsampling methods rely on random subsampling, making the `seed` parameter essential for reproducible results.

### Why Seeds Matter

Random subsampling selects a different subset of cells each time it is called unless the random number generator is explicitly seeded. Different subsets lead to different UMAP embeddings, different cluster compositions, and potentially different biological conclusions. [Non-reproducible preprocessing undermines the ability to verify results or build on prior analyses](https://doi.org/10.12688/f1000research.11622.3).

### How the Pipeline Ensures Reproducibility

The `downsample_flowframes()` function calls `set.seed(seed)` at the start, before any random sampling occurs:

```r
downsample_flowframes <- function(flow_frames, sample_names,
                                   target_cells = 10000,
                                   method = "Maximal",
                                   seed = 42, ...) {
  set.seed(seed)
  # All subsequent sample() calls are deterministic
  ...
}
```

With the same `seed`, `flow_frames`, and `target_cells`, the function produces identical output across runs, R sessions, and machines (given the same R version and platform).

### Reporting Checklist

When publishing analyses that involved downsampling, always report:

- **Method:** Consistent or Maximal
- **Target cell count:** The value of `target_cells`
- **Actual cell count:** The total number of cells after downsampling (may differ slightly from target due to rounding or capping)
- **Random seed:** The exact seed used
- **Per-sample counts:** Original and downsampled cell counts for each sample
- **Software:** R version and flowCore version

This information allows collaborators to reproduce the exact same downsampled dataset.

---

## How to Use in the Pipeline

### Proportional Downsampling (Default)

```r
source("src/preprocessing/downsampling.R")

downsampled <- downsample_flowframes(
  flow_frames  = transformed_ffs,    # Named list of flowFrame objects
  sample_names = sample_names,       # Character vector of sample identifiers
  target_cells = 50000,              # Total cells after downsampling
  method       = "Maximal",          # Proportional to original size
  seed         = 42                  # Reproducible random seed
)

# downsampled is a named list of flowFrame objects
# Each frame contains a proportional share of target_cells
```

### Equal Representation

```r
source("src/preprocessing/downsampling.R")

downsampled <- downsample_flowframes(
  flow_frames  = transformed_ffs,
  sample_names = sample_names,
  target_cells = 50000,
  method       = "Consistent",       # Equal cells per sample
  seed         = 42
)

# Each sample contributes exactly min(50000/n_samples, min(cell_counts)) cells
```

### Saving Downsampled FCS Files

To export the downsampled data as FCS files for use in external tools (e.g., FlowJo):

```r
downsampled <- downsample_flowframes(
  flow_frames  = transformed_ffs,
  sample_names = sample_names,
  target_cells = 50000,
  method       = "Maximal",
  seed         = 42,
  save_fcs     = TRUE,
  output_dir   = "output/downsampled_fcs/"
)

# Files saved as: output/downsampled_fcs/ds_{sample_name}.fcs
```

### Function Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | (required) | Named list of flowFrame objects |
| `sample_names` | character | (required) | Vector of sample identifiers, same length as `flow_frames` |
| `target_cells` | integer | 10000 | Target total number of cells after downsampling |
| `method` | character | "Maximal" | `"Consistent"` for equal per-sample counts, `"Maximal"` for proportional |
| `seed` | integer | 42 | Random seed for reproducibility |
| `save_fcs` | logical | FALSE | Whether to save downsampled flowFrames as FCS files |
| `output_dir` | character | NULL | Directory for saved FCS files (required if `save_fcs = TRUE`) |

### Return Value

A named list of flowFrame objects, one per sample, with cell counts reduced according to the selected method. Sample names are preserved as list names.

---

## Downsampling in the Pipeline

Downsampling occupies a specific position in the preprocessing pipeline. It should be applied **after** data quality and transformation steps and **before** computationally intensive analysis steps.

### Pipeline Position

```
FCS Loading
    |
    v
Quality Control (flowAI / PeacoQC)     -- see quality-control.md
    |
    v
Transformation (arcsinh / logicle)      -- see arcsinh-transform.md, logicle-transform.md
    |
    v
>>> Downsampling <<<                    -- this step
    |
    v
Panel Homogenization                    -- see harmonization docs
    |
    v
Batch Correction                        -- see harmonization docs
    |
    v
Leiden Clustering                       -- see leiden-clustering.md
    |
    v
UMAP Visualization                      -- see umap.md
    |
    v
Statistical Analysis (full data)        -- no downsampling
```

### Why This Order

- **After QC:** Quality control removes aberrant events (doublets, debris, acquisition anomalies). Downsampling before QC wastes the downsampling budget on events that will be discarded, and the QC step may unevenly remove cells across samples, producing unbalanced representation even after consistent downsampling.
- **After transformation:** Transformation does not change cell count; it only rescales expression values. Order relative to downsampling does not affect the result, but performing transformation on the full dataset ensures that transformation parameters (e.g., logicle width estimation) are [computed from the complete data distribution rather than a subset](https://doi.org/10.1002/cyto.a.22625).
- **Before clustering and UMAP:** These are the computationally expensive steps that benefit most from reduced cell count. [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) recommended downsampling before dimensionality reduction and clustering in their standard CyTOF analysis workflow.

### Related Documentation

- [Quality Control](quality-control.md) --- Upstream: remove aberrant events before downsampling
- [Arcsinh Transform](arcsinh-transform.md) --- Upstream: transform intensities before downsampling
- [Logicle Transform](logicle-transform.md) --- Upstream: transform fluorescence data before downsampling
- [Leiden Clustering](leiden-clustering.md) --- Downstream: cluster the downsampled data
- [UMAP](umap.md) --- Downstream: visualize the downsampled data

---

## Alternative Approaches

The two methods implemented in this pipeline (Consistent and Maximal) perform uniform random subsampling within each sample. Alternative strategies exist that use non-uniform sampling to address specific analytical needs.

### Density-Dependent Downsampling (SPADE)

[Qiu et al. (2011)](https://doi.org/10.1038/nbt.1991) introduced density-dependent downsampling as part of the SPADE (Spanning-tree Progression Analysis of Density-normalized Events) algorithm. Instead of sampling cells uniformly at random, SPADE [estimates the local density around each cell and preferentially retains cells in low-density regions while aggressively downsampling high-density regions](https://doi.org/10.1038/nbt.1991). This produces a dataset where rare populations are overrepresented relative to their true frequency, improving their visibility in downstream analyses.

**Advantages:**
- [Rare cell populations that would be lost under uniform random downsampling are preserved](https://doi.org/10.1038/nbt.1991)
- Produces a more uniform density landscape, which can improve clustering and visualization of rare events
- Particularly useful for [detecting rare progenitor populations in hematopoietic differentiation studies](https://doi.org/10.1038/nbt.1991)

**Disadvantages:**
- [Distorts the true frequency distribution of cell populations](https://doi.org/10.1002/cyto.a.22764), making it unsuitable for abundance analysis
- Requires kernel density estimation, which adds computational cost and introduces bandwidth parameters
- [Weber and Robinson (2016)](https://doi.org/10.1002/cyto.a.22764) found that SPADE's density-dependent downsampling, while beneficial for rare population retention, [degraded performance in cluster-frequency estimation benchmarks](https://doi.org/10.1002/cyto.a.22764)

**When to consider:** If your primary goal is detecting rare populations (< 0.1% frequency) and you are willing to sacrifice accurate frequency estimation.

### Geometric Sketching

[Hie et al. (2019)](https://doi.org/10.1016/j.cels.2019.05.003) proposed geometric sketching, a method that [selects a maximally representative subset of cells by partitioning the high-dimensional expression space into equal-volume geometric regions and sampling one cell from each region](https://doi.org/10.1016/j.cels.2019.05.003). This approach is grounded in the theory of epsilon-nets and kernel herding, ensuring that the sketch [approximates the full dataset's distribution with formal guarantees on coverage](https://doi.org/10.1016/j.cels.2019.05.003).

**Advantages:**
- [Provides theoretical guarantees that rare and common populations are both represented in the sketch](https://doi.org/10.1016/j.cels.2019.05.003)
- Outperformed uniform random sampling and density-dependent downsampling in [preserving cell-type diversity in single-cell RNA-seq benchmarks](https://doi.org/10.1016/j.cels.2019.05.003)
- Deterministic (no random seed needed) when the partitioning is fixed

**Disadvantages:**
- Computationally more expensive than uniform random sampling due to the geometric partitioning step
- Primarily developed and benchmarked for scRNA-seq; less extensively validated on protein-level cytometry data
- Implementation requires additional dependencies not currently included in this pipeline

**When to consider:** When maximum population diversity in a small sample is critical and computational overhead is acceptable.

### Stratified Downsampling

A pragmatic middle ground is to first cluster the full dataset (e.g., with [FlowSOM](https://doi.org/10.1002/cyto.a.22625), which is fast enough to run on millions of cells) and then sample cells proportionally from each cluster. This ensures that every identified population is represented in the downsampled data, regardless of its frequency in the original dataset. [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) used a similar approach in their recommended CyTOF workflow, performing initial overclustering with FlowSOM followed by downsampling for visualization.

**Advantages:**
- Preserves all detected populations, including rare ones
- Straightforward to implement with existing clustering tools
- Population frequencies can be made either proportional or equal, depending on the analytical goal

**Disadvantages:**
- Requires an initial clustering step, which adds computational cost and introduces dependence on clustering parameters
- The quality of the stratification depends on the quality of the initial clustering

### Why This Pipeline Uses Uniform Random Sampling

Uniform random sampling (both Consistent and Maximal) is the simplest approach, introduces no distributional bias, and requires no additional parameters beyond the target cell count and the method choice. For the vast majority of immunophenotyping experiments where major and minor populations (> 1% frequency) are of interest, [uniform downsampling to 50,000--100,000 cells reliably preserves population structure](https://doi.org/10.1002/cyto.a.22764). The pipeline's modular design allows density-dependent or geometric methods to be substituted in `src/preprocessing/downsampling.R` if specialized needs arise.
