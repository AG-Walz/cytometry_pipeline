# Methods Guide

This guide covers every statistical, clustering, dimensionality reduction, transformation, normalization, and quality control method available in the cytometry pipeline. Each method has a dedicated deep-dive page with peer-reviewed references, usage examples, and decision guidance.

**How to use this guide:** Start with the [decision flowcharts](#decision-flowcharts) below to identify the right method for your analysis scenario, then follow the link to the deep-dive page for implementation details and scientific background.

---

## Methods at a Glance

### Statistical Methods

| Method | When to Use | Deep Dive |
|--------|------------|-----------|
| **Fisher's Exact Test** | Cluster enrichment with 1 sample/condition or pooled data | [fisher-exact-test.md](methods/fisher-exact-test.md) |
| **Mixed-Effects Logistic Regression** | Cell-level binary cluster membership with multiple replicates | [mixed-effects-logistic.md](methods/mixed-effects-logistic.md) |
| **Negative Binomial GLMM** | Differential abundance (sample-level counts) with replicates | [negative-binomial-glmm.md](methods/negative-binomial-glmm.md) |
| **Wilcoxon Rank-Sum Test** | Comparing continuous marker expression between two groups | [wilcoxon-rank-sum.md](methods/wilcoxon-rank-sum.md) |
| **Benjamini-Hochberg Correction** | Multiple testing correction (applied automatically) | [benjamini-hochberg.md](methods/benjamini-hochberg.md) |
| **Fleiss' Kappa** | Inter-rater agreement for gating comparison (Pipeline D) | [fleiss-kappa.md](methods/fleiss-kappa.md) |

### Clustering

| Method | When to Use | Deep Dive |
|--------|------------|-----------|
| **Leiden Algorithm** | Unsupervised clustering without predefined gates | [leiden-clustering.md](methods/leiden-clustering.md) |
| **Hierarchical Clustering (Ward.D2)** | Ordering samples/clusters in heatmaps and barplots | [hierarchical-clustering.md](methods/hierarchical-clustering.md) |

### Dimensionality Reduction

| Method | When to Use | Deep Dive |
|--------|------------|-----------|
| **UMAP** | 2D visualization of cell populations | [umap.md](methods/umap.md) |
| **PCA** | Variance structure, batch effect QC, noise reduction | [pca.md](methods/pca.md) |

### Data Transformation

| Method | When to Use | Deep Dive |
|--------|------------|-----------|
| **Arcsinh Transform** | CyTOF (mass cytometry) data | [arcsinh-transform.md](methods/arcsinh-transform.md) |
| **Logicle Transform** | Compensated fluorescence flow cytometry | [logicle-transform.md](methods/logicle-transform.md) |
| **Log Transform** | Count data offsets in GLMMs, fold changes | [log-transform.md](methods/log-transform.md) |

### Normalization & Batch Correction

| Method | When to Use | Deep Dive |
|--------|------------|-----------|
| **Harmony Batch Correction** | Multi-batch data with visible batch effects | [batch-correction.md](methods/batch-correction.md) |
| **Quantile Normalization & Z-Score** | Cross-sample distribution alignment, marker scaling | [normalization.md](methods/normalization.md) |

### Preprocessing

| Method | When to Use | Deep Dive |
|--------|------------|-----------|
| **Quality Control (flowAI / PeacoQC)** | Automated removal of aberrant events from FCS data | [quality-control.md](methods/quality-control.md) |
| **Downsampling** | Reducing cell counts for visualization and clustering | [downsampling.md](methods/downsampling.md) |
| **Panel Homogenization** | Unifying marker panels across samples/batches | [panel-homogenization.md](methods/panel-homogenization.md) |

### Visualization Methods

| Method | When to Use | Deep Dive |
|--------|------------|-----------|
| **Venn / Euler / UpSet Plots** | Comparing set overlaps (gating strategies, populations) | [set-visualization.md](methods/set-visualization.md) |

---

## Decision Flowcharts

### Choosing a Statistical Test

```
What is your research question?
|
+-- "Is this cluster enriched in one condition?"
|   |
|   +-- How many samples per condition?
|       |
|       +-- 1 sample per condition (or pooled data)
|       |   --> Fisher's exact test
|       |       [Deep dive](methods/fisher-exact-test.md)
|       |
|       +-- Multiple samples per condition
|           |
|           +-- Cell-level (binary: is each cell in the cluster?)
|           |   --> Mixed-effects logistic regression
|           |       [Deep dive](methods/mixed-effects-logistic.md)
|           |
|           +-- Sample-level (counts per cluster per sample)
|               --> Negative binomial GLMM
|                   [Deep dive](methods/negative-binomial-glmm.md)
|
+-- "Does marker expression differ between conditions?"
|   |
|   +-- No / Unknown normality
|   |   --> Wilcoxon rank-sum test
|   |       [Deep dive](methods/wilcoxon-rank-sum.md)
|   |
|   +-- Yes, account for sample-level variation
|       --> Mixed-effects linear model (lmer)
|
+-- "Do multiple analysts agree on gating?"
|   --> Fleiss' kappa
|       [Deep dive](methods/fleiss-kappa.md)
|
+-- "I am testing many clusters/markers simultaneously"
    --> Apply Benjamini-Hochberg correction to all p-values
        [Deep dive](methods/benjamini-hochberg.md)
```

### Choosing a Transform

```
What type of instrument generated the data?
|
+-- CyTOF (mass cytometry)
|   --> Arcsinh transform (cofactor 5)
|       [Deep dive](methods/arcsinh-transform.md)
|
+-- Fluorescence flow cytometry
|   |
|   +-- Compensated (contains negative values)
|   |   --> Logicle / biexponential transform
|   |       [Deep dive](methods/logicle-transform.md)
|   |
|   +-- No / pre-gated positive data only
|       --> Logicle or arcsinh
|
+-- Count data (cell counts per cluster)
    --> Log transform (as GLM offset)
        [Deep dive](methods/log-transform.md)
```

### Choosing a Clustering Method

```
Do you have predefined gates?
|
+-- Yes (FlowJo workspace or manual gates)
|   --> Import gates via flowjo_import.R
|       For comparison: [Fleiss' Kappa](methods/fleiss-kappa.md),
|                       [Set Visualization](methods/set-visualization.md)
|
+-- No (unsupervised discovery)
    |
    +-- Is data transformed and batch-corrected?
    |   |
    |   +-- No --> Transform first, then correct batch if needed
    |   |         [Arcsinh](methods/arcsinh-transform.md) |
    |   |         [Logicle](methods/logicle-transform.md) |
    |   |         [Batch Correction](methods/batch-correction.md)
    |   |
    |   +-- Yes --> Leiden algorithm
    |               [Deep dive](methods/leiden-clustering.md)
    |               - Resolution 0.5-1.0 for major populations
    |               - Resolution 1.5-3.0 for subpopulations
    |
    +-- Validate with heatmaps (hierarchical clustering for ordering)
        [Deep dive](methods/hierarchical-clustering.md)
```

### Choosing a Dimensionality Reduction Method

```
What is your goal?
|
+-- 2D visualization of single-cell data
|   --> UMAP
|       [Deep dive](methods/umap.md)
|
+-- Variance structure / marker importance
|   --> PCA
|       [Deep dive](methods/pca.md)
|
+-- Batch effect QC
|   --> PCA first (fast, interpretable), then UMAP
|
+-- Pre-processing before clustering
    --> PCA to reduce noise (top N PCs), then Leiden
```

### Do I Need Batch Correction?

```
How many batches/acquisition sessions?
|
+-- Single batch --> No correction needed
|
+-- Multiple batches
    |
    +-- Check PCA/UMAP: do samples cluster by batch?
        |
        +-- No (samples intermix) --> Not needed
        |
        +-- Yes (clear batch separation)
            |
            +-- Batch confounded with biology?
                |
                +-- Yes --> Cannot correct; include batch as covariate in GLMM
                |           [NB-GLMM](methods/negative-binomial-glmm.md)
                |
                +-- No  --> Apply Harmony batch correction
                            [Deep dive](methods/batch-correction.md)
```

### Preprocessing Pipeline Order

```
Raw FCS files
|
+-- 1. Quality Control (remove aberrant events)
|      [Deep dive](methods/quality-control.md)
|
+-- 2. Transformation (arcsinh or logicle)
|      [Arcsinh](methods/arcsinh-transform.md) | [Logicle](methods/logicle-transform.md)
|
+-- 3. Panel Homogenization (if multi-panel)
|      [Deep dive](methods/panel-homogenization.md)
|
+-- 4. Normalization (optional: z-score, quantile)
|      [Deep dive](methods/normalization.md)
|
+-- 5. Batch Correction (if multi-batch)
|      [Deep dive](methods/batch-correction.md)
|
+-- 6. Downsampling (optional, for visualization/clustering)
|      [Deep dive](methods/downsampling.md)
|
+-- 7. Clustering --> 8. Statistical Testing --> 9. Visualization
```

---

## Cross-References

- **[Composability Guide](composability-guide.md)** — Data flow diagram, format contracts, 8 example workflows
- **[Pipeline Guides](pipelines/)** — Step-by-step pipeline walkthroughs (A through E)
- **[Standalone Scripts Reference](standalone-scripts.md)** — Function signatures for all scripts
