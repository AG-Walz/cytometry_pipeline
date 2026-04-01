# diffcyt: Differential Discovery Framework for Cytometry

Integrated framework that combines FlowSOM clustering with established statistical methods (edgeR, voom, limma, GLMM) to perform differential abundance (DA) and differential state (DS) analysis on high-dimensional cytometry data.

---

## When to Use

- You want a standardized, well-validated workflow for differential testing in cytometry experiments.
- You need both differential abundance (which populations change in frequency?) and differential state (which markers change in expression?) in a single framework.
- Your experimental design includes biological replicates, paired samples, or covariates.
- You prefer using established genomics statistical methods (edgeR, limma) adapted for cytometry rather than custom implementations.

If you need more control over clustering (e.g., Leiden instead of FlowSOM) or statistical models (e.g., negative binomial GLMM with custom random effects), use the individual scripts in `src/analysis/` directly.

---

## How It Works

The diffcyt framework follows a three-stage pipeline: clustering, feature extraction, and differential testing.

In the clustering stage, diffcyt uses FlowSOM to partition all cells across all samples into a fixed number of clusters (meta-clusters). FlowSOM is chosen for its speed and scalability (see [FlowSOM](flowsom.md)). The user specifies the number of meta-clusters; diffcyt defaults to 20 but this should be tuned based on panel complexity. Each cell receives a cluster assignment, and all downstream statistics operate on per-cluster, per-sample summaries rather than individual cells.

In the feature extraction stage, two types of features are computed. For differential abundance (DA), the feature is the cell count per cluster per sample, optionally normalized to proportions. For differential state (DS), the features are the median marker expression values per cluster per sample. These summary statistics reduce millions of individual cell measurements to a manageable matrix (clusters x samples) suitable for standard statistical models.

In the differential testing stage, the summarized features are passed to established statistical frameworks. For DA, the options are: edgeR (negative binomial GLM on counts), voom (variance-stabilized linear model), or a GLMM with binomial or negative binomial family. For DS, the options are: limma (moderated t-tests on median expression) or a linear mixed model (LMM). All methods support complex experimental designs including paired samples, blocking factors, and continuous covariates via design matrices. Multiple testing correction is applied across clusters using the Benjamini-Hochberg procedure.

[Weber et al. (2019)](https://doi.org/10.1038/s42003-019-0415-5) showed that diffcyt methods consistently outperformed cell-level approaches and other cluster-based methods in comprehensive benchmarks with known ground truth, particularly for detecting subtle abundance changes in the presence of biological variability.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `flow_frames` | (required) | List of flowFrame objects (one per sample) |
| `experiment_info` | (required) | data.frame with `sample_id`, `group_id`, `patient_id` columns |
| `markers_of_interest` | (required) | Character vector of markers for clustering and/or testing |
| `analysis_type` | "DA" | "DA" (differential abundance) or "DS" (differential state) |
| `method_DA` | "diffcyt-DA-edgeR" | DA method: "diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM" |
| `method_DS` | "diffcyt-DS-limma" | DS method: "diffcyt-DS-limma", "diffcyt-DS-LMM" |
| `n_clusters` | 20 | Number of FlowSOM meta-clusters |
| `seed` | 42 | Random seed for reproducibility |

**Choosing DA method**: Use edgeR for simple group comparisons with biological replicates. Use voom when you want limma's moderated statistics on count data. Use GLMM for complex designs with random effects (e.g., paired/longitudinal data).

**Choosing DS method**: Use limma for most cases (fast, robust). Use LMM when you need explicit random effects modeling.

---

## Pipeline Integration

**Script**: `src/analysis/diffcyt_analysis.R`

**Function**: `run_diffcyt(flow_frames, experiment_info, markers_of_interest, ...)`

**Input from**:
- `load_fcs_raw()` in `src/loading/load_fcs_raw.R`
- `build_combined_expr_matrix()` in `src/R_shared/flowframe_utils.R`

**Output to**:
- Standalone results CSV with per-cluster p-values, adjusted p-values, and log fold changes
- `plot_volcano()` in `src/analysis/differential_abundance.R` for visualization

**Returns**: List containing `results` (full diffcyt output), `top_table` (significant clusters), and `cluster_assignments` (per-cell cluster labels).

**Used in**: Pipelines A and B as an alternative to the manual clustering + GLMM workflow. Particularly useful when you want a self-contained analysis that handles clustering and statistics in one call.

---

## References

1. Weber, L. M., Nowicka, M., Soneson, C. & Robinson, M. D. (2019). diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. *Communications Biology*, 2, 183. [doi:10.1038/s42003-019-0415-5](https://doi.org/10.1038/s42003-019-0415-5)

2. Robinson, M. D., McCarthy, D. J. & Smyth, G. K. (2010). edgeR: A Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139-140. [doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

3. Nowicka, M., Krieg, C., Crowell, H. L. et al. (2017). CyTOF workflow: Differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*, 6, 748. [doi:10.12688/f1000research.11622.3](https://doi.org/10.12688/f1000research.11622.3)
