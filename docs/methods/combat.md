# ComBat: Empirical Bayes Batch Correction

Parametric empirical Bayes method that estimates and removes location (mean shift) and scale (variance change) batch effects per feature, borrowing strength across features to improve estimates when sample sizes are small.

---

## When to Use

- You have expression matrices from multiple acquisition batches with systematic intensity shifts.
- You want a simple, well-established correction with minimal tuning parameters.
- Your batch effects are approximately uniform across cell populations (i.e., all cell types experience a similar shift).
- You have biological covariates (e.g., condition) that you want to protect during correction.

ComBat is best suited for bulk-level normalization or situations where batch effects are global and uniform. If different cell populations experience different batch shifts (common in cytometry), use [Harmony](batch-correction.md) or [CytoNorm](cytonorm.md) instead, which provide cell-type-aware correction.

---

## How It Works

ComBat models batch effects as additive (location) and multiplicative (scale) distortions applied independently to each feature (marker), then removes them using empirical Bayes shrinkage.

For each marker, the expression in batch b is modeled as: Y_bgj = alpha_g + X * beta_g + gamma_bg + delta_bg * epsilon_bgj, where alpha_g is the overall expression level for marker g, X * beta_g captures biological covariates, gamma_bg is the additive batch effect (location shift), and delta_bg is the multiplicative batch effect (scale change). The error term epsilon_bgj is assumed to follow a normal distribution. By including biological covariates in X, ComBat can protect known biological effects from being removed during correction.

The empirical Bayes component is what distinguishes ComBat from simple per-batch standardization. Rather than estimating gamma and delta independently for each marker, ComBat assumes that the batch effect parameters across all markers follow a common prior distribution (normal for location effects, inverse gamma for scale effects). The prior hyperparameters are estimated from the data (hence "empirical" Bayes), and each marker's batch effect estimate is then shrunk toward the overall mean. This borrowing of information across markers stabilizes the estimates, particularly when the number of samples per batch is small.

[Johnson et al. (2007)](https://doi.org/10.1093/biostatistics/kxj037) originally developed ComBat for microarray gene expression data and showed that the empirical Bayes approach outperformed simple mean-centering and quantile normalization, especially with small batch sizes (3--5 samples per batch). The method has since been widely adopted across omics disciplines. For cytometry, ComBat operates on cells-as-samples, treating each cell's marker vector as an observation. The parametric version assumes normal distributions; a nonparametric variant is available for data that deviates strongly from normality.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `expr_matrix` | (required) | Numeric matrix (cells x markers) |
| `batch_vector` | (required) | Character/factor vector of batch assignments (one per row) |
| `covariates` | NULL | Optional data.frame of biological covariates to preserve (e.g., condition column) |
| `method` | "parametric" | "parametric" (assumes normal priors) or "nonparametric" (kernel density estimation) |

**Using covariates**: When a biological variable (e.g., treatment group) is partially correlated with batch, include it as a covariate to prevent ComBat from removing the biological effect along with the batch effect. Pass a data.frame with one row per cell and columns for each covariate.

**Parametric vs nonparametric**: The parametric method is faster and works well when expression distributions within each batch are approximately normal (typical for arcsinh-transformed cytometry data). Use nonparametric when distributions are heavily skewed or multimodal within batches.

---

## Pipeline Integration

**Script**: `src/harmonization/combat.R`

**Function**: `run_combat(expr_matrix, batch_vector, covariates, method)`

**Input from**:
- `build_combined_expr_matrix()` in `src/R_shared/flowframe_utils.R`
- `combine_expression_matrices()` in `src/analysis/expression_matrices.R`

**Output to**:
- `run_leiden()` in `src/clustering/leiden.R`
- `plot_umap_from_matrix()` in `src/visualization/umap.R`

**Returns**: Corrected numeric matrix with the same dimensions as the input.

**Placement in pipeline**: ComBat operates on transformed expression data. Apply after arcsinh/logicle transformation and before clustering/visualization. The sequence is: Load -> Transform -> ComBat -> Cluster -> Analyze.

**Used in**: Pipelines A and E as an alternative to Harmony batch correction. Implemented via `sva::ComBat()`.

---

## Limitations for Cytometry

ComBat applies a single global correction per marker per batch, which assumes the batch effect magnitude is the same for all cell types. In flow cytometry, this assumption is frequently violated:

- Different cell populations may have different autofluorescence levels, causing population-specific batch shifts.
- Compensation-related artifacts affect markers differently depending on co-expression patterns.
- Detector nonlinearity means a fixed voltage drift causes different intensity shifts at different expression levels.

When cell-type-specific batch effects are present, ComBat can overcorrect some populations while undercorrecting others, potentially distorting the relative expression profiles between cell types. [Harmony](batch-correction.md) and [CytoNorm](cytonorm.md) address this by performing corrections within cell-type-stratified groups.

---

## References

1. Johnson, W. E., Li, C. & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118-127. [doi:10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)

2. Leek, J. T., Johnson, W. E., Parker, H. S. et al. (2012). The sva package for removing batch effects and other unwanted variation in high-throughput experiments. *Bioinformatics*, 28(6), 882-883. [doi:10.1093/bioinformatics/bts034](https://doi.org/10.1093/bioinformatics/bts034)

3. Zhang, Y., Parmigiani, G. & Johnson, W. E. (2020). ComBat-seq: Batch effect adjustment for RNA-seq count data. *NAR Genomics and Bioinformatics*, 2(3), lqaa078. [doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078)
