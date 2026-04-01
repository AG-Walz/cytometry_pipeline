[< Back to Documentation Index](../README.md)

# src/analysis/ -- Statistical Analysis

Twelve scripts for cluster splitting, differential abundance testing, differential expression, cluster enrichment, marker statistics, effect sizes, marker correlations, diffcyt framework analysis, survival analysis, expression matrix generation, abundance plotting, and cell ratio calculation. These scripts are the statistical backbone of the pipeline.

## Summary Table

| Script | Functions | Description |
|--------|-----------|-------------|
| `cluster_splitting.R` | `split_by_cluster()`, `tracking_matrix_to_long()` | Split cells by cluster into separate CSVs |
| `differential_abundance.R` | `build_abundance_matrix()`, `run_differential_abundance()` | GLMM-based differential abundance testing |
| `cluster_enrichment.R` | `run_cluster_enrichment()` | Fisher / mixed-effects cluster enrichment |
| `marker_statistics.R` | `run_marker_statistics()` | Per-marker lmer / Wilcoxon / Kruskal-Wallis tests |
| `differential_expression.R` | `run_differential_expression()` | Per-marker differential expression (lmer/limma/wilcox) |
| `effect_sizes.R` | `compute_effect_sizes()`, `run_effect_sizes()` | Cohen's d, Cliff's delta, Hedges' g |
| `correlation_analysis.R` | `run_correlation_analysis()` | Spearman/Pearson marker correlations |
| `diffcyt_analysis.R` | `run_diffcyt()`, `run_diffcyt_and_save()` | diffcyt DA/DS framework |
| `survival_analysis.R` | `run_survival_analysis()` | Cox PH regression and Kaplan-Meier analysis |
| `expression_matrices.R` | `create_expression_matrices()`, `combine_expression_matrices()` | Grouped expression matrix generation |
| `abundance_plotting.R` | `build_tracking_matrix()`, `plot_abundance()` | Abundance barplots and heatmaps |
| `cell_ratio.R` | `calculate_cell_ratios()` | Cell type ratio calculations |

---

## cluster_splitting.R

Splits combined expression data into per-cluster CSV files based on binary gating columns. Produces a tracking matrix counting cells per cluster per sample.

### `split_by_cluster()`

```r
split_by_cluster(input_dir, output_dir,
                 cluster_regex = "Clustering annotated",
                 markers_of_interest = NULL,
                 min_cells_per_cluster = 10,
                 chunk_size = 10000,
                 overwrite_existing = TRUE)
```

Reads CSVs with binary cluster indicator columns, splits cells by population, and writes per-cluster output files. Also produces a tracking matrix CSV.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with CSVs containing expression + cluster columns |
| `output_dir` | character | *required* | Directory for per-cluster output CSVs |
| `cluster_regex` | character | `"Clustering annotated"` | Regex to identify cluster columns |
| `markers_of_interest` | character vector | `NULL` | Markers to retain; `NULL` prompts interactive selection |
| `min_cells_per_cluster` | integer | `10` | Minimum cells to output a cluster |
| `chunk_size` | integer | `10000` | Rows per chunk for memory efficiency |
| `overwrite_existing` | logical | `TRUE` | Overwrite existing output files |

**Returns:**

```r
list(
  tracking_matrix  = data.frame(),  # Sample x cluster cell counts
  cluster_summary  = data.frame(),
  output_files     = character(),
  markers_used     = character()
)
```

**Input from:** `clusters_to_gating_csv()` (clustering/leiden.R), `load_csv_gating()` (loading/load_csv_gating.R).
**Output to:** `tracking_matrix_to_long()` (this file), `run_cluster_enrichment()`, `run_differential_abundance()`, `create_expression_matrices()`.

**CLI Usage:**

```bash
Rscript cluster_splitting.R <input_dir> <output_dir>
```

### `tracking_matrix_to_long()`

```r
tracking_matrix_to_long(tracking_matrix, id_col = "file",
                        names_to = "cluster", values_to = "count")
```

Pivots a wide tracking matrix into long format with one row per sample-cluster combination.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tracking_matrix` | data.frame | *required* | Wide tracking matrix from `split_by_cluster()` |
| `id_col` | character | `"file"` | Column identifying each sample |
| `names_to` | character | `"cluster"` | Name for the cluster column |
| `values_to` | character | `"count"` | Name for the count column |

**Returns:** `data.frame` in long format with columns `file`, `cluster`, `count`.

**Dependencies:** `io_utils`, `gui_utils`, `cluster_utils`, `memory_utils`

---

## differential_abundance.R

Tests for differential cell population abundance between experimental conditions using a negative binomial GLMM.

### `build_abundance_matrix()`

```r
build_abundance_matrix(data_dir, conditions_csv = NULL,
                       cluster_regex = "Clustering annotated",
                       condition_levels = NULL,
                       min_cluster_size = 10,
                       max_cells_per_file = 100000)
```

Reads CSV files with binary cluster indicator columns, counts cells per cluster per sample, and returns a data.frame ready for statistical testing.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_dir` | character | *required* | Directory with per-cell CSV files |
| `conditions_csv` | character | `NULL` | Path to CSV mapping file_name to condition |
| `cluster_regex` | character | `"Clustering annotated"` | Regex for cluster columns |
| `condition_levels` | character vector | `NULL` | Condition labels |
| `min_cluster_size` | integer | `10` | Minimum total cells to retain a cluster |
| `max_cells_per_file` | integer | `100000` | Max cells to load per file |

**Returns:** `list(abundance_df, all_clusters, valid_clusters)`.

**Input from:** `clusters_to_gating_csv()`, `load_csv_gating()`.
**Output to:** `run_differential_abundance()`, `run_cluster_enrichment()`, `plot_abundance()`.

### `run_differential_abundance()`

```r
run_differential_abundance(data_dir, conditions_csv, output_dir,
                           cluster_regex = "Clustering annotated",
                           condition_levels = c("pCR", "npCR"),
                           markers_of_interest = NULL,
                           min_cluster_size = 10,
                           p_adjust_method = "BH",
                           max_cells_per_file = 100000,
                           seed = 123)
```

Fits a negative binomial GLMM (`count ~ condition + offset(log_total) + (1 | sample_id)`) per cluster. Falls back to Poisson GLMM if NB fails.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_dir` | character | *required* | Directory with per-sample CSVs |
| `conditions_csv` | character | *required* | Path to conditions mapping file |
| `output_dir` | character | *required* | Output directory |
| `cluster_regex` | character | `"Clustering annotated"` | Regex for cluster columns |
| `condition_levels` | character vector | `c("pCR", "npCR")` | Condition factor levels |
| `min_cluster_size` | integer | `10` | Minimum cells per cluster |
| `p_adjust_method` | character | `"BH"` | Multiple testing correction |
| `max_cells_per_file` | integer | `100000` | Max cells per file |
| `seed` | integer | `123` | Random seed |

**Output Files:** `glmm_results.csv`, `volcano_plot.png`.

**Output CSV columns:** `cluster`, `estimate`, `std_error`, `z_value`, `p_value`, `log2fc`, `odds_ratio`, `method`, `p_adjusted`, `significant`.

**CLI Usage:**

```bash
Rscript differential_abundance.R <data_dir> <conditions_csv> <output_dir>
```

**Dependencies:** `io_utils`, `cluster_utils`, `memory_utils`, `stats_utils`

---

## cluster_enrichment.R

Tests whether cell clusters are enriched in one condition versus another using Fisher's exact test and optional mixed-effects logistic regression.

### `run_cluster_enrichment()`

```r
run_cluster_enrichment(data_dir, conditions_csv, output_dir,
                       cluster_regex = "Clustering annotated",
                       condition_levels = c("pCR", "npCR"),
                       stat_method = "mixed_effects",
                       min_cluster_frac = 0.001,
                       p_adjust_method = "BH",
                       resume = FALSE, seed = 123)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_dir` | character | *required* | Directory with cluster data |
| `conditions_csv` | character | *required* | Conditions mapping file |
| `output_dir` | character | *required* | Output directory |
| `stat_method` | character | `"mixed_effects"` | `"mixed_effects"`, `"pooled"`, or `"fisher"` |
| `resume` | logical | `FALSE` | Resume from `intermediate_results.rds` |

**Methods:**

| Method | Description |
|--------|-------------|
| `"pooled"` / `"fisher"` | Pool all cells, build 2x2 contingency table, run Fisher's exact test |
| `"mixed_effects"` | Mixed-effects logistic regression with sample as random effect |

**Output File -- `enrichment_results.csv`:**

| Column | Description |
|--------|-------------|
| `in_A` / `not_in_A` | Cells in/out of cluster in condition A |
| `odds_ratio` | Fisher's odds ratio |
| `log2_enrichment` | log2(proportion_A / proportion_B) |
| `p_adjusted` | BH-adjusted p-value |

**Example:**

```r
source("cluster_enrichment.R")
results <- run_cluster_enrichment(
  data_dir         = "/data/csvs_with_gating",
  conditions_csv   = "/data/conditions.csv",
  output_dir       = "/output/enrichment",
  condition_levels = c("Control", "Treatment"),
  stat_method      = "pooled"
)
```

---

## marker_statistics.R

Computes per-marker statistics within each population using linear mixed-effects models, Wilcoxon, or Kruskal-Wallis tests.

### `run_marker_statistics()`

```r
run_marker_statistics(input_folder, output_folder,
                      filename_regex = "^(.+)_(.+)\\.csv$",
                      markers = NULL, test_type = "lmer",
                      conditions_csv = NULL, n_workers = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_folder` | character | *required* | Directory with per-population CSVs |
| `output_folder` | character | *required* | Output directory |
| `filename_regex` | character | `"^(.+)_(.+)\\.csv$"` | Regex to extract sample + population |
| `markers` | character vector | `NULL` | Markers to analyze; `NULL` = interactive |
| `test_type` | character | `"lmer"` | `"lmer"`, `"wilcox"`, `"kruskal"`, or `"both"` |
| `conditions_csv` | character | `NULL` | Conditions mapping CSV (required for wilcox/kruskal) |
| `n_workers` | integer | `NULL` | Parallel workers; `NULL` = `detectCores() - 1` |

**Statistical Models:**

| Test | Model | Package |
|------|-------|---------|
| `"lmer"` | `marker ~ 1 + (1 \| sample_id)` | lme4 |
| `"wilcox"` | Wilcoxon rank-sum between conditions | stats |
| `"kruskal"` | Kruskal-Wallis across conditions | stats |
| `"both"` | lmer + wilcox combined | lme4 + stats |

**Output -- `marker_statistics.csv`:**

| Column | Description |
|--------|-------------|
| `population` | Cell population name |
| `marker` | Marker name |
| `mean` | Mean expression |
| `between_sample_var` | Variance of sample random intercepts |
| `within_sample_var` | Residual (cell-to-cell) variance |
| `icc` | Intraclass correlation coefficient |
| `method` | `"lmer"` or `"simple_mean"` (fallback) |

**Dependencies:** `io_utils`, `gui_utils`, `stats_utils`. Requires `lme4`, `lmerTest`.

---

## differential_expression.R

Per-marker differential expression between conditions using mixed-effects models, limma pseudobulk, or Wilcoxon rank-sum tests.

### `run_differential_expression()`

```r
run_differential_expression(input_dir, conditions_csv, output_dir,
                             cluster_regex = "Clustering annotated",
                             markers = NULL, test_type = "lmer",
                             p_adjust_method = "BH", seed = 123)
```

Tests each marker within each cluster for differential expression between two conditions.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with per-cluster CSVs |
| `conditions_csv` | character | *required* | CSV mapping `sample_id` to `condition` |
| `output_dir` | character | *required* | Output directory |
| `markers` | character vector | `NULL` | Markers to test; `NULL` = auto-detect |
| `test_type` | character | `"lmer"` | `"lmer"` (mixed-effects), `"limma"` (pseudobulk), `"wilcox"` (Wilcoxon) |
| `p_adjust_method` | character | `"BH"` | Multiple testing correction |
| `seed` | integer | `123` | Random seed |

**Returns:** `data.frame(cluster, marker, log2fc, pvalue, padj, test)`.

**Output File:** `differential_expression_results.csv`.

**Input from:** `clusters_to_gating_csv()`, `load_csv_gating()`.
**Output to:** `plot_volcano()`, `plot_heatmap()`.

**CLI Usage:**

```bash
Rscript differential_expression.R <input_dir> <conditions_csv> <output_dir> [test_type]
```

---

## effect_sizes.R

Computes Cohen's d, Cliff's delta, and Hedges' g for marker expression differences between two conditions.

### `compute_effect_sizes()`

```r
compute_effect_sizes(expr_data, condition_col, value_cols,
                     methods = c("cohens_d", "cliffs_delta"))
```

Computes effect sizes with confidence intervals and standard interpretations.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_data` | data.table/data.frame | *required* | Data with marker + condition columns |
| `condition_col` | character | *required* | Column name for condition grouping |
| `value_cols` | character vector | *required* | Marker column names |
| `methods` | character vector | `c("cohens_d", "cliffs_delta")` | Methods: `"cohens_d"`, `"cliffs_delta"`, `"hedges_g"` |

**Returns:** `data.frame(marker, effect_size, ci_lower, ci_upper, method, interpretation)`.

**Interpretation thresholds for Cohen's d / Hedges' g:**

| |d| | Interpretation |
|------|----------------|
| < 0.2 | negligible |
| 0.2--0.5 | small |
| 0.5--0.8 | medium |
| >= 0.8 | large |

**Interpretation thresholds for Cliff's delta:**

| |delta| | Interpretation |
|---------|----------------|
| < 0.147 | negligible |
| 0.147--0.33 | small |
| 0.33--0.474 | medium |
| >= 0.474 | large |

### `run_effect_sizes()`

```r
run_effect_sizes(input_csv, output_dir, condition_col = "condition",
                 value_cols = NULL, methods = c("cohens_d", "cliffs_delta"))
```

Convenience wrapper that reads from CSV and saves results to `effect_sizes.csv`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_csv` | character | *required* | Path to expression CSV |
| `output_dir` | character | *required* | Output directory |
| `condition_col` | character | `"condition"` | Condition column name |
| `value_cols` | character vector | `NULL` | Marker columns; `NULL` = auto-detect |
| `methods` | character vector | `c("cohens_d", "cliffs_delta")` | Effect size methods |

**Returns:** `data.frame` with effect sizes.

**CLI Usage:**

```bash
Rscript effect_sizes.R <input_csv> <output_dir> [condition_col]
```

---

## correlation_analysis.R

Computes pairwise marker correlations with p-values, optionally stratified by cluster.

### `run_correlation_analysis()`

```r
run_correlation_analysis(expr_matrix, method = "spearman",
                          output_dir = NULL, cluster_labels = NULL,
                          p_adjust_method = "BH")
```

Computes correlation matrices, identifies significant pairs, and produces heatmap plots.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | numeric matrix | *required* | Cells x markers matrix |
| `method` | character | `"spearman"` | `"spearman"` or `"pearson"` |
| `output_dir` | character | `NULL` | Output directory (for plots + CSV); `NULL` = no saving |
| `cluster_labels` | vector | `NULL` | Optional cluster labels for per-cluster correlations |
| `p_adjust_method` | character | `"BH"` | P-value adjustment method |

**Returns:**

```r
list(
  cor_matrix        = matrix,       # Correlation matrix
  p_matrix          = matrix,       # P-value matrix
  significant_pairs = data.frame,   # Pairs with padj < 0.05
  cluster_results   = list          # Per-cluster results (when cluster_labels provided)
)
```

**Output Files:** `{prefix}_correlation_matrix.csv`, `{prefix}_pvalue_matrix.csv`, `{prefix}_significant_pairs.csv`, `{prefix}_correlation_heatmap.png`.

**Input from:** `build_combined_expr_matrix()`, `combine_expression_matrices()`.
**Output to:** `plot_heatmap()` (visualization/heatmaps.R).

**CLI Usage:**

```bash
Rscript correlation_analysis.R <expression_csv> <output_dir> [method]
```

---

## diffcyt_analysis.R

Wraps the diffcyt Bioconductor package for differential abundance (DA) and differential state (DS) analysis.

### `run_diffcyt()`

```r
run_diffcyt(flow_frames, experiment_info, markers_of_interest,
            analysis_type = "DA", method_DA = "diffcyt-DA-edgeR",
            method_DS = "diffcyt-DS-limma", n_clusters = 20, seed = 42)
```

Runs the full diffcyt pipeline: prepare data, transform, cluster with FlowSOM, compute summary statistics, and run statistical tests.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | List of flowFrame objects |
| `experiment_info` | data.frame | *required* | Must have `sample_id`, `group_id` columns (+ optional `patient_id`) |
| `markers_of_interest` | character vector | *required* | Markers for clustering and testing |
| `analysis_type` | character | `"DA"` | `"DA"` (differential abundance) or `"DS"` (differential state) |
| `method_DA` | character | `"diffcyt-DA-edgeR"` | DA method: `"diffcyt-DA-edgeR"`, `"diffcyt-DA-voom"`, `"diffcyt-DA-GLMM"` |
| `method_DS` | character | `"diffcyt-DS-limma"` | DS method: `"diffcyt-DS-limma"`, `"diffcyt-DS-LMM"` |
| `n_clusters` | integer | `20` | Number of FlowSOM meta-clusters |
| `seed` | integer | `42` | Random seed |

**Returns:**

```r
list(
  results              = <diffcyt result>,
  top_table            = data.frame,
  cluster_assignments  = data.frame,
  method_used          = character
)
```

**Input from:** flowFrame lists from `load_fcs_raw()`, `build_combined_expr_matrix()`.
**Output to:** Standalone results CSV.

### `run_diffcyt_and_save()`

```r
run_diffcyt_and_save(flow_frames, experiment_info, markers_of_interest,
                      output_dir, analysis_type = "DA",
                      method_DA = "diffcyt-DA-edgeR",
                      method_DS = "diffcyt-DS-limma",
                      n_clusters = 20, seed = 42)
```

Convenience wrapper that runs `run_diffcyt()` and saves results to `output_dir`.

**Output Files:** `diffcyt_{da|ds}_results.csv`, `diffcyt_cluster_assignments.csv`.

**Example:**

```r
source("diffcyt_analysis.R")

experiment_info <- data.frame(
  sample_id  = c("S1", "S2", "S3", "S4"),
  group_id   = c("Control", "Control", "Treatment", "Treatment"),
  patient_id = c("P1", "P2", "P3", "P4")
)

result <- run_diffcyt_and_save(
  flow_frames         = ff_list,
  experiment_info     = experiment_info,
  markers_of_interest = c("CD3", "CD4", "CD8"),
  output_dir          = "/output/diffcyt",
  analysis_type       = "DA"
)
```

---

## survival_analysis.R

Performs Cox proportional hazards regression and Kaplan-Meier survival analysis using cluster abundance proportions as predictors.

### `run_survival_analysis()`

```r
run_survival_analysis(abundance_matrix, clinical_csv, output_dir,
                       method = "both", split_method = "median",
                       p_adjust_method = "BH")
```

Tests each cluster's abundance as a predictor of survival. Supports median and maximally-selected-rank-statistics (optimal) cutpoints.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `abundance_matrix` | data.frame | *required* | Sample rows x cluster columns (proportions). First column = sample ID. |
| `clinical_csv` | character | *required* | Path to CSV with `sample_id`, `time`, `event` columns |
| `output_dir` | character | *required* | Output directory |
| `method` | character | `"both"` | `"cox"`, `"km"` (Kaplan-Meier), or `"both"` |
| `split_method` | character | `"median"` | `"median"` or `"optimal"` (maxstat via survminer) |
| `p_adjust_method` | character | `"BH"` | P-value adjustment |

**Returns:** `data.frame(cluster, hr, ci_lower, ci_upper, pvalue, padj, km_pvalue, cutpoint, n_high, n_low, cox_concordance)`.

**Output Files:** `survival_results.csv`, `forest_plot_hr.png`, `km_plots/km_{cluster}.png` (one per cluster).

**Input from:** `build_tracking_matrix()` (abundance_plotting.R), clinical CSV.
**Output to:** Standalone Kaplan-Meier plots and Cox results.

**CLI Usage:**

```bash
Rscript survival_analysis.R <abundance_csv> <clinical_csv> <output_dir> [method] [split_method]
```

**Packages:** `survival`, `survminer`

---

## expression_matrices.R

Generates grouped expression matrices from per-cluster data.

### `create_expression_matrices()`

```r
create_expression_matrices(input_dir, output_dir, markers = NULL,
                           group_by_pattern = NULL)
```

Groups cluster-split CSV files by cluster name and writes one expression matrix per group.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with cluster-split CSVs |
| `output_dir` | character | *required* | Output directory |
| `markers` | character vector | `NULL` | Markers to include; `NULL` auto-detects |
| `group_by_pattern` | character | `NULL` | Regex for cluster group extraction |

**Returns:** List of group-level expression matrices.

**Output Files:** `exprmat_{group}.csv`, `group_matrices_summary.csv`.

**Input from:** `split_by_cluster()`.
**Output to:** `plot_heatmap()`, `plot_cluster_heatmap()`.

### `combine_expression_matrices()`

```r
combine_expression_matrices(expr_data_list, markers, max_cells = NULL, seed = 42)
```

Combines a named list of data.tables into a single expression matrix, filtering to specified markers and optionally subsampling.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_data_list` | named list | *required* | List of data.tables (one per file) |
| `markers` | character vector | *required* | Marker columns to retain |
| `max_cells` | integer | `NULL` | Maximum cells in output; `NULL` keeps all |
| `seed` | integer | `42` | Random seed for subsampling |

**Returns:** Numeric matrix (cells x markers), or `NULL` if no matching markers.

**Example:**

```r
source("expression_matrices.R")

expr_data <- load_csv_expression(input_dir = "/data/csv_expression")
combined <- combine_expression_matrices(
  expr_data_list = expr_data,
  markers        = c("CD3", "CD4", "CD8", "CD19"),
  max_cells      = 50000
)
dim(combined)
# [1] 50000     4
```

---

## abundance_plotting.R

Generates stacked, proportional, and heatmap abundance plots from cluster tracking matrices.

### `build_tracking_matrix()`

```r
build_tracking_matrix(input_dir, cluster_regex, exclude_patterns, output_dir)
```

Scans CSVs in `input_dir`, counts cells per cluster per sample, saves `cluster_tracking_matrix.csv`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with per-sample CSVs |
| `cluster_regex` | character | *required* | Regex to identify cluster columns |
| `exclude_patterns` | character vector | *required* | Regex patterns to exclude (e.g., sub-clusters) |
| `output_dir` | character | *required* | Output directory |

**Returns:** `data.frame` with columns `file`, plus one column per cluster (cell counts).

**Input from:** `clusters_to_gating_csv()`, `load_csv_gating()`.
**Output to:** `plot_abundance()`.

### `plot_abundance()`

```r
plot_abundance(tracking_matrix_path = NULL, output_dir,
               tracking_matrix = NULL, metadata_file = NULL)
```

Generates multiple abundance plots from a tracking matrix.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tracking_matrix_path` | character | `NULL` | Path to tracking matrix CSV |
| `output_dir` | character | *required* | Output directory |
| `tracking_matrix` | data.frame | `NULL` | Pre-loaded tracking matrix |
| `metadata_file` | character | `NULL` | Optional metadata sheet for grouped plots |

**Generated Plots:**

| File | Description |
|------|-------------|
| `abundance_stacked.png` | Stacked barplot sorted by total cells per sample |
| `abundance_proportional.png` | Proportional barplot with Ward.D2 dendrogram |
| `total_cells_per_cluster.png` | Bar chart of total cells per cluster |
| `cluster_proportion_heatmap.png` | Biclustered heatmap (inferno palette) |

**Example:**

```r
source("abundance_plotting.R")

tm <- build_tracking_matrix(
  input_dir        = "/data/csv_gating",
  cluster_regex    = "Clustering annotated",
  exclude_patterns = c("Clustering raw"),
  output_dir       = "/output/abundance"
)

plot_abundance(output_dir = "/output/abundance", tracking_matrix = tm)
```

**Dependencies:** `io_utils`, `plot_utils`

---

## cell_ratio.R

Calculates ratios between specified cluster populations per sample.

### `calculate_cell_ratios()`

```r
calculate_cell_ratios(abundance_csv, numerator_clusters,
                      denominator_clusters, output_path = NULL)
```

Reads a cluster tracking matrix and computes per-sample ratios with log2 transformation.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `abundance_csv` | character | *required* | Path to CSV with cluster counts per sample |
| `numerator_clusters` | character vector | *required* | Cluster names for numerator |
| `denominator_clusters` | character vector | *required* | Cluster names for denominator |
| `output_path` | character | `NULL` | Output CSV path |

**Returns:** `data.frame(sample, numerator_sum, denominator_sum, ratio, log2_ratio)`.

**Example:**

```r
source("cell_ratio.R")
ratios <- calculate_cell_ratios(
  abundance_csv        = "/output/cluster_tracking_matrix.csv",
  numerator_clusters   = c("CD4+ T cells"),
  denominator_clusters = c("CD8+ T cells"),
  output_path          = "/output/cd4_cd8_ratio.csv"
)
```

**Dependencies:** `io_utils`
