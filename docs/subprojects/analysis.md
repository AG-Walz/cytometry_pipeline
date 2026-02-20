[< Back to Documentation Index](../README.md)

# src/analysis/ — Statistical Analysis

Seven scripts for cluster splitting, differential abundance testing, enrichment analysis, marker statistics, expression matrix generation, abundance plotting, and cell ratio calculation. These scripts are the statistical backbone of Pipelines A and C.

---

## cluster_splitting.R

Splits combined expression data into per-cluster CSV files based on binary gating columns. Produces a tracking matrix counting cells per cluster per sample.

### Function

```r
split_by_cluster(input_dir, output_dir,
                 cluster_regex = "Clustering annotated",
                 markers_of_interest = NULL,
                 min_cells_per_cluster = 10,
                 chunk_size = 10000,
                 overwrite_existing = TRUE)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with CSVs containing expression + cluster columns |
| `output_dir` | character | *required* | Directory for per-cluster output CSVs |
| `cluster_regex` | character | `"Clustering annotated"` | Regex to identify cluster columns |
| `markers_of_interest` | character vector | `NULL` | Markers to retain; `NULL` keeps all |
| `min_cells_per_cluster` | integer | `10` | Minimum cells to output a cluster |
| `chunk_size` | integer | `10000` | Rows per chunk for memory efficiency |

### Returns

```r
list(
  tracking_matrix  = data.frame(),  # Sample × cluster cell counts
  cluster_summary  = data.frame(),
  output_files     = character(),
  markers_used     = character()
)
```

### Example Input

CSVs with expression columns and binary cluster columns:

```
CD3;CD4;CD8;CD4+ (Clustering annotated);CD8+ (Clustering annotated)
2.46;4.01;2.15;1;0
3.89;2.12;4.56;0;1
0.98;1.56;2.78;0;1
```

### Example Output

**`cluster_tracking_matrix.csv`:**

```
file;CD4+;CD8+
sample_001;20;15
sample_002;20;15
sample_003;20;15
```

**Per-cluster files:**

```
output_dir/
├── cluster_tracking_matrix.csv
├── sample_001_CD4+.csv         # Expression data for CD4+ cells from sample_001
├── sample_001_CD8+.csv         # Expression data for CD8+ cells from sample_001
├── sample_002_CD4+.csv
└── ...
```

### CLI Usage

```bash
Rscript cluster_splitting.R <input_dir> <output_dir>
```

### Dependencies

`io_utils`, `gui_utils`, `cluster_utils`, `memory_utils`

### tracking_matrix_to_long

Pivots a wide tracking matrix (samples as rows, clusters as columns) into long format with one row per sample-cluster combination.

```r
tracking_matrix_to_long(tracking_matrix, id_col = "file",
                        names_to = "cluster", values_to = "count")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tracking_matrix` | data.frame | *required* | Wide tracking matrix from `split_by_cluster()` |
| `id_col` | character | `"file"` | Column identifying each sample |
| `names_to` | character | `"cluster"` | Name for the cluster column in long output |
| `values_to` | character | `"count"` | Name for the count column in long output |

**Returns:** `data.frame` in long format with columns `file`, `cluster`, `count`.

#### Example

```r
source("cluster_splitting.R")

result <- split_by_cluster(input_dir = "/data/csvs", output_dir = "/output/split")
long_df <- tracking_matrix_to_long(result$tracking_matrix)
head(long_df)
#        file cluster count
# 1 sample_001    CD4+    20
# 2 sample_001    CD8+    15
# 3 sample_002    CD4+    22
```

---

## differential_abundance.R

Tests for differential cell population abundance between experimental conditions using a negative binomial GLMM.

### Function

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

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_dir` | character | *required* | Directory with per-sample CSVs |
| `conditions_csv` | character | *required* | Path to conditions mapping file |
| `output_dir` | character | *required* | Output directory |
| `cluster_regex` | character | `"Clustering annotated"` | Regex for cluster columns |
| `condition_levels` | character vector | `c("pCR", "npCR")` | Condition factor levels |
| `min_cluster_size` | integer | `10` | Minimum cells per cluster to test |
| `p_adjust_method` | character | `"BH"` | Multiple testing correction |
| `seed` | integer | `123` | Random seed |

### Conditions CSV Format

```
file_name;condition
sample_001;Control
sample_002;Control
sample_003;Treatment
```

### Statistical Model

Negative binomial GLMM fitted via `lme4`:

```
count ~ condition + offset(log_total) + (1 | sample_id)
```

Falls back to Poisson GLMM if negative binomial fails to converge.

### Output

**`glmm_results.csv`:**

```
cluster;estimate;std_error;z_value;p_value;log2fc;odds_ratio;method;p_adjusted;significant
Cluster_1;-0.23;0.15;-1.53;0.126;-0.33;0.79;glmerMod;0.126;FALSE
Cluster_2;0.87;0.21;4.14;0.00003;1.26;2.39;glmerMod;0.0001;TRUE
```

| Column | Description |
|--------|-------------|
| `estimate` | Model coefficient for condition effect |
| `log2fc` | Log2 fold change between conditions |
| `odds_ratio` | Exponentiated estimate |
| `p_adjusted` | BH-adjusted p-value |
| `significant` | `TRUE` if `p_adjusted < 0.05` |

**`volcano_plot.png`:** Scatter plot of log2FC vs -log10(adjusted p-value).

### CLI Usage

```bash
Rscript differential_abundance.R <data_dir> <conditions_csv> <output_dir>
```

### Dependencies

`io_utils`, `cluster_utils`, `memory_utils`, `stats_utils` (`fit_nb_glmm`)

---

## cluster_enrichment.R

Tests whether cell clusters are enriched in one condition versus another using Fisher's exact test and optional mixed-effects logistic regression.

### Function

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

### Methods

- **`"pooled"` / `"fisher"`**: Pool all cells, build 2x2 contingency table, run Fisher's exact test
- **`"mixed_effects"`**: Mixed-effects logistic regression with sample as random effect

### Output

**`enrichment_results.csv`:**

```
cluster;in_A;not_in_A;in_B;not_in_B;odds_ratio;p_value;ci_lower;ci_upper;method;p_adjusted;significant;log2_enrichment
CD4+;40;60;20;30;1.00;1.00;0.47;2.13;fisher;1.00;FALSE;0.00
CD8+;30;70;15;35;1.00;1.00;0.45;2.27;fisher;1.00;FALSE;0.00
```

| Column | Description |
|--------|-------------|
| `in_A` / `not_in_A` | Cells in/out of cluster in condition A |
| `in_B` / `not_in_B` | Cells in/out of cluster in condition B |
| `odds_ratio` | Fisher's odds ratio |
| `log2_enrichment` | log2(proportion_A / proportion_B) |

### Example

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

### Dependencies

`io_utils`, `cluster_utils`, `memory_utils`, `stats_utils`

---

## marker_statistics.R

Computes per-marker statistics within each population using linear mixed-effects models.

### Function

```r
run_marker_statistics(input_folder, output_folder,
                      filename_regex = "^(.+)_(.+)\\.csv$",
                      markers = NULL, test_type = "lmer",
                      n_workers = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_folder` | character | *required* | Directory with per-population CSVs |
| `output_folder` | character | *required* | Output directory |
| `filename_regex` | character | `"^(.+)_(.+)\\.csv$"` | Regex to extract sample + population |
| `markers` | character vector | `NULL` | Markers to analyze |
| `test_type` | character | `"lmer"` | Statistical method |

### Statistical Model

Linear mixed-effects model via `lme4::lmer`:

```
marker ~ 1 + (1 | sample_id)
```

Decomposes variance into between-sample and within-sample components, computing the ICC.

### Output

**`marker_statistics.csv`:**

```
population;marker;mean;between_sample_var;within_sample_var;icc;n_cells;n_samples;method
CD4_T_cells;CD3;3.02;0.00;0.97;0.00;150;3;lmer
CD4_T_cells;CD4;3.08;0.003;1.10;0.003;150;3;lmer
```

| Column | Description |
|--------|-------------|
| `between_sample_var` | Variance of sample random intercepts |
| `within_sample_var` | Residual (cell-to-cell) variance |
| `icc` | Intraclass correlation coefficient |
| `method` | `"lmer"` or `"simple_mean"` (fallback) |

### Dependencies

`io_utils`, `gui_utils`, `stats_utils`. Requires `lme4`, `lmerTest`.

---

## abundance_plotting.R

Generates stacked, proportional, and heatmap abundance plots from cluster tracking matrices.

### Functions

```r
build_tracking_matrix(input_dir, cluster_regex, exclude_patterns, output_dir)
```

Scans CSVs in `input_dir`, counts cells per cluster per sample, saves `cluster_tracking_matrix.csv`.

```r
plot_abundance(tracking_matrix_path = NULL, output_dir,
               tracking_matrix = NULL, metadata_file = NULL)
```

Generates multiple abundance plots from a tracking matrix.

### Generated Plots

| File | Description |
|------|-------------|
| `abundance_stacked.png` | Stacked barplot sorted by total cells per sample |
| `abundance_proportional.png` | Proportional barplot with Ward.D2 dendrogram |
| `total_cells_per_cluster.png` | Bar chart of total cells per cluster |
| `cluster_proportion_heatmap.png` | Biclustered heatmap (inferno palette) |

If `metadata_file` is provided (Cytolution metadata sheet), additional grouped plots are generated based on metadata columns.

### Example

```r
source("abundance_plotting.R")

tm <- build_tracking_matrix(
  input_dir        = "/data/csv_gating",
  cluster_regex    = "Clustering annotated",
  exclude_patterns = c("Clustering raw"),
  output_dir       = "/output/abundance"
)

plot_abundance(
  output_dir      = "/output/abundance",
  tracking_matrix = tm
)
```

### Dependencies

`io_utils`, `plot_utils`

---

## expression_matrices.R

Generates grouped expression matrices from per-cluster data.

### Function

```r
create_expression_matrices(input_dir, output_dir, markers = NULL,
                           group_by_pattern = NULL)
```

### Output

- `exprmat_<group>.csv` — one matrix per cluster group
- `group_matrices_summary.csv` — summary of generated matrices

### combine_expression_matrices

Combines a named list of data.tables into a single expression matrix, filtering to the specified markers. Optionally subsamples to a maximum cell count.

```r
combine_expression_matrices(expr_data_list, markers,
                            max_cells = NULL, seed = 42)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_data_list` | named list | *required* | List of data.tables (one per sample/file) |
| `markers` | character vector | *required* | Marker columns to retain |
| `max_cells` | integer | `NULL` | Maximum cells in output; `NULL` keeps all |
| `seed` | integer | `42` | Random seed for subsampling |

**Returns:** Numeric matrix (cells x markers), or `NULL` if no matching markers are found.

#### Example

```r
source("expression_matrices.R")

expr_data <- load_csv_expression(input_dir = "/data/csv_expression")
combined <- combine_expression_matrices(
  expr_data_list = expr_data,
  markers        = c("CD3", "CD4", "CD8", "CD19"),
  max_cells      = 50000
)
# Combined matrix: 120000 cells x 4 markers from 10 files
# Subsampled to 50000 cells
dim(combined)
# [1] 50000     4
```

### Dependencies

`io_utils`, `cluster_utils`, `stats_utils`

---

## cell_ratio.R

Calculates ratios between specified cluster populations per sample.

### Function

```r
calculate_cell_ratios(abundance_csv, numerator_clusters,
                      denominator_clusters, output_path = NULL)
```

### Returns

```r
data.frame(sample, numerator_sum, denominator_sum, ratio, log2_ratio)
```

### Example

```r
source("cell_ratio.R")
ratios <- calculate_cell_ratios(
  abundance_csv        = "/output/cluster_tracking_matrix.csv",
  numerator_clusters   = c("CD4+ T cells"),
  denominator_clusters = c("CD8+ T cells"),
  output_path          = "/output/cd4_cd8_ratio.csv"
)
```

**Output — `cd4_cd8_ratio.csv`:**

```
sample;numerator_sum;denominator_sum;ratio;log2_ratio
sample_001;1234;567;2.18;1.12
sample_002;987;1234;0.80;-0.32
sample_003;1456;432;3.37;1.75
```

### Dependencies

`io_utils`
