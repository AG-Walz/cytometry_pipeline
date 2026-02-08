[< Back to Documentation Index](../README.md)

# src/analysis/ -- Statistical Analysis

7 scripts for cluster splitting, differential abundance testing, enrichment analysis, marker statistics, expression matrix generation, abundance plotting, and cell ratio calculation.

---

## cluster_splitting.R

Splits combined expression data into per-cluster CSV files for downstream analysis.

### Main Function

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
| `input_dir` | character | *required* | Directory with input CSVs containing expression + cluster columns |
| `output_dir` | character | *required* | Directory for per-cluster output CSVs |
| `cluster_regex` | character | `"Clustering annotated"` | Regex to identify cluster label columns |
| `markers_of_interest` | character vector | `NULL` | Markers to retain; `NULL` keeps all expression columns |
| `min_cells_per_cluster` | integer | `10` | Minimum cells required to output a cluster file |
| `chunk_size` | integer | `10000` | Rows per chunk for memory-efficient processing |
| `overwrite_existing` | logical | `TRUE` | Whether to overwrite existing output files |

### Output

```r
list(
  tracking_matrix  = data.frame(), # Sample-by-cluster cell count matrix
  cluster_summary  = data.frame(), # Per-cluster statistics
  output_files     = character(),  # Paths to generated CSV files
  markers_used     = character()   # Markers included in output
)
```

Generated files:
- Per-cluster CSVs in `output_dir`
- `cluster_tracking_matrix.csv` -- sample-by-cluster cell counts

### CLI Usage

```bash
Rscript cluster_splitting.R <input_dir> <output_dir>
```

### src/R_shared Dependencies

`io_utils`, `gui_utils`, `cluster_utils`, `memory_utils`

---

## differential_abundance.R

Tests for differential cell population abundance between experimental conditions using negative binomial GLMM.

### Main Function

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
| `data_dir` | character | *required* | Directory containing per-cluster CSVs |
| `conditions_csv` | character | *required* | Path to conditions file (see format below) |
| `output_dir` | character | *required* | Directory for results |
| `cluster_regex` | character | `"Clustering annotated"` | Regex to identify cluster columns |
| `condition_levels` | character vector | `c("pCR", "npCR")` | Factor levels for the condition variable |
| `markers_of_interest` | character vector | `NULL` | Markers to include |
| `min_cluster_size` | integer | `10` | Minimum cells per cluster to test |
| `p_adjust_method` | character | `"BH"` | Multiple testing correction method |
| `max_cells_per_file` | integer | `100000` | Maximum cells to load per file |
| `seed` | integer | `123` | Random seed |

### Conditions CSV Format

Semicolon-separated file with two columns:

```
file_name;condition
sample_001.csv;pCR
sample_002.csv;npCR
sample_003.csv;pCR
```

### Statistical Method

Negative binomial GLMM fitted via `lme4`:

```
count ~ condition + offset(log_total) + (1 | sample_id)
```

### Output

Results `data.frame`:

| Column | Description |
|--------|-------------|
| `cluster` | Cluster identifier |
| `estimate` | Model coefficient |
| `std_error` | Standard error of estimate |
| `z_value` | Wald z-statistic |
| `p_value` | Raw p-value |
| `log2fc` | Log2 fold change |
| `odds_ratio` | Exponentiated estimate |
| `method` | `"NB-GLMM"` |
| `p_adjusted` | Adjusted p-value |
| `significant` | Logical, `TRUE` if `p_adjusted < 0.05` |

Generated files:
- `glmm_results.csv`
- `volcano_plot.png`

### CLI Usage

```bash
Rscript differential_abundance.R <data_dir> <conditions_csv> <output_dir>
```

### src/R_shared Dependencies

`io_utils`, `cluster_utils`, `memory_utils`, `stats_utils` (`fit_nb_glmm`)

---

## cluster_enrichment.R

Tests whether cell clusters are enriched in one condition versus another using Fisher's exact test and optional mixed-effects logistic regression.

### Main Function

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
| `cluster_regex` | character | `"Clustering annotated"` | Regex for cluster columns |
| `condition_levels` | character vector | `c("pCR", "npCR")` | Condition factor levels |
| `stat_method` | character | `"mixed_effects"` | `"mixed_effects"` or `"fisher"` |
| `min_cluster_frac` | numeric | `0.001` | Minimum cluster fraction to test |
| `p_adjust_method` | character | `"BH"` | P-value adjustment method |
| `resume` | logical | `FALSE` | Resume from checkpoint (`intermediate_results.rds`) |
| `seed` | integer | `123` | Random seed |

### Methods

- **Fisher's exact test**: 2x2 contingency table (cluster membership x condition)
- **Mixed-effects logistic regression**: Accounts for sample-level random effects
- **Log2 enrichment**: `log2(fraction_condition1 / fraction_condition2)`

### Output

- `enrichment_results.csv`
- `intermediate_results.rds` (checkpointing for long-running analyses)

### src/R_shared Dependencies

`io_utils`, `cluster_utils`, `memory_utils`, `stats_utils` (`calculate_enrichment`, `do_mixed_effects`)

---

## marker_statistics.R

Computes per-marker statistics within each population, including between- and within-sample variance via linear mixed-effects models.

### Main Function

```r
run_marker_statistics(input_folder, output_folder,
                      filename_regex = "^(.+)_(.+)\\.csv$",
                      markers = NULL, test_type = "lmer",
                      n_workers = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_folder` | character | *required* | Directory with per-cluster CSVs |
| `output_folder` | character | *required* | Output directory |
| `filename_regex` | character | `"^(.+)_(.+)\\.csv$"` | Pattern to extract sample and cluster from filenames |
| `markers` | character vector | `NULL` | Markers to analyze; `NULL` uses all |
| `test_type` | character | `"lmer"` | Statistical method (currently `"lmer"`) |
| `n_workers` | integer | `NULL` | Number of parallel workers; `NULL` for sequential |

### Statistical Method

Linear mixed-effects model via `lme4::lmer`:

```
marker ~ 1 + (1 | sample_id)
```

Extracts the intraclass correlation coefficient (ICC) to quantify between-sample vs within-sample variance.

### Output

`marker_statistics.csv` with columns:

| Column | Description |
|--------|-------------|
| `population` | Cluster/population name |
| `marker` | Marker name |
| `mean` | Grand mean expression |
| `between_sample_var` | Between-sample variance component |
| `within_sample_var` | Within-sample (residual) variance |
| `icc` | Intraclass correlation coefficient |
| `n_cells` | Total number of cells |
| `n_samples` | Number of samples |
| `method` | `"lmer"` |

### Packages

`lme4`, `lmerTest`, `future`

### src/R_shared Dependencies

`io_utils`, `gui_utils`, `stats_utils`

---

## expression_matrices.R

Generates grouped expression matrices from per-cluster data.

### Main Function

```r
create_expression_matrices(input_dir, output_dir, markers = NULL,
                           group_by_pattern = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | character | *required* | Directory with per-cluster CSVs |
| `output_dir` | character | *required* | Output directory |
| `markers` | character vector | `NULL` | Markers to include; `NULL` uses all |
| `group_by_pattern` | character | `NULL` | Regex pattern for grouping clusters |

### Output

- `exprmat_<group>.csv` -- one expression matrix per cluster group
- `group_matrices_summary.csv` -- summary of generated matrices

### src/R_shared Dependencies

`io_utils`, `cluster_utils`, `stats_utils`

---

## abundance_plotting.R

Generates stacked, grouped, or proportional barplots from cluster abundance data.

### Main Function

```r
plot_abundance(tracking_matrix_path, conditions_csv = NULL,
               output_dir, plot_type = "stacked")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tracking_matrix_path` | character | *required* | Path to `cluster_tracking_matrix.csv` |
| `conditions_csv` | character | `NULL` | Path to conditions file for condition-level grouping |
| `output_dir` | character | *required* | Output directory for plots |
| `plot_type` | character | `"stacked"` | Plot style: `"stacked"`, `"grouped"`, or `"proportional"` |

### Output

`abundance_*.png` files in `output_dir`.

### src/R_shared Dependencies

`io_utils`, `plot_utils`

---

## cell_ratio.R

Calculates ratios between specified cluster populations per sample.

### Main Function

```r
calculate_cell_ratios(abundance_csv, numerator_clusters,
                      denominator_clusters, output_path = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `abundance_csv` | character | *required* | Path to abundance/tracking matrix CSV |
| `numerator_clusters` | character vector | *required* | Cluster names for the numerator |
| `denominator_clusters` | character vector | *required* | Cluster names for the denominator |
| `output_path` | character | `NULL` | Path to save results CSV |

### Output

`data.frame`:

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier |
| `numerator_sum` | Total cells in numerator clusters |
| `denominator_sum` | Total cells in denominator clusters |
| `ratio` | `numerator_sum / denominator_sum` |
| `log2_ratio` | `log2(ratio)` |

### Example

```r
# Calculate CD4/CD8 ratio
ratios <- calculate_cell_ratios(
  abundance_csv        = "output/cluster_tracking_matrix.csv",
  numerator_clusters   = c("CD4+ T cells"),
  denominator_clusters = c("CD8+ T cells"),
  output_path          = "output/cd4_cd8_ratio.csv"
)
```

### src/R_shared Dependencies

`io_utils`
