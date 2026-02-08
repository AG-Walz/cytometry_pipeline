[< Back to Documentation Index](../README.md)

# src/metadata_qc/ -- Metadata & Quality Assurance

4 scripts for validating study metadata, checking FCS acquisition dates, computing inter-rater agreement, and comparing gating strategies across FlowJo workspaces.

---

## check_metadata.R

Validates and visualizes study metadata, detecting missing values and temporal patterns.

### Functions

#### check_metadata

```r
check_metadata(metadata_file, output_dir, sep = ",",
               skip_columns = NULL, only_cols = NULL,
               sort_column = NULL, sort_desc = FALSE,
               n_splits = 3,
               placeholders = c("?", "na", "n/a", "null", "-"))
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `metadata_file` | character | *required* | Path to metadata CSV file |
| `output_dir` | character | *required* | Output directory for plots and summaries |
| `sep` | character | `","` | Column separator |
| `skip_columns` | character vector | `NULL` | Columns to exclude from analysis |
| `only_cols` | character vector | `NULL` | If set, analyze only these columns |
| `sort_column` | character | `NULL` | Column to sort rows by before splitting |
| `sort_desc` | logical | `FALSE` | Sort in descending order |
| `n_splits` | integer | `3` | Number of row splits for temporal analysis |
| `placeholders` | character vector | `c("?", "na", "n/a", "null", "-")` | Strings treated as missing values |

#### plot_temporal_heatmap

```r
plot_temporal_heatmap(metadata_file, output_dir, sep = ",",
                      skip_columns = NULL, date_column = "Datum",
                      placeholders = ...)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `metadata_file` | character | *required* | Path to metadata CSV file |
| `output_dir` | character | *required* | Output directory |
| `sep` | character | `","` | Column separator |
| `skip_columns` | character vector | `NULL` | Columns to exclude |
| `date_column` | character | `"Datum"` | Name of the date column |
| `placeholders` | character vector | see `check_metadata` | Missing value placeholders |

### Features

- Auto-detects date formats
- Converts placeholder strings to `"Unassigned"`
- Generates faceted barplots showing value distributions per column
- Splits data into temporal groups for trend detection

### Output

| File | Description |
|------|-------------|
| `consensus_barplots_compact.png` | Faceted barplots for all metadata columns |
| `consensus_barplots_compact_bySplit_N.png` | Barplots per temporal split (N = 1 to `n_splits`) |
| `metadata_summary_table.csv` | Summary statistics per column |
| `temporal_feature_heatmap.png` | Heatmap of feature distributions over time |

---

## check_acquisition.R

Extracts and validates acquisition metadata from FCS file headers, with duplicate file detection.

### Functions

#### get_fcs_keyword

```r
get_fcs_keyword(fcs_path, keyword_name)
```

Returns the value of a specific FCS keyword (e.g., `"$DATE"`, `"$CYT"`).

#### explore_fcs_keywords

```r
explore_fcs_keywords(fcs_path)
```

Returns a named list of all keywords in the FCS file header.

#### get_fcs_metadata

```r
get_fcs_metadata(fcs_path)
```

Returns `c(acquisition_date, file_size)` for a single FCS file.

#### check_acquisition_dates

```r
check_acquisition_dates(parent_folder, output_dir)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `parent_folder` | character | *required* | Root directory to scan recursively for FCS files |
| `output_dir` | character | *required* | Output directory for the results CSV |

### Output

`acquisition_dates.csv` with columns:

| Column | Description |
|--------|-------------|
| `file_relative_path` | Path relative to `parent_folder` |
| `acquisition_date` | Date extracted from FCS header |
| `file_size` | File size in bytes |
| `unique_file_id` | Duplicate detection ID (identical date + size = same ID) |

### Duplicate Detection

Files with identical `acquisition_date` and `file_size` are assigned the same `unique_file_id`, flagging likely duplicate acquisitions.

---

## overlap_statistics.R

Computes Fleiss' kappa for inter-rater agreement across multiple gating strategies.

### Main Function

```r
compute_fleiss_kappa(label_lists, output_dir,
                     subset_specs = NULL, subset_labels = NULL,
                     row_range = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `label_lists` | list | *required* | Named list of label vectors (one per rater/gating strategy) |
| `output_dir` | character | *required* | Output directory for the log file |
| `subset_specs` | list | `NULL` | Specifications for subsetting cells before computing kappa |
| `subset_labels` | character vector | `NULL` | Labels to include in the analysis |
| `row_range` | integer vector | `NULL` | Row indices to subset (e.g., `1:1000`) |

### Algorithm

Uses `irr::kappam.fleiss()` to compute Fleiss' kappa, a measure of inter-rater reliability for categorical data with multiple raters.

### Output

- Returns a named list of kappa results
- Writes `fleiss_kappa_log.txt` to `output_dir`

### Interpretation

| Kappa | Agreement |
|-------|-----------|
| < 0.00 | Less than chance |
| 0.01 -- 0.20 | Slight |
| 0.21 -- 0.40 | Fair |
| 0.41 -- 0.60 | Moderate |
| 0.61 -- 0.80 | Substantial |
| 0.81 -- 1.00 | Almost perfect |

---

## compare_gatings.R

Compares gating strategies from multiple FlowJo workspaces on a common FCS file, generating Venn/Euler diagrams for each population.

### Main Function

```r
compare_gatings(pop_strings_list, target_fcs, label_dict,
                custom_name = NULL, output_dir,
                idx_labels = NULL, scale_size = TRUE)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pop_strings_list` | list | *required* | Population string lists from `load_flowjo_wsp` (one per workspace) |
| `target_fcs` | character | *required* | Name of the FCS file to compare across workspaces |
| `label_dict` | named list | *required* | Dictionary translating population names to common labels |
| `custom_name` | character | `NULL` | Custom prefix for output filenames |
| `output_dir` | character | *required* | Output directory |
| `idx_labels` | character vector | `NULL` | Labels for the workspaces (e.g., `c("Rater A", "Rater B")`) |
| `scale_size` | logical | `TRUE` | If `TRUE`, use area-proportional Euler diagrams; if `FALSE`, use Venn |

### Process

1. Extract labels for `target_fcs` from each workspace's population strings
2. Translate labels via `label_dict` to a common vocabulary
3. Generate a Venn or Euler diagram for each population showing overlap between workspaces

### Output

| File | Description |
|------|-------------|
| `fj_label_lists_{target_fcs}.csv` | Raw FlowJo labels per workspace |
| `label_lists_{target_fcs}.csv` | Translated labels per workspace |
| `{population}_{target_fcs}.png` | Venn/Euler diagram per population |

### src/R_shared Dependencies

`io_utils`, `cluster_utils` (`translate_labels`)

### Additional Sources

Sources `src/visualization/venn_diagrams.R` directly for diagram generation functions.
