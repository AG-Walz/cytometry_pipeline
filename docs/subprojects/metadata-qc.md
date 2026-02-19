[< Back to Documentation Index](../README.md)

# src/metadata_qc/ — Metadata & Quality Assurance

Four scripts for validating study metadata, checking FCS acquisition dates, computing inter-rater agreement, and comparing gating strategies across FlowJo workspaces. These scripts support Pipeline D and can also be used standalone for quality checks.

---

## check_metadata.R

Validates and visualizes study metadata, detecting missing values and temporal patterns. Useful for spotting data quality issues before running analysis.

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
| `metadata_file` | character | *required* | Path to metadata CSV |
| `output_dir` | character | *required* | Output directory |
| `sep` | character | `","` | Column separator |
| `skip_columns` | character vector | `NULL` | Columns to exclude |
| `only_cols` | character vector | `NULL` | If set, analyze only these columns |
| `n_splits` | integer | `3` | Number of temporal splits |
| `placeholders` | character vector | see above | Strings treated as missing |

#### plot_temporal_heatmap

```r
plot_temporal_heatmap(metadata_file, output_dir, sep = ",",
                      skip_columns = NULL, date_column = "Datum")
```

### What It Does

1. Reads the metadata CSV
2. Converts placeholder strings (`"?"`, `"na"`, `"n/a"`, `"null"`, `"-"`) to `"Unassigned"`
3. Auto-detects date columns
4. Generates faceted barplots showing value distributions for each column
5. Splits data into temporal groups (by date) to detect trends
6. Creates a heatmap of feature distributions over time

### Example Input

```
Sample_ID;Datum;Age;Sex;Diagnosis;Treatment_Group;Response;Batch
ICK-001;15/03/2024;45;M;AML;Chemo;pCR;Batch1
ICK-002;15/03/2024;62;F;AML;Chemo;npCR;Batch1
ICK-003;20/03/2024;?;M;ALL;Immuno;?;Batch2
ICK-004;20/03/2024;38;F;AML;Chemo;pCR;Batch2
```

### Generated Files

| File | Description |
|------|-------------|
| `consensus_barplots_compact.png` | Faceted barplots for all metadata columns |
| `consensus_barplots_compact_bySplit_1.png` | Barplots for temporal split 1 |
| `consensus_barplots_compact_bySplit_2.png` | Barplots for temporal split 2 |
| `metadata_summary_table.csv` | Summary statistics per column |
| `temporal_feature_heatmap.png` | Feature distributions over time |

### Example

```r
source("check_metadata.R")
check_metadata(
  metadata_file = "/data/study_metadata.csv",
  output_dir    = "/output/metadata_qc",
  sep           = ";",
  sort_column   = "Datum",
  n_splits      = 4
)
```

---

## check_acquisition.R

Extracts and validates acquisition metadata from FCS file headers. Detects duplicate files by comparing acquisition dates and file sizes.

### Functions

```r
get_fcs_keyword(fcs_path, keyword_name)    # Get a single FCS keyword
explore_fcs_keywords(fcs_path)             # Get all FCS keywords
get_fcs_metadata(fcs_path)                 # Get acquisition date + file size

check_acquisition_dates(parent_folder, output_dir)
```

### check_acquisition_dates

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `parent_folder` | character | *required* | Root directory to scan recursively |
| `output_dir` | character | *required* | Output directory |

Recursively scans `parent_folder` for FCS files, extracts acquisition dates from the `$DATE` keyword, and detects duplicates.

### Output

**`acquisition_dates.csv`:**

```
file_relative_path;acquisition_date;file_size;unique_file_id
batch1/sample_001.fcs;2024-03-15;52428800;1
batch1/sample_002.fcs;2024-03-15;48234496;2
batch2/sample_001_copy.fcs;2024-03-15;52428800;1
```

Files with identical `acquisition_date` AND `file_size` get the same `unique_file_id`, flagging likely duplicates (e.g., `sample_001.fcs` and `sample_001_copy.fcs` in the example above both have ID 1).

### Example

```r
source("check_acquisition.R")

# Quick keyword lookup
date <- get_fcs_keyword("/data/sample_001.fcs", "$DATE")
cat("Acquired:", date, "\n")
# Acquired: 15-MAR-2024

# Full acquisition check
check_acquisition_dates(
  parent_folder = "/data/all_fcs/",
  output_dir    = "/output/acquisition"
)
```

---

## overlap_statistics.R

Computes Fleiss' kappa for inter-rater agreement across multiple gating strategies. Used by Pipeline D to quantify agreement between analysts.

### Function

```r
compute_fleiss_kappa(label_lists, output_dir,
                     subset_specs = NULL, subset_labels = NULL,
                     row_range = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `label_lists` | list | *required* | Named list of label vectors (one per rater) |
| `output_dir` | character | *required* | Output directory |
| `subset_specs` | list | `NULL` | Subsetting specifications |
| `subset_labels` | character vector | `NULL` | Labels to include |
| `row_range` | integer vector | `NULL` | Row indices to subset |

### Algorithm

Uses `irr::kappam.fleiss()` to compute Fleiss' kappa — a measure of inter-rater reliability for categorical data with two or more raters.

### Output

Returns a named list of kappa results. Also writes `fleiss_kappa_log.txt`:

```
=== Fleiss' Kappa Results ===

Population: CD4+ T cells
  Kappa = 0.85 (Almost perfect agreement)
  n_events = 50000, n_raters = 3

Population: CD8+ T cells
  Kappa = 0.72 (Substantial agreement)
  n_events = 50000, n_raters = 3

Population: B cells
  Kappa = 0.91 (Almost perfect agreement)
  n_events = 50000, n_raters = 3
```

### Kappa Interpretation

| Kappa | Agreement Level |
|-------|----------------|
| < 0.00 | Less than chance |
| 0.01 – 0.20 | Slight |
| 0.21 – 0.40 | Fair |
| 0.41 – 0.60 | Moderate |
| 0.61 – 0.80 | Substantial |
| 0.81 – 1.00 | Almost perfect |

### Example

```r
source("overlap_statistics.R")

# label_lists: one vector per rater, same length, aligned by event index
label_lists <- list(
  "Analyst A" = c("CD4+", "CD8+", "CD4+", "NK", ...),
  "Analyst B" = c("CD4+", "CD8+", "CD4+", "B cell", ...),
  "Analyst C" = c("CD4+", "CD4+", "CD4+", "NK", ...)
)

kappa_results <- compute_fleiss_kappa(
  label_lists = label_lists,
  output_dir  = "/output/overlap"
)
```

### Packages

`irr`

---

## compare_gatings.R

Compares gating strategies from multiple FlowJo workspaces on a common FCS file. Generates per-population Venn/Euler diagrams showing overlap between workspaces.

### Function

```r
compare_gatings(pop_strings_list, target_fcs, label_dict,
                custom_name = NULL, output_dir,
                idx_labels = NULL, scale_size = TRUE)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pop_strings_list` | list | *required* | Population strings from `load_flowjo_wsp` |
| `target_fcs` | character | *required* | FCS filename to compare across workspaces |
| `label_dict` | named list | *required* | Gate path → display label mapping |
| `custom_name` | character vector | `NULL` | Labels for workspaces |
| `output_dir` | character | *required* | Output directory |
| `idx_labels` | character vector | `NULL` | Additional label vector (e.g. from Pipeline C) |
| `scale_size` | logical | `TRUE` | `TRUE` = Euler, `FALSE` = Venn |

### Process

1. Extract labels for `target_fcs` from each workspace
2. Translate labels via `label_dict` to a common vocabulary
3. Build a binary comparison matrix (event × workspace × population)
4. Generate Venn/Euler diagram for each population

### Generated Files

| File | Description |
|------|-------------|
| `fj_label_lists_{target_fcs}.csv` | Raw FlowJo labels per workspace |
| `label_lists_{target_fcs}.csv` | Translated labels per workspace |
| `{population}_{target_fcs}.png` | Venn/Euler diagram per population |

### Example

```r
source("compare_gatings.R")
source("../visualization/venn_diagrams.R")

comparison <- compare_gatings(
  pop_strings_list = all_labels,
  target_fcs       = "sample_001.fcs",
  label_dict       = my_label_dict,
  custom_name      = c("Analyst A", "Analyst B"),
  output_dir       = "/output/comparison",
  scale_size       = TRUE
)
```

### Dependencies

`io_utils`, `cluster_utils` (`translate_labels`). Sources `src/visualization/venn_diagrams.R` directly.
