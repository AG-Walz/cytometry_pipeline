[< Back to Documentation Index](../README.md)

# src/preprocessing/ — QC, Transformation, Downsampling, Normalization

Five scripts for quality control, signal transformation, cell downsampling, and data normalization. The `preprocess_pipeline.R` orchestrator chains QC and transformation into a single configurable pipeline. These scripts operate on flowFrame objects and are used by Pipeline A.

---

## preprocess_pipeline.R (Orchestrator)

Chains QC, transformation, normalization, and scaling into a single configurable pipeline. This is the function called by Pipeline A's Step 2.

### Function

```r
run_preprocessing(flow_frames, params, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `params` | list | *required* | Configuration list (see keys below) |
| `output_dir` | character | `NULL` | Root output directory |

### Configuration Keys (from `params`)

| Key | Type | Description |
|-----|------|-------------|
| `pp_order` | character | `"QT"` (QC then Transform) or `"TQ"` (Transform then QC) |
| `qc_type` | character | `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `markers_of_interest` | character vector | Markers to process |
| `offset_list` | numeric vector | Per-marker offsets for transformation |
| `width_list` | numeric vector | Per-marker widths for transformation |
| `transform_type` | character | `"arcsinh"`, `"logicle"`, `"biexponential"`, or `"None"` |
| `quant_norm` | logical | Apply quantile normalization |
| `scaling` | logical | Apply z-score scaling |
| `save_pp_fcs` | logical | Save preprocessed FCS files |
| `save_pp_csv` | logical | Save preprocessed CSVs |

### Processing Steps

1. **QC + Transform** (in the order specified by `pp_order`)
2. **Quantile normalization** (optional, when `quant_norm = TRUE`)
3. **Z-score scaling** (optional, when `scaling = TRUE`)
4. **Save outputs** (FCS and/or CSV based on flags)

### Returns

Named list of preprocessed `flowFrame` objects.

### Example

```r
source("preprocess_pipeline.R")

params <- list(
  pp_order            = "QT",
  qc_type             = "None",
  transform_type      = "arcsinh",
  markers_of_interest = c("CD3", "CD4", "CD8"),
  offset_list         = rep(0, 3),
  width_list          = rep(5, 3),
  save_pp_fcs         = FALSE,
  save_pp_csv         = FALSE
)

preprocessed <- run_preprocessing(flow_frames, params, output_dir = "/output/preprocessing")
```

**Console output:**

```
  Processing order: QT (QC → Transform)
  Processing sample_001 ...
    QC: None (skipped)
    Transform: arcsinh on 3 markers
  Processing sample_002 ...
    QC: None (skipped)
    Transform: arcsinh on 3 markers
  Preprocessing complete.
```

### Output Files (when enabled)

```
output_dir/
├── preprocessed_files/
│   ├── sample_001_preprocessed.csv
│   ├── sample_001_preprocessed.fcs
│   └── ...
└── qc_results/
    ├── sample_001_qc_report.csv
    └── ...
```

---

## quality_control.R

Runs automated quality control on individual flowFrames to remove aberrant events.

### Function

```r
quality_control(flowframe, sample_name, qc_type = "PeacoQC",
                marker_indices = NULL, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `sample_name` | character | *required* | Sample identifier for reporting |
| `qc_type` | character | `"PeacoQC"` | Algorithm: `"PeacoQC"`, `"flowAI"`, `"None"` |
| `marker_indices` | integer vector | `NULL` | Column indices to evaluate; `NULL` uses all |
| `output_dir` | character | `NULL` | Directory for QC reports |

### Available Algorithms

| Algorithm | What It Removes | Package |
|-----------|----------------|---------|
| `"PeacoQC"` | Margin events, signal anomalies, monotone regions | PeacoQC |
| `"flowAI"` | Anomalous flow rate, signal acquisition, dynamic range | flowAI |
| `"None"` | Nothing — passthrough | — |

### Returns

Filtered `flowFrame` with low-quality events removed.

### Example

```r
source("quality_control.R")
filtered_ff <- quality_control(
  flowframe   = my_flowframe,
  sample_name = "sample_001",
  qc_type     = "PeacoQC",
  output_dir  = "/output/qc"
)
cat("Before QC:", nrow(flowCore::exprs(my_flowframe)), "cells\n")
cat("After QC:", nrow(flowCore::exprs(filtered_ff)), "cells\n")
```

---

## transformation.R

Applies mathematical transformations to fluorescence/mass cytometry channels.

### Function

```r
transform_flowframe(flowframe, sample_name, markers_of_interest,
                    offset_list, width_list, transform_type = "arcsinh")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `sample_name` | character | *required* | Sample identifier |
| `markers_of_interest` | character vector | *required* | Markers to transform |
| `offset_list` | numeric vector | *required* | Per-marker offset (cofactor `a`) |
| `width_list` | numeric vector | *required* | Per-marker width (cofactor `b`) |
| `transform_type` | character | `"arcsinh"` | See transformation types below |

**Important:** `offset_list` and `width_list` must have the same length as `markers_of_interest`.

### Transformation Types

| Type | flowCore Function | Use Case |
|------|-------------------|----------|
| `"arcsinh"` | `arcsinhTransform(a = -offset, b = 1/width)` | Mass cytometry (CyTOF). Standard cofactor = 5 |
| `"logicle"` | `logicleTransform()` | Fluorescence cytometry |
| `"biexponential"` | `biexponentialTransform()` | Fluorescence cytometry |
| `"None"` | — | No transformation |

### Returns

Transformed `flowFrame`.

### Example

```r
source("transformation.R")
transformed <- transform_flowframe(
  flowframe           = my_ff,
  sample_name         = "sample_001",
  markers_of_interest = c("CD3", "CD4", "CD8"),
  offset_list         = c(0, 0, 0),
  width_list          = c(5, 5, 5),
  transform_type      = "arcsinh"
)
```

### Helper Function

```r
find_optimal_b_cofactor(data, target_bimodality_coefficient = 0.9)
```

Uses `mousetrap::bimodality_coefficient` to search for the arcsinh cofactor that produces the target bimodality.

---

## downsampling.R

Reduces cell counts across samples for computational tractability while preserving relative proportions.

### Function

```r
downsample_flowframes(flow_frames, sample_names, target_cells = 10000,
                      method = "Maximal", seed = 42,
                      save_fcs = FALSE, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Corresponding sample names |
| `target_cells` | integer | `10000` | Target number of cells |
| `method` | character | `"Maximal"` | Downsampling strategy |
| `seed` | integer | `42` | Random seed |
| `save_fcs` | logical | `FALSE` | Write downsampled FCS files |
| `output_dir` | character | `NULL` | Directory for saved files |

### Methods

| Method | Behavior | Best For |
|--------|----------|----------|
| `"Consistent"` | Equal number of cells per sample | Balanced representation |
| `"Maximal"` | Proportional to original counts | Preserving proportions |

### Returns

List of downsampled `flowFrame` objects.

### Example

```r
source("downsampling.R")
downsampled <- downsample_flowframes(
  flow_frames  = my_frames,
  sample_names = c("S1", "S2", "S3"),
  target_cells = 5000,
  method       = "Consistent",
  seed         = 42
)
```

**Console output:**

```
Downsampling 3 samples (total: 150000 cells)
Method: Consistent | Target: 5000 cells
Consistent: taking 1666 cells per sample
  S1 : 50000 -> 1666
  S2 : 50000 -> 1666
  S3 : 50000 -> 1666
Total after downsampling: 4998 cells
```

---

## normalization.R

Three normalization strategies for expression data.

### Functions

#### z_score_scaling

```r
z_score_scaling(flowframe, marker_indices)
```

Scales each marker to zero mean and unit variance. Returns a scaled `flowFrame`.

#### quantile_norm

```r
quantile_norm(flowframe, marker_indices)
```

Applies quantile normalization across markers. Requires `preprocessCore` package.

#### harmony_normalize

```r
harmony_normalize(expr_matrix, batch_vector, ...)
```

Runs Harmony batch correction. Requires `harmony` package.

### Example

```r
source("normalization.R")
scaled_ff <- z_score_scaling(my_ff, marker_indices = 1:6)
```
