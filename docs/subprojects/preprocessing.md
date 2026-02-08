[< Back to Documentation Index](../README.md)

# src/preprocessing/ -- QC, Transformation, Downsampling, Normalization

5 scripts for quality control, signal transformation, cell downsampling, and data normalization. Includes a pipeline orchestrator that chains these steps together.

---

## quality_control.R

Runs automated quality control on individual flowFrames.

### Main Function

```r
quality_control(flowframe, sample_name, qc_type = "PeacoQC",
                marker_indices = NULL, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame to filter |
| `sample_name` | character | *required* | Sample identifier for reporting |
| `qc_type` | character | `"PeacoQC"` | Algorithm: `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `marker_indices` | integer vector | `NULL` | Column indices to evaluate; `NULL` uses all |
| `output_dir` | character | `NULL` | Directory for QC reports and plots |

### Algorithms

| Algorithm | Function | Package |
|-----------|----------|---------|
| `"PeacoQC"` | `PeacoQC::PeacoQC()` | PeacoQC |
| `"flowAI"` | `flowAI::flow_auto_qc()` | flowAI |
| `"None"` | Passthrough (no filtering) | -- |

### Input / Output

- **Input**: `flowFrame`
- **Output**: Filtered `flowFrame` with low-quality events removed

---

## transformation.R

Applies mathematical transformations to fluorescence/mass cytometry data.

### Main Function

```r
transform_flowframe(flowframe, sample_name, markers_of_interest,
                    offset_list, width_list, transform_type = "arcsinh")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `sample_name` | character | *required* | Sample identifier |
| `markers_of_interest` | character vector | *required* | Markers to transform |
| `offset_list` | named numeric | *required* | Per-marker offset (cofactor `a` for arcsinh) |
| `width_list` | named numeric | *required* | Per-marker width (cofactor `b` for arcsinh) |
| `transform_type` | character | `"arcsinh"` | Transformation type (see table below) |

### Transformation Types

| Type | flowCore Function | Formula / Notes |
|------|-------------------|-----------------|
| `"arcsinh"` | `arcsinhTransform(a = -offset, b = 1/width)` | Standard for mass cytometry; typical cofactor = 5 |
| `"logicle"` | `logicleTransform()` | Biexponential approximation for fluorescence data |
| `"biexponential"` | `biexponentialTransform()` | Full biexponential |
| `"None"` | -- | No transformation applied |

### Helper Function

```r
# Find optimal arcsinh cofactor by maximizing bimodality
find_optimal_b_cofactor(data, target_bimodality_coefficient = 0.9)
```

Uses `mousetrap::bimodality_coefficient` to search for the cofactor value that produces the target bimodality.

### Input / Output

- **Input**: `flowFrame`
- **Output**: Transformed `flowFrame`

---

## downsampling.R

Reduces cell counts across samples for computational tractability while preserving relative proportions.

### Main Function

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
| `method` | character | `"Maximal"` | Downsampling strategy (see below) |
| `seed` | integer | `42` | Random seed for reproducibility |
| `save_fcs` | logical | `FALSE` | Whether to write downsampled FCS files |
| `output_dir` | character | `NULL` | Directory for saved FCS files |

### Methods

| Method | Behavior |
|--------|----------|
| `"Consistent"` | Equal number of cells per sample (`target_cells` each) |
| `"Maximal"` | Proportional to original counts, scaled so total equals `target_cells * n_samples` |

### Output

List of downsampled `flowFrame` objects. When `save_fcs = TRUE`, writes `ds_*.fcs` files to `output_dir`.

---

## normalization.R

Provides three normalization strategies for expression data.

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

Applies quantile normalization across markers. Returns a normalized `flowFrame`. Requires the `preprocessCore` package.

#### harmony_normalize

```r
harmony_normalize(expr_matrix, batch_vector, ...)
```

Runs Harmony batch correction on an expression matrix. Returns a corrected matrix. Requires the `harmony` package.

---

## preprocess_pipeline.R (Orchestrator)

Chains QC, transformation, normalization, and scaling into a single configurable pipeline.

### Main Function

```r
run_preprocessing(flow_frames, params, output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `params` | list | *required* | Configuration list (see keys below) |
| `output_dir` | character | `NULL` | Root output directory |

### Configuration Keys (`params`)

| Key | Type | Description |
|-----|------|-------------|
| `pp_order` | character | Processing order: `"QT"` (QC then Transform) or `"TQ"` (Transform then QC) |
| `qc_type` | character | QC algorithm: `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `markers_of_interest` | character vector | Markers to process |
| `offset_list` | named numeric | Per-marker offsets for transformation |
| `width_list` | named numeric | Per-marker widths for transformation |
| `transform_type` | character | `"arcsinh"`, `"logicle"`, `"biexponential"`, or `"None"` |
| `quant_norm` | logical | Whether to apply quantile normalization |
| `scaling` | logical | Whether to apply z-score scaling |
| `save_pp_fcs` | logical | Save preprocessed flowFrames as FCS files |
| `save_pp_csv` | logical | Save preprocessed expression matrices as CSVs |

### Processing Steps

1. **QC + Transform** (in the order specified by `pp_order`)
2. **Quantile normalization** (optional, when `quant_norm = TRUE`)
3. **Z-score scaling** (optional, when `scaling = TRUE`)
4. **Save outputs** (FCS and/or CSV based on configuration)

### Output

- Named list of preprocessed `flowFrame` objects
- Files written to:
  - `output_dir/preprocessed_files/` -- preprocessed FCS/CSV files
  - `output_dir/qc_results/` -- QC reports and plots
