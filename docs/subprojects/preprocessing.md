[< Back to Documentation Index](../README.md)

# src/preprocessing/ -- QC, Transformation, Compensation, Gating, Downsampling, Normalization

Thirteen scripts for fluorescence compensation, spectral unmixing, bead normalization, doublet detection, live/dead gating, quality control, signal transformation, downsampling, normalization, SuperCellCyto compression, automated gating, and gate utilities. The `preprocess_pipeline.R` orchestrator chains QC and transformation into a single configurable pipeline. These scripts operate on flowFrame objects and are used primarily by Pipeline A.

## Summary Table

| Script | Functions | Description |
|--------|-----------|-------------|
| `preprocess_pipeline.R` | `run_preprocessing()` | Orchestrates QC, transform, normalize |
| `quality_control.R` | `quality_control()` | PeacoQC / flowAI quality control |
| `transformation.R` | `transform_flowframe()`, `find_optimal_b_cofactor()` | Arcsinh / logicle / biexponential transforms |
| `downsampling.R` | `downsample_flowframes()` | Consistent / Maximal downsampling |
| `normalization.R` | `z_score_scaling()`, `quantile_norm()`, `harmony_normalize()` | Scaling and normalization |
| `compensation.R` | `extract_compensation_matrix()`, `apply_compensation()`, `load_compensation_matrix()` | Fluorescence compensation |
| `spectral_unmixing.R` | `build_unmixing_matrix()`, `unmix_spectral()` | Spectral cytometry NNLS unmixing |
| `bead_normalization.R` | `normalize_beads()`, `identify_beads()` | CyTOF bead-based normalization |
| `doublet_detection.R` | `detect_doublets()` | Doublet removal via FSC-A/FSC-H ratio or model |
| `live_dead_gating.R` | `gate_live_cells()`, `find_density_valley()`, `fit_gmm_threshold()` | Viability gating |
| `automated_gating.R` | `run_auto_gating()`, `define_gate()` | flowDensity / openCyto automated gating |
| `gate_utils.R` | `apply_gate_sequence()`, `boolean_gate()` | Gate application and boolean logic |
| `supercell.R` | `run_supercell()`, `run_supercell_per_cluster()` | SuperCellCyto metacell compression |

---

## preprocess_pipeline.R (Orchestrator)

Chains QC, transformation, normalization, and scaling into a single configurable pipeline. Called by Pipeline A Step 2.

### `run_preprocessing()`

```r
run_preprocessing(flow_frames, params, output_dir = NULL)
```

Iterates over all flowFrames and applies QC, transformation, quantile normalization, and z-score scaling in the order specified by `params$pp_order`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | named list | *required* | Named list of flowFrame objects |
| `params` | list | *required* | Configuration list (see keys below) |
| `output_dir` | character | `NULL` | Root output directory |

**Configuration Keys (from `params`):**

| Key | Type | Description |
|-----|------|-------------|
| `pp_order` | character | `"QT"` (QC then Transform) or `"TQ"` (Transform then QC) |
| `qc_type` | character | `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `markers_of_interest` | character vector | Markers to process |
| `offset_list` | numeric vector | Per-marker offsets for transformation |
| `width_list` | numeric vector | Per-marker widths for transformation |
| `transform_type` | character | `"arcsinh"`, `"logicle"`, `"biexponential"`, `"log"`, or `"None"` |
| `quant_norm` | logical | Apply quantile normalization |
| `scaling` | logical | Apply z-score scaling |
| `save_pp_fcs` | logical | Save preprocessed FCS files |
| `save_pp_csv` | logical | Save preprocessed CSVs |

**Returns:** Named list of preprocessed `flowFrame` objects.

**Input from:** `load_fcs_raw()` (loading/load_fcs_raw.R).
**Output to:** `downsample_flowframes()` (preprocessing/downsampling.R), `plot_umap()` (visualization/umap.R).

**Example:**

```r
source("preprocess_pipeline.R")

params <- list(
  pp_order            = "QT",
  qc_type             = "None",
  transform_type      = "arcsinh",
  markers_of_interest = c("CD3", "CD4", "CD8"),
  offset_list         = rep(0, 3),
  width_list          = rep(5, 3),
  quant_norm          = FALSE,
  scaling             = FALSE,
  save_pp_fcs         = FALSE,
  save_pp_csv         = FALSE
)

preprocessed <- run_preprocessing(flow_frames, params, output_dir = "/output/preprocessing")
```

---

## quality_control.R

Runs automated quality control on individual flowFrames to remove aberrant events.

### `quality_control()`

```r
quality_control(flowframe, sample_name, qc_type = "PeacoQC",
                marker_indices = NULL, output_dir = NULL)
```

Dispatches to PeacoQC, flowAI, or passthrough based on `qc_type`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `sample_name` | character | *required* | Sample identifier for reporting |
| `qc_type` | character | `"PeacoQC"` | Algorithm: `"PeacoQC"`, `"flowAI"`, `"None"` |
| `marker_indices` | integer vector | `NULL` | Column indices to evaluate; `NULL` uses all |
| `output_dir` | character | `NULL` | Directory for QC reports |

**Available Algorithms:**

| Algorithm | What It Removes | Package |
|-----------|----------------|---------|
| `"PeacoQC"` | Margin events, signal anomalies, monotone regions | PeacoQC |
| `"flowAI"` | Anomalous flow rate, signal acquisition, dynamic range | flowAI |
| `"None"` | Nothing -- passthrough | -- |

**Returns:** Filtered `flowFrame` with low-quality events removed.

**Input from:** `run_preprocessing()` (this subproject).
**Output to:** `transform_flowframe()` (transformation.R).

**Example:**

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

### `transform_flowframe()`

```r
transform_flowframe(flowframe, sample_name, markers_of_interest,
                    offset_list, width_list, transform_type = "arcsinh")
```

Applies the specified transform to each marker independently. Skips markers not found in the flowFrame with a warning.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `sample_name` | character | *required* | Sample identifier |
| `markers_of_interest` | character vector | *required* | Markers to transform |
| `offset_list` | numeric vector | *required* | Per-marker offset (cofactor `a`). Same length as `markers_of_interest`. |
| `width_list` | numeric vector | *required* | Per-marker width (cofactor `b`). Same length as `markers_of_interest`. |
| `transform_type` | character | `"arcsinh"` | See transformation types below |

**Transformation Types:**

| Type | flowCore Function | Use Case |
|------|-------------------|----------|
| `"arcsinh"` | `arcsinhTransform(a = -offset, b = 1/width)` | Mass cytometry (CyTOF). Standard cofactor = 5 |
| `"logicle"` | `logicleTransform()` | Fluorescence cytometry |
| `"biexponential"` | `biexponentialTransform()` | Fluorescence cytometry |
| `"log"` | `log10(x + offset)` | Simple log transformation |
| `"None"` | -- | No transformation |

**Returns:** Transformed `flowFrame`.

**Input from:** `run_preprocessing()`, `apply_compensation()`.
**Output to:** Downstream clustering and visualization.

### `find_optimal_b_cofactor()`

```r
find_optimal_b_cofactor(data, target_bimodality_coefficient = 0.9)
```

Uses `mousetrap::bimodality_coefficient` to search for the arcsinh cofactor that produces the target bimodality via `optimize()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | numeric vector | *required* | Expression values for a single marker |
| `target_bimodality_coefficient` | numeric | `0.9` | Target bimodality coefficient |

**Returns:** Numeric scalar: the optimal `b` cofactor value.

---

## downsampling.R

Reduces cell counts across samples for computational tractability while preserving relative proportions.

### `downsample_flowframes()`

```r
downsample_flowframes(flow_frames, sample_names, target_cells = 10000,
                      method = "Maximal", seed = 42,
                      save_fcs = FALSE, output_dir = NULL)
```

Downsamples cells using either equal-per-sample or proportional strategies.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Corresponding sample names |
| `target_cells` | integer | `10000` | Target total number of cells |
| `method` | character | `"Maximal"` | `"Consistent"` (equal per sample) or `"Maximal"` (proportional) |
| `seed` | integer | `42` | Random seed |
| `save_fcs` | logical | `FALSE` | Write downsampled FCS files |
| `output_dir` | character | `NULL` | Directory for saved files |

**Returns:** Named list of downsampled `flowFrame` objects.

**Input from:** `run_preprocessing()`, `load_fcs_raw()`.
**Output to:** `run_leiden_with_umap()` (clustering/leiden.R), `plot_umap()` (visualization/umap.R).

**Example:**

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

---

## normalization.R

Three normalization strategies for expression data.

### `z_score_scaling()`

```r
z_score_scaling(flowframe, marker_indices)
```

Scales each marker to zero mean and unit variance using `flowCore::transformList()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `marker_indices` | integer vector | *required* | Column indices to scale |

**Returns:** Scaled `flowFrame`.

### `quantile_norm()`

```r
quantile_norm(flowframe, marker_indices)
```

Applies quantile normalization across markers using `preprocessCore::normalize.quantiles()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `marker_indices` | integer vector | *required* | Column indices to normalize |

**Returns:** Quantile-normalized `flowFrame`.

### `harmony_normalize()`

```r
harmony_normalize(expr_matrix, batch_vector, ...)
```

Runs Harmony batch correction via `harmony::HarmonyMatrix()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | numeric matrix | *required* | Cells x markers matrix |
| `batch_vector` | character/factor | *required* | Batch labels per cell |
| `...` | | | Additional arguments to `harmony::HarmonyMatrix()` |

**Returns:** Batch-corrected expression matrix.

---

## compensation.R

Applies spillover compensation to conventional flow cytometry data.

### `extract_compensation_matrix()`

```r
extract_compensation_matrix(flowframe)
```

Searches flowFrame keywords (`SPILL`, `SPILLOVER`, `$SPILLOVER`, etc.) for a compensation matrix. Falls back to `flowCore::spillover()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |

**Returns:** Compensation matrix, or `NULL` if not found.

**Input from:** `load_fcs_raw()`.
**Output to:** `apply_compensation()`.

### `apply_compensation()`

```r
apply_compensation(flowframe, comp_matrix = NULL, auto_extract = TRUE)
```

Applies fluorescence compensation using `flowCore::compensate()`. If `comp_matrix` is `NULL` and `auto_extract = TRUE`, extracts the matrix from the flowFrame keywords.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `comp_matrix` | matrix | `NULL` | Compensation matrix; `NULL` triggers auto-extraction |
| `auto_extract` | logical | `TRUE` | Extract matrix from flowFrame keywords if `comp_matrix` is `NULL` |

**Returns:** Compensated `flowFrame`.

**Input from:** `load_fcs_raw()`.
**Output to:** `quality_control()`, `transform_flowframe()`.

### `load_compensation_matrix()`

```r
load_compensation_matrix(csv_path, sep = ";")
```

Loads a compensation matrix from an external CSV file.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_path` | character | *required* | Path to CSV file |
| `sep` | character | `";"` | CSV separator |

**Returns:** Compensation matrix with row/column names.

---

## spectral_unmixing.R

Unmixes spectral cytometry data using reference spectra from single-stain controls via non-negative least squares (NNLS).

### `build_unmixing_matrix()`

```r
build_unmixing_matrix(control_frames, unstained_frame, channels = NULL)
```

Builds an unmixing reference matrix from single-stain control flowFrames. Subtracts autofluorescence baseline (median of unstained control), normalizes each spectrum to unit sum.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `control_frames` | named list | *required* | Single-stain flowFrames (names = fluorochrome labels) |
| `unstained_frame` | flowFrame | *required* | Unstained control for autofluorescence baseline |
| `channels` | character vector | `NULL` | Spectral channel names; `NULL` auto-detects (excludes scatter/time) |

**Returns:** Unmixing matrix (channels x fluorochromes).

**Input from:** `load_fcs_raw()` (single-stain controls).
**Output to:** `unmix_spectral()`.

### `unmix_spectral()`

```r
unmix_spectral(flowframe, unmixing_matrix, autofluorescence = NULL)
```

Unmixes spectral cytometry data cell-by-cell using `nnls::nnls()`. Processes in chunks of 10,000 cells for memory efficiency. Preserves non-spectral channels (scatter, time).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | flowFrame with spectral channels |
| `unmixing_matrix` | matrix | *required* | Reference unmixing matrix from `build_unmixing_matrix()` |
| `autofluorescence` | numeric vector | `NULL` | Autofluorescence to subtract before unmixing |

**Returns:** Unmixed `flowFrame` with virtual marker channels.

**Input from:** `load_fcs_raw()`.
**Output to:** `quality_control()`, `transform_flowframe()`.

---

## bead_normalization.R

Normalizes CyTOF/mass cytometry data using bead standards to correct for signal decay over acquisition time.

### `normalize_beads()`

```r
normalize_beads(flowframe, bead_channels = NULL, bead_type = "four_element",
                remove_beads = TRUE, baseline = "median",
                window_size = 500, output_dir = NULL)
```

Identifies bead events, computes per-channel sliding-window normalization factors, and optionally removes bead events. Generates diagnostic plots when `output_dir` is provided.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | flowFrame with bead and cell events |
| `bead_channels` | character vector | `NULL` | Bead channel names; `NULL` uses defaults for `bead_type` |
| `bead_type` | character | `"four_element"` | `"four_element"` (standard EQ beads) or `"custom"` |
| `remove_beads` | logical | `TRUE` | Remove bead events after normalization |
| `baseline` | character | `"median"` | `"median"` (normalize to global median) or `"start"` (first window) |
| `window_size` | integer | `500` | Events per smoothing window |
| `output_dir` | character | `NULL` | Directory for diagnostic plots |

**Returns:**

```r
list(
  normalized_ff        = <flowFrame>,
  bead_indices         = logical vector,
  normalization_factors = matrix
)
```

**Input from:** `load_fcs_raw()`.
**Output to:** `quality_control()`, `transform_flowframe()`.

### `identify_beads()`

```r
identify_beads(bead_expr, method = "mahalanobis")
```

Identifies bead events using Mahalanobis distance or intensity thresholding.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bead_expr` | matrix | *required* | Matrix of bead channel expression values |
| `method` | character | `"mahalanobis"` | `"mahalanobis"` or `"threshold"` |

**Returns:** Logical vector (`TRUE` = bead event).

---

## doublet_detection.R

Detects and removes doublets using forward scatter area vs. height ratios or bivariate density gating.

### `detect_doublets()`

```r
detect_doublets(flowframe, area_col = "FSC-A", height_col = "FSC-H",
                method = "ratio", threshold = 1.1, width_col = NULL)
```

Classifies cells as singlets or doublets. The "ratio" method flags events where FSC-A/FSC-H exceeds a threshold. The "model" method fits a linear model and flags residual outliers. Partial matching is used for column names.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `area_col` | character | `"FSC-A"` | Forward scatter area column |
| `height_col` | character | `"FSC-H"` | Forward scatter height column |
| `method` | character | `"ratio"` | `"ratio"` (FSC-A/FSC-H threshold) or `"model"` (bivariate density) |
| `threshold` | numeric | `1.1` | Ratio threshold for "ratio" method |
| `width_col` | character | `NULL` | Optional FSC-W column for additional filtering |

**Returns:**

```r
list(
  filtered_ff     = <flowFrame>,   # Singlets only
  doublet_indices = integer vector, # Row indices of doublets
  singlet_rate    = numeric         # Fraction of singlets
)
```

**Input from:** `load_fcs_raw()`, `apply_compensation()`.
**Output to:** `quality_control()`, `transform_flowframe()`.

---

## live_dead_gating.R

Gates on live cells using a viability marker.

### `gate_live_cells()`

```r
gate_live_cells(flowframe, viability_channel, method = "threshold",
                threshold = NULL, dead_high = TRUE)
```

Classifies cells as live or dead based on a viability marker. Supports fixed threshold and Gaussian mixture model methods. When `threshold = NULL` with the "threshold" method, auto-detects via density valley.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `viability_channel` | character | *required* | Name of live/dead marker column |
| `method` | character | `"threshold"` | `"threshold"` (fixed/auto cutoff) or `"mixture"` (GMM) |
| `threshold` | numeric | `NULL` | Manual threshold; `NULL` auto-detects |
| `dead_high` | logical | `TRUE` | If `TRUE`, dead cells are viability-HIGH (standard for amine-reactive dyes) |

**Returns:**

```r
list(
  live_ff        = <flowFrame>,   # Live cells only
  dead_indices   = integer vector, # Row indices of dead cells
  viability_rate = numeric         # Fraction of live cells
)
```

**Input from:** `load_fcs_raw()`, `detect_doublets()`.
**Output to:** `quality_control()`, `transform_flowframe()`.

### `find_density_valley()` (internal)

Finds the valley between two density peaks for auto-threshold detection. Falls back to median if no bimodal pattern is detected.

### `fit_gmm_threshold()` (internal)

Fits a 2-component Gaussian mixture model using EM and returns the threshold at the intersection of the two Gaussians.

---

## automated_gating.R

Data-driven automated gating using flowDensity or openCyto, plus programmatic gate construction.

### `run_auto_gating()`

```r
run_auto_gating(flowframe, method = "flowDensity", gating_template = NULL,
                channels, population = NULL)
```

Runs automated gating using density-based cutpoints (flowDensity) or template-based gating (openCyto).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `method` | character | `"flowDensity"` | `"flowDensity"` or `"openCyto"` |
| `gating_template` | character | `NULL` | Path to openCyto gating template CSV (required for "openCyto") |
| `channels` | character vector | *required* | Channels to gate on |
| `population` | character | `NULL` | Target population name |

**Returns:**

```r
list(
  gated_ff         = <flowFrame>,
  gate_coords      = list(...),
  population_stats = list(n_parent, n_gated, pct_parent, population, method)
)
```

**Input from:** `load_fcs_raw()`, `apply_compensation()`.
**Output to:** `apply_gate_sequence()` (gate_utils.R), downstream analysis.

### `define_gate()`

```r
define_gate(gate_type, channels, coords)
```

Constructs flowCore gate objects from specified coordinates. Supports rectangle, polygon, threshold, and quadrant gate types.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gate_type` | character | *required* | `"rectangle"`, `"polygon"`, `"threshold"`, `"quadrant"` |
| `channels` | character vector | *required* | Channel names (1 for threshold, 2 for rectangle/polygon/quadrant) |
| `coords` | varies | *required* | Gate coordinates (format depends on `gate_type`) |

**Coords format by gate_type:**

| Gate Type | Coords Format |
|-----------|--------------|
| `"rectangle"` | `list(x = c(min, max), y = c(min, max))` |
| `"polygon"` | Matrix with 2 columns (x, y vertices) |
| `"threshold"` | Single numeric value (cutpoint) |
| `"quadrant"` | `list(x = numeric, y = numeric)` (intersection point) |

**Returns:** flowCore gate object (`rectangleGate` or `polygonGate`). For `"quadrant"`, returns a named list of 4 `rectangleGate` objects (Q1-Q4).

**Input from:** User-defined coordinates.
**Output to:** `apply_gate_sequence()` (gate_utils.R).

**Example:**

```r
source("automated_gating.R")

# Rectangle gate
gate <- define_gate("rectangle", c("CD3", "CD4"),
                    coords = list(x = c(2, Inf), y = c(1.5, Inf)))

# Threshold gate
gate <- define_gate("threshold", "CD3", coords = 2.5)

# Quadrant gate
gates <- define_gate("quadrant", c("CD3", "CD4"),
                     coords = list(x = 2.0, y = 1.5))
```

---

## gate_utils.R

Utilities for applying sequences of flowCore gates and combining gate results with boolean logic.

### `apply_gate_sequence()`

```r
apply_gate_sequence(flowframe, gate_list, logical_op = "and")
```

Applies a list of flowCore gate objects to a flowFrame and combines results using AND or OR logic.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `gate_list` | list | *required* | List of flowCore gate objects |
| `logical_op` | character | `"and"` | `"and"` (intersection) or `"or"` (union) |

**Returns:**

```r
list(
  gated_ff         = <flowFrame>,
  gate_indices     = integer vector,
  population_stats = list(n_parent, n_gated, pct_parent, n_gates, logical_op)
)
```

**Input from:** `define_gate()` (automated_gating.R), `run_auto_gating()`.
**Output to:** `boolean_gate()` (this file), downstream analysis.

### `boolean_gate()`

```r
boolean_gate(flowframe, gate_results, expression)
```

Evaluates a boolean expression string referencing named gate results. Gate names in the expression are replaced with their logical vectors and the expression is evaluated.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame |
| `gate_results` | named list | *required* | Named list of logical vectors from previous gates |
| `expression` | character | *required* | Boolean expression, e.g., `"CD3+ & CD4+ & !CD8+"` |

**Returns:** Logical vector of included cells.

**Example:**

```r
source("gate_utils.R")

gate_results <- list(
  "CD3+" = cd3_logical,
  "CD4+" = cd4_logical,
  "CD8+" = cd8_logical
)
selected <- boolean_gate(ff, gate_results, "CD3+ & CD4+ & !CD8+")
```

---

## supercell.R

Compresses single-cell cytometry data into supercells (metacells) using SuperCellCyto.

### `run_supercell()`

```r
run_supercell(expr_matrix, sample_labels, markers,
              gam = 20, k_knn = 30, n_pc = 10,
              do_pca = TRUE, seed = 42)
```

Compresses cells into supercells using `SuperCellCyto::runSuperCellCyto()`. Higher `gam` values produce fewer, larger supercells.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | numeric matrix | *required* | Cells x markers matrix |
| `sample_labels` | character vector | *required* | Sample label per cell |
| `markers` | character vector | *required* | Marker columns to use |
| `gam` | numeric | `20` | Graining level (higher = fewer supercells, ~1/gam compression) |
| `k_knn` | integer | `30` | Nearest neighbors for graph construction |
| `n_pc` | integer | `10` | PCs for initial dimensionality reduction |
| `do_pca` | logical | `TRUE` | Run PCA before supercell construction |
| `seed` | integer | `42` | Random seed |

**Returns:**

```r
list(
  supercell_matrix = matrix,         # Supercell expression (supercells x markers)
  supercell_labels = character,      # Sample labels per supercell
  cell_map         = data.table,     # CellId -> SuperCellId mapping
  supercell_sizes  = data.table      # SuperCellId -> size
)
```

**Input from:** `build_combined_expr_matrix()` (R_shared/flowframe_utils.R), `combine_expression_matrices()` (analysis/expression_matrices.R).
**Output to:** `run_leiden()` / `run_flowsom()` / `run_phenograph()` (clustering/), `plot_umap_from_matrix()` (visualization/).

### `run_supercell_per_cluster()`

```r
run_supercell_per_cluster(expr_matrix, cluster_ids, sample_labels,
                           markers, gam = 20, seed = 42)
```

Applies supercell compression independently within each cluster, preserving cluster boundaries. Clusters too small for meaningful compression (fewer than `gam * 2` cells) are kept as-is.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | numeric matrix | *required* | Cells x markers matrix |
| `cluster_ids` | integer vector | *required* | Cluster assignments per cell |
| `sample_labels` | character vector | *required* | Sample labels per cell |
| `markers` | character vector | *required* | Marker columns to use |
| `gam` | numeric | `20` | Graining level per cluster |
| `seed` | integer | `42` | Random seed |

**Returns:**

```r
list(
  supercell_matrix      = matrix,     # Combined supercell expression
  supercell_cluster_ids = integer,    # Cluster assignment per supercell
  supercell_labels      = character,  # Sample labels per supercell
  cell_map              = data.table  # With cluster, CellId, SuperCellId, OriginalIndex
)
```

**Input from:** Clustering output + `build_combined_expr_matrix()`.
**Output to:** `run_differential_abundance()` / `run_differential_expression()` (analysis/), `plot_umap_from_matrix()` (visualization/).

**Example:**

```r
source("supercell.R")

result <- run_supercell(
  expr_matrix   = combined_matrix,
  sample_labels = sample_labels,
  markers       = c("CD3", "CD4", "CD8", "CD19"),
  gam           = 20
)
cat("Compression:", nrow(combined_matrix), "->", nrow(result$supercell_matrix), "\n")
```
