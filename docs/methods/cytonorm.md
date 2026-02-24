# CytoNorm: Quantile Normalization for Cytometry Batches

Cytometry-specific batch correction that uses FlowSOM-stratified quantile normalization with reference samples to align marker distributions across acquisition batches while preserving cell-type-specific expression patterns.

---

## When to Use

- You have multi-batch cytometry data with systematic intensity shifts between batches.
- You have reference samples (bead standards, aliquots of the same biological sample, or technical replicates) measured in every batch.
- You want a normalization method purpose-built for cytometry rather than adapted from genomics.
- You need batch correction before transformation or at the raw-intensity level.

CytoNorm is the gold standard when shared reference samples are available. If you lack reference samples across batches, use [Harmony](batch-correction.md) instead, which requires only batch labels.

---

## How It Works

CytoNorm normalizes marker distributions in a cell-type-aware manner by combining FlowSOM clustering with quantile mapping.

In the training phase, CytoNorm clusters the reference samples using FlowSOM to identify cell populations that should have identical marker expression across batches. Within each FlowSOM cluster, quantile functions are learned for each marker in each batch. These quantile functions describe how the marker distribution in a given batch deviates from a target distribution (either the mean across batches or a designated reference batch). The per-cluster, per-marker, per-batch mapping functions are stored as the normalization model.

In the normalization phase, the learned model is applied to all samples (both reference and experimental). Each cell is first assigned to its nearest FlowSOM cluster, then the appropriate quantile mapping is applied to each marker. Because the correction is stratified by cell population, different cell types can receive different correction magnitudes. This is critical because batch effects in cytometry are often non-uniform: a fluorescence shift that is large for dim populations may be negligible for bright populations due to the nonlinear relationship between fluorescence intensity and detector response.

The cluster-stratified approach is what distinguishes CytoNorm from simple global quantile normalization. [Van Gassen et al. (2020)](https://doi.org/10.1002/cyto.a.23904) demonstrated that global quantile normalization can distort the relative expression patterns between cell types, while CytoNorm's stratified approach preserves biological heterogeneity. In benchmarks using spike-in reference standards, CytoNorm achieved the best alignment of reference sample distributions while maintaining separation between biologically distinct populations.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `train_frames` | (required) | List of reference flowFrames used to learn the normalization model |
| `validation_frames` | (required) | List of experimental flowFrames to normalize |
| `batch_labels` | (required) | Character vector of batch assignments (one per frame, covering both sets) |
| `markers` | (required) | Character vector of markers to normalize |
| `n_clusters` | 10 | Number of FlowSOM clusters for stratified normalization |
| `goal` | "mean" | Normalization target: "mean" (average across batches) or a specific batch name |
| `seed` | 42 | Random seed for reproducibility |

**Choosing `n_clusters`**: Use 5--10 for panels with few distinct populations. Use 10--20 for complex panels. Overclustering is generally safer than underclustering because it produces more homogeneous groups for quantile mapping.

**Choosing `goal`**: Use "mean" to normalize toward the average of all batches. Specify a batch name to normalize all other batches toward one reference batch (useful when one batch is considered the "gold standard").

---

## Pipeline Integration

**Script**: `src/harmonization/cytonorm.R`

**Function**: `run_cytonorm(train_frames, validation_frames, batch_labels, markers, ...)`

**Input from**:
- `load_fcs_raw()` in `src/loading/load_fcs_raw.R`
- `apply_compensation()` in `src/preprocessing/compensation.R`

**Output to**:
- `transform_flowframe()` in `src/preprocessing/transformation.R`
- `run_leiden()` in `src/clustering/leiden.R`

**Returns**: List containing `corrected_frames` (normalized flowFrame list) and `model` (fitted CytoNorm model for reuse).

**Placement in pipeline**: CytoNorm operates on compensated but pre-transformation data. Apply before arcsinh/logicle transformation. The sequence is: Load FCS -> Compensate -> CytoNorm -> Transform -> Downstream analysis.

**Used in**: Pipeline A when reference samples are available across batches.

---

## CytoNorm vs Harmony

| Aspect | CytoNorm | Harmony |
|--------|----------|---------|
| Reference samples required | Yes (essential) | No (batch labels only) |
| Correction basis | Quantile mapping from known ground truth | Iterative soft clustering + centroid alignment |
| Cell-type awareness | FlowSOM-stratified | Soft-cluster-stratified |
| Operates on | Raw/compensated intensities (pre-transform) | Transformed expression matrices (post-transform) |
| Output | Corrected flowFrames | Corrected numeric matrix |
| Best for | Studies with bead standards or technical replicates | Studies without shared reference samples |

---

## References

1. Van Gassen, S., Gaudilliere, B., Angst, M. S. et al. (2020). CytoNorm: A normalization algorithm for cytometry data. *Cytometry Part A*, 97(3), 268-278. [doi:10.1002/cyto.a.23904](https://doi.org/10.1002/cyto.a.23904)

2. Van Gassen, S., Callebaut, B., Van Helden, M. J. et al. (2015). FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data. *Cytometry Part A*, 87(7), 636-645. [doi:10.1002/cyto.a.22625](https://doi.org/10.1002/cyto.a.22625)

3. Trussart, M., Teh, C. E., Tan, T. et al. (2020). Removing unwanted variation with CytofRUV to integrate multiple CyTOF datasets. *eLife*, 9, e59630. [doi:10.7554/eLife.59630](https://doi.org/10.7554/eLife.59630)
