# Spectral Unmixing for Spectral Flow Cytometry

Decomposes raw spectral detector signals into individual fluorochrome contributions using non-negative least squares (NNLS) regression against reference spectra from single-stain controls, with autofluorescence subtraction.

---

## When to Use

- You have spectral cytometry data (e.g., Cytek Aurora, Sony ID7000) with full-spectrum detection.
- Your instrument provides raw spectral data rather than pre-unmixed conventional channel data.
- You have single-stain control samples for each fluorochrome in your panel.
- You have an unstained control for autofluorescence estimation.

Spectral unmixing replaces traditional compensation for spectral cytometers. If you have conventional flow cytometry data with discrete detectors, use standard compensation (`src/preprocessing/compensation.R`) instead.

---

## How It Works

Spectral unmixing decomposes each cell's measured spectrum into a weighted sum of known fluorochrome reference spectra plus autofluorescence.

The process begins by building a reference spectral matrix from single-stain controls. For each fluorochrome, the median emission spectrum is extracted from the positive population of its single-stain control, and the autofluorescence baseline (median spectrum of the unstained control) is subtracted. This produces a matrix where each column represents the "pure" spectral signature of one fluorochrome across all detector channels. The autofluorescence spectrum itself is retained as an additional component to be fitted.

For each cell in the experimental data, the measured spectrum (a vector of intensities across all detector channels) is modeled as a linear combination of the reference spectra. The weights (abundance of each fluorochrome) are estimated by non-negative least squares (NNLS), which solves the system `measured = reference_matrix * weights` subject to the constraint that all weights are non-negative. The non-negativity constraint is physically meaningful: a cell cannot have a negative amount of any fluorochrome. NNLS avoids the negative spillover artifacts that can occur with ordinary least squares or matrix inversion approaches used in traditional compensation.

The quality of unmixing depends critically on the accuracy of the reference spectra. Instrument drift between control and sample acquisition, differences in cell autofluorescence between control and sample cell types, and insufficient spectral separation between fluorochromes can all degrade unmixing quality. The residual (difference between measured and reconstructed spectrum) for each cell provides a quality metric: large residuals indicate poor unmixing, possibly due to an unexpected fluorescent component not included in the reference matrix.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `control_frames` | (required) | Named list of single-stain flowFrames (names = fluorochrome labels) |
| `unstained_frame` | (required) | Unstained control flowFrame for autofluorescence baseline |
| `channels` | NULL (auto-detect) | Character vector of spectral channel names |
| `method` | "nnls" | Unmixing method: "nnls" (non-negative least squares) or "ols" (ordinary least squares) |
| `include_autofluorescence` | TRUE | Whether to include autofluorescence as a fitted component |

**Reference spectra quality**: The single-stain controls must be acquired on the same instrument with the same settings as the experimental samples. Controls should have a clearly positive population for the target fluorochrome. Dim fluorochromes benefit from brighter single-stain controls (e.g., compensation beads) to improve the signal-to-noise ratio of the reference spectrum.

---

## Pipeline Integration

**Script**: `src/preprocessing/spectral_unmixing.R`

**Functions**:
- `build_unmixing_matrix(control_frames, unstained_frame, channels)` -- constructs reference matrix
- `unmix_spectral(flow_frames, unmix_matrix, autofluorescence)` -- applies unmixing to samples

**Input from**:
- `load_fcs_raw()` in `src/loading/load_fcs_raw.R` -- loads single-stain controls and experimental FCS files

**Output to**:
- `transform_flowframe()` in `src/preprocessing/transformation.R` -- arcsinh/logicle transformation after unmixing
- `run_quality_control()` in `src/preprocessing/quality_control.R` -- QC on unmixed data

**Placement in pipeline**: Spectral unmixing replaces the compensation step and should be applied before transformation. The sequence is: Load FCS -> Spectral Unmixing -> QC -> Transform -> Downstream analysis.

**Used in**: Pipeline A when working with spectral cytometry instruments.

---

## Spectral Unmixing vs Compensation

| Aspect | Spectral Unmixing (NNLS) | Traditional Compensation |
|--------|--------------------------|--------------------------|
| Instrument type | Spectral cytometers (full-spectrum) | Conventional cytometers (discrete detectors) |
| Approach | Regression against reference spectra | Matrix inversion of spillover matrix |
| Negative values | Prevented by NNLS constraint | Common in compensated data |
| Fluorochrome capacity | Supports more fluorochromes than detectors | Limited to approximately one fluorochrome per detector |
| Autofluorescence | Explicitly modeled and subtracted | Not directly addressed |

---

## References

1. Novo, D., Griesemer, S. B. & Ambrozak, D. (2013). Full spectrum flow cytometry. In *Flow Cytometry Protocols*, Methods in Molecular Biology, vol 1678. [doi:10.1007/978-1-4939-7346-0_1](https://doi.org/10.1007/978-1-4939-7346-0_1)

2. Nolan, J. P. & Condello, D. (2013). Spectral flow cytometry. *Current Protocols in Cytometry*, 63, 1.27.1-1.27.13. [doi:10.1002/0471142956.cy0127s63](https://doi.org/10.1002/0471142956.cy0127s63)

3. Verma, A. & Bhatt, A. N. (2021). Spectral flow cytometry: Advantages and practical considerations. *Methods in Molecular Biology*, 2386, 39-55. [doi:10.1007/978-1-0716-1771-7_3](https://doi.org/10.1007/978-1-0716-1771-7_3)
