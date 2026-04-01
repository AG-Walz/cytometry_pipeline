# Logicle (Biexponential) Transform for Fluorescence Flow Cytometry

## Summary

The logicle transform is a [parameterized biexponential function](https://doi.org/10.1002/cyto.a.20258)
that maps compensated fluorescence cytometry data onto a display scale that is
logarithmic at high intensities, linear near zero, and capable of representing
negative values arising from spectral compensation. It is the
[industry-standard display and analysis transform](https://doi.org/10.1002/cyto.a.22030)
for multi-color fluorescence flow cytometry, implemented in FlowJo, flowCore,
and essentially all modern cytometry software. In this pipeline, it is applied
via `transform_flowframe()` in `src/preprocessing/transformation.R` using
`flowCore::logicleTransform()` with auto-estimated parameters.

---

## The Method

### The Problem the Logicle Solves

Traditional logarithmic display of fluorescence data fails when values approach
or cross zero. In multi-color flow cytometry, fluorescence compensation
subtracts spectral spillover between channels, routinely producing
[negative intensity values for cells with low or no true expression](https://doi.org/10.1002/cyto.a.20258).
A log scale cannot represent these values at all -- they are either clipped to
zero or discarded, creating artificial "piling up" of events along axes that
[distorts population distributions and obscures dim populations](https://doi.org/10.1002/cyto.a.20258).
Parks, Roederer, and Moore (2006) demonstrated that this display artifact leads
to incorrect gating decisions and misinterpretation of compensation quality.

### The Biexponential Function

The logicle is a specific parameterization of the
[generalized biexponential function](https://doi.org/10.1002/cyto.a.20258):

$$
S(x, T, W, M, A) = T \cdot 10^{-(M-W)} \left( 10^{(x-W)} - p^2 \cdot 10^{-(x-W)/p} + p^2 - 1 \right)
$$

where the parameter *p* is determined implicitly by the constraint that the
display function has zero slope at the negative boundary. In practice, the
transform is defined as the inverse of a display function that maps data space
to display space (typically [0, 1] for a normalized axis). The key behavior is:

- **For large positive values**: the function is approximately logarithmic
  (base 10), compressing the wide dynamic range of fluorescence intensity
  (often 4-5 decades) into a viewable scale.
- **Near zero**: the function transitions smoothly to a
  [linear mapping that passes through zero](https://doi.org/10.1002/cyto.a.20258),
  preventing the compression and distortion that log scaling causes in this
  critical region.
- **For negative values**: the function extends symmetrically below zero,
  allowing compensated negative events to be displayed in their true position
  rather than being clipped.

[Moore and Parks (2012)](https://doi.org/10.1002/cyto.a.22030) published an
updated formulation with numerically stable reference implementations,
including the root-finding algorithm used by flowCore, and clarified the
relationship between the logicle and the earlier "biexponential" parameterization
used in FlowJo.

### The Four Parameters: T, W, M, A

The logicle transform is fully determined by four parameters:

| Parameter | Name | Meaning | Typical Value |
|-----------|------|---------|---------------|
| **T** | Top of scale | Maximum data value (instrument ceiling) | 2^18 = 262,144 |
| **W** | Width | Width of the linearization region in asymptotic decades; controls how far the linear portion extends | Auto-estimated from data |
| **M** | Decades | Total number of decades of data range the display covers | 4.5 |
| **A** | Additional negative decades | Extra decades below zero for displaying negative events | 0 (auto-estimated) |

**T** is determined by the instrument digitizer. Modern digital cytometers
typically have an 18-bit range, giving T = 262,144. This value is
[stored in the FCS file $PnR keyword](https://doi.org/10.1002/cyto.a.20258) and
can be read directly from the flowFrame.

**W** is the most important parameter for the shape of the linear region. A
larger W produces a wider linear zone around zero, which is necessary when
compensation introduces [substantial negative spread](https://doi.org/10.1002/cyto.a.20258).
Parks et al. (2006) showed that W should be set based on the width of the
negative population: specifically, W is estimated as the number of decades
spanned by the negative data, ensuring that the linear region is wide enough
to display compensated negative events without distortion.

**M** sets the total number of logarithmic decades displayed. The
[standard value of 4.5](https://doi.org/10.1002/cyto.a.22030) covers the
typical dynamic range of fluorescence cytometers. Reducing M compresses the
display; increasing it spreads populations further apart at the expense of
resolution in the low-intensity region.

**A** adds extra decades below zero. When A = 0 (the default),
the negative range is determined entirely by W. Setting A > 0 extends the
display further into negative territory, which can be useful for channels with
[very broad compensation-induced spread](https://doi.org/10.1002/cyto.a.22030).

### How flowCore Auto-Estimates Parameters

When `flowCore::logicleTransform()` is called without explicit parameters (as
in this pipeline), flowCore uses the
[automatic parameter estimation procedure described by Parks et al. (2006)](https://doi.org/10.1002/cyto.a.20258):

1. **T** is read from the flowFrame's `$PnR` keyword (channel range maximum).
2. **M** defaults to 4.5 decades.
3. **W** is estimated from the data by finding the width of the negative
   population in the channel. Specifically, it uses the robust second-percentile
   negative value to determine how many decades the negative tail spans, then
   sets W to accommodate that spread within the linear portion.
4. **A** defaults to 0.

This auto-estimation works well for standard compensated fluorescence data but
can fail when data distributions are unusual (see Troubleshooting below).

---

## Why Logicle for Fluorescence

### The Compensation Problem

Multi-color fluorescence flow cytometry requires
[spectral compensation to correct for emission overlap between fluorophores](https://doi.org/10.1002/cyto.a.20258).
Compensation is a linear subtraction: for each cell, the intensity measured in a
channel is reduced by a fraction of the intensity in every overlapping channel.
When a cell has high intensity in one channel and low or no true signal in an
adjacent channel, the subtraction can produce a negative compensated value.
These negatives are not errors -- they are a
[mathematical consequence of correct compensation](https://doi.org/10.1002/cyto.a.20258)
and occur most frequently in channels with high spillover coefficients.

On a logarithmic scale, these negative events are undefined and either vanish or
pile up at zero, creating an artificial wall that
[obscures the true shape of dim populations](https://doi.org/10.1002/cyto.a.20258)
and makes it impossible to judge whether compensation has been applied correctly.
The logicle transform was specifically invented to solve this problem by
providing a continuous display from negative through zero to positive values.

### Autofluorescence Near Zero

Even in perfectly compensated data, cells have intrinsic autofluorescence that
varies from cell to cell. This creates a spread of events centered near zero
for unstained or dimly stained populations. On a log scale, this spread is
[asymmetrically compressed and distorted](https://doi.org/10.1002/cyto.a.20258),
making dim-positive populations appear to merge with the negative population.
The logicle's linear region near zero preserves the true symmetric spread of
these events, allowing proper resolution of dim populations from negative
controls.

### FlowJo Compatibility

FlowJo (the dominant commercial cytometry analysis software) adopted the
biexponential/logicle display as its default for fluorescence data based on
the [Parks et al. (2006)](https://doi.org/10.1002/cyto.a.20258) publication.
Using the logicle transform in this pipeline ensures that visualizations and
gate positions are directly comparable to FlowJo output, which is critical for
[reproducibility and validation of computational gating against manual expert gating](https://doi.org/10.1002/cyto.a.22030).

### ISAC Recommendations

The International Society for Advancement of Cytometry (ISAC) has recommended
the logicle display for compensated fluorescence data in its
[guidelines for the use of flow cytometry](https://doi.org/10.1002/cyto.a.20258).
The recommendation is based on Parks et al.'s demonstration that the logicle
avoids the display artifacts introduced by logarithmic scaling and produces
plots that more accurately represent the underlying biology.

---

## Parameters

### T (Top of Scale)

The top-of-scale parameter corresponds to the maximum channel value of the
cytometer's analog-to-digital converter. For modern digital cytometers operating
at 18-bit resolution, T = 2^18 = 262,144. Older instruments with lower
resolution may use T = 2^14 = 16,384 or T = 2^12 = 4,096. The
[FCS file standard](https://doi.org/10.1002/cyto.a.20258) stores this value in
the `$PnR` keyword for each parameter, and flowCore reads it automatically.

Setting T incorrectly (e.g., using T = 262,144 for data acquired with a 12-bit
instrument) will produce a transform that leaves most of the display range
empty, compressing all events into a small portion of the scale.

### W (Width of Linear Region)

W controls the size of the linear zone around zero. It is expressed in
[asymptotic decades](https://doi.org/10.1002/cyto.a.20258) -- the number of
decades that the linear region spans when extrapolated to where the function
would converge to zero in the log domain.

- **Small W** (e.g., 0.1-0.3): Narrow linear region. Appropriate for channels
  with little compensation spillover and minimal negative spread.
- **Large W** (e.g., 0.5-2.0): Wide linear region. Necessary for channels with
  [heavy compensation where the negative population spans a wide range](https://doi.org/10.1002/cyto.a.20258).

flowCore estimates W by examining the negative portion of the data distribution.
If the negative tail is broad (indicating heavy compensation), W is set larger
to ensure the linear region encompasses the entire negative spread.

### M (Decades)

M defines the total number of logarithmic decades covered by the display. The
[standard value is 4.5](https://doi.org/10.1002/cyto.a.22030), which covers the
typical range of fluorescence intensities from noise floor to maximum signal on
modern instruments. Adjusting M is rarely necessary unless working with unusual
instrument configurations or extremely compressed or expanded data ranges.

### A (Additional Negative Decades)

A extends the display below the range that W already accommodates. With A = 0
(the default), the negative range is determined solely by W. Increasing A
provides additional space for visualizing
[extreme negative events](https://doi.org/10.1002/cyto.a.22030) that may appear
in channels with very high spillover coefficients. In most cases, A = 0 is
sufficient because W already accounts for the negative spread seen in typical
multi-color panels.

### Summary of Auto-Estimation Logic

```
flowCore::logicleTransform() auto-estimation:

1. T  <--  $PnR keyword from FCS file (or max channel value)
2. M  <--  4.5 (default)
3. W  <--  estimated from the width of the negative data tail
           (uses robust percentile-based estimation)
4. A  <--  0 (default; can be increased if negative range is insufficient)
```

---

## When Your Data Calls For It

The logicle transform is the appropriate choice when all of the following
conditions hold:

### Compensated Fluorescence Flow Cytometry

Your data was acquired on a fluorescence cytometer (not CyTOF) and has been
compensated for spectral overlap, either during acquisition or post-acquisition.
Compensation produces the characteristic
[negative value artifacts](https://doi.org/10.1002/cyto.a.20258) that the
logicle is designed to handle.

### Data with Negative Values from Spectral Overlap

You observe negative intensity values in one or more channels after
compensation. These negatives are expected in any
[properly compensated multi-color panel](https://doi.org/10.1002/cyto.a.20258),
especially in channels with high spillover from neighboring fluorophores (e.g.,
FITC spilling into PE, or tandem dyes spilling into their parent channels).

### FlowJo-Compatible Display Desired

You want the transformed data to match the display scaling used in FlowJo,
enabling direct comparison of computational results with manual gating. FlowJo's
biexponential display is based on
[the same logicle function](https://doi.org/10.1002/cyto.a.22030) implemented
in flowCore, so results will be visually consistent.

### Preservation of Population Spacing Near Origin

You need to resolve dim-positive populations from negative populations near zero
intensity. The logicle's
[linear behavior in this region](https://doi.org/10.1002/cyto.a.20258) ensures
that the spread caused by autofluorescence and compensation is displayed
faithfully, without the compression artifacts of log scaling.

### Concrete Scenarios

1. **Standard immunophenotyping panel.** You are running a 10-color panel on a
   fluorescence cytometer, compensated with single-stain controls. Multiple
   channels show negative events. The logicle ensures proper display and
   analysis of all channels.

2. **High-spillover tandem dyes.** Your panel includes BV421, BV510, and BV605,
   which have significant spectral overlap. After compensation, BV510 has a
   broad negative tail. The logicle with auto-estimated W will accommodate this
   spread.

3. **Validating computational gating against FlowJo.** A collaborator gated the
   data in FlowJo using biexponential axes. To compare your Leiden clustering
   results, you apply the logicle transform so that overlaid plots use the same
   display scale.

---

## When to Choose Something Else

The logicle is not the right transform for every cytometry dataset. The table
below maps common situations to better alternatives.

| Situation | Problem with Logicle | Better Method |
|-----------|---------------------|---------------|
| CyTOF (mass cytometry) data | CyTOF does not use fluorescence compensation and lacks the negative artifacts the logicle was designed for. The [ion counting noise model is different](https://doi.org/10.1126/science.1198704) from photon-counting fluorescence. | [Arcsinh Transform](arcsinh-transform.md) with cofactor 5 |
| Count data (cells per cluster per sample) | Counts are non-negative integers, not continuous intensities. The logicle's negative-handling behavior is irrelevant. | [Log Transform](log-transform.md) (as offset in GLMMs) |
| Data already transformed | Applying the logicle to data that has already been transformed (e.g., FCS files exported from FlowJo with pre-applied biexponential scaling) produces nonsensical double-transformed values. | Skip transformation (`transform_type = "None"`) |
| Uncompensated fluorescence data | Without compensation, fluorescence data is strictly positive. The logicle still works but its negative-handling capability is unnecessary. | Logicle still works; arcsinh is also acceptable |
| Spectral cytometry (unmixed data) | Spectral unmixing produces negative values similar to compensation, but the noise characteristics may differ. | Logicle or arcsinh; evaluate empirically |

### Cross-Links

- For CyTOF data, see [Arcsinh Transform](arcsinh-transform.md).
- For count-based offsets in statistical models, see [Log Transform](log-transform.md).

---

## How to Use in the Pipeline

### Transform Fluorescence FCS Files (Pipeline A)

```r
# Source the transformation module
source("src/preprocessing/transformation.R")

# Apply logicle transform with auto-estimated parameters
transformed <- transform_flowframe(
  flowframe           = raw_flowframe,
  sample_name         = "patient_01",
  markers_of_interest = markers_of_interest,
  offset_list         = rep(0, length(markers_of_interest)),
  width_list          = rep(1, length(markers_of_interest)),
  transform_type      = "logicle"
)
```

The `offset_list` and `width_list` arguments are required by the function
signature but are **not used by the logicle path** -- `flowCore::logicleTransform()`
ignores them and auto-estimates its own parameters (T, W, M, A) from the data.
Pass dummy values (zeros or ones) for these arguments when using the logicle.

### Biexponential Variant

The biexponential transform (`flowCore::biexponentialTransform()`) is a closely
related generalization that uses a different parameterization. It is available
for compatibility with specific FlowJo versions:

```r
# Apply biexponential transform
transformed <- transform_flowframe(
  flowframe           = raw_flowframe,
  sample_name         = "patient_01",
  markers_of_interest = markers_of_interest,
  offset_list         = rep(0, length(markers_of_interest)),
  width_list          = rep(1, length(markers_of_interest)),
  transform_type      = "biexponential"
)
```

### Within the Full Pipeline A Orchestrator

In the standard Pipeline A workflow, transformation happens after quality
control. The config file controls the choice:

```yaml
# config/pipeline_A.yml (excerpt)
transform_type: "logicle"    # "arcsinh", "logicle", "biexponential", or "None"
pp_order: "QT"               # QC first, then transform
```

The pipeline orchestrator (`pipelines/pipeline_A.R`) reads this config and calls
`transform_flowframe()` on each flowFrame in the dataset.

### Implementation Details

The core of the implementation (from `src/preprocessing/transformation.R`):

```r
trans <- switch(transform_type,
  "arcsinh"        = flowCore::arcsinhTransform(a = -offset, b = 1 / width),
  "logicle"        = flowCore::logicleTransform(),
  "biexponential"  = flowCore::biexponentialTransform(),
  "None"           = NULL,
  stop("Invalid transform_type: '", transform_type, "'.")
)

if (!is.null(trans)) {
  tflist <- flowCore::transformList(from = marker, tfun = trans)
  tf_ff <- flowCore::transform(tf_ff, tflist, truncate_max_range = FALSE)
}
```

Key implementation choices:

- **`truncate_max_range = FALSE`**: Prevents flowCore from clipping values
  that exceed the instrument range after transformation. This is important
  because the logicle transform can produce values outside the nominal [0, 1]
  display range for extreme data values.
- **Per-marker application**: The transform is applied one marker at a time in a
  loop. Each marker gets its own auto-estimated parameters, because different
  channels have different compensation levels and therefore different optimal
  W values.
- **Missing marker handling**: If a marker is not found in the flowFrame
  column names, it is skipped with a warning rather than causing a fatal error.

---

## Effect on Downstream Analysis

### Impact on Leiden Clustering

The logicle transform affects [Leiden clustering](leiden-clustering.md) in two
important ways:

1. **Distance metric behavior.** Leiden clustering in this pipeline operates on
   a KNN graph built using Euclidean distance in marker space. On the logicle
   scale, distances between cells in the high-intensity regime are compressed
   (logarithmic behavior), while distances near zero are expanded (linear
   behavior). This means that the KNN graph gives
   [relatively more weight to differences among dim populations](https://doi.org/10.1002/cyto.a.20258)
   than a raw-intensity KNN graph would. This is generally desirable for
   immunophenotyping, where distinguishing dim-positive from negative populations
   is biologically critical.

2. **Negative value handling.** Without the logicle, negative compensated values
   would need to be either clipped (losing information) or left on a linear
   scale (where they dominate the dynamic range). The logicle's smooth handling
   of negatives ensures that [all cells contribute meaningfully to the neighbor graph](https://doi.org/10.12688/f1000research.11622.3)
   without distortion.

The [Nowicka et al. (2017) CyTOF workflow](https://doi.org/10.12688/f1000research.11622.3)
emphasizes that the choice of data transformation is a critical preprocessing
step that directly influences cluster quality and should be validated by
inspecting marker expression heatmaps after clustering.

### Impact on UMAP Visualization

[UMAP](umap.md) embeddings are computed from the same transformed marker space
used for clustering. The logicle transform therefore determines how populations
are spaced in the UMAP plot:

- Populations that differ in high-intensity markers will appear at
  log-proportional distances.
- Populations that differ in dim markers (near the linear region) will be
  [better separated than they would be on a raw or log scale](https://doi.org/10.1002/cyto.a.20258),
  because the linear region prevents compression.
- [Becht et al. (2019)](https://doi.org/10.1038/nbt.4314) noted that
  appropriate data transformation before UMAP is essential for producing
  embeddings that reflect biological rather than technical variation.

### Impact on Heatmaps and Expression Plots

Marker expression heatmaps use the transformed values directly. On the logicle
scale, expression values are approximately log-distributed for positive
populations, which produces
[visually balanced heatmaps](https://doi.org/10.12688/f1000research.11622.3)
where both dim and bright markers contribute visible contrast. Without
transformation, bright markers would dominate the color scale and dim markers
would be invisible.

---

## Comparison: Logicle vs Arcsinh

Both the logicle and the [arcsinh transform](arcsinh-transform.md) address the
same fundamental challenge: transforming data that spans a wide dynamic range
and includes values near or below zero. They differ in their design origin,
parameterization, and optimal use cases.

### Functional Comparison

| Property | Logicle | Arcsinh |
|----------|---------|---------|
| **Design origin** | [Fluorescence flow cytometry with compensation](https://doi.org/10.1002/cyto.a.20258) | [Mass cytometry (CyTOF) with ion counting](https://doi.org/10.1126/science.1198704) |
| **Behavior at high values** | Logarithmic (base 10) | Logarithmic (natural log) |
| **Behavior near zero** | Linear with controlled width (W parameter) | Linear with width controlled by cofactor |
| **Negative value handling** | [Explicitly parameterized (A parameter)](https://doi.org/10.1002/cyto.a.22030) | Symmetric by construction (arcsinh is an odd function) |
| **Parameters** | T, W, M, A (4 parameters, auto-estimated) | Cofactor only (1 parameter) |
| **Auto-estimation** | [Data-driven via flowCore](https://doi.org/10.1002/cyto.a.22030) | Fixed cofactor (5 for CyTOF) or optimized for bimodality |
| **FlowJo compatibility** | Direct match | Not used by FlowJo |
| **Computational cost** | Requires root-finding for inverse | Closed-form formula |

### When to Choose Each

**Choose logicle when:**
- You have [compensated fluorescence data](https://doi.org/10.1002/cyto.a.20258)
  with negative spillover artifacts.
- You need FlowJo-compatible display for validation or publication.
- The negative population spread varies by channel and needs per-channel
  parameter estimation.

**Choose arcsinh when:**
- You have [CyTOF data](https://doi.org/10.1126/science.1198704) where the
  noise model is based on ion counting rather than photon detection.
- You want a single, simple parameter (cofactor) rather than four parameters.
- You are following established CyTOF analysis workflows such as
  [Nowicka et al. (2017)](https://doi.org/10.12688/f1000research.11622.3) or
  [Levine et al. (2015)](https://doi.org/10.1016/j.cell.2015.05.047).
- You need a fast closed-form transform without root-finding.

**In practice**, the choice is often dictated by instrument type: logicle for
fluorescence cytometry, arcsinh for CyTOF. Using the wrong transform is not
catastrophic, but using the right one
[produces cleaner population separation and more interpretable visualizations](https://doi.org/10.1002/cyto.a.20258).

For cross-links, see [Arcsinh Transform](arcsinh-transform.md).

---

## Troubleshooting

### Auto-Estimation Fails with an Error

**Symptom:** `flowCore::logicleTransform()` throws an error such as
`"IllegalArgument: w is too large"` or `"f(a) and f(b) have the same sign"`.

**Cause:** The auto-estimation algorithm could not find valid W and A parameters
for the data in a channel. This typically happens when:

- The channel has no negative values (e.g., a scatter parameter or an
  uncompensated channel). The algorithm expects
  [negative events to estimate W](https://doi.org/10.1002/cyto.a.20258) and
  fails when none exist.
- The data range is extremely narrow or contains only a few unique values.
- The channel was pre-transformed or pre-scaled before loading.

**Solution:** Exclude non-fluorescence parameters (FSC, SSC, Time) from the
`markers_of_interest` list. Only pass compensated fluorescence channels to the
logicle transform:

```r
# Exclude scatter and time parameters
fluorescence_markers <- setdiff(
  markers_of_interest,
  c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "Time")
)

transformed <- transform_flowframe(
  flowframe           = raw_flowframe,
  sample_name         = "patient_01",
  markers_of_interest = fluorescence_markers,
  offset_list         = rep(0, length(fluorescence_markers)),
  width_list          = rep(1, length(fluorescence_markers)),
  transform_type      = "logicle"
)
```

### Transform Produces Unexpected Compression or Stretching

**Symptom:** After transformation, populations that should be well-separated
appear compressed together, or the display range seems wrong.

**Cause:** The auto-estimated T (top of scale) does not match the actual data
range. This happens when:

- The FCS file's `$PnR` keyword is set incorrectly (e.g., default 262,144 for
  a 12-bit instrument that only uses range 0-4,095).
- Data was acquired with gain settings that limit the effective range.

**Solution:** Use `flowCore::estimateLogicle()` to inspect the auto-estimated
parameters, then construct a custom transform with corrected values:

```r
# Inspect auto-estimated parameters
el <- flowCore::estimateLogicle(raw_flowframe, channels = "FL1-A")
print(el)

# Construct with manual parameters if needed
custom_logicle <- flowCore::logicleTransform(
  w = 0.5,    # Adjust width
  t = 262144, # Correct top of scale
  m = 4.5,    # Standard decades
  a = 0       # No extra negative decades
)
```

### Some Channels Have Very Wide Negative Spread

**Symptom:** A channel with high spillover shows a very broad negative
population that extends far below zero, and the auto-estimated transform
seems to linearize too much of the range.

**Cause:** The W parameter is set very large to accommodate the broad negative
tail, which extends the linear region and compresses the logarithmic region.
This is the [correct behavior for heavily compensated channels](https://doi.org/10.1002/cyto.a.20258),
but it can make the display look unusual compared to channels with less
spillover.

**Solution:** This is usually not a problem -- the transform is doing what it
should. If the display seems too compressed, verify that compensation values
are correct. Over-compensation (spillover coefficients too high) artificially
broadens the negative population. If compensation is correct and the spread is
genuine, consider increasing A to extend the negative display range:

```r
custom_logicle <- flowCore::logicleTransform(
  w = 1.0,  # Wide linear region for heavy compensation
  t = 262144,
  m = 4.5,
  a = 0.5   # Extra negative decades
)
```

### Logicle vs Biexponential: Which to Use

**Symptom:** Uncertainty about whether to use `transform_type = "logicle"` or
`transform_type = "biexponential"`.

**Explanation:** In this pipeline, both call flowCore functions with default
parameters. The
[logicle is a specific constrained parameterization of the general biexponential](https://doi.org/10.1002/cyto.a.22030).
The key differences:

- `flowCore::logicleTransform()` uses the [Parks et al. (2006)](https://doi.org/10.1002/cyto.a.20258)
  parameterization with auto-estimation. This is the recommended default.
- `flowCore::biexponentialTransform()` uses a more general parameterization with
  different defaults. It may not auto-estimate from data as robustly.

**Recommendation:** Use `"logicle"` unless you have a specific reason to prefer
the biexponential parameterization (e.g., matching a particular FlowJo version's
output).

### Performance with Large Datasets

**Symptom:** Transformation is slow when applied to millions of events.

**Cause:** The logicle transform involves numerical root-finding (inverting the
biexponential function) for each data point, which is slower than the
closed-form arcsinh transform.

**Solution:** The per-event computation is still fast (microseconds), but for
datasets with millions of events across many markers, total time can add up.
Consider:

1. **Downsampling before transformation** using
   `src/preprocessing/downsampling.R` if the full dataset is not needed.
2. **Parallel processing** across samples (each flowFrame can be transformed
   independently).
3. For exploration only, using the arcsinh transform for speed and switching to
   logicle for final analysis and publication-quality figures.

---

## Source Code Location

- **Transform function:** [`transform_flowframe()`](/src/preprocessing/transformation.R)
  in `src/preprocessing/transformation.R`
- **Pipeline orchestrator:** `pipelines/pipeline_A.R` (calls transform as part
  of the preprocessing step)
- **Config template:** `config/pipeline_A.yml` (set `transform_type: "logicle"`)
- **Shared utilities:** `src/R_shared/flowframe_utils.R` (flowFrame helper
  functions used during transformation)
