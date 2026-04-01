# Automated Quality Control for Flow Cytometry Data

## Summary

Automated quality control removes aberrant events caused by instrument
instabilities -- flow rate fluctuations, air bubbles, electronic noise, and
detector saturation -- before any downstream analysis. This pipeline supports
two complementary algorithms:
[flowAI](https://doi.org/10.1093/bioinformatics/btw191), which independently
evaluates three acquisition dimensions (flow rate, signal acquisition, and
dynamic range), and
[PeacoQC](https://doi.org/10.1002/cyto.a.24501), which uses
isolation trees and kernel density estimation to detect outlier events across
multiple channels simultaneously. Both are implemented in
`src/preprocessing/quality_control.R` and are called automatically by the
preprocessing orchestrator.

---

## Why QC Matters

Flow cytometers are complex fluidic and optical instruments. During a typical
acquisition session -- which may run for hours across dozens of samples -- a
variety of technical artifacts can corrupt the recorded data.

### Sources of Artifact

**Flow rate instabilities.** The fluidic system that delivers cells to the
interrogation point is susceptible to
[pressure fluctuations, partial clogs, and air bubble introduction](https://doi.org/10.1093/bioinformatics/btw191).
When the flow rate changes abruptly, cells pass through the laser at
irregular intervals, producing time-dependent shifts in scatter and
fluorescence intensities that do not reflect biology. Monaco et al. (2016)
demonstrated that
[flow rate anomalies are the most common acquisition artifact in large-scale cytometry studies](https://doi.org/10.1093/bioinformatics/btw191).

**Air bubbles and clogs.** Introduction of air into the sample line produces
transient bursts of high-scatter events followed by periods of abnormally low
event rates. These
[manifest as sharp spikes or gaps in the time-versus-event-rate profile](https://doi.org/10.1093/bioinformatics/btw191)
and contaminate both scatter and fluorescence channels.

**Electronic noise and detector saturation.** Photomultiplier tube (PMT)
voltage fluctuations, laser power drift, and detector saturation produce
signal anomalies that appear as
[sudden shifts in channel medians or increases in coefficient of variation over time](https://doi.org/10.1093/bioinformatics/btw191).
These events carry no biological information but can distort population
boundaries.

**Aggregate events (doublets and debris).** Cell aggregates and debris
particles produce anomalous scatter profiles. While doublet discrimination is
typically handled by scatter-based gating (FSC-H vs. FSC-A), automated QC
can flag time windows where
[aggregate frequency increases due to sample degradation or concentration changes](https://doi.org/10.1002/cyto.a.24501).

### Impact on Downstream Analysis

Aberrant events propagate through every step of the analysis pipeline.
[Transformation algorithms](arcsinh-transform.md) applied to corrupted data
produce distorted distributions.
[Leiden clustering](leiden-clustering.md) can assign artifact events to their
own spurious clusters or absorb them into genuine populations, shifting
cluster centroids. Differential abundance tests then report false positives
driven by technical variation rather than biology.
[Aghaeepour et al. (2013)](https://doi.org/10.1038/nmeth.2365) showed in
the FlowCAP-III challenge that
[automated gating and clustering performance degrades substantially when QC is omitted](https://doi.org/10.1038/nmeth.2365),
and recommended automated QC as a mandatory preprocessing step.

The
[Minimum Information about a Flow Cytometry Experiment (MIFlowCyt) guidelines](https://doi.org/10.1002/cyto.a.20623)
explicitly require reporting of quality control steps applied to cytometry
data, underscoring the community consensus that QC is not optional.

---

## flowAI

### Overview

[flowAI](https://doi.org/10.1093/bioinformatics/btw191) (Monaco et al.,
2016) is an R/Bioconductor package that performs fully automated quality
control by independently evaluating three dimensions of the acquisition
process: flow rate stability, signal acquisition consistency, and dynamic
range conformity. Events flagged as anomalous in any dimension are removed
from the flowFrame.

### Three-Dimension Anomaly Detection

flowAI decomposes the QC problem into three orthogonal checks, each targeting
a different class of instrument artifact.

#### 1. Flow Rate Check

The flow rate check detects time intervals where the event acquisition rate
deviates from the expected stable rate. flowAI
[divides the acquisition timeline into fixed-width bins and counts events per bin](https://doi.org/10.1093/bioinformatics/btw191).
It then applies a changepoint detection algorithm to identify abrupt shifts
in the event rate time series. Bins containing
[statistically significant rate changes (indicating clogs, bubbles, or pressure drops)](https://doi.org/10.1093/bioinformatics/btw191)
are flagged, and all events within those bins are marked for removal.

The changepoint detection uses the
[PELT (Pruned Exact Linear Time) algorithm](https://doi.org/10.1093/bioinformatics/btw191),
which identifies the optimal number and location of changepoints by
minimizing a penalized cost function. This approach is computationally
efficient (linear in the number of time bins) and
[avoids the need for user-specified thresholds](https://doi.org/10.1093/bioinformatics/btw191).

#### 2. Signal Acquisition Check

The signal acquisition check evaluates whether the fluorescence and scatter
signals remain stable over the course of acquisition. flowAI
[computes a moving median for each channel across time](https://doi.org/10.1093/bioinformatics/btw191)
and flags intervals where the local median deviates significantly from the
global channel median. This detects
[PMT voltage drift, laser power fluctuations, and gradual clogging](https://doi.org/10.1093/bioinformatics/btw191)
that cause progressive shifts in signal intensity.

The detection uses
[kernel density estimation on the time-binned channel medians](https://doi.org/10.1093/bioinformatics/btw191)
to identify outlier time intervals. Events within flagged intervals are
removed regardless of their individual intensity values, since the entire
interval is considered compromised.

#### 3. Dynamic Range Check

The dynamic range check identifies events that fall outside the expected
measurement boundaries for each channel. flowAI
[flags events at or near the instrument detection limits](https://doi.org/10.1093/bioinformatics/btw191)
-- both floor (negative saturation) and ceiling (positive saturation of the
analog-to-digital converter). These margin events carry no quantitative
information because the true signal has been clipped by the detector.

### Implementation in This Pipeline

```r
# From src/preprocessing/quality_control.R
qc_ff <- flowAI::flow_auto_qc(
  flowframe,
  folder_results = file.path(output_dir, paste0("FlowAI_", sample_name))
)
```

`flow_auto_qc()` runs all three checks sequentially, produces a diagnostic
report in `folder_results`, and returns a filtered flowFrame containing only
the events that passed all three checks. The diagnostic report includes
time-series plots of flow rate and channel medians with flagged intervals
highlighted.

### Strengths

- **Fully automated**: No user parameters required; thresholds are determined
  from the data via changepoint detection and density estimation.
- **Three independent dimensions**: Each check targets a different artifact
  class, providing comprehensive coverage.
- **Detailed diagnostics**: Produces per-sample HTML reports showing exactly
  which time intervals were flagged and why.
- [Validated on diverse instrument platforms](https://doi.org/10.1093/bioinformatics/btw191)
  including BD FACSCanto, BD LSRFortessa, and Beckman Coulter CytoFLEX.

### Limitations

- **Time-centric approach**: flowAI's core strategy is temporal -- it
  identifies bad time windows. Events with aberrant signal values that occur
  during otherwise normal time windows may be missed.
- **Single-channel evaluation**: Signal acquisition and dynamic range checks
  are performed per channel independently. Multivariate anomalies (events that
  are individually normal in each channel but jointly implausible) are not
  detected.
- **Runtime**: For large FCS files (> 1 million events), the per-channel
  density estimation can be
  [computationally expensive](https://doi.org/10.1093/bioinformatics/btw191),
  though still feasible on modern hardware.

---

## PeacoQC

### Overview

[PeacoQC](https://doi.org/10.1002/cyto.a.24501) (Emmaneel et al., 2022) is
an R/Bioconductor package that implements peak-based quality control using
isolation trees and multivariate kernel density estimation. Unlike flowAI's
three-dimension sequential approach, PeacoQC evaluates events across
multiple channels simultaneously, identifying outliers that deviate from the
expected multivariate signal distribution.

### Algorithm Details

PeacoQC operates in multiple stages to detect and remove anomalous events.

#### Peak-Based Framework

The "peak-based" approach in PeacoQC refers to its
[identification of stable signal peaks in each channel over time](https://doi.org/10.1002/cyto.a.24501).
For each channel, PeacoQC divides the acquisition timeline into consecutive
bins and estimates the density distribution within each bin. The mode (peak)
of the distribution is tracked over time. When the peak position
[shifts beyond an expected tolerance](https://doi.org/10.1002/cyto.a.24501),
the corresponding time interval is flagged as containing anomalous events.
This peak-tracking approach is particularly sensitive to
[gradual drift artifacts that may not trigger abrupt changepoint detection](https://doi.org/10.1002/cyto.a.24501).

#### Isolation Trees

PeacoQC incorporates
[isolation forests (Liu, Ting, and Zhou, 2008)](https://doi.org/10.1109/ICDM.2008.17)
as a complementary anomaly detection mechanism. Isolation trees work by
recursively partitioning the data space with random splits. Anomalous events,
being sparse and distinct from the bulk population,
[require fewer random splits to isolate](https://doi.org/10.1002/cyto.a.24501)
and therefore have shorter average path lengths in the tree. PeacoQC uses
the isolation tree anomaly score to
[flag events that are multivariate outliers](https://doi.org/10.1002/cyto.a.24501)
-- events that may appear normal in any single channel but are jointly
implausible when all channels are considered together.

This multivariate capability is a key advantage over univariate approaches:
[doublet events, for instance, may have plausible scatter values and plausible fluorescence values individually, but the combination is anomalous](https://doi.org/10.1002/cyto.a.24501).

#### Multivariate Approach and the Channels Parameter

PeacoQC evaluates anomalies across a user-specified set of channels, passed
via the `channels` parameter. This allows the user to focus QC on the marker
channels most relevant to the analysis, excluding channels (such as time or
viability dyes) that may add noise to the outlier detection.

In this pipeline, the `marker_indices` argument to `quality_control()` maps
directly to PeacoQC's `channels` parameter:

```r
# From src/preprocessing/quality_control.R
qc_result <- PeacoQC::PeacoQC(
  flowframe,
  channels = marker_indices,
  plot = TRUE,
  save_fcs = !is.null(output_dir),
  determine_good_cells = "all",
  output_directory = output_dir %||% tempdir()
)
```

The `determine_good_cells = "all"` setting instructs PeacoQC to apply all
available detection methods (peak-based, isolation tree, and margin event
checks) and remove events flagged by any method.

#### Margin Event Detection

In addition to the peak-based and isolation tree methods, PeacoQC includes
a
[margin event filter that removes events at detector saturation boundaries](https://doi.org/10.1002/cyto.a.24501).
These are events where one or more channels record the maximum or minimum
possible value for the detector, indicating that the true signal was clipped.

### Diagnostic Outputs

When `plot = TRUE`, PeacoQC generates per-sample diagnostic plots showing:

- **Channel-by-time heatmaps**: Signal intensity distributions over
  acquisition time, with flagged intervals highlighted in red.
- **Removal summary**: Percentage of events removed by each detection method
  (peak-based, isolation tree, margin events).
- **Before/after density overlays**: Kernel density plots for each channel
  before and after QC, allowing visual verification that the QC did not
  remove genuine biological populations.

### Strengths

- **Multivariate detection**: Isolation trees evaluate events across all
  specified channels simultaneously, catching anomalies that
  [univariate methods miss](https://doi.org/10.1002/cyto.a.24501).
- **Peak-based sensitivity**: The peak-tracking approach detects gradual drift
  artifacts that
  [changepoint-based methods may overlook because the change is too slow to trigger a changepoint](https://doi.org/10.1002/cyto.a.24501).
- **User-configurable channels**: The `channels` parameter allows focusing
  QC on biologically relevant markers.
- **Combined detection**: The `determine_good_cells = "all"` mode applies
  multiple complementary methods, providing comprehensive artifact removal.
- [Designed specifically for high-dimensional cytometry panels](https://doi.org/10.1002/cyto.a.24501),
  including spectral flow cytometry with 30+ channels.

### Limitations

- **Channel selection matters**: If `marker_indices` is poorly chosen
  (e.g., including irrelevant channels or omitting important ones), the
  multivariate outlier detection may be suboptimal.
- **Computational cost**: The isolation tree component adds computational
  overhead compared to purely temporal methods, particularly for very large
  files.
- **Sensitivity tuning**: While PeacoQC's defaults are well-calibrated for
  most cytometry data, datasets with unusual noise characteristics may require
  adjustment of internal thresholds.

---

## flowAI vs PeacoQC

| Aspect | flowAI | PeacoQC |
|--------|--------|---------|
| **Publication** | [Monaco et al. 2016](https://doi.org/10.1093/bioinformatics/btw191) | [Emmaneel et al. 2022](https://doi.org/10.1002/cyto.a.24501) |
| **Detection strategy** | Three independent temporal checks (flow rate, signal, dynamic range) | Peak-based tracking + isolation trees + margin events |
| **Multivariate** | No -- each channel evaluated independently | Yes -- isolation trees evaluate channels jointly |
| **Temporal sensitivity** | Changepoint detection (abrupt shifts) | Peak tracking (gradual and abrupt drift) |
| **User parameters** | None required (fully automated) | `channels` (which markers to evaluate) |
| **Output** | Filtered flowFrame + HTML diagnostic report | Filtered flowFrame + channel-by-time heatmaps + density overlays |
| **Speed** | Moderate (per-channel density estimation) | Moderate (isolation tree + density estimation) |
| **Best for** | Detecting discrete acquisition failures (clogs, bubbles, laser drops) | Detecting multivariate anomalies and gradual instrument drift |
| **Margin event removal** | Yes (dynamic range check) | Yes (explicit margin event filter) |
| **Handles high-dimensional panels** | Yes, but per-channel only | [Yes, designed for 30+ channel panels](https://doi.org/10.1002/cyto.a.24501) |
| **Default in this pipeline** | No | **Yes** (recommended) |

### When to Prefer flowAI

- Your primary concern is **discrete acquisition failures** -- sudden clogs,
  air bubbles, or laser shutoffs that produce sharp temporal discontinuities.
- You want **fully automated QC** with zero parameter specification.
- Your panel has few channels (< 10) and multivariate outlier detection
  provides limited additional benefit.

### When to Prefer PeacoQC

- Your panel is **high-dimensional** (15+ channels) and you want
  multivariate outlier detection.
- You suspect **gradual instrument drift** that may not trigger abrupt
  changepoints.
- You want **fine-grained control** over which channels drive the QC
  decision via the `channels` parameter.
- Your data contains **multivariate anomalies** (events normal in each
  channel individually but jointly implausible).

---

## When Your Data Calls For It

Automated QC is recommended in the following situations.

### Every FCS Acquisition (Default Recommendation)

Instrument artifacts are
[unpredictable and may occur in any sample during any acquisition session](https://doi.org/10.1093/bioinformatics/btw191).
Even well-maintained instruments produce occasional flow rate disturbances or
transient signal anomalies. Running automated QC on every sample is
computationally inexpensive relative to the downstream analysis and provides
a safety net against artifacts that manual inspection might miss.

### Before Transformation and Clustering

QC should precede transformation (arcsinh, logicle) and clustering. Aberrant
events distort the
[density distributions that transformation algorithms rely on for parameter estimation](https://doi.org/10.1002/cyto.a.20258).
If an entire time interval has shifted signal intensities, transforming those
events alongside valid data produces a
[bimodal artifact in the transformed space that can generate spurious clusters](https://doi.org/10.12688/f1000research.11622.3).

### Data with Suspected Instrument Issues

If acquisition logs or operator notes indicate flow rate warnings, laser
realignment during the run, or sample tube changes, targeted QC on the
affected samples is essential. Both flowAI and PeacoQC can isolate the
specific time windows affected by these events.

### Large Multi-Center Studies

In multi-center studies where data is acquired on different instruments by
different operators,
[instrument-specific artifacts are more likely and more variable](https://doi.org/10.1038/nmeth.2365).
Automated QC standardizes the data cleaning process across sites, reducing
inter-center technical variation that could confound biological comparisons.

---

## When to Skip QC

There are legitimate scenarios where automated QC is unnecessary or
counterproductive.

### Pre-Cleaned or Pre-Gated Data

If the data has already been gated to remove debris, doublets, and dead cells
in FlowJo or another gating tool, the most egregious artifact events have
already been excluded. Running automated QC on top of pre-gated data is
unlikely to remove many additional events and may flag genuine rare
populations as outliers. Use `qc_type = "None"` in this case.

### Very Small Files Where Manual Inspection Is Feasible

For FCS files with fewer than a few thousand events, manual inspection of
time-vs-channel plots is straightforward and may be more reliable than
automated algorithms whose
[statistical power depends on having sufficient events to estimate distributions](https://doi.org/10.1093/bioinformatics/btw191).
Automated QC on very small files risks removing a disproportionate fraction
of events.

### Data Already QC'd by a Core Facility

If the cytometry core facility has already applied their own QC pipeline
(many core facilities run flowAI or equivalent checks before releasing data),
re-running QC is redundant. Confirm with the core facility what QC steps
were applied and skip accordingly.

### Pre-Processed FCS Files (Pipeline E)

[Pipeline E](../pipelines/pipeline-E-preprocessed-umap.md) is designed for
FCS files that have already been QC'd and transformed. It skips the entire
preprocessing stage, including QC. If your FCS files are Pipeline E inputs,
QC has already been handled upstream.

---

## How to Use in the Pipeline

### PeacoQC (Default, Recommended)

```r
source("src/preprocessing/quality_control.R")

qc_ff <- quality_control(
  flowframe = raw_ff,
  sample_name = "Sample_001",
  qc_type = "PeacoQC",
  marker_indices = c(1:15),
  output_dir = "output/01_qc/"
)

cat(sprintf("Removed %d of %d events (%.1f%%)\n",
            nrow(raw_ff) - nrow(qc_ff),
            nrow(raw_ff),
            100 * (1 - nrow(qc_ff) / nrow(raw_ff))))
```

The `marker_indices` parameter controls which channels PeacoQC uses for
multivariate outlier detection. Pass the indices of your phenotypic markers
(excluding time, event number, and viability channels). If `NULL`, all
columns are used.

### flowAI Alternative

```r
source("src/preprocessing/quality_control.R")

qc_ff <- quality_control(
  flowframe = raw_ff,
  sample_name = "Sample_001",
  qc_type = "flowAI",
  output_dir = "output/01_qc/"
)
```

flowAI does not use the `marker_indices` parameter -- it evaluates all
channels automatically. The diagnostic report is saved to
`output/01_qc/FlowAI_Sample_001/`.

### Skip QC

```r
source("src/preprocessing/quality_control.R")

qc_ff <- quality_control(
  flowframe = raw_ff,
  sample_name = "Sample_001",
  qc_type = "None"
)
# Returns the input flowFrame unchanged
```

### Within the Pipeline Orchestrator

When using `run_preprocessing()` from `src/preprocessing/preprocess_pipeline.R`,
QC is configured via the `params` list:

```r
source("src/preprocessing/preprocess_pipeline.R")

params <- list(
  pp_order            = "QT",
  qc_type             = "PeacoQC",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56"),
  offset_list         = rep(0, 5),
  width_list          = rep(5, 5),
  transform_type      = "arcsinh",
  quant_norm          = FALSE,
  scaling             = FALSE,
  save_pp_fcs         = TRUE,
  save_pp_csv         = FALSE
)

preprocessed <- run_preprocessing(flow_frames, params, output_dir = "output/")
```

---

## QC in the Pipeline Context

Quality control is the first processing step in Pipeline A. After QC, the
data flows through transformation, downsampling, and clustering. Each
downstream step assumes that the input has been cleaned of instrument
artifacts.

### Processing Order

```
Raw FCS files
  |
  v
quality_control()          <-- You are here
  |
  v
transform_flowframe()      --> Arcsinh Transform / Logicle Transform
  |
  v
downsample_flowframes()    --> Downsampling
  |
  v
run_leiden()               --> Leiden Clustering
  |
  v
Differential abundance, visualization
```

### Related Methods

- [Arcsinh Transform](arcsinh-transform.md) -- The standard CyTOF
  transformation, applied after QC in the default `pp_order = "QT"` mode.
- [Logicle Transform](logicle-transform.md) -- The standard fluorescence
  transform, also applied after QC.
- [Downsampling](downsampling.md) -- Reduces cell counts for computational
  tractability, applied after QC and transformation.
- [Leiden Clustering](leiden-clustering.md) -- Unsupervised population
  identification, applied to QC'd and transformed data.

---

## Interpreting QC Results

### Percentage of Events Removed

Both flowAI and PeacoQC report the number of events retained after QC. The
pipeline prints a summary for each sample:

```
Sample Sample_001 QC with PeacoQC.
After QC: 48523/50000 cells remained.
```

**Expected removal rates.** For well-behaved acquisitions on properly
maintained instruments, automated QC typically removes
[1--5% of events](https://doi.org/10.1093/bioinformatics/btw191). A removal
rate under 1% suggests the acquisition was very clean. A rate above 10%
warrants investigation -- either the instrument had significant issues during
that sample, or the QC algorithm may be overly aggressive for that particular
data.

**Alarm thresholds.** Consider the following guidelines:

| Removal Rate | Interpretation | Action |
|--------------|----------------|--------|
| < 1% | Very clean acquisition | No action needed |
| 1--5% | Normal; minor artifacts removed | No action needed |
| 5--10% | Moderate artifact load | Review QC diagnostic plots |
| 10--20% | Significant instrument issues | Review acquisition logs; consider re-acquiring |
| > 20% | Severe problems | Data may be unreliable; consult instrument operator |

### QC Diagnostic Reports

**flowAI** produces an HTML report in the `folder_results` directory for each
sample. The report contains:

- Time-vs-event-rate plot with changepoints marked
- Per-channel time-vs-median plots with flagged intervals shaded
- Summary table of events removed by each check (flow rate, signal, dynamic
  range)

**PeacoQC** produces diagnostic plots when `plot = TRUE`:

- Channel-by-time heatmaps showing signal stability
- Red-shaded regions indicating flagged time intervals
- Summary of events removed by each detection method

### Visual Verification

After QC, overlay the density distributions of key channels before and after
removal. The post-QC distributions should show:

- Tighter, more symmetric peaks for scatter channels (FSC, SSC)
- Removal of events at detector saturation boundaries (margin events)
- Preserved population structure -- the QC should not have removed an entire
  biological population

If a genuine population has been removed, the QC may be too aggressive.
Consider:

- Adjusting the `marker_indices` in PeacoQC to exclude channels that
  define the rare population
- Switching to flowAI, which is less likely to flag multivariate outliers
- Using `qc_type = "None"` and performing manual gating of artifacts

---

## QC Order: QT vs TQ

The `pp_order` parameter in the preprocessing pipeline controls whether QC
runs before transformation (`"QT"`) or after it (`"TQ"`).

### QT: Quality Control First, Then Transform (Default)

```r
params$pp_order <- "QT"
```

In QT mode, QC operates on raw, untransformed intensity values. This is the
default and recommended order because:

- **Margin events are defined on the raw scale.** Detector saturation
  boundaries (0 and 2^18 for an 18-bit ADC) are
  [meaningful on the raw scale but meaningless after transformation](https://doi.org/10.1093/bioinformatics/btw191).
  QC on the raw scale correctly identifies margin events; QC after
  transformation may miss them.
- **Flow rate and signal checks are scale-invariant.** Temporal anomaly
  detection works on relative changes, so it produces identical results
  regardless of whether the data has been transformed.
- **Fewer transformed artifacts.** If aberrant events are removed before
  transformation, the transformation algorithm estimates its parameters
  (e.g., logicle W parameter or arcsinh bimodality-optimized cofactor) from
  clean data, producing
  [more accurate parameter estimates](https://doi.org/10.1002/cyto.a.20258).

### TQ: Transform First, Then Quality Control

```r
params$pp_order <- "TQ"
```

In TQ mode, the data is first transformed (arcsinh, logicle, etc.) and then
passed to QC. This order is appropriate in specific scenarios:

- **PeacoQC with isolation trees on transformed data.** If you are using
  PeacoQC's multivariate isolation tree detection, the transformed scale may
  provide better separation between normal and anomalous events. On the raw
  scale, high-intensity markers dominate the distance calculations; after
  transformation, all markers contribute more equally to the isolation tree
  partitioning.
- **Logicle-transformed fluorescence data.** The logicle transform maps
  compensated fluorescence data (which may include negative values) onto a
  scale where
  [population structure is more clearly resolved](https://doi.org/10.1002/cyto.a.20258).
  Running PeacoQC on this resolved structure can improve its ability to
  distinguish genuine rare populations from artifact events.
- **Consistency with upstream analysis tools.** Some analysis platforms (e.g.,
  FlowJo) apply transformations before any further processing. If you want
  your QC results to be comparable with QC performed in such tools, TQ mode
  reproduces that order.

### Decision Guide

| Scenario | Recommended Order | Rationale |
|----------|-------------------|-----------|
| CyTOF with flowAI | QT | flowAI's dynamic range check needs raw scale |
| CyTOF with PeacoQC | QT | Standard; margin events detected correctly |
| Fluorescence with flowAI | QT | Dynamic range check needs raw scale |
| Fluorescence with PeacoQC | QT or TQ | TQ can improve multivariate detection on logicle-transformed data |
| High-dimensional panel (30+ channels) with PeacoQC | TQ | Transformation equalizes channel contributions to isolation trees |

---

## Implementation Details

### Source Code Location

- **QC function:** [`quality_control()`](/src/preprocessing/quality_control.R)
  in `src/preprocessing/quality_control.R`
- **Pipeline orchestrator:** [`run_preprocessing()`](/src/preprocessing/preprocess_pipeline.R)
  in `src/preprocessing/preprocess_pipeline.R`

### Function Signature

```r
quality_control(flowframe, sample_name,
                qc_type = "PeacoQC",
                marker_indices = NULL,
                output_dir = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | *required* | Input flowFrame object |
| `sample_name` | character | *required* | Sample identifier for reporting and output file naming |
| `qc_type` | character | `"PeacoQC"` | Algorithm: `"PeacoQC"`, `"flowAI"`, or `"None"` |
| `marker_indices` | integer vector | `NULL` | Column indices for PeacoQC channels; `NULL` uses all columns |
| `output_dir` | character | `NULL` | Directory for diagnostic reports; `NULL` uses `tempdir()` |

### Return Value

A filtered `flowFrame` object containing only the events that passed quality
control. The returned flowFrame has the same column structure as the input
but fewer rows.

### Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `flowCore` | >= 2.0 | flowFrame I/O and manipulation |
| `flowAI` | >= 1.18 | `flow_auto_qc()` function |
| `PeacoQC` | >= 1.4 | `PeacoQC()` function |

---

## References by Topic

All references are linked inline throughout this document. The key
publications are:

- **flowAI (automated QC):** [Monaco, G. et al. (2016)](https://doi.org/10.1093/bioinformatics/btw191). flowAI: automatic and interactive anomaly discerning tools for flow cytometry data. *Bioinformatics*, 32(16), 2473--2480.
- **PeacoQC (peak-based QC):** [Emmaneel, A. et al. (2022)](https://doi.org/10.1002/cyto.a.24501). PeacoQC: Peak-based selection of high quality cytometry data. *Cytometry Part A*, 101(4), 325--338.
- **FlowCAP QC benchmark:** [Aghaeepour, N. et al. (2013)](https://doi.org/10.1038/nmeth.2365). Critical assessment of automated flow cytometry data analysis techniques. *Nat. Methods*, 10(3), 228--238.
- **MIFlowCyt reporting guidelines:** [Lee, J. A. et al. (2008)](https://doi.org/10.1002/cyto.a.20623). MIFlowCyt: The minimum information about a flow cytometry experiment. *Cytometry Part A*, 73A(10), 926--930.
- **Logicle transform (referenced):** [Parks, D. R., Roederer, M., & Moore, W. A. (2006)](https://doi.org/10.1002/cyto.a.20258). A new "logicle" display method avoids deceptive effects of logarithmic scaling for low signals and compensated data. *Cytometry Part A*, 69(6), 541--551.
- **CyTOF analysis workflow (referenced):** [Nowicka, M. et al. (2017)](https://doi.org/10.12688/f1000research.11622.3). CyTOF workflow: Differential discovery in high-throughput high-dimensional cytometry datasets. *F1000Research*, 6, 748.
- **Isolation forests (algorithm):** [Liu, F. T., Ting, K. M., & Zhou, Z.-H. (2008)](https://doi.org/10.1109/ICDM.2008.17). Isolation forest. In *Proceedings of the 2008 Eighth IEEE International Conference on Data Mining*, 413--422.
