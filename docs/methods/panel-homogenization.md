# Panel Comparison and Homogenization

## Summary

Panel homogenization is the process of aligning marker panels across flow cytometry samples so that every sample presents the same set of markers in the same column order, enabling valid cross-sample comparison. In multi-panel and multi-site studies, different antibody configurations, naming conventions, and panel evolution over time produce samples with [heterogeneous marker sets that cannot be directly concatenated or jointly analyzed](https://doi.org/10.12688/f1000research.11622.3). This pipeline provides `compare_panels_detailed()` to diagnose panel differences and `homogenize_flowframes()` to resolve them through marker renaming, reordering, and NA-filling, producing a unified set of flowFrame objects ready for downstream [batch correction](batch-correction.md), [clustering](leiden-clustering.md), and visualization.

---

## The Panel Mismatch Problem

### Why Panels Differ

A "panel" in flow cytometry is the specific combination of antibody-fluorochrome conjugates used to stain a sample. In practice, panels vary across samples for several reasons:

- **Multi-site studies**: Consortium studies such as the [Human Immunology Project Consortium (HIPC)](https://doi.org/10.1038/ni.2549) and the [ONE Study](https://doi.org/10.1016/j.trim.2012.11.001) deliberately or inadvertently use different antibody clones, fluorochrome conjugates, or instrument configurations across sites. [Maecker et al. (2012)](https://doi.org/10.1038/nri3158) documented that even when sites attempt to use the same panel, differences in instrument setup, reagent lots, and operator protocols produce effectively different panels.
- **Panel evolution**: Longitudinal studies that span months or years frequently modify panels as [new reagents become available, markers are added for emerging hypotheses, or discontinued antibodies are replaced](https://doi.org/10.1002/cyto.a.22624). This creates a time-dependent panel structure where early samples lack markers that later samples include, and vice versa.
- **Instrument differences**: Different cytometer models support different numbers of detectors, laser lines, and fluorochrome-detector assignments. A panel designed for a [5-laser instrument cannot be run identically on a 3-laser instrument](https://doi.org/10.1002/cyto.a.22624), forcing modifications that alter the marker set.
- **Naming conventions**: The same biological marker may be recorded under different names across files. A CD3 channel might appear as `CD3`, `CD3-FITC`, `CD3_FITC`, `FITC-A`, or `<FITC-A> CD3` depending on the instrument software, export settings, and operator conventions. [Spidlen et al. (2021)](https://doi.org/10.1002/cyto.a.24225) highlighted that inconsistent marker naming in FCS files is one of the most common obstacles to automated multi-file analysis.

### Consequences of Ignoring Panel Differences

Attempting to concatenate or jointly analyze samples with mismatched panels produces several failure modes:

1. **Column misalignment**: If sample A has markers [CD3, CD4, CD8] and sample B has [CD3, CD8, CD19], naive row-binding produces a matrix where columns represent different markers for different rows, corrupting all downstream analysis.
2. **Silent data loss**: Subsetting to the intersection of markers without verifying that the intersection is biologically meaningful can discard markers that define key populations.
3. **Dimension errors**: Algorithms that expect a fixed number of features (KNN graph construction for [Leiden clustering](leiden-clustering.md), [UMAP](umap.md) embedding, [PCA](pca.md)) will fail or produce meaningless results if the input matrix has inconsistent column counts across samples.
4. **False biological conclusions**: If panel differences correlate with experimental conditions (e.g., all treatment samples were run with panel version 2, all controls with version 1), unresolved panel effects can [masquerade as biological differences in differential abundance analysis](https://doi.org/10.12688/f1000research.11622.3).

### The Scale of the Problem

In a [systematic evaluation of multi-center cytometry studies, Finak et al. (2016)](https://doi.org/10.1002/cyto.a.22859) found that center-to-center variation in panel composition, instrument settings, and gating strategies was a dominant source of variability, often exceeding the biological signal of interest. The [EuroFlow Consortium](https://doi.org/10.1038/leu.2012.122) addressed this by mandating standardized 8-color panels with fixed antibody-fluorochrome combinations across all participating sites, demonstrating that [enforced panel standardization dramatically improved cross-site reproducibility](https://doi.org/10.1038/leu.2012.122). When such prospective standardization is not possible -- which is the common situation in retrospective analyses, meta-analyses, and studies combining legacy data -- computational panel homogenization is the necessary alternative.

---

## Panel Comparison

### Purpose

Before homogenizing panels, you must first understand how they differ. `compare_panels_detailed()` provides a complete diagnostic report of panel heterogeneity across a collection of flowFrame objects. This is the essential first step: without knowing which markers are shared, which are unique, and which samples belong to the same panel group, you cannot make informed decisions about how to homogenize.

### How `compare_panels_detailed()` Works

The function accepts a named list of flowFrame objects and performs four operations:

#### 1. Extract Markers from Each Sample

For each flowFrame, the function calls `explore_data()` (from `src/R_shared/flowframe_utils.R`) to retrieve the marker annotations stored in the FCS file's `$PnS` (parameter short name) keywords. These are the human-readable marker names (e.g., "CD3", "CD4") as distinct from the channel/dye names (e.g., "FITC-A", "PE-A") stored in `$PnN`. Only columns with non-empty marker annotations are included -- scatter parameters (FSC, SSC) and time channels are excluded because they lack `$PnS` entries.

```r
all_markers <- lapply(flow_frames, function(ff) {
  data <- explore_data(ff)
  data$markers[data$mcols]
})
```

#### 2. Compute Common and Unique Markers

The function computes the set intersection across all samples to find markers present in every sample (`common_markers`), and the set union to find all markers present in any sample (`all_markers`). The difference between these two sets (`not_common`) identifies markers that are present in some samples but not others.

```r
common <- Reduce(intersect, all_markers)
all_unique <- unique(unlist(all_markers))
not_common <- setdiff(all_unique, common)
```

#### 3. Build Presence Matrix

A Boolean matrix is constructed with markers as rows and samples as columns. Each cell indicates whether that marker is present in that sample. This matrix is the primary diagnostic output: it reveals the exact pattern of panel heterogeneity.

```r
presence <- matrix(FALSE, nrow = length(all_unique), ncol = length(flow_frames),
                   dimnames = list(all_unique, sample_names))
for (i in seq_along(all_markers)) {
  presence[all_markers[[i]], i] <- TRUE
}
```

For example, with 20 samples and 35 unique markers, the presence matrix might reveal that 28 markers are shared by all samples, 4 markers are present in only the first 10 samples (panel version 1), and 3 markers are present in only the last 10 samples (panel version 2).

#### 4. Group Samples by Panel

The function calls `compare_panels()` (from `src/R_shared/flowframe_utils.R`) to group samples that share an identical marker set. Two samples belong to the same panel group if and only if their marker vectors are identical (same markers in the same order). This grouping helps identify how many distinct panels exist in the dataset and which samples belong to each.

### Return Value

`compare_panels_detailed()` returns a list with six components:

| Component | Type | Description |
|-----------|------|-------------|
| `common_markers` | character vector | Markers present in every sample |
| `all_markers` | character vector | Union of all markers across all samples |
| `not_common` | character vector | Markers present in some but not all samples |
| `presence_matrix` | logical matrix | Samples (columns) x markers (rows) presence/absence |
| `panel_groups` | list | Groups of samples sharing identical panels |
| `per_sample_markers` | list | Per-sample marker vectors |

### Interpreting the Output

The presence matrix and common markers list inform the homogenization strategy:

| Observation | Interpretation | Action |
|-------------|----------------|--------|
| All markers common | Panels are identical | Homogenization may be unnecessary; verify column order |
| Most markers common, few unique | Minor panel evolution or site-specific additions | Homogenize to the common set or the full union (filling missing with NA) |
| Few markers common, many unique | Substantially different panels | Homogenize to common set only; consider whether the common set is sufficient for your analysis |
| Markers appear similar but not matching | Naming convention differences | Build a rename dictionary before re-running comparison |

---

## Panel Homogenization

### Purpose

Once panel differences are understood, `homogenize_flowframes()` resolves them by transforming every sample to a common target panel. After homogenization, all samples have identical marker columns in identical order, enabling valid concatenation and joint analysis.

### The Three-Step Process

#### Step 1: Rename Markers (Synonym Resolution)

Different naming conventions for the same biological marker are resolved using a user-provided dictionary (`markername_replacing`). The `replace_markernames()` function iterates through the flowFrame's marker annotations and substitutes any marker name found in the dictionary with its canonical replacement.

```r
replace_markernames <- function(flowframe, replace_map) {
  marker_list <- flowCore::markernames(flowframe)
  marker_pos <- explore_data(flowframe)$mcols
  for (i in seq_along(marker_list)) {
    if (marker_list[i] %in% names(replace_map)) {
      new_marker <- replace_map[[marker_list[i]]]
      names(new_marker) <- flowCore::featureNames(flowframe)[marker_pos[i]]
      flowCore::markernames(flowframe)[i] <- new_marker
    }
  }
  flowframe
}
```

This step addresses the naming heterogeneity problem. Common renaming scenarios include:

| Original Name | Canonical Name | Reason |
|--------------|----------------|--------|
| `CD3-FITC` | `CD3` | Fluorochrome suffix included in marker name |
| `CD3_FITC` | `CD3` | Underscore separator variant |
| `TCRab` | `TCRa/b` | Alternative abbreviation |
| `HLA-DR` | `HLADR` | Hyphenation difference |
| `CD279 (PD-1)` | `PD-1` | Parenthetical alternate name |
| `<PE-Cy7-A> CD8` | `CD8` | Channel prefix in marker field |

#### Step 2: Reorder to Target Panel

The `reorder_exprmtx()` function creates a new expression matrix with columns in the exact order specified by the target marker list (`markers_of_interest`). For each target marker, it searches the flowFrame's marker annotations for a match. If found, the corresponding expression column is copied into the correct position. If not found, the column is filled with `NA`.

```r
reorder_exprmtx <- function(flowframe, marker_name_list) {
  data <- explore_data(flowframe)
  result <- matrix(NA, nrow = nrow(flowCore::exprs(flowframe)),
                   ncol = length(marker_name_list))
  colnames(result) <- marker_name_list
  for (j in seq_along(marker_name_list)) {
    col_idx <- which(data$markers == marker_name_list[j])
    if (length(col_idx) > 0) {
      result[, j] <- flowCore::exprs(flowframe)[, col_idx[1]]
    }
  }
  result
}
```

This step simultaneously accomplishes two things:

1. **Reordering**: Markers that are present but in a different column position are moved to the correct position.
2. **NA-filling**: Markers that are absent from this sample are represented as columns of `NA` values, maintaining the matrix dimensions.

#### Step 3: Convert Back to flowFrame

The reordered expression matrix is converted back to a flowFrame object using `from_mtx_to_ffr()` (from `src/R_shared/flowframe_utils.R`), which constructs the necessary `AnnotatedDataFrame` metadata (parameter names, descriptions, ranges) from the matrix. The resulting flowFrame has the target panel's markers as both feature names and marker descriptions.

### Handling Missing Markers

When a sample lacks a marker that is in the target panel, `reorder_exprmtx()` fills the corresponding column with `NA`. This is a deliberate design choice with specific implications:

- **Downstream clustering**: Algorithms such as [Leiden clustering](leiden-clustering.md) that compute cell-to-cell distances will need to handle `NA` values. In practice, you should either (a) restrict the target panel to `common_markers` so no NAs are introduced, or (b) use imputation or NA-tolerant distance metrics.
- **Visualization**: [UMAP](umap.md) and [heatmaps](../subprojects/visualization-ref.md) can handle NA columns by excluding them from computation, but the analyst must be aware of which markers are imputed versus measured.
- **Statistical analysis**: Differential abundance analysis via [NB-GLMM](negative-binomial-glmm.md) or [mixed-effects logistic regression](mixed-effects-logistic.md) operates on cluster proportions, not raw marker values, so NA markers in the expression matrix do not directly affect abundance statistics -- provided the clustering step handled them correctly.

The safest approach for most analyses is to set `markers_of_interest` to the `common_markers` returned by `compare_panels_detailed()`, ensuring that every sample has measured values for every marker in the homogenized panel.

### Output Options

`homogenize_flowframes()` can optionally save the homogenized data as CSV files and/or FCS files:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `save_csv` | `FALSE` | Write reordered expression matrix as CSV per sample |
| `save_fcs` | `FALSE` | Write reordered flowFrame as FCS per sample |
| `output_dir` | `NULL` | Directory for saved files (required if saving) |

Output filenames follow the pattern `hg_<sample_name>.csv` and `hg_<sample_name>.fcs`.

---

## When Your Data Calls For It

Panel homogenization is necessary whenever samples that will be analyzed together have non-identical marker panels. Specific scenarios include:

### Multi-Site Studies with Different Panels

Consortium studies and multi-center clinical trials frequently use [site-specific panel configurations that share a core set of markers but differ in supplementary markers](https://doi.org/10.1038/nri3158). For example, the [HIPC consortium](https://doi.org/10.1038/ni.2549) established standardized panels for immune monitoring, yet participating sites often supplemented the core panel with site-specific markers of interest. Homogenization to the shared core enables cross-site analysis while preserving the option to analyze site-specific markers separately.

### Panels That Evolved Over Time

Longitudinal studies spanning months or years are particularly susceptible to panel evolution. [Reagent discontinuation, new antibody clone availability, and evolving biological questions](https://doi.org/10.1002/cyto.a.22624) all motivate panel changes. A study that begins with a 15-marker panel and adds 3 markers at the 6-month timepoint creates two panel versions that must be harmonized before joint analysis.

### Different Naming Conventions

Even when the biological panel is identical, [inconsistent annotation of FCS file parameters](https://doi.org/10.1002/cyto.a.24225) across instruments, software versions, or operators can prevent automated matching. The marker renaming step in `homogenize_flowframes()` resolves this by mapping all naming variants to a canonical form.

Common sources of naming inconsistency:

| Source | Example | Resolution |
|--------|---------|------------|
| Fluorochrome included in marker name | `CD3-FITC` vs. `CD3` | Strip fluorochrome suffix |
| Separator variation | `CD3_FITC` vs. `CD3-FITC` | Normalize separators |
| Case differences | `cd3` vs. `CD3` | Standardize to uppercase |
| Channel name leakage | `<FITC-A> CD3` vs. `CD3` | Strip channel prefix |
| Alternative nomenclature | `PD-1` vs. `CD279` | Map to preferred name |

### Before Any Cross-Sample Analysis

Any analysis that combines data from multiple samples -- [pooled clustering, joint UMAP embedding, cross-sample differential abundance testing](https://doi.org/10.12688/f1000research.11622.3) -- requires that all samples present the same features. Panel homogenization is the step that establishes this requirement. Without it, the combined expression matrix will have inconsistent column semantics.

---

## When to Skip

Panel homogenization is unnecessary and should be skipped when:

### Single-Panel Study

All samples were stained with the same antibody panel, acquired on the same instrument with the same settings, and exported with the same software. In this case, `compare_panels_detailed()` will report that all markers are common and all samples belong to a single panel group. Homogenization would be a no-op.

### All Samples Already Have Identical Markers

Even in multi-site or longitudinal studies, it is possible that all sites used the same panel and naming conventions. Run `compare_panels_detailed()` as a verification step: if `not_common` is empty and all samples fall into a single panel group, skip homogenization.

### After FlowJo Gating Export

If your workflow uses FlowJo-exported gating results (Pipeline D), the analysis operates on gate labels and population frequencies rather than raw expression matrices. Panel differences at the expression level are irrelevant because the gating software has already resolved marker identities internally.

---

## How to Use in the Pipeline

### Step 1: Compare Panels

Always begin by diagnosing panel heterogeneity:

```r
source("src/harmonization/panel_comparison.R")

panel_info <- compare_panels_detailed(flow_frames_list)

# Review the diagnostic output
cat("Common markers:", length(panel_info$common_markers), "\n")
cat("Total unique markers:", length(panel_info$all_markers), "\n")
cat("Non-common markers:", length(panel_info$not_common), "\n")
cat("Number of distinct panels:", length(panel_info$panel_groups), "\n")

# Inspect the presence matrix
print(panel_info$presence_matrix)
```

### Step 2: Build a Rename Dictionary (if Needed)

If the comparison reveals naming inconsistencies (e.g., `CD3-FITC` and `CD3` both appear as distinct markers), construct a rename dictionary mapping non-canonical names to their canonical forms:

```r
rename_dict <- list(
  "CD3-FITC"  = "CD3",
  "CD4-PE"    = "CD4",
  "CD8-APC"   = "CD8",
  "CD19-BV421" = "CD19",
  "CD56-PE-Cy7" = "CD56"
)
```

### Step 3: Choose the Target Panel

Decide which markers to include in the homogenized panel:

```r
# Option A: Use only common markers (safest, no NAs)
target_markers <- panel_info$common_markers

# Option B: Use all markers (fills missing with NA)
target_markers <- panel_info$all_markers

# Option C: Use a curated subset
target_markers <- c("CD3", "CD4", "CD8", "CD19", "CD14", "CD56",
                     "CD45", "CD25", "CD127", "HLADR")
```

### Step 4: Homogenize

```r
source("src/harmonization/panel_homogenization.R")

homogenized <- homogenize_flowframes(
  flow_frames          = flow_frames_list,
  sample_names         = sample_names,
  markers_of_interest  = target_markers,
  markername_replacing = rename_dict,
  save_fcs             = TRUE,
  save_csv             = TRUE,
  output_dir           = "output/02_homogenized/"
)

# homogenized is a named list of flowFrame objects
# All have identical markers in identical order
```

### Step 5: Verify

Re-run panel comparison on the homogenized output to confirm that all samples now have identical panels:

```r
verify <- compare_panels_detailed(homogenized)
stopifnot(length(verify$not_common) == 0)
cat("Homogenization verified: all samples have identical panels.\n")
```

---

## Homogenization in the Pipeline

Panel homogenization sits at a specific point in the processing sequence. It must occur after data loading and before any step that assumes a shared feature space.

### Pipeline A (Raw FCS)

```
Load FCS -> QC -> Transform -> [Homogenize panels] -> Batch correct -> Downsample -> Cluster -> Analyze -> Visualize
```

### Pipeline E (Preprocessed FCS)

```
Load preprocessed FCS -> [Homogenize panels] -> Batch correct -> Cluster -> Visualize
```

### Relationship to Adjacent Steps

| Step | Relationship | Link |
|------|-------------|------|
| **Loading** | Produces the list of flowFrame objects that `compare_panels_detailed()` and `homogenize_flowframes()` consume | `src/loading/` |
| **Batch Correction** | Requires homogenized input: all samples must have the same markers before [Harmony can compute inter-batch corrections](https://doi.org/10.1038/s41592-019-0619-0) | [Batch Correction](batch-correction.md) |
| **Downsampling** | Operates on the homogenized flowFrames to reduce cell counts for computational efficiency | `src/preprocessing/downsampling.R` |
| **Leiden Clustering** | Builds a KNN graph from the combined expression matrix, which requires [consistent marker columns across all samples](https://doi.org/10.1038/s41598-019-41695-z) | [Leiden Clustering](leiden-clustering.md) |
| **UMAP** | Computes a [joint embedding across all samples](https://doi.org/10.48550/arXiv.1802.03426), which requires consistent feature dimensions | [UMAP](umap.md) |

### Why Homogenization Must Come Before Batch Correction

[Batch correction algorithms such as Harmony](https://doi.org/10.1038/s41592-019-0619-0) operate on a concatenated matrix where every row has the same columns. If samples have different marker sets, concatenation either fails (different column counts) or silently misaligns columns (same count but different markers). Homogenization ensures that the matrix passed to batch correction has consistent semantics in every column.

---

## Best Practices

### Creating Marker Dictionaries

A well-constructed rename dictionary is the most important manual input to the homogenization process. Guidelines:

1. **Start from the data**: Use the `per_sample_markers` output of `compare_panels_detailed()` to identify all unique marker names across your dataset.
2. **Choose canonical names**: Select a single canonical name for each biological marker. Prefer short, unambiguous names without fluorochrome suffixes (e.g., `CD3` rather than `CD3-FITC`). The [OMIP series of publications](https://doi.org/10.1002/cyto.a.22916) provides standardized marker nomenclature that can serve as a reference.
3. **Document the dictionary**: Save the rename dictionary as a separate R file or CSV in `config/study_params/` so it can be reviewed, version-controlled, and reused.
4. **Handle ambiguity explicitly**: If two different biological markers share the same fluorochrome channel (e.g., two different panels use PE-Cy7 for different markers), the FCS `$PnS` annotation should distinguish them. If it does not, manual inspection is required.

Example of a documented dictionary:

```r
# config/study_params/marker_rename_dict.R
# Rename dictionary for Study XYZ
# Maps site-specific and instrument-specific marker names to canonical forms
marker_rename_dict <- list(
  # Site A uses fluorochrome suffixes
  "CD3-FITC"     = "CD3",
  "CD4-PE"       = "CD4",
  "CD8-APC"      = "CD8",
  "CD19-BV421"   = "CD19",

  # Site B uses underscores
  "CD3_FITC"     = "CD3",
  "CD4_PE"       = "CD4",

  # Alternative nomenclature
  "TCRab"        = "TCRa/b",
  "CD279"        = "PD-1",
  "CD274"        = "PD-L1"
)
```

### Choosing the Target Panel

The choice of target panel determines what biological questions can be addressed:

| Strategy | Target | Pros | Cons |
|----------|--------|------|------|
| **Common markers only** | `panel_info$common_markers` | No NA values; all data is measured | May lose informative markers present in most but not all samples |
| **Full union** | `panel_info$all_markers` | Retains all available information | Introduces NA columns; requires NA-tolerant downstream methods |
| **Curated subset** | Manually selected markers | Focused on biological question; avoids irrelevant markers | Requires domain knowledge to select |

For most immunophenotyping studies, restricting to common markers is the safest choice. If common markers exclude a marker critical for distinguishing a key population (e.g., CD56 for NK cells), consider adding it to the target panel and handling the resulting NAs explicitly.

### Handling Missing Markers

When the target panel includes markers not present in all samples:

1. **Document which samples lack which markers**: The presence matrix from `compare_panels_detailed()` provides this information.
2. **Exclude NA markers from clustering**: When running [Leiden clustering](leiden-clustering.md) or [UMAP](umap.md), use only the subset of target markers that are present in all samples. The full target panel can still be used for per-sample analyses where NAs are handled individually.
3. **Consider separate analyses**: If two groups of samples have substantially different panels (e.g., fewer than 50% overlap), it may be more appropriate to analyze them separately and compare results at the population level rather than forcing homogenization.

### OMIP Panels as Reference Standards

The [Optimized Multicolor Immunofluorescence Panel (OMIP) series](https://doi.org/10.1002/cyto.a.22916) published in *Cytometry Part A* provides peer-reviewed, standardized panel designs for specific biological questions. OMIP panels are valuable references when designing rename dictionaries and choosing target markers because they establish [community-consensus marker combinations for defined immunophenotyping goals](https://doi.org/10.1002/cyto.a.22916).

Relevant OMIP publications for common panel design scenarios:

| OMIP | Focus | Markers | Citation |
|------|-------|---------|----------|
| OMIP-069 | [40-color full-spectrum panel for deep human immunophenotyping](https://doi.org/10.1002/cyto.a.24213) | 40 surface markers covering all major immune lineages | [Liechti et al. (2021)](https://doi.org/10.1002/cyto.a.24213) |
| OMIP-009 | [Characterization of antigen-responsive human T cells](https://doi.org/10.1002/cyto.a.22269) | T cell subset markers including activation and memory | [Chattopadhyay et al. (2012)](https://doi.org/10.1002/cyto.a.22269) |
| OMIP-030 | [Assessment of human T cell and NK cell immunodeficiency](https://doi.org/10.1002/cyto.a.22792) | T and NK cell lineage and maturation markers | [Hernandez et al. (2016)](https://doi.org/10.1002/cyto.a.22792) |
| OMIP-044 | [28-color immunophenotyping of the human dendritic cell compartment](https://doi.org/10.1002/cyto.a.23331) | Dendritic cell subsets and activation markers | [Jahn et al. (2018)](https://doi.org/10.1002/cyto.a.23331) |

When your study's panel overlaps with a published OMIP, adopting the OMIP's marker nomenclature in your rename dictionary improves reproducibility and facilitates comparison with other studies using the same standard.

### Quality Control After Homogenization

After homogenization, verify the result before proceeding:

1. **Re-run `compare_panels_detailed()`** on the homogenized output. All markers should be common and all samples should belong to a single panel group.
2. **Check expression ranges**: For each marker, compare the range (min, max, median) across samples. Large discrepancies may indicate incorrect renaming (e.g., mapping two biologically different markers to the same canonical name) or pre-existing [batch effects that require correction](batch-correction.md).
3. **Inspect NA patterns**: If your target panel includes markers not present in all samples, verify that the NA columns correspond to the expected missing markers using the presence matrix.

---

## Implementation Details

### Source Code Locations

- **Panel comparison:** [`compare_panels_detailed()`](/src/harmonization/panel_comparison.R) in `src/harmonization/panel_comparison.R`
- **Panel homogenization:** [`homogenize_flowframes()`](/src/harmonization/panel_homogenization.R) in `src/harmonization/panel_homogenization.R`
- **Shared utilities:** [`explore_data()`](/src/R_shared/flowframe_utils.R), [`compare_panels()`](/src/R_shared/flowframe_utils.R), [`from_mtx_to_ffr()`](/src/R_shared/flowframe_utils.R) in `src/R_shared/flowframe_utils.R`

### Dependencies

- `flowCore` (for `exprs()`, `markernames()`, `featureNames()`, `write.FCS()`)
- `src/R_shared/io_utils.R` (for `make_outdir()`)
- `src/R_shared/flowframe_utils.R` (for `explore_data()`, `compare_panels()`, `from_mtx_to_ffr()`)

### How `homogenize_flowframes()` Works Internally

1. **Creates output directory** if `output_dir` is provided, using `make_outdir()`.
2. **Generates cell IDs** across all input flowFrames for tracking provenance.
3. **Iterates** over each flowFrame:
   a. Applies marker renaming via `replace_markernames()` if `markername_replacing` is non-empty.
   b. Reorders the expression matrix to match `markers_of_interest` via `reorder_exprmtx()`, filling missing markers with `NA`.
   c. Converts the reordered matrix back to a flowFrame via `from_mtx_to_ffr()`.
   d. Optionally writes CSV and/or FCS output.
4. **Returns** a named list of homogenized flowFrame objects.

---

## Related Documentation

- [Batch Correction](batch-correction.md) -- Remove technical variation between acquisition batches after homogenization
- [Leiden Clustering](leiden-clustering.md) -- Unsupervised population identification requiring homogenized input
- [UMAP](umap.md) -- Joint embedding requiring consistent marker dimensions
- [PCA](pca.md) -- Dimensionality reduction on the homogenized expression matrix
- [Arcsinh Transform](arcsinh-transform.md) -- Transformation applied before homogenization (CyTOF)
- [Logicle Transform](logicle-transform.md) -- Transformation applied before homogenization (fluorescence)
