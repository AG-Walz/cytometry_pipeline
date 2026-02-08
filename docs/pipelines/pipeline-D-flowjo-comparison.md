[< Back to Documentation Index](../README.md)

# Pipeline D: FlowJo Workspace Comparison

## Purpose

Compare gating strategies from multiple FlowJo workspaces. This pipeline is designed for inter-rater or inter-strategy comparisons: given two or more `.wsp` files that gate the same FCS data, it extracts gate labels, computes overlap statistics (Fleiss' kappa), and generates Venn/Euler and UpSet plots showing agreement and disagreement between gating strategies.

## Pipeline Flow

```
Load WSP -> Compare gatings -> Overlap statistics -> Venn / UpSet
```

```
  .wsp files + matching .fcs files
      |
      v
  [1] Loop: load_flowjo_wsp()          # Per-workspace gate extraction
        + get_population_strings_list()
      |
      v
  [2] compare_gatings()                # Cross-workspace comparison
      |
      v
  [3] compute_fleiss_kappa()           # Inter-rater agreement
      + Venn/Euler + UpSet plots
```

## Input

- Multiple `.wsp` files (FlowJo workspace files), each representing a different gating strategy or expert.
- Matching `.fcs` files referenced by the workspaces.

All workspaces must reference the same FCS file(s) for comparison to be meaningful.

## Configuration

```bash
Rscript pipelines/pipeline_D_flowjo_comparison.R --config config/study_params/Compare_gatings_my_study.R
```

**The `--config` flag is mandatory for this pipeline.** The pipeline calls `stop()` if no config file is provided, because workspace paths and label mappings are inherently study-specific.

Example config files are provided as templates:

- `config/study_params/Compare_gatings_Endo_CD45_Julia.R`
- `config/study_params/Compare_gatings_ICK_broad_narrow.R`

## Key Config Variables

| Variable | Type | Description |
|----------|------|-------------|
| `wsp_paths` | character vector | Paths to the `.wsp` files to compare |
| `custom_name` | character vector | Short nicknames for each workspace (used in plot legends) |
| `target_fcs` | character | FCS filename to compare across workspaces |
| `label_dict` | named character vector | Maps full FlowJo gate paths to simplified display labels |
| `idx_labels` | character or `NULL` | Optional path to a CSV file with pre-computed labels (alternative to WSP-derived labels) |
| `output_dir` | character | Root directory for output subdirectories |

### Example `label_dict`

The `label_dict` maps verbose FlowJo hierarchy paths to concise labels suitable for plots and tables:

```r
label_dict <- c(
  "/no erys/Single Cells/lymphocytes/CD3+"                              = "CD3+",
  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+"                 = "CD4+ T cells",
  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+"                 = "CD8+ T cells",
  "/no erys/Single Cells/lymphocytes/CD3-/CD19+"                        = "B cells",
  "/no erys/Single Cells/lymphocytes/CD3-/CD56+"                        = "NK cells"
)
```

Only gates present in `label_dict` are included in the comparison. Gates not listed are silently skipped.

## Steps in Detail

### Step 1: Load FlowJo Workspaces

```r
# src/loading/load_flowjo_wsp.R
all_labels <- list()
for (i in seq_along(params$wsp_paths)) {
  ws <- load_flowjo_wsp(params$wsp_paths[i])
  pops <- get_population_strings_list(ws, target_fcs = params$target_fcs)
  all_labels[[params$custom_name[i]]] <- pops
}
```

Each workspace is opened with `flowWorkspace::open_flowjo_xml()`, and gate memberships are extracted for the target FCS file. The `get_population_strings_list()` function returns a per-event data frame of gate paths, which are then mapped through `label_dict` to produce simplified labels.

### Step 2: Compare Gatings

```r
# src/metadata_qc/compare_gatings.R
comparison <- compare_gatings(
  all_labels,
  label_dict  = params$label_dict,
  custom_name = params$custom_name
)
```

Aligns events across workspaces (by event index) and produces a comparison matrix: for each event and each population, a binary vector indicates which workspaces assigned that event to that population. Generates Venn/Euler diagrams per population and an UpSet plot showing all intersection patterns.

### Step 3: Overlap Statistics

```r
# src/metadata_qc/overlap_statistics.R
kappa <- compute_fleiss_kappa(
  comparison,
  output_dir = file.path(params$output_dir, "02_overlap_stats")
)
```

Computes Fleiss' kappa for each population across all raters (workspaces). Kappa values are written to `fleiss_kappa_log.txt` with interpretation guidelines:

| Kappa | Interpretation |
|-------|----------------|
| < 0.20 | Poor agreement |
| 0.21 -- 0.40 | Fair agreement |
| 0.41 -- 0.60 | Moderate agreement |
| 0.61 -- 0.80 | Substantial agreement |
| 0.81 -- 1.00 | Almost perfect agreement |

## Output Directory Structure

```
output_dir/
  01_gating_comparison/    # Venn/Euler plots (per population), UpSet plot, label CSVs
  02_overlap_stats/        # fleiss_kappa_log.txt, per-population kappa values
```

## Notes

- All workspaces must gate the same FCS file (specified by `target_fcs`). If a workspace does not contain the target file, the pipeline stops with an error.
- The `idx_labels` parameter provides an alternative input path: if a CSV with pre-computed per-event labels is available (e.g., from a previous Pipeline C run), it can be included as an additional "rater" in the comparison.
- Venn/Euler diagrams are limited to 2--5 sets; for more workspaces, only the UpSet plot is generated.
- This pipeline does not perform any expression analysis or clustering -- it is purely a gating agreement tool.
