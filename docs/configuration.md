[< Back to Documentation Index](README.md)

# Configuration Guide

All pipelines are driven by a `params` list defined in a plain R configuration file. Each pipeline has a template in `config/`.

---

## Creating a Config

Copy a pipeline template and edit it for your study:

```bash
cp config/pipeline_A_config.R config/study_params/my_study_A.R
# Edit paths and parameters, then run:
Rscript pipelines/pipeline_A_raw_fcs.R --config config/study_params/my_study_A.R
```

Config files automatically compute `base_dir` from their own file location, so all relative paths work regardless of where you run from.

---

## Key Parameters Across Pipelines

| Parameter | Used By | Description |
|-----------|---------|-------------|
| `fcs_folder` | A, E | Directory containing `.fcs` files |
| `csv_folder` | B, C | Directory containing CSV files |
| `output_dir` | All | Root output directory (subdirs created automatically) |
| `markers_of_interest` | A, E | Character vector of marker names to analyze |
| `markers` | B | Character vector of marker column names |
| `conditions_csv` | A, C | Path to file_name -> condition mapping CSV |
| `condition_levels` | A, C | Two-element vector: `c("Control", "Treatment")` |
| `truncation` | A, E | Truncate FCS values to `$PnR` max range |
| `qc_type` | A | Quality control: `"PeacoQC"`, `"flowAI"`, `"None"` |
| `transform_type` | A | Transformation: `"arcsinh"`, `"logicle"`, `"biexponential"`, `"None"` |
| `pp_order` | A | Preprocessing order: `"QT"` (QC first) or `"TQ"` (Transform first) |
| `offset_list` / `width_list` | A | Per-marker numeric vectors for transformation (length must match markers) |
| `downsample_method` | A, E | `"Consistent"` (fixed N per sample) or `"Maximal"` (proportional) |
| `downsample_n` | A, E | Target cell count for downsampling |
| `n_neighbors` | A, E | KNN graph neighbors for UMAP/Leiden |
| `resolution` | A, E | Leiden clustering resolution (higher = more clusters) |
| `cluster_regex` | A, C | Regex identifying gating/cluster columns |
| `analysis_method` | C | `"enrichment"` (Fisher) or `"glmm"` (NB-GLMM) |
| `stat_method` | C | `"pooled"` or `"mixed_effects"` |
| `seed` | All | Random seed for reproducibility |
| `output_formats` | All | `list(pdf=TRUE, png=TRUE, svg=FALSE)` |

---

## Standalone Script Configs

For running individual scripts outside pipelines, see `config/standalone/` for ready-to-edit config templates.

---

## Example: Pipeline A Config

```r
# config/pipeline_A_config.R (key parameters)
params <- list(
  fcs_folder          = "/path/to/fcs/files",
  output_dir          = "/path/to/output",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),
  conditions_csv      = "/path/to/conditions.csv",
  condition_levels    = c("Control", "Treatment"),
  truncation          = FALSE,
  pp_order            = "QT",
  qc_type             = "None",
  transform_type      = "arcsinh",
  offset_list         = rep(0, 6),
  width_list          = rep(5, 6),
  downsample_method   = "Consistent",
  downsample_n        = 10000,
  n_neighbors         = 15,
  resolution          = 0.5,
  seed                = 42
)
```

See each [pipeline guide](pipelines/) for pipeline-specific config examples and parameter details.
