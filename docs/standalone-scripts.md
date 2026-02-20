[< Back to Documentation Index](README.md)

# Standalone Scripts Reference

Every script in `src/` can be used independently outside of a pipeline. Each script defines one or more functions that you call after sourcing.

---

## Loading Scripts (`src/loading/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `load_fcs_raw.R` | `load_fcs_raw()` | `(input_dir, sample_names=NULL, truncation=FALSE, min_limit=NULL, extract_csv=FALSE, output_dir=NULL)` | `list(flow_frames, sample_names, input_files)` |
| `load_csv_expression.R` | `load_csv_expression()` | `(input_dir, sep=NULL, max_memory_gb=4)` | Named list of data.tables |
| `load_csv_gating.R` | `load_csv_gating()` | `(input_fcs_path, input_csv_path, output_path, gating_pattern, csv_label_dict=NULL, gating_tree=NULL, processing_mode="automatic", output_format="wide", chunk_size=100000)` | `list(output_dir, results, n_processed, n_total)` |
| `load_flowjo_wsp.R` | `load_flowjo_wsp()` | `(wsp_path, fcs_dir, ...)` | GatingSet objects |
| `load_manual_gating.R` | `load_manual_gating()` | `(csv_dir, fcs_dir, ...)` | Labelled data.tables |
| `create_flowframes.R` | `create_flowframes_from_matrices()` | `(expr_list, ...)` | List of flowFrames |

**Example — Loading FCS files standalone:**

```r
setwd("/path/to/cytometry/src/loading")
source("load_fcs_raw.R")

result <- load_fcs_raw(
  input_dir  = "/path/to/fcs/files",
  truncation = FALSE
)
# result$flow_frames  -> list of flowFrame objects
# result$sample_names -> c("sample_001", "sample_002", ...)
```

---

## Preprocessing Scripts (`src/preprocessing/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `preprocess_pipeline.R` | `run_preprocessing()` | `(flow_frames, params, output_dir=NULL)` | Named list of preprocessed flowFrames |
| `quality_control.R` | `quality_control()` | `(flowframe, sample_name, qc_type, marker_indices, output_dir)` | Filtered flowFrame |
| `transformation.R` | `transform_flowframe()` | `(flowframe, sample_name, markers_of_interest, offset_list, width_list, transform_type)` | Transformed flowFrame |
| `downsampling.R` | `downsample_flowframes()` | `(flow_frames, sample_names, target_cells=10000, method="Maximal", seed=42, save_fcs=FALSE, output_dir=NULL)` | Named list of downsampled flowFrames |
| `normalization.R` | `z_score_scaling()` / `quantile_norm()` | `(flowframe, marker_indices)` | Normalized flowFrame |

**Example — Downsampling standalone:**

```r
setwd("/path/to/cytometry/src/preprocessing")
source("downsampling.R")

downsampled <- downsample_flowframes(
  flow_frames  = my_flow_frames,
  sample_names = c("S1", "S2", "S3"),
  target_cells = 5000,
  method       = "Consistent",
  seed         = 42
)
```

---

## Harmonization Scripts (`src/harmonization/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `panel_homogenization.R` | `homogenize_flowframes()` | `(flow_frames, sample_names, markers_of_interest, markername_replacing=list(), save_csv=FALSE, save_fcs=FALSE, output_dir=NULL)` | Named list of homogenized flowFrames |
| `panel_comparison.R` | `compare_panels_detailed()` | `(flow_frames)` | Panel comparison report |
| `batch_correction.R` | `run_batch_correction()` | `(flow_frames, batch_labels, ...)` | Batch-corrected flowFrames |

**Example — Homogenizing panels standalone:**

```r
setwd("/path/to/cytometry/src/harmonization")
source("panel_homogenization.R")

homogenized <- homogenize_flowframes(
  flow_frames         = my_flow_frames,
  sample_names        = c("S1", "S2", "S3"),
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19"),
  markername_replacing = list(`PD1 (CD279)` = "PD1")
)
```

---

## Clustering Scripts (`src/clustering/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `leiden.R` | `run_leiden()` | `(expr_matrix, n_neighbors=30, resolution=1.0, seed=42)` | Integer vector of cluster IDs |
| `leiden.R` | `run_leiden_with_umap()` | `(expr_matrix, sample_labels=NULL, output_dir=NULL, n_neighbors=30, resolution=1.0, seed=42)` | data.frame(UMAP1, UMAP2, cluster, sample) |

**Example — Leiden clustering standalone:**

```r
setwd("/path/to/cytometry/src/clustering")
source("leiden.R")

result <- run_leiden_with_umap(
  expr_matrix   = combined_expression_matrix,
  sample_labels = sample_labels_vector,
  output_dir    = "/path/to/output",
  n_neighbors   = 15,
  resolution    = 0.8,
  seed          = 42
)
```

---

## Analysis Scripts (`src/analysis/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `cluster_splitting.R` | `split_by_cluster()` | `(input_dir, output_dir, cluster_regex="Clustering annotated", markers_of_interest=NULL, min_cells_per_cluster=10)` | `list(tracking_matrix, cluster_summary, output_files, markers_used)` |
| `differential_abundance.R` | `run_differential_abundance()` | `(data_dir, conditions_csv, output_dir, cluster_regex, condition_levels, markers_of_interest=NULL, min_cluster_size=10, p_adjust_method="BH", seed=123)` | data.frame of GLMM results |
| `cluster_enrichment.R` | `run_cluster_enrichment()` | `(data_dir, conditions_csv, output_dir, cluster_regex, condition_levels, stat_method="mixed_effects", min_cluster_frac=0.001, p_adjust_method="BH", seed=123)` | data.frame of enrichment results |
| `marker_statistics.R` | `run_marker_statistics()` | `(input_folder, output_folder, filename_regex="^(.+)_(.+)\\.csv$", markers=NULL, test_type="lmer", n_workers=NULL)` | data.frame of per-marker stats |
| `abundance_plotting.R` | `build_tracking_matrix()` | `(input_dir, cluster_regex, exclude_patterns, output_dir)` | data.frame (tracking matrix) |
| `abundance_plotting.R` | `plot_abundance()` | `(tracking_matrix_path=NULL, output_dir, tracking_matrix=NULL, metadata_file=NULL)` | invisible(NULL) |

**Example — Cluster enrichment standalone:**

```r
setwd("/path/to/cytometry/src/analysis")
source("cluster_enrichment.R")

results <- run_cluster_enrichment(
  data_dir         = "/path/to/csvs/with/gating/columns",
  conditions_csv   = "/path/to/conditions.csv",
  output_dir       = "/path/to/output",
  condition_levels = c("Control", "Treatment"),
  stat_method      = "pooled"
)
```

---

## Visualization Scripts (`src/visualization/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `umap.R` | `plot_umap()` | `(flow_frames, sample_names, markers_of_interest=NULL, output_dir, n_neighbors=15, min_dist=0.1, output_formats=list(pdf=F, png=T, svg=T), seed=42)` | list(umap_coords, sample_labels, plots) |
| `heatmaps.R` | `plot_heatmap()` | `(expr_matrix, output_path, title, scale="row", width=10, height=8)` | invisible(NULL) |
| `heatmaps.R` | `plot_cluster_heatmap()` | `(expr_data, cluster_col, marker_cols, output_path, title)` | invisible(NULL) |
| `density_plots.R` | `plot_density()` | `(flow_frames, sample_names, markers, output_dir, plot_type="density", group_by="sample")` | invisible(NULL) |
| `barplots.R` | `plot_composition_barplot()` | `(data, output_dir, x_col="cluster", fill_col="sample", position="stack", title)` | invisible(ggplot) |

---

## Metadata QC Scripts (`src/metadata_qc/`)

| Script | Function | Signature |
|--------|----------|-----------|
| `check_metadata.R` | `check_metadata()` | `(metadata_file, ...)` |
| `check_acquisition.R` | `check_acquisition_dates()` | `(parent_folder, ...)` |
| `overlap_statistics.R` | `compute_fleiss_kappa()` | `(label_lists, ...)` |
| `compare_gatings.R` | `compare_gatings()` | `(pop_strings_list, ...)` |
