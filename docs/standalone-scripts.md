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
| `load_flowjo_wsp.R` | `load_multiple_workspaces()` | `(wsp_paths, names, fcs_dirs=NULL)` | Named list of population_strings_lists |
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
| `compensation.R` | `extract_compensation_matrix()` | `(flowframe)` | Compensation matrix or NULL |
| `compensation.R` | `apply_compensation()` | `(flowframe, comp_matrix=NULL, auto_extract=TRUE)` | Compensated flowFrame |
| `compensation.R` | `load_compensation_matrix()` | `(csv_path, sep=";")` | Compensation matrix |
| `spectral_unmixing.R` | `build_unmixing_matrix()` | `(control_frames, unstained_frame, channels=NULL)` | Unmixing matrix (channels x fluorochromes) |
| `spectral_unmixing.R` | `unmix_spectral()` | `(flowframe, unmixing_matrix, autofluorescence=NULL)` | Unmixed flowFrame with virtual marker channels |
| `bead_normalization.R` | `normalize_beads()` | `(flowframe, bead_channels=NULL, bead_type="four_element", remove_beads=TRUE, baseline="median", window_size=500, output_dir=NULL)` | `list(normalized_ff, bead_indices, normalization_factors)` |
| `bead_normalization.R` | `identify_beads()` | `(bead_expr, method="mahalanobis")` | Logical vector (TRUE = bead event) |
| `doublet_detection.R` | `detect_doublets()` | `(flowframe, area_col="FSC-A", height_col="FSC-H", method="ratio", threshold=1.1, width_col=NULL)` | `list(filtered_ff, doublet_indices, singlet_rate)` |
| `live_dead_gating.R` | `gate_live_cells()` | `(flowframe, viability_channel, method="threshold", threshold=NULL, dead_high=TRUE)` | `list(live_ff, dead_indices, viability_rate)` |
| `automated_gating.R` | `run_auto_gating()` | `(flowframe, method="flowDensity", gating_template=NULL, channels, population=NULL)` | `list(gated_ff, gate_coords, population_stats)` |
| `automated_gating.R` | `define_gate()` | `(gate_type, channels, coords)` | flowCore gate object (or list of 4 gates for "quadrant") |
| `gate_utils.R` | `apply_gate_sequence()` | `(flowframe, gate_list, logical_op="and")` | `list(gated_ff, gate_indices, population_stats)` |
| `gate_utils.R` | `boolean_gate()` | `(flowframe, gate_results, expression)` | Logical vector of included cells |
| `supercell.R` | `run_supercell()` | `(expr_matrix, sample_labels, markers, gam=20, k_knn=30, n_pc=10, do_pca=TRUE, seed=42)` | `list(supercell_matrix, supercell_labels, cell_map, supercell_sizes)` |
| `supercell.R` | `run_supercell_per_cluster()` | `(expr_matrix, cluster_ids, sample_labels, markers, gam=20, seed=42)` | `list(supercell_matrix, supercell_cluster_ids, supercell_labels, cell_map)` |

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

**Example — Compensation standalone:**

```r
setwd("/path/to/cytometry/src/preprocessing")
source("compensation.R")

# Auto-extract compensation matrix from FCS keywords and apply
compensated_ff <- apply_compensation(my_flowframe)

# Or load from external CSV
comp_matrix <- load_compensation_matrix("/path/to/spillover.csv")
compensated_ff <- apply_compensation(my_flowframe, comp_matrix = comp_matrix)
```

**Example — Doublet detection standalone:**

```r
setwd("/path/to/cytometry/src/preprocessing")
source("doublet_detection.R")

result <- detect_doublets(
  my_flowframe,
  method    = "ratio",
  threshold = 1.1
)
# result$filtered_ff   -> singlet flowFrame
# result$singlet_rate  -> 0.95 (95% singlets)
```

**Example — SuperCellCyto compression standalone:**

```r
setwd("/path/to/cytometry/src/preprocessing")
source("supercell.R")

sc_result <- run_supercell(
  expr_matrix   = combined_matrix,
  sample_labels = sample_labels,
  markers       = c("CD3", "CD4", "CD8", "CD19"),
  gam           = 20
)
# sc_result$supercell_matrix -> compressed expression matrix
```

---

## Harmonization Scripts (`src/harmonization/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `panel_homogenization.R` | `homogenize_flowframes()` | `(flow_frames, sample_names, markers_of_interest, markername_replacing=list(), save_csv=FALSE, save_fcs=FALSE, output_dir=NULL)` | Named list of homogenized flowFrames |
| `panel_comparison.R` | `compare_panels_detailed()` | `(flow_frames)` | Panel comparison report |
| `batch_correction.R` | `run_batch_correction()` | `(flow_frames, batch_labels, ...)` | Batch-corrected flowFrames |
| `cytonorm.R` | `run_cytonorm()` | `(train_frames, validation_frames, batch_labels, markers, n_clusters=10, goal="mean", seed=42)` | `list(corrected_frames, model)` |
| `combat.R` | `run_combat()` | `(expr_matrix, batch_vector, covariates=NULL, method="parametric")` | Corrected numeric matrix (same dimensions) |

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

**Example — CytoNorm batch correction standalone:**

```r
setwd("/path/to/cytometry/src/harmonization")
source("cytonorm.R")

result <- run_cytonorm(
  train_frames      = reference_flowframes,
  validation_frames = target_flowframes,
  batch_labels      = c("batch1", "batch1", "batch2", "batch2"),
  markers           = c("CD3", "CD4", "CD8"),
  n_clusters        = 10
)
# result$corrected_frames -> list of corrected flowFrames
```

**Example — ComBat batch correction standalone:**

```r
setwd("/path/to/cytometry/src/harmonization")
source("combat.R")

corrected_matrix <- run_combat(
  expr_matrix  = my_expression_matrix,
  batch_vector = c(rep("batch1", 500), rep("batch2", 500)),
  method       = "parametric"
)
```

---

## Clustering Scripts (`src/clustering/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `leiden.R` | `run_leiden()` | `(expr_matrix, n_neighbors=30, resolution=1.0, seed=42)` | Integer vector of cluster IDs |
| `leiden.R` | `run_leiden_with_umap()` | `(expr_matrix, sample_labels=NULL, output_dir=NULL, n_neighbors=30, resolution=1.0, seed=42)` | data.frame(UMAP1, UMAP2, cluster, sample) |
| `cluster_utils.R` | `clusters_to_gating_csv()` | `(combined_matrix, sample_labels, cluster_ids, output_dir, sep=";")` | Character vector of written file paths (invisibly) |
| `meta_clustering.R` | `run_resolution_sweep()` | `(expr_matrix, method="leiden", resolutions=seq(0.2,2.0,by=0.2), metric="silhouette", n_neighbors=30, seed=42)` | `list(best_clusters, best_param, scores_df, all_clusterings)` |
| `meta_clustering.R` | `consensus_clusters()` | `(cluster_list, method="majority", n_consensus_clusters=NULL)` | Integer vector of consensus cluster assignments |
| `flowjo_import.R` | `import_flowjo_gates()` | `(wsp_path, fcs_dir=dirname(wsp_path), target_populations=NULL)` | `list(gating_set, assignments, populations, sample_names)` |
| `cluster_assignment.R` | `get_cluster_assignments()` | `(source=c("leiden","flowjo","csv","flowsom"), ...)` | data.frame with cell_id and cluster columns |
| `cluster_assignment.R` | `get_clusters_from_leiden()` | `(csv_path)` | data.frame(cell_id, cluster, cluster_name) |
| `cluster_assignment.R` | `get_clusters_from_csv()` | `(csv_path, cluster_col="cluster")` | data.frame(cell_id, cluster, cluster_name) |
| `cluster_assignment.R` | `get_clusters_from_flowjo()` | `(wsp_path, fcs_dir=dirname(wsp_path))` | data.frame with cell_id and cluster columns |
| `cluster_assignment.R` | `get_clusters_from_flowsom()` | `(fsom_object)` | data.frame(cell_id, cluster, cluster_name) |

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

**Example — Resolution sweep standalone:**

```r
setwd("/path/to/cytometry/src/clustering")
source("meta_clustering.R")

sweep <- run_resolution_sweep(
  expr_matrix = combined_expression_matrix,
  resolutions = seq(0.2, 2.0, by = 0.2),
  metric      = "silhouette"
)
# sweep$best_clusters -> optimal cluster assignments
# sweep$best_param    -> best resolution value
# sweep$scores_df     -> data.frame of all scores
```

**Example — Unified cluster assignment standalone:**

```r
setwd("/path/to/cytometry/src/clustering")
source("cluster_assignment.R")

assignments <- get_cluster_assignments("csv", csv_path = "/path/to/clusters.csv")
# assignments$cell_id -> integer IDs
# assignments$cluster -> cluster labels
```

---

## Analysis Scripts (`src/analysis/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `cluster_splitting.R` | `split_by_cluster()` | `(input_dir, output_dir, cluster_regex="Clustering annotated", markers_of_interest=NULL, min_cells_per_cluster=10)` | `list(tracking_matrix, cluster_summary, output_files, markers_used)` |
| `cluster_splitting.R` | `tracking_matrix_to_long()` | `(tracking_matrix, id_col="file", names_to="cluster", values_to="count")` | data.frame with columns: file, cluster, count |
| `differential_abundance.R` | `run_differential_abundance()` | `(data_dir, conditions_csv, output_dir, cluster_regex, condition_levels, markers_of_interest=NULL, min_cluster_size=10, p_adjust_method="BH", seed=123)` | data.frame of GLMM results |
| `cluster_enrichment.R` | `run_cluster_enrichment()` | `(data_dir, conditions_csv, output_dir, cluster_regex, condition_levels, stat_method="mixed_effects", min_cluster_frac=0.001, p_adjust_method="BH", seed=123)` | data.frame of enrichment results |
| `marker_statistics.R` | `run_marker_statistics()` | `(input_folder, output_folder, filename_regex="^(.+)_(.+)\\.csv$", markers=NULL, test_type="lmer", n_workers=NULL)` | data.frame of per-marker stats |
| `abundance_plotting.R` | `build_tracking_matrix()` | `(input_dir, cluster_regex, exclude_patterns, output_dir)` | data.frame (tracking matrix) |
| `abundance_plotting.R` | `plot_abundance()` | `(tracking_matrix_path=NULL, output_dir, tracking_matrix=NULL, metadata_file=NULL)` | invisible(NULL) |
| `abundance_plotting.R` | `plot_abundance_stacked()` | `(tracking_matrix, output_dir)` | invisible(NULL) |
| `abundance_plotting.R` | `plot_abundance_proportional()` | `(tracking_matrix, output_dir)` | invisible(NULL) |
| `abundance_plotting.R` | `plot_abundance_totals()` | `(tracking_matrix, output_dir)` | invisible(NULL) |
| `abundance_plotting.R` | `plot_abundance_heatmap()` | `(tracking_matrix, output_dir)` | invisible(NULL) |
| `abundance_plotting.R` | `plot_abundance_by_metadata()` | `(tracking_matrix, metadata_file, output_dir, ...)` | invisible(NULL) |
| `expression_matrices.R` | `create_expression_matrices()` | `(input_dir, output_dir, markers_of_interest=NULL, cluster_regex, conditions_csv=NULL, ...)` | list(expression_matrices, summary) |
| `expression_matrices.R` | `combine_expression_matrices()` | `(expr_data_list, markers, max_cells=NULL, seed=42)` | Numeric matrix (cells x markers) |
| `differential_abundance.R` | `build_abundance_matrix()` | `(data_dir, cluster_regex, conditions_csv, condition_levels, min_cluster_size=10)` | list(counts, metadata, clusters) |
| `differential_abundance.R` | `plot_volcano()` | `(results_df, output_dir, title="Volcano Plot", p_col="padj", fc_col="log2fc", label_col="cluster")` | invisible(NULL) |
| `correlation_analysis.R` | `run_correlation_analysis()` | `(expr_matrix, method="spearman", output_dir=NULL, cluster_labels=NULL, p_adjust_method="BH")` | `list(cor_matrix, p_matrix, significant_pairs)` |
| `diffcyt_analysis.R` | `run_diffcyt()` | `(flow_frames, experiment_info, markers_of_interest, analysis_type="DA", method_DA="diffcyt-DA-edgeR", method_DS="diffcyt-DS-limma", n_clusters=20, seed=42)` | `list(results, top_table, cluster_assignments, method_used)` |
| `diffcyt_analysis.R` | `run_diffcyt_and_save()` | `(flow_frames, experiment_info, markers_of_interest, output_dir, analysis_type="DA", method_DA="diffcyt-DA-edgeR", method_DS="diffcyt-DS-limma", n_clusters=20, seed=42)` | `list(results, top_table, cluster_assignments, method_used)` |
| `differential_expression.R` | `run_differential_expression()` | `(input_dir, conditions_csv, output_dir, cluster_regex="Clustering annotated", markers=NULL, test_type="lmer", p_adjust_method="BH", seed=123)` | data.frame(cluster, marker, log2fc, pvalue, padj, test) |
| `effect_sizes.R` | `compute_effect_sizes()` | `(expr_data, condition_col, value_cols, methods=c("cohens_d","cliffs_delta"))` | data.frame(marker, effect_size, ci_lower, ci_upper, method, interpretation) |
| `effect_sizes.R` | `run_effect_sizes()` | `(input_csv, output_dir, condition_col="condition", value_cols=NULL, methods=c("cohens_d","cliffs_delta"))` | data.frame with effect sizes |
| `survival_analysis.R` | `run_survival_analysis()` | `(abundance_matrix, clinical_csv, output_dir, method="both", split_method="median", p_adjust_method="BH")` | data.frame(cluster, hr, ci_lower, ci_upper, pvalue, padj) |
| `cell_ratio.R` | `calculate_cell_ratios()` | `(abundance_csv, numerator_clusters, denominator_clusters, output_path=NULL)` | data.frame with ratios per sample |

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

**Example — Correlation analysis standalone:**

```r
setwd("/path/to/cytometry/src/analysis")
source("correlation_analysis.R")

result <- run_correlation_analysis(
  expr_matrix = combined_matrix,
  method      = "spearman",
  output_dir  = "/path/to/output"
)
# result$cor_matrix        -> correlation matrix
# result$significant_pairs -> data.frame of significant pairs (padj < 0.05)
```

**Example — diffcyt analysis standalone:**

```r
setwd("/path/to/cytometry/src/analysis")
source("diffcyt_analysis.R")

experiment_info <- data.frame(
  sample_id = c("S1", "S2", "S3", "S4"),
  group_id  = c("Control", "Control", "Treatment", "Treatment"),
  patient_id = c("P1", "P2", "P3", "P4")
)

result <- run_diffcyt(
  flow_frames         = my_flow_frames,
  experiment_info     = experiment_info,
  markers_of_interest = c("CD3", "CD4", "CD8"),
  analysis_type       = "DA",
  method_DA           = "diffcyt-DA-edgeR"
)
# result$top_table -> data.frame of differential results
```

**Example — Differential expression standalone:**

```r
setwd("/path/to/cytometry/src/analysis")
source("differential_expression.R")

results <- run_differential_expression(
  input_dir      = "/path/to/cluster/csvs",
  conditions_csv = "/path/to/conditions.csv",
  output_dir     = "/path/to/output",
  test_type      = "lmer"
)
# results -> data.frame(cluster, marker, log2fc, pvalue, padj, test)
```

**Example — Effect sizes standalone:**

```r
setwd("/path/to/cytometry/src/analysis")
source("effect_sizes.R")

results <- run_effect_sizes(
  input_csv     = "/path/to/expression_data.csv",
  output_dir    = "/path/to/output",
  condition_col = "condition",
  methods       = c("cohens_d", "cliffs_delta", "hedges_g")
)
```

**Example — Survival analysis standalone:**

```r
setwd("/path/to/cytometry/src/analysis")
source("survival_analysis.R")

results <- run_survival_analysis(
  abundance_matrix = abundance_df,
  clinical_csv     = "/path/to/clinical.csv",
  output_dir       = "/path/to/output",
  method           = "both",
  split_method     = "optimal"
)
# results -> data.frame(cluster, hr, ci_lower, ci_upper, pvalue, padj)
```

**Example — Cell ratio standalone:**

```r
setwd("/path/to/cytometry/src/analysis")
source("cell_ratio.R")

ratios <- calculate_cell_ratios(
  abundance_csv        = "/path/to/tracking_matrix.csv",
  numerator_clusters   = c("CD4_Tcells", "CD8_Tcells"),
  denominator_clusters = c("Tregs"),
  output_path          = "/path/to/output/ratios.csv"
)
```

---

## Visualization Scripts (`src/visualization/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `umap.R` | `plot_umap()` | `(flow_frames, sample_names, markers_of_interest=NULL, output_dir, n_neighbors=15, min_dist=0.1, output_formats=list(pdf=F, png=T, svg=T), seed=42)` | list(umap_coords, sample_labels, plots) |
| `umap.R` | `plot_umap_from_matrix()` | `(expr_matrix, sample_labels=NULL, cluster_labels=NULL, output_dir, n_neighbors=15, min_dist=0.1, seed=42, output_formats=list(pdf=F, png=T, svg=T))` | list(umap_coords, sample_labels, plots) |
| `heatmaps.R` | `plot_heatmap()` | `(expr_matrix, output_path, title, scale="row", width=10, height=8)` | invisible(NULL) |
| `heatmaps.R` | `plot_cluster_heatmap()` | `(expr_data, cluster_col, marker_cols, output_path, title)` | invisible(NULL) |
| `density_plots.R` | `plot_density()` | `(flow_frames, sample_names, markers, output_dir, plot_type="density", group_by="sample")` | invisible(NULL) |
| `density_plots.R` | `plot_density_from_matrix()` | `(expr_matrix, sample_labels=NULL, markers=NULL, output_dir, plot_type="density", group_by="sample")` | invisible(NULL) |
| `barplots.R` | `plot_composition_barplot()` | `(data, output_dir, x_col="cluster", fill_col="sample", position="stack", title)` | invisible(ggplot) |
| `tsne.R` | `plot_tsne()` | `(flow_frames, sample_names, markers_of_interest=NULL, output_dir, perplexity=30, max_iter=1000, marker_selection_mode="all", output_formats=list(pdf=F, png=T, svg=T), plot_width=10, plot_height=8, seed=42)` | list(tsne_coords, sample_labels, plots) |
| `tsne.R` | `plot_tsne_from_matrix()` | `(expr_matrix, sample_labels=NULL, cluster_labels=NULL, output_dir, perplexity=30, max_iter=1000, tsne_coords=NULL, output_formats=list(pdf=F, png=T, svg=T), plot_width=10, plot_height=8, seed=42)` | list(tsne_coords, sample_labels, plots) |
| `diffusion_map.R` | `plot_diffusion_map()` | `(expr_matrix, sample_labels=NULL, cluster_labels=NULL, output_dir, n_components=3, sigma="local", seed=42, output_formats=list(pdf=F, png=T, svg=T), plot_width=10, plot_height=8)` | list(dm_coords, eigenvalues, plots) |
| `violin_plots.R` | `plot_violin()` | `(expr_data, markers, group_col="cluster", plot_type="violin_box", output_dir, output_formats=list(pdf=F, png=T, svg=T), plot_width=12, plot_height=8)` | invisible(NULL) |
| `radar_plots.R` | `plot_radar()` | `(cluster_medians, markers, output_dir, normalize=TRUE, output_formats=list(pdf=F, png=T, svg=T), plot_width=10, plot_height=10)` | invisible(NULL) |
| `contour_plots.R` | `plot_contour()` | `(coords, group_labels=NULL, x_col="x", y_col="y", method="kde2d", output_dir, output_formats=list(pdf=F, png=T, svg=T), plot_width=10, plot_height=8)` | invisible(ggplot) |
| `biaxial_plots.R` | `plot_biaxial()` | `(flowframe=NULL, expr_matrix=NULL, x_marker, y_marker, color_by="density", sample_label=NULL, cluster_labels=NULL, gate_coords=NULL, output_dir, output_formats=list(pdf=F, png=T, svg=T), plot_width=10, plot_height=8)` | invisible(ggplot) |

**Example — t-SNE from matrix standalone:**

```r
setwd("/path/to/cytometry/src/visualization")
source("tsne.R")

result <- plot_tsne_from_matrix(
  expr_matrix    = combined_matrix,
  sample_labels  = sample_labels,
  cluster_labels = cluster_ids,
  output_dir     = "/path/to/output",
  perplexity     = 30,
  seed           = 42
)
# result$tsne_coords -> 2-column matrix of t-SNE coordinates
```

**Example — Diffusion map standalone:**

```r
setwd("/path/to/cytometry/src/visualization")
source("diffusion_map.R")

result <- plot_diffusion_map(
  expr_matrix    = combined_matrix,
  sample_labels  = sample_labels,
  cluster_labels = cluster_ids,
  output_dir     = "/path/to/output",
  n_components   = 3
)
# result$dm_coords   -> matrix with DC1, DC2, DC3 columns
# result$eigenvalues -> diffusion eigenvalues
```

**Example — Violin plots standalone:**

```r
setwd("/path/to/cytometry/src/visualization")
source("violin_plots.R")

plot_violin(
  expr_data  = my_data_frame,
  markers    = c("CD3", "CD4", "CD8", "CD19"),
  group_col  = "cluster",
  plot_type  = "violin_box",
  output_dir = "/path/to/output"
)
```

**Example — Biaxial plot standalone:**

```r
setwd("/path/to/cytometry/src/visualization")
source("biaxial_plots.R")

plot_biaxial(
  expr_matrix = combined_matrix,
  x_marker    = "CD3",
  y_marker    = "CD4",
  color_by    = "density",
  output_dir  = "/path/to/output"
)
```

**Example — Radar plot standalone:**

```r
setwd("/path/to/cytometry/src/visualization")
source("radar_plots.R")

# cluster_medians: data.frame with cluster rows and marker columns
plot_radar(
  cluster_medians = cluster_medians_df,
  markers         = c("CD3", "CD4", "CD8", "CD19", "CD56"),
  output_dir      = "/path/to/output",
  normalize       = TRUE
)
```

**Example — Contour plot standalone:**

```r
setwd("/path/to/cytometry/src/visualization")
source("contour_plots.R")

coords_df <- data.frame(x = umap_result[, 1], y = umap_result[, 2])
plot_contour(
  coords       = coords_df,
  group_labels = sample_labels,
  method       = "kde2d",
  output_dir   = "/path/to/output"
)
```

---

## Metadata QC Scripts (`src/metadata_qc/`)

| Script | Function | Signature |
|--------|----------|-----------|
| `check_metadata.R` | `check_metadata()` | `(metadata_file, ...)` |
| `check_acquisition.R` | `check_acquisition_dates()` | `(parent_folder, ...)` |
| `overlap_statistics.R` | `compute_fleiss_kappa()` | `(label_lists, ...)` |
| `compare_gatings.R` | `compare_gatings()` | `(pop_strings_list, ...)` |
