[< Back to Documentation Index](README.md)

# Input and Output File Formats

This reference describes every file format accepted as input and produced as output by the cytometry pipeline.

---

## Input Formats

### FCS Files (Pipelines A, E)

Standard FCS 3.0/3.1 files from flow cytometers. Marker names are read from `$PnS` keywords. All `.fcs` files in the input directory are loaded.

Key keywords the pipeline uses:

| Keyword | Example | Usage |
|---------|---------|-------|
| `$PnS` | `CD3`, `CD4`, `HLA-DR` | Marker (stain) names — used as column headers |
| `$PnN` | `BV421-A`, `APC-A` | Channel (fluorochrome) names |
| `$TOT` | `50000` | Total events — validated by `check_tot()` |
| `$DATE` | `15-MAR-2024` | Acquisition date — extracted by `check_acquisition.R` |
| `$FIL` | `Specimen_001.fcs` | The name the file was acquired under, and **what FlowJo displays by default**. Renaming a file on disk does not change it — see [Sync FCS Embedded Names](scripts/sync-fcs-names.md) |
| `ORIGINALFIL` | `Specimen_001.fcs` | Written by `sync_fcs_names.R` when it replaces `$FIL`, so the acquisition name survives the rename. Write-once. Not an FCS-standard keyword, and deliberately un-prefixed: `$` is reserved for the standard's own vocabulary (FCS 3.1 §3.2.12) |

### Expression CSVs (Pipeline B)

Semicolon-separated CSVs with one row per cell and one column per marker. No gating columns, no metadata columns — pure expression data. The separator is auto-detected (semicolon, comma, or tab).

Filename convention: `{sample_id}_{population}.csv`

```
CD3;CD4;CD8;CD19;CD56;HLA-DR
3.94;4.31;5.33;3.23;2.35;2.78
2.75;4.80;3.52;3.42;2.00;4.06
```

### Gating CSVs (Pipeline C)

Semicolon-separated CSVs with expression columns **and** binary gating columns (0 = not in gate, 1 = in gate). Typically exported from Cytolution or similar tools. Each CSV should have a matching `.fcs` file in the same directory for row-count validation.

```
CD3;CD4;CD8;CD19;CD56;HLA-DR;CD4+ (Clustering annotated);CD8+ (Clustering annotated);Clustering raw
2.46;4.01;2.15;1.29;2.37;1.22;0;0;1
3.89;2.12;4.56;1.87;3.21;3.45;1;0;1
```

The column classifier automatically detects six column types:

| Type | Pattern | Example |
|------|---------|---------|
| Index | `index`, `Index` | `index` |
| Scatter/Time | FSC, SSC, Time keywords | `SSC.H \| SSC.H`, `TIME \| Time` |
| Expression | Marker \| Channel format | `CD4 \| BV421.A` |
| Gating | Ends with `(Clustering annotated)`, `(Manual gating)`, or is `Clustering raw` | `CD4+ T cells (Clustering annotated)` |
| Dim. reduction | `cyto_2d_*` | `cyto_2d_Samples_All_X` |
| Autofluorescence | Marker name = channel name | `AF.P1 \| AF.P1.A` |

### FlowJo Workspaces (Pipeline D)

`.wsp` files with their associated FCS files. The pipeline uses `CytoML::open_flowjo_xml()` and `flowWorkspace::parseWorkspace()` to extract gating hierarchies and population assignments.

### Conditions CSV (Pipelines A, C)

Maps sample filenames to experimental groups for differential testing. Semicolon-separated, two required columns:

```
file_name;condition
sample_001;Control
sample_002;Control
sample_003;Treatment
```

### Metadata CSV (Pipeline D, abundance plotting)

Study metadata for grouping plots. Semicolon-separated, with a `Sample_ID` column for matching:

```
Sample_ID;Datum;Age;Sex;Diagnosis;Treatment_Group;Response;Batch
ICK-001;15/03/2024;45;M;AML;Chemo;pCR;Batch1
ICK-002;15/03/2024;62;F;AML;Chemo;npCR;Batch1
```

### Cluster Tracking Matrix

Output of `split_by_cluster()` and `build_tracking_matrix()`. Used as input for abundance plotting. One row per sample, one column per cluster:

```
file;CD4+ T cells;CD8+ T cells;B cells;NK cells;Monocytes
sample_001.csv;1234;567;890;234;456
sample_002.csv;987;1234;345;567;234
```

---

## Output Formats

### CSV Outputs

All CSV outputs use **semicolon** (`;`) as separator, consistent with the lab's conventions.

| File | Produced By | Columns |
|------|-------------|---------|
| `leiden_clusters.csv` | `run_leiden_with_umap()` | UMAP1, UMAP2, cluster, sample |
| `glmm_results.csv` | `run_differential_abundance()` | cluster, estimate, std_error, z_value, p_value, log2fc, odds_ratio, method, p_adjusted, significant |
| `enrichment_results.csv` | `run_cluster_enrichment()` | cluster, in_A, not_in_A, in_B, not_in_B, odds_ratio, p_value, ci_lower, ci_upper, method, p_adjusted, significant, log2_enrichment |
| `cluster_tracking_matrix.csv` | `split_by_cluster()` / `build_tracking_matrix()` | file, {cluster1}, {cluster2}, ... |
| `marker_statistics.csv` | `run_marker_statistics()` | population, marker, mean, between_sample_var, within_sample_var, icc, n_cells, n_samples, method |

### Plot Outputs

Plots are saved in multiple formats controlled by the `output_formats` config parameter:

```r
output_formats = list(pdf = TRUE, png = TRUE, svg = FALSE)
```

| Plot | Produced By | Description |
|------|-------------|-------------|
| `umap_sample.png` | `plot_umap()` | UMAP scatter plot colored by sample origin |
| `umap_{marker}.png` | `plot_umap()` | UMAP scatter plot colored by marker expression (viridis palette) |
| `volcano_plot.png` | `run_differential_abundance()` | log2 fold-change vs -log10 adjusted p-value |
| `expression_heatmap.png` | `plot_heatmap()` | Clustered heatmap of expression (pheatmap) |
| `barplot_*.png` | `plot_composition_barplot()` | Cell composition stacked/dodge barplot |
| `abundance_stacked.png` | `plot_abundance_stacked()` | Cluster abundance per sample, sorted by total |
| `abundance_proportional.png` | `plot_abundance_proportional()` | Proportional barplot with Ward.D2 dendrogram |
| `total_cells_per_cluster.png` | `plot_abundance_totals()` | Total cells per cluster across all samples |
| `cluster_proportion_heatmap.png` | `plot_abundance_heatmap()` | Biclustered proportion heatmap (inferno palette) |
