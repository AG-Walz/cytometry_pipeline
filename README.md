# Cytometry Data Analysis Pipeline

A modular flow cytometry analysis framework in R, supporting raw FCS files, preprocessed FCS, CSV expression matrices, FlowJo workspace output, and CSV gating data. Designed for immune cell phenotyping with Leiden clustering, differential abundance analysis (NB-GLMM), and enrichment statistics.

---

## Table of Contents

1. [Quick Start](#1-quick-start)
2. [Choosing a Pipeline](#2-choosing-a-pipeline)
3. [Pipeline A — Raw FCS End-to-End](#3-pipeline-a--raw-fcs-end-to-end)
4. [Pipeline B — CSV Expression Statistics](#4-pipeline-b--csv-expression-statistics)
5. [Pipeline C — CSV Gating Cluster Analysis](#5-pipeline-c--csv-gating-cluster-analysis)
6. [Pipeline D — FlowJo Workspace Comparison](#6-pipeline-d--flowjo-workspace-comparison)
7. [Pipeline E — Preprocessed FCS Quick Clustering](#7-pipeline-e--preprocessed-fcs-quick-clustering)
8. [Input File Formats](#8-input-file-formats)
9. [Output File Formats](#9-output-file-formats)
10. [Standalone Scripts Reference](#10-standalone-scripts-reference)
11. [Shared Utility Library (R_shared)](#11-shared-utility-library-r_shared)
12. [Configuration Guide](#12-configuration-guide)
13. [Environment Setup](#13-environment-setup)
14. [Project Structure](#14-project-structure)
15. [Testing](#15-testing)
16. [Key Packages](#16-key-packages)
17. [Troubleshooting](#17-troubleshooting)

---

## 1. Quick Start

```bash
# 1. Clone and enter the project
git clone https://github.com/AG-Walz/cytometry_pipeline.git
cd cytometry_pipeline

# 2. Set up the environment
conda env create -f environment.yml
conda activate cytometry_env
Rscript install_extra_packages.R

# 3. Generate test data (creates synthetic FCS + CSV files)
Rscript test_data/generate_test_data.R

# 4. Run a pipeline (from the pipelines/ directory)
cd pipelines
Rscript pipeline_A_raw_fcs.R --config ../config/pipeline_A_config.R
```

After the run completes, explore the output:

```
output/pipeline_A_results/
├── 03_clustering/
│   └── leiden_clusters.csv          # UMAP coordinates + cluster labels
├── 04_analysis/
│   ├── cluster_data/                # Per-sample CSVs with cluster columns
│   ├── split_clusters/              # Per-cluster split files + tracking matrix
│   └── differential_abundance/      # GLMM results + volcano plot
└── 05_visualization/
    ├── umap_sample.png              # UMAP colored by sample
    ├── umap_CD3.png                 # UMAP colored by CD3 expression
    ├── umap_CD4.png                 # ... one per marker
    └── ...
```

---

## 2. Choosing a Pipeline

The framework provides five pipelines. Choose based on your input data format:

```
                         ┌──────────────────────────────────────────────┐
                         │          What input data do you have?        │
                         └─────────────────────┬────────────────────────┘
                                               │
                 ┌─────────────────┬───────────┼───────────┬──────────────────┐
                 │                 │           │           │                  │
           Raw .fcs files   Preprocessed   CSV with     CSV with       FlowJo .wsp
           (untransformed)  .fcs files     expression   gating         workspaces
                 │                │        columns only  columns            │
                 ▼                ▼           ▼           ▼                  ▼
            Pipeline A      Pipeline E   Pipeline B   Pipeline C       Pipeline D
            (full analysis) (skip QC)    (statistics) (cluster DA)     (comparison)
```

| Pipeline | Input | What It Does | When to Use |
|----------|-------|-------------|-------------|
| **A** | Raw `.fcs` files | QC → Transform → Harmonize → Downsample → Leiden → DA → UMAP | You have raw FCS from the cytometer and want full end-to-end analysis |
| **B** | Expression CSVs | Load → Marker statistics → Heatmap | You exported expression matrices from FlowJo/Cytolution and want per-marker stats |
| **C** | CSVs with gating columns | Load → Split clusters → Enrichment/GLMM → Abundance plots | You have Cytolution or manual gating exports with binary gate columns |
| **D** | FlowJo `.wsp` files | Compare gating strategies → Kappa → Venn/UpSet | You want to compare multiple analysts' gating of the same samples |
| **E** | Preprocessed `.fcs` | Harmonize → Downsample → Leiden → UMAP | Your FCS files are already QC'd and transformed; you just want clustering |

### Running a Pipeline

All pipelines follow the same pattern:

```bash
cd pipelines
Rscript pipeline_X_name.R --config ../config/pipeline_X_config.R
```

Or from R:

```r
source("config/pipeline_A_config.R")   # loads params list
source("pipelines/pipeline_A_raw_fcs.R")
```

---

## 3. Pipeline A — Raw FCS End-to-End

**Full guide:** [docs/pipelines/pipeline-A-raw-fcs.md](docs/pipelines/pipeline-A-raw-fcs.md)

### What It Does

Takes raw `.fcs` files from the cytometer through every processing step: quality control, mathematical transformation, panel harmonization, downsampling, Leiden clustering, differential abundance testing, and UMAP visualization.

### Flow Diagram

```
  .fcs files (one per sample)
       │
       ▼
  ┌─────────────────────┐
  │ 1. load_fcs_raw()   │  Read FCS files into flowFrame objects
  └─────────┬───────────┘
            ▼
  ┌─────────────────────┐
  │ 2. run_preprocessing│  QC (PeacoQC/flowAI) + arcsinh/logicle transform
  └─────────┬───────────┘
            ▼
  ┌─────────────────────────────┐
  │ 3. homogenize_flowframes()  │  Rename markers, reorder to common panel
  └─────────┬───────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 4. downsample_flowframes()   │  Reduce to N cells per sample
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 5. run_leiden_with_umap()    │  KNN graph → Leiden clusters + UMAP
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────────────────┐
  │ 6. split_by_cluster() +                  │  Per-cluster CSVs + volcano plot
  │    run_differential_abundance()           │
  └─────────┬────────────────────────────────┘
            ▼
  ┌─────────────────────┐
  │ 7. plot_umap()      │  UMAP plots per sample + per marker
  └─────────────────────┘
```

### Example Config

```r
# config/pipeline_A_config.R (key parameters)
params <- list(
  fcs_folder          = "/path/to/fcs/files",
  output_dir          = "/path/to/output",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),
  conditions_csv      = "/path/to/conditions.csv",
  condition_levels    = c("Control", "Treatment"),
  truncation          = FALSE,
  pp_order            = "QT",           # QC before Transform
  qc_type             = "None",         # "PeacoQC", "flowAI", or "None"
  transform_type      = "arcsinh",
  offset_list         = rep(0, 6),      # one per marker
  width_list          = rep(5, 6),      # cofactor per marker
  downsample_method   = "Consistent",
  downsample_n        = 10000,
  n_neighbors         = 15,
  resolution          = 0.5,
  seed                = 42
)
```

### Example Output

After running with the test data (3 samples, 50 cells each, 6 markers):

**Leiden cluster assignments** (`03_clustering/leiden_clusters.csv`):

```
"UMAP1","UMAP2","cluster","sample"
-0.134,-1.288,1,"sample_001"
-0.039,0.814,1,"sample_001"
1.706,-0.682,1,"sample_001"
...
```

**GLMM differential abundance results** (`04_analysis/differential_abundance/glmm_results.csv`):

```
cluster;estimate;std_error;z_value;p_value;log2fc;odds_ratio;method;p_adjusted;significant
Cluster_1;-0.23;0.15;-1.53;0.126;-0.33;0.79;glmerMod;0.126;FALSE
Cluster_2;0.87;0.21;4.14;0.00003;1.26;2.39;glmerMod;0.0001;TRUE
```

**Cluster tracking matrix** (`04_analysis/split_clusters/cluster_tracking_matrix.csv`):

```
file;Cluster_1;Cluster_2;Cluster_3
sample_001;33;12;5
sample_002;28;18;4
sample_003;35;8;7
```

**UMAP plots** (`05_visualization/`): one PNG per marker showing expression intensity overlaid on the UMAP embedding, plus `umap_sample.png` colored by sample origin.

### Output Directory Structure

```
output/pipeline_A_results/
├── 01_preprocessing/                    # (if save_pp_fcs/csv enabled)
├── 02_harmonization/                    # (if save_csv/fcs enabled)
├── 03_clustering/
│   └── leiden_clusters.csv              # UMAP1, UMAP2, cluster, sample
├── 04_analysis/
│   ├── cluster_data/                    # Per-sample CSVs with binary cluster columns
│   │   ├── sample_001.csv
│   │   ├── sample_002.csv
│   │   └── sample_003.csv
│   ├── split_clusters/
│   │   ├── cluster_tracking_matrix.csv  # Sample-by-cluster cell counts
│   │   ├── sample_001_Cluster_1.csv     # Expression data for cells in cluster 1
│   │   └── ...
│   └── differential_abundance/
│       ├── glmm_results.csv             # Per-cluster GLMM test results
│       └── volcano_plot.png             # log2FC vs -log10(p-adjusted)
└── 05_visualization/
    ├── umap_sample.png                  # UMAP colored by sample
    ├── umap_CD3.png                     # UMAP colored by CD3 expression
    ├── umap_CD4.png
    └── ...                              # One plot per marker
```

---

## 4. Pipeline B — CSV Expression Statistics

**Full guide:** [docs/pipelines/pipeline-B-csv-expression.md](docs/pipelines/pipeline-B-csv-expression.md)

### What It Does

The simplest pipeline. Loads pre-exported expression CSV files (e.g. from FlowJo or Cytolution), computes per-marker statistics using linear mixed-effects models, and generates a heatmap.

### Flow Diagram

```
  CSV files (one per sample × population)
       │
       ▼
  ┌────────────────────────────┐
  │ 1. load_csv_expression()   │  Read CSVs, auto-detect separator
  └─────────┬──────────────────┘
            ▼
  ┌────────────────────────────┐
  │ 2. run_marker_statistics() │  Per-marker lmer: mean, variance, ICC
  └─────────┬──────────────────┘
            ▼
  ┌────────────────────────────┐
  │ 3. plot_heatmap()          │  Clustered expression heatmap
  └────────────────────────────┘
```

### Expected Input

A directory of semicolon-separated CSVs. Each file contains expression columns only (no gating columns). Filenames encode sample and population: `{sample}_{population}.csv`.

```
csv_expression/
├── sample_001_CD4_T_cells.csv
├── sample_001_CD8_T_cells.csv
├── sample_002_CD4_T_cells.csv
├── sample_002_CD8_T_cells.csv
├── sample_003_CD4_T_cells.csv
└── sample_003_CD8_T_cells.csv
```

Each CSV looks like:

```
CD3;CD4;CD8;CD19;CD56;HLA-DR
3.94;4.31;5.33;3.23;2.35;2.78
2.75;4.80;3.52;3.42;2.00;4.06
1.89;3.12;2.67;4.01;3.55;1.92
...
```

### Example Config

```r
# config/pipeline_B_config.R
params <- list(
  csv_folder = "/path/to/csv/expression/files",
  output_dir = "/path/to/output",
  markers    = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),
  test_type  = "lmer",    # "lmer" (mixed-effects), "wilcox", or "both"
  seed       = 42
)
```

### Example Output

**Marker statistics** (`01_statistics/marker_statistics.csv`):

```
population;marker;mean;between_sample_var;within_sample_var;icc;n_cells;n_samples;method
CD4_T_cells;CD3;3.02;0.00;0.97;0.00;300;6;lmer
CD4_T_cells;CD4;3.08;0.003;1.10;0.003;300;6;lmer
CD4_T_cells;CD8;2.95;0.024;1.04;0.022;300;6;lmer
...
```

Columns explained:
- **population**: Cell population (parsed from filename)
- **mean**: Grand mean expression across all cells and samples
- **between_sample_var / within_sample_var**: Variance components from the mixed model
- **icc**: Intraclass correlation coefficient (how much variance is between-sample vs within)
- **method**: `lmer` (mixed-effects model succeeded) or `simple_mean` (fallback)

**Expression heatmap** (`02_visualization/expression_heatmap.png`): A clustered heatmap (pheatmap) showing expression patterns across markers.

### Output Directory Structure

```
output/pipeline_B_results/
├── 01_statistics/
│   └── marker_statistics.csv     # Per-population, per-marker statistics
└── 02_visualization/
    ├── expression_heatmap.png    # Clustered heatmap
    └── expression_heatmap.svg
```

---

## 5. Pipeline C — CSV Gating Cluster Analysis

**Full guide:** [docs/pipelines/pipeline-C-csv-gating.md](docs/pipelines/pipeline-C-csv-gating.md)

### What It Does

Processes CSV exports that contain both expression columns and binary gating columns (e.g. from Cytolution). Splits data by cluster, tests for enrichment between conditions, and generates abundance plots.

### Flow Diagram

```
  CSV + FCS file pairs (from Cytolution/manual gating)
       │
       ▼
  ┌──────────────────────────────┐
  │ 1. load_csv_gating()         │  Classify columns, split per population
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 2. split_by_cluster()        │  Count cells per cluster per sample
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 3. run_cluster_enrichment()  │  Fisher's exact test between conditions
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 4. plot_composition_barplot()│  Cluster composition barplot
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 5. build_tracking_matrix()   │  Abundance stacked/proportional/heatmap
  │    plot_abundance()          │
  └──────────────────────────────┘
```

### Expected Input

A directory containing paired CSV + FCS files. CSVs have expression columns plus binary gating columns (0/1 membership). Gating columns are identified by a regex pattern.

```
csv_gating/
├── sample_001.csv        # Expression + gating columns
├── sample_001.fcs        # Matching FCS file (for row-count validation)
├── sample_002.csv
├── sample_002.fcs
├── sample_003.csv
└── sample_003.fcs
```

Each CSV looks like:

```
CD3;CD4;CD8;CD19;CD56;HLA-DR;CD4+ (Clustering annotated);CD8+ (Clustering annotated);Clustering raw
2.46;4.01;2.15;1.29;2.37;1.22;0;0;1
2.57;4.40;3.12;4.38;4.07;2.33;0;0;0
3.89;2.12;4.56;1.87;3.21;3.45;1;0;1
0.98;1.56;2.78;3.45;1.23;2.67;0;1;0
...
```

The column auto-classifier recognizes these patterns:
- `(Clustering annotated)` → gating column
- `(Manual gating)` → gating column
- `Clustering raw` → gating column
- Everything else numeric → expression marker

A separate **conditions file** maps samples to experimental groups:

```
file_name;condition
sample_001;Control
sample_002;Control
sample_003;Treatment
```

### Example Config

```r
# config/pipeline_C_config.R
params <- list(
  csv_folder       = "/path/to/csv/gating/files",
  output_dir       = "/path/to/output",
  conditions_csv   = "/path/to/conditions.csv",
  condition_levels = c("Control", "Treatment"),
  analysis_method  = "enrichment",     # "enrichment" (Fisher) or "glmm" (NB-GLMM)
  stat_method      = "pooled",         # "pooled" or "mixed_effects"
  cluster_regex    = "Clustering annotated",
  gating_pattern   = "\\(Clustering annotated\\)$|^Clustering raw$",
  processing_mode  = "manual",         # "automatic" (GUI), "manual", "expression_only"
  exclude_patterns = c(),
  seed             = 42
)
```

### Example Output

**Enrichment results** (`03_analysis/enrichment_results.csv`):

```
cluster;in_A;not_in_A;in_B;not_in_B;odds_ratio;p_value;ci_lower;ci_upper;method;p_adjusted;significant;log2_enrichment
CD4+;40;60;20;30;1.00;1.00;0.47;2.13;fisher;1.00;FALSE;0.00
CD8+;30;70;15;35;1.00;1.00;0.45;2.27;fisher;1.00;FALSE;0.00
```

Columns explained:
- **in_A / not_in_A**: Cells in/out of cluster in condition A (Control)
- **in_B / not_in_B**: Cells in/out of cluster in condition B (Treatment)
- **odds_ratio**: Fisher's exact test odds ratio
- **log2_enrichment**: log2 of the proportion ratio between conditions

**Cluster tracking matrix** (`05_abundance/cluster_tracking_matrix.csv`):

```
file;CD4+;CD8+
sample_001;20;15
sample_002;20;15
sample_003;20;15
```

**Abundance plots** (`05_abundance/`): Stacked barplot (sorted by total cells), proportional barplot with dendrogram, total cells per cluster, proportion heatmap.

### Output Directory Structure

```
output/pipeline_C_results/
├── 01_loaded/                           # Per-population expression CSVs
│   ├── sample_001_CD4_.csv
│   ├── sample_001_CD8_.csv
│   ├── sample_001_Clustering_raw.csv
│   └── ...
├── 02_clusters/
│   ├── cluster_tracking_matrix.csv      # Sample × cluster cell counts
│   ├── sample_001_CD4+.csv             # Expression data for CD4+ cells
│   └── ...
├── 03_analysis/
│   └── enrichment_results.csv           # Fisher/mixed-effects results
├── 04_visualization/
│   └── barplot_cluster_file.png         # Composition barplot
└── 05_abundance/
    ├── cluster_tracking_matrix.csv
    ├── abundance_stacked.png            # Stacked barplot by total count
    ├── abundance_proportional.png       # Proportional barplot + dendrogram
    ├── total_cells_per_cluster.png      # Bar chart of cluster sizes
    └── cluster_proportion_heatmap.png   # Biclustered proportion heatmap
```

---

## 6. Pipeline D — FlowJo Workspace Comparison

**Full guide:** [docs/pipelines/pipeline-D-flowjo-comparison.md](docs/pipelines/pipeline-D-flowjo-comparison.md)

### What It Does

Compares gating strategies across multiple FlowJo workspaces (e.g. from different analysts). Extracts population assignments, computes inter-rater agreement (Fleiss' kappa), and visualizes overlap using Venn/Euler diagrams and UpSet plots.

### Flow Diagram

```
  .wsp files (one per analyst)  +  shared .fcs files
       │
       ▼
  ┌────────────────────────────┐
  │ 1. load_flowjo_wsp()       │  Extract gating sets from each WSP
  └─────────┬──────────────────┘
            ▼
  ┌────────────────────────────┐
  │ 2. compare_gatings()       │  Align populations across workspaces
  └─────────┬──────────────────┘
            ▼
  ┌────────────────────────────┐
  │ 3. compute_fleiss_kappa()  │  Inter-rater agreement statistics
  └─────────┬──────────────────┘
            ▼
  ┌────────────────────────────┐
  │ 4. plot_venn() / plot_euler│  Overlap visualizations
  │    plot_upset()            │
  └────────────────────────────┘
```

### Expected Input

- Two or more FlowJo `.wsp` workspace files (same samples, different gating strategies)
- The `.fcs` files referenced by the workspaces

### Note

Pipeline D requires real FlowJo workspace files and cannot be tested with synthetic data. See the [full guide](docs/pipelines/pipeline-D-flowjo-comparison.md) for details.

---

## 7. Pipeline E — Preprocessed FCS Quick Clustering

**Full guide:** [docs/pipelines/pipeline-E-preprocessed-umap.md](docs/pipelines/pipeline-E-preprocessed-umap.md)

### What It Does

A streamlined version of Pipeline A for FCS files that have already been QC'd and transformed. Skips preprocessing and goes directly to panel harmonization, downsampling, Leiden clustering, and UMAP visualization.

### Flow Diagram

```
  Preprocessed .fcs files
       │
       ▼
  ┌─────────────────────┐
  │ 1. load_fcs_raw()   │  Load FCS (truncation=FALSE, no QC)
  └─────────┬───────────┘
            ▼
  ┌─────────────────────────────┐
  │ 2. homogenize_flowframes()  │  Unify marker panels
  └─────────┬───────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 3. downsample_flowframes()   │  Reduce cell count
  └─────────┬────────────────────┘
            ▼
  ┌──────────────────────────────┐
  │ 4. run_leiden_with_umap()    │  Cluster + embed
  └─────────┬────────────────────┘
            ▼
  ┌─────────────────────┐
  │ 5. plot_umap()      │  UMAP per sample + per marker
  └─────────────────────┘
```

### Example Config

```r
# config/pipeline_E_config.R
params <- list(
  fcs_folder          = "/path/to/preprocessed/fcs",
  output_dir          = "/path/to/output",
  markers_of_interest = c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR"),
  markername_replacing = list(),     # e.g. list(`PD1 (CD279)` = "PD1")
  downsample_method   = "Consistent",
  downsample_n        = 10000,
  n_neighbors         = 15,
  resolution          = 0.5,
  seed                = 42
)
```

### Output Directory Structure

```
output/pipeline_E_results/
├── 01_harmonized/           # (if save_csv/fcs enabled)
├── 02_clustering/
│   └── leiden_clusters.csv  # UMAP1, UMAP2, cluster, sample
└── 03_visualization/
    ├── umap_sample.png      # UMAP colored by sample
    ├── umap_CD3.png         # UMAP colored by CD3 expression
    └── ...                  # One plot per marker
```

---

## 8. Input File Formats

### FCS Files (Pipelines A, E)

Standard FCS 3.0/3.1 files from flow cytometers. Marker names are read from `$PnS` keywords. All `.fcs` files in the input directory are loaded.

### Expression CSVs (Pipeline B)

Semicolon-separated CSVs with one row per cell and one column per marker. No gating columns, no metadata columns — pure expression data. The separator is auto-detected (semicolon, comma, or tab).

Filename convention for marker_statistics: `{sample_id}_{population}.csv`

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

### Conditions CSV (Pipelines A, C)

Maps sample filenames to experimental groups for differential testing. Semicolon-separated, two required columns:

```
file_name;condition
sample_001;Control
sample_002;Control
sample_003;Treatment
```

### Metadata CSV (Pipeline D, abundance plotting)

Study metadata for grouping plots. Semicolon-separated, with a `*file_name` column for matching:

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
sample_003.csv;1456;432;1023;189;567
```

---

## 9. Output File Formats

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

---

## 10. Standalone Scripts Reference

Every script in `src/` can be used independently outside of a pipeline. Each script defines one or more functions that you call after sourcing.

### Loading Scripts (`src/loading/`)

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

### Preprocessing Scripts (`src/preprocessing/`)

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
  method       = "Consistent",   # Same number per sample
  seed         = 42
)
# Console output:
#   Downsampling 3 samples (total: 150000 cells)
#   Method: Consistent | Target: 5000 cells
#   Consistent: taking 1666 cells per sample
#     S1 : 50000 -> 1666
#     S2 : 50000 -> 1666
#     S3 : 50000 -> 1666
#   Total after downsampling: 4998 cells
```

### Harmonization Scripts (`src/harmonization/`)

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
# Console output:
#   Homogenizing 3 samples to 4 markers
#     Homogenizing S1 ...  50000 cells x 4 markers
#     Homogenizing S2 ...  48000 cells x 4 markers
#     Homogenizing S3 ...  52000 cells x 4 markers
```

### Clustering Scripts (`src/clustering/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `leiden.R` | `run_leiden()` | `(expr_matrix, n_neighbors=30, resolution=1.0, seed=42)` | Integer vector of cluster IDs |
| `leiden.R` | `run_leiden_with_umap()` | `(expr_matrix, sample_labels=NULL, output_dir=NULL, n_neighbors=30, resolution=1.0, seed=42)` | data.frame(UMAP1, UMAP2, cluster, sample) |

**Example — Leiden clustering standalone:**

```r
setwd("/path/to/cytometry/src/clustering")
source("leiden.R")

result <- run_leiden_with_umap(
  expr_matrix   = combined_expression_matrix,  # cells × markers
  sample_labels = sample_labels_vector,
  output_dir    = "/path/to/output",
  n_neighbors   = 15,
  resolution    = 0.8,
  seed          = 42
)
# Console output:
#   Computing UMAP...
#   Building KNN graph (k = 15)...
#   Running Leiden clustering (resolution = 0.8)...
#   Found 12 clusters
#   Cluster sizes: 1234, 987, 876, 765, ...
#   Results saved to: /path/to/output

# result is a data.frame:
#   UMAP1       UMAP2      cluster  sample
#   -0.134      -1.288     1        sample_001
#   -0.039       0.814     3        sample_001
#   ...
```

### Analysis Scripts (`src/analysis/`)

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
# Console output:
#   Loading 42 files...
#   Total cells with conditions: 125000
#   Cluster columns: 8
#   Clusters passing filter: 6
#     Analyzing: CD4+ ... p = 0.023
#     Analyzing: CD8+ ... p = 0.456
#     ...
#   Enrichment analysis complete.
#   Significant clusters: 2
```

### Visualization Scripts (`src/visualization/`)

| Script | Function | Signature | Returns |
|--------|----------|-----------|---------|
| `umap.R` | `plot_umap()` | `(flow_frames, sample_names, markers_of_interest=NULL, output_dir, n_neighbors=15, min_dist=0.1, output_formats=list(pdf=F, png=T, svg=T), seed=42)` | list(umap_coords, sample_labels, plots) |
| `heatmaps.R` | `plot_heatmap()` | `(expr_matrix, output_path, title, scale="row", width=10, height=8)` | invisible(NULL) |
| `heatmaps.R` | `plot_cluster_heatmap()` | `(expr_data, cluster_col, marker_cols, output_path, title)` | invisible(NULL) |
| `density_plots.R` | `plot_density()` | `(flow_frames, sample_names, markers, output_dir, plot_type="density", group_by="sample")` | invisible(NULL) |
| `barplots.R` | `plot_composition_barplot()` | `(data, output_dir, x_col="cluster", fill_col="sample", position="stack", title)` | invisible(ggplot) |

**Example — UMAP standalone:**

```r
setwd("/path/to/cytometry/src/visualization")
source("umap.R")

result <- plot_umap(
  flow_frames         = my_flow_frames,
  sample_names        = c("S1", "S2", "S3"),
  markers_of_interest = c("CD3", "CD4", "CD8"),
  output_dir          = "/path/to/output/umap",
  n_neighbors         = 15,
  seed                = 42
)
# Creates: umap_sample.png, umap_CD3.png, umap_CD4.png, umap_CD8.png
```

### Metadata QC Scripts (`src/metadata_qc/`)

| Script | Function | Signature |
|--------|----------|-----------|
| `check_metadata.R` | `check_metadata()` | `(metadata_file, ...)` |
| `check_acquisition.R` | `check_acquisition_dates()` | `(parent_folder, ...)` |
| `overlap_statistics.R` | `compute_fleiss_kappa()` | `(label_lists, ...)` |
| `compare_gatings.R` | `compare_gatings()` | `(pop_strings_list, ...)` |

---

## 11. Shared Utility Library (R_shared)

All subproject scripts source shared utilities from `src/R_shared/`. These 8 modules provide the building blocks used throughout the framework.

**Full API reference:** [docs/r-shared-reference.md](docs/r-shared-reference.md)

| Module | Key Functions | Purpose |
|--------|--------------|---------|
| **io_utils.R** | `detect_csv_sep()`, `classify_csv_columns()`, `list_csvs()`, `make_outdir()`, `die()` | File I/O, directory creation, CSV separator detection, column classification |
| **gui_utils.R** | `select_markers()`, `select_columns_gui()`, `is_gui_available()` | Interactive marker/column selection with automatic CLI fallback for headless servers |
| **memory_utils.R** | `check_memory()`, `read_csv_chunked()`, `get_available_memory_gb()` | Memory monitoring (triggers GC at 85%), chunked CSV reading for large files |
| **cluster_utils.R** | `clean_cluster_name()`, `labels_from_binary_df()`, `translate_labels()` | Parse cluster names from column headers, convert binary gating to categorical labels |
| **plot_utils.R** | `save_plot_formats()`, `save_base_plot_formats()`, `create_safe_filename()` | Save ggplot/base plots to PDF/PNG/SVG, sanitize filenames |
| **flowframe_utils.R** | `explore_data()`, `extract_exprmtx_list()`, `from_mtx_to_ffr()`, `compare_panels()` | Inspect flowFrame markers, extract expression matrices, convert between matrix and flowFrame |
| **fcs_utils.R** | `check_tot()`, `create_csv_fcs_mapping()`, `get_fcs_event_count()` | FCS file validation, CSV↔FCS file matching |
| **stats_utils.R** | `run_umap()`, `run_pca()`, `fit_nb_glmm()`, `calculate_enrichment()`, `do_mixed_effects()` | Dimensionality reduction, statistical models (NB-GLMM, Fisher, mixed-effects logistic) |

### Sourcing Pattern

Every script begins by sourcing only the shared modules it needs:

```r
shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "flowframe_utils.R"))
# ... source only what you need
```

---

## 12. Configuration Guide

### Pipeline Configs

Each pipeline has a dedicated config in `config/`. Copy and edit for your study:

```bash
cp config/pipeline_A_config.R config/study_params/my_study_A.R
# Edit paths and parameters, then run:
cd pipelines
Rscript pipeline_A_raw_fcs.R --config ../config/study_params/my_study_A.R
```

Config files automatically compute `base_dir` from their own file location, so all relative paths "just work" regardless of where you run from.

### Key Parameters Across Pipelines

| Parameter | Used By | Description |
|-----------|---------|-------------|
| `fcs_folder` | A, E | Directory containing `.fcs` files |
| `csv_folder` | B, C | Directory containing CSV files |
| `output_dir` | All | Root output directory (subdirs created automatically) |
| `markers_of_interest` | A, E | Character vector of marker names to analyze |
| `markers` | B, C | Character vector of marker column names |
| `conditions_csv` | A, C | Path to file_name → condition mapping CSV |
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

### Standalone Script Configs

For running individual scripts outside pipelines, see `config/standalone/` for ready-to-edit config templates.

---

## 13. Environment Setup

### Prerequisites

- Linux system (memory monitoring uses `/proc/meminfo`)
- Conda (Miniconda or Mambaforge)
- At least 8 GB RAM (16 GB recommended for large FCS datasets)

### Automated Setup

```bash
bash setup_env.sh
```

### Manual Setup

```bash
# Create conda environment with R and all packages
conda env create -f environment.yml
conda activate cytometry_env

# Install additional R packages not available via conda
Rscript install_extra_packages.R

# Install Python dependencies for Leiden clustering
pip install leidenalg igraph numpy
```

### Verify Installation

```bash
conda activate cytometry_env
Rscript -e "library(flowCore); library(data.table); library(lme4); library(uwot); library(leiden); cat('All OK\n')"
```

---

## 14. Project Structure

```
cytometry/
├── README.md                         # This file
├── CLAUDE.md                         # Project instructions for AI assistants
├── environment.yml                   # Conda environment specification
├── install_extra_packages.R          # Post-conda package installer
├── setup_env.sh                      # Automated setup script
│
├── src/                              # Source code (9 subprojects)
│   ├── R_shared/                     #   8 shared utility modules
│   │   ├── io_utils.R
│   │   ├── gui_utils.R
│   │   ├── memory_utils.R
│   │   ├── cluster_utils.R
│   │   ├── plot_utils.R
│   │   ├── flowframe_utils.R
│   │   ├── fcs_utils.R
│   │   └── stats_utils.R
│   ├── loading/                      #   6 data loaders
│   │   ├── load_fcs_raw.R
│   │   ├── load_csv_expression.R
│   │   ├── load_csv_gating.R
│   │   ├── load_flowjo_wsp.R
│   │   ├── load_manual_gating.R
│   │   └── create_flowframes.R
│   ├── preprocessing/                #   5 preprocessing scripts
│   │   ├── quality_control.R
│   │   ├── transformation.R
│   │   ├── downsampling.R
│   │   ├── normalization.R
│   │   └── preprocess_pipeline.R
│   ├── harmonization/                #   3 harmonization scripts
│   │   ├── panel_comparison.R
│   │   ├── panel_homogenization.R
│   │   └── batch_correction.R
│   ├── clustering/                   #   3 clustering scripts
│   │   ├── leiden.R
│   │   ├── flowjo_import.R
│   │   └── cluster_assignment.R
│   ├── analysis/                     #   7 analysis scripts
│   │   ├── cluster_splitting.R
│   │   ├── differential_abundance.R
│   │   ├── cluster_enrichment.R
│   │   ├── marker_statistics.R
│   │   ├── expression_matrices.R
│   │   ├── abundance_plotting.R
│   │   └── cell_ratio.R
│   ├── visualization/                #   7 visualization scripts
│   │   ├── umap.R
│   │   ├── heatmaps.R
│   │   ├── density_plots.R
│   │   ├── barplots.R
│   │   ├── upset_plots.R
│   │   ├── venn_diagrams.R
│   │   └── pca.R
│   └── metadata_qc/                 #   4 metadata QC scripts
│       ├── check_metadata.R
│       ├── check_acquisition.R
│       ├── overlap_statistics.R
│       └── compare_gatings.R
│
├── pipelines/                        # Pipeline orchestrators
│   ├── pipeline_A_raw_fcs.R
│   ├── pipeline_B_csv_expression.R
│   ├── pipeline_C_csv_gating.R
│   ├── pipeline_D_flowjo_comparison.R
│   └── pipeline_E_preprocessed_umap.R
│
├── config/                           # Configuration files
│   ├── pipeline_A_config.R           #   Per-pipeline configs
│   ├── pipeline_B_config.R
│   ├── pipeline_C_config.R
│   ├── pipeline_D_config.R
│   ├── pipeline_E_config.R
│   ├── Parameters_template.R         #   Master template (~80 params)
│   ├── standalone/                   #   Per-script standalone configs
│   └── study_params/                 #   Your study-specific configs
│
├── docs/                             # Documentation
│   ├── getting-started.md
│   ├── r-shared-reference.md
│   ├── pipelines/                    #   Per-pipeline guides
│   └── subprojects/                  #   Per-module references
│
├── test_data/                        # Test fixtures
│   ├── generate_test_data.R          #   Generate synthetic FCS + CSV
│   ├── fcs_samples/                  #   Synthetic FCS files + conditions
│   ├── csv_expression/               #   Per-sample expression CSVs
│   ├── csv_gating/                   #   Gating CSVs + FCS pairs
│   └── examples/                     #   Format reference files
│
└── output/                           # Pipeline output (gitignored)
    ├── pipeline_A_results/
    ├── pipeline_B_results/
    ├── pipeline_C_results/
    └── pipeline_E_results/
```

---

## 15. Testing

### Generate Test Data

```bash
Rscript test_data/generate_test_data.R
```

Creates:
- `test_data/fcs_samples/`: 3 synthetic FCS files (50 cells x 6 markers) + conditions.csv
- `test_data/csv_expression/`: 6 expression CSVs (3 samples x 2 populations)
- `test_data/csv_gating/`: 3 gating CSVs + 3 matching FCS files

### Run All Pipelines with Test Data

```bash
cd pipelines

# Pipeline A: Raw FCS → Full analysis
Rscript pipeline_A_raw_fcs.R --config ../config/pipeline_A_config.R

# Pipeline B: CSV Expression → Statistics
Rscript pipeline_B_csv_expression.R --config ../config/pipeline_B_config.R

# Pipeline C: CSV Gating → Cluster analysis
Rscript pipeline_C_csv_gating.R --config ../config/pipeline_C_config.R

# Pipeline E: Preprocessed FCS → UMAP
Rscript pipeline_E_preprocessed_umap.R --config ../config/pipeline_E_config.R
```

### Unit Tests

Each subproject has a `tests/` directory using testthat:

```bash
# Test a single subproject
cd src/loading && Rscript -e "testthat::test_dir('tests/')"

# Test all subprojects
for dir in src/R_shared src/loading src/preprocessing src/harmonization \
           src/clustering src/analysis src/visualization src/metadata_qc; do
  echo "=== Testing $dir ==="
  (cd $dir && Rscript -e "testthat::test_dir('tests/')" 2>&1) || true
done
```

---

## 16. Key Packages

| Package | Purpose |
|---------|---------|
| [flowCore](https://bioconductor.org/packages/flowCore/) | FCS file I/O, flowFrame/flowSet objects |
| [flowWorkspace](https://bioconductor.org/packages/flowWorkspace/) | FlowJo workspace parsing, gating sets |
| [data.table](https://rdatatable.gitlab.io/data.table/) | Fast CSV I/O (semicolon separator default) |
| [lme4](https://github.com/lme4/lme4) | Mixed-effects models (GLMM for differential abundance) |
| [uwot](https://github.com/jlmelville/uwot) | UMAP dimensionality reduction |
| [leiden](https://github.com/TomKellyGenetics/leiden) | Leiden community detection clustering (via reticulate → Python leidenalg) |
| [ggplot2](https://ggplot2.tidyverse.org/) | Core plotting engine |
| [pheatmap](https://github.com/raivokolde/pheatmap) | Clustered heatmaps |
| [patchwork](https://patchwork.data-imaginist.com/) | Combine multiple ggplots |
| [viridis](https://sjmgarnier.github.io/viridis/) | Color-blind-friendly palettes |
| [flowAI](https://bioconductor.org/packages/flowAI/) / [PeacoQC](https://bioconductor.org/packages/PeacoQC/) | Automated quality control |

---

## 17. Troubleshooting

### Memory Issues

Large FCS/CSV files can exhaust RAM. Mitigations:
- Use `downsample_method = "Consistent"` with a moderate `downsample_n` (e.g. 10000)
- Set `max_memory_gb` in config (default 4)
- The framework auto-monitors `/proc/meminfo` and triggers `gc()` at 85% utilization
- Use `read_csv_chunked()` from `memory_utils.R` for large CSVs

### Panel Mismatch

Different FCS files may have different marker sets. The `homogenize_flowframes()` step resolves this by subsetting to a common panel. Use `markername_replacing` to handle naming inconsistencies:

```r
markername_replacing = list(
  `PD1 (CD279)` = "PD1",
  `Fas (CD95)`  = "CD95"
)
```

### tcltk Not Available

GUI dialogs for marker selection require tcltk (X11 display). On headless servers, scripts automatically fall back to CLI prompts. To suppress warnings, ensure `Sys.getenv("DISPLAY")` is empty.

### Leiden Python Dependencies

The R `leiden` package uses `reticulate` to call Python's `leidenalg`. Install the Python dependencies:

```bash
pip install leidenalg igraph numpy
```

If reticulate uses the wrong Python, set:

```bash
export RETICULATE_PYTHON=/path/to/conda/env/bin/python
```

### FCS $TOT Mismatch

Some FCS files report an incorrect `$TOT` keyword (total events). Validate with:

```r
source("src/R_shared/fcs_utils.R")
check_tot("/path/to/file.fcs")
```

---

## Design Conventions

- **Shared utilities first**: All scripts source `src/R_shared/` modules — never duplicate functions across subprojects
- **Function arguments, not globals**: Parameters are passed explicitly, not read from a global CONFIG object
- **GUI with CLI fallback**: Interactive dialogs (tcltk) automatically fall back to `readline()` in headless environments
- **data.table for I/O**: CSV read/write via `fread`/`fwrite` with auto-detected separator (semicolon default)
- **Memory-safe**: Large files use chunked reading; memory checks trigger GC at 85% utilization
- **Reproducible**: All stochastic operations accept a `seed` parameter
- **Multi-format output**: Plots saved as PDF/PNG/SVG via `save_plot_formats()` based on config toggles
