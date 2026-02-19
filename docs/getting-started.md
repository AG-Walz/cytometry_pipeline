[< Back to Documentation Index](README.md)

# Getting Started

This guide walks you through everything you need to run the cytometry pipeline: installing the environment, understanding the project layout, writing your first configuration file, running a pipeline end-to-end, and interpreting the results. By the end, you will have a working setup and a clear mental model of how the pieces fit together.

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Environment Setup](#2-environment-setup)
3. [Project Structure Overview](#3-project-structure-overview)
4. [How Shared Utilities Work](#4-how-shared-utilities-work)
5. [Configuration Files](#5-configuration-files)
6. [Running Your First Pipeline](#6-running-your-first-pipeline)
7. [Input File Formats](#7-input-file-formats)
8. [Output Structure](#8-output-structure)
9. [Troubleshooting](#9-troubleshooting)
10. [GUI vs Headless Mode](#10-gui-vs-headless-mode)
11. [Memory Management](#11-memory-management)
12. [Running Tests](#12-running-tests)
13. [Next Steps](#13-next-steps)

---

## 1. Prerequisites

Before you begin, make sure you have:

- **Linux system** — Memory monitoring uses `/proc/meminfo` (Ubuntu, Debian, CentOS, etc.)
- **Miniconda or Anaconda** — [Install Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't have it
- **At least 8 GB RAM** — 16 GB recommended for datasets with >50,000 cells per sample
- **Python 3.8+** — Required for Leiden clustering (the R `leiden` package calls Python's `leidenalg` via `reticulate`)

---

## 2. Environment Setup

### Option A: Automated Setup (Recommended)

The project includes a setup script that handles everything:

```bash
cd /path/to/cytometry
bash setup_env.sh
```

This script will:

1. Check that conda is installed
2. Create (or update) the `cytometry_env` conda environment from `environment.yml`
3. Install all R packages from CRAN, Bioconductor, and GitHub via `install_extra_packages.R`
4. Install the Python `leidenalg` and `igraph` packages needed by the Leiden clustering module

### Option B: Manual Setup

If you prefer to set things up step by step:

```bash
# Step 1: Create the conda environment
conda env create -f environment.yml

# Step 2: Activate it
conda activate cytometry_env

# Step 3: Install additional R packages not available in conda
Rscript install_extra_packages.R
```

### Verify Your Installation

Run this quick smoke test to confirm all critical packages are available:

```bash
conda activate cytometry_env

Rscript -e "
  library(flowCore)
  library(data.table)
  library(lme4)
  library(uwot)
  library(ggplot2)
  library(pheatmap)
  cat('All core packages loaded successfully.\n')
"
```

Expected output:

```
All core packages loaded successfully.
```

If any package fails to load, see [Troubleshooting](#9-troubleshooting) below.

### Generate Test Data

The project includes a test data generator that creates synthetic FCS files and CSV datasets for testing every pipeline:

```bash
cd /path/to/cytometry
Rscript test_data/generate_test_data.R
```

This creates:

```
test_data/
├── fcs_samples/                    # For Pipelines A and E
│   ├── sample_001.fcs              #   50 cells × 6 markers
│   ├── sample_002.fcs              #   50 cells × 6 markers
│   ├── sample_003.fcs              #   50 cells × 6 markers
│   └── conditions.csv              #   Sample-to-condition mapping
├── csv_expression/                 # For Pipeline B
│   ├── sample_001_CD4_T_cells.csv  #   Expression values only
│   ├── sample_001_CD8_T_cells.csv
│   ├── sample_002_CD4_T_cells.csv
│   ├── sample_002_CD8_T_cells.csv
│   ├── sample_003_CD4_T_cells.csv
│   └── sample_003_CD8_T_cells.csv
└── csv_gating/                     # For Pipeline C
    ├── sample_001.csv              #   Expression + binary gating columns
    ├── sample_001.fcs              #   Matching FCS file
    ├── sample_002.csv
    ├── sample_002.fcs
    ├── sample_003.csv
    └── sample_003.fcs
```

---

## 3. Project Structure Overview

The project is organized into **9 subprojects** under `src/`, with pipeline orchestrators in `pipelines/` and configuration templates in `config/`:

```
cytometry/
├── src/                      Source code
│   ├── R_shared/                 Shared utility library (8 modules, 50+ functions)
│   │   ├── io_utils.R                CSV I/O, column classification, path helpers
│   │   ├── gui_utils.R               tcltk dialogs with CLI fallback
│   │   ├── memory_utils.R            RAM monitoring, chunked CSV reading
│   │   ├── cluster_utils.R           Cluster name parsing, label translation
│   │   ├── plot_utils.R              Multi-format plot saving, barplot helpers
│   │   ├── flowframe_utils.R         FlowFrame manipulation, expression extraction
│   │   ├── fcs_utils.R               FCS parsing, CSV-FCS mapping, batch loading
│   │   └── stats_utils.R             Enrichment, GLMM, PCA, UMAP, summary stats
│   ├── loading/                  Data loaders
│   │   ├── load_fcs_raw.R            Raw FCS file loading with truncation
│   │   ├── load_fcs_preprocessed.R   Preprocessed FCS loading
│   │   ├── load_csv_expression.R     Expression CSV loading
│   │   ├── load_csv_gating.R         Gating CSV loading with column classification
│   │   ├── load_flowjo_wsp.R         FlowJo workspace parsing
│   │   └── load_manual_gating.R      Manual gating label import
│   ├── preprocessing/            QC and transformation
│   │   ├── preprocessing.R           Orchestrator (dispatches to QC + transform)
│   │   ├── qc.R                      Quality control (PeacoQC, flowAI)
│   │   ├── transformation.R          arcsinh, logicle, biexponential transforms
│   │   ├── downsampling.R            Maximal/consistent downsampling
│   │   └── normalization.R           Quantile normalization
│   ├── harmonization/            Panel unification and batch correction
│   │   ├── panel_comparison.R        Compare marker panels across samples
│   │   ├── panel_homogenization.R    Unify to common marker set
│   │   └── batch_correction.R        Harmony batch correction
│   ├── clustering/               Cell clustering
│   │   ├── leiden.R                  Leiden community detection + UMAP
│   │   ├── flowjo_import.R           FlowJo gate import
│   │   └── cluster_assignment.R      Unified cluster assignment interface
│   ├── analysis/                 Statistical analysis
│   │   ├── cluster_splitting.R       Split data by cluster
│   │   ├── differential_abundance.R  NB-GLMM testing
│   │   ├── cluster_enrichment.R      Fisher's exact / mixed-effects
│   │   ├── marker_statistics.R       Per-marker lmer + ICC
│   │   ├── abundance_plotting.R      Abundance bar/heatmap plots
│   │   ├── expression_matrices.R     Grouped expression matrices
│   │   └── cell_ratio.R             Population ratio calculation
│   ├── visualization/            Plotting
│   │   ├── umap.R                    UMAP scatter plots
│   │   ├── heatmaps.R               Clustered heatmaps (pheatmap)
│   │   ├── density_plots.R           Density and ridge overlays
│   │   ├── barplots.R                Composition barplots
│   │   ├── pca.R                     PCA scatter + loadings
│   │   ├── upset_plots.R             UpSet intersection plots
│   │   └── venn_diagrams.R           Venn and Euler diagrams
│   └── metadata_qc/             Quality assurance
│       ├── check_metadata.R          Metadata validation barplots
│       ├── check_acquisition.R       FCS acquisition date extraction
│       ├── overlap_statistics.R      Fleiss' kappa agreement
│       └── compare_gatings.R         Cross-workspace gating comparison
├── pipelines/                Config-driven orchestrators
│   ├── pipeline_A_raw_fcs.R          Raw FCS → full analysis
│   ├── pipeline_B_csv_expression.R   CSV expression → marker stats
│   ├── pipeline_C_csv_gating.R       CSV gating → cluster analysis
│   ├── pipeline_D_flowjo_comparison.R FlowJo → gating comparison
│   └── pipeline_E_preprocessed_umap.R Preprocessed FCS → UMAP
├── config/                   Configuration templates
│   ├── pipeline_A_config.R
│   ├── pipeline_B_config.R
│   ├── pipeline_C_config.R
│   ├── pipeline_D_config.R
│   └── pipeline_E_config.R
└── test_data/                Synthetic test data
```

Each subproject directory also contains its own `CLAUDE.md` (project-specific instructions) and `tests/` directory (testthat tests).

---

## 4. How Shared Utilities Work

Every script in every subproject begins by sourcing the utility modules it needs from `src/R_shared/`. This is the foundational pattern of the entire codebase:

```r
# At the top of every script in src/loading/, src/analysis/, etc.
shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
source(file.path(shared_dir, "memory_utils.R"))
# ... source only the modules you need
```

The `shared_dir` path is computed relative to the script's working directory, so it works regardless of where the project is installed. The convention is:

- If your working directory is `src/loading/`, then `dirname(getwd())` gives `src/`, and `file.path(dirname(getwd()), "R_shared")` gives `src/R_shared/`.

The 8 shared modules are:

| Module | Key Functions | Purpose |
|--------|---------------|---------|
| `io_utils.R` | `detect_csv_sep()`, `classify_csv_columns()`, `make_outdir()`, `die()` | CSV I/O, column classification, file management |
| `gui_utils.R` | `tk_multiselect()`, `select_markers()` | Interactive dialogs with CLI fallback |
| `memory_utils.R` | `check_memory()`, `read_csv_chunked()` | RAM monitoring, safe large-file loading |
| `cluster_utils.R` | `clean_cluster_name()`, `translate_labels()`, `labels_from_binary_df()` | Cluster name parsing, label manipulation |
| `plot_utils.R` | `save_plot_formats()`, `create_safe_filename()` | Multi-format plot saving |
| `flowframe_utils.R` | `explore_data()`, `extract_exprmtx_list()`, `from_mtx_to_ffr()` | FlowFrame manipulation |
| `fcs_utils.R` | `check_tot()`, `load_fcs_folder()`, `create_csv_fcs_mapping()` | FCS parsing and loading |
| `stats_utils.R` | `run_umap()`, `run_pca()`, `fit_nb_glmm()`, `calculate_enrichment()` | Statistical analysis and dimensionality reduction |

For complete function signatures and usage examples, see the [R_shared API Reference](r-shared-reference.md).

---

## 5. Configuration Files

All pipelines are driven by a `params` list defined in a configuration file. The config file is a plain R script that creates the `params` list with all necessary settings. Each pipeline has a template in `config/`.

### Anatomy of a Config File

Here is a minimal config for Pipeline A:

```r
# config/study_params/my_study.R

params <- list()

# ─── Paths ──────────────────────────────────────────
params$fcs_folder     <- "/data/raw_fcs/"          # Where your FCS files are
params$output_dir     <- "/output/my_study/"        # Where results go
params$conditions_csv <- "/data/conditions.csv"     # Sample → condition mapping

# ─── Markers ────────────────────────────────────────
params$markers_of_interest <- c("CD3", "CD4", "CD8", "CD19", "CD56", "HLA-DR")
params$markername_replacing <- list()               # Rename map (empty = no renames)

# ─── Preprocessing ──────────────────────────────────
params$qc_type        <- "None"                     # "PeacoQC", "flowAI", or "None"
params$transform_type <- "arcsinh"                  # "arcsinh", "logicle", "biexponential", "None"
params$pp_order       <- "QT"                       # "QT" = QC first, "TQ" = Transform first
params$truncation     <- TRUE                       # Truncate expression to max range

# ─── Downsampling ───────────────────────────────────
params$target_ncell_subsampling <- 10000            # Target cells per sample
params$subsampling_method       <- "Maximal"        # "Maximal" or "Consistent"

# ─── Clustering ─────────────────────────────────────
params$n_neighbors <- 15                            # KNN graph neighbors
params$resolution  <- 0.8                           # Leiden resolution (higher = more clusters)

# ─── Analysis ───────────────────────────────────────
params$condition_levels <- c("Control", "Treatment")
params$cluster_regex    <- "Clustering annotated"
params$p_adjust_method  <- "BH"                     # Benjamini-Hochberg

# ─── Output ─────────────────────────────────────────
params$output_formats <- list(pdf = FALSE, png = TRUE, svg = FALSE)

# ─── Reproducibility ───────────────────────────────
params$seed <- 42
```

### Creating Your Own Config

The easiest way to get started is to copy a pipeline template and edit it:

```bash
# For Pipeline A:
cp config/pipeline_A_config.R config/study_params/my_study.R

# For Pipeline C:
cp config/pipeline_C_config.R config/study_params/my_gating_study.R
```

Then edit the copied file to point to your actual data paths and set your analysis parameters.

### Key Parameters by Pipeline

| Parameter | Pipeline A | Pipeline B | Pipeline C | Pipeline E |
|-----------|-----------|-----------|-----------|-----------|
| `fcs_folder` | Required | — | — | Required |
| `csv_folder` | — | Required | Required | — |
| `markers_of_interest` | Required | Optional | Optional | Required |
| `conditions_csv` | Required | — | Required | — |
| `qc_type` | Required | — | — | — |
| `transform_type` | Required | — | — | — |
| `n_neighbors` | Required | — | — | Required |
| `resolution` | Required | — | — | Required |

---

## 6. Running Your First Pipeline

### Step 1: Activate the Environment

```bash
conda activate cytometry_env
```

### Step 2: Choose a Pipeline and Config

| Your Data | Pipeline | Command |
|-----------|----------|---------|
| Raw `.fcs` files | A | `Rscript pipelines/pipeline_A_raw_fcs.R --config config/pipeline_A_config.R` |
| Expression CSVs | B | `Rscript pipelines/pipeline_B_csv_expression.R --config config/pipeline_B_config.R` |
| CSVs with gating columns | C | `Rscript pipelines/pipeline_C_csv_gating.R --config config/pipeline_C_config.R` |
| FlowJo `.wsp` files | D | `Rscript pipelines/pipeline_D_flowjo_comparison.R --config config/pipeline_D_config.R` |
| Preprocessed `.fcs` | E | `Rscript pipelines/pipeline_E_preprocessed_umap.R --config config/pipeline_E_config.R` |

All pipeline commands are run from the project root directory.

### Step 3: Run It

Let's run Pipeline C as an example, using the bundled test data:

```bash
cd /path/to/cytometry
Rscript pipelines/pipeline_C_csv_gating.R --config config/pipeline_C_config.R
```

You will see console output like this:

```
══════════════════════════════════════════════════════════
  Pipeline C — CSV Gating Cluster Analysis
══════════════════════════════════════════════════════════

── Step 1: Loading CSV gating data ──
  Input FCS: test_data/csv_gating/
  Input CSV: test_data/csv_gating/
  Gating pattern: (Clustering annotated)|(Manual gating)

  Loading sample_001...
    Columns classified:
      Index:       1 (index)
      Scatter:     3 (FSC.A|FSC.A, SSC.A|SSC.A, TIME|Time)
      Markers:     6 (CD3|BV421.A, CD4|BV510.A, CD8|APC.A, ...)
      Gating:      4 (CD4+ T cells (Clustering annotated), ...)
      DimReduc:    0
      AutoFluor:   0
  Loading sample_002...  done
  Loading sample_003...  done
  Loaded 3 samples (150 total cells)

── Step 2: Splitting by cluster ──
  Cluster regex: Clustering annotated
  Found 4 clusters: CD4+ T cells, CD8+ T cells, B cells, NK cells
  Writing per-cluster files...
    sample_001_CD4+_T_cells.csv (12 cells)
    sample_001_CD8+_T_cells.csv (8 cells)
    ...
  Saved cluster_tracking_matrix.csv

── Step 3: Cluster enrichment analysis ──
  Method: pooled (Fisher's exact)
  Conditions: Control vs Treatment
  Testing 4 clusters...
  Results saved to enrichment_results.csv

── Step 4: Abundance plots ──
  Generated abundance_stacked.png
  Generated abundance_proportional.png
  Generated cluster_proportion_heatmap.png
  Generated total_cells_per_cluster.png

══════════════════════════════════════════════════════════
  Pipeline C complete. Output: output/pipeline_C_results/
══════════════════════════════════════════════════════════
```

### Step 4: Explore the Output

```bash
ls -la output/pipeline_C_results/
```

```
output/pipeline_C_results/
├── 01_loaded/
│   └── (loaded gating data)
├── 02_clusters/
│   ├── cluster_tracking_matrix.csv
│   ├── sample_001_CD4+_T_cells.csv
│   ├── sample_001_CD8+_T_cells.csv
│   └── ...
├── 03_analysis/
│   └── enrichment_results.csv
└── 04_abundance/
    ├── abundance_stacked.png
    ├── abundance_proportional.png
    ├── cluster_proportion_heatmap.png
    └── total_cells_per_cluster.png
```

### Running Individual Scripts Interactively

You do not have to use a pipeline orchestrator. You can call any function directly from an R session:

```r
# Start R from the subproject directory
setwd("/path/to/cytometry/src/loading")

# Source the script (loads functions into the environment)
source("load_csv_gating.R")

# Call the function you need
result <- load_csv_gating(
  input_fcs_path  = "/data/fcs/",
  input_csv_path  = "/data/csvs/",
  output_path     = "/output/loaded/",
  gating_pattern  = "(Clustering annotated)",
  processing_mode = "automatic"
)
```

For function signatures and parameter documentation, see the relevant [subproject reference](README.md#subproject-references).

---

## 7. Input File Formats

### FCS Files (Pipelines A, E)

Standard FCS 3.0/3.1 files from any flow cytometer. The pipeline reads marker names from the `$PnS` keywords in the FCS header.

```
                    ┌──────────────────────────────────────┐
                    │         FCS File Structure           │
                    ├──────────────────────────────────────┤
                    │ HEADER  ($BEGINDATA, $ENDDATA, ...)  │
                    │ TEXT    ($PnN, $PnS, $TOT, $DATE)    │
                    │ DATA    (expression matrix)          │
                    └──────────────────────────────────────┘
```

Key keywords the pipeline uses:

| Keyword | Example | Usage |
|---------|---------|-------|
| `$PnS` | `CD3`, `CD4`, `HLA-DR` | Marker (stain) names — used as column headers |
| `$PnN` | `BV421-A`, `APC-A` | Channel (fluorochrome) names |
| `$TOT` | `50000` | Total events — validated by `check_tot()` |
| `$DATE` | `15-MAR-2024` | Acquisition date — extracted by `check_acquisition.R` |

### CSV Expression Files (Pipeline B)

Pure expression data with markers as columns and cells as rows. No gating or index columns. Auto-detected separator (`;` or `,`):

```
CD3;CD4;CD8;CD19;CD56;HLA-DR
2.46;4.01;0.15;0.23;0.45;1.89
3.89;0.12;4.56;0.08;0.34;0.67
0.98;1.56;2.78;3.45;0.12;2.34
```

**Filename convention for Pipeline B:** Files should follow the pattern `{sample}_{population}.csv` so that `run_marker_statistics()` can parse the sample and population from the filename. For example: `sample_001_CD4_T_cells.csv`.

### CSV Gating Files — Cytolution Export (Pipeline C)

Expression data with additional columns for index, scatter/time, binary gating, dimensionality reduction, and autofluorescence. The `classify_csv_columns()` function auto-detects all six column types:

```
index;FSC.A | FSC.A;SSC.A | SSC.A;CD3 | BV421.A;CD4 | BV510.A;CD8 | APC.A;CD4+ T cells (Clustering annotated);CD8+ T cells (Clustering annotated);cyto_2d_Samples_All_X;cyto_2d_Samples_All_Y
1;45000;32000;2.46;4.01;0.15;1;0;-5.2;3.1
2;52000;28000;3.89;0.12;4.56;0;1;2.8;-1.5
3;38000;41000;0.98;1.56;2.78;0;1;0.3;4.7
```

Column type classification:

| Column Type | Detection Rule | Example |
|-------------|---------------|---------|
| **Index** | Column named `index` or `Index` | `index` |
| **Scatter/Time** | Contains `SSC`, `FSC`, or `TIME` | `FSC.A \| FSC.A` |
| **Expression** | Format `MarkerName \| ChannelName` where names differ | `CD4 \| BV510.A` |
| **Gating** | Ends with `(Clustering annotated)`, `(Manual gating)`, or starts with `Clustering raw` | `CD4+ T cells (Clustering annotated)` |
| **Dim. Reduction** | Starts with `cyto_2d_` | `cyto_2d_Samples_All_X` |
| **Autofluorescence** | Format `Name \| Name.A` where marker matches channel | `AF.P1 \| AF.P1.A` |

### FlowJo Workspaces (Pipeline D)

`.wsp` files with their associated FCS files in the same (or a specified) directory. The pipeline uses `CytoML::open_flowjo_xml()` and `flowWorkspace::parseWorkspace()` to extract gating hierarchies and population assignments.

### Conditions CSV (Pipelines A, C)

A two-column CSV mapping sample filenames to experimental conditions. Used by `run_differential_abundance()` and `run_cluster_enrichment()`:

```
file_name;condition
sample_001;Control
sample_002;Control
sample_003;Treatment
```

**Important:** The `file_name` column should match the filenames in your data directory (without path, with or without extension depending on the script). Place this file **outside** your data directory so the loading function doesn't try to process it as a data file.

### Metadata CSV (Metadata QC)

Clinical or experimental metadata with samples as rows. Used by `check_metadata.R` for quality checking:

```
Sample_ID;Datum;Age;Sex;Diagnosis;Treatment_Group;Response;Batch
ICK-001;15/03/2024;45;M;AML;Chemo;pCR;Batch1
ICK-002;15/03/2024;62;F;AML;Chemo;npCR;Batch1
ICK-003;20/03/2024;?;M;ALL;Immuno;?;Batch2
```

Placeholder values (`?`, `na`, `n/a`, `null`, `-`) are automatically detected and treated as missing.

---

## 8. Output Structure

Each pipeline creates numbered subdirectories under `params$output_dir`. The exact structure varies by pipeline:

### Pipeline A Output

```
output/pipeline_A_results/
├── 01_preprocessing/           # QC'd, transformed FCS files
├── 02_harmonization/           # Panel-unified flowFrames
├── 03_clustering/
│   └── leiden_clusters.csv     # UMAP1, UMAP2, cluster, sample
├── 04_analysis/
│   ├── cluster_data/           # Per-sample CSVs with cluster assignments
│   ├── split_clusters/         # Per-cluster files + tracking matrix
│   └── differential_abundance/ # glmm_results.csv + volcano_plot.png
└── 05_visualization/
    ├── umap_sample.png         # UMAP colored by sample
    ├── umap_CD3.png            # UMAP colored by marker expression
    └── ...
```

### Pipeline B Output

```
output/pipeline_B_results/
├── 01_statistics/
│   └── marker_statistics.csv   # Per-marker summary with ICC
└── 02_visualization/
    └── (density/heatmap plots)
```

### Pipeline C Output

```
output/pipeline_C_results/
├── 01_loaded/
├── 02_clusters/
│   ├── cluster_tracking_matrix.csv
│   └── (per-sample per-cluster CSVs)
├── 03_analysis/
│   └── enrichment_results.csv
└── 04_abundance/
    ├── abundance_stacked.png
    ├── abundance_proportional.png
    ├── cluster_proportion_heatmap.png
    └── total_cells_per_cluster.png
```

### Pipeline E Output

```
output/pipeline_E_results/
├── 01_harmonized/
├── 02_clustering/
│   └── leiden_clusters.csv
└── 03_visualization/
    ├── umap_sample.png
    └── umap_*.png
```

### Key Output Files

These are the most important output files you will use for downstream analysis:

| File | Pipeline | Format | What It Contains |
|------|----------|--------|------------------|
| `leiden_clusters.csv` | A, E | CSV | UMAP1, UMAP2, cluster (integer), sample — one row per cell |
| `cluster_tracking_matrix.csv` | A, C | CSV | Rows = samples, columns = clusters, values = cell counts |
| `glmm_results.csv` | A | CSV | Per-cluster statistical test results: estimate, log2FC, p_adjusted, significant |
| `enrichment_results.csv` | C | CSV | Per-cluster Fisher's exact or mixed-effects results: odds_ratio, p_adjusted, log2_enrichment |
| `marker_statistics.csv` | B | CSV | Per-population per-marker summary: mean, between/within variance, ICC |
| `umap_sample.png` | A, E | PNG | UMAP plot colored by sample identity |
| `volcano_plot.png` | A | PNG | Log2FC vs -log10(p_adjusted) for all clusters |

---

## 9. Troubleshooting

### Package Installation Failures

**Problem:** `install_extra_packages.R` fails on a Bioconductor package.

**Solution:** Make sure BiocManager is installed and up to date:

```r
install.packages("BiocManager")
BiocManager::install(c("flowCore", "flowWorkspace", "CytoML"))
```

**Problem:** The `leiden` R package can't find Python's `leidenalg`.

**Solution:** Install the Python dependencies and set the `RETICULATE_PYTHON` environment variable:

```bash
conda activate cytometry_env
pip install leidenalg igraph numpy

# If reticulate picks the wrong Python:
export RETICULATE_PYTHON=$(which python)
```

### Memory Errors

**Problem:** R crashes with "cannot allocate vector of size X Gb".

**Solution:** Use chunked reading and set memory limits in your config:

```r
params$chunk_size       <- 50000    # Rows per chunk (smaller = less memory)
params$max_cells_per_file <- 200000 # Cap cells per file
params$max_memory_gb    <- 4        # Memory limit for chunked reading
```

For very large datasets, also consider downsampling more aggressively:

```r
params$target_ncell_subsampling <- 5000  # Fewer cells per sample
```

### tcltk / GUI Errors

**Problem:** `Error in structure(.External(.C_dotTclObjv, ...): [tcl] invalid command name "toplevel"`

**Solution:** You're on a headless server without a display. The scripts fall back to CLI mode automatically in most cases. If not, make sure `DISPLAY` is unset:

```bash
unset DISPLAY
```

### FCS `$TOT` Mismatch

**Problem:** Warning about `$TOT` not matching event count.

**Solution:** This means the FCS header disagrees with the actual data matrix. Use `check_tot()` to diagnose:

```r
source("src/R_shared/io_utils.R")
source("src/R_shared/fcs_utils.R")
ff <- flowCore::read.FCS("data/problem_file.fcs")
check_tot(ff, index = "problem_file")
```

The file may be corrupted or truncated. Try re-exporting from the cytometer software.

### Pipeline Exits Immediately

**Problem:** The pipeline exits after printing `"Config file not found"`.

**Solution:** Provide the config path relative to your current working directory, or use an absolute path:

```bash
# From project root:
Rscript pipelines/pipeline_A_raw_fcs.R --config config/pipeline_A_config.R

# Or with absolute path:
Rscript pipelines/pipeline_A_raw_fcs.R --config /path/to/cytometry/config/pipeline_A_config.R
```

---

## 10. GUI vs Headless Mode

Many scripts present interactive dialogs for marker or column selection using tcltk. This is convenient when running interactively on a desktop, but on headless servers (no X11 display), the scripts fall back to CLI prompts using `readline()`.

### How It Works

```r
# Inside gui_utils.R:
is_gui_available <- function() {
  requireNamespace("tcltk", quietly = TRUE) &&
    nchar(Sys.getenv("DISPLAY")) > 0
}
```

If `is_gui_available()` returns `FALSE`, all dialog functions automatically switch to numbered text menus:

```
Select markers (enter numbers separated by commas, or 'all'):
  1) CD3
  2) CD4
  3) CD8
  4) CD19
  5) CD56
  6) HLA-DR
> 1,2,3
```

### Bypassing Selection Entirely

In pipeline configs, you can set `marker_selection_mode = "all"` to skip the selection dialog and use all available markers:

```r
params$marker_selection_mode <- "all"
```

---

## 11. Memory Management

The pipeline is designed to handle large datasets (hundreds of thousands of cells across many samples) on machines with limited RAM.

### Automatic Memory Monitoring

The `memory_utils.R` module monitors `/proc/meminfo` and triggers garbage collection when memory usage exceeds a threshold:

```r
# Called internally by long-running functions
check_memory(warn_threshold = 0.85)
```

When triggered, you'll see:

```
Warning: Memory usage at 87%. Running garbage collection...
  Freed 1.2 GB. Available: 8.3 GB
```

### Chunked CSV Reading

For very large CSV files (>1 GB), use `read_csv_chunked()` instead of `data.table::fread()`:

```r
source(file.path(shared_dir, "memory_utils.R"))

big_data <- read_csv_chunked(
  file_path    = "data/merged_all_samples.csv",
  chunk_size   = 50000,     # Rows per chunk
  max_memory_gb = 4         # Pause and GC if memory exceeds this
)
```

### Memory Tips

| Dataset Size | Recommendation |
|-------------|----------------|
| < 100K cells total | No special handling needed |
| 100K–1M cells | Set `params$max_cells_per_file = 500000`, downsample to 10K/sample |
| > 1M cells | Use chunked reading, downsample to 5K/sample, increase `max_memory_gb` |

---

## 12. Running Tests

Each subproject has a `tests/` directory using `testthat`. Tests verify that functions produce correct output for known inputs.

### Test a Single Subproject

```bash
cd /path/to/cytometry/src/loading
Rscript -e "testthat::test_dir('tests/')"
```

### Test All Subprojects

```bash
cd /path/to/cytometry
for dir in src/R_shared src/loading src/preprocessing src/harmonization \
           src/clustering src/analysis src/visualization src/metadata_qc; do
  echo "═══ Testing $dir ═══"
  (cd "$dir" && Rscript -e "testthat::test_dir('tests/')" 2>&1) || true
  echo ""
done
```

Expected output for a passing test suite:

```
═══ Testing src/R_shared ═══
✔ | F W  S  OK | Context
✔ |          12 | io_utils
✔ |           8 | cluster_utils
✔ |           5 | memory_utils

── Results ──────────────────────────────────
[ PASS: 25 | FAIL: 0 | WARN: 0 | SKIP: 0 ]
```

---

## 13. Next Steps

Now that your environment is set up and you've run a test pipeline, here's where to go next:

1. **Read the pipeline guide for your data type.** Each guide walks through every step with example input, console output, and generated files:
   - [Pipeline A — Raw FCS](pipelines/pipeline-A-raw-fcs.md)
   - [Pipeline B — CSV Expression](pipelines/pipeline-B-csv-expression.md)
   - [Pipeline C — CSV Gating](pipelines/pipeline-C-csv-gating.md)
   - [Pipeline D — FlowJo Comparison](pipelines/pipeline-D-flowjo-comparison.md)
   - [Pipeline E — Preprocessed UMAP](pipelines/pipeline-E-preprocessed-umap.md)

2. **Create a config for your study.** Copy the relevant template from `config/` to `config/study_params/` and edit it with your paths and parameters.

3. **Browse the R_shared API.** The [R_shared API Reference](r-shared-reference.md) documents every shared utility function. This is especially useful if you're calling functions directly from an R session rather than using a pipeline orchestrator.

4. **Explore the subproject references.** If you want to understand exactly what a specific function does — its parameters, return value, and edge cases — see the relevant [subproject reference](README.md#subproject-references).

5. **Check the example data.** The `test_data/examples/` directory contains synthetic files demonstrating every input format. Use these as templates when preparing your own data.
