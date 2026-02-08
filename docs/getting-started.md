# Getting Started

This guide walks you through setting up the environment, understanding the project layout, and running your first pipeline.

## Prerequisites

- Linux system (memory monitoring uses `/proc/meminfo`)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda installed
- At least 8 GB RAM (16 GB recommended for large FCS datasets)

## 1. Environment Setup

### Option A: Automated (recommended)

```bash
cd /path/to/cytometry
bash setup_env.sh
```

This will:
1. Check for conda
2. Create (or update) the `cytometry_env` conda environment
3. Install all R packages (CRAN, Bioconductor, GitHub)

### Option B: Manual

```bash
# Create conda environment
conda env create -f environment.yml
conda activate cytometry_env

# Install additional packages not in conda
Rscript install_extra_packages.R
```

### Verify Installation

```r
# Quick smoke test
Rscript -e "library(flowCore); library(data.table); library(lme4); library(uwot); cat('All OK\n')"
```

## 2. Project Structure Overview

The project is organized into **9 subprojects** + supporting directories:

```
cytometry/
├── src/               Source code (9 subprojects)
│   ├── R_shared/          8 utility modules shared by all scripts
│   ├── loading/           Data loaders (FCS, CSV, FlowJo WSP)
│   ├── preprocessing/     QC, transformation, downsampling, normalization
│   ├── harmonization/     Panel comparison & unification, batch correction
│   ├── clustering/        Leiden clustering, FlowJo import, assignment
│   ├── analysis/          GLMM, enrichment, marker statistics
│   ├── visualization/     UMAP, heatmaps, density, barplots, Venn/Euler
│   └── metadata_qc/      Metadata checks, acquisition QC, gating comparison
├── pipelines/         5 orchestrator scripts (A through E)
├── config/            Parameter templates and study-specific configs
└── test_data/         Example data for testing
```

Every subproject script sources shared utilities from `src/R_shared/` using this pattern:

```r
shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "gui_utils.R"))
# ... source only what you need
```

## 3. Configuration

All pipelines are driven by a **params list** defined in a configuration file. Start from the template:

```bash
cp config/Parameters_template.R config/study_params/my_study.R
```

Edit `my_study.R` to set your paths, markers, and analysis options. Key parameters:

```r
params <- list()

# --- Paths ---
params$input_dir     <- "/path/to/your/data/"
params$project_name  <- "My_Study"
params$output_dir    <- "/path/to/output/"

# --- Markers ---
params$markers_of_interest <- c("CD4", "CD8", "CD19", "CD56", "HLA-DR")

# --- Preprocessing ---
params$qc_type        <- "PeacoQC"          # "PeacoQC", "flowAI", or "None"
params$transform_type <- "arcsinh"           # "arcsinh", "logicle", "biexponential", "None"
params$pp_order       <- "QT"                # "QT" (QC first) or "TQ" (Transform first)

# --- Downsampling ---
params$target_ncell_subsampling <- 10000
params$subsampling_method       <- "Maximal" # "Maximal" or "Consistent"

# --- Analysis ---
params$condition_levels <- c("Control", "Treatment")
params$stat_method      <- "mixed_effects"   # "mixed_effects" or "pooled"
params$p_adjust_method  <- "BH"              # Benjamini-Hochberg

# --- Reproducibility ---
params$seed <- 42
```

See [config/Parameters_template.R](../config/Parameters_template.R) for the full list of ~80 configurable parameters.

## 4. Running Your First Pipeline

### Choose a pipeline

| Your data | Pipeline | Command |
|-----------|----------|---------|
| Raw `.fcs` files | A | `Rscript pipelines/pipeline_A_raw_fcs.R --config config/study_params/my_study.R` |
| Expression CSVs | B | `Rscript pipelines/pipeline_B_csv_expression.R --config config/study_params/my_study.R` |
| CSVs with gating columns | C | `Rscript pipelines/pipeline_C_csv_gating.R --config config/study_params/my_study.R` |
| FlowJo `.wsp` files | D | `Rscript pipelines/pipeline_D_flowjo_comparison.R --config config/study_params/my_compare.R` |
| Preprocessed `.fcs` | E | `Rscript pipelines/pipeline_E_preprocessed_umap.R --config config/study_params/my_study.R` |

### Example: Pipeline C with CSV gating data

```bash
cd /path/to/cytometry

# 1. Activate environment
conda activate cytometry_env

# 2. Run pipeline
Rscript pipelines/pipeline_C_csv_gating.R --config config/study_params/my_study.R
```

The pipeline will:
1. Load CSV files from `params$csv_folder`
2. Detect gating columns (binary 0/1 membership)
3. Classify columns (scatter, markers, gating, etc.)
4. Present marker selection GUI (or CLI prompt)
5. Split cells by cluster
6. Run differential abundance (GLMM) or enrichment (Fisher's exact)
7. Generate UMAP, barplots, and heatmaps

### Example: Running individual scripts interactively

You can also run scripts one at a time from an R session:

```r
# Set working directory to the subproject
setwd("/path/to/cytometry/src/loading")

# Source the script (this loads functions but doesn't execute)
source("load_csv_gating.R")

# Call the main function
result <- load_csv_gating(
  input_fcs_path = "/path/to/fcs/",
  input_csv_path = "/path/to/csvs/",
  output_path    = "/path/to/output/",
  processing_mode = "automatic"
)
```

## 5. Input File Formats

### FCS Files
Standard FCS 3.0/3.1 files from flow cytometers. The pipeline reads marker names from FCS `$PnS` keywords.

### CSV Expression Files
Semicolon-separated (`;`) by default, auto-detected. Columns are markers, rows are cells:

```
CD4;CD8;CD19;CD56;HLA-DR
1.23;0.45;2.34;0.12;3.45
0.87;1.56;0.23;2.67;1.89
```

### CSV Gating Files (Cytolution export)
Same format but with additional binary gating columns:

```
index;CD4 | BV421.A;CD8 | APC.A;CD4+ T cells (Clustering annotated);CD8+ T cells (Clustering annotated)
1;1.23;0.45;1;0
2;0.87;1.56;0;1
3;2.10;0.33;1;0
```

Column types are auto-classified:
- **Index**: `index`, `Index`
- **Scatter/Time**: `SSC.H | SSC.H`, `FSC.A | FSC.A`, `TIME | Time`
- **Expression markers**: `CD4 | BV421.A` (MarkerName | ChannelName)
- **Gating**: columns ending in `(Clustering annotated)`, `(Manual gating)`, or `Clustering raw`
- **Dim. reduction**: `cyto_2d_Samples_All_X`, `cyto_2d_Samples_All_Y`
- **Autofluorescence**: `AF.P1 | AF.P1.A` or where marker name equals channel name

### FlowJo Workspaces
`.wsp` files with associated FCS files in the same (or specified) directory. The pipeline reads gating hierarchies and extracts population assignments.

### Conditions CSV
Two-column CSV mapping sample filenames to experimental conditions:

```
file_name;condition
sample_001.csv;Control
sample_002.csv;Treatment
sample_003.csv;Control
```

## 6. Output Structure

Each pipeline creates numbered subdirectories under `params$output_dir`:

```
output/
├── 01_preprocessing/     QC'd and transformed FCS/CSV files
├── 02_harmonization/     Panel-unified files
├── 03_downsampled/       Downsampled data
├── 04_clustering/        Leiden cluster assignments + UMAP coords
├── 05_analysis/          GLMM results, enrichment tables, volcano plots
└── 06_visualization/     UMAP plots, heatmaps, barplots, density plots
```

(Exact subdirectories vary by pipeline. See each pipeline's documentation for details.)

## 7. GUI vs Headless Mode

Many scripts present interactive dialogs for marker/column selection using tcltk. If running on a headless server (no display), the scripts automatically fall back to CLI prompts via `readline()`.

To force CLI mode, ensure `Sys.getenv("DISPLAY")` is empty, or the tcltk package is not installed.

## 8. Memory Management

For large datasets (>1 GB CSVs, >50k cells per file):

- The `memory_utils.R` module monitors `/proc/meminfo` and triggers `gc()` when memory usage exceeds 85%
- Use `read_csv_chunked()` for memory-efficient loading of large CSVs
- Set `params$chunk_size` (default: 50000) and `params$max_cells_per_file` (default: 500000) in your config
- Set `params$max_memory_gb` (default: 4) to cap memory usage

## 9. Running Tests

Each subproject has a `tests/` directory using `testthat`:

```bash
cd /path/to/cytometry/src/loading
Rscript -e "testthat::test_dir('tests/')"
```

Or test all subprojects:

```bash
cd /path/to/cytometry
for dir in src/R_shared src/loading src/preprocessing src/harmonization src/clustering src/analysis src/visualization src/metadata_qc; do
  echo "=== Testing $dir ==="
  (cd $dir && Rscript -e "testthat::test_dir('tests/')" 2>&1) || true
done
```

## Next Steps

- Read the [pipeline guide](pipelines/) for your data type
- Browse the [R_shared API reference](r-shared-reference.md) for available utility functions
- Check the [example data](../test_data/examples/) for input format samples
