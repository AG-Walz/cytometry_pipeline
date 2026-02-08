# Cytometry Data Analysis Pipeline -- Documentation

Multi-pipeline flow cytometry analysis framework for immune cell phenotyping, supporting raw FCS files, preprocessed FCS, CSV expression matrices, FlowJo workspace output, and CSV gating data.

## Quick Navigation

| Document | Description |
|----------|-------------|
| [Getting Started](getting-started.md) | Environment setup, first run, project overview |
| [R_shared API Reference](r-shared-reference.md) | All shared utility functions across 8 modules |

### Pipeline Guides

| Pipeline | Input | Guide |
|----------|-------|-------|
| [A: Raw FCS](pipelines/pipeline-A-raw-fcs.md) | `.fcs` files | Full pipeline: Load, QC, Transform, Harmonize, Cluster, Analyze, Visualize |
| [B: CSV Expression](pipelines/pipeline-B-csv-expression.md) | Expression CSVs | Lightweight: Load CSVs, marker statistics, density/heatmap |
| [C: CSV Gating](pipelines/pipeline-C-csv-gating.md) | CSVs with gates | Load gating labels, split clusters, GLMM/enrichment, UMAP |
| [D: FlowJo Comparison](pipelines/pipeline-D-flowjo-comparison.md) | `.wsp` files | Compare gating strategies across FlowJo workspaces |
| [E: Preprocessed UMAP](pipelines/pipeline-E-preprocessed-umap.md) | Preprocessed `.fcs` | Harmonize, downsample, Leiden cluster, UMAP |

### Subproject Reference

| Module | Guide |
|--------|-------|
| [loading/](subprojects/loading.md) | FCS, CSV, FlowJo workspace loaders |
| [preprocessing/](subprojects/preprocessing.md) | QC, transformation, downsampling, normalization |
| [harmonization/](subprojects/harmonization.md) | Panel comparison, homogenization, batch correction |
| [clustering/](subprojects/clustering.md) | Leiden clustering, FlowJo gate import, unified assignment |
| [analysis/](subprojects/analysis.md) | Differential abundance, enrichment, marker statistics |
| [visualization/](subprojects/visualization.md) | UMAP, heatmaps, density, barplots, UpSet, Venn/Euler, PCA |
| [metadata_qc/](subprojects/metadata-qc.md) | Metadata barplots, acquisition QC, Fleiss' kappa, gating comparison |

### Example Data

| File | Purpose |
|------|---------|
| [`test_data/examples/`](../test_data/examples/) | Synthetic CSV examples for each pipeline input format |
| [`test_data/examples/README.md`](../test_data/examples/README.md) | Description of all example files |

---

## Project Architecture

```
cytometry/
├── docs/                  <-- You are here
│   ├── README.md              Main index (this file)
│   ├── getting-started.md     First-time setup guide
│   ├── r-shared-reference.md  Shared utility API reference
│   ├── pipelines/             Per-pipeline deep dives
│   └── subprojects/           Per-module reference docs
├── R_shared/              Shared utility library (8 modules)
├── loading/               Data loading (FCS, CSV, WSP)
├── preprocessing/         QC, transformation, downsampling
├── harmonization/         Panel unification, batch correction
├── clustering/            Leiden, FlowJo import, assignment
├── analysis/              GLMM, enrichment, marker stats
├── visualization/         Plots (UMAP, heatmap, density, etc.)
├── metadata_qc/           Metadata checks, gating comparison
├── pipelines/             Config-driven orchestrators (A-E)
├── config/                Parameter templates & study configs
├── test_data/             Example data and synthetic samples
├── environment.yml        Conda environment specification
├── install_extra_packages.R  Post-conda package installation
└── setup_env.sh           Automated environment setup script
```

## Pipeline Selection Guide

```
What is your input data?
│
├─ Raw .fcs files (unprocessed)
│  └─> Pipeline A (full analysis)
│
├─ Preprocessed .fcs files (already QC'd/transformed)
│  └─> Pipeline E (skip QC, go straight to clustering)
│
├─ CSV files with expression columns only
│  └─> Pipeline B (marker statistics + visualization)
│
├─ CSV files with binary gating/clustering columns
│  └─> Pipeline C (cluster splitting + differential abundance)
│
└─ FlowJo .wsp workspace files
   └─> Pipeline D (compare gating strategies across experts)
```

## Key Conventions

- **Shared utilities**: All scripts `source()` from `R_shared/` -- never duplicate functions
- **Parameters**: Passed as function arguments via a `params` list, not as global variables
- **GUI with CLI fallback**: Interactive dialogs use tcltk when available, falling back to `readline()` in headless environments
- **CSV format**: `data.table::fread`/`fwrite` with semicolon separator by default
- **FlowFrame objects**: Built on `flowCore`/`flowWorkspace` Bioconductor packages
- **Statistical tests**: Fisher's exact, mixed-effects logistic regression, negative binomial GLMM
- **Testing**: Each subproject has a `tests/` directory using `testthat`
