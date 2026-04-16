# Cytometry Data Analysis Pipeline — Documentation

Welcome to the documentation for the Cytometry Data Analysis Pipeline, a modular R framework for flow cytometry data processing and immune cell phenotyping. This documentation is organized as a reference manual: start with the Getting Started guide, then dive into the pipeline or module that matches your data.

---

## How to Read This Documentation

This documentation is organized in three layers, from high-level to low-level:

1. **Getting Started** — Environment setup, project layout, running your first pipeline
2. **Pipeline Guides** — End-to-end walkthroughs for each of the five pipelines, with example input, console output, and generated files
3. **Module References** — Detailed function-level documentation for every script in every subproject, with signatures, parameter tables, and example code

If you are **new to the project**, start with [Getting Started](getting-started.md), then read the pipeline guide that matches your input data format.

If you are **calling functions directly** from an R session (not using the Nextflow orchestrator), jump to the subproject reference for the module you need.

If you are **looking for a specific utility function** (CSV I/O, memory management, FlowFrame manipulation, etc.), see the [R_shared API Reference](r-shared-reference.md).

---

## Quick Navigation

### Setup & Overview

| Document | What You'll Learn |
|----------|-------------------|
| [Getting Started](getting-started.md) | How to install the environment, understand the project layout, write your first config, and run a pipeline end-to-end |
| [Configuration Guide](configuration.md) | How to create config files, key parameters for each pipeline, standalone script configs |
| [File Formats](file-formats.md) | Input and output format specifications for all pipelines (FCS, CSV, WSP, plots) |
| [Standalone Scripts](standalone-scripts.md) | Using individual module functions outside of pipelines, with signatures and examples |
| [R_shared API Reference](r-shared-reference.md) | Complete function reference for all 12 shared utility modules (io, gui, memory, cluster, plot, flowframe, fcs, stats, clustering_dispatch, dimred_dispatch, batch_correction_dispatch) |
| [Composability Guide](composability-guide.md) | Data flow contracts, format specs, and 8 example cross-pipeline workflows |
| [Methods Guide](methods-guide.md) | Decision flowcharts and 19 deep-dive pages on every statistical, clustering, transformation, and QC method |

### Pipeline Guides

Each pipeline guide is a self-contained walkthrough. It shows you what input the pipeline expects, what it does at each step, what it prints to the console, and what files it generates.

| Pipeline | Input Data | Guide | When to Use |
|----------|------------|-------|-------------|
| **A: Raw FCS** | Unprocessed `.fcs` files | [Pipeline A Guide](pipelines/pipeline-A-raw-fcs.md) | You have raw FCS files straight from the cytometer and want the full analysis: QC, transformation, harmonization, Leiden clustering, differential abundance, and UMAP visualization |
| **B: CSV Expression** | Expression CSVs (markers only) | [Pipeline B Guide](pipelines/pipeline-B-csv-expression.md) | You have per-population CSV files with expression values and want marker statistics (mean, variance, ICC) plus density and heatmap plots |
| **C: CSV Gating** | CSVs with binary gating columns | [Pipeline C Guide](pipelines/pipeline-C-csv-gating.md) | You have Cytolution-exported CSVs with binary gating columns and want cluster splitting, differential abundance (GLMM), enrichment (Fisher), and abundance plots |
| **D: FlowJo Comparison** | `.wsp` workspace files | [Pipeline D Guide](pipelines/pipeline-D-flowjo-comparison.md) | You have multiple FlowJo workspaces analyzing the same FCS files and want to compare gating strategies with Venn/Euler diagrams and Fleiss' kappa |
| **E: Preprocessed UMAP** | Preprocessed `.fcs` files | [Pipeline E Guide](pipelines/pipeline-E-preprocessed-umap.md) | You have already QC'd and transformed FCS files and want to jump straight to harmonization, Leiden clustering, and UMAP |

### Deciding Which Pipeline to Use

| Your Input Data | Pipeline |
|-----------------|----------|
| Raw `.fcs` files (unprocessed, straight from the cytometer) | **Pipeline A** — full analysis from scratch |
| Preprocessed `.fcs` files (already QC'd and transformed) | **Pipeline E** — skip preprocessing, go straight to clustering |
| CSV files with expression columns only (no gating) | **Pipeline B** — marker statistics + visualization |
| CSV files with binary gating/clustering columns (Cytolution export) | **Pipeline C** — cluster splitting + differential abundance |
| FlowJo `.wsp` workspace files (multiple analysts) | **Pipeline D** — inter-rater gating comparison |

### Subproject References

Each subproject reference documents every function in that module: its signature, parameters, return value, and example usage. These are the building blocks that the pipelines call internally.

| Module | Scripts | Reference | What It Does |
|--------|---------|-----------|--------------|
| `src/loading/` | 6 scripts | [Loading Reference](subprojects/loading.md) | Load data from FCS files, CSV files, FlowJo workspaces, and manual gating exports |
| `src/preprocessing/` | 13 scripts | [Preprocessing Reference](subprojects/preprocessing.md) | Compensation, unmixing, beads, doublets, live/dead, QC (PeacoQC, flowAI), transformation (arcsinh, logicle), downsampling, normalization, SuperCellCyto, auto-gating |
| `src/harmonization/` | 5 scripts | [Harmonization Reference](subprojects/harmonization.md) | Compare marker panels across samples, unify to a common marker set, correct batch effects (Harmony, CytoNorm, ComBat) |
| `src/clustering/` | 9 scripts | [Clustering Reference](subprojects/clustering.md) | Leiden, FlowSOM, PhenoGraph, DBSCAN, meta-clustering, validation, annotation, FlowJo gate import, unified cluster assignment interface |
| `src/analysis/` | 13 scripts | [Analysis Reference](subprojects/analysis.md) | Cluster splitting, NB-GLMM differential abundance, Fisher's enrichment, marker statistics (lmer/ICC), differential expression, effect sizes, correlations, diffcyt, survival analysis, abundance plotting, expression matrices, cell ratios |
| `src/visualization/` | 14 scripts | [Visualization Reference](subprojects/visualization.md) | UMAP, t-SNE, diffusion maps, clustered heatmaps, density/ridge overlays, biaxial, violin, radar, contour, composition barplots, PCA, UpSet intersection plots, Venn/Euler diagrams |
| `src/metadata_qc/` | 4 scripts | [Metadata QC Reference](subprojects/metadata-qc.md) | Metadata validation barplots, FCS acquisition date extraction, Fleiss' kappa inter-rater agreement, gating comparison across FlowJo workspaces |

### Example Data

The `test_data/` directory contains synthetic data for testing every pipeline without needing real patient samples:

| Directory | Contents | Used By |
|-----------|----------|---------|
| [`test_data/fcs_samples/`](../test_data/fcs_samples/) | 3 synthetic FCS files (50 cells, 6 markers each) + `conditions.csv` | Pipelines A and E |
| [`test_data/csv_expression/`](../test_data/csv_expression/) | 6 per-sample/per-population expression CSVs | Pipeline B |
| [`test_data/csv_gating/`](../test_data/csv_gating/) | 3 CSV + 3 FCS file pairs with binary gating columns | Pipeline C |
| [`test_data/examples/`](../test_data/examples/) | Standalone example files (expression matrix, Cytolution export, conditions, metadata, tracking matrix, manual gating, config) | Individual scripts |

See [`test_data/examples/README.md`](../test_data/examples/README.md) for detailed descriptions of each example file.

---

## Project Architecture

```
cytometry/
├── README.md                 Main project README (comprehensive overview)
├── docs/                     Documentation (you are here)
│   ├── README.md                 This index page
│   ├── getting-started.md        First-time setup and walkthrough
│   ├── configuration.md          Config creation and parameter reference
│   ├── file-formats.md           Input/output format specifications
│   ├── standalone-scripts.md     Using modules outside of pipelines
│   ├── r-shared-reference.md     Shared utility API (12 modules, 50+ functions)
│   ├── composability-guide.md   Cross-pipeline data flow and example workflows
│   ├── methods-guide.md         Decision flowcharts + 19 method deep-dives
│   ├── methods/                 Per-method deep-dive pages (peer-reviewed refs)
│   ├── pipelines/                Per-pipeline deep-dive guides
│   │   ├── pipeline-A-raw-fcs.md
│   │   ├── pipeline-B-csv-expression.md
│   │   ├── pipeline-C-csv-gating.md
│   │   ├── pipeline-D-flowjo-comparison.md
│   │   └── pipeline-E-preprocessed-umap.md
│   └── subprojects/              Per-module function references
│       ├── loading.md
│       ├── preprocessing.md
│       ├── harmonization.md
│       ├── clustering.md
│       ├── analysis.md
│       ├── visualization.md
│       └── metadata-qc.md
├── src/                      Source code (8 subprojects, 76 scripts)
│   ├── R_shared/                 Shared utility library (12 modules)
│   ├── loading/                  Data loaders (FCS, CSV, FlowJo WSP)
│   ├── preprocessing/            QC, transformation, downsampling, normalization
│   ├── harmonization/            Panel unification, batch correction
│   ├── clustering/               Leiden, FlowJo import, cluster assignment
│   ├── analysis/                 GLMM, enrichment, marker stats, abundance
│   ├── visualization/            UMAP, heatmaps, density, barplots, Venn/Euler
│   └── metadata_qc/             Metadata checks, gating comparison
├── nextflow/                 Nextflow pipeline orchestration
│   ├── main.nf                   Entry point (--pipeline A|B|C|D|E)
│   ├── modules/                  Per-step Nextflow process definitions
│   └── bin/                      CLI wrappers that call src/ R functions
├── config/                   YAML parameter files
│   ├── pipeline_A_params.yaml    Template for Pipeline A
│   ├── pipeline_B_params.yaml    Template for Pipeline B
│   ├── pipeline_C_params.yaml    Template for Pipeline C
│   ├── pipeline_D_params.yaml    Template for Pipeline D
│   ├── pipeline_E_params.yaml    Template for Pipeline E
│   └── study_params/             Your study-specific YAML configs
├── test_data/                Synthetic test data
│   ├── generate_test_data.R      Script to regenerate all test fixtures
│   ├── fcs_samples/              Synthetic FCS files for Pipelines A/E
│   ├── csv_expression/           Expression CSVs for Pipeline B
│   ├── csv_gating/               Gating CSVs for Pipeline C
│   └── examples/                 Standalone example files
├── environment.yml           Conda environment specification
├── install_extra_packages.R  Post-conda R package installation
└── setup_env.sh              Automated environment setup script
```

---

## Key Conventions

These conventions apply across all pipelines and modules:

- **Shared utilities first.** Every script sources functions from `src/R_shared/` before doing anything else. Utilities are never duplicated across scripts.

- **Config-driven parameters.** Pipelines read all settings from a `params` list defined in a config file. No global variables, no hardcoded paths.

- **GUI with CLI fallback.** Interactive marker/column selection dialogs use tcltk when available. On headless servers, they fall back to `readline()` prompts automatically.

- **Semicolon-separated CSVs.** The lab's standard export format uses `;` as a separator. `detect_csv_sep()` auto-detects the separator at load time, so both `;` and `,` work.

- **FlowFrame objects.** Flow cytometry data is represented as `flowCore::flowFrame` objects throughout the pipeline. Expression matrices are extracted with `extract_exprmtx_list()` when needed for statistics or clustering.

- **Multi-format plot output.** All visualization functions save plots via `save_plot_formats()`, which writes PDF, PNG, and/or SVG based on the `output_formats` parameter in your config.

- **testthat testing.** Each subproject has a `tests/` directory. Run `testthat::test_dir('tests/')` from any subproject to verify correctness.

---

## Getting Help

- **First-time setup problems?** See the [Troubleshooting section](getting-started.md#9-troubleshooting) in the Getting Started guide.
- **Function signature questions?** Check the [R_shared API Reference](r-shared-reference.md) or the relevant subproject reference.
- **Pipeline not producing expected output?** Read the pipeline guide for your data type — each guide includes expected console output and file listings.
