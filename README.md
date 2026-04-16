# Cytometry Data Analysis Pipeline

A modular flow cytometry analysis framework in R, supporting raw FCS files, preprocessed FCS, CSV expression matrices, FlowJo workspace output, and CSV gating data. Designed for immune cell phenotyping with Leiden clustering, differential abundance analysis (NB-GLMM), and enrichment statistics.

---

## Installation

### Prerequisites

- Nextflow (>= 23.04.0) and Docker (recommended for reproducible runs)
- Alternatively: Conda (Miniconda or Mambaforge) for local R development
- At least 8 GB RAM (16 GB recommended for large FCS datasets)

### Setup

```bash
git clone https://github.com/AG-Walz/cytometry_pipeline.git
cd cytometry_pipeline

# Option A: Automated
bash setup_env.sh

# Option B: Manual
conda env create -f environment.yml
conda activate cytometry_env
Rscript install_extra_packages.R
pip install leidenalg igraph numpy
```

### Verify Installation

```bash
conda activate cytometry_env
Rscript -e "library(flowCore); library(data.table); library(lme4); library(uwot); library(leiden); cat('All OK\n')"
```

---

## Quick Start

```bash
# Generate synthetic test data
Rscript test_data/generate_test_data.R

# Run a pipeline via Nextflow (from project root)
nextflow run nextflow/main.nf --pipeline A --input_dir /data/fcs \
    --markers_of_interest "CD3,CD4,CD8" --output_dir ./results -profile docker

# Or with a param file:
nextflow run nextflow/main.nf -params-file config/pipeline_A_params.yaml -profile docker
```

---

## Choosing a Pipeline

| Pipeline | Input | What It Does | When to Use |
|----------|-------|--------------|-------------|
| [**A**](docs/pipelines/pipeline-A-raw-fcs.md) | Raw `.fcs` files | QC, transform, harmonize, Leiden clustering, differential abundance, UMAP | Raw FCS from the cytometer, full end-to-end analysis |
| [**B**](docs/pipelines/pipeline-B-csv-expression.md) | Expression CSVs | Marker statistics (lmer, ICC), heatmap | Exported expression matrices, per-marker stats |
| [**C**](docs/pipelines/pipeline-C-csv-gating.md) | CSVs with gating columns | Split clusters, Fisher/GLMM enrichment, abundance plots | Cytolution or manual gating exports with binary gate columns |
| [**D**](docs/pipelines/pipeline-D-flowjo-comparison.md) | FlowJo `.wsp` files | Compare gating strategies, Fleiss' kappa, Venn/UpSet | Comparing multiple analysts' gating of the same samples |
| [**E**](docs/pipelines/pipeline-E-preprocessed-umap.md) | Preprocessed `.fcs` | Harmonize, Leiden clustering, UMAP | Already QC'd and transformed FCS, quick clustering |

All pipelines follow the same Nextflow invocation pattern:

```bash
nextflow run nextflow/main.nf --pipeline X --input_dir /path/to/data --output_dir ./results -profile docker
# Or with a param file:
nextflow run nextflow/main.nf -params-file config/pipeline_X_params.yaml -profile docker
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [Getting Started](docs/getting-started.md) | Environment setup, first pipeline run, troubleshooting |
| [Configuration Guide](docs/configuration.md) | Writing config files, parameter reference |
| [File Formats](docs/file-formats.md) | Input and output file format specifications |
| [Standalone Scripts](docs/standalone-scripts.md) | Using individual scripts outside pipelines |
| [R_shared API Reference](docs/r-shared-reference.md) | Shared utility library (12 modules, 50+ functions) |
| [Composability Guide](docs/composability-guide.md) | Cross-pipeline data flow, format contracts, example workflows |
| [Methods Guide](docs/methods-guide.md) | Decision flowcharts and 19 deep-dive pages with peer-reviewed references |

### Pipeline Guides

- [Pipeline A — Raw FCS](docs/pipelines/pipeline-A-raw-fcs.md)
- [Pipeline B — CSV Expression](docs/pipelines/pipeline-B-csv-expression.md)
- [Pipeline C — CSV Gating](docs/pipelines/pipeline-C-csv-gating.md)
- [Pipeline D — FlowJo Comparison](docs/pipelines/pipeline-D-flowjo-comparison.md)
- [Pipeline E — Preprocessed UMAP](docs/pipelines/pipeline-E-preprocessed-umap.md)

### Module References

- [Loading](docs/subprojects/loading.md) | [Preprocessing](docs/subprojects/preprocessing.md) | [Harmonization](docs/subprojects/harmonization.md) | [Clustering](docs/subprojects/clustering.md)
- [Analysis](docs/subprojects/analysis.md) | [Visualization](docs/subprojects/visualization.md) | [Metadata QC](docs/subprojects/metadata-qc.md)

---

## Project Structure

```
cytometry/
├── src/                      Source code (8 subprojects)
│   ├── R_shared/                 Shared utility library (12 modules)
│   ├── loading/                  Data loaders (FCS, CSV, FlowJo WSP)
│   ├── preprocessing/            QC, transformation, downsampling
│   ├── harmonization/            Panel unification, batch correction
│   ├── clustering/               Leiden, FlowJo import
│   ├── analysis/                 GLMM, enrichment, marker stats, abundance
│   ├── visualization/            UMAP, heatmaps, density, barplots
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
├── docs/                     Documentation
├── test_data/                Synthetic test data
└── output/                   Pipeline output (gitignored)
```

---

## Key Packages

| Package | Purpose |
|---------|---------|
| [flowCore](https://bioconductor.org/packages/flowCore/) | FCS file I/O, flowFrame objects |
| [flowWorkspace](https://bioconductor.org/packages/flowWorkspace/) | FlowJo workspace parsing |
| [data.table](https://rdatatable.gitlab.io/data.table/) | Fast CSV I/O |
| [lme4](https://github.com/lme4/lme4) | Mixed-effects models (GLMM) |
| [uwot](https://github.com/jlmelville/uwot) | UMAP dimensionality reduction |
| [leiden](https://github.com/TomKellyGenetics/leiden) | Leiden community detection |
| [ggplot2](https://ggplot2.tidyverse.org/) | Plotting engine |
