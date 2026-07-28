# Cytometry Data Analysis Pipeline

A modular flow cytometry analysis framework in R, supporting raw FCS files, preprocessed FCS, CSV expression matrices, FlowJo workspace output, and CSV gating data. Designed for immune cell phenotyping with Leiden clustering, differential abundance analysis (NB-GLMM), and enrichment statistics.

---

## Installation

### Prerequisites

- Linux system (memory monitoring uses `/proc/meminfo`)
- Conda (Miniconda or Mambaforge)
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

# Run a pipeline (from project root)
Rscript pipelines/pipeline_A_raw_fcs.R --config config/pipeline_A_config.R
```

### UMAP from a FlowJo workspace

To go straight from a FlowJo `.wsp` + raw FCS files to a shared UMAP coloured by
sample, by FlowJo subgate, and by every marker, see the
[UMAP from a FlowJo Population](docs/scripts/umap-from-flowjo.md) guide. After
[setting up the environment](#setup):

```bash
conda activate cytometry_env
# List the gating tree, then pick a population:
Rscript scripts/umap_from_flowjo.R --wsp data/Killing_AML.wsp --fcs_dir data --list_populations
# Build the UMAP for one population:
Rscript scripts/umap_from_flowjo.R --population "Lymphocytes/Single Cells/Viable cells" --output_dir results/umap_viable
```

Prefer a GUI? [`app/`](app/README.md) is a Dash front-end for exactly this
script — pick the paths in a file browser, click a gate in the workspace's own
gating tree, watch the log, and browse the plots it wrote. It runs the same
command line and shows it to you before running it.

```bash
conda activate cytometry_env
python app/app.py    # -> http://127.0.0.1:8050
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

All pipelines follow the same invocation pattern:

```bash
Rscript pipelines/pipeline_X_name.R --config config/pipeline_X_config.R
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [Getting Started](docs/getting-started.md) | Environment setup, first pipeline run, troubleshooting |
| [Configuration Guide](docs/configuration.md) | Writing config files, parameter reference |
| [File Formats](docs/file-formats.md) | Input and output file format specifications |
| [Standalone Scripts](docs/standalone-scripts.md) | Using individual scripts outside pipelines |
| [UMAP from a FlowJo Population](docs/scripts/umap-from-flowjo.md) | Extract a gated population from a `.wsp` + FCS, build a shared UMAP, colour by sample, FlowJo subgate, and every marker (conjoint + per file) |
| [FlowJo → UMAP Dash app](app/README.md) | Browser front-end for that script: file/gate pickers, live log, results gallery |
| [R_shared API Reference](docs/r-shared-reference.md) | Shared utility library (8 modules, 50+ functions) |
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
│   ├── R_shared/                 Shared utility library (8 modules)
│   ├── loading/                  Data loaders (FCS, CSV, FlowJo WSP)
│   ├── preprocessing/            QC, transformation, downsampling
│   ├── harmonization/            Panel unification, batch correction
│   ├── clustering/               Leiden, FlowJo import
│   ├── analysis/                 GLMM, enrichment, marker stats, abundance
│   ├── visualization/            UMAP, heatmaps, density, barplots
│   └── metadata_qc/             Metadata checks, gating comparison
├── pipelines/                Pipeline orchestrators (A through E)
├── config/                   Parameter templates and study-specific configs
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
