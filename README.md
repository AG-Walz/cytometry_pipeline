# Cytometry Pipeline

A modular flow cytometry analysis framework in R, supporting raw FCS files, preprocessed FCS, CSV expression matrices, FlowJo workspace output, and CSV gating data. Designed for immune cell phenotyping with Leiden clustering, differential abundance analysis (NB-GLMM), and enrichment statistics.

---

## Table of Contents

- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Pipelines](#pipelines)
- [Modules](#modules)
  - [R\_shared — Shared Utility Library](#r_shared--shared-utility-library)
  - [loading](#loading)
  - [preprocessing](#preprocessing)
  - [harmonization](#harmonization)
  - [clustering](#clustering)
  - [analysis](#analysis)
  - [visualization](#visualization)
  - [metadata\_qc](#metadata_qc)
- [Configuration](#configuration)
- [Documentation](#documentation)
- [Test Data](#test-data)
- [Key Packages](#key-packages)

---

## Quick Start

### Prerequisites

- R 4.3+ (RStudio recommended but not required)
- Conda (Miniconda or Mambaforge)

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/AG-Walz/cytometry_pipeline.git
cd cytometry_pipeline

# 2. Create the conda environment
conda env create -f environment.yml
conda activate cytometry

# 3. Install additional R packages not available via conda
Rscript install_extra_packages.R
```

Or use the automated setup script:

```bash
bash setup_env.sh
```

### Running a Pipeline

1. Copy and edit the parameter template for your study:

```r
file.copy("config/Parameters_template.R", "config/study_params/my_study.R")
# Edit my_study.R with your paths, markers, and analysis options
```

2. Source the appropriate pipeline script in R:

```r
source("config/study_params/my_study.R")  # loads `params` list
source("pipelines/pipeline_A_raw_fcs.R")  # or B, C, D, E
```

See the [Getting Started Guide](docs/getting-started.md) for detailed instructions.

---

## Project Structure

```
cytometry_pipeline/
├── R_shared/            Shared utility library (sourced by all subprojects)
├── loading/             Data loading — FCS, CSV, FlowJo WSP
├── preprocessing/       QC, transformation, downsampling, normalization
├── harmonization/       Panel comparison, homogenization, batch correction
├── clustering/          Leiden clustering, FlowJo import, cluster assignment
├── analysis/            Differential abundance, enrichment, marker stats
├── visualization/       UMAP, heatmaps, density, barplots, PCA, Venn/Euler
├── metadata_qc/         Metadata validation, acquisition QC, gating comparison
├── pipelines/           Config-driven pipeline orchestrators (A–E)
├── config/
│   ├── pipeline_A_config.R   Ready-to-edit config for each pipeline
│   ├── pipeline_B_config.R
│   ├── pipeline_C_config.R
│   ├── pipeline_D_config.R
│   ├── pipeline_E_config.R
│   ├── standalone/           Configs for running individual scripts
│   └── study_params/         Your study-specific parameter files
├── docs/                In-depth documentation
└── test_data/           Example input files for testing
```

---

## Pipelines

Each pipeline chains together loading, processing, analysis, and visualization steps. Choose the one that matches your input data.

> Detailed documentation for each pipeline is linked below.

| Pipeline | Input Type | What It Does | Docs |
|----------|-----------|--------------|------|
| **[A — Raw FCS](pipelines/pipeline_A_raw_fcs.R)** | `.fcs` files | Full end-to-end: Load → QC → Transform → Harmonize → Downsample → Cluster → Differential Abundance → Visualize | [Guide](docs/pipelines/pipeline-A-raw-fcs.md) |
| **[B — CSV Expression](pipelines/pipeline_B_csv_expression.R)** | Expression CSVs | Lightweight: Load CSVs → Marker statistics → Density plots / Heatmaps | [Guide](docs/pipelines/pipeline-B-csv-expression.md) |
| **[C — CSV Gating](pipelines/pipeline_C_csv_gating.R)** | CSVs with gating columns | Load gating CSVs → Split clusters → DA / Enrichment → UMAP / Barplots | [Guide](docs/pipelines/pipeline-C-csv-gating.md) |
| **[D — FlowJo Comparison](pipelines/pipeline_D_flowjo_comparison.R)** | `.wsp` workspace files | Compare FlowJo gating strategies → Fleiss' kappa → Venn/Euler/UpSet | [Guide](docs/pipelines/pipeline-D-flowjo-comparison.md) |
| **[E — Preprocessed UMAP](pipelines/pipeline_E_preprocessed_umap.R)** | Preprocessed `.fcs` | Quick clustering: Load (skip QC) → Harmonize → Downsample → Leiden → UMAP | [Guide](docs/pipelines/pipeline-E-preprocessed-umap.md) |

### Pipeline Flow Diagrams

**Pipeline A** (complete):
```
Raw FCS → QC → Transform → Harmonize → Downsample → Leiden → Split Clusters
    → Differential Abundance → Enrichment → UMAP / Heatmaps / Barplots
```

**Pipeline D** (comparison):
```
WSP₁ + WSP₂ + ... → Extract Gates → Compare Populations
    → Fleiss' Kappa → Venn / Euler / UpSet Plots
```

---

## Modules

### R\_shared — Shared Utility Library

All subproject scripts source from `R_shared/` for common functionality. Never duplicate utilities — always add shared functions here.

> Full API reference: [docs/r-shared-reference.md](docs/r-shared-reference.md)

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| [io\_utils.R](R_shared/io_utils.R) | File I/O, CSV handling, column classification, directory management | `detect_csv_sep`, `classify_csv_columns`, `make_outdir`, `die`, `save_and_cleanup` |
| [gui\_utils.R](R_shared/gui_utils.R) | Interactive GUI dialogs (tcltk) with automatic CLI fallback | `select_markers`, `tk_multiselect`, `is_gui_available` |
| [memory\_utils.R](R_shared/memory_utils.R) | Memory monitoring, chunked CSV reading, garbage collection | `check_memory`, `read_csv_chunked`, `get_available_memory_gb` |
| [cluster\_utils.R](R_shared/cluster_utils.R) | Cluster name parsing, label manipulation, binary-to-categorical conversion | `clean_cluster_name`, `labels_from_binary_df`, `convert_binary_to_single_cluster` |
| [plot\_utils.R](R_shared/plot_utils.R) | Multi-format plot saving (PDF/PNG/SVG), barplots, filename sanitization | `save_plot_formats`, `population_barplot`, `create_safe_filename` |
| [flowframe\_utils.R](R_shared/flowframe_utils.R) | FlowFrame/FlowSet manipulation, expression matrix extraction | `explore_data`, `compare_panels`, `extract_exprmtx_list`, `from_mtx_to_ffr` |
| [fcs\_utils.R](R_shared/fcs_utils.R) | Low-level FCS parsing, CSV-FCS mapping, batch loading | `readFCStext`, `check_tot`, `load_fcs_folder`, `create_csv_fcs_mapping` |
| [stats\_utils.R](R_shared/stats_utils.R) | Statistical tests, mixed-effects models, PCA, UMAP | `calculate_enrichment`, `fit_nb_glmm`, `run_pca`, `run_umap`, `calculate_statistics` |

---

### loading

Load cytometry data from various file formats into R objects.

> Detailed reference: [docs/subprojects/loading.md](docs/subprojects/loading.md)

| Script | Purpose | Main Function |
|--------|---------|---------------|
| [load\_fcs\_raw.R](loading/load_fcs_raw.R) | Load raw FCS files into flowFrame objects, optionally extract CSVs | `load_fcs_raw(input_dir, sample_names, ...)` |
| [load\_csv\_expression.R](loading/load_csv_expression.R) | Load CSV expression matrices, auto-detect separator, validate markers | `load_csv_expression(input_dir, sep, ...)` |
| [load\_csv\_gating.R](loading/load_csv_gating.R) | Load CSVs with binary gating columns, classify columns, match to FCS | `load_csv_gating(input_fcs_path, input_csv_path, ...)` |
| [load\_flowjo\_wsp.R](loading/load_flowjo_wsp.R) | Load FlowJo workspaces, extract gating sets and population assignments | `load_flowjo_wsp(wsp_path, fcs_dir, ...)` |
| [load\_manual\_gating.R](loading/load_manual_gating.R) | Load manual gating labels from CSV, supports label translation | `load_manual_gating(csv_dir, fcs_dir, ...)` |
| [create\_flowframes.R](loading/create_flowframes.R) | Create flowFrame objects from expression matrices, merge into flowSet | `create_flowframes_from_matrices(expr_list, ...)` |

---

### preprocessing

Quality control, mathematical transformation, downsampling, and normalization of flow cytometry data.

> Detailed reference: [docs/subprojects/preprocessing.md](docs/subprojects/preprocessing.md)

| Script | Purpose | Main Function |
|--------|---------|---------------|
| [quality\_control.R](preprocessing/quality_control.R) | Automated QC via flowAI or PeacoQC — removes margin events, doublets, anomalies | `quality_control(flowframe, sample_name, qc_type, ...)` |
| [transformation.R](preprocessing/transformation.R) | Arcsinh, logicle, or biexponential transformation of fluorescence channels | `transform_flowframe(flowframe, transform_type, ...)` |
| [downsampling.R](preprocessing/downsampling.R) | Reduce event counts: proportional ("Maximal") or fixed ("Consistent") | `downsample_flowframes(flow_frames, target_cells, method, ...)` |
| [normalization.R](preprocessing/normalization.R) | Z-score scaling, quantile normalization, Harmony batch correction | `z_score_scaling()`, `quantile_norm()`, `harmony_normalize()` |
| [preprocess\_pipeline.R](preprocessing/preprocess_pipeline.R) | Orchestrates QC → Transform → Normalize → Scale with configurable order | `run_preprocessing(flow_frames, params, ...)` |

---

### harmonization

Unify marker panels across heterogeneous FCS files and correct batch effects.

> Detailed reference: [docs/subprojects/harmonization.md](docs/subprojects/harmonization.md)

| Script | Purpose | Main Function |
|--------|---------|---------------|
| [panel\_comparison.R](harmonization/panel_comparison.R) | Compare marker panels across flowFrames — identify shared/divergent markers | `compare_panels_detailed(flow_frames)` |
| [panel\_homogenization.R](harmonization/panel_homogenization.R) | Unify flowFrames to common marker set, fill missing with NA | `homogenize_flowframes(flow_frames, markers_of_interest, ...)` |
| [batch\_correction.R](harmonization/batch_correction.R) | Harmony-based batch effect correction preserving biological signal | `run_batch_correction(flow_frames, batch_labels, ...)` |

---

### clustering

Unsupervised and supervised cell population identification.

> Detailed reference: [docs/subprojects/clustering.md](docs/subprojects/clustering.md)

| Script | Purpose | Main Function |
|--------|---------|---------------|
| [leiden.R](clustering/leiden.R) | Leiden community detection with optional UMAP embedding (KNN graph via uwot) | `run_leiden(expr_matrix, ...)`, `run_leiden_with_umap(...)` |
| [flowjo\_import.R](clustering/flowjo_import.R) | Import gating hierarchies from FlowJo workspaces, extract per-cell assignments | `import_flowjo_gates(wsp_path, fcs_dir, ...)` |
| [cluster\_assignment.R](clustering/cluster_assignment.R) | Unified interface for cluster assignments from any source (Leiden, FlowJo, CSV, FlowSOM) | `get_cluster_assignments(source, ...)` |

---

### analysis

Statistical analysis of clustered cytometry data.

> Detailed reference: [docs/subprojects/analysis.md](docs/subprojects/analysis.md)

| Script | Purpose | Main Function |
|--------|---------|---------------|
| [cluster\_splitting.R](analysis/cluster_splitting.R) | Split combined expression data into per-cluster CSV files + tracking matrix | `split_by_cluster(input_dir, output_dir, ...)` |
| [differential\_abundance.R](analysis/differential_abundance.R) | Test differential abundance between conditions via negative binomial GLMM | `run_differential_abundance(data_dir, conditions_csv, ...)` |
| [cluster\_enrichment.R](analysis/cluster_enrichment.R) | Fisher's exact test or mixed-effects logistic regression for cluster enrichment | `run_cluster_enrichment(data_dir, conditions_csv, ...)` |
| [marker\_statistics.R](analysis/marker_statistics.R) | Per-marker summary statistics (mean, SD, ICC) using linear mixed-effects models | `run_marker_statistics(input_folder, ...)` |
| [expression\_matrices.R](analysis/expression_matrices.R) | Generate grouped expression matrices from per-cluster data | `create_expression_matrices(input_dir, ...)` |
| [abundance\_plotting.R](analysis/abundance_plotting.R) | Stacked, grouped, or proportional barplots from cluster abundance data | `plot_abundance(tracking_matrix_path, ...)` |
| [cell\_ratio.R](analysis/cell_ratio.R) | Calculate ratios between cluster populations per sample (e.g., CD4/CD8) | `calculate_cell_ratios(abundance_csv, ...)` |

---

### visualization

Publication-quality plots for cytometry data.

> Detailed reference: [docs/subprojects/visualization.md](docs/subprojects/visualization.md)

| Script | Purpose | Main Function |
|--------|---------|---------------|
| [umap.R](visualization/umap.R) | UMAP embeddings with per-sample, per-marker, and per-cluster overlays | `plot_umap(flow_frames, markers_of_interest, ...)` |
| [heatmaps.R](visualization/heatmaps.R) | Clustered heatmaps (pheatmap) — expression and median-expression per cluster | `plot_heatmap(...)`, `plot_cluster_heatmap(...)` |
| [density\_plots.R](visualization/density_plots.R) | Per-marker density or ridge plots, grouped by sample or condition | `plot_density(flow_frames, markers, ...)` |
| [barplots.R](visualization/barplots.R) | Composition barplots — stacked, grouped, proportional layouts | `plot_composition_barplot(data, ...)` |
| [upset\_plots.R](visualization/upset_plots.R) | UpSet plots for set intersection visualization | `plot_upset(data, set_columns, ...)` |
| [venn\_diagrams.R](visualization/venn_diagrams.R) | Venn diagrams and area-proportional Euler diagrams | `plot_euler(sets, ...)`, `plot_venn(sets, ...)` |
| [pca.R](visualization/pca.R) | PCA scatter plots and loadings plots | `plot_pca(expr_matrix, ...)` |

---

### metadata\_qc

Metadata validation, acquisition quality checks, and inter-rater agreement.

> Detailed reference: [docs/subprojects/metadata-qc.md](docs/subprojects/metadata-qc.md)

| Script | Purpose | Main Function |
|--------|---------|---------------|
| [check\_metadata.R](metadata_qc/check_metadata.R) | Validate study metadata, detect missing values and temporal patterns | `check_metadata(metadata_file, ...)` |
| [check\_acquisition.R](metadata_qc/check_acquisition.R) | Extract FCS acquisition dates, detect likely duplicate files | `check_acquisition_dates(parent_folder, ...)` |
| [overlap\_statistics.R](metadata_qc/overlap_statistics.R) | Fleiss' kappa for inter-rater agreement across gating strategies | `compute_fleiss_kappa(label_lists, ...)` |
| [compare\_gatings.R](metadata_qc/compare_gatings.R) | Compare multiple FlowJo gating strategies with Venn/Euler diagrams | `compare_gatings(pop_strings_list, ...)` |

---

## Configuration

Each pipeline has a dedicated config file with only the parameters it needs. Copy the relevant config, edit the paths and markers, and run.

### Pipeline Configs

| Config | Pipeline | What to Edit |
|--------|----------|--------------|
| [pipeline\_A\_config.R](config/pipeline_A_config.R) | A — Raw FCS | `fcs_folder`, `markers`, `conditions_csv`, QC/transform/clustering options |
| [pipeline\_B\_config.R](config/pipeline_B_config.R) | B — CSV Expression | `csv_folder`, `markers`, `condition_col`, `sample_col` |
| [pipeline\_C\_config.R](config/pipeline_C_config.R) | C — CSV Gating | `csv_folder`, `conditions_csv`, analysis method (`glmm` or `enrichment`) |
| [pipeline\_D\_config.R](config/pipeline_D_config.R) | D — FlowJo Comparison | `wsp_paths`, `custom_name`, `target_fcs`, `label_dict` |
| [pipeline\_E\_config.R](config/pipeline_E_config.R) | E — Preprocessed UMAP | `fcs_folder`, `markers`, downsampling/clustering options |

Usage:
```bash
Rscript pipelines/pipeline_A_raw_fcs.R --config config/pipeline_A_config.R
```

### Standalone Script Configs

For running individual scripts outside the pipelines. Each config shows the exact function call with all available parameters.

**Analysis:**

| Config | Script | Purpose |
|--------|--------|---------|
| [cluster\_splitting\_config.R](config/standalone/cluster_splitting_config.R) | [cluster\_splitting.R](analysis/cluster_splitting.R) | Split data by cluster, generate tracking matrix |
| [differential\_abundance\_config.R](config/standalone/differential_abundance_config.R) | [differential\_abundance.R](analysis/differential_abundance.R) | NB-GLMM differential abundance testing |
| [cluster\_enrichment\_config.R](config/standalone/cluster_enrichment_config.R) | [cluster\_enrichment.R](analysis/cluster_enrichment.R) | Fisher's exact / mixed-effects enrichment |
| [marker\_statistics\_config.R](config/standalone/marker_statistics_config.R) | [marker\_statistics.R](analysis/marker_statistics.R) | Per-marker summary stats with mixed models |
| [expression\_matrices\_config.R](config/standalone/expression_matrices_config.R) | [expression\_matrices.R](analysis/expression_matrices.R) | Grouped expression matrices |
| [abundance\_plotting\_config.R](config/standalone/abundance_plotting_config.R) | [abundance\_plotting.R](analysis/abundance_plotting.R) | Abundance barplots from tracking matrix |
| [cell\_ratio\_config.R](config/standalone/cell_ratio_config.R) | [cell\_ratio.R](analysis/cell_ratio.R) | Cell population ratios (e.g., CD4/CD8) |

**Visualization:**

| Config | Script | Purpose |
|--------|--------|---------|
| [umap\_config.R](config/standalone/umap_config.R) | [umap.R](visualization/umap.R) | UMAP embeddings and overlay plots |
| [heatmap\_config.R](config/standalone/heatmap_config.R) | [heatmaps.R](visualization/heatmaps.R) | Clustered expression heatmaps |
| [density\_plots\_config.R](config/standalone/density_plots_config.R) | [density\_plots.R](visualization/density_plots.R) | Per-marker density/ridge plots |
| [barplots\_config.R](config/standalone/barplots_config.R) | [barplots.R](visualization/barplots.R) | Composition barplots |
| [pca\_config.R](config/standalone/pca_config.R) | [pca.R](visualization/pca.R) | PCA scatter and loadings plots |
| [upset\_plots\_config.R](config/standalone/upset_plots_config.R) | [upset\_plots.R](visualization/upset_plots.R) | UpSet set intersection plots |
| [venn\_diagrams\_config.R](config/standalone/venn_diagrams_config.R) | [venn\_diagrams.R](visualization/venn_diagrams.R) | Venn and Euler diagrams |

**Loading & Metadata QC:**

| Config | Script | Purpose |
|--------|--------|---------|
| [load\_fcs\_config.R](config/standalone/load_fcs_config.R) | [load\_fcs\_raw.R](loading/load_fcs_raw.R) | Load raw FCS files |
| [load\_csv\_gating\_config.R](config/standalone/load_csv_gating_config.R) | [load\_csv\_gating.R](loading/load_csv_gating.R) | Load CSVs with binary gating columns |
| [load\_flowjo\_wsp\_config.R](config/standalone/load_flowjo_wsp_config.R) | [load\_flowjo\_wsp.R](loading/load_flowjo_wsp.R) | Load FlowJo workspace |
| [check\_metadata\_config.R](config/standalone/check_metadata_config.R) | [check\_metadata.R](metadata_qc/check_metadata.R) | Validate study metadata |
| [check\_acquisition\_config.R](config/standalone/check_acquisition_config.R) | [check\_acquisition.R](metadata_qc/check_acquisition.R) | Extract FCS acquisition dates |
| [compare\_gatings\_config.R](config/standalone/compare_gatings_config.R) | [compare\_gatings.R](metadata_qc/compare_gatings.R) | Compare FlowJo gating strategies |
| [overlap\_statistics\_config.R](config/standalone/overlap_statistics_config.R) | [overlap\_statistics.R](metadata_qc/overlap_statistics.R) | Fleiss' kappa inter-rater agreement |

### Other Config Files

| File | Purpose |
|------|---------|
| [config/Parameters\_template.R](config/Parameters_template.R) | Master template with all ~80 parameters (for reference) |
| [config/Packages.R](config/Packages.R) | Explicit package loading |
| [config/study\_params/](config/study_params/) | Study-specific configs (copy from pipeline configs above) |

### Key Parameters

The `params` list in pipeline configs controls everything. The most important ones:

```r
params <- list(
  fcs_folder      = "/path/to/data",          # or csv_folder for pipelines B/C
  output_dir      = "/path/to/output",
  markers         = c("CD3", "CD4", "CD8", "CD19", "CD56"),
  qc_method       = "PeacoQC",               # Pipeline A only
  transform_method = "arcsinh",               # Pipeline A only
  downsample_method = "Consistent",           # Pipelines A, E
  downsample_n    = 10000,
  resolution      = 0.8,                      # Leiden clustering
  n_neighbors     = 30,                       # KNN for Leiden/UMAP
  condition_levels = c("Control", "Treatment"),
  seed            = 42
)
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/getting-started.md](docs/getting-started.md) | First-time setup, environment installation, running your first pipeline |
| [docs/r-shared-reference.md](docs/r-shared-reference.md) | Complete API reference for all 8 shared utility modules |
| **Pipeline guides** | |
| [docs/pipelines/pipeline-A-raw-fcs.md](docs/pipelines/pipeline-A-raw-fcs.md) | Pipeline A — raw FCS end-to-end |
| [docs/pipelines/pipeline-B-csv-expression.md](docs/pipelines/pipeline-B-csv-expression.md) | Pipeline B — CSV expression matrices |
| [docs/pipelines/pipeline-C-csv-gating.md](docs/pipelines/pipeline-C-csv-gating.md) | Pipeline C — CSV gating data |
| [docs/pipelines/pipeline-D-flowjo-comparison.md](docs/pipelines/pipeline-D-flowjo-comparison.md) | Pipeline D — FlowJo workspace comparison |
| [docs/pipelines/pipeline-E-preprocessed-umap.md](docs/pipelines/pipeline-E-preprocessed-umap.md) | Pipeline E — preprocessed FCS quick clustering |
| **Subproject references** | |
| [docs/subprojects/loading.md](docs/subprojects/loading.md) | Loading module reference |
| [docs/subprojects/preprocessing.md](docs/subprojects/preprocessing.md) | Preprocessing module reference |
| [docs/subprojects/harmonization.md](docs/subprojects/harmonization.md) | Harmonization module reference |
| [docs/subprojects/clustering.md](docs/subprojects/clustering.md) | Clustering module reference |
| [docs/subprojects/analysis.md](docs/subprojects/analysis.md) | Analysis module reference |
| [docs/subprojects/visualization.md](docs/subprojects/visualization.md) | Visualization module reference |
| [docs/subprojects/metadata-qc.md](docs/subprojects/metadata-qc.md) | Metadata QC module reference |

---

## Test Data

Example input files are provided for testing and understanding expected formats.

| File | Description |
|------|-------------|
| [test\_data/examples/expression\_matrix.csv](test_data/examples/expression_matrix.csv) | Sample expression matrix with marker columns |
| [test\_data/examples/metadata.csv](test_data/examples/metadata.csv) | Sample metadata with conditions and sample IDs |
| [test\_data/examples/conditions.csv](test_data/examples/conditions.csv) | Condition assignments for differential analysis |
| [test\_data/examples/cytolution\_gating\_export.csv](test_data/examples/cytolution_gating_export.csv) | Cytolution-format CSV with binary gating columns |
| [test\_data/examples/manual\_gating\_labels.csv](test_data/examples/manual_gating_labels.csv) | Manual gating labels per cell |
| [test\_data/examples/cluster\_tracking\_matrix.csv](test_data/examples/cluster_tracking_matrix.csv) | Sample-by-cluster cell count matrix |
| [test\_data/examples/config\_example.R](test_data/examples/config_example.R) | Example parameter configuration |

See [test\_data/examples/README.md](test_data/examples/README.md) for column format specifications.

---

## Key Packages

| Package | Purpose |
|---------|---------|
| [flowCore](https://bioconductor.org/packages/flowCore/) | FCS file I/O, flowFrame/flowSet objects |
| [flowWorkspace](https://bioconductor.org/packages/flowWorkspace/) | FlowJo workspace parsing, gating sets |
| [data.table](https://rdatatable.gitlab.io/data.table/) | Fast CSV I/O (semicolon separator default) |
| [lme4](https://github.com/lme4/lme4) | Mixed-effects models (GLMM for differential abundance) |
| [uwot](https://github.com/jlmelville/uwot) | UMAP dimensionality reduction |
| [leiden](https://github.com/TomKellyGenetics/leiden) | Leiden community detection clustering |
| [ggplot2](https://ggplot2.tidyverse.org/) | Core plotting engine |
| [pheatmap](https://github.com/raivokolde/pheatmap) | Clustered heatmaps |
| [harmony](https://github.com/immunogenomics/harmony) | Batch effect correction |
| [flowAI](https://bioconductor.org/packages/flowAI/) / [PeacoQC](https://bioconductor.org/packages/PeacoQC/) | Automated quality control |
| [eulerr](https://github.com/jolars/eulerr) | Area-proportional Euler diagrams |
| [UpSetR](https://github.com/hms-dbmi/UpSetR) | UpSet plots for set intersections |

---

## Environment Files

| File | Purpose |
|------|---------|
| [environment.yml](environment.yml) | Conda environment spec — R 4.3+, 50+ packages from conda-forge/bioconda |
| [install\_extra\_packages.R](install_extra_packages.R) | Post-conda installer for GitHub-only packages (Harmony, CytoExploreR, Spectre, etc.) |
| [setup\_env.sh](setup_env.sh) | Automated setup script — creates conda env + installs extras |

---

## Design Conventions

- **Shared utilities first**: All scripts source `R_shared/` modules — never duplicate functions across subprojects
- **Function arguments, not globals**: Parameters are passed explicitly, not read from global `CONFIG` objects
- **GUI with CLI fallback**: Interactive dialogs (tcltk) automatically fall back to `readline()` in headless environments
- **data.table for I/O**: CSV read/write via `fread`/`fwrite` with auto-detected separator (semicolon default)
- **Memory-safe**: Large files use chunked reading; memory checks trigger GC at 85% utilization
- **Reproducible**: All stochastic operations accept a `seed` parameter
- **Multi-format output**: Plots saved as PDF/PNG/SVG via `save_plot_formats()` based on config toggles
