[< Back to Documentation Index](README.md)

# Configuration Guide

All pipelines are configured via YAML parameter files and Nextflow CLI overrides. Each pipeline has a template in `config/`.

---

## Creating a Config

Copy a pipeline template and edit it for your study:

```bash
cp config/pipeline_A_params.yaml config/study_params/my_study_A.yaml
# Edit paths and parameters, then run:
nextflow run nextflow/main.nf -params-file config/study_params/my_study_A.yaml -profile docker
```

You can also pass individual parameters directly on the command line:

```bash
nextflow run nextflow/main.nf --pipeline A --input_dir /data/fcs \
    --markers_of_interest "CD3,CD4,CD8" --output_dir ./results -profile docker
```

CLI `--param value` overrides take precedence over values in the YAML file.

---

## Nextflow Configuration

Pipeline execution settings (executor, container images, resource limits) live in `nextflow.config` at the project root. Study-specific analysis parameters live in YAML files under `config/`.

### `nextflow.config` (execution settings)

Controls how and where processes run (local, SLURM, AWS Batch), container images, and resource defaults. You typically do not need to edit this unless changing the execution environment.

### YAML param files (analysis settings)

All analysis parameters (input paths, markers, methods, thresholds) are specified in YAML:

```yaml
# config/pipeline_A_params.yaml
pipeline: "A"
input_dir: "/path/to/fcs/files"
output_dir: "/path/to/output"
markers_of_interest: "CD3,CD4,CD8,CD19,CD56,HLA-DR"
conditions_csv: "/path/to/conditions.csv"
condition_levels: "Control,Treatment"
truncation: false
pp_order: "QT"
qc_type: "None"
transform_type: "arcsinh"
offset_list: "0,0,0,0,0,0"
width_list: "5,5,5,5,5,5"
downsample_method: "Consistent"
downsample_n: 10000
n_neighbors: 15
resolution: 0.5
clustering_method: "leiden"
dimred_method: "umap"
batch_method: "none"
seed: 42
```

---

## Complete Parameter Reference

All parameters below are defined in `nextflow.config` with their defaults. Override them via YAML config files or `--param value` on the command line.

### Pipeline Selection

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `pipeline` | String | `null` | All | **Required.** Pipeline to run: `A`, `B`, `C`, `D`, or `E` |

### Paths

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `input_dir` | Path | `null` | A, B, C, E | **Required.** Directory containing input files (`.fcs` or `.csv`) |
| `output_dir` | Path | `./results` | All | Root output directory (subdirectories created automatically) |
| `conditions_csv` | Path | `null` | A, C | Path to CSV mapping file names to conditions |

### Markers

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `markers_of_interest` | String | `null` | A, B, C, E | **Required.** Comma-separated marker names to analyze: `"CD3,CD4,CD8"` |
| `markername_replacing` | String | `null` | A, E | JSON string for renaming markers: `'{"old_name":"new_name"}'` |

### Methods (swappable)

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `clustering_method` | String | `leiden` | A, E | Clustering algorithm: `leiden`, `flowsom`, `phenograph`, `dbscan` |
| `batch_method` | String | `none` | A, E | Batch correction method: `harmony`, `cytonorm`, `combat`, `none` |
| `dimred_method` | String | `umap` | A, E | Dimensionality reduction: `umap`, `tsne`, `diffusion_map` |
| `transform_type` | String | `arcsinh` | A | Data transformation: `arcsinh`, `logicle`, `biexponential`, `log` |
| `analysis_method` | String | `enrichment` | C | Statistical analysis approach: `enrichment` (Fisher's exact) or `abundance` (NB-GLMM) |

### Conditions / Metadata

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `condition_levels` | String | `Control,Treatment` | A, C | Comma-separated condition group names |
| `sample_names` | String | `null` | A, E | Comma-separated sample identifiers (auto-detected from filenames if null) |
| `batch_names` | String | `null` | A, E | Comma-separated batch labels for batch correction |
| `patient_names` | String | `null` | A | Comma-separated patient identifiers for paired analyses |

### Loading

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `truncation` | Boolean | `false` | A, E | Truncate FCS expression values to the `$PnR` maximum range |

### Preprocessing

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `pp_order` | String | `QT` | A | Preprocessing step order: `QT` (QC then Transform) or `TQ` (Transform then QC) |
| `qc_type` | String | `None` | A | Quality control method: `None`, `flowAI`, `PeacoQC` |
| `quant_norm` | Boolean | `false` | A | Apply quantile normalization after transformation |
| `scaling` | Boolean | `false` | A | Apply z-score scaling to expression values |

### Downsampling

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `downsample_method` | String | `Consistent` | A, E | Downsampling strategy: `Consistent` (fixed N per sample) or `Maximal` (proportional to sample size) |
| `downsample_n` | Integer | `10000` | A, E | Target number of cells per sample after downsampling |

### Clustering

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `n_neighbors` | Integer | `15` | A, E | Number of nearest neighbors for KNN graph construction (used by UMAP and Leiden) |
| `resolution` | Float | `0.5` | A, E | Leiden clustering resolution; higher values produce more clusters |
| `n_clusters` | Integer | `20` | A, E | Target number of meta-clusters for FlowSOM |
| `grid_xdim` | Integer | `10` | A, E | FlowSOM SOM grid X dimension |
| `grid_ydim` | Integer | `10` | A, E | FlowSOM SOM grid Y dimension |
| `dbscan_eps` | Float | `null` | A, E | DBSCAN epsilon (neighborhood radius); `null` triggers automatic estimation via k-distance plot |
| `dbscan_min_pts` | Integer | `5` | A, E | DBSCAN minimum number of points to form a dense region |

### Batch Correction

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `batch_do_pca` | Boolean | `false` | A, E | Run PCA before Harmony batch correction |
| `combat_method` | String | `parametric` | A, E | ComBat estimation method: `parametric` or `nonparametric` |
| `cytonorm_n_clusters` | Integer | `10` | A, E | Number of FlowSOM clusters used internally by CytoNorm for reference alignment |
| `cytonorm_goal` | String | `mean` | A, E | CytoNorm alignment target: `mean` (align to batch mean) or a specific batch name |

### Dimensionality Reduction

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `umap_min_dist` | Float | `0.1` | A, E | UMAP minimum distance; smaller values produce tighter clusters in the embedding |
| `tsne_perplexity` | Integer | `30` | A, E | t-SNE perplexity; roughly the expected number of neighbors per point (typical range: 5-50) |
| `tsne_max_iter` | Integer | `1000` | A, E | Maximum number of t-SNE optimization iterations |

### Analysis

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `cluster_regex` | String | `Clustering annotated` | A, C | Regex pattern identifying gating/cluster columns in CSV files |
| `gating_pattern` | String | `\\(Clustering annotated\\)$\|^Clustering raw$` | C | Extended regex for matching gating hierarchy paths |
| `processing_mode` | String | `manual` | C | Gating data processing mode: `manual` |
| `min_cluster_size` | Integer | `5` | A, B, C | Minimum number of cells in a cluster to include it in statistical analysis |
| `p_adjust_method` | String | `BH` | A, B, C | Multiple-testing correction: `BH` (Benjamini-Hochberg), `bonferroni`, `holm`, `fdr`, `BY` |
| `alpha_level` | Float | `0.05` | A, B, C | Significance threshold for adjusted p-values |
| `stat_method` | String | `pooled` | A, C | Abundance analysis mode: `pooled` (aggregate across samples) or `mixed_effects` (per-sample random effects) |
| `test_type` | String | `lmer` | A, B, C | Statistical test for differential analysis: `lmer` (mixed-effects), `wilcox` (Wilcoxon), `kruskal` (Kruskal-Wallis) |
| `exclude_patterns` | String | _(empty)_ | C | Comma-separated patterns of gate/cluster names to exclude from analysis |
| `cytolution_metadata_file` | Path | `null` | A | Path to Cytolution-format metadata file (optional external metadata source) |

### Pipeline D Specific

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `wsp_paths` | String | `null` | D | **Required for D.** Comma-separated paths to FlowJo `.wsp` workspace files |
| `wsp_names` | String | `null` | D | **Required for D.** Comma-separated analyst/workspace display names (must match `wsp_paths` count) |
| `target_fcs` | Path | `null` | D | **Required for D.** Path to the target `.fcs` file to compare gating across workspaces |
| `label_dict` | String | `null` | D | JSON string mapping gate hierarchy paths to display names: `'{"gate/path":"Display Name"}'` |

### Output

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `output_formats` | String | `png` | All | Comma-separated output formats for plots: `png`, `pdf`, `svg`, or combinations like `png,pdf` |
| `plot_width` | Integer | `10` | All | Default plot width in inches |
| `plot_height` | Integer | `8` | All | Default plot height in inches |
| `save_intermediates` | Boolean | `false` | A, E | Save intermediate RDS files (flowFrame objects, matrices) for debugging or downstream reuse |

### Memory

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `max_cells_per_file` | Integer | `500000` | A, E | Maximum cells to read per FCS file before triggering chunked processing |
| `max_cells_heatmap` | Integer | `5000` | A, B, E | Maximum cells to include in heatmap visualizations (random subsample applied if exceeded) |

### Reproducibility

| Parameter | Type | Default | Pipelines | Description |
|-----------|------|---------|-----------|-------------|
| `seed` | Integer | `42` | All | Random seed for reproducible results across runs |

---

## Example: Pipeline A YAML Config

```yaml
# config/pipeline_A_params.yaml
pipeline: "A"
input_dir: "/path/to/fcs/files"
output_dir: "/path/to/output"
markers_of_interest: "CD3,CD4,CD8,CD19,CD56,HLA-DR"
conditions_csv: "/path/to/conditions.csv"
condition_levels: "Control,Treatment"
truncation: false
pp_order: "QT"
qc_type: "None"
transform_type: "arcsinh"
offset_list: "0,0,0,0,0,0"
width_list: "5,5,5,5,5,5"
downsample_method: "Consistent"
downsample_n: 10000
n_neighbors: 15
resolution: 0.5
seed: 42
```

See each [pipeline guide](pipelines/) for pipeline-specific config examples and parameter details.
