[< Back to Documentation Index](../README.md)

# Pipeline D: FlowJo Workspace Comparison

Pipeline D is designed for comparing gating strategies across multiple FlowJo workspaces. Given two or more `.wsp` files that gate the same FCS data (e.g. from different analysts or different gating approaches), it extracts population assignments, computes inter-rater agreement using Fleiss' kappa, and visualizes the overlap between strategies using Venn/Euler diagrams and UpSet plots.

This is the only pipeline in the framework that does not perform any expression analysis or clustering. It is purely a gating agreement and quality assessment tool.

---

## When to Use Pipeline D

Use this pipeline when:

- You have **multiple FlowJo `.wsp` files** that gate the same FCS data
- You want to assess **inter-rater agreement** between different analysts
- You want to compare **manual vs automated gating** strategies
- You need **Venn/Euler diagrams** showing overlap between gating sets

Pipeline D requires **real FlowJo workspace files** and matching FCS data. It cannot be tested with synthetic data.

---

## Pipeline Steps

1. **`load_multiple_workspaces()`** — Load all WSP files and extract per-event population labels
2. **`compare_gatings()`** — Align populations across workspaces, generate Venn/Euler per population
3. **`compute_fleiss_kappa()`** — Inter-rater agreement statistics

---

## Input

### FlowJo Workspaces

Two or more `.wsp` files, each representing a different gating strategy or analyst. All workspaces must reference the same FCS file(s).

### FCS Files

The `.fcs` files referenced by the workspaces. These must be accessible at the paths stored in the workspace, or in a directory you specify.

### Directory Layout

```
study_data/
├── analyst_A.wsp           # Workspace from analyst A
├── analyst_B.wsp           # Workspace from analyst B
├── sample_001.fcs          # Shared FCS file
├── sample_002.fcs
└── sample_003.fcs
```

---

## Configuration

### Running the Pipeline

```bash
# Minimal invocation with inline parameters
nextflow run nextflow/main.nf --pipeline D \
    --wsp_paths "/data/analyst_A.wsp,/data/analyst_B.wsp" \
    --wsp_names "Analyst_A,Analyst_B" --output_dir ./results -profile docker

# Using a YAML parameter file (recommended for reproducibility)
nextflow run nextflow/main.nf -params-file config/pipeline_D_params.yaml -profile docker
```

**A params file is strongly recommended for this pipeline**, because workspace paths and label mappings are inherently study-specific.

Pipeline-level defaults are defined in `nextflow.config`. Any parameter can be overridden on the command line or via a YAML params file.

### Example Params File

```yaml
# config/pipeline_D_params.yaml

pipeline:      "D"

# --- Workspace paths ---
wsp_paths:     "/path/to/analyst_A.wsp,/path/to/analyst_B.wsp"
wsp_names:     "Analyst_A,Analyst_B"    # Short labels for legends

# --- Target FCS ---
target_fcs:    "sample_001.fcs"          # Compare this FCS across workspaces
fcs_dir:       "/path/to/fcs/files"      # Where FCS files are located

# --- Label mapping ---
label_dict:
  "/no erys/Single Cells/lymphocytes/CD3+":                      "CD3+"
  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+":         "CD4+ T cells"
  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+":         "CD8+ T cells"
  "/no erys/Single Cells/lymphocytes/CD3-/CD19+":                "B cells"
  "/no erys/Single Cells/lymphocytes/CD3-/CD56+":                "NK cells"

# --- Optional: pre-computed labels CSV ---
idx_labels:    null     # or "/path/to/labels.csv"

# --- Output ---
output_dir:    "/path/to/output/pipeline_D_results"

# --- Options ---
scale_size:    true     # true = Euler (proportional), false = Venn
```

### Key Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `wsp_paths` | comma-separated string | Paths to `.wsp` files to compare |
| `wsp_names` | comma-separated string | Short nicknames for each workspace (used in legends) |
| `target_fcs` | string | FCS filename to compare across workspaces |
| `fcs_dir` | string | Directory containing FCS files |
| `label_dict` | YAML map | Maps FlowJo gate hierarchy paths to simplified display labels |
| `idx_labels` | string or `null` | Optional CSV with pre-computed labels (alternative to WSP labels) |
| `output_dir` | string | Output directory |
| `scale_size` | boolean | `true` for Euler (area-proportional), `false` for Venn |

### Understanding `label_dict`

The `label_dict` maps verbose FlowJo hierarchy paths to concise labels:

```r
label_dict <- c(
  "/no erys/Single Cells/lymphocytes/CD3+"                              = "CD3+",
  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+"                 = "CD4+ T cells",
  "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+"                 = "CD8+ T cells",
  "/no erys/Single Cells/lymphocytes/CD3-/CD19+"                        = "B cells",
  "/no erys/Single Cells/lymphocytes/CD3-/CD56+"                        = "NK cells"
)
```

**Only gates listed in `label_dict` are included in the comparison.** Gates not listed are silently skipped. This lets you focus on the populations you care about.

To find the gate paths in your workspace, use:

```r
source("src/loading/load_flowjo_wsp.R")
ws <- load_flowjo_wsp("analyst_A.wsp")
print(ws$population_names)
# [1] "/no erys" "/no erys/Single Cells" "/no erys/Single Cells/lymphocytes" ...
```

---

## Steps in Detail

### Step 1: Load FlowJo Workspaces

```r
source(file.path(base_dir, "src", "loading", "load_flowjo_wsp.R"))

pop_strings_list <- load_multiple_workspaces(wsp_paths = wsp_paths, names = custom_name)
```

**What it does:** Iterates over all workspace paths, and for each one:

1. Opens the `.wsp` file with `CytoML::open_flowjo_xml()`
2. Builds a GatingSet linking the workspace gates to the FCS data
3. Extracts per-event population assignments using `gh_pop_get_indices()`
4. Maps gate paths through `label_dict` to produce simplified labels

Returns a named list (keyed by `custom_name`) of per-sample population string lists.

**`load_multiple_workspaces()` signature:**

```r
load_multiple_workspaces(wsp_paths, names = NULL, fcs_dir = NULL)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `wsp_paths` | character vector | Paths to `.wsp` files |
| `names` | character vector | Short labels for each workspace (defaults to filenames) |
| `fcs_dir` | character | Directory containing FCS files (defaults to each WSP's parent directory) |

Returns a named list of population string lists, one entry per workspace. Internally calls `load_flowjo_wsp()` and `get_population_strings_list()` for each workspace.

---

### Step 2: Compare Gatings

```r
source(file.path(base_dir, "src", "metadata_qc", "compare_gatings.R"))
comparison <- compare_gatings(
  pop_strings_list = all_labels,
  target_fcs       = params$target_fcs,
  label_dict       = params$label_dict,
  custom_name      = params$custom_name,
  output_dir       = file.path(params$output_dir, "01_gating_comparison"),
  scale_size       = params$scale_size
)
```

**What it does:**

1. Extracts labels for `target_fcs` from each workspace's population strings
2. Translates labels via `label_dict` to a common vocabulary
3. Aligns events by index across workspaces
4. For each population, creates a binary matrix (event × workspace) indicating membership
5. Generates a **Venn or Euler diagram** for each population showing overlap

**Function signature:**

```r
compare_gatings(pop_strings_list, target_fcs, label_dict,
                custom_name = NULL, output_dir,
                idx_labels = NULL, scale_size = TRUE)
```

**Output files:**

| File | Description |
|------|-------------|
| `fj_label_lists_{target_fcs}.csv` | Raw FlowJo labels per workspace |
| `label_lists_{target_fcs}.csv` | Translated labels per workspace |
| `CD4_T_cells_{target_fcs}.png` | Venn/Euler diagram for CD4+ T cells |
| `CD8_T_cells_{target_fcs}.png` | Venn/Euler diagram for CD8+ T cells |
| `B_cells_{target_fcs}.png` | Venn/Euler diagram for B cells |
| `...` | One diagram per population in `label_dict` |

**Venn vs Euler:**

- **Venn diagram** (`scale_size = FALSE`): Fixed-size circles with overlaps. Works with any number of sets.
- **Euler diagram** (`scale_size = TRUE`): Area-proportional circles. More visually informative but limited to 2–5 sets.

---

### Step 3: Overlap Statistics (Fleiss' Kappa)

```r
source(file.path(base_dir, "src", "metadata_qc", "overlap_statistics.R"))
kappa <- compute_fleiss_kappa(
  label_lists = comparison$aligned_labels,
  output_dir  = file.path(params$output_dir, "02_overlap_stats")
)
```

**What it does:** Computes Fleiss' kappa for each population across all workspaces. Kappa measures inter-rater reliability — how well the different analysts agree beyond what would be expected by chance.

**Function signature:**

```r
compute_fleiss_kappa(label_lists, output_dir,
                     subset_specs = NULL, subset_labels = NULL,
                     row_range = NULL)
```

**Output** — `fleiss_kappa_log.txt`:

```
=== Fleiss' Kappa Results ===

Population: CD4+ T cells
  Kappa = 0.85 (Almost perfect agreement)
  n_events = 50000, n_raters = 2

Population: CD8+ T cells
  Kappa = 0.72 (Substantial agreement)
  n_events = 50000, n_raters = 2

Population: B cells
  Kappa = 0.91 (Almost perfect agreement)
  n_events = 50000, n_raters = 2
```

### Kappa Interpretation Guide

| Kappa | Agreement Level |
|-------|----------------|
| < 0.00 | Less than chance |
| 0.01 – 0.20 | Slight |
| 0.21 – 0.40 | Fair |
| 0.41 – 0.60 | Moderate |
| 0.61 – 0.80 | Substantial |
| 0.81 – 1.00 | Almost perfect |

---

## Output Directory Structure

```
output/pipeline_D_results/
├── 01_gating_comparison/
│   ├── fj_label_lists_sample_001.csv     # Raw FlowJo labels
│   ├── label_lists_sample_001.csv        # Translated labels
│   ├── CD4_T_cells_sample_001.png        # Venn/Euler per population
│   ├── CD8_T_cells_sample_001.png
│   ├── B_cells_sample_001.png
│   ├── NK_cells_sample_001.png
│   └── upset_plot.png                    # UpSet intersection plot
└── 02_overlap_stats/
    └── fleiss_kappa_log.txt              # Per-population kappa values
```

---

## Tips

- **All workspaces must gate the same FCS file** (specified by `target_fcs`). If a workspace does not contain that file, the pipeline stops with an error.
- **`idx_labels` alternative input:** If you have pre-computed per-event labels from a previous Pipeline C run (or any external source), you can include them as an additional "rater" via the `idx_labels` parameter. This lets you compare automated clustering (e.g. Cytolution) against manual FlowJo gating.
- **Euler diagram limits:** Area-proportional Euler diagrams work best with 2–5 sets. For more workspaces, only the UpSet plot is generated.
- **Finding gate paths:** Use `get_population_stats()` or print `ws$population_names` after loading a workspace to see all available gate paths.
- Pipeline D does **not** perform expression analysis, clustering, or differential testing. It is strictly a gating comparison tool.

---

## Platform Notes

### Windows Limitation

Pipeline D depends on `flowWorkspace` and `CytoML`, which require HDF5 libraries that are not natively available on Windows. **On Windows, this pipeline must be run through WSL2 or Docker:**

```bash
# Option 1: Docker via Nextflow (recommended)
nextflow run nextflow/main.nf -params-file config/pipeline_D_params.yaml -profile docker

# Option 2: WSL2 (Ubuntu)
wsl -d Ubuntu
# Then run Nextflow or Rscript as normal inside WSL2
```

On Linux and macOS, the pipeline runs natively without any special configuration.

### Available Parameters Reference

In addition to the key parameters documented above, the following options are available:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `scale_size` | boolean | `true` | `true` = Euler (area-proportional diagrams), `false` = Venn (fixed-size circles) |
| `idx_labels` | string | `null` | Path to a CSV with pre-computed per-event labels from Pipeline C or external clustering. Included as an additional "rater" alongside WSP gatings. |
| `label_dict` | YAML map | *required* | Maps FlowJo gate hierarchy paths to simplified display labels. Only gates listed here are compared. |
| `output_format_pdf` | boolean | `false` | Save plots as PDF |
| `output_format_png` | boolean | `true` | Save plots as PNG |
| `output_format_svg` | boolean | `false` | Save plots as SVG |

---

## Required Packages

This pipeline requires packages not needed by others:

| Package | Purpose |
|---------|---------|
| `flowWorkspace` | Parse FlowJo workspaces, build GatingSet |
| `CytoML` | FlowJo XML parsing (`open_flowjo_xml()`) |
| `irr` | Fleiss' kappa computation (`kappam.fleiss()`) |
| `eulerr` | Area-proportional Euler diagrams |
| `UpSetR` | UpSet intersection plots |
| `venn` | Non-proportional Venn diagrams |
