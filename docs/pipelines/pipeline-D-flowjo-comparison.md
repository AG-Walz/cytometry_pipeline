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

1. **`load_flowjo_wsp()`** + **`get_population_strings_list()`** — Extract gating sets from each WSP, get per-event population labels
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
cd pipelines
Rscript pipeline_D_flowjo_comparison.R --config ../config/study_params/my_comparison.R
```

**The `--config` flag is mandatory for this pipeline.** The pipeline calls `stop()` if no config file is provided, because workspace paths and label mappings are inherently study-specific. Example config templates:

- `config/study_params/Compare_gatings_Endo_CD45_Julia.R`
- `config/study_params/Compare_gatings_ICK_broad_narrow.R`

### Example Config

```r
# config/study_params/my_comparison.R

params <- list(
  # --- Workspace paths ---
  wsp_paths    = c(
    "/path/to/analyst_A.wsp",
    "/path/to/analyst_B.wsp"
  ),
  custom_name  = c("Analyst_A", "Analyst_B"),    # Short labels for legends

  # --- Target FCS ---
  target_fcs   = "sample_001.fcs",               # Compare this FCS across workspaces
  fcs_dir      = "/path/to/fcs/files",           # Where FCS files are located

  # --- Label mapping ---
  label_dict   = c(
    "/no erys/Single Cells/lymphocytes/CD3+"                      = "CD3+",
    "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD4+"         = "CD4+ T cells",
    "/no erys/Single Cells/lymphocytes/CD3+/T cells/CD8+"         = "CD8+ T cells",
    "/no erys/Single Cells/lymphocytes/CD3-/CD19+"                = "B cells",
    "/no erys/Single Cells/lymphocytes/CD3-/CD56+"                = "NK cells"
  ),

  # --- Optional: pre-computed labels CSV ---
  idx_labels   = NULL,     # or "/path/to/labels.csv"

  # --- Output ---
  output_dir   = "/path/to/output/pipeline_D_results",

  # --- Options ---
  scale_size   = TRUE      # TRUE = Euler (proportional), FALSE = Venn
)
```

### Key Config Variables

| Variable | Type | Description |
|----------|------|-------------|
| `wsp_paths` | character vector | Paths to `.wsp` files to compare |
| `custom_name` | character vector | Short nicknames for each workspace (used in legends) |
| `target_fcs` | character | FCS filename to compare across workspaces |
| `fcs_dir` | character | Directory containing FCS files |
| `label_dict` | named character vector | Maps FlowJo gate hierarchy paths to simplified display labels |
| `idx_labels` | character or `NULL` | Optional CSV with pre-computed labels (alternative to WSP labels) |
| `output_dir` | character | Output directory |
| `scale_size` | logical | `TRUE` for Euler (area-proportional), `FALSE` for Venn |

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

all_labels <- list()
for (i in seq_along(params$wsp_paths)) {
  ws <- load_flowjo_wsp(params$wsp_paths[i], fcs_dir = params$fcs_dir)
  pops <- get_population_strings_list(ws$gating_set, ws$sample_names)
  all_labels[[params$custom_name[i]]] <- pops
}
```

**What it does:** For each workspace:

1. Opens the `.wsp` file with `CytoML::open_flowjo_xml()`
2. Builds a GatingSet linking the workspace gates to the FCS data
3. Extracts per-event population assignments using `gh_pop_get_indices()`
4. Maps gate paths through `label_dict` to produce simplified labels

**Function signatures:**

```r
load_flowjo_wsp(wsp_path, fcs_dir = dirname(wsp_path), name = NULL)
```

Returns:

```r
list(
  gating_set       = <GatingSet>,
  sample_names     = c("sample_001", "sample_002"),
  population_names = c("/no erys", "/no erys/Single Cells", ...),
  workspace_name   = "Analyst_A"
)
```

```r
get_population_strings_list(gating_set, sample_names)
```

Returns a named list of data.frames, one per sample, with a `labels` column containing the gate path for each event.

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
