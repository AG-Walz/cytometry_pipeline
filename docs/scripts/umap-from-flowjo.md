[< Back to Documentation Index](../README.md)

# UMAP from a FlowJo Workspace Population

`scripts/umap_from_flowjo.R` turns one **or many** FlowJo workspaces plus their
raw FCS files into a single shared UMAP, coloured by sample, by **workspace**,
by **FlowJo subgate**, and by every marker.

The script itself is a thin CLI wrapper: it parses options, locates the repo,
and calls **`umap_from_flowjo()`** in
[`src/visualization/umap_from_flowjo.R`](../../src/visualization/umap_from_flowjo.R),
where all the analysis and plotting logic lives. See
[Calling it from R](#calling-it-from-r) to use that function directly.

It:

1. Opens every FlowJo `.wsp` it is pointed at and builds a `GatingSet` per
   workspace over the FCS files that workspace references.
2. Extracts the cells of **one selected gated population** from each sample.
   The data is taken **as-is** from FlowJo's working space (already compensated
   and transformed), so **no further transform is applied**.
3. Discovers the **subgates underneath** that population and assigns each cell a
   subpopulation label (see [Subpopulation colouring](#subpopulation-colouring)).
4. Concatenates the per-sample expression matrices, optionally normalising each
   marker (`--normalize`).
5. Computes a **single shared UMAP embedding** over the concatenated data —
   one embedding for everything, however many workspaces contributed.
6. Plots that embedding coloured by **sample**, by **workspace** (when there is
   more than one), by **subpopulation**, and by **every marker** — markers both
   **conjointly** (all files together) and **per file** (each file on the same
   coordinates, with a shared colour scale for direct comparison).

It reuses the repo's shared utilities (`run_umap`, `save_plot_formats`,
`create_safe_filename`, `make_outdir`, `explore_data`) rather than
reimplementing them.

---

## Prerequisites

The script runs inside the `cytometry_env` conda environment (it needs
`CytoML`, `flowWorkspace`, `flowCore`, `uwot`, `ggplot2`, `viridis`,
`optparse`). If you have not set the environment up yet, follow the
[Getting Started](../getting-started.md) guide first, or on a fresh machine:

```bash
git clone https://github.com/AG-Walz/cytometry_pipeline.git
cd cytometry_pipeline
bash setup_env.sh
conda activate cytometry_env
```

> **Note:** raw data (`*.fcs`, `*.wsp`) is gitignored and is **not** included in
> a clone. Copy your data into `data/` (or point `--wsp` / `--fcs_dir` elsewhere)
> before running.

> **`undefined symbol` errors (e.g. `SETLENGTH`, `NO_REFERENCES`)?** A global
> `~/.Rprofile` or `~/.Renviron` that pins `.libPaths()` / `R_LIBS_USER` to a
> different R version's library (such as `~/R/x86_64-pc-linux-gnu-library/4.4`)
> shadows the conda environment and loads ABI-incompatible packages. `setup_env.sh`
> now neutralises those files and installs a conda activation hook so the env is
> self-contained — just (re)run `bash setup_env.sh`. If you invoke R outside the
> activated env, prefix it with `R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null`.

---

## Usage

Always run from the repository root with the environment activated.

**1. Inspect the workspace so you can pick a population:**

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp data/Killing_AML.wsp \
  --fcs_dir data \
  --list_populations
```

This loads the workspace once and prints everything you need to configure a
run — the samples in the group, every population path, and the marker panel
common to all samples — then exits without analysing anything:

```
Samples in loaded group (2):
  - Killing_AML_Pool_001.fcs_0
  - Killing_AML_Pool_002.fcs_1

Population tree (37 nodes):
  /Lymphocytes
  /Lymphocytes/Single Cells
  ...

Common markers (12):
  CD3
  CD4
  ...
```

Markers are a property of the channels rather than of a gate, so they are
evaluated at the root node and are the same set a real run would use. These
three block headers are a stable format that [the Dash front-end](../../app/README.md)
parses — if you change them, update `app/services/rscript.py` too. With several
workspaces two further blocks follow, indented four spaces so those parsers
cannot swallow them; see [Gating trees that differ](#gating-trees-that-differ).

**2. Run the analysis** (defaults already target the bundled example data and
the `Lymphocytes/Single Cells/Viable cells` population):

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp data/Killing_AML.wsp \
  --fcs_dir data \
  --population "Lymphocytes/Single Cells/Viable cells" \
  --output_dir results/umap_viable
```

A leading `/` on `--population` is optional, and a unique suffix
(e.g. `"Single Cells/Viable cells"`) is matched against the tree.

**Common variations:**

```bash
# Restrict to / reorder a subset of markers (must be common to both files)
Rscript scripts/umap_from_flowjo.R --markers "CD3,CD4,CD8,CD19" \
  --output_dir results/umap_tcells

# Colour the subpopulation plot by the CD4+/CD8+ lineage split only
Rscript scripts/umap_from_flowjo.R --subpop_mode cd4cd8

# Cap cells per sample AFTER gating (faster UMAP on large populations)
Rscript scripts/umap_from_flowjo.R --downsample 20000

# z-score each marker before the UMAP, also emit PDF + PNG
Rscript scripts/umap_from_flowjo.R --normalize zscore \
  --output_formats png,pdf,svg --n_neighbors 30 --min_dist 0.3
```

---

## Many workspaces at once

`--wsp` also takes a **folder**, searched recursively for `.wsp`. Every
workspace it finds contributes its cells to **one shared embedding**, so their
coordinates are directly comparable, and an extra `..._by_workspace` plot
colours the embedding by which workspace each cell came from.

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp  data/workspaces \
  --fcs_dir data/fcs \
  --group "All Samples" \
  --population "Single Cells/Viable cells" \
  --output_dir results/umap_cohort
```

> **Batch effects.** Pooling workspaces from different acquisition runs into one
> embedding does **not** correct for batch — nothing here does. If the
> `..._by_workspace` plot shows each workspace on its own island, the embedding
> is partly driven by acquisition rather than biology; see
> [`src/harmonization/`](../../src/harmonization/) (ComBat, CytoNorm) for
> correction, applied before this step.

### How workspaces are paired with their FCS files

A `.wsp` records a full `uri` per sample, but it is the path on the Windows
machine that built the workspace, so it does not exist on the analysis host.
What a workspace reliably gives is the **list of FCS file names** — so those are
read out of the XML first, `--fcs_dir` is indexed **recursively** by file name,
and the two are matched.

That means **no directory convention is imposed on your data**. A flat dump,
mirrored subfolders and a messy mix all work, and one workspace's files may even
live in several directories:

```
--wsp data/workspaces          --fcs_dir data/fcs
  workspaces/runA.wsp            fcs/2024/runA/Pool_001.fcs
  workspaces/2025/runB.wsp       fcs/2024/runA/Pool_002.fcs
                                 fcs/2025/Pool_003.fcs
```

Name matching is tolerant: case-insensitive, with `.fcs` added if the workspace
omitted it, and Windows backslash paths reduced to their base name.

Everything is resolved **before any GatingSet is built**, and every problem in
every workspace is reported together — a GatingSet load is slow, and finding one
mislaid file per load is not a workflow:

```
Error: could not build the workspace set.

  runB (/data/workspaces/2025/runB.wsp)
    missing:   Pool_007.fcs, Pool_008.fcs
  runC (/data/workspaces/runC.wsp)
    ambiguous: Pool_001.fcs matches
                 /data/fcs/2024/Pool_001.fcs
                 /data/fcs/2025/Pool_001.fcs

FCS root searched (recursively): /data/fcs
Nothing was loaded. Fix the FCS root (or --group) and re-run.
```

### Sample groups across workspaces

The group of interest can sit at a different index in each workspace, so
`--group` accepts a **name** as well as an index:

```bash
--group 1               # index, as always
--group "All Samples"   # name, resolved separately in each workspace
```

`--list_populations` prints each workspace's groups, so a name that does not
resolve tells you what the valid ones are.

### Gating trees that differ

`--population` is resolved **independently in each workspace**, using the same
suffix matching as always. Two workspaces may therefore gate the same biology
under differently-named ancestors — `/Lymphocytes/Single Cells/Viable cells` in
one and `/Lymph/Singlets/Viable cells` in the other — and still line up.

A population that resolves in *some* workspaces but not all is a **hard error
naming the workspace**. Pooling cells from gates that are not the same gate
would be a silent scientific error, so it is never done quietly.

Subgate colouring works the same way: the colour classes are the gate paths
**relative to `--population`**, which is the only key comparable across
workspaces, and they are unioned. A gate that only some workspaces have simply
has no cells from the others.

`--list_populations` reflects all of this. The `Population tree` block lists
only the gates **shared by every workspace** — the ones a run can actually
use — and a separate block names the ones that are not:

```
Workspaces resolved (2):
    runA | /data/workspaces/runA.wsp | group 1 "All Samples" | 2 sample(s)
    runB | /data/workspaces/2025/runB.wsp | group 1 "All Samples" | 2 sample(s)

Populations not shared by all workspaces (2):
    /Lymphocytes/Single Cells/Viable cells/CD8+   [missing in: runB]
    /Lymphocytes/Single Cells/Viable cells/NK     [missing in: runA]
```

### Sample names

A sample keeps its own name whenever it is unique across the whole run. Only
names that genuinely **collide** between workspaces are qualified, so a
single-workspace run produces exactly the legends and file names it always did:

```
no collision:            Pool_001, Pool_002, Pool_003
Pool_001 in two wsp:     runA: Pool_001, runB: Pool_001, Pool_002
```

This matters on disk as well as in legends: two identically-named samples would
otherwise share one set of per-file plots, the second silently overwriting the
first.

---

## Subpopulation colouring

Beyond colouring by sample and marker, the script produces a
`..._by_subpopulation` plot that colours each cell on the **same embedding** by
the FlowJo subgate it belongs to. The subgate set is chosen by `--subpop_mode`
(or overridden with `--subpop_select`).

**The overlap problem.** Under a typical population the child gates are **not
mutually exclusive** — e.g. a cell sits in `CD107a-`, `IFNg-`, `IL2-`, `TNF-`,
*and* the combined `CD107a-IFNg-IL2-TNF-` gate all at once. Colouring "by every
child gate" therefore can't give each cell one colour. The modes below pick a
set that *does* yield a single label per cell.

| `--subpop_mode` | Each cell is labelled… | Notes |
|-----------------|------------------------|-------|
| `furthest` *(default)* | by the **most specific / furthest gating-tree node it reaches** | A deeper node is a subset of its parents (fewer cells), so this uniquely picks e.g. `CD4+/CD107a-IFNg-IL2-TNF-` over the overlapping single-marker gates. Falls back to a parent gate when no deeper node is reached. Lineage (CD4+/CD8+) is encoded in the node path. No `(multiple)`. |
| `cd4cd8` | `CD4+`, `CD8+`, or `(none)` | The clean lineage split — fewest colours, most readable. |
| `direct` | by its single immediate child gate, else `(multiple)` | Immediate children overlap heavily, so most cells become `(multiple)`. |
| `deepest` | by the deepest gate **by path depth**, ties → `(multiple)` | Sibling single-marker and Boolean gates share a path depth, so this also tends to collapse to `(multiple)`. |
| `none` | — | Skips the subpopulation plot entirely. |

Cells in **no** selected subgate are labelled `(none)`; under `direct`/`deepest`,
cells in several equally-deep gates are labelled `(multiple)`. These catch-all
levels are drawn in grey so real subgates keep distinguishable colours, and the
legend reports the cell count per label.

**Explicit gate set.** To colour by a specific mutually-exclusive set (e.g. the
16 Boolean functional combos), pass them with `--subpop_select`, overriding
`--subpop_mode`. Entries match either the path **relative** to `--population` or
the bare gate name:

```bash
Rscript scripts/umap_from_flowjo.R --subpop_select \
  "CD107a+IFNg+IL2+TNF+,CD107a+IFNg+IL2+TNF-,...,CD107a-IFNg-IL2-TNF-"
```

> **Reading the plot.** Real data is usually dominated by a few "all-negative"
> nodes while functional-positive combos are rare (often a handful of cells), so
> the embedding is mostly two or three colours with tiny sprinkles of the rest —
> that reflects the biology, not a plotting fault. For spotting rare responders,
> a per-marker view (one clean `+/-` plot per functional marker) is more
> effective than one many-colour plot.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--wsp` | `data/Killing_AML.wsp` | A `.wsp` file, **or a folder** searched recursively for `.wsp` — see [Many workspaces at once](#many-workspaces-at-once) |
| `--fcs_dir` | `data` | Folder holding the FCS files, searched **recursively**; each workspace is paired with its own files by name |
| `--population` | `Lymphocytes/Single Cells/Viable cells` | Population path to extract (leading `/` optional; suffix-matched). Resolved independently in each workspace |
| `--list_populations` | off | Load the workspace(s), print their samples, population tree, common markers and workspace/divergence blocks, then exit (no analysis) |
| `--group` | `1` | FlowJo sample group to load: an index, **or a group name** such as `"All Samples"` |
| `--expected_samples` | `0` | If `> 0`, warn when the run does not total this many samples across all workspaces |
| `--markers` | all common | Comma-separated subset of marker names to use/plot |
| `--subpop_mode` | `furthest` | Subgate set for the subpopulation plot: `furthest`, `cd4cd8`, `direct`, `deepest`, `none` (see [above](#subpopulation-colouring)) |
| `--subpop_select` | — | Comma-separated explicit gate list (relative path or bare name) used as the colour partition, overriding `--subpop_mode` |
| `--normalize` | `none` | Per-marker normalisation from the pooled gated cells: `none`, `zscore`, `minmax`. Applied before UMAP and used for the marker colour scales |
| `--output_dir` | `umap_flowjo` | Output directory for plots |
| `--output_formats` | `svg` | Comma-separated formats: `png,pdf,svg` |
| `--downsample` | `0` | Max cells per sample after gating (`0` = no downsampling) |
| `--n_neighbors` | `15` | UMAP `n_neighbors` |
| `--min_dist` | `0.1` | UMAP `min_dist` |
| `--metric` | `euclidean` | UMAP distance metric |
| `--seed` | `42` | Random seed (UMAP + downsampling) |
| `--point_size` | `0.4` | ggplot point size |
| `--point_alpha` | `0.6` | ggplot point alpha |
| `--plot_width` | `10` | Plot width (inches) |
| `--plot_height` | `8` | Plot height (inches) |
| `--src_dir` | auto | Path to the repo `src` directory (auto-detected from the script location) |

Run `Rscript scripts/umap_from_flowjo.R --help` for the same list at the
command line.

---

## Output

All files are written to `--output_dir` (created if needed). Plots have no
background grid. File names are sanitised from the population and marker names.

| File | Content |
|------|---------|
| `umap_<pop>_by_file.<fmt>` | Shared embedding, points coloured categorically by source file |
| `umap_<pop>_by_workspace.<fmt>` | Shared embedding, points coloured by source **workspace**. Only written when more than one workspace was loaded — see [Many workspaces at once](#many-workspaces-at-once) |
| `umap_<pop>_by_subpopulation.<fmt>` | Shared embedding, points coloured by FlowJo subgate (per `--subpop_mode`); omitted when `--subpop_mode none` finds no subgates |
| `umap_<pop>_marker_<M>_ALL.<fmt>` | **Conjoint** plot per marker `M` — all cells, viridis-coloured by expression |
| `umap_<pop>_marker_<M>_file_<sample>.<fmt>` | **Per-file** plot per marker `M` — that file's cells on the shared coordinates and colour scale |

---

## Plot layout

Every plot shares the same geometry so panels are directly comparable:

- **Equal aspect (no distortion).** All plots use `coord_fixed(ratio = 1)` over a
  shared axis range, so UMAP1 and UMAP2 render at the same scale and the
  embedding shape is never stretched. The data panel is the same size in every
  plot regardless of its legend.
- **Legends never squish the panel.** Categorical figures (by-file,
  by-subpopulation) grow in **width** to fit their legend rather than shrinking
  the panel; the figure **height** stays at `--plot_height`.
- **Multi-column legends.** When a categorical legend has more labels than fit in
  one column at the plot height, it automatically wraps into 2, 3, … columns
  (≈23 labels per column at the default 8-inch height) and the figure widens
  accordingly — the panel size is unaffected.
- **Many-label palette.** The subpopulation plot uses a maximally-distinct
  palette (hues spread by the golden angle, cycled across several lightness and
  chroma levels) so dozens of labels stay distinguishable, where the default
  hue-only palette would not.

---

## Calling it from R

The logic lives in a shared module, so you can bypass the CLI and call it
directly (e.g. from another script or an interactive session):

```r
shared_dir <- "src/R_shared"
source(file.path(shared_dir, "io_utils.R"))        # make_outdir
source(file.path(shared_dir, "plot_utils.R"))      # save_plot_formats, create_safe_filename
source(file.path(shared_dir, "stats_utils.R"))     # run_umap
source(file.path(shared_dir, "flowframe_utils.R")) # explore_data
source("src/loading/discover_flowjo_workspaces.R") # workspace discovery/pairing
source("src/visualization/umap_from_flowjo.R")     # umap_from_flowjo

umap_from_flowjo(
  wsp_path   = "data/Killing_AML.wsp",   # or a folder of .wsp files
  fcs_dir    = "data",
  population = "Lymphocytes/Single Cells/Viable cells",
  output_dir = "results/umap_viable",
  subpop_mode = "furthest",
  output_formats = list(pdf = FALSE, png = FALSE, svg = TRUE)
)
```

`src/loading/discover_flowjo_workspaces.R` must be sourced too: it holds the
workspace discovery and FCS pairing (`find_workspace_files`, `index_fcs_dir`,
`resolve_workspace_specs`, …) that `umap_from_flowjo()` calls.

`umap_from_flowjo()` parameters mirror the CLI options (e.g. `markers` and
`subpop_select` take character vectors rather than comma-separated strings).
With `list_populations = TRUE` it prints the inspect blocks and returns
(invisibly) a list of `populations`, `samples` and `markers`; otherwise it
returns (invisibly) the `umap_coords`, the master `umap_df`, and the resolved
`output_dir`.

---

## Notes & troubleshooting

- **Population not found / ambiguous.** Run with `--list_populations` and copy a
  path verbatim. The script matches exact paths first, then a unique suffix.
- **Subpopulation plot is one colour.** With `--subpop_mode direct` or `deepest`
  the overlapping child gates collapse most cells to `(multiple)`. Use the
  default `furthest`, or `cd4cd8`, or an explicit `--subpop_select` set.
- **Wrong / too few samples loaded.** A workspace can contain several sample
  groups; only `--group` is loaded. The script logs the loaded sample list;
  set `--expected_samples 2` to get a warning if the count is unexpected, or
  change `--group`. With several workspaces the count is the **total** across
  them, and `--list_populations` prints each workspace's groups by name.
- **"could not build the workspace set."** One or more FCS files a workspace
  references could not be found under `--fcs_dir`, or the same file name exists
  in two places under it. The message names every workspace and every file; fix
  the FCS root (it is searched recursively) and re-run. Nothing is loaded until
  all of them resolve.
- **"--population could not be resolved in every workspace."** The gate exists
  in some workspaces but not all. Run `--list_populations` to see the gates they
  share — those are the only ones a run can use — and the block listing which
  workspace lacks which gate. Pooling cells from gates that are not the same
  gate is never done silently.
- **Markers silently missing.** Only channels with a non-empty marker name are
  kept, and only markers **common to all samples** are used. Scatter/`Time`
  channels are excluded by channel name. Excluded channels are logged per
  sample.
- **No transform is applied** by design — the embedding uses FlowJo's
  compensated/transformed working-space values. If a channel is still on a raw
  linear scale it can dominate the UMAP; restrict with `--markers` or try
  `--normalize zscore` if a plot looks off.
- **Non-finite values** (`NA`/`NaN`/`Inf`) in any marker abort the run before
  UMAP, naming the offending columns, to protect row alignment between
  coordinates, file labels, subpopulation labels, and colours.

---

## Tests

The workspace discovery and pairing logic is pure base R, so it can be checked
on any machine — no conda environment, no CytoML, no `.wsp` file needed:

```bash
Rscript test_data/test_umap_multi_wsp.R
```

It builds a fake directory tree and covers recursive `.wsp` discovery, the FCS
index, name matching (including missing and ambiguous files), group resolution
by index and by name, sample-label collision handling, population resolution
across differing trees, and subgate merging. What it cannot cover is CytoML's
own behaviour — that needs a real workspace on an analysis host.

---

## Related

- [FlowJo → UMAP Dash app](../../app/README.md) — a GUI over this script
- [`src/loading/discover_flowjo_workspaces.R`](../subprojects/loading.md#discover_flowjo_workspacesr) — the workspace discovery and FCS pairing functions
- [Getting Started](../getting-started.md) — environment setup
- [Pipeline D — FlowJo Comparison](../pipelines/pipeline-D-flowjo-comparison.md) — inter-rater gating comparison from `.wsp` files
- [Pipeline E — Preprocessed UMAP](../pipelines/pipeline-E-preprocessed-umap.md) — UMAP/clustering from preprocessed FCS (no gating)
- [Standalone Scripts](../standalone-scripts.md) — using module functions directly
