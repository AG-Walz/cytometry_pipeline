[< Back to Documentation Index](../README.md)

# UMAP from a FlowJo Workspace Population

`scripts/umap_from_flowjo.R` is a standalone script (no Nextflow) that turns a
FlowJo workspace plus its raw FCS files into a single shared UMAP, coloured by
every marker.

It:

1. Opens a FlowJo `.wsp` and builds a `GatingSet` over a directory of FCS files.
2. Extracts the cells of **one selected gated population** from each sample.
   The data is taken **as-is** from FlowJo's working space (already compensated
   and transformed), so **no further transform is applied**.
3. Concatenates the per-sample expression matrices.
4. Computes a **single shared UMAP embedding** over the concatenated data.
5. Plots that embedding coloured by sample and by every marker — both
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

**1. List the gating tree so you can pick a population:**

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp data/Killing_AML.wsp \
  --fcs_dir data \
  --list_populations
```

This prints every population path in the workspace and exits.

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

# Cap cells per sample AFTER gating (faster UMAP on large populations)
Rscript scripts/umap_from_flowjo.R --downsample 20000

# Also emit PDF and SVG, tweak the embedding
Rscript scripts/umap_from_flowjo.R --output_formats png,pdf,svg \
  --n_neighbors 30 --min_dist 0.3
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--wsp` | `data/Killing_AML.wsp` | Path to the FlowJo `.wsp` workspace file |
| `--fcs_dir` | `data` | Directory containing the FCS files referenced by the wsp |
| `--population` | `Lymphocytes/Single Cells/Viable cells` | Population path to extract (leading `/` optional; suffix-matched) |
| `--list_populations` | off | Print the population tree and exit (no analysis) |
| `--group` | `1` | FlowJo workspace sample-group index to load |
| `--expected_samples` | `0` | If `> 0`, warn when the loaded group does not contain this many samples |
| `--markers` | all common | Comma-separated subset of marker names to use/plot |
| `--output_dir` | `umap_flowjo` | Output directory for plots |
| `--output_formats` | `png` | Comma-separated formats: `png,pdf,svg` |
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
| `umap_<pop>_marker_<M>_ALL.<fmt>` | **Conjoint** plot per marker `M` — all cells, viridis-coloured by expression |
| `umap_<pop>_marker_<M>_file_<sample>.<fmt>` | **Per-file** plot per marker `M` — that file's cells on the shared coordinates and colour scale |

---

## Notes & troubleshooting

- **Population not found / ambiguous.** Run with `--list_populations` and copy a
  path verbatim. The script matches exact paths first, then a unique suffix.
- **Wrong / too few samples loaded.** A workspace can contain several sample
  groups; only `--group` is loaded. The script logs the loaded sample list;
  set `--expected_samples 2` to get a warning if the count is unexpected, or
  change `--group`.
- **Markers silently missing.** Only channels with a non-empty marker name are
  kept, and only markers **common to all samples** are used. Scatter/`Time`
  channels are excluded by channel name. Excluded channels are logged per
  sample.
- **No transform is applied** by design — the embedding uses FlowJo's
  compensated/transformed working-space values. If a channel is still on a raw
  linear scale it can dominate the UMAP; restrict with `--markers` if a plot
  looks off.
- **Non-finite values** (`NA`/`NaN`/`Inf`) in any marker abort the run before
  UMAP, naming the offending columns, to protect row alignment between
  coordinates, file labels, and colours.

---

## Related

- [Getting Started](../getting-started.md) — environment setup
- [Pipeline D — FlowJo Comparison](../pipelines/pipeline-D-flowjo-comparison.md) — inter-rater gating comparison from `.wsp` files
- [Pipeline E — Preprocessed UMAP](../pipelines/pipeline-E-preprocessed-umap.md) — UMAP/clustering from preprocessed FCS (no gating)
- [Standalone Scripts](../standalone-scripts.md) — using module functions directly
