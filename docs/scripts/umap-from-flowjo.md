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

With `--layout gates`, steps 3 to 5 change: the features are rebuilt from the
FlowJo gates themselves, so that every combination of gate calls gets its own
island, and the combination, split and conjoint marker plots are drawn in a
muted, axis-free style. See
[Gate-separated layout](#gate-separated-layout). A sample metadata table can
select which samples go in and split every plot by a column such as timepoint;
see [Sample metadata, filters and split plots](#sample-metadata-filters-and-split-plots).

---

## Prerequisites

The script runs inside the `cytometry_env` conda environment (it needs
`CytoML`, `flowWorkspace`, `flowCore`, `uwot`, `ggplot2`, `viridis`,
`optparse`). `--layout gates` also uses `RANN` (or `FNN`) for its separation
metrics, which are skipped with a message when neither is installed, and
`ggrastr`, when installed, to rasterise the points in SVG and PDF output so a
plot of many cells stays a small file. Both are declared in `environment.yml`.
If you have not set the environment up yet, follow the
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
evaluated at the root node and are the same set a real run would use. With
several workspaces two further blocks follow, indented four spaces; see
[Gating trees that differ](#gating-trees-that-differ).

This text is for people. [The Dash front-end](../../app/README.md) does not
parse it: it runs the same command with `--report_json <file>`, which writes
the whole load (workspaces, sample pairing, gating trees, marker panels and
what diverges between them) as JSON. `src/visualization/umap_load_report.R`
writes that file, stamped with a `schema` version, and
`app/services/rscript.py` reads it and refuses a major version it does not
know, so change them together. Of what a run prints, the app uses
the `::progress::` lines for its progress bar, the closing
`Done. Plots written to:` line to tell that a detached run from before exit
codes were recorded finished, and the last line of the log as the reason
when a run fails.

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

The report is printed in full to stderr and the error itself is a one-line
summary. That is deliberate: R truncates the error message it *prints* at
`getOption("warning.length")` — 1000 characters by default — which was enough to
swallow most of a 33-workspace report.

### Loading anyway when FCS files are missing

A collection of real workspaces routinely names files that were never copied to
the analysis host — negative controls, re-exports, a folder someone forgot. One
of them should not cost the other 32 workspaces, so `--allow_missing_samples`
drops the samples that cannot be resolved and loads the rest:

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp /data/workspaces --fcs_dir /data/fcs \
  --population "Lymphocytes/Single Cells/Viable cells" \
  --allow_missing_samples
```

Every dropped sample is named, before the slow load and again after the plots
are written, because the embedding is then missing cells the workspaces asked
for and that must not be lost in the scrollback:

```
WARNING: 3 sample(s) across 2 workspace(s) were NOT FOUND and have been dropped:

  MTB039_078  (11 of 12 sample(s) kept)
    Neg39_A10_A10_038.fcs
  MTB082  (0 of 21 sample(s) kept)  <- WORKSPACE DROPPED ENTIRELY
    MTB082_V0W0_Negativ.fcs, MTB082_V0W0_245024.fcs, ... and 19 more

1 workspace(s) contributed NO cells at all: MTB082
The embedding below is built WITHOUT those samples. Drop --allow_missing_samples to make this fatal again.
```

It relaxes **only** the missing-file case. An **ambiguous** name stays fatal:
a missing sample is a known absence, while an ambiguous one is a sample that
two files could be, and picking either would put the wrong cells in the
embedding. A workspace that will not open, or a bad `--group`, also stays fatal
— that is a broken input, not a gap in the data. If every workspace loses every
file the run still stops, since that means `--fcs_dir` is simply wrong.

In the app the same switch is **Advanced → Load anyway when FCS files are
missing**. It applies to loading the workspaces as well as to the run, because
the population tree cannot be printed until the FCS files resolve; while it is
on, the collapsed Advanced header carries a `missing FCS allowed` badge.

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

### Workspaces whose samples have different channels

FlowJo holds every sample as its own object and resolves each gate's dimension
against that sample's own parameter list, so one workspace can happily contain
acquisitions from two instruments. CytoML cannot: it materialises a sample group
as a single cytoset with **one** column schema, and refuses a group whose FCS
files do not all carry exactly the same `$PnN` set:

```
Found channel inconsistency across samples.
'Comp-Aqua LiveDead-A' is missing from MTB017_1_V4W279_neg.fcs_15850
```

That is not corruption. A cohort acquired over two years typically has:

* the same viability dye on a detector named `Aqua LiveDead-A` on one cytometer
  and `BV480-A` on another (identical `$PnS`, different configuration name);
* the scatter-width parameters `FSC-W`/`SSC-W` recorded under some acquisition
  templates and not others (`$PAR` 14 vs 12);
* the odd extra detector (`SP SSC-A`) from a third template.

Such a workspace used to be lost entirely — `refused, 0/0`. It is now split into
channel-homogeneous **parts**, one GatingSet each, and reported:

```
MTB017_1: its 20 sample(s) were acquired with 2 different channel sets, so CytoML
refused the workspace whole. Loading it as 2 parts, each with its own
channels; they rejoin in the UMAP on marker names.
    part 1 — 10 sample(s), e.g. MTB017_1_V4W279_neg.fcs
    part 2 — 10 sample(s), e.g. MTB017_1_V4W228_neg.fcs
        lacks: Aqua LiveDead-A, FSC-W, SSC-W
        adds:  BV480-A
```

Nothing is renamed and no parameter is dropped: each sample keeps the channel
names its gates were drawn against, which is the only way the gates stay the
gates FlowJo applied. The parts rejoin downstream because cells are pooled on
**marker**, not on channel — marker recovery (on unless `--no_recover_markers`) already maps
`Comp-Aqua LiveDead-A` and `Comp-BV480-A` onto the same marker. Counts, the
population map, plot legends and the load report all remain **per workspace**;
the parts show up only as a `2 channel sets` badge in the app's data-health
panel and in the log above.

A sample whose FCS file was never resolved by name has no channels to read and
joins the largest part. If a part still refuses after the split, only that
part's samples are lost — the rest of the workspace still loads.

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
| `combinations` | by its **combination of gate calls**, e.g. `CD4+ \| IFNy+ TNFa+` | The labels of the [gate-separated layout](#combination-labels), in its muted palette, with its combination, split and fold-change plots. Colours only: the embedding still uses the marker values. Markers layout only |
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

To give every combination of gate calls, rare ones included, its own island,
use the gate-separated layout below.

---

## Fit set of the markers layout

Fitting UMAP on every cell of a cohort of millions takes hours: the neighbour
graph and the initialisation grow faster than the cell count, the projection of
further cells onto a fitted model does not. `--fit_subset` chooses what the
markers layout fits on; every cell is then projected onto that model and
plotted.

| `--fit_subset` | Fitted on | Settings |
|----------------|-----------|----------|
| `random` *(default)* | `--umap_subsample` randomly chosen cells (default: 300,000 when the run has more, else every cell; `0` = every cell) | `--umap_subsample` |
| `supercell` | the mean marker values of [SuperCellCyto](https://github.com/phipsonlab/SuperCellCyto) supercells: each sample's cells grouped into supercells of about `--supercell_gamma` similar cells (never across samples), so the fit set is about the cell count / gamma and every phenotype keeps a share in proportion to its cells | `--supercell_gamma` (default 20, above 1), `--supercell_k_knn` (default 5), `--supercell_max` (default 100000) |

**Large cohorts.** Gamma is raised when the cohort would give more than
`--supercell_max` supercells: 8.1 million MTB cells at gamma 20 give 405,521,
so gamma becomes 82 and the fit about 99,000 supercells. The run says so. Two
reasons. The fit stays at the size the gate layout fits on (about a minute),
and it stays under the 300,000 points above which a spectral start becomes slow
and the run warns. A PCA start (`--umap_init pca`, never chosen automatically)
is scaled to the spectral start's 0-10
box: unscaled, on FlowJo values, it put the points about 100 units apart, the
optimisation never brought them together, and the 2026-10-14 MTB run came out
as specks over a range of 200. A fit whose layout still spans more than 100
units now warns. The projection of the cells onto the model runs on every
core (`--n_cores`); with a threaded embedding its result does not depend on
the core count.

Supercells need the R package SuperCellCyto: `install_extra_packages.R`
installs it from GitHub (bioconda has `bioconductor-supercellcyto` for R 4.5
only). Without it a supercell run stops before the workspaces load and says
how to install it. Samples with fewer than 2 x gamma cells are not grouped;
their cells join the fit set as they are. The gates layout has its own
balanced fit set and ignores `--fit_subset supercell` with a warning.

Example, the MTB cohort of the [gate-separated layout](#mtb-example) on the
plain marker values, every stimulated cell plotted:

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp /mnt/volume/patrick/UMAP/wsp --fcs_dir /mnt/volume/patrick/UMAP/fcs \
  --population "/Lymphocytes/Single Cells/Live cells" \
  --marker_aliases "INF=IFNy;TNF=TNFa;IFN=IFNy;INFy=IFNy" \
  --markers "CD4,CD8,IFNy,IL-2,TNFa,CD107a" \
  --sample_metadata /mnt/volume/patrick/UMAP/mtb_umap_metadata.csv \
  --filter stimulation=stimulated --split_by timepoint --annotate timepoint \
  --layout markers --subpop_mode combinations \
  --fit_subset supercell --supercell_gamma 20 \
  --n_neighbors 80 --min_dist 0.1 --n_cores 24 --output_formats png
```

### Background subtraction and an existing map

Cells cannot be subtracted, only their shares per tube. `--background
stimulation=unstimulated --pair_by patient,visit` pairs every sample with the
sample of the same patient and visit whose `stimulation` is `unstimulated`
(`background_pairs.csv` lists the pairs) and, in every legend, split plot,
fold-change map and frequency table, replaces a combination's percentage by
the sample's minus its background's, negatives set to 0. Only cytokine
responses are subtracted (combinations with a positive call below the
lineage, such as `CD4+ | IFNy+ TNFa+`). For `CD4+ | none`, `CD4- CD8-` and the
like, stimulated minus unstimulated is a near-zero difference of two large
numbers, and those are shown grey, "not a cytokine response". The background
samples drop out of the frequencies and the split plots; the whole-embedding
plots still show every cell. The embedding is not affected.

`--umap_model <run>/umap_model.uwot` skips the fit and projects every cell
onto that run's map, so a cohort with its unstimulated tubes lands on the map
of the stimulated-only run (same `--markers`, in the same order; the run
checks this against the model's `run_info.json`).

`scripts/compare_umap_runs.R --run_a A --run_b B --column timepoint
--output_dir C` then puts two such runs side by side: each level's
fold-change maps next to each other (`compare_<col>_<level>.png`), every
combination's fold change in A against B (`compare_<col>_scatter.png`) and
the numbers (`compare_<col>.csv`).

---

## Gate-separated layout

`--layout gates` builds the embedding from the gates under `--population`
instead of from the raw marker values, so that every combination of gate calls
(`CD4+ | IFNy+ TNFa+`, `CD8+ | none`, `CD4- CD8-`, ...) lands in its own
island. This section covers running it and reading the output; the method and
the reasoning behind each step are in
[Gate-anchored UMAP layout](../methods/gate-anchored-layout.md).

### Why

In the default layout a cytokine-positive cell differs from its negative
neighbours in one channel, often by a few units on FlowJo's 0 to 256 scale,
and shares every other channel with them, so it is easily drawn among them
instead of in an island of its own. Gates are also moved from visit to visit
as analysts follow the detector voltages, so a single global threshold would
disagree with the gater in many samples. The gate layout uses each sample's own
gates instead.

### How calls and cuts are read from the gates

For every sample, every non-Boolean gate below the population is checked for
the markers it decides:

* A rectangle bound on a marker channel that on its own excludes at least 0.5 %
  of the parent's cells binds. A gate that binds only at its lower bound is a
  positive decision (cells above the cut), one that binds only at its upper
  bound a negative decision. When both bounds bind, the gate is a band (a
  "CD4 dim" gate) and decides nothing on that marker, unless the smaller side
  excludes under 5 % of the parent: that is a far edge clipping a tail, and the
  larger side is the decision.
* Scatter, `Time` and viability-dye channels never decide, so the wide FSC-A
  range of a cytokine-versus-FSC gate drops out. Viability dyes are recognised
  by marker name (`--viability_pattern`).
* A polygon or other gate decides a marker when membership separates well along
  it (AUC of at least 0.8 or at most 0.2, with at least 10 members and 10
  non-members); its cut is the 1-D threshold that best reproduces membership.
* A gate never decides a marker that one of its ancestor gates below the
  population already decides. A CD4 x IFNy box under `CD4+` that clips the
  CD4-dim tail trims cells inside the side the `CD4+` gate chose; it is not a
  second CD4 threshold.

Each cell then gets the gater's call per marker. Only the deepest gates that
apply to a cell count, so a cytokine gate under `CD4+` wins over one hung
directly under the population for a CD4+ cell.

* Inside one or more of those gates: positive if any of them keeps cells above
  its cut, otherwise negative.
* Outside all of them: the opposite side when every one of them binds on this
  marker alone (outside the IFNy+ gate is IFNy-negative, whatever its FSC-A
  edge did). When a gate also binds on another marker the cell may be outside
  it because of that other marker, and when the gates have opposite
  polarities they disagree; in both cases the cell's own value against the
  median cut of those gates decides (`x >= cut` is positive).
* No gate applies: the marker is not decided (cytokines in CD4-CD8- cells, when
  the cytokine gates sit under CD4 and CD8 only).

Every decision goes to `gate_decisions.csv` with its gate, cut, polarity,
method and `n_axes`, the number of marker axes the gate binds on.

### The features

Each decided marker `m` becomes

```
f     = tanh((x - anchor) / (gate_softness * sigma)) + gap_m * s
s     = +1 for a positive call, -1 for a negative or undecided one
gap_m = gate_gap * gate_tier_ratio ^ (deepest tier - tier of m)
```

with `x` the FlowJo value. The `anchor` is the marker's median over the
sample's own cells with the same combination of calls, so in each tube the
CD4+ IFNy+ cells are centred on their own medians and the CD8+ cells on
theirs. That removes per-visit voltage shifts without depending on where the
gater drew the line in that tube, which can differ between tubes even when the
cells do not, or on which other populations share a side of the gate. It also
means a whole combination that is brighter in one tube than in another sits in
the same place in both; the marker plots keep the raw values. `sigma` is the
spread (MAD) of the negative cells around their anchor, pooled over samples
and at least 1 % of the 0 to 256 scale.

The tanh term always lies between -1 and 1, so two cells whose calls on `m`
differ as positive against negative or undecided are more than
`2 * (gap_m - 1)` apart on that coordinate, whatever the voltages of their
samples were. Negative and undecided share a side. With the defaults
(`--gate_gap 2`, `--gate_tier_ratio 2`) the cytokines of the MTB tree get a
jump of 2 and CD4 and CD8 one of 4, so the lineage split is the largest
distance and the islands of one lineage gather in one region.
`--gate_tier_ratio` must be 1 or more. `--gate_softness` sets how much of the
spread inside each island survives.

Cells that no gate decides on any marker (with a population above the live
gate, the cells that gate excluded) would otherwise sit on the negative side of
every marker, on top of the all-negative combination. They get one more
coordinate, `no_gate_decision`, at plus or minus the largest jump, and the
label `(no decision)`, and the log says how many there are.

By default the features are the markers that some gate below the population
decides and that every sample measures, which leaves out scatter, `Time`,
live/dead and anything gated only above the population. The marker plots cover
the same markers. `--markers` gives the exact list instead: decided markers
get the transform above, the others enter as a bounded robust z-score without a
jump, shaping the islands without splitting them.

### Combination labels

A marker's tier is how deep below the population its gates sit: in the MTB
tree under `Live cells`, CD4 and CD8 (the quadrant gates) are tier 1 and the
cytokines, gated under CD4 and CD8, are tier 2. A cell's label lists the
positive markers of each tier in panel order, tiers joined by ` | `, for
example `CD4+ | IFNy+ TNFa+`. A tier with no positive marker reads `none`,
except the first, which names its markers as negative (`CD4- CD8-`); tiers in
which nothing was decided for the cell are left out.

Labels come from marker names, not gate names, so gates spelled `IFNy` in one
workspace and `INFy` in another give one label once `--marker_aliases` makes
the marker names agree. A cell with no decision on any marker is labelled
`(no decision)`.

Combinations with fewer than `--min_combination_cells` cells in the whole run
are pooled as `other (rare)` everywhere: in the colours, the fit balancing,
supervision, the separation metrics, `combination_counts.csv` and the frequency
and change tables. They keep their own features, and their calls stay in the
`call_<marker>` columns of `cells.csv.gz`; lower `--min_combination_cells`
(`0` turns pooling off) to count them separately.

When a marker is decided for more than half of a lineage's cells (counted with
the sampling weights, so a share of cells rather than of tubes), the cells of
that lineage left without a decision on it, typically those of a sample whose
workspace lacks the gate, are treated as negative for the marker. The log
reports it, and `run_info.json` lists each fill per sample.

### Sampling and the balanced fit

1. **Counts first.** For every sample, the cells of each combination are
   counted before anything is sampled (`combination_counts.csv`). Every
   frequency in a legend or a CSV comes from these counts.
2. **Stratified sampling.** Per sample, at most `--cells_per_combination` cells
   of each call pattern are kept; rarer patterns are kept whole. Each kept cell
   carries a weight (cells of its pattern in that sample / cells kept), which
   restores the full counts wherever frequencies are rebuilt from the plotted
   cells. A sample in which no gate decides anything is sampled as one pattern,
   not dropped.
3. **Balanced fit.** The UMAP is fitted on a set in which no combination has
   more than `--umap_subsample` / (number of combinations) cells (100000 /
   number of combinations when it is `0`). Smaller combinations are fitted
   whole, and what they leave unused is not handed to the large ones, so the
   fit set can be much smaller than `--umap_subsample`, and a small run is
   balanced too. The other kept cells are projected onto the fitted embedding
   (`uwot::umap_transform`). This gives rare and common combinations comparable
   room; it also means that **island area does not reflect frequency**.

The projection's neighbour search runs on every plotting core (`--n_cores`).
It finds the same neighbours at any thread count, so the coordinates do not
depend on it.

Some options do not apply in this layout. `--downsample` and
`--normalize_data` are ignored with a warning (`--normalize_umap` still sets
the marker colours). The subgate plots of the markers layout
(`by_subpopulation`, `by_level<N>`) are not written, and `--subpop_mode`,
`--subpop_select`, `--subpop_levels_select` and `--no_subpop_levels` have no
effect.

### Supervision

`--supervise w`, with `w` between 0 and 1, fits the UMAP with the combination
label as target (uwot's `target_weight = w`); the default `0` fits without one.
It pulls each combination together more tightly. Projected cells are placed
from their features only, since uwot's projection takes no target.

### Separation metrics

`separation_metrics.csv` scores the embedding against the combination labels,
on up to 50000 cells balanced across labels, using each cell's 15 nearest
neighbours in the 2-D embedding. The first row, `(macro average)`, weights
every combination equally; one row per combination follows.

| Column | Meaning |
|--------|---------|
| `label` | The combination, or `(macro average)` |
| `cells_scored` | Cells of the label among the scored cells (not the run's total; see `combination_counts.csv` for those) |
| `knn_purity` | Share of a cell's 15 nearest neighbours that carry its label, averaged over the label's cells |
| `isolated` | Share of the label's cells whose 15 nearest neighbours all carry its label |
| `silhouette` | Distance to the nearest other label's centroid against distance to its own: 1 = far from every other island, 0 or below = overlapping |
| `group_mixing` | Within the label, how often neighbours come from a different workspace relative to chance (NA for a label found in one workspace only): 1 = as mixed as chance, 0 = every neighbour from the cell's own workspace |
| `sample_mixing` | The same over samples (NA for a label found in one sample only): 1 = the label's samples are mixed as by chance, 0 = each sample sits apart inside the island |

The metrics need the `RANN` or `FNN` package; without either the file is not
written. The two mixing columns are the most useful in practice: a low value
means the same combination from different workspaces or samples lands in
different places, which points at acquisition or gating differences, or a real
per-sample shift, rather than at the combination itself. The log lists every
island of at least 200 scored cells with either mixing below 0.5, and every
island with `isolated` below 0.9: islands apart in feature space can still be
drawn touching, and two touching colours are then not one population.

> **Read the islands as the gating, not as a finding.** The layout is
> gate-driven by construction: the gap puts cells with different calls apart
> before UMAP sees them, and `--supervise` adds the labels themselves. Separate
> islands, high purity and a high silhouette restate the gates; they are not
> independent evidence that the combinations differ biologically, and a
> misplaced gate gives a clean island of wrongly called cells. What the layout
> adds is the structure inside each island and how the composition of
> combinations changes between samples and metadata levels.
> [Reading the plots](../methods/gate-anchored-layout.md#reading-the-plots)
> says in plain terms what each figure shows and what it does not.

### Plot style

Combinations are coloured by their first-tier group: one muted hue family per
group (CD4+, CD8+, ...), the group's `| none` combination in its lightest tint
and darker shades for more positive markers, an all-negative first tier
(`CD4- CD8-`) in neutral greys, `(no decision)` in a near-black slate and
`other (rare)` in light grey. A group with a single combination gets an accent
colour between the families. The positive combinations of a family, the
accents, and every colour after the twelve reference colours when there are no
families are each the muted candidate furthest in CIE-Lab from every colour
already used, so even a large family stays tellable apart. Points are
small and translucent, sized and faded per combination to its density unless
`--point_size` and `--point_alpha` are given. There are no axes, grid or box,
only a short
UMAP1/UMAP2 arrow in the lower-left corner, and no text labels or leader lines
on the islands; the legend names the combinations.

The file, workspace, annotation and per-file marker plots keep their usual
style on the same coordinates. The conjoint marker plots
(`umap_<pop>_marker_<M>_ALL`) are redrawn in the clean style, still in viridis
and with the same colour limits as the per-file plots of that marker, and
their subtitle says whether `--normalize_umap` scaled the values. This also
happens with `--split_by` in the markers layout.

### MTB example

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp /mnt/volume/data/MTB_FACS/wsp \
  --fcs_dir /mnt/volume/data/MTB_FACS/fcs \
  --population "Lymphocytes/Single Cells/Live cells" \
  --marker_aliases "INFy=IFNy;IFN=IFNy;TNF=TNFa" \
  --layout gates \
  --sample_metadata mtb_umap_metadata.csv \
  --filter "stimulation=stimulated" \
  --split_by timepoint \
  --output_formats png,svg \
  --output_dir results/umap_mtb_gates
```

This pools the Live cells of all 20 workspaces. The files spell IFNg as `INFy`
or `IFN` and TNFa as `TNF` or `TNFa`, so the aliases give each one name. The
decision markers found under `Live cells` are CD4 and CD8 (the quadrant gates)
and CD107a, IL-2, IFNy and TNFa (the gates under CD4 and CD8). Only stimulated
tubes go in, and the embedding is drawn once per vaccination level, with the
fold-change maps drawn against `V0 pre-vaccination`, the first level in natural
order. The metadata file is described in the next section.

---

## Sample metadata, filters and split plots

One table describes the samples. Any of its columns can then colour the
embedding, select which samples go in, or split every plot into one per level.
Nothing here is specific to patient or timepoint: a column is a column.

### The metadata table

`--sample_metadata` takes a CSV, comma- or semicolon-separated (detected from
the header line). Every column is read as text, so `01` stays `01`.

* **The sample column.** `--sample_id_column` names it. Otherwise the first of
  `sample_id`, `sample`, `sampleID`, `sample_name`, `filename` and `file` that
  the table has is used, in that order of preference rather than the order of
  the columns (so a `sample` column wins over `filename`), and failing that the
  first column.
* **Matching.** A row matches a sample by FCS file name, ignoring case, any
  folder in front of the name, and the `.fcs` extension (with or without the
  `.fcs_0`-style suffix of FlowJo sample names). `MTB011_1_V3W315_neg.fcs` and
  `MTB011_1_V3W315_neg` are the same sample. The run's display name and the
  workspace name are tried after the file name, so a table keyed on workspace
  annotates every sample of that workspace.
* **Duplicates** give a warning, and the first row wins.
* **`Unassigned`.** A sample with no row is kept, with `Unassigned` in every
  column and a warning naming it. Blank cells and `NA`, `n/a`, `null`, `?` and
  `-` also read as `Unassigned`. It is drawn grey in the annotation plots, never
  becomes a split level, and fails every `=` filter clause that does not name it.

Write the `filename` column from the FCS folder rather than typing it. In MTB,
the `MTB16_1_...`, `MTB57_...` and `MTB59_...` files are not zero-padded
although their workspaces are, and a retyped `MTB016_1_...` leaves those samples
`Unassigned`. The app can write a template for the loaded samples: a `filename`
column, one column per named group of the pattern box, and one row per FCS file.
In the app the table can also be uploaded as TSV or XLSX (first sheet); it is
stored as CSV.

**MTB example** (`mtb_umap_metadata.csv`, one row per FCS file, 328 rows):

```
filename,patient,prior_peptide_vaccination,visit,timepoint,week,stimulation,cocktail
MTB011_1_V3W315_234216.fcs,MTB011,no,V3W315,V3 after 3rd vaccination,W315,stimulated,234216
MTB011_1_V3W315_neg.fcs,MTB011,no,V3W315,V3 after 3rd vaccination,W315,unstimulated,neg
MTB012_2_V4W70_215008.fcs,MTB012,yes,V4W70,V4 after 4th vaccination,W070,stimulated,215008
MTB59_1_V0W0_224065.fcs,MTB059,no,V0W0,V0 pre-vaccination,W000,stimulated,224065
MTB59_1_V1W5_214154.fcs,MTB059,no,V1W5,V1 after 1st vaccination,W005,stimulated,214154
```

| Column | Levels |
|--------|--------|
| `patient` | `MTB011` ... `MTB088`, zero-padded even where the file name is not |
| `prior_peptide_vaccination` | `no`, `yes` |
| `visit` | the `V..W..` token of the file name, e.g. `V3W315` |
| `timepoint` | `V0 pre-vaccination`, `V1 after 1st vaccination` ... `V4 after 4th vaccination` |
| `week` | `W000` ... `W315`, zero-padded |
| `stimulation` | `stimulated`, `unstimulated` (the neg tube) |
| `cocktail` | the peptide pool, `neg` for the unstimulated tube |

**Without a table**, `--sample_pattern` takes a regular expression with named
capture groups, applied to each sample's display name; the group names become
the columns:

```bash
--sample_pattern "(?<patient>MTB[0-9]+)_[0-9]_(?<visit>V[0-9]W[0-9]+)"
```

It is ignored when `--sample_metadata` is given. A pattern only knows what the
names contain: in MTB it would make `MTB016` and `MTB16` two patients.

`--annotate` picks which annotation columns get plots, in that order (default:
all). Each chosen column gives `umap_<pop>_annot_<col>_ALL`, all cells coloured
by the column, and one `umap_<pop>_annot_<col>_level_<level>` per level, that
level's cells in one colour over all cells in grey. The `--split_by` columns
are added to the list, since a split needs their values per cell; `--filter`
reads every column of the metadata, whatever `--annotate` selects.

### Filters

`--filter` keeps samples by their metadata before anything is extracted, so a
filtered-out sample is never read. It works in both layouts and filters on the
annotation columns, so it needs `--sample_metadata` or `--sample_pattern`.

```bash
--filter "stimulation=stimulated"                 # keep samples with this value
--filter "stimulation!=unstimulated"              # drop samples with this value
--filter "stimulation=stimulated;timepoint!=V4 after 4th vaccination"
```

* Clauses are separated by `;` and must all hold.
* `column=a,b` keeps samples whose value is `a` or `b`; `column!=a,b` drops them.
* Values are compared exactly, case included, after trimming spaces around them.
  A value containing `,` or `;` cannot be named.
* A column the metadata does not have is an error that lists the columns it has.

Quote the whole value: an unquoted `;` ends the shell command (see
[Quoting on the command line](#quoting-on-the-command-line)).

### Split plots

`--split_by col1,col2` draws, for each column, the same embedding once per
level, that level's cells in colour over a light grey silhouette of all cells.
It works in both layouts. In the `gates` layout cells are coloured by gate
combination and every percentage comes from the full counts. In the `markers`
layout they are coloured by their subgate label (`--subpop_mode`); with
`--subpop_mode combinations` by gate combination, the percentages again from
the full counts. With the other modes the percentages are shares of the
plotted cells, so of the downsampled ones with `--downsample`.

* **Level order** is natural: runs of digits compare as numbers, so
  `W5 < W12 < W110` and `V0 ... < V1 ... < V4 ...`. `--level_order` overrides
  it per column, `"timepoint=V0 pre-vaccination,V2 after 2nd vaccination;week=W000,W005"`:
  the listed levels come first in that order, the rest follow in natural order,
  and a listed level that does not occur gives a warning.
* **The reference level** for the fold-change maps is the first level in that
  order. `--split_reference "timepoint=V0 pre-vaccination"` picks another.

Per column the run writes:

* `umap_<pop>_split_<col>_<level>`: that level's cells coloured by combination
  over the grey silhouette of all cells. The legend gives each combination's
  mean per-sample percentage of the population in that level.
* `umap_<pop>_split_<col>_panel`: all levels side by side, with shared colours
  and one legend.
* `umap_<pop>_diff_<col>_<level>`: every island coloured by how the frequency
  of its combination changed at this level against the reference, as
  `log2((mean % at level + p) / (mean % at reference + p))`. A mean % is the
  combination's percentage of the population averaged over the level's
  samples, so every **sample** counts equally. `p` is one cell in the median
  sample of the two levels (100 / its population cells, in %), which keeps a
  change of well under one cell per sample pale. Colours stop at 8-fold either
  way, and a combination absent at both levels is grey. There is no map for
  the reference itself.
* `combination_frequencies_<col>.csv`: per level and combination, the mean and
  SD over the level's samples of the combination's percentage of the
  population, and the number of samples.
* `combination_change_<col>.csv`: the numbers behind the fold-change maps, one
  row per level other than the reference and combination: `column`, `level`,
  `reference`, `combination`, `mean_pct_level`, `mean_pct_reference`,
  `log2_fold_change` (NA when the combination is absent at both levels),
  `n_samples_level`, `n_samples_reference` and `pseudo_pct` (the `p` above).

A fold-change map is descriptive: it compares two means over samples, without
a test or a confidence interval, and the samples of two levels need not come
from the same patients. In MTB, `V0 pre-vaccination` has 8 stimulated tubes
from 4 patients, none of them previously peptide-vaccinated, while
`V1 after 1st vaccination` has 60 from 12,
so a change against V0 mixes time with differences between patients, and one
patient can drive the change of a rare combination. For tests, start from the
per-sample counts in `combination_counts.csv`.

### Polyfunctionality maps

A run whose cells are labelled by gate combination (`--subpop_mode
combinations` or `--layout gates`) also draws, unless `--no_polyfunctionality`,
the maps of the MTB cohort (`scripts/plot_polyfunctionality.R
--size_by_polyfunctionality`): every cell coloured by how many of the gates
below its lineage it is in, `CD4+ | IFNy+ TNFa+` counting 2. The gates are the
last tier of the run's gate decisions (IFNy, IL-2, TNFa and CD107a in MTB), so
the scale runs 0 to 4 even when no cell is positive for all of them. The
colours are the five MTB blues, light (0) to dark (4; the same ramp for
another number of gates), and a cell is drawn `1 + 2 * n / k` times the base
size and on top of those with fewer, so 1x to 3x for four gates. The base size
is `--point_size`, else 0.4 or the automatic size if that is smaller; the
alpha is `--point_alpha`, else 0.85. Cells whose label has no tier below the
lineage, and `other (rare)`, have no count and are not drawn.

Written to `polyfunctionality/`, on the frame of the whole embedding so the
maps line up:

* `umap_<pop>_polyfunctionality`: every cell of the run. The legend gives the
  cells and percentage with each count, of the cells shown (weighted in the
  gates layout).
* `polyfunctionality/<col>/umap_<pop>_polyfunctionality_<col>_<level>`: one
  per level of each `--split_by` column, in its level order.
* `polyfunctionality_per_sample.csv`: per sample and count, `sample`,
  `n_positive`, `cells`, `percent` (of the sample's counted cells) and the
  annotation columns, for tests or plots of your own.

Unlike the MTB script, the run does not cross columns (patient x timepoint x
stimulation) or draw equal cells per map: use `scripts/plot_polyfunctionality.R`
or `scripts/plot_polyfunctionality_timepoints.R` on the finished run for those.

---

## Quoting on the command line

Every option is passed through your shell first, so the rule is a shell rule,
not a rule of this script: **quote a value if it contains a space, a `;`, a
`*`, a `(`, a `$`, or a `+`.** Nothing is harmed by quoting a value that did not
need it, so when in doubt, quote.

In practice that means:

```bash
# Spaces — the usual case. Gate paths and group names almost always have them.
--population "Lymphocytes/Single Cells/Viable cells"
--group "All Samples"

# Semicolons — these separate entries, and an unquoted ';' ENDS THE COMMAND.
# Without the quotes the shell runs everything after the first ';' as a
# separate command, and the script silently receives only the first pair.
--marker_aliases "TNF=TNFa;IFN=INFy"
--population_map "MTB016=/Lymphocytes/Live cells;MTB023=/Lymphocytes/Viable"
--filter "stimulation=stimulated;timepoint!=V4 after 4th vaccination"
--split_reference "timepoint=V0 pre-vaccination"

# Parentheses and '+' — real FlowJo workspaces contain both.
--workspaces "MTB061_ (2),MTB074"
--subpop_select "CD4+/IFN,CD8+/IFN"

# No spaces or shell characters — quotes optional, and usually omitted.
--n_neighbors 15
--match_mode content
--output_formats png,svg
```

Two things that are NOT about quoting and catch people anyway:

* **Commas separate list items; do not put spaces after them.** `--markers
  "CD4, CD8"` asks for a marker literally called `" CD8"`. Write
  `--markers CD4,CD8`.
* **A gate name containing a comma** cannot go in a comma-separated list. That
  is why `--population_map` and `--marker_aliases` use `;` between entries — a
  FlowJo gate may contain a comma but, in practice, never a semicolon.

Inside a params YAML (`config/pipeline_F_params.yaml`) there is no shell, so
these rules do not apply — quote only what YAML itself requires:

```yaml
population: Lymphocytes/Single Cells/Viable cells   # no quotes needed
marker_aliases: "TNF=TNFa;IFN=INFy"                 # quoted: starts a YAML flow
group: "1"                                          # quoted: keep it a string
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--wsp` | `data/Killing_AML.wsp` | A `.wsp` file, **or a folder** searched recursively for `.wsp` — see [Many workspaces at once](#many-workspaces-at-once) |
| `--fcs_dir` | `data` | Folder holding the FCS files, searched **recursively**; each workspace is paired with its own files by name |
| `--population` | `Lymphocytes/Single Cells/Viable cells` | Population path to extract (leading `/` optional; suffix-matched). Resolved independently in each workspace |
| `--list_populations` | off | Load the workspace(s), print their samples, population tree, common markers and workspace/divergence blocks, then exit (no analysis) |
| `--group` | `1` | FlowJo sample group to load: an index, **or a group name** such as `"All Samples"` |
| `--allow_missing_samples` | off | Drop samples whose FCS file is not found instead of failing, reporting each one (see [above](#loading-anyway-when-fcs-files-are-missing)). An *ambiguous* name is still fatal |
| `--expected_samples` | `0` | If `> 0`, warn when the run does not total this many samples across all workspaces |
| `--report_json` | none | With `--list_populations`, also write the whole load to this path as JSON. This is what the Dash app reads |
| `--sample_metadata` | none | CSV keyed on sample (see [The metadata table](#the-metadata-table)). Every other column becomes a per-sample annotation that can be plotted, filtered on and split by |
| `--sample_id_column` | see description | The `--sample_metadata` column holding the sample identifier. Default: the first of `sample_id`, `sample`, `sampleID`, `sample_name`, `filename`, `file` present, else the first column |
| `--sample_pattern` | none | Alternative to `--sample_metadata`: a regex with named capture groups applied to the sample name, e.g. `"(?<patient>MTB[0-9]+)_[0-9]_(?<visit>V[0-9]W[0-9]+)"`. Group names become the columns. Ignored when `--sample_metadata` is given |
| `--annotate` | all columns | Comma-separated subset of annotation columns to plot, in this order. The `--split_by` columns are always added |
| `--filter` | none | Keep only matching samples, before anything is extracted: `"col=a,b;col2!=c"` (see [Filters](#filters)). Works in both layouts and on any metadata column, whatever `--annotate` selects |
| `--split_by` | none | Comma-separated metadata columns. Per column: one plot per level, a panel of all levels, fold-change maps against a reference level, `combination_frequencies_<col>.csv` and `combination_change_<col>.csv` (see [Split plots](#split-plots)). Works in both layouts |
| `--level_order` | natural sort | Display order of levels: `"col=lvA,lvB;col2=..."`. Levels not listed follow in natural order |
| `--split_reference` | first level | Reference level per column for the fold-change maps: `"col=level;col2=level"` |
| `--markers` | all common | Comma-separated subset of marker names to use/plot. With `--layout gates`: exactly these markers become features; gate-decided ones get the gap transform, the others enter without a jump |
| `--subpop_mode` | `furthest` | Subgate set for the subpopulation plot: `furthest`, `cd4cd8`, `direct`, `deepest`, `combinations`, `none` (see [above](#subpopulation-colouring)). No effect with `--layout gates` |
| `--subpop_select` | — | Comma-separated explicit gate list (relative path or bare name) used as the colour partition, overriding `--subpop_mode`. No effect with `--layout gates` |
| `--no_polyfunctionality` | off | Do not draw the [polyfunctionality maps](#polyfunctionality-maps), which a run labelled by gate combination draws by default |
| `--gate_separation` | `0` | Markers layout: a middle ground between one continuum and the gate islands. Each cell's gate label (the subpopulation colouring, `--subpop_mode`) is added to the embedding as a column worth this many typical marker spreads (median MAD), so cells of different gates are pushed apart while the marker structure inside each stays. `0.5`-`1`: gated populations become compact clusters next to their parents; `2` or more: islands. Colours and saved values are unaffected. Not with `--layout gates` or `--umap_model` |
| `--normalize_data` | `none` | Per-marker normalisation of the values the **embedding** is computed from: `none`, `zscore`, `minmax`. Changing it moves the points; it is a different analysis. Ignored with `--layout gates`, with a warning |
| `--normalize_umap` | `same` | Scale of the **marker colours**. `same`: the scale the embedding was computed on, so a z-scored embedding is coloured by z-scores and a plain one by FlowJo values. `none` (FlowJo values), `zscore`, `minmax` choose one explicitly. Recolours only; nothing moves. With `--replot` it recolours a finished run (default: the run's own scale) |
| `--normalize` | none | Deprecated alias of `--normalize_data`; passing both is an error |
| `--output_dir` | `umap_flowjo` | Output directory for plots |
| `--output_formats` | `png` | Comma-separated formats: `png,pdf,svg`. SVG/PDF hold one vector point per cell (hundreds of MB per plot on millions of cells) |
| `--downsample` | `0` | Max cells per sample after gating (`0` = no downsampling). Ignored with `--layout gates`, which samples per combination instead. With `--subpop_mode combinations` a sample also keeps one cell of every gate-call pattern the random draw missed, so a few more than this may be kept and every pattern's full count stays labelled |
| `--layout` | `markers` | `markers`: embed the FlowJo values of the markers, as described above. `gates`: gate-anchored features, combination labels, stratified sampling, balanced fit and the muted plot style (see [Gate-separated layout](#gate-separated-layout)) |
| `--gate_gap` | `2` | `gates` layout: jump added at each gate cut of the deepest tier, in feature units (shallower tiers get it times `--gate_tier_ratio` per tier). The within-side term is bounded by 1, so the two sides of a gate are more than `2 * (gate_gap - 1)` apart |
| `--gate_tier_ratio` | `2` | `gates` layout: how much larger the jump is one tier up the gating tree; with 2, CD4 and CD8 get twice the cytokines' jump and the islands of one lineage gather in one region. Must be 1 or more |
| `--gate_softness` | `4` | `gates` layout: tanh softness, in units of the marker's negative-population MAD around its anchor |
| `--viability_pattern` | built-in list | `gates` layout: case-insensitive regular expression of marker names that are viability dyes and never decide a combination. The built-in list covers names such as `live`, `dead`, `viab`, `zombie`, `fixable`, `ALD`, `7-AAD`, `DAPI`, `PI`, `FVD`, `FVS`, `SYTOX`, `DRAQ7`, `Aqua`, `L/D` and `LD`; `""` turns it off |
| `--supervise` | `0` | `gates` layout: uwot `target_weight` with the combination label as target; `0` = unsupervised |
| `--cells_per_combination` | `2000` | `gates` layout: per sample, the most cells kept of each call pattern; rarer ones are kept whole |
| `--min_combination_cells` | `50` | `gates` layout: combinations with fewer cells in total (weighted back to the full counts) are pooled as `other (rare)`, in the plots and in every count and frequency table; `0` turns pooling off |
| `--n_neighbors` | `15` | UMAP `n_neighbors` |
| `--min_dist` | `0.1`; `0.25` with `--layout gates` | UMAP `min_dist` |
| `--metric` | `euclidean` | UMAP distance metric |
| `--umap_init` | `auto` | Initial layout: `auto` (= `spectral`), `spectral`, `pca`, `random`. Spectral is the best start but single-threaded; it is never swapped automatically. Past 300k fitted cells the run warns at once and suggests fitting on fewer cells (`--umap_subsample`) or passing `pca`. A PCA start is scaled to a 0-10 range; unscaled, it does not converge on FlowJo values |
| `--umap_subsample` | `-1` | Fit the embedding on this many randomly chosen cells and project the rest onto it; `0` = fit on every cell; `-1` = 300,000 when the run has more cells, else every cell (the fit then keeps the spectral start, see `--umap_init`). With `--layout gates`: the fit budget, each combination contributing at most this / (number of combinations) cells, and `0` or `-1` means 100000 (see [Sampling and the balanced fit](#sampling-and-the-balanced-fit)). Ignored with `--fit_subset supercell` |
| `--background` | — | `column=level` of the background samples, e.g. `stimulation=unstimulated`; frequencies become sample minus background for cytokine responses (see [Background subtraction](#background-subtraction-and-an-existing-map)). Needs `--pair_by` |
| `--pair_by` | — | Comma-separated metadata columns a sample shares with its background sample, e.g. `patient,visit` |
| `--umap_model` | — | A saved `umap_model.uwot` (markers layout): nothing is fitted, every cell is projected onto that map |
| `--no_resume` | off | Recompute every step instead of reusing the checkpoints an earlier run into the same `--output_dir` saved (see [Restarting a run](#restarting-a-run)) |
| `--replot` | — | Run folder to redraw from its saved files (see [Output](#output)); nothing is read or recomputed |
| `--fit_subset` | `random` | Markers layout: fit on `random` cells (`--umap_subsample`) or on SuperCellCyto `supercell` means, then project every cell (see [Fit set of the markers layout](#fit-set-of-the-markers-layout)) |
| `--supercell_gamma` | `20` | Cells per supercell with `--fit_subset supercell`; above 1 |
| `--supercell_k_knn` | `5` | Neighbours of SuperCellCyto's within-sample kNN graph with `--fit_subset supercell` |
| `--supercell_max` | `100000` | Most supercells UMAP is fitted on with `--fit_subset supercell`; gamma is raised to cells / this on a larger cohort (see [Fit set of the markers layout](#fit-set-of-the-markers-layout)) |
| `--n_cores` | `0` | Worker processes for the per-file plots: `0` = automatic, leaving one core free; `1` = no parallelism. Capped at the cgroup CPU quota. An explicit value above 1 also runs the embedding multi-threaded, which is reproducible but gives different coordinates from a single-threaded run. With `--layout gates` the projection of the non-fitted cells always uses the plotting cores, which changes its speed but not the coordinates |
| `--seed` | `42` | Random seed (UMAP, downsampling, and the sampling and fit set of `--layout gates`) |
| `--point_size` | `0.4`; automatic in the styled plots | ggplot point size. Unset, the plots in the usual style (file, workspace, subpopulation, annotation and marker plots) use 0.4, and the combination, split, panel and styled conjoint marker plots size the points of each combination to its density. A value applies to every plot |
| `--point_alpha` | `0.6`; automatic in the styled plots | ggplot point alpha, set the same way as `--point_size` |
| `--plot_width` | `10` | Plot width (inches) |
| `--plot_height` | `8` | Plot height (inches) |
| `--src_dir` | auto | Path to the repo `src` directory (auto-detected from the script location) |

Run `Rscript scripts/umap_from_flowjo.R --help` for the same list at the
command line.

---

## Output

All files are written to `--output_dir` (created if needed), one subfolder
per kind of output, so a figure can be found without scrolling past two
thousand per-sample plots:

```
<output_dir>/
  run_info.json, plot_index.json    what was run, and every plot written
  overview/                         by combination, file, workspace, subgate
  markers/                          one map per marker, all cells
  markers/per_sample/<marker>/      the per-file marker plots
  metadata/<column>/                coloured by an annotation, one plot per level
  splits/<column>/                  per-level plots and the side-by-side panel
  changes/<column>/                 fold-change maps, combination_frequencies_ and
                                    combination_change_<column>.csv
  polyfunctionality/[<column>/]     positive gates per cell: the whole run, each
                                    split level, polyfunctionality_per_sample.csv
  data/                             cells.csv.gz, supercells.csv.gz,
                                    combination_counts.csv, gate_decisions.csv,
                                    separation_metrics.csv, background_pairs.csv,
                                    sample_metadata.csv
  model/                            umap_model.uwot
  checkpoints/                      for resuming (see Restarting a run)
```

File names are the same as before and say what a plot is wherever it sits.
Runs before 2026-10-15 wrote everything flat into `--output_dir`;
`Rscript scripts/organise_umap_run.R --run <folder> [--dry_run]` moves such a
run into this layout (and rewrites its `plot_index.json`), and `--replot`, the
app and `scripts/compare_umap_runs.R` read both. Plots have no background
grid. File names are sanitised from the population and marker names.

| File | Content |
|------|---------|
| `umap_<pop>_by_file.<fmt>` | Shared embedding, points coloured categorically by source file |
| `umap_<pop>_by_workspace.<fmt>` | Shared embedding, points coloured by source **workspace**. Only written when more than one workspace was loaded — see [Many workspaces at once](#many-workspaces-at-once) |
| `umap_<pop>_by_subpopulation.<fmt>` | Shared embedding, points coloured by FlowJo subgate (per `--subpop_mode`); omitted when `--subpop_mode none` finds no subgates, and not written with `--layout gates`, which colours by combination instead |
| `umap_<pop>_marker_<M>_ALL.<fmt>` | **Conjoint** plot per marker `M` — all cells, viridis-coloured by expression. With `--layout gates` or `--split_by` it is redrawn in the clean style on the same colour scale (see [Plot style](#plot-style)) |
| `umap_<pop>_marker_<M>_file_<sample>.<fmt>` | **Per-file** plot per marker `M` — that file's cells on the shared coordinates and colour scale |
| `umap_<pop>_annot_<col>_ALL.<fmt>` | All cells coloured by annotation column `col`, with `Unassigned` in grey. Written per column chosen by `--annotate` when `--sample_metadata` or `--sample_pattern` is given |
| `umap_<pop>_annot_<col>_level_<level>.<fmt>` | One per level of `col`: that level's cells in one colour over all cells in grey |

With `--layout gates`, or with `--split_by` in either layout, the run also
writes the files below. Those marked *gates* are written only in the `gates`
layout; the split, change and frequency files only with `--split_by`.

| File | Content |
|------|---------|
| `umap_<pop>_by_combination.<fmt>` | All plotted cells coloured by combination label (by subgate label in the markers layout), in the muted style (see [Gate-separated layout](#gate-separated-layout)) |
| `umap_<pop>_split_<col>_<level>.<fmt>` | One per level of each `--split_by` column: that level's cells coloured by combination over a light grey silhouette of all cells. The legend gives each combination's mean per-sample percentage in that level |
| `umap_<pop>_split_<col>_panel.<fmt>` | All levels of `col` side by side, shared colours, one legend |
| `umap_<pop>_diff_<col>_<level>.<fmt>` | Each island coloured by the log2 fold change of its combination's mean per-sample percentage at the level against the reference level, capped at 8-fold, grey where the combination is absent at both, with a grey outline of the embedding (see [Split plots](#split-plots)). Not drawn for the reference itself |
| `combination_frequencies_<col>.csv` | Per `--split_by` column: `column`, `level`, `combination`, `mean_pct` (mean over the level's samples of 100 x cells / population cells), `sd_pct`, `n_samples`. From the full counts in the gates layout, from the plotted cells in the markers layout |
| `combination_change_<col>.csv` | Per `--split_by` column, the numbers behind the fold-change maps: `column`, `level`, `reference`, `combination`, `mean_pct_level`, `mean_pct_reference`, `log2_fold_change`, `n_samples_level`, `n_samples_reference`, `pseudo_pct` |
| `combination_counts.csv` | *gates*. Per sample and combination, from all cells before sampling and zero-filled: `sample`, `combination`, `cells`, `population_cells` |
| `gate_decisions.csv` | *gates*. Per sample, every decision read from the gates: `gate`, `parent`, `depth`, `channel`, `marker`, `polarity`, `cut`, `method`, `excluded` (share of the parent on the excluded side), `n_parent`, `n_member`, `n_axes` (marker axes the gate binds on, bands included, markers an ancestor gate decides left out), `sample`, `workspace` |
| `separation_metrics.csv` | *gates*. How well the embedding separates the combinations: `label`, `cells_scored`, `knn_purity`, `isolated`, `silhouette`, `group_mixing`, `sample_mixing` (see [Separation metrics](#separation-metrics)) |
| `cells.csv.gz` | *gates*. Every plotted cell: `UMAP1`, `UMAP2`, `sample`, `workspace`, `combination`, one `call_<marker>` column per decided marker (1, 0 or NA), `weight`, `fit` (in the fit set or projected) and one `annot_<col>` column per annotation |
| `run_info.json` | *gates*. Creation time, command line, git commit and dirty state, R and package versions, seed, thread count, the parameters that shape the result (population, gate, sampling and UMAP settings, filter, splits, metadata source, markers and aliases), the features, feature scales, a summary of the decisions and fills, cell, fit-set and sample counts, samples left out by the filter, and the weighted total of each combination |
| `sample_metadata.csv` | *gates*. Copy of the `--sample_metadata` table (not written with `--sample_pattern`) |
| `plot_index.json` | Every plot written, with its kind (`combination`, `split`, `split_panel`, `diff`, `marker`, ...), column, level and reference. The app's gallery sorts plots by it |

Every markers-layout run also keeps what it computed, written before any plot
so a plotting failure cannot cost it, and enough to redraw or re-analyse the
embedding without reading a workspace or recomputing anything:

| File | Content |
|------|---------|
| `cells.csv.gz` | *markers*. Every cell: `UMAP1`, `UMAP2`, `sample`, `workspace`, `label` (subgate or, with `--subpop_mode combinations`, gate combination), `supercell` (with `--fit_subset supercell`: its supercell's number within the sample, NA for a sample too small to group), one `annot_<col>` column per annotation and one `value_<marker>` column per embedded marker: the FlowJo values, before any normalisation, so `--replot --normalize_umap` can show any scale exactly (runs before 2026-11-03 saved the values UMAP saw; `run_info.json`'s `values_scale` says which) |
| `supercells.csv.gz` | *markers*, `--fit_subset supercell`. The fit set, one row per supercell: `sample`, `supercell`, `cells` (its size) and the mean of each marker. Written as soon as SuperCellCyto finishes |
| `umap_model.uwot` | The fitted UMAP model, whenever the embedding was fitted on a subset and the rest projected (supercells, `--umap_subsample`, the gates layout's balanced fit). `uwot::load_uwot()` reads it and `uwot::umap_transform()` places new cells on the same map. Written right after the embedding |
| `run_info.json` | *markers*. As in the gates layout, with the fit set (`fit_subset`, gamma and k_knn or `umap_subsample`), the embedded markers, the supercell count and how long SuperCellCyto took, and the cell count per label |
| `gate_decisions.csv`, `combination_counts.csv`, `sample_metadata.csv` | *markers*, as in the gates layout: the decisions and full counts with `--subpop_mode combinations`, the metadata copy with `--sample_metadata` |

### Settings record: what a run used, and how to repeat it

Before it reads anything, every run writes what it was asked to do:

| File | What it holds |
|---|---|
| `_run/command.sh` | The exact command, with **every** option written out, defaults included, as `--option='value'`. `bash _run/command.sh` repeats the run. |
| `_run/settings.json` | Every option's effective value, which ones were typed (`options_given`), the git commit, branch and whether the checkout had local changes, R and package versions (uwot, CytoML, flowWorkspace, ...), host, user, working directory, `TMPDIR`, the inputs (the workspace's size, date and md5; the FCS folder's file count, with its full listing in `fcs_dir_files.tsv`), and when the run ends: `status` (`completed` / `failed`), `minutes` and the `error`. |
| `_run/settings/<YYYYmmdd-HHMMSS>/` | The same files for every run into this folder, kept: a re-run does not overwrite the record of the one before. Holds `uncommitted.diff` when the code had local changes. |

Defaults are written out because defaults change: on 2026-11-03 the fit set,
the output format and the UMAP start all did. A command that only lists what
was typed would repeat a different run today. A run whose `status` is still
`running` long after it stopped was killed before R could record the outcome
(the app's `_run/exit_code` has the exit code). `run_info.json` points to the
record (`"settings"`); the app shows it under **Settings used** in the results.

To repeat a run into another folder, edit `--output_dir` in a copy of
`command.sh`; to repeat it on the same code, check out the commit its header
names first.

### Restarting a run

Every step saves a checkpoint in `checkpoints/` inside the run folder as soon
as it is done: `extract.rds` (the cells, labels and values read from the
workspaces), `supercells.rds`, `model.rds` (with `umap_model.uwot`) and
`embedding.rds` (the coordinates). Each carries a key of the settings it
depends on, and of the inputs' file sizes and modification times. Running the
same command into the same folder again reuses every step whose key still
matches and continues from the first that does not:

| Rerun after | Recomputed |
|-------------|------------|
| a crash, cancel or power cut | from the step that was running |
| changing only plot settings (point size, formats, ...) | the plots only |
| changing `--min_dist`, `--n_neighbors`, `--seed`, ... | the fit, projection and plots |
| changing `--supercell_gamma` or `--supercell_max` | from the supercells on |
| changing the markers, filter, population or an input file | everything |

The projection is seeded per chunk, so a run resumed from a saved model gives
the coordinates the uninterrupted run would have. On the two-workspace MTB
test a rerun with the same settings took 66 s instead of 185 s and wrote a
byte-identical `cells.csv.gz`, as did a run "crashed" after its model was
saved. `--no_resume` recomputes everything. The checkpoints of an 8-million
cell run take a few GB; delete `checkpoints/` once a run is final.

To redraw every plot of a run from these files, for example with other point
sizes, plot dimensions or more cores, without reading a workspace or
recomputing anything:

```bash
Rscript scripts/umap_from_flowjo.R --replot /path/to/run_folder --n_cores 24
```

`--replot` reads `cells.csv.gz`, `run_info.json` and `combination_counts.csv`,
takes the population, splits, level order and colouring from `run_info.json`,
and writes the plots into the run folder, replacing the old ones. Only
`--n_cores`, `--point_size`, `--point_alpha`, `--plot_width`, `--plot_height`
and `--output_formats` apply. On the two-workspace MTB test it gave the same
132 PNG files, byte for byte, as the run that computed them, in 51 s.

The overview plots (by file, by workspace, per annotation level, per marker)
are drawn in parallel, at most 12 at a time; each one draws every cell, about
1.3 minutes at 8 million cells.

To reuse a run in R:

```r
cells <- data.table::fread(cmd = "zcat cells.csv.gz")    # or read.csv(gzfile(...))
model <- uwot::load_uwot("umap_model.uwot")
x <- as.matrix(cells[, grep("^value_", names(cells)), with = FALSE])
colnames(x) <- sub("^value_", "", colnames(x))
new_xy <- uwot::umap_transform(x[1:1000, ], model)       # same map, no refit
```

Names go through the same sanitising as every other plot: characters other than
letters, digits, `.`, `_` and `-` become `_`, runs of `_` or `-` shrink to one,
and both are stripped from the ends, so the level `V0 pre-vaccination`
gives `umap_Lymphocytes_Single_Cells_Live_cells_split_timepoint_V0_pre-vaccination.png`.
Names are cut to 100 characters, and a name that would collide with one
already written gets `_2`, `_3`, and so on.

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
source(file.path(shared_dir, "memory_utils.R"))    # resolve_cores, parallel_lapply
source(file.path(shared_dir, "stats_utils.R"))     # run_umap
source(file.path(shared_dir, "dimred_dispatch.R")) # run_dimred
source(file.path(shared_dir, "flowframe_utils.R")) # explore_data
source("src/loading/flowjo_workspace.R")           # workspace discovery/pairing
source("src/loading/sample_annotations.R")         # --sample_metadata, --sample_pattern
source("src/harmonization/panel_homogenization.R") # marker name recovery, aliases
source("src/harmonization/panel_comparison.R")     # panel comparison
source("src/visualization/umap_load_report.R")     # --report_json
source("src/visualization/umap_from_flowjo.R")     # umap_from_flowjo
source("src/visualization/umap_gate_layout.R")     # layout = "gates", filter, split_by
source("src/visualization/umap_style.R")           # look of the styled plots
source("src/visualization/umap_split_plots.R")     # combination, split, fold-change plots

umap_from_flowjo(
  wsp_path   = "data/Killing_AML.wsp",   # or a folder of .wsp files
  fcs_dir    = "data",
  population = "Lymphocytes/Single Cells/Viable cells",
  output_dir = "results/umap_viable",
  subpop_mode = "furthest",
  output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE)
)
```

This is the list `scripts/umap_from_flowjo.R` sources. `src/loading/flowjo_workspace.R`
holds the workspace discovery and FCS pairing (`find_workspace_files`,
`index_fcs_dir`, `resolve_workspace_specs`, …) that `umap_from_flowjo()` calls.
Without `umap_gate_layout.R`, `layout = "gates"`, `filter`, `split_by`,
`level_order` and `split_reference` stop with a message naming it; without
`umap_style.R` and `umap_split_plots.R` the combination, split and fold-change
plots are skipped with a note.

`umap_from_flowjo()` parameters mirror the CLI options (e.g. `markers`,
`subpop_select`, `annotate` and `split_by` take character vectors rather than
comma-separated strings). With `list_populations = TRUE` it prints the inspect
blocks and returns (invisibly) a list of `populations`, `samples` and
`markers`; otherwise it returns (invisibly) the `umap_coords`, the master
`umap_df`, `workspaces`, `dropped`, the resolved `output_dir` and `layout`, and
for a gates or split run also `cells` (one row per plotted cell),
`gate_layout` (calls, weights, decisions and counts; gates layout only) and
`plot_index`.

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

```bash
Rscript test_data/test_umap_gate_layout.R
```

covers the pure functions of `--layout gates` the same way: decisions, calls,
anchors, features, labels, sampling, filters, level order, the change table and
the palette, on hand-built gates that mirror the MTB tree.

```bash
Rscript test_data/test_run_record.R
```

covers the settings record: every option in the command, the command parsing
back to the same options (negative numbers, empty strings, spaces and quotes),
the input fingerprints, the outcome, and one record kept per run.

---

## Related

- [FlowJo → UMAP Dash app](../../app/README.md) — a GUI over this script
- [Gate-anchored UMAP layout](../methods/gate-anchored-layout.md): the method behind `--layout gates`
- [`src/loading/flowjo_workspace.R`](../subprojects/loading.md#flowjo_workspacer) — the workspace discovery and FCS pairing functions
- [Getting Started](../getting-started.md) — environment setup
- [Pipeline D — FlowJo Comparison](../pipelines/pipeline-D-flowjo-comparison.md) — inter-rater gating comparison from `.wsp` files
- [Pipeline E — Preprocessed UMAP](../pipelines/pipeline-E-preprocessed-umap.md) — UMAP/clustering from preprocessed FCS (no gating)
- [Standalone Scripts](../standalone-scripts.md) — using module functions directly
