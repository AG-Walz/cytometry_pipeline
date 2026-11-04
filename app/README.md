[< Back to Documentation Index](../docs/README.md)

# FlowJo → UMAP — Dash front-end

A GUI for two of this repo's standalone scripts, on two tabs down the left:

| Tab | Script | What it does |
|---|---|---|
| **UMAP** | [`umap_from_flowjo.R`](../scripts/umap_from_flowjo.R) | FlowJo workspaces + FCS → one shared UMAP |
| **FCS names** | [`sync_fcs_names.R`](../scripts/sync_fcs_names.R) | make each FCS file's `$FIL` match its name on disk, and show what every file has stored inside it |

Both build the same command line you would type by hand, show it to you before
running it, execute it, and stream the output into the shared **Log** pane on
the right. There is no separate analysis path: if the GUI produces a figure or
a rewritten file, the CLI produces the identical one from the command it
printed.

---

## Running it

```bash
conda activate cytometry_env
python app/app.py
# → http://127.0.0.1:8050
```

Working from another machine? Do **not** pass `--host 0.0.0.0`; tunnel instead:

```bash
ssh -L 8050:localhost:8050 user@host
```

The app binds to `127.0.0.1` deliberately. Its file browser is rooted at `/`
with no confinement and no authentication, so anyone who can reach the port
gets a filesystem browser and a job runner.

| Flag | Default | Purpose |
|------|---------|---------|
| `--host` | `127.0.0.1` | Interface to bind |
| `--port` | `8050` | Port to bind |
| `--debug` | off | Dash hot-reload and error pages |

| Environment variable | Purpose |
|----------------------|---------|
| `UMAP_RSCRIPT` | How to invoke R. Defaults to `Rscript`. Set to `conda run -n cytometry_env Rscript` if Python and R are in different environments. |
| `UMAP_APP_CACHE` | Where the remembered form state is stored (default `~/.cache/umap_flowjo_app`). |
| `UMAP_APP_TMPDIR` | Temporary folder for runs (Annoy index, HDF5 cytoframes; gigabytes on a large cohort). Default: `<output_dir>/_run/tmp`, i.e. on the results volume rather than `/tmp`. |

---

## Serving it from PBIAppLauncher

`app/run.sh <port>` activates the environment and starts gunicorn on
`0.0.0.0:<port>`:

```json
{
  "name": "FlowJo UMAP",
  "start_cmd": "cd /mnt/volume/apps/cytometry_pipeline && git pull origin main && exec bash app/run.sh {port}",
  "hosts": ["highmem1", "highmem2"]
}
```

Two constraints that are not negotiable:

* **One gunicorn worker.** The active run lives in a module-level object, so a
  second worker would own a second, unrelated runner and the browser would poll
  whichever process its request reached. `run.sh` sets `--workers 1`;
  concurrency comes from threads.
* **The environment must be *activated*, not just referenced.** Calling
  `$CONDA_PREFIX/bin/gunicorn` directly skips the activation hook that
  `setup_env.sh` installs, and the `Rscript` subprocess then finds the system R
  (no CytoML) or loads ABI-incompatible packages from a global `~/R`.

**Serving it this way makes it network-reachable, which the app is not designed
for.** The file browser is rooted at `/` with no confinement and no
authentication, and anyone who can reach the port can start jobs. It is also
single-user: two people using one instance share one runner and one log. Only
put it behind the launcher on a trusted network, and preferably one at a time
per host.

---

## The FCS names tab

For the case where a whole run was acquired under one file name, you renamed
the files on disk, and FlowJo still shows a column of identical samples —
because it displays the `$FIL` keyword stored *inside* each file, not the file
name. Full background in
[docs/scripts/sync-fcs-names.md](../docs/scripts/sync-fcs-names.md).

The panel is a folder, a mode, and a button. The five modes are listed
read-only first, and that ordering is the safety design:

| Mode | Flags | Writes? |
|---|---|---|
| Preview — list what would change | — | no |
| Apply — write `$FIL` into the original files | `--apply --in_place` | **yes** |
| Undo — restore the names from the last run's ledger | `--revert` | **yes** |

Exactly one of them writes, it is the only red button, and Preview is what the
page loads on. The mode is deliberately **not** remembered between sessions — a
refresh must never bring the app back armed to write. Everything a write does
is recorded in `fcs_name_sync_ledger.csv` beside the data, which is what Undo
replays, and the name being replaced is kept in `ORIGINALFIL`.

Preview is the dry run rather than the CLI's `--probe`, because the dry run's
output is a strict superset: it prints the whole census *and* the per-file
plan. Three of the script's flags have a UI (`--fcs_dir`, the mode, `--also`);
the rest stay on the command line on purpose, and the command preview is
copy-pasteable. `--allow` especially: it downgrades a refusal that exists to
stop the tool writing to a file it does not fully understand, and that should
be a deliberate typed act, not a checkbox.

### Stored metadata

The TEXT segment of an FCS file is a delimiter-packed blob you cannot read
without a parser, so **Show metadata** decodes it for the whole folder at once:
one row per keyword, one column per file, with the keyword column frozen while
you scroll sideways through the cohort. Identity keywords come first; the
switch expands it to all of them.

A cell is highlighted when its value is shared by more than one file **and** the
keyword is a name — `$FIL`, `GUID`, `$SRC`, `TUBE NAME` and friends. A cohort
sharing `$TOT` or `$DATE` is ordinary and is not flagged, and neither is
`ORIGINALFIL`, which is *supposed* to be identical everywhere: it records the
one name the whole cohort used to have.

This runs R synchronously rather than through the shared runner, because it
only reads each file's header and TEXT segment (a few KB) and you are waiting
for it. It refuses while a run is in progress, since reading a file mid-rewrite
would give a torn view.

---

## The four steps

The order is enforced, because it has to be: the gating tree and the marker
panel do not exist until the workspaces have been read.

**1 — Data.** Pick the **folder** holding your workspaces and the **folder**
holding the FCS files, either by typing a path or through the Browse tree, then
**Load workspaces**. Both are searched recursively: the workspace field shows
how many `.wsp` it found before you commit to a load, and each workspace is
paired with its own FCS files by file name, so the folders need no particular
layout. A single `.wsp` file still works. Loading runs `--list_populations`,
which builds a full `GatingSet` per workspace; on large workspaces it takes a
while, and the Log tab shows it happening.

**Sample group** sits under **Advanced**, collapsed, because group 1 is "All
Samples" in almost every workspace and nobody should have to think about it. It
stays in step 1 rather than moving to Parameters for a reason: the group decides
which samples are loaded *at all*, so it has to be set before the load that
builds the gating tree. It accepts an index or a group name — a name is what
covers several workspaces that order their groups differently. A value other
than the default (or a blank one) is shown as a badge on the collapsed header,
so a setting that changes which samples get analysed can never hide behind a
closed panel.

**2 — Population & markers.** Click a node in the gating tree to set
`--population`. Below it, either pick a `--subpop_mode`, or turn on **explicit
gate set** to get checkboxes and tick the gates you want as the colour
partition (`--subpop_select`). The badges under the tree always show exactly
what will be passed. Markers come from the same load; leaving the field empty
means "every marker common to all samples", which is the script's own default.

With more than one workspace, the tree and the marker list are the
**intersection** across them — gates and markers every workspace has. That is
not a display choice: the R side refuses a `--population` it cannot resolve in
every workspace, so anything outside the intersection could not be run anyway.
The Log tab names the gates that were left out and which workspace lacks them.

**Sample metadata.** Further down step 2, upload a CSV, TSV or XLSX with one
row per FCS file (`--sample_metadata`), or type a regular expression with named
groups that reads the columns out of the sample names (`--sample_pattern`). The
file names usually go in a column called `filename`. R keys the table on the
first of `sample_id`, `sample`, `sampleID`, `sample_name`, `filename` and `file`
that it has, in that order of preference (so a `sample` column wins over
`filename`), otherwise on the first column. The app applies the same rule, and
the check under the upload names the column it keyed on and any FCS file
without a row. R itself resolves the table for the preview
(`--annotation_preview`), so the columns and levels the form offers are the
ones a run will see.

**Download metadata template** writes a CSV with one row per loaded sample,
keyed on its FCS file name as it is on disk. With a name pattern in the box,
its groups become pre-filled columns (R's `(?<name>...)` syntax is accepted).
A group named `sample_id`, `sample`, `sampleID` or `sample_name` is written as
`<name>_group`, since R would otherwise key the table on it, and the note under
the button says so. Fill in the rest in a spreadsheet and upload it. Generating
the names rather than typing them matters: a retyped `MTB016` for a file called
`MTB16_...` silently fails to match.

**Polyfunctionality maps** (on by default; off passes
`--no_polyfunctionality`) apply when the cells are labelled by gate
combination, i.e. with the `combinations` mode or the gate layout: every cell
coloured by how many gates below its lineage it is in, more gates larger and
on top, in the blues of the MTB cohort maps. See
[Polyfunctionality maps](../docs/scripts/umap-from-flowjo.md#polyfunctionality-maps).

With a metadata source, three more controls apply:

| Control | Flag | What it does |
|---|---|---|
| Split plots by | `--split_by` | per column: one plot per level on the same embedding, all levels side by side, and a map of each level against the reference that colours every combination by its fold change in frequency |
| Level order | `--level_order` | `column=a,b;column2=c`; the first level is the reference of the fold-change maps. Empty means natural order (`W5` before `W12`). The hint below the box lists the detected levels in that order, ready to copy and rearrange |
| Keep only samples where | `--filter` | one column, its values, and include or exclude, giving `col=a,b` or `col!=a,b`; applied before any cell is read |

Without a source these are disabled and none of them is passed. A split or
filter column does not have to be among **Annotations to plot**; R reads it
whatever `--annotate` selects, and adds a split column to the annotated ones. The run is refused when a split or filter names
a column the metadata lacks, a filter value is not one of its levels, or the
level order is malformed or names a column that is not split.

Every setting explains itself on hover: point at the small ⓘ next to its label
(or tab to it) for what it does. The form itself only shows the settings.

**3 — Parameters.** The layout comes first:

* **Marker values** *(default)*, `--layout markers`: the embedding of the
  marker values themselves. **Fit the embedding on** chooses what UMAP is fitted
  on before every cell is projected onto it and plotted: **Random cells**
  (**Random cells to fit on**, sent as `--umap_subsample`; default 300000, so
  a larger run is fitted on that many and the rest projected; 0 = every cell) or
  **SuperCellCyto supercells** (`--fit_subset supercell`, with **Cells per
  supercell** as `--supercell_gamma`, **Supercell neighbours** as
  `--supercell_k_knn` and **Max supercells** as `--supercell_max`, which
  raises gamma on a large cohort so the fit stays at most that size; needs
  the R package SuperCellCyto, see `install_extra_packages.R`). **CPU cores**
  does not change the coordinates: it threads the plots and the neighbour
  search, and the optimiser stays serial. Step 2's `combinations` colouring labels the
  cells by gate combination as the gate layout does, without the gates
  entering the embedding.
* **Gate-separated islands**, `--layout gates`: every
  combination of gate calls under the population is forced onto its own
  island, so the map shows the gating rather than how similar the cells are. Values
  are re-centred within each sample on that sample's own cells of the same
  combination, rare combinations are kept rather than downsampled away, and
  the embedding is fitted on a set balanced across combinations. Its settings
  sit in **Gate layout settings**
  (`--gate_gap`, `--gate_tier_ratio` as **Lineage gap ratio**,
  `--gate_softness`, `--supervise`, `--cells_per_combination`,
  `--min_combination_cells`, and the fit set size `--umap_subsample`).
  Normalisation and downsampling do not apply and are disabled.
  `--viability_pattern`, the regular expression for marker names that are
  never gate decisions (viability dyes by default; empty disables it), is
  CLI-only: add it to the printed command by hand.

**Per-marker normalisation** (`--normalize_data`) changes what the embedding
is computed from; **Marker colour scale** (`--normalize_umap`) only changes
what the marker plots are coloured by. Its default, **Same as the embedding**,
shows the values the map was built from: z-scores for a z-scored run, FlowJo
values otherwise. Choose FlowJo values, z-score or min-max to override it; the
points do not move.

The form starts on the marker layout, as R does, and `--layout` is always
written out, and so is `--min_dist`, whose R default
depends on the layout (0.1 for markers, 0.25 for gates; switching the layout
swaps the field between the two unless you typed your own value). Everything
else left at its R default is omitted from the command line, so the printed
command stays readable, and a layout's own settings are only passed with that
layout. Point size and alpha left empty ("auto") are omitted too: R then sizes
the points for the number of cells drawn (about 1.1 for 10k cells, 0.3 for
1M, 0.2 for 3M; per-sample plots for their own count, styled plots per
combination). A value typed in applies to every plot.

**4 — Run.** Choose a parent directory and a folder name (pre-filled with the
population and a timestamp, so consecutive runs never overwrite each other).
To continue an interrupted or changed run, choose its folder instead: the
step lists the steps saved there, and with **Reuse saved steps** on (the
default) the run continues from the first step whose settings changed; off
passes `--no_resume`.
The exact command appears above the Run button. Output formats are exactly what
you ticked in step 3 — nothing is added behind your back, which does mean a
PDF-only run has no inline preview.

While a run is active the Run button is disabled and **Cancel** kills the
process group. Only one thing runs at a time, inspect included.

---

## Results

The gallery reads the **Results folder** at the top of the tab rather than
tracking what the script said it wrote, so it shows the truth on disk. A
finished run fills the folder in; typing or browsing to an earlier run shows
that run instead, and the value is remembered with the rest of the form, so a
reload or a server restart does not lose it. It reads the run's subfolders
(`overview/`, `markers/`, `splits/<column>/`, `changes/<column>/`, `data/`,
...; see the output layout in `docs/scripts/umap-from-flowjo.md`) as well as
the flat folders of older runs, and skips `checkpoints/`. Files are grouped
into:

* **Overview** — `by_combination` first (gate layout, or any run split by
  metadata), then the polyfunctionality map (gate combinations), then `by_file`,
  `by_workspace` (only when several workspaces were loaded) and
  `by_subpopulation`
* **Split by metadata**: one accordion per split column, holding the panel with all
  levels, then each level, then the fold-change maps against the reference
  level, then the polyfunctionality map of each level, in the order R drew them
* **Markers — all samples** — one conjoint plot per marker
* **Markers — per sample** — one accordion per marker, whose images are only
  loaded once you open it. A full panel across several samples is hundreds of
  files, and Dash exposes no lazy-loading attribute on images, so the gallery
  has to control that itself.

Marker and sample dropdowns filter the marker sections; a sample filter
also hides the annotation and split sections, which are not per sample. The
run's tables (`cells.csv.gz`, `combination_frequencies_<column>.csv`,
`separation_metrics.csv`, `run_info.json` and the rest) are listed under
**Data files** for download. Clicking an image opens it full size. PNG is preferred for preview and is the default output format, in the form
and in the CLI: an SVG of a large embedding carries one
node per cell and can be tens of megabytes, so an SVG over 8 MB is offered
for download instead of being embedded.

---

## Layout

The tab bar wraps the **left column only**, so the Log and Results panes stay
visible for both tools. Tabs around the whole page would hide the log behind
the UMAP tab — and the log is where the FCS tool's entire output goes. Its id is
`tool-tabs`, not `tabs`: `tabs` belongs to the right column and is both a State
and an Output of the poll callback.


```
app/
  app.py                    entrypoint and CLI flags
  config.py                 defaults mirroring the R CLI, paths, R invocation
  layout.py                 the page
  components/filebrowser.py the reusable dmc.Tree path picker
  callbacks/
    common.py               form collection, validation, parameter mapping
    steps.py                navigation, gating tree, metadata, command preview
    execution.py            start/cancel/poll the subprocess
    results.py              the gallery
    fcsnames.py             the FCS names tab
  components/health.py      the Data health panel
  services/
    fs_tree.py              lazy filesystem listing for the browser
    gating.py               population paths -> tree data
    rscript.py              argv construction, load report and annotation preview reading
    runner.py               single-slot subprocess runner with a live log
    detached.py             runs that outlive the app (log and exit code in <out>/_run)
    report_cache.py         remembered workspace loads
    results_scan.py         output directory -> gallery structure
    serving.py              Flask route for plot files outside assets/
    session_state.py        remembered form state
  assets/dmcfuncs.js        custom dmc.Tree node renderer (file/folder icons)
```

### Why the FCS names tab needed almost no new code

`RUNNER` takes a free-form `kind`, and the poll callback in
`callbacks/execution.py` branches on `kind == "inspect"` / `kind == "run"` and
falls through for anything else — while the tuple it always returns already
drives the status badge, the elapsed clock, the log and the Cancel button. So
`kind="fcs-names"` gets all of that for **zero** lines in `execution.py`. The
folder pickers are `browser_field`s, whose callbacks are pattern-matching and
already cover them. `RUNNER` being single-slot is what stops a sync and a UMAP
run overlapping. `callbacks/fcsnames.py` is the whole tab.

The one structured thing the tab needs — the keyword table — comes back from R
as a CSV rather than being parsed out of stdout, because R owns the only
correct FCS parser in this repo and a second one in Python would drift from it.

### What the app reads from the R side

Nothing is parsed out of R's stdout any more. Everything structured comes back
as a file:

| File | Written by | Read by |
|---|---|---|
| load report (`--report_json`, schema 1.x) | `write_load_report()` in `src/visualization/umap_load_report.R` | `rscript.read_report` / `inspect_from_report`; a different major schema is refused rather than half-read |
| annotation preview (`--annotation_preview`) | `preview_annotations()`, same file | `rscript.describe_annotations`: columns, levels per column, match counts |
| `plot_index.json` in the output folder | the run | `results_scan.read_plot_index`: kind, column, level and reference of each split plot |

The plot index exists because split file names cannot be read back reliably:
`<column>_<level>` is ambiguous when a column name contains `_`, and names are
cut at 100 characters. Without an index (an older run) the gallery falls back
to the name tokens `_by_combination`, `_split_`, `_diff_`.

Two things are still read from the log: `::progress::` lines (from
`emit_progress()` in `src/R_shared/io_utils.R`), which drive the progress bar
and are hidden from the log, and, for detached runs started before exit codes
were recorded, the closing `Done. Plots written to:` line.

The progress bar reads only the LAST run in `run.log` (each run into the same
output folder adds to that log, after a `── started` line). The time left
comes from a reference run: the phase times of the last complete run with the
same settings, kept in `~/.cache/umap_flowjo_app/phase_history.json`, or else
the last complete run in the same log. For settings that have no complete run
yet, the time is predicted from the size of the job, which R states in a
`::runinfo::` line before the embedding (cells, fitted cells, neighbours,
threads; see `app/services/run_time_model.py`), and the bar says "roughly".
Each host corrects that prediction with its own complete runs. `Resuming: '<step>'` lines tell it which phases a run reuses from
`checkpoints/`. See `app/services/progress_estimate.py`.

---

## Notes

- **Not multi-user.** The active run lives in a module-level object. Two
  browsers pointed at one instance share one runner and one log.
- **R buffers stdout when it is a pipe**, which would make the log arrive only
  at the end. The app prefixes the command with `stdbuf -oL -eL` where that
  exists; without it (macOS without coreutils) the log updates in chunks.
- **Marker and sample names in the gallery are the sanitised ones**, since they
  are recovered from filenames that went through `create_safe_filename()`.
- **The form is remembered** between launches in `~/.cache/umap_flowjo_app/`.
  **Reset form** clears it. A form saved before the gate layout existed opens
  on that layout with its `min_dist` of 0.25; the other remembered values are
  kept.

---

## Related

- [UMAP from a FlowJo Workspace Population](../docs/scripts/umap-from-flowjo.md) — the script itself
- [Getting Started](../docs/getting-started.md) — environment setup
