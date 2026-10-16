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

**3 — Parameters.** UMAP and plot options. Anything left at its default is
omitted from the command line, so the printed command stays readable.

**4 — Run.** Choose a parent directory and a folder name (pre-filled with the
population and a timestamp, so consecutive runs never overwrite each other).
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
reload or a server restart does not lose it. Files are grouped into:

* **Overview** — `by_file`, `by_workspace` (only when several workspaces were
  loaded) and `by_subpopulation`
* **Markers — all samples** — one conjoint plot per marker
* **Markers — per sample** — one accordion per marker, whose images are only
  loaded once you open it. A full panel across several samples is hundreds of
  files, and Dash exposes no lazy-loading attribute on images, so the gallery
  has to control that itself.

Marker and sample dropdowns filter all three. Clicking an image opens it full
size. PNG is preferred for preview and is the form's default output format
(the CLI's own default stays `svg`): an SVG of a large embedding carries one
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
    steps.py                navigation, gating tree, command preview
    execution.py            start/cancel/poll the subprocess
    results.py              the gallery
  services/
    fs_tree.py              lazy filesystem listing for the browser
    gating.py               population paths -> tree data
    rscript.py              argv construction and inspect-output parsing
    runner.py               single-slot subprocess runner with a live log
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

### The one dependency on the R side

Step 2 needs the marker panel, which only R can compute. `--list_populations`
therefore prints four parseable blocks — `Samples in loaded group (N):`,
`Population tree (N nodes):`, `Common markers (N):` and `Workspaces resolved
(N):` — and `services/rscript.py` reads them. Both sides carry a comment saying
the format is a shared contract. If you change those headers in
`src/visualization/umap_from_flowjo.R`, change the parser too.

The last block, and the `Populations not shared by all workspaces` block after
it, are indented **four** spaces rather than two. That is load-bearing: the
first three blocks are collected by two-space item patterns, and a two-space
indent here would let the population parser swallow the divergence list.

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
  **Reset form** clears it.

---

## Related

- [UMAP from a FlowJo Workspace Population](../docs/scripts/umap-from-flowjo.md) — the script itself
- [Getting Started](../docs/getting-started.md) — environment setup
