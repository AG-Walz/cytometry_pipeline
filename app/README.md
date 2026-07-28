[< Back to Documentation Index](../docs/README.md)

# FlowJo → UMAP — Dash front-end

A GUI for [`scripts/umap_from_flowjo.R`](../scripts/umap_from_flowjo.R), and
nothing else. It builds the same command line you would type by hand, shows it
to you before running it, executes it, streams its output, and displays the
plots the script wrote. There is no separate analysis path: if the GUI produces
a figure, the CLI produces the identical figure from the command it printed.

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

## The four steps

The order is enforced, because it has to be: the gating tree and the marker
panel do not exist until the workspace has been read.

**1 — Data.** Pick the `.wsp` file and the FCS directory, either by typing a
path or through the Browse tree, then **Load workspace**. That runs
`--list_populations`, which builds the full `GatingSet`; on a large workspace it
takes a while, and the Log tab shows it happening.

**2 — Population & markers.** Click a node in the gating tree to set
`--population`. Below it, either pick a `--subpop_mode`, or turn on **explicit
gate set** to get checkboxes and tick the gates you want as the colour
partition (`--subpop_select`). The badges under the tree always show exactly
what will be passed. Markers come from the same load; leaving the field empty
means "every marker common to all samples", which is the script's own default.

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

The gallery reads the output directory rather than tracking what the script
said it wrote, so it shows the truth on disk. Files are grouped into:

* **Overview** — `by_file` and `by_subpopulation`
* **Markers — all samples** — one conjoint plot per marker
* **Markers — per sample** — one accordion per marker, whose images are only
  loaded once you open it. A full panel across several samples is hundreds of
  files, and Dash exposes no lazy-loading attribute on images, so the gallery
  has to control that itself.

Marker and sample dropdowns filter all three. Clicking an image opens it full
size. PNG is preferred for preview: an SVG of a large embedding carries one node
per cell and can be tens of megabytes.

---

## Layout

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

### The one dependency on the R side

Step 2 needs the marker panel, which only R can compute. `--list_populations`
therefore prints three parseable blocks — `Samples in loaded group (N):`,
`Population tree (N nodes):` and `Common markers (N):` — and
`services/rscript.py` reads them. Both sides carry a comment saying the format
is a shared contract. If you change those headers in
`src/visualization/umap_from_flowjo.R`, change the parser too.

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
