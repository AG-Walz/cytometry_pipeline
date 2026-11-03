"""Shared plumbing for the callbacks: form collection, validation, params.

The form spans four stepper panes and two dozen components, so the mapping from
component ids to CLI parameter names lives here once and every callback that
needs the whole form reuses it.
"""

from __future__ import annotations

import csv
import io
import re
from pathlib import Path

from dash import Input, State

from app import config
from app.services import fs_tree, gating

#: (parameter name, component id, component prop), in a fixed order that
#: ``as_dict`` relies on.
FORM_FIELDS: list[tuple[str, object, str]] = [
    ("wsp", {"comp": "fb-input", "bid": "wsp"}, "value"),
    ("fcs_dir", {"comp": "fb-input", "bid": "fcs"}, "value"),
    ("output_parent", {"comp": "fb-input", "bid": "outparent"}, "value"),
    ("output_name", "in-outname", "value"),
    ("resume", "sw-resume", "checked"),
    ("group", "in-group", "value"),
    ("allow_missing_samples", "sw-allow-missing", "checked"),
    ("match_mode", "in-match-mode", "value"),
    ("marker_aliases", "in-marker-aliases", "value"),
    ("recover_markers", "sw-recover-markers", "checked"),
    ("sample_metadata", {"comp": "fb-input", "bid": "meta"}, "value"),
    ("sample_pattern", "in-pattern", "value"),
    ("annotate", "in-annotate", "value"),
    ("split_by", "in-split-by", "value"),
    ("level_order", "in-level-order", "value"),
    ("filter_column", "in-filter-col", "value"),
    ("filter_values", "in-filter-values", "value"),
    ("filter_mode", "in-filter-mode", "value"),
    ("markers", "in-markers", "value"),
    ("explicit_subpop", "sw-explicit", "checked"),
    ("subpop_mode", "subpop-mode", "value"),
    ("layout", "in-layout", "value"),
    ("gate_gap", "in-gate-gap", "value"),
    ("gate_tier_ratio", "in-gate-tier-ratio", "value"),
    ("gate_softness", "in-gate-softness", "value"),
    ("supervise", "in-supervise", "value"),
    ("cells_per_combination", "in-cells-per-comb", "value"),
    ("min_combination_cells", "in-min-comb-cells", "value"),
    ("umap_subsample", "in-fit-cells", "value"),
    ("fit_subset", "in-fit-subset", "value"),
    ("random_fit_cells", "in-random-fit-cells", "value"),
    ("supercell_gamma", "in-sc-gamma", "value"),
    ("supercell_k_knn", "in-sc-knn", "value"),
    ("supercell_max", "in-sc-max", "value"),
    ("normalize", "in-normalize", "value"),
    ("normalize_umap", "in-colour-scale", "value"),
    ("clip_quantile", "in-clip", "value"),
    ("gate_separation", "in-gate-separation", "value"),
    ("downsample", "in-downsample", "value"),
    ("n_cores", "in-cores", "value"),
    ("n_neighbors", "in-neighbors", "value"),
    ("min_dist", "in-mindist", "value"),
    ("metric", "in-metric", "value"),
    ("seed", "in-seed", "value"),
    ("point_size", "in-pointsize", "value"),
    ("point_alpha", "in-pointalpha", "value"),
    ("plot_width", "in-width", "value"),
    ("plot_height", "in-height", "value"),
    ("output_formats", "in-formats", "value"),
    ("results_dir", {"comp": "fb-input", "bid": "results"}, "value"),
    ("population_map", "store-popmap", "data"),
    ("excluded_workspaces", "store-wsp-exclude", "data"),
    ("selected_gates", "gate-tree", "selected"),
    ("checked_gates", "gate-tree", "checked"),
    # What R's --annotation_preview said about the metadata: its columns and
    # their levels. Read by validation only; never persisted (not in DEFAULTS)
    # and never passed to R (not in run_params).
    ("annotation_preview", "store-annot-preview", "data"),
]

AUTO_NAME_RE = re.compile(r"^umap_.*_\d{8}_\d{4}$")

#: Columns R tries, in this order, as the sample identifier of a metadata
#: table before falling back to the first column. Mirrors
#: metadata_id_column() in src/loading/sample_annotations.R.
METADATA_ID_COLUMNS = ("sample_id", "sample", "sampleID", "sample_name",
                       "filename", "file")

#: The level R's annotation code gives a sample it could not annotate.
UNASSIGNED = "Unassigned"

#: Gate-layout settings, the form key and how it reads in a message.
GATE_PARAMS = (
    ("gate_gap", "Gap at each gate cut"),
    ("gate_tier_ratio", "Lineage gap ratio"),
    ("gate_softness", "Softness"),
    ("supervise", "Supervision weight"),
    ("cells_per_combination", "Cells per combination per sample"),
    ("min_combination_cells", "Minimum cells per combination"),
    ("umap_subsample", "Fit cells"),
)
#: Their component ids, in the same order.
GATE_COMPONENTS = [cid for key, _ in GATE_PARAMS
                   for name, cid, _ in FORM_FIELDS if name == key]

#: The markers layout's fit-set settings, the form key and how it reads.
FIT_PARAMS = (
    ("fit_subset", "Fit on"),
    ("random_fit_cells", "Random cells to fit on"),
    ("supercell_gamma", "Cells per supercell"),
    ("supercell_k_knn", "Supercell neighbours"),
    ("supercell_max", "Max supercells"),
)
#: Their component ids, in the same order.
FIT_COMPONENTS = [cid for key, _ in FIT_PARAMS
                  for name, cid, _ in FORM_FIELDS if name == key]


def metadata_id_column(columns: list[str]) -> str:
    """The column R keys a metadata table on, or "" for a table with none."""
    for cand in METADATA_ID_COLUMNS:
        if cand in columns:
            return cand
    return columns[0] if columns else ""


def has_metadata_source(form: dict) -> bool:
    return bool((form.get("sample_metadata") or "").strip()
                or (form.get("sample_pattern") or "").strip())


def filter_spec(form: dict) -> str:
    """The ``--filter`` value for the one clause the form offers, or "".

    Empty until both a column and at least one value are chosen: a clause
    with no values would keep (or drop) nothing, which is never what a
    half-filled form means.
    """
    column = (form.get("filter_column") or "").strip()
    values = [v for v in (form.get("filter_values") or []) if str(v).strip()]
    if not column or not values:
        return ""
    op = "!=" if form.get("filter_mode") == "exclude" else "="
    return f"{column}{op}{','.join(values)}"


def metadata_template(inspect: dict | None, pattern: str = "") -> tuple[str, str]:
    """``(csv_text, note)`` for "Download metadata template".

    One row per loaded sample, keyed on its FCS file name: the GatingSet name
    minus the "_<n>" a cytoset appends after ".fcs". Generated from the load
    rather than typed, because a retyped name ("MTB016" for a file called
    "MTB16") is exactly the row that silently fails to match.

    With a name pattern, its named groups become pre-filled columns, applied
    to the display name as R's --sample_pattern is. R writes groups as
    (?<name>...), Python as (?P<name>...); lookbehinds (?<= and (?<! are left
    alone. A group named like an ID column R prefers to "filename" (sample,
    sample_id, ...) gets a "_group" suffix, or R would key the table on it.
    ``note`` is "" unless the pattern could not be used or a group was renamed.
    """
    report = (inspect or {}).get("report") or {}
    samples = report.get("samples") or [
        {"display": s} for s in (inspect or {}).get("samples") or []]

    rx, groups, headers, note = None, [], [], ""
    pattern = (pattern or "").strip()
    if pattern:
        try:
            rx = re.compile(re.sub(r"\(\?<(?![=!])", "(?P<", pattern))
            groups = [g for g in sorted(rx.groupindex, key=rx.groupindex.get)
                      if g != "filename"]
        except re.error as exc:
            rx, note = None, (f"The name pattern could not be read ({exc}), so "
                              "the template holds the file names only.")
    outranks = METADATA_ID_COLUMNS[:METADATA_ID_COLUMNS.index("filename")]
    renamed = []
    for g in groups:
        header = g
        if g in outranks:
            while header in outranks or header in groups or header in headers:
                header += "_group"
            renamed.append(f"{g} to {header}")
        headers.append(header)
    if renamed:
        note = ("Renamed " + ", ".join(renamed) + " in the template: R keys a "
                "metadata table on a column of that name before filename, so "
                "no row would match.")

    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(["filename", *headers])
    seen: set[str] = set()
    for s in samples:
        raw = str(s.get("raw") or s.get("display") or "").strip()
        name = re.sub(r"(\.fcs)_\d+$", r"\1", raw, flags=re.IGNORECASE)
        if not name or name in seen:
            continue
        seen.add(name)
        match = rx.search(str(s.get("display") or raw)) if rx else None
        writer.writerow([name, *[(match.group(g) or "") if match else ""
                                 for g in groups]])
    return out.getvalue(), note


def parse_level_order(text: str) -> tuple[dict[str, list[str]], list[str]]:
    """``"col=a,b;col2=c"`` -> ({col: [a, b], col2: [c]}, problems).

    The same split rules as parse_level_order() in umap_gate_layout.R: ";"
    between columns, the first "=" between column and levels, "," between
    levels, whitespace trimmed.
    """
    order: dict[str, list[str]] = {}
    problems: list[str] = []
    for part in [p.strip() for p in (text or "").split(";") if p.strip()]:
        column, eq, levels = part.partition("=")
        if not eq or not column.strip():
            problems.append(f"Level order entry '{part}' is not of the form "
                            "column=a,b,c (step 2).")
            continue
        order[column.strip()] = [v.strip() for v in levels.split(",")
                                 if v.strip()]
    return order, problems


def workspace_problem(wsp: str) -> str:
    """Why ``wsp`` is not a usable --wsp value, or "" when it is.

    The field takes a folder (searched recursively for .wsp, which is the whole
    point of the multi-workspace mode) but a single .wsp file still works, so
    both are accepted and the folder case is checked for actually containing
    something.
    """
    wsp = (wsp or "").strip()
    if not wsp:
        return "Select the folder holding your .wsp workspace(s)."
    path = Path(wsp)
    if path.is_dir():
        found, _ = fs_tree.find_files(wsp, config.WSP_EXTENSIONS)
        if not found:
            return f"No .wsp workspace anywhere under: {wsp}"
        return ""
    if path.is_file():
        if path.suffix.lower() not in config.WSP_EXTENSIONS:
            return f"Not a .wsp workspace: {wsp}"
        return ""
    return f"Workspace path does not exist: {wsp}"


def form_states() -> list[State]:
    return [State(cid, prop) for _, cid, prop in FORM_FIELDS]


def form_inputs() -> list[Input]:
    return [Input(cid, prop) for _, cid, prop in FORM_FIELDS]


def as_dict(values: tuple | list) -> dict:
    """Zip raw callback values back into a named form dict."""
    return {name: values[i] for i, (name, _, _) in enumerate(FORM_FIELDS)}


def population_of(form: dict) -> str:
    selected = form.get("selected_gates") or []
    return selected[0] if selected else ""


def subpop_selection(form: dict, populations: list[str]) -> list[str]:
    """The ``--subpop_select`` list: checked gates, relative to the population.

    Only strict descendants of the chosen population qualify — the R side
    matches entries against the paths under ``--population``, so a checked node
    elsewhere in the tree would simply fail to resolve.
    """
    if not form.get("explicit_subpop"):
        return []
    population = population_of(form)
    if not population:
        return []
    allowed = set(gating.descendants(populations, population))
    checked = [c for c in (form.get("checked_gates") or []) if c in allowed]
    return [gating.relative(c, population) for c in sorted(checked)]


def output_dir_of(form: dict) -> str:
    parent = (form.get("output_parent") or "").strip()
    name = (form.get("output_name") or "").strip()
    if not parent or not name:
        return ""
    return str(Path(parent) / name)


def workspace_selection(form: dict, inspect: dict | None) -> list[str]:
    """The ``--workspaces`` INCLUDE list, or [] meaning "all of them".

    The UI stores what is switched OFF; the CLI takes what to switch ON. The
    conversion needs to know what loaded, which only the report does — so an
    empty exclusion set returns [] and the flag is omitted entirely rather than
    spelling out all thirty-three names for a run that excludes nothing.
    """
    excluded = set(form.get("excluded_workspaces") or [])
    if not excluded:
        return []
    loaded = list((inspect or {}).get("populations_by_workspace") or {})
    return [w for w in loaded if w not in excluded]


def population_map_of(form: dict) -> dict:
    """The ``--population_map``, minus workspaces this run does not load.

    Two filters, both of which R would otherwise turn into a hard stop:

    An entry with no gate is dropped because "workspace X, no gate" and "no
    entry for workspace X" mean the same thing — X contributes no cells — and
    only the second is a thing the flag can express.

    An entry naming an EXCLUDED workspace is dropped because ``--workspaces``
    is applied at discovery, before any workspace is read, so a workspace
    switched off in the health panel is never loaded at all — and the map
    naming one is exactly the "population_map names workspace(s) that were not
    loaded" error, raised at the end of a load that took minutes. Dropping it
    costs nothing: excluded and unmapped already mean the same thing.
    """
    excluded = set(form.get("excluded_workspaces") or [])
    return {w: p for w, p in (form.get("population_map") or {}).items()
            if p and w not in excluded}


def run_params(form: dict, populations: list[str],
               inspect: dict | None = None) -> dict:
    """Everything ``rscript.run_argv`` needs, derived from the form.

    An explicit whitelist, not a copy of ``form``: gate selections and other
    UI-only state must not reach the command line. The cost is that a NEW FORM
    FIELD IS INVISIBLE TO THE RUN UNTIL IT IS LISTED HERE — it reaches the form
    dict, is silently dropped on the way to the command, and the only symptom is
    a flag that does nothing. Add the key here whenever you add one to
    ``FORM_FIELDS``.
    """
    return {
        "wsp": (form.get("wsp") or "").strip(),
        "fcs_dir": (form.get("fcs_dir") or "").strip(),
        "population": population_of(form),
        "output_dir": output_dir_of(form),
        "group": form.get("group"),
        "population_map": population_map_of(form),
        "allow_missing_samples": bool(form.get("allow_missing_samples")),
        "match_mode": form.get("match_mode") or config.DEFAULTS["match_mode"],
        "workspaces": workspace_selection(form, inspect),
        "marker_aliases": form.get("marker_aliases") or "",
        "recover_markers": bool(form.get("recover_markers", True)),
        "sample_metadata": (form.get("sample_metadata") or "").strip(),
        "sample_pattern": (form.get("sample_pattern") or "").strip(),
        "annotate": form.get("annotate") or [],
        "split_by": form.get("split_by") or [],
        "level_order": (form.get("level_order") or "").strip(),
        "filter": filter_spec(form),
        "markers": form.get("markers") or [],
        "subpop_mode": form.get("subpop_mode"),
        "subpop_select": subpop_selection(form, populations),
        "layout": form.get("layout") or config.FORM_START["layout"],
        **{key: form.get(key) for key, _ in GATE_PARAMS},
        **{key: form.get(key) for key, _ in FIT_PARAMS},
        "normalize": form.get("normalize"),
        "normalize_umap": form.get("normalize_umap"),
        "clip_quantile": form.get("clip_quantile"),
        "gate_separation": form.get("gate_separation"),
        "downsample": form.get("downsample"),
        "n_cores": form.get("n_cores"),
        "n_neighbors": form.get("n_neighbors"),
        "min_dist": form.get("min_dist"),
        "metric": form.get("metric"),
        "seed": form.get("seed"),
        "point_size": form.get("point_size"),
        "point_alpha": form.get("point_alpha"),
        "plot_width": form.get("plot_width"),
        "plot_height": form.get("plot_height"),
        "output_formats": form.get("output_formats") or [],
        "resume": bool(form.get("resume", True)),
    }


def validate_run(form: dict, inspect: dict) -> list[str]:
    """Problems that would make a run fail or be meaningless, in report order."""
    problems: list[str] = []
    fcs_dir = (form.get("fcs_dir") or "").strip()

    wsp_problem = workspace_problem(form.get("wsp"))
    if wsp_problem:
        problems.append(wsp_problem)

    if not fcs_dir:
        problems.append("No FCS directory selected (step 1).")
    elif not Path(fcs_dir).is_dir():
        problems.append(f"FCS directory does not exist: {fcs_dir}")
    # Same rule as validate_inspect. The field's placeholder reads "1", so a
    # cleared field looks set, and without this the run went out with
    # --group '' and R rejected it in every workspace.
    if not str(form.get("group") or "").strip():
        problems.append("No sample group set (step 1, Advanced): enter an "
                        "index such as 1, or a name.")

    picked = {w: p for w, p in (form.get("population_map") or {}).items() if p}
    pop_map = population_map_of(form)
    loaded = inspect.get("populations_by_workspace") or {}
    if not inspect.get("populations") and not loaded:
        problems.append("Workspaces not loaded yet — run 'Load workspaces' in step 1.")
    elif picked and not pop_map:
        # Every workspace with a gate is switched off, so the run has no cells
        # at all. R says so too, but only after a load that takes minutes.
        problems.append(
            "Every workspace with a per-workspace population is switched off "
            "in step 1, so nothing would go into the UMAP. Switch one back on, "
            "or pick a gate for a workspace that is on (step 2)."
        )
    elif pop_map:
        # A gate chosen for a workspace that is no longer loaded is a hard
        # error on the R side. The mapping is pruned on every load, so this
        # should be unreachable from the UI — it survives as a backstop for a
        # hand-edited session file, and for a load that reported no workspaces.
        stale = [w for w in pop_map if loaded and w not in loaded]
        if stale:
            problems.append(
                "The per-workspace populations name workspace(s) that are not "
                "loaded: " + ", ".join(sorted(stale)) +
                ". Load the workspaces again in step 1 — that drops them."
            )
    elif not population_of(form):
        problems.append("No population selected (step 2).")

    if form.get("explicit_subpop") and not subpop_selection(
        form, inspect.get("populations") or []
    ):
        problems.append(
            "'Explicit gate set' is on but no gate under the selected "
            "population is ticked (step 2)."
        )

    problems += _metadata_problems(form)
    problems += _gate_layout_problems(form)

    parent = (form.get("output_parent") or "").strip()
    if not parent:
        problems.append("No output parent directory selected (step 4).")
    elif not Path(parent).is_dir():
        problems.append(f"Output parent directory does not exist: {parent}")
    if not (form.get("output_name") or "").strip():
        problems.append("No output folder name (step 4).")

    if not (form.get("output_formats") or []):
        problems.append("No output format ticked (step 3) — the run would write nothing.")

    return problems


def _metadata_problems(form: dict) -> list[str]:
    """Split, level order and filter settings R would reject.

    Checked against the columns and levels R's own annotation preview found,
    and only when it found some: while the preview is empty or failed, step 2
    already says why, and repeating it here as "unknown column" would point at
    the wrong cause. Without a metadata source none of these is passed at all.
    """
    if not has_metadata_source(form):
        return []
    preview = form.get("annotation_preview") or {}
    columns = list(preview.get("columns") or [])
    levels = preview.get("level_values") or {}
    problems: list[str] = []

    split_by = list(form.get("split_by") or [])
    if columns:
        unknown = [c for c in split_by if c not in columns]
        if unknown:
            problems.append("'Split plots by' names column(s) the metadata "
                            "does not have: " + ", ".join(unknown) + " (step 2).")
    if split_by:
        order, syntax = parse_level_order(form.get("level_order") or "")
        problems += syntax
        stray = [c for c in order if c not in split_by]
        if stray:
            problems.append("Level order names column(s) the plots are not "
                            "split by: " + ", ".join(stray) + " (step 2).")

    column = (form.get("filter_column") or "").strip()
    values = [str(v) for v in (form.get("filter_values") or [])]
    if column:
        if columns and column not in columns:
            problems.append(f"'Keep only samples where' uses column '{column}', "
                            "which the metadata does not have (step 2).")
        if not values:
            problems.append(f"'Keep only samples where': pick at least one "
                            f"value of '{column}', or clear the column (step 2).")
        # --filter separates clauses with ";" and values with ",", and has no
        # way to quote either.
        awkward = [v for v in values if "," in v or ";" in v]
        if awkward:
            problems.append("A filter value containing ',' or ';' cannot be "
                            "passed to --filter: " + ", ".join(awkward)
                            + ". Rename the level in the metadata (step 2).")
        known = levels.get(column) or []
        missing = [v for v in values if known and v not in known]
        if missing:
            problems.append(f"Filter value(s) not found in '{column}': "
                            + ", ".join(missing) + " (step 2).")
    return problems


def _number(value) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _gate_layout_problems(form: dict) -> list[str]:
    """Settings of the gate layout R would reject, or that make no sense."""
    layout = form.get("layout")
    if layout not in {value for value, _ in config.LAYOUTS}:
        return ["Choose a layout (step 3)."]
    if layout != "gates":
        return _fit_problems(form)
    problems: list[str] = []
    for key, label in GATE_PARAMS:
        value = _number(form.get(key))
        if value is None:
            problems.append(f"{label} is empty (step 3).")
        elif key == "gate_tier_ratio" and value < 1:
            problems.append(f"{label} must be 1 or more (step 3).")
        elif key == "gate_softness" and value <= 0:
            problems.append(f"{label} must be above 0 (step 3).")
        elif key == "supervise" and not 0 <= value <= 1:
            problems.append(f"{label} must be between 0 and 1 (step 3).")
        elif key == "cells_per_combination" and value < 1:
            problems.append(f"{label} must be at least 1 (step 3).")
        elif value < 0:
            problems.append(f"{label} must be 0 or more (step 3).")
        elif key in ("cells_per_combination", "min_combination_cells",
                     "umap_subsample") and value != int(value):
            problems.append(f"{label} must be a whole number (step 3).")
    return problems


def _fit_problems(form: dict) -> list[str]:
    """Fit-set settings of the markers layout R would reject."""
    sep = _number(form.get("gate_separation"))
    if form.get("gate_separation") not in (None, "") and (sep is None or sep < 0):
        return ["Gate separation must be 0 or more (step 3)."]
    fit_subset = form.get("fit_subset") or config.DEFAULTS["fit_subset"]
    if fit_subset not in ("random", "supercell"):
        return ["Choose what UMAP is fitted on (step 3)."]
    keys = (("supercell_gamma", "supercell_k_knn", "supercell_max")
            if fit_subset == "supercell" else ("random_fit_cells",))
    labels = dict(FIT_PARAMS)
    problems: list[str] = []
    for key in keys:
        label = labels[key]
        value = _number(form.get(key))
        if value is None:
            problems.append(f"{label} is empty (step 3).")
        elif key == "supercell_gamma" and value <= 1:
            problems.append(f"{label} must be above 1 (step 3).")
        elif key == "supercell_k_knn" and value < 1:
            problems.append(f"{label} must be at least 1 (step 3).")
        elif key == "supercell_max" and value < 1000:
            problems.append(f"{label} must be at least 1,000 (step 3).")
        elif value < 0:
            problems.append(f"{label} must be 0 or more (step 3).")
        elif key in ("random_fit_cells", "supercell_k_knn", "supercell_max") \
                and value != int(value):
            problems.append(f"{label} must be a whole number (step 3).")
    return problems


#: R's checkpoint stages (umap_from_flowjo.R, ckpt_save), in run order.
CHECKPOINT_STAGES = (("extract", "extraction"), ("supercells", "supercells"),
                     ("model", "UMAP model"), ("embedding", "coordinates"))


def saved_steps(output_dir: str) -> list[str]:
    """The steps an earlier run saved in ``output_dir/checkpoints``, in order."""
    if not output_dir:
        return []
    ckpt = Path(output_dir) / "checkpoints"
    return [label for stage, label in CHECKPOINT_STAGES
            if (ckpt / f"{stage}.rds").is_file()]


def validate_inspect(form: dict) -> list[str]:
    """Problems that would stop the workspaces from loading."""
    problems: list[str] = []
    fcs_dir = (form.get("fcs_dir") or "").strip()
    wsp_problem = workspace_problem(form.get("wsp"))
    if wsp_problem:
        problems.append(wsp_problem)
    if not fcs_dir:
        problems.append("Select the FCS directory.")
    elif not Path(fcs_dir).is_dir():
        problems.append(f"FCS directory does not exist: {fcs_dir}")
    elif not fs_tree.find_files(fcs_dir, (".fcs",))[0]:
        problems.append(f"No .fcs file anywhere under: {fcs_dir}")
    if not str(form.get("group") or "").strip():
        problems.append("Enter a sample group (an index such as 1, or a name).")
    return problems


def session_snapshot(form: dict, populations: list[str]) -> dict:
    """The subset of the form worth remembering between launches.

    The population and its gate set are omitted while nothing is selected, so
    a page load — where they are necessarily empty until the workspace has been
    re-read — does not erase what was remembered. ``session_state.save``
    merges, so an omitted key keeps its stored value.
    """
    state = {k: form.get(k) for k in config.DEFAULTS if k in form}
    population = population_of(form)
    if population:
        state["population"] = population
        state["subpop_select"] = subpop_selection(form, populations)
    return state
