"""Shared plumbing for the callbacks: form collection, validation, params.

The form spans four stepper panes and two dozen components, so the mapping from
component ids to CLI parameter names lives here once and every callback that
needs the whole form reuses it.
"""

from __future__ import annotations

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
    ("group", "in-group", "value"),
    ("allow_missing_samples", "sw-allow-missing", "checked"),
    ("match_mode", "in-match-mode", "value"),
    ("marker_aliases", "in-marker-aliases", "value"),
    ("recover_markers", "sw-recover-markers", "checked"),
    ("sample_metadata", {"comp": "fb-input", "bid": "meta"}, "value"),
    ("sample_pattern", "in-pattern", "value"),
    ("annotate", "in-annotate", "value"),
    ("markers", "in-markers", "value"),
    ("explicit_subpop", "sw-explicit", "checked"),
    ("subpop_mode", "subpop-mode", "value"),
    ("normalize", "in-normalize", "value"),
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
]

AUTO_NAME_RE = re.compile(r"^umap_.*_\d{8}_\d{4}$")


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
        "markers": form.get("markers") or [],
        "subpop_mode": form.get("subpop_mode"),
        "subpop_select": subpop_selection(form, populations),
        "normalize": form.get("normalize"),
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
