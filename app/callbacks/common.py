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
    ("sample_metadata", {"comp": "fb-input", "bid": "meta"}, "value"),
    ("sample_pattern", "in-pattern", "value"),
    ("annotate", "in-annotate", "value"),
    ("markers", "in-markers", "value"),
    ("explicit_subpop", "sw-explicit", "checked"),
    ("subpop_mode", "subpop-mode", "value"),
    ("normalize", "in-normalize", "value"),
    ("downsample", "in-downsample", "value"),
    ("n_neighbors", "in-neighbors", "value"),
    ("min_dist", "in-mindist", "value"),
    ("metric", "in-metric", "value"),
    ("seed", "in-seed", "value"),
    ("point_size", "in-pointsize", "value"),
    ("point_alpha", "in-pointalpha", "value"),
    ("plot_width", "in-width", "value"),
    ("plot_height", "in-height", "value"),
    ("output_formats", "in-formats", "value"),
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


def run_params(form: dict, populations: list[str]) -> dict:
    """Everything ``rscript.run_argv`` needs, derived from the form."""
    return {
        "wsp": (form.get("wsp") or "").strip(),
        "fcs_dir": (form.get("fcs_dir") or "").strip(),
        "population": population_of(form),
        "output_dir": output_dir_of(form),
        "group": form.get("group"),
        "sample_metadata": (form.get("sample_metadata") or "").strip(),
        "sample_pattern": (form.get("sample_pattern") or "").strip(),
        "annotate": form.get("annotate") or [],
        "markers": form.get("markers") or [],
        "subpop_mode": form.get("subpop_mode"),
        "subpop_select": subpop_selection(form, populations),
        "normalize": form.get("normalize"),
        "downsample": form.get("downsample"),
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

    if not inspect.get("populations"):
        problems.append("Workspaces not loaded yet — run 'Load workspaces' in step 1.")
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
