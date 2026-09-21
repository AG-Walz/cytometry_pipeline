"""Callbacks for the form itself: navigation, gating tree, and the command."""

from __future__ import annotations

import datetime as _dt
from pathlib import Path

import dash_mantine_components as dmc
from dash import ALL, Input, Output, State, callback, ctx, html, no_update

from app import config
from app.callbacks import common
from app.components import health
from app.services import fs_tree, gating, rscript, session_state

#: Last state written to the session file, so it is only rewritten on a change.
_last_saved: dict | None = None

SUBPOP_HELP = {
    "furthest": (
        "Each cell takes the most specific gate it reaches — a unique label per "
        "cell, no '(multiple)'."
    ),
    "cd4cd8": "Only the CD4+/CD8+ lineage split. Fewest colours, most readable.",
    "direct": (
        "Every immediate child gate. These usually overlap, so most cells "
        "collapse to '(multiple)'."
    ),
    "deepest": (
        "Deepest gate by path depth; siblings at equal depth collapse to "
        "'(multiple)'."
    ),
    "none": "Skip the subpopulation plot entirely.",
}


def _alert(children, color="red", title=None):
    return dmc.Alert(children, color=color, title=title, variant="light")


def _problem_list(problems: list[str]):
    return html.Ul(
        [html.Li(p) for p in problems],
        style={"margin": "4px 0 0 0", "paddingLeft": "18px"},
    )


def register_callbacks() -> None:
    # -- step 1: what is in the workspace folder ---------------------------

    @callback(
        Output("wsp-scan-preview", "children"),
        Input({"comp": "fb-input", "bid": "wsp"}, "value"),
    )
    def preview_workspaces(wsp):
        """Show what a --wsp folder actually contains, before loading it.

        A pure Python glob, so it is instant: it costs nothing to be told "4
        workspaces found" here rather than after a multi-minute GatingSet build.
        """
        wsp = (wsp or "").strip()
        if not wsp:
            return None
        problem = common.workspace_problem(wsp)
        if problem:
            return _alert(problem, color="yellow")

        found, truncated = fs_tree.find_files(wsp, config.WSP_EXTENSIONS)
        root = Path(wsp)
        base = root if root.is_dir() else root.parent
        names = [str(Path(f).relative_to(base)) for f in found]
        if len(names) == 1:
            return dmc.Text(f"1 workspace: {names[0]}", size="xs", c="dimmed")
        return dmc.Alert(
            [
                dmc.Text(
                    f"{len(names)} workspaces found"
                    + (" (list truncated)" if truncated else "")
                    + " — all of them go into one shared UMAP.",
                    size="sm",
                    fw=600,
                ),
                dmc.Space(h=4),
                dmc.Group(
                    [dmc.Badge(n, variant="light", size="sm", tt="none")
                     for n in names],
                    gap=4,
                ),
            ],
            color="blue",
            variant="light",
        )

    @callback(
        Output("advanced-badge", "children"),
        Input("in-group", "value"),
        Input("sw-allow-missing", "checked"),
    )
    def advanced_badge(group, allow_missing):
        """Surface non-default sample loading on the collapsed Advanced header.

        A setting that changes which samples are loaded must not be able to
        hide behind a closed panel — that is the one way this section could
        cost someone a run. "Missing FCS allowed" earns a badge for the
        opposite reason to the others: it cannot fail a run, it can quietly
        shrink one, and a UMAP that silently lost a subject looks exactly like
        a UMAP that did not.
        """
        badges = []
        group = str(group or "").strip()
        if not group:
            badges.append(dmc.Badge("no group set", color="red",
                                    variant="light", size="sm", tt="none"))
        elif group != str(config.DEFAULTS["group"]):
            badges.append(dmc.Badge(f"group: {group}", color="blue",
                                    variant="light", size="sm", tt="none"))
        if allow_missing:
            badges.append(dmc.Badge("missing FCS allowed", color="orange",
                                    variant="light", size="sm", tt="none"))
        if not badges:
            return None
        return dmc.Group(badges, gap="xs")

    @callback(
        Output("params-badge", "children"),
        Input("in-metric", "value"),
        Input("in-seed", "value"),
        Input("in-pointsize", "value"),
        Input("in-pointalpha", "value"),
        Input("in-width", "value"),
        Input("in-height", "value"),
    )
    def params_badge(*values):
        """Name anything non-default that is hidden behind the closed panel.

        Same rule as advanced_badge: collapsing a control is only acceptable if
        a changed one still announces itself. Without this, a run seeded
        differently — or plotted at a size that makes the legend unreadable —
        would look identical to a default run from the outside.
        """
        keys = ("metric", "seed", "point_size", "point_alpha",
                "plot_width", "plot_height")
        changed = [
            f"{k}: {v}" for k, v in zip(keys, values)
            if v is not None and v != config.DEFAULTS[k]
        ]
        if not changed:
            return None
        return dmc.Group(
            [dmc.Badge(c, color="blue", variant="light", size="sm", tt="none")
             for c in changed],
            gap="xs",
        )

    # -- stepper navigation ------------------------------------------------

    @callback(
        Output("stepper", "active"),
        Input("btn-back", "n_clicks"),
        Input("btn-next", "n_clicks"),
        State("stepper", "active"),
        prevent_initial_call=True,
    )
    def navigate(_back, _next, active):
        active = active or 0
        if ctx.triggered_id == "btn-back":
            return max(0, active - 1)
        return min(3, active + 1)

    @callback(
        Output("btn-back", "disabled"),
        Output("btn-next", "disabled"),
        Input("stepper", "active"),
        Input("store-inspect", "data"),
        Input("gate-tree", "selected"),
    )
    def nav_availability(active, inspect, selected):
        active = active or 0
        # "Loaded" is "at least one workspace came in", NOT "there is a shared
        # gate". Those differ exactly when the workspaces were gated
        # independently — which is when the per-workspace population picker in
        # step 2 is the thing you need, so blocking the way to it was backwards.
        inspect = inspect or {}
        loaded = bool(inspect.get("populations")
                      or inspect.get("populations_by_workspace"))
        has_population = bool(selected)
        if active == 0:
            forward_blocked = not loaded
        elif active == 1:
            forward_blocked = not has_population
        elif active == 3:
            forward_blocked = True
        else:
            forward_blocked = False
        return active == 0, forward_blocked

    # -- data health -------------------------------------------------------
    #
    # Pure rendering: the panel shows what R reported and works nothing out.
    # See app/components/health.py.

    @callback(
        Output("health-panel", "children"),
        Output("health-badge", "children"),
        Input("store-inspect", "data"),
        Input("store-wsp-exclude", "data"),
    )
    def fill_health(inspect, excluded):
        return (health.panel(inspect, excluded),
                health.headline(inspect, excluded))

    @callback(
        Output("store-wsp-exclude", "data"),
        Input({"comp": "wsp-toggle", "bid": ALL}, "checked"),
        Input({"comp": "wsp-all"}, "n_clicks"),
        State({"comp": "wsp-toggle", "bid": ALL}, "id"),
        State("store-wsp-exclude", "data"),
        State("store-inspect", "data"),
        prevent_initial_call=True,
    )
    def collect_workspace_selection(checked, _all_clicks, ids, current,
                                    inspect):
        """Which workspaces are switched off.

        Stored as the EXCLUDED set: the default — use everything that loaded —
        is then the empty list, and a workspace that appears in the folder later
        is included rather than silently dropped by an inclusion list written
        before it existed.
        """
        # ctx.triggered_id is a dict for a pattern id; compare on the part
        # that identifies the control rather than the whole dict.
        if isinstance(ctx.triggered_id, dict) and \
                ctx.triggered_id.get("comp") == "wsp-all":
            return []
        # Only workspaces that LOADED. One that failed draws its switch off and
        # disabled (health._workspace_row) because it cannot be included — so
        # reading that switch as "the user switched this off" records a choice
        # nobody made, and puts a workspace in the excluded set that the user
        # cannot see a reason for and cannot take back out.
        loaded = (inspect or {}).get("populations_by_workspace") or {}
        excluded = [i["bid"] for i, on in zip(ids or [], checked or [])
                    if not on and (not loaded or i["bid"] in loaded)]
        # Re-rendering the panel re-fires this with the values it just drew, so
        # returning an equal-but-new list would loop between the two callbacks.
        if sorted(excluded) == sorted(current or []):
            return no_update
        return excluded

    @callback(
        Output("acc-health", "value"),
        Input("store-inspect", "data"),
        prevent_initial_call=True,
    )
    def open_health_when_it_matters(inspect):
        """Open the panel on its own only when something needs attention.

        A load where every workspace came in and every panel agrees has nothing
        in here worth reading, and popping it open each time would train people
        to close it without looking — which is precisely when they would miss
        the load that lost eleven workspaces.
        """
        report = (inspect or {}).get("report")
        if not report:
            return no_update
        totals = report.get("totals") or {}
        markers = report.get("markers") or {}
        trouble = (
            totals.get("workspaces_loaded") != totals.get("workspaces_found")
            or not markers.get("common")
            or bool(markers.get("not_common"))
        )
        return "health" if trouble else None

    @callback(
        Output("in-marker-aliases", "value"),
        Input({"comp": "alias-suggestions"}, "value"),
        State("in-marker-aliases", "value"),
        prevent_initial_call=True,
    )
    def apply_alias_suggestions(picked, current):
        """Fold the ticked suggestions into the alias field.

        Merged rather than replaced, and de-duplicated on the left-hand side, so
        an alias typed by hand is not thrown away by ticking a box — and so
        ticking two suggestions that name the same spelling cannot produce a map
        R will reject as a duplicate.
        """
        chosen: list[str] = []
        for group in picked or []:
            for value in group or []:
                chosen.extend(v for v in value.split(";") if v)

        merged: dict[str, str] = {}
        for pair in [p for p in (current or "").split(";") if p] + chosen:
            written, _, used = pair.partition("=")
            if written.strip() and used.strip():
                merged[written.strip()] = used.strip()
        return rscript.alias_string(merged)

    # -- gating tree -------------------------------------------------------

    @callback(
        Output("gate-tree", "data"),
        Output("gate-tree", "expanded"),
        Output("gate-tree", "selected"),
        Output("gate-tree", "checked"),
        Output("gate-tree", "checkboxes"),
        Output("gate-tree", "selectOnClick"),
        Output("step2-summary", "children"),
        Output("subpop-mode-wrap", "style"),
        Output("gate-tree-mode", "children"),
        Input("store-inspect", "data"),
        Input("sw-explicit", "checked"),
        Input("sw-per-wsp", "checked"),
        Input("sel-wsp", "value"),
        State("gate-tree", "selected"),
        State("store-popmap", "data"),
    )
    def fill_gate_tree(inspect, explicit, per_wsp, wsp, current_selected,
                       pop_map):
        """Sole owner of the gating tree's contents and mode.

        The tree does two jobs — pick the population, or tick a gate set — and
        both the workspace load and the explicit switch can change which. They
        are handled together so neither can overwrite the other's tree.
        """
        inspect = inspect or {}
        by_wsp = inspect.get("populations_by_workspace") or {}
        explicit = bool(explicit)
        per_wsp = bool(per_wsp) and bool(by_wsp)
        pop_map = pop_map or {}

        # In per-workspace mode the tree shows ONE workspace's own gates —
        # otherwise the intersection, which is what a single --population must
        # live within.
        if per_wsp:
            if wsp not in by_wsp:
                wsp = next(iter(by_wsp))
            populations = by_wsp.get(wsp) or []
            gate_labels = (inspect.get("gate_labels_by_workspace")
                           or {}).get(wsp) or {}
        else:
            populations = inspect.get("populations") or []
            gate_labels = inspect.get("population_labels") or {}

        if not populations:
            return (
                [], [], [], [], False, True,
                _alert("Load a workspace in step 1 to see its gating tree.",
                       color="gray"),
                {"display": "none"} if explicit else {},
                None,
            )

        # Keep the current population if it still exists (e.g. after reloading
        # the same workspace), else fall back to the remembered one.
        keep = [
            p
            for p in ([pop_map.get(wsp) or ""] if per_wsp else
                      list(current_selected or [])
                      + [inspect.get("remembered_population") or ""])
            if p in set(populations)
        ]
        selected = keep[:1]
        population = selected[0] if selected else ""

        if explicit:
            # A flat list of the gates under the population: Mantine reports
            # only leaf nodes in `checked`, so nesting would make an
            # intermediate gate impossible to tick. selectOnClick is off so a
            # click cannot silently change the population underneath.
            data = gating.flat_subgate_tree(populations, population)
            expanded: list[str] = []
            mode_note = dmc.Text(
                "Gates under the selected population.", size="xs", c="dimmed",
                mb=6,
            )
        else:
            data = gating.build_tree(populations)
            expanded = gating.default_expanded(populations, depth=2)
            mode_note = None

        samples = inspect.get("samples") or []
        markers = inspect.get("markers") or []
        workspaces = inspect.get("workspaces") or []
        multi = len(workspaces) > 1
        body = [
            dmc.Text(
                (f"{len(workspaces)} workspace(s), " if multi else "")
                + f"{len(samples)} sample(s), {len(populations)} gate(s), "
                f"{len(markers)} common marker(s).",
                size="sm",
                fw=600,
            ),
        ]
        if multi:
            # With several workspaces the gate list is their INTERSECTION, and
            # saying so saves you hunting for a gate only one of them has. In
            # per-workspace mode that is no longer true.
            body += [
                dmc.Space(h=4),
                dmc.Text(
                    (f"Tree: {wsp}'s own gates." if per_wsp else
                     "Tree: gates shared by every workspace."),
                    size="xs", c="dimmed",
                ),
            ]
        summary = dmc.Alert(body, color="green", variant="light")
        return (
            data,
            expanded,
            selected,
            [],
            explicit,
            not explicit,
            summary,
            {"display": "none"} if explicit else {},
            mode_note,
        )

    # -- per-workspace population picker -----------------------------------

    @callback(
        Output("sel-wsp", "data"),
        Output("sel-wsp", "value"),
        Output("per-wsp-select-wrap", "style"),
        Input("sw-per-wsp", "checked"),
        Input("store-inspect", "data"),
        Input("store-popmap", "data"),
        State("sel-wsp", "value"),
    )
    def fill_workspace_select(per_wsp, inspect, pop_map, current):
        by_wsp = (inspect or {}).get("populations_by_workspace") or {}
        if not per_wsp or not by_wsp:
            return [], "", {"display": "none"}   # "" not None — see layout.py
        pop_map = pop_map or {}
        labels = list(by_wsp)
        if current not in by_wsp:
            current = labels[0]
        data = [
            {"value": w,
             "label": f"{w}  ({'set' if pop_map.get(w) else 'not set'})"}
            for w in labels
        ]
        return data, current, {}

    @callback(
        Output("per-wsp-picker", "children"),
        Input("sw-per-wsp", "checked"),
        Input("store-inspect", "data"),
        Input("store-popmap", "data"),
    )
    def per_workspace_picker(per_wsp, inspect, pop_map):
        """A running list of what each workspace contributes.

        The safeguard, not decoration: this mode pools cells from gates with
        different names, so which gate every workspace contributes has to be
        readable at a glance BEFORE the run, not inferred afterwards from the
        plots.
        """
        by_wsp = (inspect or {}).get("populations_by_workspace") or {}
        if not per_wsp:
            return None
        if not by_wsp:
            return _alert("Load the workspaces in step 1 first.", color="gray")

        pop_map = pop_map or {}
        labels = list(by_wsp)
        chosen = sum(1 for w in labels if pop_map.get(w))
        rows = [
            dmc.Group(
                [
                    dmc.Badge("set" if pop_map.get(w) else "not set",
                              color="green" if pop_map.get(w) else "gray",
                              variant="light", size="sm", tt="none",
                              style={"minWidth": "72px"}),
                    dmc.Text(w, size="xs", fw=600,
                             style={"minWidth": "190px"}),
                    # The prop is OMITTED rather than set to None for the
                    # default colour: None crosses to the browser as JSON null,
                    # `typeof null` is "object", and Mantine rejects it with
                    # "Failed to parse color ... instead got object" — which
                    # takes down the whole step, not just this line.
                    dmc.Text(pop_map.get(w) or "not in the UMAP", size="xs",
                             style={"flex": 1, "wordBreak": "break-all"},
                             **({} if pop_map.get(w) else {"c": "dimmed"})),
                    dmc.Button("clear", id={"comp": "popmap-clear", "wsp": w},
                               size="compact-xs", variant="subtle",
                               disabled=not pop_map.get(w)),
                ],
                gap="xs", wrap="nowrap",
            )
            for w in labels
        ]
        return [
            dmc.Alert(
                [
                    dmc.Text(f"{chosen} of {len(labels)} workspace(s) will "
                             f"contribute cells.", size="sm", fw=600),
                    dmc.Space(h=6),
                    dmc.Stack(rows, gap=4),
                ],
                color="blue" if chosen else "yellow",
                variant="light",
                title="What each workspace contributes",
            ),
            dmc.Space(h=10),
        ]

    @callback(
        Output("store-popmap", "data"),
        Input("gate-tree", "selected"),
        Input({"comp": "popmap-clear", "wsp": ALL}, "n_clicks"),
        Input("sw-per-wsp", "checked"),
        Input("store-inspect", "data"),
        State("sel-wsp", "value"),
        State("store-popmap", "data"),
        prevent_initial_call=True,
    )
    def record_population_choice(selected, _clears, per_wsp, inspect, wsp,
                                 pop_map):
        """Sole writer of the per-workspace mapping."""
        pop_map = dict(pop_map or {})
        by_wsp = (inspect or {}).get("populations_by_workspace") or {}
        trigger = ctx.triggered_id

        # A LOAD is the only moment the real workspace names are known, so it is
        # where entries naming a workspace that is no longer loaded are dropped.
        # Without this the mapping could only ever GROW — the mapping persists
        # with the rest of the form, so pointing the app at a different folder
        # carried the previous study's names into the new one forever.
        #
        # It was unfixable from the UI, which is what made it worse than
        # untidy: the picker draws one row per LOADED workspace, so a stale
        # entry had no row and no "clear" button, and the run step's advice to
        # re-pick in step 2 could not be acted on — picking added a correct
        # entry BESIDE the stale one rather than replacing it.
        #
        # Only ever against a load that produced something. An empty inspect
        # means "nothing loaded yet" — the app has just started and the mapping
        # was restored from the session file — and pruning against nothing
        # would throw away exactly the choices that were remembered.
        if trigger == "store-inspect":
            if not by_wsp:
                return no_update
            kept = {w: g for w, g in pop_map.items() if w in by_wsp}
            return kept if len(kept) != len(pop_map) else no_update

        # Turning the mode OFF clears the mapping rather than leaving it to
        # apply invisibly: with it set, the R side ignores --population
        # entirely, so a stale mapping would silently decide the whole run.
        if trigger == "sw-per-wsp" and not per_wsp:
            return {}
        if not per_wsp:
            return no_update

        if isinstance(trigger, dict) and trigger.get("comp") == "popmap-clear":
            pop_map.pop(trigger.get("wsp"), None)
            return pop_map

        if trigger == "gate-tree":
            if not wsp or wsp not in by_wsp:
                return no_update
            gate = (selected or [None])[0]
            if gate:
                pop_map[wsp] = gate
            else:
                pop_map.pop(wsp, None)
            return pop_map
        return no_update

    @callback(
        Output("gate-tree", "expanded", allow_duplicate=True),
        Input("btn-expand-all", "n_clicks"),
        Input("btn-collapse", "n_clicks"),
        State("store-inspect", "data"),
        State("sw-per-wsp", "checked"),
        State("sel-wsp", "value"),
        prevent_initial_call=True,
    )
    def expand_collapse(_expand, _collapse, inspect, per_wsp, wsp):
        if ctx.triggered_id == "btn-collapse":
            return []
        # Must expand the tree that is actually SHOWN. Expanding the shared
        # list while the tree shows one workspace's own gates unfolds nothing
        # below the gates they happen to have in common — which is exactly the
        # part you switched modes to get at.
        inspect = inspect or {}
        by_wsp = inspect.get("populations_by_workspace") or {}
        if per_wsp and by_wsp:
            paths = by_wsp.get(wsp) or by_wsp[next(iter(by_wsp))]
        else:
            paths = inspect.get("populations") or []
        return gating.all_values(paths)

    @callback(
        Output("population-echo", "children"),
        Input("gate-tree", "selected"),
        Input("sw-per-wsp", "checked"),
        Input("sel-wsp", "value"),
    )
    def echo_population(selected, per_wsp, wsp):
        if per_wsp and wsp:
            if not selected:
                return dmc.Text(
                    f"{wsp} has no gate yet — it will be left out of the UMAP.",
                    size="sm", c="dimmed")
            return dmc.Text(
                [dmc.Text(f"{wsp} \u2192 ", span=True, c="dimmed", size="sm"),
                 dmc.Text(selected[0], span=True, size="sm", fw=600)],
                size="sm")
        return _echo_single_population(selected)

    def _echo_single_population(selected):
        if not selected:
            return _alert("Click a gate above to choose the population.",
                          color="yellow")
        return dmc.Text(
            [dmc.Text("--population ", span=True, c="dimmed", size="sm"),
             dmc.Text(selected[0], span=True, size="sm", fw=600)],
            style={"wordBreak": "break-all"},
        )

    @callback(
        Output("subpop-mode-help", "children"),
        Input("subpop-mode", "value"),
    )
    def mode_help(mode):
        return dmc.Text(SUBPOP_HELP.get(mode, ""), size="xs", c="dimmed")

    @callback(
        Output("subpop-chips", "children"),
        Input("gate-tree", "checked"),
        Input("gate-tree", "selected"),
        Input("sw-explicit", "checked"),
        State("store-inspect", "data"),
    )
    def show_subpop_selection(checked, selected, explicit, inspect):
        if not explicit:
            return None
        populations = (inspect or {}).get("populations") or []
        form = {
            "explicit_subpop": True,
            "selected_gates": selected,
            "checked_gates": checked,
        }
        chosen = common.subpop_selection(form, populations)
        if not chosen:
            return _alert("Tick gates under the selected population.",
                          color="yellow")
        return dmc.Box(
            [
                dmc.Text(f"--subpop_select ({len(chosen)} gate(s))",
                         size="xs", c="dimmed"),
                dmc.Space(h=4),
                dmc.Group(
                    # tt="none": Badge uppercases by default, which would
                    # misrepresent case-sensitive gate names like CD107a.
                    [dmc.Badge(c, variant="light", size="sm", tt="none")
                     for c in chosen],
                    gap=4,
                ),
            ],
            mt=8,
        )

    # -- sample metadata upload --------------------------------------------

    @callback(
        Output({"comp": "fb-input", "bid": "meta"}, "value"),
        Input("up-meta", "contents"),
        State("up-meta", "filename"),
        prevent_initial_call=True,
    )
    def store_metadata(contents, filename):
        """Save the upload as a comma CSV the R side can read, and hand back its path.

        XLSX is converted here because R only reads CSV; every column is kept
        as text so a timepoint like "01" stays "01".
        """
        import base64
        import io
        import pandas as pd

        raw = base64.b64decode(contents.split(",", 1)[1])
        if filename.lower().endswith(".xlsx"):
            df = pd.read_excel(io.BytesIO(raw), dtype=str)
        else:
            # Same rule as the R reader: ";" if the header has one, else ",".
            first = raw.split(b"\n", 1)[0]
            sep = ";" if b";" in first else "\t" if b"\t" in first else ","
            df = pd.read_csv(io.BytesIO(raw), dtype=str, sep=sep)
        dest = config.CACHE_DIR / "metadata" / (Path(filename).stem + ".csv")
        dest.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(dest, index=False)
        return str(dest)

    @callback(
        Output("meta-check", "children"),
        Output("meta-title", "children"),
        Output("meta-table", "columns"),
        Output("meta-table", "data"),
        Input({"comp": "fb-input", "bid": "meta"}, "value"),
        Input({"comp": "fb-input", "bid": "fcs"}, "value"),
    )
    def check_metadata(metadata, fcs_dir):
        """Validate the table and show it: first column = FCS file name, and
        every FCS file in the folder must have exactly one row."""
        import pandas as pd

        path = (metadata or "").strip()
        if not path:
            return None, "No metadata uploaded.", [], []
        try:
            df = pd.read_csv(path, dtype=str, keep_default_na=False)
        except Exception as exc:  # noqa: BLE001 - shown to the user as is
            return _alert(f"Could not read the metadata: {exc}"), "No metadata uploaded.", [], []

        key = lambda n: Path(str(n)).name.lower().removesuffix(".fcs")  # noqa: E731
        problems = []
        if df.shape[1] < 2:
            problems.append("At least two columns are needed: the file name plus one annotation.")
        if any(c.startswith("Unnamed:") or not c.strip() for c in df.columns):
            problems.append("Every column needs a header.")
        keys = [key(v) for v in df.iloc[:, 0]]
        if "" in keys:
            problems.append("The first column has empty file names.")
        dups = sorted({k for k in keys if keys.count(k) > 1} - {""})
        if dups:
            problems.append(f"Duplicate file names: {', '.join(dups[:5])}.")
        fcs, _ = fs_tree.find_files((fcs_dir or "").strip(), (".fcs",))
        missing = sorted(Path(f).name for f in fcs if key(f) not in set(keys))
        if missing:
            problems.append(f"{len(missing)} of {len(fcs)} FCS file(s) have no metadata row: "
                            + ", ".join(missing[:5]) + ("…" if len(missing) > 5 else ""))
        title = f"{Path(path).name}: {len(df)} row(s), {df.shape[1] - 1} annotation column(s)."
        note = (_alert(_problem_list(problems), color="red", title="Metadata problems")
                if problems else
                _alert(title + (f" All {len(fcs)} FCS file(s) covered." if fcs else
                                " Pick an FCS folder to check coverage."), color="green"))
        return note, title, [{"name": c, "id": c} for c in df.columns], df.to_dict("records")

    # -- sample annotations ------------------------------------------------

    @callback(
        Output("in-annotate", "data"),
        Output("in-annotate", "value"),
        Output("annotate-preview", "children"),
        Input({"comp": "fb-input", "bid": "meta"}, "value"),
        Input("in-pattern", "value"),
        Input("store-inspect", "data"),
        State("in-annotate", "value"),
    )
    def fill_annotations(metadata, pattern, inspect, current):
        """Offer the annotation columns, and show what they will resolve to.

        Done in Python rather than by asking R: a CSV header is a single line
        and a regex is a regex, so the picker can be populated and the match
        rate shown instantly instead of after a workspace load.
        """
        samples = (inspect or {}).get("samples") or []
        columns, preview = rscript.describe_annotations(
            (metadata or "").strip(), (pattern or "").strip(), samples)

        kept = [c for c in (current or []) if c in set(columns)]
        note = None
        if preview.get("error"):
            note = _alert(preview["error"], color="red")
        elif not columns:
            note = None
        elif preview.get("unmatched"):
            note = _alert(
                [
                    dmc.Text(
                        f"{preview['matched']} of {preview['total']} loaded "
                        f"sample(s) resolve; the rest become 'Unassigned'.",
                        size="sm", fw=600),
                    dmc.Space(h=4),
                    dmc.Text("Not resolved: "
                             + ", ".join(preview["unmatched"][:5])
                             + ("…" if len(preview["unmatched"]) > 5 else ""),
                             size="xs"),
                ],
                color="yellow",
            )
        else:
            detail = ", ".join(
                f"{c} ({preview['levels'].get(c, 0)} level(s))"
                for c in columns) if preview.get("levels") else ", ".join(columns)
            note = _alert(
                f"{preview['source']}: {detail}."
                + (f" All {preview['total']} loaded sample(s) resolve."
                   if preview.get("total") else ""),
                color="green",
            )
        return [{"value": c, "label": c} for c in columns], kept, note

    # -- markers -----------------------------------------------------------

    @callback(
        Output("in-markers", "data"),
        Output("in-markers", "value"),
        Input("store-inspect", "data"),
        State("in-markers", "value"),
        prevent_initial_call=True,
    )
    def fill_markers(inspect, current):
        markers = (inspect or {}).get("markers") or []
        # Drop anything that is not in this workspace, so a remembered
        # selection from another panel cannot silently fail the run.
        kept = [m for m in (current or []) if m in set(markers)]
        return [{"value": m, "label": m} for m in markers], kept

    # -- output directory --------------------------------------------------

    @callback(
        Output("in-outname", "value"),
        Input("gate-tree", "selected"),
        State("in-outname", "value"),
        prevent_initial_call=True,
    )
    def suggest_output_name(selected, current):
        if not selected:
            return no_update
        # Only overwrite an empty field or a previously auto-generated name, so
        # a name the user typed themselves is never clobbered.
        if current and not common.AUTO_NAME_RE.match(current):
            return no_update
        leaf = selected[0].rstrip("/").split("/")[-1]
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in leaf).strip("_")
        stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M")
        return f"umap_{safe.lower() or 'population'}_{stamp}"

    @callback(
        Output("output-echo", "children"),
        Input({"comp": "fb-input", "bid": "outparent"}, "value"),
        Input("in-outname", "value"),
    )
    def echo_output(parent, name):
        form = {"output_parent": parent, "output_name": name}
        path = common.output_dir_of(form)
        if not path:
            return _alert("Choose a parent directory and a folder name.",
                          color="yellow")
        exists = Path(path).exists()
        return dmc.Box(
            [
                dmc.Text(path, size="sm", fw=600,
                         style={"wordBreak": "break-all"}),
                dmc.Text(
                    "Directory already exists — files with matching names will "
                    "be overwritten."
                    if exists
                    else "Will be created by the run.",
                    size="xs",
                    c="orange" if exists else "dimmed",
                ),
            ]
        )

    @callback(
        Output("format-warning", "children"),
        Input("in-formats", "value"),
    )
    def format_warning(formats):
        formats = formats or []
        if not formats:
            return _alert("Nothing ticked — the run would produce no files.")
        if not set(formats) & {"png", "svg"}:
            return _alert(
                "PDF only: the results gallery cannot embed PDFs, so it will "
                "show download links instead.",
                color="blue",
            )
        if "svg" in formats and "png" not in formats:
            return _alert(
                "Large SVGs are not previewed in the gallery; tick PNG too.",
                color="blue",
            )
        return None

    # -- command preview and run validation --------------------------------

    @callback(
        Output("cmd-preview", "children"),
        Output("run-validation", "children"),
        Output("btn-run", "disabled"),
        *common.form_inputs(),
        Input("store-inspect", "data"),
        Input("poll", "n_intervals"),
    )
    def refresh_command(*values):
        form = common.as_dict(values)
        inspect = values[len(common.FORM_FIELDS)] or {}
        populations = inspect.get("populations") or []

        params = common.run_params(form, populations, inspect)
        problems = common.validate_run(form, inspect)

        # The command is shown even when incomplete, so it is obvious which
        # piece is still missing.
        try:
            command = rscript.command_string(rscript.run_argv(params))
        except (KeyError, TypeError) as exc:
            command = f"(cannot build command yet: {exc})"

        if problems:
            feedback = _alert(
                [dmc.Text("Not ready to run:", size="sm", fw=600),
                 _problem_list(problems)],
                color="yellow",
            )
        else:
            feedback = _alert("Ready to run.", color="green")

        # This callback already sees every form field on every change, so it is
        # also where the remembered state is kept up to date. The poll ticks it
        # too, hence the comparison — otherwise this would rewrite the file
        # more than once a second for no reason.
        global _last_saved
        snapshot = common.session_snapshot(form, populations)
        if snapshot != _last_saved:
            _last_saved = snapshot
            session_state.save(snapshot)

        from app.services.runner import RUNNER

        return command, feedback, bool(problems) or RUNNER.busy

    # -- theme and reset ---------------------------------------------------

    @callback(
        Output("mantine-provider", "forceColorScheme"),
        Input("sw-theme", "checked"),
    )
    def theme(dark):
        return "dark" if dark else "light"

    @callback(
        Output({"comp": "fb-input", "bid": ALL}, "value"),
        Output("in-group", "value"),
        Output("sw-allow-missing", "checked"),
        Output("in-normalize", "value"),
        Output("in-downsample", "value"),
        Output("in-neighbors", "value"),
        Output("in-mindist", "value"),
        Output("in-metric", "value"),
        Output("in-seed", "value"),
        Output("in-pointsize", "value"),
        Output("in-pointalpha", "value"),
        Output("in-width", "value"),
        Output("in-height", "value"),
        Output("in-formats", "value"),
        Output("in-markers", "value", allow_duplicate=True),
        Output("sw-explicit", "checked"),
        Output("sw-per-wsp", "checked"),
        Output("store-popmap", "data", allow_duplicate=True),
        Output("subpop-mode", "value"),
        Output("in-outname", "value", allow_duplicate=True),
        Input("btn-reset", "n_clicks"),
        State({"comp": "fb-input", "bid": ALL}, "value"),
        prevent_initial_call=True,
    )
    def reset_form(_clicks, path_values):
        d = config.DEFAULTS
        return (
            [""] * len(path_values),
            d["group"],
            d["allow_missing_samples"],
            d["normalize"],
            d["downsample"],
            d["n_neighbors"],
            d["min_dist"],
            d["metric"],
            d["seed"],
            d["point_size"],
            d["point_alpha"],
            d["plot_width"],
            d["plot_height"],
            list(d["output_formats"]),
            [],
            d["explicit_subpop"],
            bool(d["population_map"]),
            dict(d["population_map"]),
            d["subpop_mode"],
            "",
        )
