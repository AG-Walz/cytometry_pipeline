"""Callbacks for the form itself: navigation, gating tree, and the command."""

from __future__ import annotations

import datetime as _dt
from pathlib import Path

import dash_mantine_components as dmc
from dash import ALL, Input, Output, State, callback, ctx, html, no_update

from app import config
from app.callbacks import common
from app.services import annotations, fs_tree, gating, rscript, session_state

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
    )
    def advanced_badge(group):
        """Surface a non-default sample group on the collapsed Advanced header.

        A setting that changes which samples are loaded must not be able to
        hide behind a closed panel — that is the one way this section could
        cost someone a run.
        """
        group = str(group or "").strip()
        if not group:
            return dmc.Badge("no group set", color="red", variant="light",
                             size="sm", tt="none")
        if group == str(config.DEFAULTS["group"]):
            return None
        return dmc.Badge(f"group: {group}", color="blue", variant="light",
                         size="sm", tt="none")

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
        loaded = bool((inspect or {}).get("populations"))
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
        State("gate-tree", "selected"),
    )
    def fill_gate_tree(inspect, explicit, current_selected):
        """Sole owner of the gating tree's contents and mode.

        The tree does two jobs — pick the population, or tick a gate set — and
        both the workspace load and the explicit switch can change which. They
        are handled together so neither can overwrite the other's tree.
        """
        inspect = inspect or {}
        populations = inspect.get("populations") or []
        explicit = bool(explicit)

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
            for p in list(current_selected or [])
            + [inspect.get("remembered_population") or ""]
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
                "Every gate under the selected population. Turn this off to "
                "change the population.",
                size="xs",
                c="dimmed",
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
            # saying so here saves you hunting for a gate that only one of them
            # has. The R side refuses a population it cannot resolve in all.
            body += [
                dmc.Space(h=4),
                dmc.Text(
                    "The tree below shows only the gates shared by every "
                    "workspace — those are the ones a run can use.",
                    size="xs",
                ),
                dmc.Space(h=6),
                dmc.Group(
                    [dmc.Badge(w.get("label", "?"), variant="light", size="sm",
                               tt="none")
                     for w in workspaces],
                    gap=4,
                ),
            ]
        body += [
            dmc.Space(h=4),
            dmc.Text(", ".join(samples) or "—", size="xs", c="dimmed"),
        ]
        summary = dmc.Alert(
            body,
            color="green",
            variant="light",
            title="Workspaces loaded" if multi else "Workspace loaded",
        )
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

    @callback(
        Output("gate-tree", "expanded", allow_duplicate=True),
        Input("btn-expand-all", "n_clicks"),
        Input("btn-collapse", "n_clicks"),
        State("store-inspect", "data"),
        prevent_initial_call=True,
    )
    def expand_collapse(_expand, _collapse, inspect):
        if ctx.triggered_id == "btn-collapse":
            return []
        return gating.all_values((inspect or {}).get("populations") or [])

    @callback(
        Output("population-echo", "children"),
        Input("gate-tree", "selected"),
    )
    def echo_population(selected):
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
            return _alert(
                "Tick gates underneath the selected population to build the "
                "colour partition.",
                color="yellow",
            )
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
        columns, preview = annotations.describe(
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
                "SVG previews of a large embedding carry one node per cell and "
                "can be slow to render in the browser. Tick PNG too for a fast "
                "gallery.",
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

        params = common.run_params(form, populations)
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
            d["subpop_mode"],
            "",
        )
