"""Callbacks that start, stop and report on the R subprocess.

One ``dcc.Interval`` drives everything the UI needs while something is running:
the status badge, the elapsed clock, the log, and — when the process exits — the
follow-up work (parsing an inspect result, pointing the gallery at a finished
run). Keeping that in a single callback means there is exactly one place where a
completed run is turned into UI state, so a completion can never be handled
twice or missed.
"""

from __future__ import annotations

import dash_mantine_components as dmc
from dash import Input, Output, State, callback, html, no_update

from app import config
from app.callbacks import common
from app.services import rscript, session_state
from app.services.runner import (
    CANCELLED,
    DONE,
    FAILED,
    RUNNING,
    RUNNER,
    format_elapsed,
)

STATUS_COLOURS = {
    "idle": "gray",
    RUNNING: "blue",
    DONE: "green",
    FAILED: "red",
    CANCELLED: "orange",
}


def _alert(children, color="red", title=None):
    return dmc.Alert(children, color=color, title=title, variant="light")


def register_callbacks() -> None:
    # -- start: inspect ----------------------------------------------------

    @callback(
        Output("inspect-feedback", "children", allow_duplicate=True),
        Input("btn-inspect", "n_clicks"),
        State({"comp": "fb-input", "bid": "wsp"}, "value"),
        State({"comp": "fb-input", "bid": "fcs"}, "value"),
        State("in-group", "value"),
        State("sw-allow-missing", "checked"),
        prevent_initial_call=True,
    )
    def start_inspect(_clicks, wsp, fcs_dir, group, allow_missing):
        form = {"wsp": wsp, "fcs_dir": fcs_dir, "group": group}
        problems = common.validate_inspect(form)
        if problems:
            return _alert(
                html.Ul([html.Li(p) for p in problems],
                        style={"margin": 0, "paddingLeft": "18px"})
            )
        if RUNNER.busy:
            return _alert(
                "Something is already running — wait for it or cancel it first.",
                color="yellow",
            )

        argv = rscript.inspect_argv(wsp.strip(), fcs_dir.strip(), group,
                                    allow_missing_samples=bool(allow_missing))
        ok, message = RUNNER.start(
            argv, cwd=str(config.REPO_ROOT), kind="inspect",
            meta={"wsp": wsp, "fcs_dir": fcs_dir, "group": group},
        )
        if not ok:
            return _alert(f"Could not start R: {message}")
        return _alert("Loading — watch the Log tab.", color="blue")

    # -- start: run --------------------------------------------------------

    @callback(
        Output("run-start-sink", "children"),
        Output({"comp": "fb-input", "bid": "results"}, "value",
               allow_duplicate=True),
        Input("btn-run", "n_clicks"),
        *common.form_states(),
        State("store-inspect", "data"),
        prevent_initial_call=True,
    )
    def start_run(_clicks, *values):
        form = common.as_dict(values)
        inspect = values[len(common.FORM_FIELDS)] or {}
        populations = inspect.get("populations") or []

        if common.validate_run(form, inspect) or RUNNER.busy:
            return no_update, no_update

        params = common.run_params(form, populations)
        session_state.save(common.session_snapshot(form, populations))

        argv = rscript.run_argv(params)
        RUNNER.start(
            argv, cwd=str(config.REPO_ROOT), kind="run",
            meta={"output_dir": params["output_dir"]},
        )
        return "", params["output_dir"]

    @callback(
        Output("run-start-sink", "children", allow_duplicate=True),
        Input("btn-cancel", "n_clicks"),
        prevent_initial_call=True,
    )
    def cancel(_clicks):
        RUNNER.cancel()
        return ""

    # -- poll --------------------------------------------------------------

    @callback(
        Output("status-badge", "children"),
        Output("status-badge", "color"),
        Output("elapsed-text", "children"),
        Output("log-text", "children"),
        Output("btn-cancel", "disabled"),
        Output("btn-inspect", "loading"),
        Output("store-run-generation", "data"),
        Output("store-inspect", "data"),
        Output({"comp": "fb-input", "bid": "results"}, "value",
               allow_duplicate=True),
        Output("inspect-feedback", "children"),
        Output("tabs", "value"),
        Input("poll", "n_intervals"),
        State("store-run-generation", "data"),
        State("store-inspect", "data"),
        State("tabs", "value"),
        # The results-folder output is shared with start_run (and the file
        # browser); this is what lets a duplicate output fire on page load.
        prevent_initial_call="initial_duplicate",
    )
    def poll(_n, seen_generation, inspect_store, current_tab):
        snap = RUNNER.snapshot()
        status = snap["status"]
        kind = snap["kind"]
        running = status == RUNNING

        label = f"{status}" + (f" · {kind}" if kind and status != "idle" else "")
        elapsed = format_elapsed(snap["elapsed"])
        log = snap["log"] or "Nothing has run yet."

        base = (
            label,
            STATUS_COLOURS.get(status, "gray"),
            elapsed,
            log,
            not running,
            running and kind == "inspect",
        )

        # Nothing finished since the last poll: only the live fields change.
        if snap["generation"] == seen_generation or status == "idle":
            return (*base, no_update, no_update, no_update, no_update, no_update)

        generation = snap["generation"]
        new_inspect = no_update
        new_results = no_update
        feedback = no_update
        tab = no_update

        if kind == "inspect":
            if status == DONE:
                parsed = rscript.parse_inspect(snap["log"])
                if not parsed["populations"]:
                    n_wsp = len(parsed.get("workspaces") or [])
                    feedback = _alert(
                        "The workspaces loaded but no usable gate was "
                        "reported — with several workspaces the tree is the "
                        "gates they all share, so this means their gating "
                        "trees have nothing in common. See the Log tab."
                        if n_wsp > 1 else
                        "The workspace loaded but no gating populations were "
                        "reported. Check the Log tab and the sample group."
                    )
                else:
                    remembered = session_state.load().get("population") or ""
                    new_inspect = {
                        **parsed,
                        "wsp": snap["meta"].get("wsp"),
                        "fcs_dir": snap["meta"].get("fcs_dir"),
                        "group": snap["meta"].get("group"),
                        "remembered_population": remembered,
                    }
                    note = ""
                    if not parsed["markers"]:
                        note = (
                            " No marker is common to all samples — the marker "
                            "picker will be empty."
                        )
                    n_wsp = len(parsed.get("workspaces") or [])
                    where = f" across {n_wsp} workspaces" if n_wsp > 1 else ""
                    feedback = _alert(
                        f"Loaded {len(parsed['samples'])} sample(s){where} and "
                        f"{len(parsed['populations'])} gate(s)."
                        f"{note} Continue to step 2.",
                        color="green",
                    )
            elif status == CANCELLED:
                feedback = _alert("Workspace load cancelled.", color="orange")
            else:
                feedback = _alert(
                    f"Loading the workspace failed (exit code "
                    f"{snap['returncode']}). See the Log tab.",
                )

        elif kind == "run":
            output_dir = snap["meta"].get("output_dir") or ""
            new_results = output_dir
            if status == DONE:
                tab = "results"
            elif current_tab == "results":
                # A failed run leaves the gallery empty; the log explains why.
                tab = "log"

        return (*base, generation, new_inspect, new_results, feedback, tab)
