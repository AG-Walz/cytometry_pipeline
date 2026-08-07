"""Callbacks that start, stop and report on the R subprocess.

One ``dcc.Interval`` drives everything the UI needs while something is running:
the status badge, the elapsed clock, the log, and — when the process exits — the
follow-up work (parsing an inspect result, pointing the gallery at a finished
run). Keeping that in a single callback means there is exactly one place where a
completed run is turned into UI state, so a completion can never be handled
twice or missed.
"""

from __future__ import annotations

import time
from pathlib import Path

import dash_mantine_components as dmc
from dash import Input, Output, State, callback, ctx, html, no_update

from app import config
from app.callbacks import common
from app.services import detached, report_cache, rscript, session_state
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


def _failure_message(what: str, snap: dict) -> str:
    """Why a run is being reported as failed, in terms a person can act on.

    The old text printed the raw return code, which for a detached run is None
    — nobody wait()s for it — so it read "failed (exit code None)" and said
    nothing at all. What is actually useful is the last thing the run printed,
    because that is where it stopped.
    """
    code = snap.get("returncode")
    tail = [ln for ln in (snap.get("log") or "").splitlines() if ln.strip()]
    last = tail[-1].strip() if tail else ""
    reason = (f"{what} failed (exit code {code})." if code is not None
              else f"{what} stopped before it finished.")
    return f"{reason} Last line: {last}" if last else f"{reason} See the Log tab."


def _format_remaining(seconds: float | None) -> str:
    """A duration a person would say out loud, not a stopwatch reading.

    Deliberately coarse. The estimate is a linear extrapolation over uneven
    work, so "about 25 min" is the most it can honestly claim; printing 24:37
    would dress a guess up as a measurement.
    """
    if seconds is None:
        return ""
    if seconds < 60:
        return "under a minute left"
    minutes = int(round(seconds / 60))
    if minutes < 60:
        return f"about {minutes} min left"
    hours, minutes = divmod(minutes, 60)
    return f"about {hours} h {minutes:02d} min left" if minutes else \
        f"about {hours} h left"


def _phase_bar(progress: dict):
    """One bar per phase, side by side, each as wide as that phase costs.

    A single bar over a six-phase run answers the wrong question. It says how
    far through the CURRENT phase you are, so the 20-minute neighbour search and
    the 20-second initialisation look the same size, and the bar resets to empty
    five times before the run ends — which reads as five restarts. Laying the
    phases out at their real relative cost makes one glance answer both
    questions people actually have: what is happening, and how much is left.

    Every number comes from Runner._phase_view(); this only arranges them.
    """
    segments = progress.get("phases") or []
    if not segments:
        return dmc.Progress(value=progress.get("percent", 0), size="sm",
                            animated=True, striped=True)

    # A done phase is solid, the running one is striped and animated, a pending
    # one is an empty track. The distinction has to survive a screenshot, so it
    # is carried by fill and colour rather than by motion alone.
    bars = []
    for seg in segments:
        state = seg["state"]
        if state == "done":
            inner = dmc.Progress(value=100, size="sm", color="blue")
        elif state == "running":
            inner = dmc.Progress(value=round(seg["fill"] * 100), size="sm",
                                 color="blue", animated=True, striped=True)
        elif state == "skipped":
            inner = dmc.Progress(value=0, size="sm", color="gray")
        else:
            inner = dmc.Progress(value=0, size="sm", color="gray")
        bars.append(
            dmc.Tooltip(
                dmc.Box(inner, style={"flex": f"{seg['weight']} 0 0",
                                      "minWidth": "6px"}),
                label=f"{seg['label']} — {round(seg['fill'] * 100)}%",
                position="top", withArrow=True, openDelay=200,
            )
        )
    return dmc.Stack(
        [
            dmc.Group(bars, gap=3, wrap="nowrap", align="stretch"),
            # The whole-run figure, said once. Without it the segments show the
            # shape of the run but never its total, and "how far along am I"
            # goes back to being a guess.
            dmc.Text(f"{progress.get('overall_percent', 0)}% of the run",
                     size="xs", c="dimmed"),
        ],
        gap=4,
    )


def _progress_bar(progress: dict, running: bool):
    """The live progress bar, or nothing.

    Renders only while something is running: a bar left at 40% after a run
    finished reads as a stalled job, and after a FAILED one it reads as a lie.
    Everything shown here was counted by R (see emit_progress in
    src/R_shared/io_utils.R); this arranges it and works nothing out except the
    wording of the estimate.
    """
    if not running or not progress:
        return None
    remaining = _format_remaining(progress.get("remaining"))
    return dmc.Stack(
        [
            dmc.Space(h=8),
            dmc.Group(
                [
                    dmc.Text(progress["phase_label"], size="xs", fw=600),
                    dmc.Text(
                        f"{progress['done']} / {progress['total']}"
                        + (f" · {remaining}" if remaining else ""),
                        size="xs", c="dimmed"),
                ],
                justify="space-between", gap="xs",
            ),
            _phase_bar(progress),
            # The item in flight, so a run that has stopped moving is obvious
            # from WHICH item it stopped on.
            dmc.Text(progress.get("label") or "", size="xs", c="dimmed",
                     truncate="end"),
        ],
        gap=4,
    )


def _cached_feedback(report: dict, parsed: dict, cache_key: str):
    """The same summary as a fresh load, plus where the answer came from.

    Said explicitly, because a load that returns instantly is otherwise
    indistinguishable from one that silently did nothing — and because a user
    who has just changed something outside the key (added an FCS file) needs to
    know this did not look at it.
    """
    base = _inspect_feedback(report, parsed)
    age = report_cache.age_text(cache_key)
    note = dmc.Text(
        (f"Remembered from a load {age} — use Re-read to look again."
         if age else "Remembered from an earlier load — use Re-read to look "
                     "again."),
        size="xs", c="dimmed")
    return dmc.Stack([base, note], gap=4)


def _inspect_feedback(report: dict, parsed: dict):
    """One line about the load, pointing at whatever needs attention.

    Reports the load as it actually happened rather than as pass/fail: a load
    where twenty of thirty-three workspaces came in is neither, and calling it a
    success hides the ten that did not while calling it a failure hides the
    twenty that did.
    """
    totals = report.get("totals") or {}
    found = totals.get("workspaces_found") or 0
    loaded = totals.get("workspaces_loaded") or 0
    samples = totals.get("samples") or 0
    gates = len(parsed.get("populations") or [])
    markers = len(parsed.get("markers") or [])

    if not loaded:
        return _alert(
            "No workspace loaded. Open Data health below — it names each one "
            "and why."
        )

    lines = [
        dmc.Text(f"Loaded {samples} sample(s) from {loaded} of {found} "
                 f"workspace(s).", fw=600, size="sm"),
    ]
    colour = "green"
    if loaded < found:
        colour = "yellow"
        lines.append(dmc.Text(
            f"{found - loaded} workspace(s) are not in the analysis — Data "
            f"health names each one and why.", size="xs"))
    if not gates:
        colour = "yellow"
        lines.append(dmc.Text(
            "No gate is shared by every workspace, so there is nothing to pick "
            "in step 2. Choose the gate per workspace there instead.",
            size="xs"))
    if not markers:
        colour = "yellow"
        lines.append(dmc.Text(
            "No marker is shared by every workspace. Data health shows which "
            "spellings differ and can suggest aliases.", size="xs"))
    if gates and markers and loaded == found:
        lines.append(dmc.Text(f"{gates} shared gate(s), {markers} shared "
                              f"marker(s). Continue to step 2.", size="xs"))
    return _alert(dmc.Stack(lines, gap=4), color=colour)


def register_callbacks() -> None:
    # -- start: inspect ----------------------------------------------------

    @callback(
        Output("inspect-feedback", "children", allow_duplicate=True),
        Output("store-inspect", "data", allow_duplicate=True),
        Input("btn-inspect", "n_clicks"),
        Input("btn-reload", "n_clicks"),
        State({"comp": "fb-input", "bid": "wsp"}, "value"),
        State({"comp": "fb-input", "bid": "fcs"}, "value"),
        State("in-group", "value"),
        State("sw-allow-missing", "checked"),
        State("in-match-mode", "value"),
        State("in-marker-aliases", "value"),
        State("sw-recover-markers", "checked"),
        prevent_initial_call=True,
    )
    def start_inspect(_clicks, _reload_clicks, wsp, fcs_dir, group,
                      allow_missing, match_mode, marker_aliases,
                      recover_markers):
        """Load the workspaces — from the remembered result where possible.

        "Reload" is the same load with the cache bypassed, handled here rather
        than in a second callback so there is one place that decides what a
        load is.

        A load reads EVERY workspace under the folder; the health panel's on/off
        switches are not consulted here. They used to be, and backwards: the
        stored set is what is switched OFF while R's --workspaces is what to
        switch ON, so one workspace switched off made the next load read that
        workspace and nothing else. Everything downstream is keyed on the
        workspaces the report contains, so the whole form — the shared gating
        tree, the panel, the per-workspace populations — collapsed onto the
        rejects, and it was self-sustaining: the panel could then only redraw
        the handful it had loaded. The selection belongs to the run, which
        converts it properly (see common.workspace_selection).
        """
        force = ctx.triggered_id == "btn-reload"
        form = {"wsp": wsp, "fcs_dir": fcs_dir, "group": group}
        problems = common.validate_inspect(form)
        if problems:
            return _alert(
                html.Ul([html.Li(p) for p in problems],
                        style={"margin": 0, "paddingLeft": "18px"})
            ), no_update
        if RUNNER.busy:
            return _alert(
                "Something is already running — wait for it or cancel it first.",
                color="yellow",
            ), no_update

        request = {
            "wsp": wsp.strip(), "fcs_dir": fcs_dir.strip(),
            "group": str(group).strip(),
            "match_mode": match_mode or config.DEFAULTS["match_mode"],
            "marker_aliases": marker_aliases or "",
            "recover_markers": bool(recover_markers),
            "allow_missing_samples": bool(allow_missing),
        }
        cache_key = report_cache.key_for(request)

        # A previous load of exactly this, still valid. Answered without
        # starting R at all: ten minutes of GatingSet building is the single
        # worst thing to repeat for no reason.
        if not force:
            cached = report_cache.load(cache_key)
            if cached is not None:
                try:
                    parsed = rscript.inspect_from_report(cached)
                except rscript.ReportError:
                    parsed = None
                if parsed is not None:
                    store = {
                        **parsed,
                        "wsp": wsp, "fcs_dir": fcs_dir, "group": group,
                        "remembered_population":
                            session_state.load().get("population") or "",
                        "cache_key": cache_key,
                    }
                    return _cached_feedback(cached, parsed, cache_key), store

        report_json = report_cache.reserve(cache_key)
        argv = rscript.inspect_argv(
            wsp.strip(), fcs_dir.strip(), group,
            allow_missing_samples=bool(allow_missing),
            match_mode=match_mode or config.DEFAULTS["match_mode"],
            marker_aliases=marker_aliases or "",
            recover_markers=bool(recover_markers),
            report_json=report_json,
        )
        ok, message = RUNNER.start(
            argv, cwd=str(config.REPO_ROOT), kind="inspect",
            meta={"wsp": wsp, "fcs_dir": fcs_dir, "group": group,
                  "report_json": report_json, "cache_key": cache_key},
        )
        if not ok:
            return _alert(f"Could not start R: {message}"), no_update
        return _alert(
            "Reading the workspaces. The result is remembered, so this is a "
            "one-off for these folders.",
            color="blue",
        ), no_update

    # -- start: run --------------------------------------------------------

    @callback(
        Output("run-start-sink", "children"),
        Output("store-results-dir", "data", allow_duplicate=True),
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

        params = common.run_params(form, populations, inspect)
        session_state.save(common.session_snapshot(form, populations))

        argv = rscript.run_argv(params)
        # Detached, always. A run over millions of cells is an overnight job:
        # it has to survive the browser closing, the app being restarted and the
        # app being redeployed. Its log goes into the run's own output directory
        # so it can be re-read afterwards — by this app, by the next one, or by
        # someone with tail -f. See app/services/detached.py.
        ok, message = RUNNER.start_detached(
            argv, cwd=str(config.REPO_ROOT),
            output_dir=params["output_dir"], kind="run",
            meta={"output_dir": params["output_dir"]},
        )
        if not ok:
            return _alert(f"Could not start the run: {message}"), no_update
        return "", params["output_dir"]

    @callback(
        Output("run-start-sink", "children", allow_duplicate=True),
        Input("btn-cancel", "n_clicks"),
        prevent_initial_call=True,
    )
    def cancel(_clicks):
        RUNNER.cancel()
        return ""

    # -- reattach ----------------------------------------------------------

    @callback(
        Output("reattach-note", "children"),
        Output("store-results-dir", "data", allow_duplicate=True),
        Input("poll", "n_intervals"),
        State({"comp": "fb-input", "bid": "outparent"}, "value"),
        prevent_initial_call=True,
    )
    def reattach(n_intervals, output_parent):
        """Pick a still-running detached job back up after the app restarts.

        Runs once, on the first poll. An overnight job that survives a restart
        is only useful if the app can find it again in the morning — otherwise
        it is a process nobody can see, writing into a directory nobody is
        watching.

        Only the newest live run is attached, and only when nothing else is
        going, so this can never displace what the user is currently watching.
        """
        if n_intervals != 1 or RUNNER.busy or not (output_parent or "").strip():
            return no_update, no_update

        live = [s for s in detached.find_recent([output_parent.strip()])
                if s.get("running")]
        if not live:
            return no_update, no_update

        newest = live[0]
        if not RUNNER.attach(newest["output_dir"]):
            return no_update, no_update
        started = time.strftime("%H:%M", time.localtime(newest["started_at"]))
        return _alert(
            f"Picked up a run started at {started} that is still going — "
            f"{Path(newest['output_dir']).name}.",
            color="blue",
        ), newest["output_dir"]

    # -- poll --------------------------------------------------------------

    @callback(
        Output("status-badge", "children"),
        Output("status-badge", "color"),
        Output("elapsed-text", "children"),
        Output("log-text", "children"),
        Output("btn-cancel", "disabled"),
        Output("btn-inspect", "loading"),
        Output("progress-wrap", "children"),
        Output("store-run-generation", "data"),
        Output("store-inspect", "data"),
        Output("store-results-dir", "data"),
        Output("inspect-feedback", "children"),
        Output("tabs", "value"),
        Input("poll", "n_intervals"),
        State("store-run-generation", "data"),
        State("store-inspect", "data"),
        State("tabs", "value"),
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
            _progress_bar(snap["progress"], running),
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
                try:
                    report = rscript.read_report(
                        snap["meta"].get("report_json") or "")
                except rscript.ReportError as exc:
                    feedback = _alert(str(exc))
                    report = None

                if report is not None:
                    parsed = rscript.inspect_from_report(report)
                    remembered = session_state.load().get("population") or ""
                    # Stored even when nothing is usable: the health panel is
                    # most needed exactly then, and the old behaviour — throw
                    # the result away and print one line about the Log tab —
                    # left the user with 862 lines of prose and no next step.
                    new_inspect = {
                        **parsed,
                        "wsp": snap["meta"].get("wsp"),
                        "fcs_dir": snap["meta"].get("fcs_dir"),
                        "group": snap["meta"].get("group"),
                        "remembered_population": remembered,
                    }
                    feedback = _inspect_feedback(report, parsed)
            elif status == CANCELLED:
                feedback = _alert("Workspace load cancelled.", color="orange")
            else:
                feedback = _alert(_failure_message("Loading the workspaces",
                                                   snap))

        elif kind == "run":
            output_dir = snap["meta"].get("output_dir") or ""
            new_results = output_dir
            if status == DONE:
                tab = "results"
            else:
                # A failed run leaves the gallery empty, and an empty gallery on
                # its own reads as "nothing happened". Say what stopped it.
                if status == CANCELLED:
                    feedback = _alert("Run cancelled.", color="orange")
                else:
                    feedback = _alert(_failure_message("The run", snap))
                tab = "log"

        return (*base, generation, new_inspect, new_results, feedback, tab)
