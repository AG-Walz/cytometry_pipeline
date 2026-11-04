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


def _loading(text: str):
    """A spinner and a few words. The log tab has the detail; this only has
    to say that something is happening."""
    return dmc.Group([dmc.Loader(size="xs"), dmc.Text(text, size="sm", c="dimmed")],
                     gap="xs")


def _failure_message(what: str, snap: dict) -> str:
    """Why a run is being reported as failed, in terms a person can act on.

    The old text printed the raw return code, which for a detached run is None
    — nobody wait()s for it — so it read "failed (exit code None)" and said
    nothing at all. What is actually useful is R's own error message. The last
    line alone is not it: after stop() R always ends on "Execution halted",
    with the error and its warnings above, so the message is taken from the
    last "Error" line up to that point.
    """
    code = snap.get("returncode")
    tail = [ln.rstrip() for ln in (snap.get("log") or "").splitlines()
            if ln.strip()]
    if tail and tail[-1].strip() == "Execution halted":
        tail = tail[:-1]
    starts = [i for i, ln in enumerate(tail) if ln.startswith("Error")]
    # R's "In addition: Warning message:" block after the error is noise here;
    # the Log tab still has it.
    error = tail[starts[-1]:] if starts else tail[-1:]
    for i, ln in enumerate(error):
        if ln.startswith("In addition:"):
            error = error[:i]
            break
    reason = (f"{what} failed (exit code {code})." if code is not None
              else f"{what} stopped before it finished.")
    if not error:
        return f"{reason} See the Log tab."
    return f"{reason}\n" + "\n".join(error[:25])


def _left_out_warning(log: str) -> str | None:
    """The "sample(s) are LEFT OUT" block R prints before extraction, or None.

    Contract with umap_from_flowjo.R (section 2a): the block starts with
    "WARNING: <n> sample(s) are LEFT OUT" and ends with the "To include them"
    sentence. A run that succeeds without some samples must not read as a run
    that used them all.
    """
    lines = log.splitlines()
    starts = [i for i, ln in enumerate(lines)
              if "sample(s) are LEFT OUT" in ln]
    if not starts:
        return None
    block = []
    for ln in lines[starts[-1]:]:
        block.append(ln.rstrip())
        if ln.startswith("To include them"):
            break
    block = block[:40]
    block[0] = block[0].replace("WARNING: ", "", 1)
    return "\n".join(block)


def _left_out_alert(text: str):
    return _alert(html.Div(text, style={"whiteSpace": "pre-wrap"}),
                  color="orange", title="Samples left out of this run")


def _format_remaining(progress: dict) -> str:
    """The time left until the WHOLE run is done, as a person would say it.

    Always the whole run, never the current step: one time, with one meaning.
    Rounded to 5 min above half an hour, because the estimate comes from one
    earlier run and a finer step would claim more than is known.
    """
    seconds = progress.get("remaining")
    if seconds is None:
        return ""
    if seconds < 90:
        return "less than 2 min left"
    step = 300 if seconds > 1800 else 60
    minutes = int(round(seconds / step) * step / 60)
    if minutes < 60:
        text = f"about {minutes} min left"
    else:
        hours, minutes = divmod(minutes, 60)
        text = (f"about {hours} h {minutes:02d} min left" if minutes
                else f"about {hours} h left")
    return ("at least " + text.removeprefix("about ")) if progress.get("overrun") \
        else text


#: One colour per step, so the steps stay apart on the bar. The done part of
#: a step is the full shade, the part still to do is a light shade of it.
_PHASE_COLOURS = {"load": "grape", "extract": "violet", "embed_nn": "indigo",
                  "embed_init": "blue", "embed": "cyan", "plot": "teal"}


def _phase_bar(progress: dict):
    """One Mantine progress bar, one section per step, as wide as the step.

    Each step fills up in its own colour, so one look tells which step runs
    and how much of the run is done. A step that is not part of this run (or
    reuses a saved result) has no width and is left out.
    """
    sections = []
    for seg in progress.get("phases") or []:
        width = seg["weight"]
        if width < 0.5:
            continue
        colour = _PHASE_COLOURS.get(seg["phase"], "blue")
        done = width * seg["fill"]
        running = seg["state"] == "running"
        label = (seg.get("short") or seg["label"]) if width >= 10 else ""
        # The name goes in the larger of the two parts, so it has space.
        in_done = done >= width - done
        if done > 0:
            sections.append(dmc.ProgressSection(
                dmc.ProgressLabel(label if in_done else ""),
                value=done, color=f"{colour}.6",
                striped=running, animated=running))
        if width - done > 0:
            sections.append(dmc.ProgressSection(
                dmc.ProgressLabel(label if not in_done else "",
                                  style={"color": "var(--mantine-color-dark-4)"}),
                value=width - done, color=f"{colour}.1"))
    if not sections:
        return dmc.Progress(value=progress.get("percent", 0), size="lg",
                            striped=True, animated=True)
    return dmc.ProgressRoot(sections, size=18)


def _progress_bar(progress: dict, running: bool):
    """The live progress bar with one line above it, or nothing.

    Renders only while something is running: a bar left at 40% after a run
    finished reads as a stalled job, and after a FAILED one it reads as a lie.
    The numbers come from app/services/progress_estimate.py.
    """
    if not running or not progress:
        return None
    step = progress.get("step")
    title = progress["phase_label"]
    if step and progress.get("steps"):
        title = f"Step {step} of {progress['steps']}: {title}"
    status = [f"{progress.get('overall_percent', 0)}%"]
    remaining = _format_remaining(progress)
    if remaining:
        status.append(remaining)
    elif not progress.get("has_reference"):
        status.append("time left known after the first complete run")
    return dmc.Stack(
        [
            dmc.Space(h=8),
            dmc.Group(
                [dmc.Text(title, size="xs", fw=600),
                 dmc.Text(" · ".join(status), size="xs", c="dimmed")],
                justify="space-between", gap="xs",
            ),
            _phase_bar(progress),
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
        f"From a load {age}. Re-read to look again." if age
        else "From an earlier load. Re-read to look again.",
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
                              f"marker(s).", size="xs"))
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
        return _loading("Loading the workspaces…"), no_update

    # -- start: run --------------------------------------------------------

    @callback(
        Output("run-feedback", "children", allow_duplicate=True),
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
        # Clears the previous run's failure, which would otherwise sit next to
        # the Run button while the new run goes.
        return None, params["output_dir"]

    @callback(
        Output("run-feedback", "children", allow_duplicate=True),
        Output({"comp": "fb-input", "bid": "results"}, "value",
               allow_duplicate=True),
        Input("btn-preview", "n_clicks"),
        *common.form_states(),
        State("store-inspect", "data"),
        prevent_initial_call=True,
    )
    def start_preview(_clicks, *values):
        """Quick preview: the same settings on ~1,000 cells per sample."""
        form = common.as_dict(values)
        inspect = values[len(common.FORM_FIELDS)] or {}
        populations = inspect.get("populations") or []
        if common.validate_run(form, inspect) or RUNNER.busy:
            return no_update, no_update
        params = common.preview_params(
            common.run_params(form, populations, inspect))
        session_state.save(common.session_snapshot(form, populations))
        ok, message = RUNNER.start_detached(
            rscript.run_argv(params), cwd=str(config.REPO_ROOT),
            output_dir=params["output_dir"], kind="run",
            meta={"output_dir": params["output_dir"], "preview": True},
        )
        if not ok:
            return _alert(f"Could not start the preview: {message}"), no_update
        return _alert(f"Preview running into {params['output_dir']} "
                      f"(about {common.PREVIEW_CELLS_PER_SAMPLE:,} cells per "
                      "sample).", color="blue"), params["output_dir"]

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
        Output({"comp": "fb-input", "bid": "results"}, "value",
               allow_duplicate=True),
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
        Output({"comp": "fb-input", "bid": "results"}, "value",
               allow_duplicate=True),
        Output("inspect-feedback", "children"),
        Output("run-feedback", "children"),
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
            _progress_bar(snap["progress"], running),
        )

        # Nothing finished since the last poll: only the live fields change.
        # A RUNNING job has not finished either, whatever the generation says:
        # starting a job does not bump it, so the first job after the app starts
        # (generation 0 against the page's -1) used to land below and be
        # reported as "stopped before it finished" for as long as it ran.
        if (snap["generation"] == seen_generation or status == "idle"
                or running):
            # Samples left out are reported by R right after the load, so the
            # warning is up while the run is still going, not only at the end.
            live = no_update
            if running and kind == "run":
                left_out = _left_out_warning(snap.get("log") or "")
                if left_out:
                    live = _left_out_alert(left_out)
            return (*base, no_update, no_update, no_update, no_update,
                    live, no_update)

        generation = snap["generation"]
        new_inspect = no_update
        new_results = no_update
        feedback = no_update
        run_feedback = no_update
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
                feedback = _alert(html.Div(
                    _failure_message("Loading the workspaces", snap),
                    style={"whiteSpace": "pre-wrap"}))

        elif kind == "run":
            output_dir = snap["meta"].get("output_dir") or ""
            new_results = output_dir
            if status == DONE:
                tab = "results"
                # The runner can take back a FAILED verdict it reached from a
                # log read too early (see Runner._refresh_detached); the alert
                # that verdict raised must go with it, or a finished run keeps
                # saying it stopped. Samples the run had to leave out are the
                # one thing a finished run still has to say.
                left_out = _left_out_warning(snap.get("log") or "")
                run_feedback = _left_out_alert(left_out) if left_out else None
            else:
                # A failed run leaves the gallery empty, and an empty gallery on
                # its own reads as "nothing happened". Say what stopped it, in
                # step 4 next to the Run button: in step 1 nobody at step 4
                # sees it.
                if status == CANCELLED:
                    run_feedback = _alert("Run cancelled.", color="orange")
                else:
                    run_feedback = _alert(html.Div(
                        _failure_message("The run", snap),
                        style={"whiteSpace": "pre-wrap"}), title="Run failed")
                tab = "log"

        return (*base, generation, new_inspect, new_results, feedback,
                run_feedback, tab)
