"""Callbacks for the "FCS names" tab.

This tab is deliberately thin. It owns a form and a keyword table, and hands
everything else to machinery that already exists:

* the run itself goes through the shared ``RUNNER`` with ``kind="fcs-names"``,
  and the poll callback in ``execution.py`` already drives the status badge, the
  elapsed clock, the log and the Cancel button for *any* kind — so the report
  streams into the Log pane with no changes there at all;
* the folder pickers are ``browser_field``s, whose callbacks are
  pattern-matching and therefore already cover them;
* ``RUNNER`` being single-slot is what stops a sync and a UMAP run overlapping.

There is no parser here. ``sync_fcs_names.R``'s output *is* the deliverable and
is shown verbatim; the one structured thing the tab needs — the keyword table —
comes from R as a CSV, because R already owns the only correct FCS parser in
this repo and a second one in Python would drift from it.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path

import dash_mantine_components as dmc
from dash import Input, Output, State, callback, html, no_update

from app import config
from app.services import rscript, runner, session_state

#: The keyword read is synchronous, unlike a sync run: it only touches each
#: file's header and TEXT segment (a few KB), so even a large cohort returns in
#: seconds, and doing it through RUNNER would mean a new poll branch and a
#: completion handler for something the user is waiting on anyway.
KEYWORD_TIMEOUT_S = 180

#: Sticky first column. The table is keywords x files, and a cohort is wide, so
#: the keyword name has to stay put while you scroll sideways through samples.
_STICKY = {
    "position": "sticky",
    "left": 0,
    "zIndex": 2,
    "background": "var(--mantine-color-body)",
    "borderRight": "1px solid var(--mantine-color-default-border)",
    "whiteSpace": "nowrap",
}


def _alert(children, color="red", title=None):
    return dmc.Alert(children, color=color, title=title, variant="light")


def _folder_problem(folder: str) -> str:
    folder = (folder or "").strip()
    if not folder:
        return "Choose the folder holding your renamed FCS files."
    if not Path(folder).is_dir():
        return f"Not a folder: {folder}"
    return ""


def register_callbacks() -> None:
    # -- the command preview, and what the button says ---------------------
    #
    # One callback for all three, because they are three views of the same
    # thing: the mode decides the flags, the label and the colour. Deliberately
    # NOT driven by the poll interval — it would re-render a dmc.Code twice a
    # second forever to grey out a button that RUNNER already refuses.

    @callback(
        Output("fcs-cmd", "children"),
        Output("fcs-note", "children"),
        Output("btn-fcs", "children"),
        Output("btn-fcs", "color"),
        Input({"comp": "fb-input", "bid": "fcssync"}, "value"),
        Input("fcs-mode", "value"),
        Input("fcs-also", "value"),
    )
    def preview(folder, mode, also):
        spec = rscript.SYNC_MODES.get(mode) or rscript.SYNC_MODES["dry"]
        argv = rscript.sync_argv((folder or "").strip(), mode or "dry", also or "")
        problem = _folder_problem(folder)
        if problem:
            note = _alert(problem, color="yellow")
        elif spec.get("writes") and mode != "revert":
            note = _alert(
                "This writes to your files. The name being replaced is kept in "
                "ORIGINALFIL, and every change is logged to "
                "fcs_name_sync_ledger.csv next to the data — Undo replays it.",
                color="orange",
            )
        else:
            note = None
        return rscript.command_string(argv), note, spec["button"], spec["color"]

    # -- start a run -------------------------------------------------------

    @callback(
        Output("fcs-feedback", "children"),
        Output("tabs", "value", allow_duplicate=True),
        Input("btn-fcs", "n_clicks"),
        State({"comp": "fb-input", "bid": "fcssync"}, "value"),
        State("fcs-mode", "value"),
        State("fcs-also", "value"),
        prevent_initial_call=True,
    )
    def start_sync(_clicks, folder, mode, also):
        problem = _folder_problem(folder)
        if problem:
            return _alert(problem), no_update

        # The folder only. The mode is deliberately not remembered: a refresh
        # that restored "Apply" would turn one stray click into a cohort-wide
        # write.
        session_state.save({"fcs_sync_dir": folder.strip()})

        argv = rscript.sync_argv(folder.strip(), mode, also or "")
        ok, message = runner.RUNNER.start(
            argv, cwd=str(config.REPO_ROOT), kind="fcs-names",
            meta={"fcs_dir": folder.strip(), "mode": mode},
        )
        if not ok:
            return _alert(message, color="yellow"), no_update
        # Switch the right-hand pane to the log: the report is the whole output.
        return _alert("Running — watch the Log.", color="blue"), "log"

    # -- the stored metadata table -----------------------------------------

    @callback(
        Output("store-fcs-kw", "data"),
        Output("fcs-feedback", "children", allow_duplicate=True),
        Input("btn-fcs-kw", "n_clicks"),
        State({"comp": "fb-input", "bid": "fcssync"}, "value"),
        prevent_initial_call=True,
    )
    def load_keywords(_clicks, folder):
        problem = _folder_problem(folder)
        if problem:
            return no_update, _alert(problem)
        # Reading a file that is being rewritten would give a torn view. The
        # writer is the only thing that can be doing that, and it is the runner.
        if runner.RUNNER.busy:
            return no_update, _alert(
                "Something is running — wait for it to finish, or cancel it.",
                color="yellow",
            )

        with tempfile.TemporaryDirectory() as tmp:
            csv_path = str(Path(tmp) / "keywords.csv")
            argv = rscript.keyword_table_argv(folder.strip(), csv_path)
            try:
                proc = subprocess.run(
                    argv, cwd=str(config.REPO_ROOT), capture_output=True,
                    text=True, timeout=KEYWORD_TIMEOUT_S,
                )
            except subprocess.TimeoutExpired:
                return no_update, _alert(
                    f"Reading the folder took longer than {KEYWORD_TIMEOUT_S}s "
                    "and was stopped."
                )
            except OSError as exc:
                return no_update, _alert(f"Could not start R: {exc}")

            if not Path(csv_path).exists():
                tail = "\n".join((proc.stderr or proc.stdout or "").splitlines()[-6:])
                return no_update, _alert(
                    html.Div([
                        html.Div("R did not produce a keyword table."),
                        html.Pre(tail, style={"whiteSpace": "pre-wrap",
                                              "fontSize": "11px", "margin": 0}),
                    ])
                )
            with open(csv_path, newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))

        if not rows:
            return no_update, _alert("No readable FCS files in that folder.",
                                     color="yellow")
        return rows, None

    @callback(
        Output("fcs-kw-table", "children"),
        Input("store-fcs-kw", "data"),
        Input("sw-fcs-allkw", "checked"),
    )
    def render_keywords(rows, show_all):
        if not rows:
            return None
        return _keyword_table(rows, bool(show_all))


def _keyword_table(rows: list[dict], show_all: bool):
    """Keywords down, files across, keyword column frozen.

    Rows first so a value shared by the whole cohort reads as one horizontal
    band — which is exactly the symptom this tab exists to fix. R supplies
    ``keyword_rank`` (1..n for the identity keywords, 999 for the rest) so the
    ordering does not have to be duplicated here.
    """
    files = sorted({r["file"] for r in rows})
    ranks: dict[str, int] = {}
    name_like: dict[str, bool] = {}
    for r in rows:
        try:
            ranks[r["keyword"]] = int(r["keyword_rank"])
        except (TypeError, ValueError):
            ranks[r["keyword"]] = 999
        name_like[r["keyword"]] = str(r.get("name_like", "")).upper() == "TRUE"
    keywords = [k for k, rank in ranks.items() if show_all or rank < 999]
    keywords.sort(key=lambda k: (ranks[k], k))
    lookup = {(r["file"], r["keyword"]): r["value"] for r in rows}

    header = dmc.TableTr(
        [dmc.TableTh("keyword", style={**_STICKY, "zIndex": 3})]
        + [dmc.TableTh(Path(f).name, style={"whiteSpace": "nowrap"})
           for f in files]
    )

    body = []
    for kw in keywords:
        values = [lookup.get((f, kw)) for f in files]
        # Only flag duplicates among the NAME-LIKE keywords, which R marks for
        # us. A cohort sharing $TOT, $PAR or $DATE is ordinary and highlighting
        # it would drown the signal; a cohort sharing $FIL or GUID is the exact
        # thing this tab exists to fix.
        shared: set[str] = set()
        if name_like.get(kw):
            seen: dict[str, int] = {}
            for v in values:
                if v:
                    seen[v] = seen.get(v, 0) + 1
            shared = {v for v, n in seen.items() if n > 1}

        cells = [dmc.TableTd(kw, style={**_STICKY, "fontWeight": 500})]
        for v in values:
            style = {"whiteSpace": "nowrap"}
            if v is None:
                v = "—"
                style["color"] = "var(--mantine-color-dimmed)"
            elif v in shared:
                style["background"] = "var(--mantine-color-orange-light)"
            cells.append(dmc.TableTd(v, style=style))
        body.append(dmc.TableTr(cells))

    return html.Div([
        dmc.Group(
            [
                dmc.Text(f"{len(keywords)} keyword(s) x {len(files)} file(s)",
                         size="xs", c="dimmed"),
                dmc.Text("shared name — the problem", size="xs", c="dimmed",
                         style={"background": "var(--mantine-color-orange-light)",
                                "padding": "0 6px", "borderRadius": "3px"}),
            ],
            justify="space-between",
        ),
        dmc.Space(h=6),
        dmc.ScrollArea(
            dmc.Table(
                [dmc.TableThead(header), dmc.TableTbody(body)],
                striped=True, highlightOnHover=True, withTableBorder=True,
                horizontalSpacing="sm", verticalSpacing=4, fz="xs",
            ),
            h=360, type="auto", offsetScrollbars=True,
        ),
    ])
