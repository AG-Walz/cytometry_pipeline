"""The version badge in the header: what code this app runs, and whether it
is the newest. See app/services/version.py."""

from __future__ import annotations

import dash_mantine_components as dmc
from dash import Input, Output, callback

from app.services import version

_COLOUR = {"current": "green", "restart": "orange", "behind": "orange",
           "unknown": "gray"}


def badge(info: dict | None = None):
    """The badge and its tooltip for ``version.status()``."""
    info = info or version.status()
    run = info["running"]
    short = run.get("short") or "?"
    date = (run.get("date") or "")[:10]
    state = info["state"]
    text = {
        "current": f"Version {short} · {date} · up to date",
        "restart": f"Version {short} · newer code on the server: restart the app",
        "behind": f"Version {short} · a newer version exists",
        "unknown": f"Version {short} · {date}",
    }[state]
    lines = [
        f"Running: {short} ({run.get('date', '?')}) {run.get('subject', '')}",
        f"App started: {run.get('started', '?')}",
    ]
    if run.get("dirty"):
        lines.append("The checkout had local changes when the app started.")
    if state == "restart":
        disk = info["on_disk"]
        lines.append(f"On the server now: {disk.get('short')} ({disk.get('date')}) "
                     f"{disk.get('subject', '')}. Close the app and open it "
                     "again to load it.")
    elif state == "behind":
        lines.append(f"GitHub main is at {(info['remote_commit'] or '')[:7]}. "
                     "Opening the app again from the launcher pulls it.")
    elif info.get("remote_error"):
        lines.append(f"Not compared with GitHub: {info['remote_error']}.")
    return dmc.Tooltip(
        dmc.Badge(text, color=_COLOUR[state], variant="light", size="sm",
                  style={"textTransform": "none"}),
        label=dmc.Stack([dmc.Text(l, size="xs") for l in lines], gap=2),
        multiline=True, w=380, withArrow=True, position="bottom-start")


def register_callbacks() -> None:
    @callback(Output("version-badge", "children"),
              Input("version-poll", "n_intervals"))
    def refresh(_n):
        return badge()
