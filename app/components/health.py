"""Render the R load report as a per-workspace data-health table.

Presentation only. Every number, verdict and diagnosis shown here was computed
by ``src/visualization/umap_load_report.R``; this module decides what is worth
looking at first and nothing else. If you find yourself wanting to work
something out in here, it belongs in R — the app cannot see a GatingSet, so
anything it derives is a second opinion that will eventually disagree with the
first.

Why a table at all: a 33-workspace load prints an 862-line report. Prose is the
right shape for one workspace and the wrong shape for thirty. The questions
someone actually has — *is my data in there, and if not, which part isn't* — are
per-workspace, and a table answers them at a glance.
"""

from __future__ import annotations

import dash_mantine_components as dmc
from dash import html

#: status -> (badge colour, what it means in one phrase)
_STATUS = {
    "loaded": ("green", "loaded"),
    "build_failed": ("red", "CytoML refused it"),
    "unusable": ("red", "could not be read"),
    "no_files": ("orange", "no files found"),
    "not_loaded": ("gray", "not loaded"),
}


def _pill(label: str, value, color: str = "gray") -> dmc.Badge:
    return dmc.Badge(f"{label} {value}", color=color, variant="light", size="sm")


def _sample_list(title: str, names: list[str], color: str, limit: int = 40):
    """A collapsible list of sample names, or nothing when there are none.

    Capped, and the cap is stated. A truncation that does not say it truncated
    reads as a complete list, and someone will act on it as one.
    """
    if not names:
        return None
    shown = names[:limit]
    extra = len(names) - len(shown)
    return dmc.Accordion(
        chevronPosition="left",
        variant="contained",
        children=[
            dmc.AccordionItem(
                value=title,
                children=[
                    dmc.AccordionControl(
                        dmc.Text(f"{title} ({len(names)})", size="xs",
                                 c=color, fw=600)),
                    dmc.AccordionPanel(
                        dmc.Stack(
                            [dmc.Text(n, size="xs", ff="monospace")
                             for n in shown]
                            + ([dmc.Text(f"… and {extra} more — see the Log tab",
                                         size="xs", c="dimmed", fs="italic")]
                               if extra else []),
                            gap=2)),
                ],
            )
        ],
    )


def _workspace_row(ws: dict) -> dmc.Paper:
    color, meaning = _STATUS.get(ws.get("status", ""), ("gray", ws.get("status", "")))
    routes = ws.get("routes") or {}
    n_loaded = ws.get("n_loaded") or 0
    n_listed = ws.get("n_listed") or 0

    head = dmc.Group(
        [
            dmc.Text(ws["label"], fw=600, size="sm"),
            dmc.Badge(meaning, color=color, variant="filled", size="sm"),
            dmc.Text(f"{n_loaded} of {n_listed} sample(s)", size="xs",
                     c="dimmed"),
        ],
        gap="xs", wrap="wrap",
    )

    pills = []
    # Only the routes that actually contributed. A row of zeroes is noise.
    if routes.get("name"):
        pills.append(_pill("by name", routes["name"], "green"))
    if routes.get("recorded"):
        pills.append(_pill("by recorded name", routes["recorded"], "teal"))
    if routes.get("content"):
        pills.append(_pill("by $FIL content", routes["content"], "blue"))
    if ws.get("n_gates"):
        pills.append(_pill("gates", ws["n_gates"]))
    if ws.get("markers"):
        pills.append(_pill("markers", len(ws["markers"])))

    body = [head]
    if pills:
        body.append(dmc.Group(pills, gap=4, wrap="wrap"))

    if ws.get("reason"):
        body.append(
            dmc.Text(ws["reason"], size="xs", c="red",
                     style={"whiteSpace": "pre-wrap"}))

    if ws.get("markers"):
        body.append(dmc.Text("Panel: " + ", ".join(ws["markers"]),
                             size="xs", c="dimmed"))
    if ws.get("markers_not_shared"):
        body.append(dmc.Text(
            "Not shared with every workspace: "
            + ", ".join(ws["markers_not_shared"]), size="xs", c="orange"))

    for title, key, colour in (
        ("Paired by content, not by name", "rescued", "blue"),
        ("Listed but did not load", "rejected", "red"),
        ("Loaded but never listed", "unlisted", "orange"),
    ):
        block = _sample_list(title, ws.get(key) or [], colour)
        if block is not None:
            body.append(block)

    return dmc.Paper(dmc.Stack(body, gap=6), p="xs", withBorder=True,
                     radius="sm")


def _headline(report: dict) -> dmc.Alert:
    """The one paragraph someone reads before deciding whether to continue."""
    totals = report.get("totals") or {}
    found = totals.get("workspaces_found") or 0
    loaded = totals.get("workspaces_loaded") or 0
    samples = totals.get("samples") or 0
    listed = totals.get("samples_listed") or 0
    routes = totals.get("routes") or {}

    lines = [
        dmc.Text(f"{loaded} of {found} workspace(s) loaded — "
                 f"{samples} of {listed} sample(s) are in the analysis.",
                 fw=600, size="sm"),
    ]
    if routes.get("content"):
        lines.append(dmc.Text(
            f"{routes['content']} of them were paired by the $FIL keyword "
            f"inside the file rather than by file name. Those would be missing "
            f"under strict matching.", size="xs"))
    if totals.get("samples_without_marker_names"):
        lines.append(dmc.Text(
            f"{totals['samples_without_marker_names']} sample(s) carry no "
            f"marker names of their own; their markers were recovered from the "
            f"channel each was measured on.", size="xs"))

    colour = "green" if loaded == found else ("yellow" if loaded else "red")
    return dmc.Alert(dmc.Stack(lines, gap=4), color=colour, variant="light")


def _marker_section(report: dict):
    """Common panel, what diverges, and the aliases R suggests for it."""
    markers = report.get("markers") or {}
    common = markers.get("common") or []
    not_common = markers.get("not_common") or []
    suggestions = markers.get("suggestions") or []
    conflicts = markers.get("conflicts") or []

    body = [
        dmc.Text(f"Markers shared by every workspace ({len(common)})",
                 fw=600, size="sm"),
        dmc.Text(", ".join(common) if common
                 else "None — no joint UMAP is possible until this is fixed.",
                 size="xs", c="dimmed" if common else "red"),
    ]

    if not_common:
        body.append(dmc.Space(h=6))
        body.append(dmc.Text(f"Not shared ({len(not_common)})", fw=600,
                             size="sm"))
        body.append(dmc.Stack(
            [dmc.Text(f"{m['marker']} — in {len(m['present_in'])}: "
                      + ", ".join(m["present_in"][:6])
                      + ("…" if len(m["present_in"]) > 6 else ""),
                      size="xs")
             for m in not_common],
            gap=2))

    if conflicts:
        body.append(dmc.Space(h=6))
        body.append(dmc.Alert(
            dmc.Stack(
                [dmc.Text("A detector is named differently in different files. "
                          "The most common name was used:", size="xs")]
                + [dmc.Text(f"{c['channel']}: {c['counts']} → {c['chosen']}",
                            size="xs", ff="monospace")
                   for c in conflicts],
                gap=2),
            color="yellow", variant="light"))

    if suggestions:
        body.append(dmc.Space(h=6))
        body.append(dmc.Text("Probably the same marker, spelled differently",
                             fw=600, size="sm"))
        body.append(dmc.Text(
            "Nothing is merged unless you say so. Ticking these fills in the "
            "Marker aliases field under Advanced; load the workspaces again to "
            "apply them.", size="xs", c="dimmed"))
        # Pattern-matching id, because this group only exists once a report has
        # been rendered. A callback keyed on a plain id would fail on every page
        # load before the first inspect; a pattern with zero matches is fine.
        body.append(dmc.CheckboxGroup(
            id={"comp": "alias-suggestions"},
            value=[],
            children=dmc.Stack(
                [dmc.Checkbox(
                    label=f"{' = '.join(s['names'])}  →  {s['suggested']}"
                          f"   ({s['confidence']})",
                    # The value IS the CLI fragment, so nothing has to rebuild
                    # it later: every non-chosen spelling maps to the chosen one.
                    value=";".join(f"{n}={s['suggested']}"
                                   for n in s["names"]
                                   if n != s["suggested"]),
                    size="sm")
                 for s in suggestions],
                gap=4),
        ))
    return dmc.Paper(dmc.Stack(body, gap=4), p="sm", withBorder=True,
                     radius="sm")


def badge(inspect: dict | None):
    """The one-glance verdict shown on the collapsed accordion header."""
    report = (inspect or {}).get("report")
    if not report:
        return None
    totals = report.get("totals") or {}
    found = totals.get("workspaces_found") or 0
    loaded = totals.get("workspaces_loaded") or 0
    if loaded == found and found:
        return dmc.Badge(f"{found}/{found} workspaces", color="green",
                         variant="light", size="sm")
    return dmc.Badge(f"{loaded}/{found} workspaces",
                     color="red" if not loaded else "yellow",
                     variant="filled", size="sm")


def panel(inspect: dict | None) -> object:
    """The whole health panel, or a placeholder before anything is loaded."""
    report = (inspect or {}).get("report")
    if not report:
        return dmc.Text(
            "Load the workspaces to see what came in and what did not.",
            size="sm", c="dimmed")

    workspaces = report.get("workspaces") or []
    # Problems first. Someone opening this panel is looking for what went wrong,
    # and making them scroll past twenty successes to find it is the whole
    # failure mode of the text log this replaces.
    workspaces = sorted(workspaces,
                        key=lambda w: (w.get("status") == "loaded",
                                       w.get("label", "")))

    return dmc.Stack(
        [
            _headline(report),
            _marker_section(report),
            dmc.Text(f"Workspaces ({len(workspaces)})", fw=600, size="sm"),
            dmc.Stack([_workspace_row(w) for w in workspaces], gap=6),
        ],
        gap="sm",
    )
